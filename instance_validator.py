import pandas as pd
import numpy as np
import json
import jsonschema
import instance


class Validator:

    __inst = None
    __verbose = None
    __input_trips = None
    __output = None
    __distances = None
    __distances_rounded = None
    __durations = None
    __durations_rounded = dict()
    __rail_trips = None
    __backward_compatibility = False

    def __init__(self, inst):
        self.__inst = inst

    # Validate the correctness of the instance
    # eps: epsilon used for comparison, automatically scaled up for aggregated values
    # output_from_file: True - read output from file defined in instance parameters
    #                   False - skip reloading output, use instance.output instead
    # schema_file: json schema file used for validation
    # schema_validation: True - use jsonschema to validate the json file
    #                           skip validating the trip_splitting data with the schema (which is very slow)
    #                    'extensive' - include trip_splitting validation
    #                    False - no validation
    # ignore_rounding: True - ignore off by 1 errors in distances and durations due to inconsistent rounding
    # verbose: True - print information
    # backward_compatibility: skip score validation for backward compatibility design files (pre June 2020)
    def validate(self,
                 eps=1E-9,
                 output_from_file=True,
                 schema_file="odmts_schema.json",
                 schema_validation=True,
                 ignore_rounding=False,
                 verbose=True,
                 backward_compatibility=False):
                
        self.__verbose = verbose
        if backward_compatibility:
            self.__backward_compatibility = True
            self.__message('Running validator in backward compatibility mode.')

        self.__message('Loading output into validator...')
        if not output_from_file:
            self.__message('Using output provided by the instance.')
            if self.__inst.output == None:
                print('Exception: instance does not contain output.')
                assert False
            else:
                self.__output = self.__inst.output
        else:
            self.__message('Reading output from filepath in the instance parameters...')
            outputFile = self.getParameterValue('outputFile')
            with open(outputFile) as output_file:
                self.__output = json.load(output_file)
            output_file.close()

        self.__validate_schema(schema_file, schema_validation)

        self.__load_from_file()

        self.__validate_parameters(eps=eps)
        self.__validate_design(eps=eps, ignore_rounding=ignore_rounding)
        self.__validate_routes(eps=eps)
        self.__validate_trips(eps=eps, ignore_rounding=ignore_rounding)
        self.__validate_trip_splitting()
        if self.__backward_compatibility:
            self.__message('Skip validating scores in backward compatibility mode.')
        else:
            self.__validate_scores(eps=eps)

        self.__message('Validation completed.')

    def __validate_schema(self, schema_file, schema_validation):

        if not schema_validation:
            return

        with open(schema_file) as file:
            schema = json.load(file)

        trip_splitting = self.__output['trip_splitting']
        if schema_validation == True:
            self.__message('Validating schema without trip_splitting...')
            self.__output['trip_splitting'] = None
        else:
            self.__message('Validating schema with trip_splitting...')

        jsonschema.Draft7Validator(schema).validate(self.__output)

        self.__output['trip_splitting'] = trip_splitting

    def __validate_parameters(self, eps):

        self.__message("Comparing instance parameters with values reported by the output...")
        def approxEqual(x, y, eps):
            return np.max(np.abs(np.array(x) - np.array(y))) < eps
        
        def compare(output_key, param_flag=None):
            if param_flag == None:
                param_flag = output_key
            return approxEqual(self.__output['parameters'][output_key],
                               self.getParameterValue(param_flag),
                               eps)

        assert compare('alpha')
        assert compare('shuttleCostPerKm')
        assert approxEqual( self.__output['parameters']['shuttleCostPerMile'],
                            self.getParameterValue('shuttleCostPerKm') * 1.609344,
                            eps)
        assert compare('busCostPerHour')
        assert compare('fixedTransferTime')
        assert compare('passengerFactor')
        assert compare('frequencies_per_hour', 'busesPerHour')
        assert compare('maximum_number_of_transfers', 'maximumNumberOfTransfers')
        assert approxEqual( self.__output['parameters']['time_horizon'],
                            self.getParameterValue('timeHorizonInHours') * 3600,
                            eps)
        assert self.__output['parameters']['time_units'] == 'seconds'
        assert self.__output['parameters']['distance_units'] == 'meters'
        assert compare('travelTimeFactorShuttle')
        assert compare('travelTimeFactorBus')

        assert 'shuttle' in self.__output['parameters']['modes']
        assert 'bus' in self.__output['parameters']['modes']
        assert 'rail' in self.__output['parameters']['modes']
        assert len(self.__output['parameters']['modes']) == 3

    def __validate_design(self, eps, ignore_rounding):
        
        self.__message("Validating the design...")

        output = self.__output
        hubs = output['design']['hubs']
        legs = pd.DataFrame(output['design']['legs'])
        rail_trips = self.__rail_trips
        durations_rounded = self.__durations_rounded
        distances_rounded = self.__distances_rounded

        # Check for duplicate information
        assert len(legs) == len(legs.drop_duplicates(subset=['board_stop_id','alight_stop_id','mode']))

        # Check for NaNs
        columns = [column for column in legs.columns if column != 'frequency_per_hour_index']
        assert len(legs) == len(legs.dropna(subset=columns))

        # Check that all rail legs are in rail trips
        rail = output['parameters']['modes'].index('rail')
        rail_legs = legs[legs['mode'] == rail]
        rail_legs = rail_legs.merge(rail_trips,
                                    how='left',
                                    left_on=['board_stop_id', 'alight_stop_id'],
                                    right_on=['from_id','to_id'],
                                    validate='1:1'
        )
        columns = [column for column in legs.columns if column != 'frequency_per_hour_index']
        assert len(rail_legs) == len(rail_legs.dropna(subset=columns)) # all lookups succesfull
    
        # Compare rail leg travel time and frequency information
        assert (rail_legs['frequency_per_hour'] - rail_legs['rail_frequency']).abs().max() < eps
        assert (rail_legs['travel_time'] - rail_legs['time'].apply(lambda x: 60*x)).abs().max() < eps

        # Verify no shuttle legs in the design
        shuttle = output['parameters']['modes'].index('shuttle')
        assert (legs['mode'] == shuttle).sum() == 0

        # Verify bus frequency lookup is correct
        bus = output['parameters']['modes'].index('bus')
        bus_legs = legs[legs['mode'] == bus]
        assert(len(bus_legs) == len(bus_legs.dropna()))
        if len(bus_legs) > 0:
            bus_legs['freq_check'] = bus_legs['frequency_per_hour_index'].apply(
                lambda x: output['parameters']['frequencies_per_hour'][int(np.round(x))])
            assert (bus_legs['frequency_per_hour'] - bus_legs['freq_check']).abs().sum() < eps

        # Verify bus travel times
        if len(bus_legs) > 0:
            bus_legs['travel_time_check'] = bus_legs.apply(lambda row:
                    durations_rounded[bus].loc[row['board_stop_id'], row['alight_stop_id']],
                axis=1)
            violation = (bus_legs['travel_time_check'] - bus_legs['travel_time']).abs().max()
            if ignore_rounding:
                assert violation <= 1
            else:
                if violation == 1:
                    assert False, "Exception: rounding error in travel times, set ignore_rounding=True to ignore."
                else:
                    assert violation == 0

        # Verify bus distances
        if len(bus_legs) > 0:
            bus_legs['distance_check'] = bus_legs.apply(lambda row:
                distances_rounded.loc[row['board_stop_id'], row['alight_stop_id']],
                axis=1)
            violation = (bus_legs['distance_check'] - bus_legs['distance']).abs().max()
            if ignore_rounding:
                assert violation <= 1
            else:
                if violation == 1:
                    assert False, "Exception: rounding error in distances, set ignore_rounding=True to ignore."
                else:
                    assert violation == 0    

        # Verify bus frequency balance
        if len(bus_legs) > 0:
            freq_out = bus_legs[['board_stop_id', 'frequency_per_hour']].groupby(by=['board_stop_id']).sum()
            freq_in = bus_legs[['alight_stop_id', 'frequency_per_hour']].groupby(by=['alight_stop_id']).sum()
            balance = freq_in.merge(freq_out, how='outer',
                                    left_index=True, right_index=True,
                                    validate='1:1')
            assert (balance.iloc[:,0] - balance.iloc[:,1]).abs().max() < eps

        # Verify used hubs are open
        def verify_hub_open(hub):
            assert hub in hubs, "Hub " + str(hub) + " does not exist."
            assert hubs[hub]['opened'], "Hub " + str(hub) + " is not open."
            return True
        legs['board_stop_id'].apply(verify_hub_open)
        legs['alight_stop_id'].apply(verify_hub_open)

        # Verify total_time and waiting_time
        assert (legs['total_time'] == legs['travel_time'] + legs['waiting_time']).all()
        assert (legs['frequency_per_hour'].apply(lambda x: 60 * (60/2)/x) +
                output['parameters']['fixedTransferTime'] -
                legs['waiting_time']).abs().max() < eps

    def __validate_routes(self, eps):
        
        self.__message("Validating routes...")
        
        output = self.__output
        durations = self.__durations
        durations_rounded = self.__durations_rounded
        legs = pd.DataFrame(output['design']['legs'])

        legs = pd.DataFrame(output['design']['legs'])
        routes = output['routes']

        # Convert routes to legs
        leg_freq_names = dict()
        for route_name in routes:

            route = routes[route_name]
            stop_ids = route['stop_ids']
            
            for i in range(len(stop_ids)):
                
                from_id = stop_ids[i]
                to_id = stop_ids[(i+1) % len(stop_ids)]
                
                key = (from_id, to_id, route['mode'])
                leg_freq_names[key] = (leg_freq_names.get(key, [0, 0])[0] + route['frequency_per_hour'],
                                    leg_freq_names.get(key, [[],[]])[1] + [route_name])

        # Add route information for checking
        legs[['freq_check','name_check']] = legs.apply(lambda row:
            leg_freq_names.get((row['board_stop_id'], row['alight_stop_id'], row['mode']), (None, None)),
            axis=1, result_type='expand')

        # Check that all bus legs are covered
        bus = output['parameters']['modes'].index('bus')
        bus_legs = legs[legs['mode'] == bus]
        assert len(bus_legs) == sum([1 for odm in leg_freq_names if odm[2] == bus])
        assert bus_legs['name_check'].isnull().sum() == 0

        # Check leg and route consistency
        assert (legs['frequency_per_hour'] - legs['freq_check']).abs().sum() < eps
        nonnan = legs['name_check'].isnull() == False
        assert (legs.loc[nonnan,'routes'].apply(lambda x: set(x)) == \
                    legs.loc[nonnan,'name_check'].apply(lambda x: set(x))).all()

    def __validate_trips(self, eps, ignore_rounding):
        
        self.__message("Validating trips...")
        
        output = self.__output
        alpha = output['parameters']['alpha']
        trips = output['trips']
        shuttle = output['parameters']['modes'].index('shuttle')
        durations_rounded = self.__durations_rounded
        distances_rounded = self.__distances_rounded

        approxEqual = lambda x, y, eps: np.abs(x - y) < eps
        def compare(x, y):
            violation = np.abs(x-y)
            if ignore_rounding:
                assert violation <= 1
            else:
                if violation == 1:
                    assert False, "Exception: rounding error in travel times, set ignore_rounding=True to ignore."
                else:
                    assert violation == 0

        for trip in trips:
            
            assert trip['origin_stop_id'] == trip['legs'][0]['board_stop_id']
            assert trip['destination_stop_id'] == trip['legs'][-1]['alight_stop_id']
            assert trip['passengers'] == len(trip['departure_times'])
            assert trip['travel_time'] == \
                sum([leg['travel_time'] for leg in trip['legs']])
            assert trip['waiting_time'] == \
                sum([leg['waiting_time'] for leg in trip['legs']])
            assert trip['total_time'] == \
                sum([leg['total_time'] for leg in trip['legs']])
            assert trip['distance'] == \
                sum([leg['distance'] for leg in trip['legs']])
            assert approxEqual(trip['objective'],
                            sum([leg['objective'] for leg in trip['legs']]),
                            eps)
            assert len(trip['legs']) <= output['parameters']['maximum_number_of_transfers'] + 1
            for i in range(len(trip['legs'])-1):
                assert trip['legs'][i]['alight_stop_id'] == trip['legs'][i+1]['board_stop_id']
            
            for leg in trip['legs']:
                
                if leg['mode'] == shuttle:
                    compare(leg['travel_time'],
                        durations_rounded[shuttle].loc[leg['board_stop_id'], leg['alight_stop_id']])
                    assert leg['waiting_time'] == 0
                    assert leg['total_time'] == leg['travel_time'] + leg['waiting_time']
                    compare(leg['distance'],
                        distances_rounded.loc[leg['board_stop_id'], leg['alight_stop_id']])        
                else:
                    leg = output['design']['legs'][leg['design_leg_index']]
                    
                    equals_keys = ['board_stop_id',
                                'alight_stop_id',
                                'mode',                          
                                'frequency_per_hour_index',
                                'travel_time',
                                'waiting_time',
                                'total_time',
                                'distance'
                                ]
                    for key in equals_keys:
                        assert leg[key] == leg[key]
                        
                    assert approxEqual(leg['frequency_per_hour'],
                                    leg['frequency_per_hour'], eps)
                    assert set(leg['routes']) == set(leg['routes'])

    def __validate_trip_splitting(self):
                       
        output = self.__output
        trip_splitting = output['trip_splitting']
        trips = output['trips']
        
        if trip_splitting == None:
            self.__message("No trip splitting file provided for validation.")
            return
        else:
            self.__message("Validating trip splitting...")
        
        # Verify trip_splitting is a square matrix
        keys_1 = set(trip_splitting.keys())
        for key in trip_splitting:
            keys_2 = set(trip_splitting[key].keys())
            assert keys_1 == keys_2
            
        # Verify that trip information is consistent with trip_splitting
        for trip in trips:
            
            origin = trip['origin_stop_id']
            destination = trip['destination_stop_id']
            
            trip_split = trip_splitting[origin][destination]
            assert len(trip['legs']) == len(trip_split['legs'])
            
            for i in range(len(trip['legs'])):
            
                leg = trip['legs'][i]
            
                if leg['design_leg_index'] == None:
                    assert trip_split['legs'][i]['b'] == leg['board_stop_id']
                    assert trip_split['legs'][i]['a'] == leg['alight_stop_id']
                    assert trip_split['legs'][i]['m'] == leg['mode']
                else:
                    assert trip_split['legs'][i]['i'] == leg['design_leg_index']
        
    def __validate_scores(self, eps):

        self.__message("Validating scores...")
        
        output = self.__output
        scores = output['scores']
        alpha = output['parameters']['alpha']
        trips = output['trips']
        shuttle = output['parameters']['modes'].index('shuttle')
        bus = output['parameters']['modes'].index('bus')
        rail = output['parameters']['modes'].index('rail')
        rail_trips = self.__rail_trips
        distances = self.__distances
        durations = self.__durations

        approxEqual = lambda x, y, eps: np.abs(x - y) < eps

        # Validate passenger scores in the trips, and add up the totals
        passenger_objective = 0
        passenger_cost = 0
        passenger_convenience= 0
        original_number_of_passengers = 0
    
        for trip in trips:
            for leg in trip['legs']:
                
                if leg['mode'] == shuttle:
                    cost = output['parameters']['shuttleCostPerKm'] * \
                            distances.loc[leg['board_stop_id'], leg['alight_stop_id']]
                    convenience = output['parameters']['travelTimeFactorShuttle'] * \
                            durations.loc[leg['board_stop_id'], leg['alight_stop_id']]
                    objective = (1-alpha) * cost + alpha * convenience              
                else:
                    cost = 0
                    if leg['mode'] == bus:
                        convenience = durations.loc[leg['board_stop_id'], leg['alight_stop_id']]
                        convenience *= output['parameters']['travelTimeFactorBus']
                    else:
                        convenience = rail_trips.loc[
                            (rail_trips['from_id'] == leg['board_stop_id']) &
                            (rail_trips['to_id'] == leg['alight_stop_id']),
                            'time'].iloc[0]
                    convenience += output['parameters']['fixedTransferTime']/60.0
                    convenience += 60/(2.0 * leg['frequency_per_hour'])
                    objective = (1-alpha) * cost + alpha * convenience
                    
                assert approxEqual(leg['cost'], cost, eps)
                assert approxEqual(leg['convenience'], convenience, eps)
                assert approxEqual(leg['objective'], objective, eps)
                
                passenger_cost += cost * trip['passengers']
                passenger_convenience += convenience * trip['passengers']
                passenger_objective += objective * trip['passengers']
                
            original_number_of_passengers += trip['passengers']
        
        passenger_cost *= output['parameters']['passengerFactor']
        passenger_convenience *= output['parameters']['passengerFactor']
        passenger_objective *= output['parameters']['passengerFactor']
        
        # Validate design scores, and add up the totals
        design_objective = 0
        design_cost = 0
        design_convenience = 0
        
        for leg in output['design']['legs']:
            
            if leg['fixed']:
                assert leg['distance'] == 0, "Check if distance is set to 0 for fixed arcs, not strictly required."
                assert approxEqual(leg['cost'], 0, eps)
                assert approxEqual(leg['objective'], 0, eps)
            else:
                cost = durations.loc[leg['board_stop_id'], leg['alight_stop_id']]/60.0 #driving hours
                if leg['mode'] == bus:
                    cost *= output['parameters']['travelTimeFactorBus'] #time correction
                cost *= output['parameters']['time_horizon']/3600.0 #extend to time horizon
                cost *= leg['frequency_per_hour'] #correct by frequency
                cost *= output['parameters']['busCostPerHour'] #cost
                convenience = 0
                objective = (1-alpha) * cost + alpha * convenience
                
                assert approxEqual(leg['cost'], cost, eps)
                assert approxEqual(leg['objective'], objective, eps)
                
                design_cost += cost
                design_convenience += convenience
                design_objective += objective
            
        # Validate aggregate scores
        assert scores['original_number_of_passengers'] == original_number_of_passengers
        assert approxEqual(scores['scaled_number_of_passengers'],
                           original_number_of_passengers * output['parameters']['passengerFactor'], eps)
        
        
        eps_design = eps * len(output['design']['legs'])
        eps_passenger = eps * scores['scaled_number_of_passengers']
        
        assert approxEqual(scores['total_objective'], design_objective + passenger_objective, eps_design + eps_passenger)
        assert approxEqual(scores['total_cost'], design_cost + passenger_cost, eps_design + eps_passenger)
        assert approxEqual(scores['total_convenience'], design_convenience + passenger_convenience, eps_design + eps_passenger)
        assert approxEqual(scores['design_objective'], design_objective, eps_design)
        assert approxEqual(scores['design_cost'], design_cost, eps_design)
        assert approxEqual(scores['design_convenience'], design_convenience, eps_design)
        assert approxEqual(scores['passenger_objective'], passenger_objective, eps_passenger)
        assert approxEqual(scores['passenger_cost'], passenger_cost, eps_passenger)
        assert approxEqual(scores['passenger_convenience'], passenger_convenience, eps_passenger)
        assert approxEqual(scores['total_cost_per_passenger'],
              (design_cost + passenger_cost)/scores['scaled_number_of_passengers'], eps)
        assert approxEqual(scores['passenger_cost_per_passenger'],
              passenger_cost/scores['scaled_number_of_passengers'], eps)
        assert approxEqual(scores['passenger_convenience_per_passenger'],
              passenger_convenience/scores['scaled_number_of_passengers'], eps)
        
    def __load_from_file(self):

        # Read input trips file
        self.__message('Reading input trips from file...')
        self.__input_trips = pd.read_csv(self.getParameterValue('tripsFile'),
                                    dtype={'start_stop': int, 'end_stop': int, 'count': 'int'})
        self.__input_trips.rename(columns={'count': 'nb_people'}, inplace=True)

        # Read distance and duration matrices
        self.__message('Reading distances from file...')
        self.__distances = pd.read_csv(self.getParameterValue('roadTravelDistancesFile'))
        self.__distances.rename(columns={self.__distances.columns[0]: 'stop_id'}, inplace=True)
        self.__distances['stop_id'] = self.__distances['stop_id'].astype(int).astype(str)
        self.__distances.set_index('stop_id', inplace=True)
        self.__distances.columns = [str(int(float(i))) for i in self.__distances.columns]
        assert (self.__distances.index == self.__distances.columns).all()  # consistent row and column labels
        self.__distances_rounded = self.__distances.applymap(lambda x: int(np.ceil(1000*x)))

        self.__message('Reading durations from file...')
        self.__durations = pd.read_csv(self.getParameterValue('roadTravelTimesFile'))
        self.__durations.rename(columns={self.__durations.columns[0]: 'stop_id'}, inplace=True)
        self.__durations['stop_id'] = self.__durations['stop_id'].astype(int).astype(str)
        self.__durations.set_index('stop_id', inplace=True)
        self.__durations.columns = [str(int(float(i))) for i in self.__durations.columns]
        assert (self.__durations.index == self.__durations.columns).all()  # consistent row and column labels
        
        shuttle = self.__output['parameters']['modes'].index('shuttle')
        bus = self.__output['parameters']['modes'].index('bus')
        self.__durations_rounded[shuttle] = self.__durations.applymap(lambda x:
                        int(np.ceil(self.__output['parameters']['travelTimeFactorShuttle']*60*x)))
        self.__durations_rounded[bus] = self.__durations.applymap(lambda x:
                        int(np.ceil(self.__output['parameters']['travelTimeFactorBus']*60*x)))

        assert (self.__distances.index == self.__durations.index).all()  # consistency between distance and time matrices

        # Read and process rail line data
        self.__message('Reading rail lines from file and determining rail trips...')
        railLine = pd.read_csv(self.getParameterValue('railLineFile'),
                                dtype={'NAME': str, 'NB_PER_HOUR': int})

        freq_time = dict()

        for index, row in railLine.iterrows():
            nbPerHour = row['NB_PER_HOUR']
            if nbPerHour > 0:

                assert railLine.iloc[index + 1]['NAME'] == 'schedule_left_right'
                assert railLine.iloc[index + 2]['NAME'] == 'schedule_right_left'

                for i in range(2, len(row)):
                    for j in range(2, len(row)):

                        if i == j or pd.isnull(row.iloc[i]) or pd.isnull(row.iloc[j]):
                            continue

                        if i < j:
                            time = railLine.iloc[index + 1, j] - railLine.iloc[index + 1, i]
                        else:
                            time = railLine.iloc[index + 2, j] - railLine.iloc[index + 2, i]

                        key = (row.iloc[i], row.iloc[j])
                        if key in freq_time:
                            assert(freq_time[key][1] == time) # travel time is consistent
                            freq_time[key] = (freq_time[key][0] + nbPerHour, time)
                        else:
                            freq_time[key] = (nbPerHour, time)
        
        self.__rail_trips = pd.DataFrame([  { 'from_id': str(int(key[0])),
                                            'to_id': str(int(key[1])),
                                            'rail_frequency': freq_time[key][0],
                                            'time': freq_time[key][1]
                                          }
                                             for key in freq_time])

    def getParameterValue(self, flag):
        return self.__inst.getParameterValue(flag)

    def __message(self, text):
        if self.__verbose:
            print(text)