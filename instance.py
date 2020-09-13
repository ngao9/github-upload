import subprocess
import pandas as pd
import numpy as np
import ast
import json
import copy
import folium
import re
import matplotlib.pyplot as plt

import plotdesign
import instance_validator
import taz_processor

# Models an ODMTS instance.
# Keeps track of all parameters, allows for running the optimization, and for plotting the results.
class Instance:

    __runFile = None     # file path to executable
    __parameters = None  # dictionary of command line flags and values
    __logFile = None     # file path for writing the optimization log
    __plotter = None     # plotdesign object (see plotdesign.py)
    __hist = None        # matplotlib figure with histogram
    output = None        # dictionary based on the json output of the optimization
    log = None           # contents of the optimization log

    def __init__(self, runFile, parameters, logFile):
        self.__runFile = runFile
        self.__parameters = copy.deepcopy(parameters)
        self.__logFile = logFile

    # Run the optimization, write the log, and load the optimization output
    def run(self, wsl=False):
  
        command = self.__getCommand(wsl)
        try:
            log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(str(e.output, encoding="utf-8").replace('\\n', '\n').replace('\\t', '\t'))
            raise

        with open(self.__logFile, 'wb') as log_file:
            log_file.write(log)
        log_file.close()
        self.log = log.decode('UTF-8')

        self.load_outputFile()
        
    # Prepare command for running the executable
    # If wsl=True, run on Windows Subsystem Linux
    def __getCommand(self, wsl):
        
        command = [self.__runFile]
        command.extend(self.__getParametersAsList())
           
        if wsl:
            # Remove .exe extension
            if command[0][-4:] == '.exe':
                command[0] = command[0][:-4]
            
            # Map to WSL folders (change D:/ to /mnt/d/)
            for i in range(len(command)):
                string = command[i]
                if type(string) == str:
                    if re.search("^([a-zA-Z]):", string) != None: #match D:/
                        string = "/mnt/" + string[0:1].lower() + string[2:]
                        command[i] = string
                        
            # Add wsl command
            command = ['wsl'] + command
            
        return command

    # Load outputFile from file
    # By default, the file in the parameters is loaded (the result of run())
    def load_outputFile(self, outputFile=None):

        if(outputFile == None):
            outputFile = self.getParameterValue('outputFile')

        with open(outputFile) as output_file:
            self.output = json.load(output_file)
        output_file.close()

    # Write outputFile to file
    # By default, the file in the parameters is overwritten (the result of run())
    # If ignore_trip_splitting=True, then trip_splitting is ignored in the output
    def write_outputFile(self, outputFile=None, ignore_trip_splitting=False):

        if(outputFile == None):
            outputFile = self.getParameterValue('outputFile')

        has_trip_splitting = 'trip_splitting' in self.output

        if has_trip_splitting and ignore_trip_splitting:
            trip_splitting = self.output['trip_splitting']
            self.output['trip_splitting'] = None

        with open(outputFile, 'w') as file:
            if has_trip_splitting and not ignore_trip_splitting:
                json.dump(self.output, file)
            else:
                json.dump(self.output, file, indent=4)

        if has_trip_splitting and ignore_trip_splitting:
            self.output['trip_splitting'] = trip_splitting

    # Replace output['trips'] by the od-pairs in the sample, routed using the trip splittings
    # Note that output['scores'], output['parameters'], etc. are not updated
    def replace_trips_by_sample(self, sampleFile):

        if (not 'trip_splitting' in self.output) or (self.output['trip_splitting'] == None):
            print("Warning: cannot replace trips by sample, as output has no trip splitting information. " +\
                  "Make sure the output is generated with generateTripSplittings=True.")
            return

        # Read distance and duration matrices
        distances = pd.read_csv(self.getParameterValue('roadTravelDistancesFile'))
        distances.rename(columns={distances.columns[0]: 'stop_id'}, inplace=True)
        distances['stop_id'] = distances['stop_id'].astype(int).astype(str)
        distances.set_index('stop_id', inplace=True)
        distances.columns = [str(int(float(i))) for i in distances.columns]
        assert (distances.index == distances.columns).all()  # consistent row and column labels
        distances_rounded = distances.applymap(lambda x: int(np.ceil(1000*x)))
        
        durations = pd.read_csv(self.getParameterValue('roadTravelTimesFile'))
        durations.rename(columns={durations.columns[0]: 'stop_id'}, inplace=True)
        durations['stop_id'] = durations['stop_id'].astype(int).astype(str)  
        durations.set_index('stop_id', inplace=True)
        durations.columns = [str(int(float(i))) for i in durations.columns]
        assert (durations.index == durations.columns).all()  # consistent row and column labels
        
        shuttle = self.output['parameters']['modes'].index('shuttle')
        durations_rounded = dict()
        durations_rounded[shuttle] = durations.applymap(lambda x:
                        int(np.ceil(self.getParameterValue('travelTimeFactorShuttle')*60*x)))

        # Use trip splitting to replace trips by sampled trips
        samples = pd.read_csv(sampleFile,
                              dtype={'start_cluster': str, 'end_cluster': str, 'count': int, 'start_times': str},
                              index_col=0).reset_index(drop=True)
        samples['start_times'] = samples['start_times'].apply(ast.literal_eval)

        def passenger_cost(leg):
            param = self.output['parameters']
            if param['modes'][leg['mode']] == 'shuttle':
                cost = distances.loc[leg['board_stop_id'], leg['alight_stop_id']]
                cost *= param['shuttleCostPerKm']
            else:
                cost = 0
            return cost

        def passenger_convenience(leg):
            param = self.output['parameters']
            if param['modes'][leg['mode']] == 'shuttle':
                convenience = durations.loc[leg['board_stop_id'], leg['alight_stop_id']]
                convenience *= self.getParameterValue('travelTimeFactorShuttle')
            else:
                if param['modes'][leg['mode']] == 'bus':
                    convenience = durations.loc[leg['board_stop_id'], leg['alight_stop_id']]
                    convenience *= self.getParameterValue('travelTimeFactorBus')
                else:
                    convenience = leg['travel_time']/60 #assuming rail travel times have not been rounded
                convenience += param['fixedTransferTime']/60.0
                convenience += 60/(2.0 * leg['frequency_per_hour'])
            return convenience

        def passenger_objective(leg):
            alpha = self.output['parameters']['alpha']
            return (1 - alpha) * passenger_cost(leg) + alpha * passenger_convenience(leg)

        trips = list()

        for index, sample in samples.iterrows():

            trip = dict()
            trip['passengers'] = int(sample['count'])
            trip['origin_stop_id'] = sample['start_cluster']
            trip['destination_stop_id'] = sample['end_cluster']
            trip['departure_times'] = sample['start_times']
            trip['travel_time'] = 0
            trip['waiting_time'] = 0
            trip['total_time'] = 0
            trip['distance'] = 0
            trip['objective'] = 0
            trip['cost'] = 0
            trip['convenience'] = 0

            trip['legs'] = list()
            for leg in self.output['trip_splitting'][trip['origin_stop_id']][trip['destination_stop_id']]['legs']:
                if 'i' in leg:
                    trip['legs'].append(self.output['design']['legs'][leg['i']].copy())
                    trip['legs'][-1]['objective'] = None
                    trip['legs'][-1]['design_leg_index'] = leg['i']
                    del trip['legs'][-1]['fixed']
                else:
                    trip['legs'].append(
                        {'board_stop_id': leg['b'],
                         'alight_stop_id': leg['a'],
                         'mode': 0,
                         'frequency_per_hour': 0,
                         'frequency_per_hour_index': None,
                         'routes': None,
                         'travel_time': int(durations_rounded[shuttle].loc[leg['b'], leg['a']]),
                         'waiting_time': 0,
                         'total_time': int(durations_rounded[shuttle].loc[leg['b'], leg['a']]),
                         'distance': int(distances_rounded.loc[leg['b'], leg['a']]),
                         'objective': None,
                         'design_leg_index': None
                         }
                    )

                trip['legs'][-1]['objective'] = passenger_objective(trip['legs'][-1])
                trip['legs'][-1]['cost'] = passenger_cost(trip['legs'][-1])
                trip['legs'][-1]['convenience'] = passenger_convenience(trip['legs'][-1])

                trip['travel_time'] += trip['legs'][-1]['travel_time']
                trip['waiting_time'] += trip['legs'][-1]['waiting_time']
                trip['total_time'] += trip['legs'][-1]['total_time']
                trip['distance'] += trip['legs'][-1]['distance']
                trip['objective'] += trip['legs'][-1]['objective']
                trip['cost'] += trip['legs'][-1]['cost']
                trip['convenience'] += trip['legs'][-1]['convenience']

            trips.append(trip)

        self.output['trips'] = trips
        
        # Update scores
        scores = self.output['scores']
        scores['passenger_objective'] = 0
        scores['passenger_cost'] = 0
        scores['passenger_convenience']= 0
        scores['original_number_of_passengers'] = 0
    
        for trip in trips:
            
            for leg in trip['legs']:
                scores['passenger_objective'] += (leg['objective'] * trip['passengers'])
                scores['passenger_cost'] += (leg['cost'] * trip['passengers'])
                scores['passenger_convenience'] += (leg['convenience'] * trip['passengers'])
            
            scores['original_number_of_passengers'] += trip['passengers']

        scores['scaled_number_of_passengers'] = scores['original_number_of_passengers']
        scores['original_number_of_passengers'] = int(np.round(scores['original_number_of_passengers']))       
        
        scores['total_objective'] = scores['design_objective'] + scores['passenger_objective']
        scores['total_cost'] = scores['design_cost'] + scores['passenger_cost']
        scores['total_convenience'] = scores['design_convenience'] + scores['passenger_convenience']
        scores['total_cost_per_passenger'] = \
            (scores['design_cost'] + scores['passenger_cost'])/scores['scaled_number_of_passengers']
        scores['passenger_cost_per_passenger'] = \
              scores['passenger_cost']/scores['scaled_number_of_passengers']
        scores['passenger_convenience_per_passenger'] = \
              scores['passenger_convenience']/scores['scaled_number_of_passengers']
        
        self.__parameters['passengerFactor'] = 1.0
        self.output['parameters']['passengerFactor'] = 1.0
        
    def __getParametersAsList(self):

        def convert_key_to_param_list(key):

            value = self.__parameters[key]

            if type(value) == str:
                return ['--' + key, value]

            elif type(value) == int or type(value) == float:
                return ['--' + key, str(value)]

            elif type(value) == bool:
                if value == True:
                    return ['--' + key]
                else:
                    return []

            elif type(value) == tuple or type(value) == list:
                return ['--' + key] + [str(i) for i in value]

            else:
                assert (False)  # Unknown parameter value type

        param_list = list()
        for key in self.__parameters:
            param_list.extend(convert_key_to_param_list(key))

        return param_list

    # Print the command that is executed on run
    def printCommand(self, wsl=False):

        command = self.__getCommand(wsl)
        command_str = '"' + '" "'.join(command) + '"'
        print(command_str)

    # Get the parameter value for a given flag
    def getParameterValue(self, flag):
        return self.__parameters[flag]

    # Create a plot and store in plotter (see plotdesign.createPlot())
    def createPlot(self, choropleth=None, bus_weights=None, rail_weight=None, point_weight=None):
        self.__plotter = plotdesign.PlotDesign(self.getParameterValue('stopFile'), self.output)
        return self.__plotter.createPlot(choropleth, bus_weights, rail_weight, point_weight)

    def prepareTripChoropleth(self,
                              trip_value = lambda trip: 0,
                              trip_stop_id = lambda trip: trip['origin_stop_id'],
                              trip_weight = lambda trip: trip['passengers'],
                              processor = None,
                              taz_shape_file='taz-2010.shp'):

        # Initialize TAZ processor if necessary
        if processor == None:
            processor = taz_processor.TazProcessor(taz_shape_file)

        # Get latlons for stop_ids
        locations = pd.read_csv(self.getParameterValue('stopFile'),
                                dtype={'stop_id': 'int',
                                        'stop_lat': 'float',
                                        'stop_lon': 'float'})
        locations['stop_id'] = locations['stop_id'].astype(str)
        locations.set_index('stop_id', inplace=True)

        # Create trips dataframe with values per trip
        trips = pd.DataFrame(
            [{'stop_id': trip_stop_id(trip),
              'value': trip_value(trip),
              'weight': trip_weight(trip)} for trip in self.output['trips']]
        )
        trips['weighted_value'] = trips['weight'] * trips['value']

        # Efficiently add TAZs to every trip
        stop_id_to_taz = {stop_id: processor.location_to_taz(locations.loc[stop_id, 'stop_lat'],
                                                            locations.loc[stop_id, 'stop_lon'])
                        for stop_id in trips['stop_id'].unique()}
        trips['taz'] = trips['stop_id'].apply(lambda x: stop_id_to_taz[x])

        # Calculate score per TAZ
        grouped_trips = trips.groupby('taz').sum()
        grouped_trips['score'] = grouped_trips['weighted_value']/grouped_trips['weight']
        grouped_trips['taz'] = grouped_trips.index #needed for choropleth

        # Convert TAZ shapes to geo_data format
        geo_data = []
        for row in grouped_trips.itertuples():
            taz_oid = row.Index
            geo_data.append({'type': 'Feature',
                            'id': taz_oid,
                            'geometry': processor.taz[taz_oid].shape.__geo_interface__})
        draw_geo_data = {'type': 'FeatureCollection', 'features': geo_data}

        # Prepare choropleth data
        choropleth = folium.Choropleth( geo_data = draw_geo_data,
                                        name = 'choropleth',
                                        key_on = 'id',
                                        data = grouped_trips,
                                        columns = ['taz','score'],
                                        fill_color = 'BuPu',
                                        fill_opacity = 0.7,
                                        line_opacity = 0.2 )
        return choropleth

    # Export created plot to html and png (see plotdesign.exportPlot())
    def exportPlot(self, htmlFile, pngFile=None):
        return self.__plotter.exportPlot(htmlFile, pngFile)

    def createHistogram(self, hist_type='transfers'):
        
        if hist_type == 'transfers':

            trips = pd.DataFrame([trip for trip in self.output['trips']])
            trips = trips.explode(column='departure_times').reset_index(drop=True)
            nb_transfers = trips['legs'].apply(lambda x: len(x)-1)
            max_nb_transfers = self.output['parameters']['maximum_number_of_transfers']

            self.__hist, ax = plt.subplots()
            nb_transfers.hist(bins=np.linspace(-0.5, max_nb_transfers + 0.5, max_nb_transfers + 2),
                            rwidth=0.5, figsize=(2*max_nb_transfers,8), ax=ax)
            ax.set_xticks(range(max_nb_transfers + 1))
            ax.set_xlabel('Number of transfers')
            ax.set_ylabel('Number of passengers')
            ax.xaxis.grid(False)

        else:
            print("Error createHistogram: Unknown hist_type.")
            return

    def exportHist(self, pngFile):
        self.__hist.savefig(pngFile)

    # Verify that the ODMTS output is feasible
    # See Validator.validate() in instanc_validator.py for options
    def verifyOutput(self,
                     eps=1E-9,
                     schema_file="odmts_schema.json",
                     verbose=False,
                     schema_validation=True,
                     throw_exception=True,
                     backward_compatibility=False):

        if backward_compatibility:
            self.output['parameters']['travelTimeFactorBus'] = self.getParameterValue('travelTimeFactorShuttle')
            self.output['parameters']['travelTimeFactorShuttle'] = self.getParameterValue('travelTimeFactorShuttle')

            for leg in self.output['design']['legs']:
                leg['cost'] = 0
            
            for trip in self.output['trips']:
                trip['distance'] = int(np.round(trip['distance']))
                for leg in trip['legs']:
                    leg['cost'] = 0
                    leg['convenience'] = 0

        validator = instance_validator.Validator(self)

        try:
            validator.validate(eps=eps,
                               output_from_file=False,
                               schema_file=schema_file,
                               schema_validation=schema_validation,
                               ignore_rounding=True,
                               verbose=verbose,
                               backward_compatibility=backward_compatibility
                               )
            return True
        except AssertionError:
            if throw_exception:
                raise
            else:
                return False
