def add_descriptions(schema):

    schema['properties']['parameters']['description'] = "Overview of the most important parameters used in the optimization and in this file"
    parameters = schema['properties']['parameters']['properties'];
    parameters['alpha']['description'] = "Weight for cost and convenience (higher is more convenient)"
    parameters['shuttleCostPerKm']['description'] = "Shuttle cost per kilometer"
    parameters['shuttleCostPerMile']['description'] = "Shuttle cost per mile"
    parameters['busCostPerHour']['description'] = "Bus cost per hour"
    parameters['fixedTransferTime']['description'] = "Fixed transfer time in seconds at the start of each non-shuttle arc"
    parameters['passengerFactor']['description'] = "Multiply the number of passengers with this factor"
    parameters['frequencies_per_hour']['description'] = "Possible numbers of non-fixed vehicles (typically buses) per hour, only guaranteed to be integer over the time horizon (not per hour)"
    parameters['maximum_number_of_transfers']['description'] = "Maximum number of transfers"
    parameters['time_horizon']['description'] = "Time horizon"
    parameters['modes']['description'] = "Possible modes, typically [shuttle, bus, rail]"
    parameters['time_units']['description'] = "Time units used throughout this file, typically seconds"
    parameters['distance_units']['description'] = "Distance units used throughout this file, typically meters"
    parameters['travelTimeFactorShuttle']['description'] = "Multiplicationfactor to obtain shuttle travel times"
    parameters['travelTimeFactorBus']['description'] = "Multiplicationfactor to obtain bus travel times"

    schema['properties']['design']['description'] = "Details on the generated ODMTS design"
    design = schema['properties']['design']['properties']
    design['hubs']['description'] = "List of all hubs, with the stop id (string) as their key. Hubs are stops that may be used as waypoints to transport passengers"
    design_hub = design['hubs']['patternProperties'][r'^.*$']['properties']
    design_hub['fixed']['description'] = "Hub was fixed in the design, typically the case for rail stations"
    design_hub['opened']['description'] = "Hub is opened and can be used, typically true. It may be that the hub is exclusively used by shuttles and not by bus and rail"

    design['legs']['description'] = (
        "List of all legs in the design, from origin to destination. "
        "A leg gives a connection that can be traveled without any transfers, but there may be other stops inbetween. "
        "Typically, for a fixed rail line of length n, there will be n x (n-1) legs, so one between each pair of stops.")
    design_leg = design['legs']['items']['properties']
    design_leg['board_stop_id']['description'] = "Boarding stop id (string)"
    design_leg['alight_stop_id']['description'] = "Alighting stop id (string)"
    design_leg['mode']['description'] = "Mode of this leg (index, lookup mode name in /parameters/modes)"
    design_leg['frequency_per_hour']['description'] = "Numbers of vehicles per hour, only guaranteed to be integer over the time horizon (not per hour)"
    design_leg['frequency_per_hour_index']['description'] = "For non-fixed legs (typically buses), give the index of the frequency in /parameters/frequency_per_hour, otherwise null"
    design_leg['routes']['description'] = "List of possible routes to travel this leg, can be looked up in /routes"
    design_leg['fixed']['description'] = "Leg was fixed in the design, typically the case for rail legs"
    design_leg['travel_time']['description'] = "Time between boarding and alighting, rounded up to the next second"
    design_leg['waiting_time']['description'] = "Expected waiting time (including fixed waiting time), rounded up to the next second"
    design_leg['total_time']['description'] = "Sum of travel time and total time"
    design_leg['distance']['description'] = "Distance between boarding and alighting, may be set to zero for fixed legs, rounded up to the next meter"
    design_leg['objective']['description'] = "Contribution of design leg to the optimization objective value, may be set to zero for fixed legs"
    design_leg['cost']['description'] = "Unweighted contribution of design leg to the cost component of the objective value (using non-rounded distances and durations)"

    schema['properties']['routes']['description'] = (
        "List of routes that together cover all the design legs with the appropriate frequencies. "
        "Routes are cyclical: after the last stop, the vehicle returns to the first stop. "
        "Note that stops may occur multiple times on the same route, this is typical for rail routes that go back and forth."
    )
    routes = schema['properties']['routes']['patternProperties'][r'^.*$']['properties']
    routes['mode']['description'] = design_leg['mode']['description']
    routes['frequency_per_hour']['description'] = design_leg['frequency_per_hour']['description']
    routes['stop_ids']['description'] = "List of stop ids on the route (strings)"

    schema['properties']['trip_splitting']['anyOf'][0]['description'] = "If trip splitting generation was disabled (generateTripSplittings setting), then /tripsplitting is set to null"
    trip_splitting = schema['properties']['trip_splitting']['anyOf'][1]
    trip_splitting['description'] = (
        "Optimal passenger paths for all origins and destinations. "
        "The next two keys ('^.*$' and '^.*$') are first the origin stop id (string), and then the destination stop id (string). "
        "Both sets of keys are the same, but the ordering may differ. "
        "origin = destination is included, but these entries are set to null. "
    )
    trip_splitting_trip = trip_splitting['patternProperties'][r'^.*$']['patternProperties'][r'^.*$']['anyOf'][1]['properties']
    trip_splitting_trip['legs']['description'] = (
        "List of all legs traveled by the passenger, from origin to destination. "
        "Either 'b', 'a', and 'm' are provided (typically for shuttles), or 'i' is provided (typically for buses and rail)."
    )
    trip_splitting_trip['legs']['items']['properties']['b']['description'] = design_leg['board_stop_id']['description']
    trip_splitting_trip['legs']['items']['properties']['a']['description'] = design_leg['alight_stop_id']['description']
    trip_splitting_trip['legs']['items']['properties']['m']['description'] = design_leg['mode']['description']
    trip_splitting_trip['legs']['items']['properties']['i']['description'] = "Leg index, lookup in /design/legs"

    schema['properties']['trips']['description'] = (
        "Detailed trip information for all passenger flows on which the design is based. "
        "Information is presented both aggregate and per leg."
    )
    trips = schema['properties']['trips']['items']['properties']
    trips['passengers']['description'] = (
        "Number of passengers that take this trip. "
        "Note that all trip information is for a SINGLE passenger. "
        "This number is the number of real passengers, so ignoring /parameters/passengerFactor"
    )
    trips['origin_stop_id']['description'] = "Origin stop id (string)"
    trips['destination_stop_id']['description'] = "Destination stop id (string)"
    trips['travel_time']['description'] = "Time between origin and destination"
    trips['waiting_time']['description'] = "Expected waiting time (including fixed waiting time) between origin and destination"
    trips['total_time']['description'] = design_leg['total_time']['description']
    trips['distance']['description'] = "Distance between origin and destination, rounded up to the next meter"
    trips['objective']['description'] = "Contribution of passenger flow to the optimization objective value"
    trips['legs']['description'] = "List of all legs traveled by the passenger, from origin to destination"
    trips_leg = trips['legs']['items']['properties']
    trips_leg['board_stop_id']['description'] = design_leg['board_stop_id']['description']
    trips_leg['alight_stop_id']['description'] = design_leg['alight_stop_id']['description']
    trips_leg['mode']['description'] = design_leg['mode']['description']
    trips_leg['frequency_per_hour']['description'] = design_leg['frequency_per_hour']['description']
    trips_leg['frequency_per_hour_index']['description'] = design_leg['frequency_per_hour_index']['description']
    trips_leg['routes']['description'] = design_leg['routes']['description']
    trips_leg['travel_time']['description'] = design_leg['travel_time']['description']
    trips_leg['waiting_time']['description'] = design_leg['waiting_time']['description']
    trips_leg['total_time']['description'] = design_leg['total_time']['description']
    trips_leg['distance']['description'] = design_leg['distance']['description']
    trips_leg['objective']['description'] = "Contribution of a single passenger on this leg to the optimization objective value, may be set to zero for fixed legs"
    trips_leg['cost']['description'] = "Unweighted contribution of a single passenger on this leg to the cost component of the objective value (using non-rounded distances and durations)"
    trips_leg['convenience']['description'] = "Unweighted contribution of a single passenger on this leg to the convenience component of the objective value (using non-rounded distances and durations)"
    trips_leg['design_leg_index']['description'] = "Leg index, if leg exists (lookup in /design/legs), null otherwise"

    schema['properties']['scores']['description'] = "Overview of aggregate scores resulting from the optimization. These scores may use different units than those defined in the parameters. Convenience is typically based on minutes."
    scores = schema['properties']['scores']['properties']
    scores['total_cost_per_passenger']['description'] = "total_cost/scaled_number_of_passengers"
    scores['passenger_cost_per_passenger']['description'] = "passenger_cost/scaled_number_of_passengers"
    scores['passenger_convenience_per_passenger']['description'] = "passenger_convenience/scaled_number_of_passengers"
    scores['original_number_of_passengers']['description'] = "passenger_cost/scaled_number_of_passengers"
    scores['scaled_number_of_passengers']['description'] = "passenger_cost/scaled_number_of_passengers"