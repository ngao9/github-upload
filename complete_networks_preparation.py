import os
import pandas as pd
import numpy as np
import copy
import instance
import folium
import json

# Determine folder paths
# In the parameters, absolute paths are required
currentFolder = os.path.abspath('').replace('\\', '/')
dataFolder = currentFolder + "/data/"     #folder with input data
outputFolder = "D:/Dropbox (GaTech)/LeapHi/Pandemic routing/ODMTS designs/"
outputFileName = "shuttlemile_1_0_transf_3_alpha_0_05_multipl_5_0_busfreqs_2_3_4_output.json"
roadTravelTimesFileName = "duration_matrix_jan1_oct31.csv"
roadTravelDistancesFileName = "distance_matrix_jan1_oct31.csv"
stopFileName = "clustered_stops_1500feet_jan1_oct31.csv"
odPairsFileName = "sample_n_33000_odx.csv"

# Set parameters
param = dict()
param['outputFile'] = outputFolder + outputFileName
param['roadTravelTimesFile'] = dataFolder + roadTravelTimesFileName
param['travelTimeFactorShuttle'] = 1
param['roadTravelDistancesFile'] = dataFolder + roadTravelDistancesFileName
param['stopFile'] = dataFolder + stopFileName

# Create a network design instance and load the existing output file
inst = instance.Instance(runFile=None, parameters=param, logFile=None)
inst.load_outputFile()

# Replace the trips that were used to create the design by other od pairs
# Paths are calculated with the trip splitting information in the instance
inst.replace_trips_by_sample(dataFolder + odPairsFileName)
output = inst.output

# Include location information from stop file
locations = pd.read_csv(param['stopFile'],
                        dtype={'stop_id': 'int',
                                'stop_lat': 'float',
                                'stop_lon': 'float'})
locations['stop_id'] = locations['stop_id'].astype(str)
locations.set_index('stop_id', inplace=True)

# Output['trips'] now contains trip information
# Next, isolate shuttle legs

# Add start_times and end_times to trip legs
for trip in output['trips']:
    times = [time for time in pd.to_datetime(trip['departure_times'])]

    for leg in trip['legs']:
        leg['start_times'] = times.copy()
        times = [i + pd.Timedelta(seconds=leg['total_time']) for i in times]
        leg['end_times'] = times.copy()

# Prepare shuttle_Legs DataFrame
shuttle_legs = pd.DataFrame([leg for trip in output['trips'] for leg in trip['legs'] if leg['mode'] == 0])
shuttle_legs = shuttle_legs.explode(column='start_times').reset_index(drop=True)
shuttle_legs['end_times'] = shuttle_legs['start_times'] + shuttle_legs['total_time'].apply(lambda x: pd.Timedelta(seconds=x)) #explode only explodes start_times, so end_times needs to be recalculated
shuttle_legs = shuttle_legs[['board_stop_id','alight_stop_id','start_times','end_times']]
shuttle_legs.rename(columns={'start_times': 'start_time', 'end_times': 'end_time'}, inplace=True)

# Create complete networks data
cn_data = dict()

mode_bus = output['parameters']['modes'].index('bus')
mode_rail = output['parameters']['modes'].index('rail')

def total_freq(hub, mode):
    x = sum([route['frequency_per_hour'] * int(hub in route['stop_ids'])
             for key in output['routes']
             for route in [output['routes'][key]]
             if route['mode'] == mode])
    if abs(x - np.round(x)) < 0.00001:
        x = int(np.round(x))
    return x

def arriving_shuttles(hub):
    return int((shuttle_legs['alight_stop_id'] == hub).sum())

def departing_shuttles(hub):
    return int((shuttle_legs['board_stop_id'] == hub).sum())

cn_data['hubs'] = {hub: { 'lat': locations.loc[hub,'stop_lat'],
                          'lon': locations.loc[hub,'stop_lon'],
                          'total_bus_route_freq': total_freq(hub, mode_bus),
                          'total_rail_route_freq': total_freq(hub, mode_rail),
                          'arriving_shuttles': arriving_shuttles(hub),
                          'departing_shuttles': departing_shuttles(hub) }
                   for hub in output['design']['hubs']}

cn_data['routes'] = output['routes'].copy()

with open(dataFolder + 'cn_data.json', 'w') as file:
    json.dump(cn_data, file, indent=4)

# Create plots
inst.createPlot(point_weight=0)

htmlFile = dataFolder + "design.html"
pngFile = dataFolder + "design.png"
inst.exportPlot(htmlFile, pngFile)

hub_properties = ['total_bus_route_freq', 'total_rail_route_freq', 'arriving_shuttles', 'departing_shuttles']
factors = [50, 50, 1, 1]

for hub_property, factor in zip(hub_properties, factors):

    m = inst.createPlot(point_weight=0)

    for key in cn_data['hubs']:
        hub = cn_data['hubs'][key]
        val = hub[hub_property]
        folium.Circle(location=(hub['lat'], hub['lon']),
                    radius=factor*val,
                    color='black',
                    weight=0.2,
                    fill=True,
                    fill_color='white',
                    fill_opacity=1,
                    popup=hub_property + "=" + str(val)
                    ).add_to(m)

    htmlFile = dataFolder + "design_" + hub_property + ".html"
    pngFile = dataFolder + "design_" + hub_property + ".png"
    inst.exportPlot(htmlFile, pngFile)
