# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
import os
import pandas as pd
import numpy as np
import copy
import instance

# Determine folder paths
# In the parameters, absolute paths are required
currentFolder = os.path.abspath('').replace('\\', '/')
dataFolder = currentFolder + "/data/"     #folder with input data
outputFolder = currentFolder + "/data/"   #folder that contains output.json
outputFileName = "shuttlemile_1_0_transf_3_alpha_0_05_multipl_5_0_busfreqs_2_3_4_output.json"
roadTravelTimesFileName = "duration_matrix_jan1_oct31.csv"
roadTravelDistancesFileName = "distance_matrix_jan1_oct31.csv"
odPairsFileName = "sample_n_33000_odx.csv"

# Set parameters
param = dict()
param['outputFile'] = outputFolder + outputFileName
param['roadTravelTimesFile'] = dataFolder + roadTravelTimesFileName
param['travelTimeFactorShuttle'] = 1
param['travelTimeFactorBus'] = 1
param['roadTravelDistancesFile'] = dataFolder + roadTravelDistancesFileName

# Create a network design instance and load the existing output file
inst = instance.Instance(runFile=None, parameters=param, logFile=None)
inst.load_outputFile()

# Replace the trips that were used to create the design by other od pairs
# Paths are calculated with the trip splitting information in the instance
inst.replace_trips_by_sample(dataFolder + odPairsFileName)
output = inst.output

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

print(shuttle_legs)

