#%%
import os
import pandas as pd
import numpy as np
import copy
import itertools
import importlib

import instance
import taz_processor

# Determine folder paths
# In the parameters, absolute paths are required
currentFolder = os.path.abspath('').replace('\\', '/')
runFile = "C:/Users/ningg/LeapHi-Benders/bin/leaphi.exe" 
dataFolder = currentFolder + "/data/"       #folder must exist
outputFolder = currentFolder + "/output/"   #folder must exist

# Lambda functions to calculate values from trips for choropleth maps
choro_lambdas = dict()
choro_lambdas['transfers'] = lambda trip: len(trip['legs'])-1
choro_lambdas['cost'] = lambda trip: sum([leg['cost'] for leg in trip['legs']])
choro_lambdas['convenience'] = lambda trip: sum([leg['convenience'] for leg in trip['legs']])
choro_lambdas['waiting_time'] = lambda trip: trip['waiting_time']
choro_lambdas['travel_time'] = lambda trip: trip['travel_time']

# Prepare TAZ processor
processor = taz_processor.TazProcessor()

# Fixed parameters in all experiments
fixed_param = dict()
fixed_param['stopFile'] = dataFolder + "clustered_stops_1500feet.csv"
fixed_param['railLineFile'] = dataFolder + "rail_lines.csv"
fixed_param['tripsFile'] = dataFolder + "odx_clustered_april3_jan1_oct31_with_departues.csv"
fixed_param['timeHorizonInHours'] = 4
fixed_param['roadTravelTimesFile'] = dataFolder + "dur_mins.csv"
fixed_param['roadTravelDistancesFile'] = dataFolder + "dist_km.csv"
fixed_param['travelTimeFactorShuttle'] = 1.2
fixed_param['travelTimeFactorBus'] = 1.2

#original with only rail  (38)
railLats = [33.753899, 33.758082, 33.640553, 33.750161, 33.789669, 33.771472, 
            33.781247, 33.774699, 33.92113, 33.823501, 33.736564, 33.650282, 
            33.847944, 33.677424, 33.766111, 33.754517, 33.845388, 33.910757, 
            33.772647, 33.756179, 33.931671, 33.748696, 33.887929, 33.94495, 
            33.7173, 33.756478, 33.749951, 33.769856, 33.902787, 33.70088, 
            33.775554, 33.772159, 33.860329, 33.753247, 33.757451, 33.756613, 
            33.765241, 33.761803]
railLons = [-84.39156, -84.387595, -84.446198, -84.385915, -84.387414, -84.387258, 
            -84.386342, -84.295417, -84.344268, -84.369313, -84.413653, -84.448612, 
            -84.367716, -84.440542, -84.387504, -84.469824, -84.358235, -84.35189, 
            -84.251607, -84.397215, -84.351069, -84.395741, -84.305556, -84.357275, 
            -84.42503, -84.41723, -84.375675, -84.228906, -84.28061, -84.428768, 
            -84.281487, -84.428873, -84.339245, -84.445568, -84.352762, -84.403902, 
            -84.312937, -84.340825]

#region 1 [38]
potentialHubLats_1 = []
potentialHubLons_1 = []

#region 2 [38]
potentialHubLats_2 = []
potentialHubLons_2 = [] 

#region 3 (6) [44/21]
potentialHubLats_3 = [34.040397, 34.021794, 34.069693, 34.073296, 34.040324, 34.064854]
potentialHubLons_3 = [-84.319276, -84.325413, -84.288296, -84.257132, -84.319267, -84.253526]

#region 4 (1) [39/rest]
potentialHubLats_4 = [34.095857]
potentialHubLons_4 = [-84.240906]

# #region 5 [38]
potentialHubLats_5 = []
potentialHubLons_5 = []

# #region 6 (5) [43/rest]
potentialHubLats_6 = [33.911130, 33.804432, 33.876609, 33.842910, 33.819182]
potentialHubLons_6 = [-84.424681, -84.476740, -84.464372, -84.426632, -84.453326]

# #region 7 (10) [48/rest]
potentialHubLats_7 = [33.801821, 33.917168, 33.840310, 33.828414, 33.985826, 33.927261, 33.882017, 33.944950, 33.931671, 33.834237]
potentialHubLons_7 = [-84.346760, -84.349101, -84.264194, -84.333299, -84.349743, -84.269468, -84.251113, -84.357275, -84.351069, -84.271889]

# #region 8 (5) [43/rest]
potentialHubLats_8 = [33.804499, 33.849433, 33.850511, 33.801023, 33.819005]
potentialHubLons_8 = [-84.191053, -84.206255, -84.216162, -84.208339, -84.193778]

# #region 9 (4) [42/rest]
potentialHubLats_9 = [33.725034, 33.735991, 33.771776, 33.725821]
potentialHubLons_9 = [-84.583102, -84.565001, -84.552421, -84.582613]

# #region 10 (20) [58/rest]
potentialHubLats_10 = [33.699048,33.717300, 33.709742, 33.603977, 33.677424, 33.700880, 33.736564, 33.648978, 33.657361, 33.654336, 33.788669, 33.684248, 33.745768, 33.659014, 33.763512, 33.767062, 33.729619, 33.614435, 33.786399, 33.706005]
potentialHubLons_10 = [-84.425965, -84.425030, -84.506663, -84.475831, -84.440542, -84.428768, -84.413653, -84.448288, -84.506595, -84.477822, -84.407668, -84.500844, -84.406163, -84.426166, -84.437984, -84.401653, -84.510554, -84.449344, -84.494598, -84.506485]

# #region 11 (15) [53/rest]
potentialHubLats_11 = [33.716633, 33.765124, 33.716731, 33.688895, 33.730776, 33.717266, 33.727363, 33.773840, 33.793652, 33.626633, 33.631202, 33.704926, 33.715668, 33.644332, 33.704363]
potentialHubLons_11 = [-84.254904, -84.312886, -84.325732, -84.337031, -84.393172, -84.309760, -84.362304, -84.352985, -84.280721, -84.367963, -84.289667, -84.277541, -84.392104, -84.325288, -84.308340]

# #region 12 (11) [49/rest]
potentialHubLats_12 =  [33.736814, 33.745833, 33.706990, 33.783276, 33.747820, 33.776248, 33.700374, 33.714038, 33.759863, 33.797554, 33.790184]
potentialHubLons_12 = [-84.190468, -84.136800, -84.113397, -84.162154, -84.166202, -84.194969, -84.168424, -84.213085, -84.174916, -84.223930, -84.176568]

# #region 12 (11) [49/rest]
potentialHubLats_12 =  [33.736814, 33.745833, 33.706990, 33.783276, 33.747820, 33.776248, 33.700374, 33.714038, 33.759863, 33.797554, 33.790184]
potentialHubLons_12 = [-84.190468, -84.136800, -84.113397, -84.162154, -84.166202, -84.194969, -84.168424, -84.213085, -84.174916, -84.223930, -84.176568]

# #region 13 (2) [40/rest]
potentialHubLats_13 = [33.567570, 33.523702]
potentialHubLons_13 = [-84.579634, -84.664692]

# #region 14 (6) [44/rest]
potentialHubLats_14 = [33.587253, 33.562847, 33.570023, 33.587319, 33.585073, 33.534351]
potentialHubLons_14 = [-84.549811, -84.517395, -84.538593, -84.548837, -84.469175, -84.418769]

# #region 15 (3) [41/rest]
potentialHubLats_15 = [33.573365, 33.520005, 33.571945]
potentialHubLons_15 =  [-84.389090, -84.361797, -84.354724]

#region 16 [38]
potentialHubLats_16 = []
potentialHubLons_16 =  []

#region 3 (6) [44/21] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_2 + potentialHubLats_3 + potentialHubLats_4 + potentialHubLats_6 + potentialHubLats_7 + potentialHubLats_8
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_2 + potentialHubLons_3 + potentialHubLons_4 + potentialHubLons_6 + potentialHubLons_7 + potentialHubLons_8
# print(len(fixed_param['desiredHubLats']))
# rail_region_3_weight = [1 for i in range(len(railLats) + len(potentialHubLats_3))]
# neighbor_region_3_weight = [0 for i in range(len(potentialHubLats_2) + len(potentialHubLats_4) + len(potentialHubLats_6) + len(potentialHubLats_7) + len(potentialHubLats_8))]
# fixed_param['hubWeights'] = rail_region_3_weight + neighbor_region_3_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_3
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_3
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_3))]

#region 4 (1) [39/21] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_4 + potentialHubLats_3 + potentialHubLats_7 + potentialHubLats_8
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_4 + potentialHubLons_3 + potentialHubLons_7 + potentialHubLons_8

# rail_region_4_weight = [1 for i in range(len(railLats) + len(potentialHubLats_4))]
# neighbor_region_4_weight = [0 for i in range(len(potentialHubLats_3) + len(potentialHubLats_7) + len(potentialHubLats_8))]
# fixed_param['hubWeights'] = rail_region_4_weight + neighbor_region_4_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_4
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_4
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_4))]

# #region 5

# #region 6 (5) [43/55] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_6 + potentialHubLats_3 + potentialHubLats_7 + potentialHubLats_9 + potentialHubLats_10 + potentialHubLats_11
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_6 + potentialHubLons_3 + potentialHubLons_7 + potentialHubLons_9 + potentialHubLons_10 + potentialHubLons_11

# rail_region_6_weight = [1 for i in range(len(railLats) + len(potentialHubLats_6))]
# neighbor_region_6_weight = [0 for i in range(len(potentialHubLats_3) + len(potentialHubLats_7) + len(potentialHubLats_9) + len(potentialHubLats_10) + len(potentialHubLats_11))]
# fixed_param['hubWeights'] = rail_region_6_weight + neighbor_region_6_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_6
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_6
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_6))]

# #region 7 (10) [48/63] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_7 + potentialHubLats_2 + potentialHubLats_3 + potentialHubLats_4 + potentialHubLats_6 + potentialHubLats_8 + potentialHubLats_10 + potentialHubLats_11 + potentialHubLats_12
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_7 + potentialHubLons_2 + potentialHubLons_3 + potentialHubLons_4 + potentialHubLons_6 + potentialHubLons_8 + potentialHubLons_10 + potentialHubLons_11 + potentialHubLons_12

# rail_region_7_weight = [1 for i in range(len(railLats) + len(potentialHubLats_7))]
# neighbor_region_7_weight = [0 for i in range(len(potentialHubLats_2) + len(potentialHubLats_3) + len(potentialHubLats_4) + len(potentialHubLats_6) + len(potentialHubLats_8) + len(potentialHubLats_10) + len(potentialHubLats_11) + len(potentialHubLats_12))]
# fixed_param['hubWeights'] = rail_region_7_weight + neighbor_region_7_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_7
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_7
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_7))]

# #region 8 (5) [43/43] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_8 + potentialHubLats_3 + potentialHubLats_4 + potentialHubLats_7 + potentialHubLats_11 + potentialHubLats_12
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_8 + potentialHubLons_3 + potentialHubLons_4 + potentialHubLons_7 + potentialHubLons_11 + potentialHubLons_12

# rail_region_8_weight = [1 for i in range(len(railLats) + len(potentialHubLats_8))]
# neighbor_region_8_weight = [0 for i in range(len(potentialHubLats_3) + len(potentialHubLats_4) + len(potentialHubLats_7) + len(potentialHubLats_11) + len(potentialHubLats_12))]
# fixed_param['hubWeights'] = rail_region_8_weight + neighbor_region_8_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_8
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_8
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_8))]

# #region 9 (4) [42/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_9 + potentialHubLats_5 + potentialHubLats_6 + potentialHubLats_10 + potentialHubLats_13 + potentialHubLats_14
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_9 + potentialHubLons_5 + potentialHubLons_6 + potentialHubLons_10 + potentialHubLons_13 + potentialHubLons_14

# rail_region_9_weight = [1 for i in range(len(railLats) + len(potentialHubLats_9))]
# neighbor_region_9_weight = [0 for i in range(len(potentialHubLats_5) + len(potentialHubLats_6) + len(potentialHubLats_10) + len(potentialHubLats_13) + len(potentialHubLats_14))]
# fixed_param['hubWeights'] = rail_region_9_weight + neighbor_region_9_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_9
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_9
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_9))]

# #region 10 (20) [58/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_10 + potentialHubLats_5 + potentialHubLats_6 + potentialHubLats_7 + potentialHubLats_9 + potentialHubLats_11 + potentialHubLats_13 + potentialHubLats_14 + potentialHubLats_15
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_10 + potentialHubLons_5 + potentialHubLons_6 + potentialHubLons_7 + potentialHubLons_9 + potentialHubLons_11 + potentialHubLons_13 + potentialHubLons_14 + potentialHubLons_15 

# rail_region_10_weight = [1 for i in range(len(railLats) + len(potentialHubLats_10))]
# neighbor_region_10_weight = [0 for i in range(len(potentialHubLats_5) + len(potentialHubLats_6) + len(potentialHubLats_7) + len(potentialHubLats_9) + len(potentialHubLats_11) + len(potentialHubLats_13) + len(potentialHubLats_14) + len(potentialHubLats_15))]
# fixed_param['hubWeights'] = rail_region_10_weight + neighbor_region_10_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_10
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_10
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_10))]

# #region 11 (15) [53/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_11 + potentialHubLats_6 + potentialHubLats_7 + potentialHubLats_8 + potentialHubLats_10 + potentialHubLats_12 + potentialHubLats_14 + potentialHubLats_15 + potentialHubLats_16
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_11 + potentialHubLons_6 + potentialHubLons_7 + potentialHubLons_8 + potentialHubLons_10 + potentialHubLons_12 + potentialHubLons_14 + potentialHubLons_15 + potentialHubLons_16 

# rail_region_11_weight = [1 for i in range(len(railLats) + len(potentialHubLats_11))]
# neighbor_region_11_weight = [0 for i in range(len(potentialHubLats_6) + len(potentialHubLats_7) + len(potentialHubLats_8) + len(potentialHubLats_10) + len(potentialHubLats_12) + len(potentialHubLats_14) + len(potentialHubLats_15) + len(potentialHubLats_16))]
# fixed_param['hubWeights'] = rail_region_11_weight + neighbor_region_11_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_11
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_11
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_11))]

# #region 12 (11) [49/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_12 + potentialHubLats_7 + potentialHubLats_8 + potentialHubLats_11 + potentialHubLats_15 + potentialHubLats_16
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_12 + potentialHubLons_7 + potentialHubLons_8 + potentialHubLons_11 + potentialHubLons_15 + potentialHubLons_16

# rail_region_12_weight = [1 for i in range(len(railLats) + len(potentialHubLats_12))]
# neighbor_region_12_weight = [0 for i in range(len(potentialHubLats_7) + len(potentialHubLats_8) + len(potentialHubLats_11) + len(potentialHubLats_15) + len(potentialHubLats_16))]
# fixed_param['hubWeights'] = rail_region_12_weight + neighbor_region_12_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_12
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_12
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_12))]

# #region 13 (2) [40/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_13 + potentialHubLats_9 + potentialHubLats_10 + potentialHubLats_14
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_13 + potentialHubLons_9 + potentialHubLons_10 + potentialHubLons_14

# rail_region_13_weight = [1 for i in range(len(railLats) + len(potentialHubLats_13))]
# neighbor_region_13_weight = [0 for i in range(len(potentialHubLats_9) + len(potentialHubLats_10) + len(potentialHubLats_14))]
# fixed_param['hubWeights'] = rail_region_13_weight + neighbor_region_13_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_13
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_13
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_13))]

# #region 14 (6) [44/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_14 + potentialHubLats_9 + potentialHubLats_10 + potentialHubLats_11 + potentialHubLats_13 + potentialHubLats_15
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_14 + potentialHubLons_9 + potentialHubLons_10 + potentialHubLons_11 + potentialHubLons_13 + potentialHubLons_15

# rail_region_14_weight = [1 for i in range(len(railLats) + len(potentialHubLats_14))]
# neighbor_region_14_weight = [0 for i in range(len(potentialHubLats_9) + len(potentialHubLats_10) + len(potentialHubLats_11) + len(potentialHubLats_13) + len(potentialHubLats_15))]
# fixed_param['hubWeights'] = rail_region_14_weight + neighbor_region_14_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_14
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_14
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_14))]

# #region 15 (3) [41/rest] parameters
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_15 + potentialHubLats_10 + potentialHubLats_11 + potentialHubLats_12 + potentialHubLats_14 + potentialHubLats_16
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_15 + potentialHubLons_10 + potentialHubLons_11 + potentialHubLons_12 + potentialHubLons_14 + potentialHubLons_16

# rail_region_15_weight = [1 for i in range(len(railLats) + len(potentialHubLats_15))]
# neighbor_region_15_weight = [0 for i in range(len(potentialHubLats_10) + len(potentialHubLats_11) + len(potentialHubLats_12) + len(potentialHubLats_14) + len(potentialHubLats_16))]
# fixed_param['hubWeights'] = rail_region_15_weight + neighbor_region_15_weight

# fixed_param['desiredHubLats'] = railLats + potentialHubLats_15
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_15
# fixed_param['hubWeights'] = [1 for i in range(len(railLats) + len(potentialHubLats_15))]

#region 16


#selected hub candidates
'''
region 1: none
region 2: none
region 3: [34.07329619983435, -84.2571320398017] [34.040396502203244, -84.31927583700414] [34.021793807692326, -84.3254132692307]
region 4: [34.09585738345859, -84.24090608270691]
region 5: none
region 6: none
region 7: [33.98582619029845, -84.34974328731334]
region 8: none
region 9: none
region 10: [33.729618999999964, -84.51055400000004] [33.603976548616984, -84.47583116193249] [33.614434799999984, -84.44934359999993]
region 11: [33.62663331944446, -84.36796263690462] [33.70492561666675, -84.27754104166675]
region 12: [33.79755400000009, -84.2239300000001] [33.71403767452817, -84.21308545754698] [33.70699007142855, -84.11339726428584]
region 13: none
region 14: [33.585073, -84.469175] [33.534351, -84.418769]
region 15: [33.57336455995267, -84.38909049200626] [33.571945, -84.354724]
region 16: none
'''
#selected hub candidates lat & lon

hub_candidates_lats = [34.07329619983435, 34.040396502203244, 34.021793807692326, 34.09585738345859, 33.98582619029845, 33.729618999999964, 33.603976548616984,
33.614434799999984, 33.62663331944446, 33.70492561666675, 33.79755400000009, 33.71403767452817, 33.70699007142855, 33.585073, 33.534351, 33.57336455995267,
 33.571945]

hub_candidates_lons = [-84.2571320398017, -84.31927583700414, -84.3254132692307, -84.24090608270691, -84.34974328731334, -84.51055400000004, -84.47583116193249,
-84.44934359999993, -84.36796263690462, -84.27754104166675, -84.2239300000001, -84.21308545754698, -84.11339726428584, -84.469175, -84.418769, -84.38909049200626, 
-84.354724]

#selected hub candidates + rail 
fixed_param['desiredHubLats'] = railLats + hub_candidates_lats
fixed_param['desiredHubLons'] = railLons + hub_candidates_lons

#naive hub candidates candidates + rail
# fixed_param['desiredHubLats'] = railLats + potentialHubLats_1 + potentialHubLats_2 + potentialHubLats_3 + potentialHubLats_4 + potentialHubLats_5 + potentialHubLats_6 + potentialHubLats_7 + potentialHubLats_8 + potentialHubLats_9 + potentialHubLats_10 + potentialHubLats_11 + potentialHubLats_12 + potentialHubLats_13 + potentialHubLats_14 + potentialHubLats_15 + potentialHubLats_16
# fixed_param['desiredHubLons'] = railLons + potentialHubLons_1 + potentialHubLons_2 + potentialHubLons_3 + potentialHubLons_4 + potentialHubLons_5 + potentialHubLons_6 + potentialHubLons_7 + potentialHubLons_8 + potentialHubLons_9 + potentialHubLons_10 + potentialHubLons_11 + potentialHubLons_12 + potentialHubLons_13 + potentialHubLons_14 + potentialHubLons_15 + potentialHubLons_16
                         
fixed_param['fixedTransferTime'] = 0
fixed_param['mipGap'] = 0.01
fixed_param['timeLimit'] = 43200 #seconds
fixed_param['busCostPerHour'] = 100
fixed_param['ignoreFrequencyCorrection'] = False
fixed_param['useGurobi'] = False
fixed_param['generateTripSplittings'] = False
fixed_param['enableGNUplot'] = False
fixed_param['outputFile'] = outputFolder + "output.json" #Set for each experiment

#Parameters that vary between experiments
#For every parameter, a list is given, and an experiment is conducted for each combination of values
# vary_param = dict()
# vary_param['shuttleCostPerKm'] = [costPerMile / 1.609344 for costPerMile in (0.75, 1.0, 1.25)]
# vary_param['maximumNumberOfTransfers'] = [2, 3]
# vary_param['alpha'] = [0.05, 0.1]
# vary_param['passengerFactor'] = [1, 2, 8]
# vary_param['busesPerHour'] = [(2, 3, 4)]

vary_param = dict()
vary_param['shuttleCostPerKm'] = [1.0 / 1.609344]
vary_param['maximumNumberOfTransfers'] = [3]
vary_param['alpha'] = [0.05]
vary_param['passengerFactor'] = [5.0]
vary_param['busesPerHour'] = [(2, 3, 4)]

# Turn a dictionary into a prefix for a file name
def dict_str(x):
    result = ""

    for key in x.keys():

        if result != "":
            result += "_"

        if key == "shuttleCostPerKm":
            result += 'shuttlemile_' + str(x[key] * 1.609344).replace('.', '_')
        elif key == "maximumNumberOfTransfers":
            result += 'transf_' + str(x[key])
        elif key == "alpha":
            result += 'alpha_' + str(x[key]).replace('.', '_')
        elif key == "passengerFactor":
            result += 'multipl_' + str(x[key]).replace('.', '_')
        elif key == "busesPerHour":
            result += 'busfreqs_' + '_'.join(map(str, x[key]))
        else:
            value_str = key + "_" + str(x[key])
            value_str = value_str.replace(",", "_")
            keepcharacters = ('_')
            value_str = "".join(c for c in value_str if c.isalnum() or c in keepcharacters).rstrip()
            result += key + value_str

    return result


# Prepare list of experiments
# Every experiment is a dictionary with keys
#   string: suggested file name prefix based on parameters
#   log: file path to write the log file
#   param: complete set of parameters for this experiment
experiments = list()
for value in itertools.product(*vary_param.values()):
    experiment = dict()

    exp_param = dict(zip(vary_param.keys(), value))

    experiment['string'] = dict_str(exp_param)
    experiment['log'] = outputFolder + experiment['string'] + "_log.txt"
    experiment['param'] = copy.deepcopy(fixed_param)
    experiment['param'].update(exp_param)
    experiment['param']['outputFile'] = outputFolder + experiment['string'] + "_output.json"

    experiments.append(experiment)

for experiment in experiments:

	# Run the experiments and save output and log to disk
    inst = instance.Instance(runFile, experiment['param'], experiment['log'])
    inst.run()

    # Verify correctness
    inst.verifyOutput() #See Instance.verifyOutput() in instance.py for possible options

# Create maps and histograms
for experiment in experiments:
    inst = instance.Instance(runFile, experiment['param'], experiment['log'])
    inst.load_outputFile()

    htmlFile = outputFolder + experiment['string'] + "_design.html"
    pngFile = outputFolder + experiment['string'] + "_design.png"
    inst.createPlot()
    inst.exportPlot(htmlFile, pngFile)

    hist_type = 'transfers'
    inst.createHistogram(hist_type)
    pngFile = outputFolder + experiment['string'] + "_hist_" + hist_type + ".png"
    inst.exportHist(pngFile)

# Create maps, histograms, and choropleths for sampled demand
for experiment in experiments:
    inst = instance.Instance(runFile, experiment['param'], experiment['log'])
    inst.load_outputFile()

    #inst.replace_trips_by_sample(dataFolder + "sample_n_33000_odx.csv")
    #inst.verifyOutput(backward_compatibility=True) #See Instance.verifyOutput() in instance.py for possible options
                                                   #Backward compatibility needed for pre June 2020 designs
    outputFile_sample = outputFolder + experiment['string'] + "_sample_output.json"
    inst.write_outputFile(outputFile_sample, ignore_trip_splitting=True)

    htmlFile = outputFolder + experiment['string'] + "_sample_design.html"
    pngFile = outputFolder + experiment['string'] + "_sample_design.png"
    inst.createPlot()
    inst.exportPlot(htmlFile, pngFile)

    hist_type = 'transfers'
    inst.createHistogram(hist_type)
    pngFile = outputFolder + experiment['string'] + "_sample_hist_" + hist_type + ".png"
    inst.exportHist(pngFile)
    
    for key in choro_lambdas:
        htmlFile = outputFolder + experiment['string'] + "_sample_" + key +  "_map.html"
        pngFile = outputFolder + experiment['string'] + "_sample_" + key +  "_map.png"
        choropleth = inst.prepareTripChoropleth(trip_value=choro_lambdas[key],
                                                trip_stop_id=lambda trip: trip['origin_stop_id'],
                                                processor=processor)
        inst.createPlot(choropleth=choropleth,
                        bus_weights=(2, 3, 4, 5),
                        rail_weight=3,
                        point_weight=0.1)
        inst.exportPlot(htmlFile, pngFile)