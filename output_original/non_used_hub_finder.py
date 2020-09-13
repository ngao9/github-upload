import json 
import pandas as pd 
  
# open output JSON file 
with open('shuttlemile_1_0_transf_3_alpha_0_05_multipl_5_0_busfreqs_2_3_4_output.json') as f:
    data = json.load(f)

# open location file
geo_location = pd.read_csv('clustered_stops_1500feet.csv') 

# create used stop_id list
stop_id_used = []

# find stop_id selected
for key,values in data['routes'].items():
    stop_ids = values['stop_ids']
    for stop_id in stop_ids:
        if stop_id not in stop_id_used:
            stop_id_used.append(float(stop_id))

# find selected hub location 
stop_id_used_df = pd.DataFrame(stop_id_used,columns=['stop_ids'])
hub_selected_df = geo_location[geo_location['stop_id'].isin(stop_id_used_df['stop_ids'])]

# output selected hub location + rail location
#print(hub_selected_df)

'''
[33.753899, 33.758082, 33.640553, 33.750161, 33.789669,
33.771472, 33.781247, 33.774699, 33.921130, 33.823501,
33.736564, 33.650282, 33.847944, 33.677424, 33.766111,
33.754517, 33.845388, 33.910757, 33.772647, 33.756179,
33.931671, 33.748696, 33.887929, 33.944950, 33.717300,
33.756478, 33.749951, 33.769856, 33.902787, 33.700880,
33.775554, 33.772159, 33.860329, 33.753247, 33.757451,
33.756613, 33.765241, 33.761803, 33.705515, 33.706005,
33.704363, 33.571945, 34.040324, 33.585073, 34.064854]

[-84.391560, -84.387595, -84.446198, -84.385915, -84.387414,
-84.387258, -84.386342, -84.295417,-84.344268, -84.369313,
-84.413653, -84.448612, -84.367716, -84.440542, -84.387504,
-84.469824, -84.358235, -84.351890, -84.251607, -84.397215,
-84.351069, -84.395741, -84.305556, -84.357275, -84.425030,
-84.417230, -84.375675, -84.228906, -84.280610, -84.428768,
-84.281487, -84.428873, -84.339245, -84.445568, -84.352762,
-84.403902, -84.312937, -84.340825, -84.135234, -84.506485,
-84.308340, -84.354724, -84.319267, -84.469175, -84.253526]
'''

#find rail_id selected
rail_id = []
all_hubs = data['design']['hubs']
for key, values in all_hubs.items():
    if values['fixed'] == True:
        rail_id.append(float(key))

#find rail hub location
rail_id_df = pd.DataFrame(rail_id,columns=['rail_ids'])
rail_selected_df = geo_location[geo_location['stop_id'].isin(rail_id_df['rail_ids'])]

#output rail location
# print(rail_selected_df)
rail_lat = rail_selected_df['stop_lat'].tolist()
rail_lon = rail_selected_df['stop_lon'].tolist()
print(rail_lat)
print(rail_lon)

'''
[33.753899, 33.758082, 33.640553000000004, 33.750161, 33.789669, 33.771471999999996, 33.781247, 
33.774699, 33.92113, 33.823501, 33.736564, 33.650282000000004, 33.847944, 33.677424, 33.766110999999995, 
33.754517, 33.845388, 33.910757000000004, 33.772647, 33.756178999999996, 33.931671, 33.748696, 33.887929, 
33.94495, 33.7173, 33.756478, 33.749951, 33.769856, 33.902787, 33.70088, 33.775554, 33.772159, 33.860329, 
33.753246999999995, 33.757451, 33.756613, 33.765240999999996, 33.761803]

[-84.39156, -84.38759499999998, -84.446198, -84.385915, -84.387414, -84.387258, -84.386342, -84.295417, 
-84.344268, -84.369313, -84.413653, -84.448612, -84.367716, -84.44054200000002, -84.387504, -84.469824, 
-84.358235, -84.35189, -84.25160699999998, -84.397215, -84.351069, -84.395741, -84.30555600000002, -84.357275, 
-84.42503, -84.41723, -84.375675, -84.22890600000002, -84.28061, -84.428768, -84.281487, -84.428873, 
-84.33924499999998, -84.44556800000002, -84.352762, -84.403902, -84.312937, -84.340825]