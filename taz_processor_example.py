import pandas as pd
import taz_processor

taz_shape_file = "taz-2010.shp"

# Example location data
locations = pd.DataFrame([[908976,33.753899,-84.391560],
                          [907842,33.758082,-84.387595],
                          [906647,33.640553,-84.446198]],
                         columns=['stop_id','lat','lon'])
locations.set_index('stop_id', inplace=True)

# Initialize TazProcessor
processor = taz_processor.TazProcessor(taz_shape_file)

# Find TAZ object id (oid) for all locations
locations['oid'] = locations.apply(lambda row:
    processor.location_to_taz(row['lat'], row['lon']), axis=1)

# Access TAZ information through taz_processor.taz[oid]
# See taz_processor.py for more information

print(locations)