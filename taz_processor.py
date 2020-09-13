import os.path
import shapefile as shp #pyshp package
import shapely.geometry
import collections

class TazProcessor:
    
    # taz dictionary
    # key: OBJECTID value from the TAZ record
    # taz[oid].record: record information
    # taz[oid].shape: shape information
    # taz[oid].geometry: geometry of the shape
    # taz[oid].bounds: bounding box around the taz
    taz = None
    
    def __init__(self, taz_shape_file='taz-2010.shp'):
        
        # Perform checks
        assert taz_shape_file[-4:] == '.shp', \
            "Exception: taz_shape_file does not have extension 'shp'."
        assert os.path.exists(taz_shape_file), \
            "Exception: " + taz_shape_file + " not found."
        assert os.path.exists(taz_shape_file[:-3] + 'dbf'), \
            "Exception: " + taz_shape_file[:-3] + "dbf not found."

        # Read TAZ information
        taz_read = shp.Reader(taz_shape_file)
        assert len(taz_read.records()) == len(taz_read.shapes()), \
            "Exception: shp and dbf files inconsistent."

        # Construct taz dictionary
        TAZ = collections.namedtuple('TAZ', 'record shape geometry bounds')
        self.taz = {record['OBJECTID']: TAZ(record=record,
                                    shape=shape,
                                    geometry=geometry,
                                    bounds=geometry.bounds
                                    )
            for (record, shape) in zip(taz_read.records(), taz_read.shapes())
            for geometry in [shapely.geometry.asShape(shape)]}
        assert len(self.taz) == len(taz_read), "Exception: taz OBJECTIDs not unique."
        del taz_read #remove to avoid confusion between list index and object id

    # Return -1 if no matching taz
    def location_to_taz(self, lat, lon):
        for oid in self.taz:
            taz_oid = self.taz[oid]
            if lon >= taz_oid.bounds[0] and \
               lon <= taz_oid.bounds[2] and \
               lat >= taz_oid.bounds[1] and \
               lat <= taz_oid.bounds[3] and \
               shapely.geometry.Point(lon, lat).within(taz_oid.geometry):
                return oid
        return -1