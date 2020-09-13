import pandas as pd
import numpy as np
import folium
import time
import shutil
from collections import namedtuple
from selenium import webdriver

class PlotDesign:

    __m = None                      # map object
    __bus_colors = ("#FCB913", "#D2662B", "#D1222D", "#000000") # colors corresponding to bus_frequency_index
    __bus_weights = (4, 6, 8, 10)   # bus arrow widths
    __rail_weight = 6               # rail line widths
    __point_weight = 1              # hub size
    __stop_file = None              # file path to stops information file
    __output = None                 # dictionary based on the json output of the optimization

    # Initialize map
    def __make_map(self):
        self.__m = folium.Map(
            location=[33.789705, -84.387789],
            tiles='cartodbpositron',
            zoom_start=11
        )

    def __init__(self, stop_file, output):
        self.__make_map()
        self.__stop_file = stop_file
        self.__output = output

    # Create plot
    # choropleth: choropleth information for coloring map tiles, no coloring if None
    # bus_weights: overwrite __bus_weights, use defaults if None
    # rail_weight: overwrite __rail_weight, use default if None
    # point_weight: overwrite __point_weight, use default if None
    def createPlot(self, choropleth=None, bus_weights=None, rail_weight=None, point_weight=None):

        if not bus_weights == None:
            self.__bus_weights = bus_weights
        if not rail_weight == None:
            self.__rail_weight = rail_weight
        if not point_weight == None:
            self.__point_weight = point_weight

        if not choropleth == None:
            choropleth.add_to(self.__m)

        # Create data DataFrame based on output
        data = pd.DataFrame(columns=['from_id', 'to_id', 'mode', 'bus_frequency_index'])

        rail = self.__output['parameters']['modes'].index('rail')
        for route_key in self.__output['routes']:
            route = self.__output['routes'][route_key]
            if route['mode'] == rail:
                stop_ids = route['stop_ids']
                for i in range(len(stop_ids)):
                    data.loc[len(data)] = [stop_ids[i], stop_ids[(i + 1) % len(stop_ids)], rail, -1]

        bus = self.__output['parameters']['modes'].index('bus')
        for leg in self.__output['design']['legs']:
            if leg['mode'] == bus:
                data.loc[len(data)] = [leg['board_stop_id'], leg['alight_stop_id'], bus, leg['frequency_per_hour_index']]

        shuttle = self.__output['parameters']['modes'].index('shuttle')
        if choropleth == None:
            shuttle_data = [(leg['board_stop_id'], leg['alight_stop_id'], shuttle, -1)
                            for trip in self.__output['trips'] for leg in trip['legs'] if leg['mode'] == shuttle]
            data = pd.concat([data, pd.DataFrame(shuttle_data, columns=['from_id','to_id','mode','bus_frequency_index'])])

        # Include location information from stop file
        locations = pd.read_csv(self.__stop_file,
                                dtype={'stop_id': 'int',
                                       'stop_lat': 'float',
                                       'stop_lon': 'float'})

        locations['stop_id'] = locations['stop_id'].astype(str)
        locations.set_index('stop_id', inplace=True)

        data = data.merge(locations, how='left', left_on='from_id', right_index=True)
        data.rename(columns={'stop_lat': 'from_lat', 'stop_lon': 'from_lon'}, inplace=True)
        data = data.merge(locations, how='left', left_on='to_id', right_index=True)
        data.rename(columns={'stop_lat': 'to_lat', 'stop_lon': 'to_lon'}, inplace=True)

        data['points'] = data.apply(lambda x: [(x['from_lat'], x['from_lon']), (x['to_lat'], x['to_lon'])], axis=1)

        # Create separate DataFrames from the three modes
        shuttleData = data.loc[data['mode'] == 0, ['points']]
        busData = data.loc[data['mode'] == 1, ['points', 'bus_frequency_index']]
        railData = data.loc[data['mode'] == 2, ['points']]

        # Select all hubs for which opened = True in the output
        usedHubs = [(locations.loc[key,'stop_lat'], locations.loc[key,'stop_lon'])
                    for key in self.__output['design']['hubs'] if self.__output['design']['hubs'][key]['opened']]

        # Create the plot
        for index, row in shuttleData.iterrows():
            if choropleth == None:
                self.__plot_shuttle(row['points'])

        for index, row in railData.iterrows():
            self.__plot_rail(row['points'])

        for index, row in busData.iterrows():
            self.__plot_bus(row['points'], row['bus_frequency_index'])

        self.__plot_hubs(usedHubs)

        for index, row in busData.iterrows():
            self.__plot_bus_arrow_marker(row['points'], row['bus_frequency_index'])

        return self.__m

    # Export created plot to html and png
    # htmlFile: file path for html export
    # pngFile: file path for png export
    #   if pngFile=None, do not export png
    #   png export requires the selenium python package, Firefox browser, and the geckodriver to be in the path
    #   if geckodriver cannot be found, png export is automatically skipped
    def exportPlot(self, htmlFile, pngFile=None):

        self.__m.save(htmlFile)

        # Beautiful hack by Kevin to remove the gradient on the arrow heads (which look very bad)
        with open(htmlFile, "a") as html:
            html.write(
                "<script> var elements = document.getElementsByTagName(\"linearGradient\"); for(i = 0; i < elements.length; i++){ elements[i].setAttribute(\"x1\", \"100%\"); elements[i].setAttribute(\"y1\", \"100%\"); } </script>")

        # Use selenium to open the html and take a screenshot
        if pngFile != None and shutil.which("geckodriver") != None: #geckodriver for Firefoxs needs to be accessible in path
            browser = webdriver.Firefox()
            browser.set_window_position(0, 0)
            browser.set_window_size(1920, 1200)
            browser.get('file://' + htmlFile)
            time.sleep(5) # Give the map tiles some time to load
            browser.save_screenshot(pngFile)
            browser.quit()

        return self.__m

    # Adapted from https://medium.com/@bobhaffner/folium-lines-with-arrows-25a0fe88e4e
    def __get_bearing(self, p1, p2):

        long_diff = np.radians(p2.lon - p1.lon)

        lat1 = np.radians(p1.lat)
        lat2 = np.radians(p2.lat)

        x = np.sin(long_diff) * np.cos(lat2)
        y = (np.cos(lat1) * np.sin(lat2)
             - (np.sin(lat1) * np.cos(lat2)
                * np.cos(long_diff)))

        bearing = np.degrees(np.arctan2(x, y))

        # adjusting for compass bearing
        if bearing < 0:
            return bearing + 360
        return bearing

    # Adapted from https://medium.com/@bobhaffner/folium-lines-with-arrows-25a0fe88e4e
    def __get_arrow(self, locations, color='blue', size=6, fraction=0.6):

        Point = namedtuple('Point', field_names=['lat', 'lon'])

        # creating point from our Point named tuple
        p1 = Point(locations[0][0], locations[0][1])
        p2 = Point(locations[1][0], locations[1][1])

        # getting the rotation needed for our marker.
        # Subtracting 90 to account for the marker's orientation
        # of due East(get_bearing returns North)
        rotation = self.__get_bearing(p1, p2) - 90

        arrow_lats = np.array([p1.lat + fraction * (p2.lat - p1.lat)])
        arrow_lons = np.array([p1.lon + fraction * (p2.lon - p1.lon)])

        arrows = []

        # creating each "arrow" and appending them to our arrows list
        for points in zip(arrow_lats, arrow_lons):
            arrows.append(folium.RegularPolygonMarker(location=points,
                                                      line_cap='square',
                                                      line_join='miter',
                                                      color='black', opacity=1,
                                                      weight=0.2,
                                                      fill_color=color, fill_opacity=1,
                                                      number_of_sides=3,
                                                      radius=size, rotation=rotation))
        return arrows[0]

    def __plot_shuttle(self, points, color="#66C2A5", opacity=0.1):
        folium.PolyLine(points, color=color, opacity=opacity, weight=2).add_to(self.__m)

    def __plot_bus(self, points, frequency_index):
        folium.PolyLine(points, color=self.__bus_colors[frequency_index], weight=self.__bus_weights[frequency_index]).add_to(self.__m)

    def __plot_bus_arrow_marker(self, points, frequency_index):
        self.__get_arrow(points, color=self.__bus_colors[frequency_index], size=2.5 * self.__bus_weights[0], fraction=0.65).add_to(self.__m)

    def __plot_rail(self, points):
        folium.PolyLine(points, color="#0B75B3", weight=self.__rail_weight).add_to(self.__m)

    def __plot_hubs(self, points):

        for point in points:
            folium.Circle(location=point,
                          radius=300,
                          color='black',
                          weight=self.__point_weight,
                          fill=True,
                          fill_color='black',
                          fill_opacity=1,
                          ).add_to(self.__m)