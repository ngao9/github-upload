#include "kevintools_functions.hpp"
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <algorithm>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double kevintools::euclideanDistance(const std::vector<double>& x, const std::vector<double>& y){

    if(x.size() != y.size()){
        throw std::runtime_error("Exception: x.size() != y.size().");
    }

    double result = 0;

    for(int i = 0; i < x.size(); ++i){
        result += (x.at(i) - y.at(i)) * (x.at(i) - y.at(i));
    }

    result = std::sqrt(result);

    return result;

}

//https://stackoverflow.com/a/1120224/11428765
std::vector<std::string> kevintools::getNextLineAndSplitIntoTokens(std::istream &str) {

    std::vector<std::string>   result;
    std::string                line;
    std::getline(str,line);

    std::stringstream          lineStream(line);
    std::string                cell;

    while(std::getline(lineStream,cell, ','))
    {
        result.push_back(cell);
    }

    //remove line endings
    std::string& lastCell = result.at(result.size()-1);
    lastCell.erase(std::remove(lastCell.begin(), lastCell.end(), '\n'), lastCell.end());
    lastCell.erase(std::remove(lastCell.begin(), lastCell.end(), '\r'), lastCell.end());

    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;

}

// This function converts decimal degrees to radians
double kevintools::deg2rad(double deg) {
    return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double kevintools::rad2deg(double rad) {
    return (rad * 180 / M_PI);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 * //https://stackoverflow.com/a/10205532/11428765
 */
double kevintools::distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r)/2);
    v = sin((lon2r - lon1r)/2);
    return 2.0 * 6371.0 * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}
