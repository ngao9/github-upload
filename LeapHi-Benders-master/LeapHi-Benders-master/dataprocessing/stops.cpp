#include "stops.hpp"

#include "kevintools_functions.hpp"

#include <string>
#include <tuple>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace dataprocessing;

VStop Stops::readStopInfo(const std::string& stopsPath) const {

    std::ifstream stops;
    stops.open(stopsPath);
    if(!stops.is_open()){
        throw std::runtime_error("Exception: cannot open file: " + stopsPath);
    }

    intStopMap uniqueStops;

    Vstring row;
    while(stops.peek() != EOF){

         row = kevintools::getNextLineAndSplitIntoTokens(stops);

        if(row.at(0) == "stop_id"){
            continue; //skip header row
        }

        int id = std::stoi(row.at(0));
        Adouble2 latlong = Adouble2{std::stod(row.at(1)), std::stod(row.at(2))};
        std::string route = "";

        auto iter = uniqueStops.find(id);
        if(iter == uniqueStops.end()){ //not found
            uniqueStops.insert(std::make_pair(id, Stop{id, latlong, Vstring{route}}));
        }
        else{ //duplicate
            Stop& stop = iter->second;
            if(latlong != std::get<1>(stop)){
                std::cerr << "Warning: id " + std::to_string(id) + " has different lat longs." << std::endl;
                //throw std::runtime_error("Exception: id " + std::to_string(id) + " has different lat longs.");
            }
            std::get<2>(stop).push_back(route);
        }

    }

    stops.close();

    VStop _allStops;
    _allStops.reserve(uniqueStops.size());
    for(auto it = uniqueStops.begin(); it != uniqueStops.end(); ++it ) {
        _allStops.push_back( it->second );
    }

    return _allStops;
}

Stops::Stops(const std::string& stopsPath) :
    allStops(readStopInfo(stopsPath))
{

}

int Stops::idToIndex(const int id) const {
    auto index = std::find_if(allStops.begin(), allStops.end(),
                              [&id](const Stop& stop) { return std::get<0>(stop) == id; });
    if(index == allStops.end()){
        throw std::runtime_error("Exception: id not valid.");
    }
    return std::distance(allStops.begin(), index);
}