#include "trips.hpp"

#include "kevintools_functions.hpp"
#include "kevintools_strings.hpp"

#include <string>
#include <fstream>

using namespace dataprocessing;

VTrip Trips::readTripData(const std::string& tripDataPath) const {

    std::ifstream file;
    file.open(tripDataPath);
    if(!file.is_open()){
        throw std::runtime_error("Exception: cannot open file: " + tripDataPath);
    }

    VAint2 allOD = VAint2();
    VTrip _allTrips;

    //Recognize headers
    Vstring headerRow = kevintools::getNextLineAndSplitIntoTokens(file);
    std::map<std::string, int> indices;
    for(int i = 0; i < headerRow.size(); ++i){
        if(headerRow.at(i) == "start_cluster" || headerRow.at(i) == "start_cluster_id" ||
                headerRow.at(i) == "start_id" || headerRow.at(i) == "start_stop" || headerRow.at(i) == "start_stop_id"){
            indices["origin"] = i;
        }
        if(headerRow.at(i) == "end_cluster" || headerRow.at(i) == "end_cluster_id" ||
           headerRow.at(i) == "end_id" || headerRow.at(i) == "end_stop" || headerRow.at(i) == "end_stop_id"){
            indices["destination"] = i;
        }
        if(headerRow.at(i) == "count" || headerRow.at(i) == "passengers"){
            indices["count"] = i;
        }
        if(headerRow.at(i) == "start_times" || headerRow.at(i) == "times"){
            indices["times"] = i;
        }
    }
    if(indices.size() != 4){
        throw std::runtime_error("Exception: not all headers can be recognized.");
    }
    if(indices["times"] != headerRow.size()-1){
        throw std::runtime_error("Exception: the departure times column is not the last columns.");
    }
    
    //Read trips
    int nbPeople = 0;
    while(file.peek() != EOF){

        Vstring row = kevintools::getNextLineAndSplitIntoTokens(file);

        int origin = std::stoi(row.at(indices.at("origin")));
        int destination = std::stoi(row.at(indices.at("destination")));
        int count = std::stoi(row.at(indices.at("count")));

        Vstring timestamps = Vstring(row.size()-indices["times"]);
        for(int i = 0; i < timestamps.size(); ++i){
            timestamps.at(i) = row.at(i+indices["times"]);
            kevintools::removeCharsFromString(timestamps.at(i), "\"'[]");
            kevintools::trimWhitespace(timestamps.at(i));
        }

        allOD.push_back(Aint2{origin, destination});
        _allTrips.push_back(Trip{origin, destination, count, timestamps});

        nbPeople += count;

    }
    file.close();

    Aint2USet uniqueOD = Aint2USet(allOD.begin(), allOD.end());
    if(uniqueOD.size() != allOD.size()){
        throw std::runtime_error("Exception: origin and destination pairs not unique.");
    }

    int sum = 0;
    for(Trip trip : _allTrips){
        sum += std::get<2>(trip);
    }

    if(sum != nbPeople){
        throw std::runtime_error("Exception: total number of people in _allTrips incorrect.");
    }

	return _allTrips;

}

Trips::Trips(const std::string& tripDataPath) :
        allTrips(readTripData(tripDataPath))
{

}