#include "distances.hpp"
#include "stops.hpp"

#include <fstream>
#include <limits>
#include <cmath>
#include <string>
#include <iostream>

#include "kevintools_primitives.hpp"

dataprocessing::VVdouble dataprocessing::Distances::readDistanceData(const std::string &path) const {

    const VStop& allStops = stops.allStops;

    VVdouble rawData = VVdouble();

    std::ifstream file;
    file.open(path);
    if(!file.is_open()){
        throw std::runtime_error("Exception: cannot open file: " + path);
    }

    while(file.peek() != EOF){

        Vstring stringRow = kevintools::getNextLineAndSplitIntoTokens(file);
        const int& m = stringRow.size();

        rawData.push_back(Vdouble(m));
        Vdouble& doubleRow = rawData.at(rawData.size() - 1);

        for(int j = 0; j < m; ++j){

            const std::string& stringVal = stringRow.at(j);

            if(stringVal == "" || stringVal == "stop_id"){
                doubleRow.at(j) = std::numeric_limits<double>::quiet_NaN();
            }
            else{
                doubleRow.at(j) = std::stod(stringVal);
            }

        }

    }

    int nRaw = rawData.size();
    int mRaw = rawData.at(0).size();

    if(nRaw != mRaw){
        throw std::runtime_error("nRaw != mRaw");
    }

    std::map<int, int> idToIndexMap = std::map<int, int>();

    for(int i = 1; i < nRaw; ++i) {

        int rowId = (int) std::round(rawData.at(i).at(0));
        int colId = (int) std::round(rawData.at(0).at(i));

        if(rowId != colId){
            throw std::runtime_error("rowId != colId");
        }

        idToIndexMap.insert(std::make_pair(rowId, i));

    }

    const int& nbStops = allStops.size();
    VVdouble distances = VVdouble(nbStops, Vdouble(nbStops, 0));

    for(int i = 0; i < nbStops; ++i){

        const Stop& fromStop = allStops.at(i);
        const int& fromId = std::get<0>(fromStop);
        int fromIndex = idToIndexMap.at(fromId);

        for(int j = 0; j < nbStops; ++j){

            const Stop& toStop = allStops.at(j);
            const int& toId = std::get<0>(toStop);
            int toIndex = idToIndexMap.at(toId);

            distances.at(i).at(j) = rawData.at(fromIndex).at(toIndex);

            if(distances.at(i).at(j) < 0){
                std::cerr << "Warning: missing value for fromId = " + std::to_string(fromId) + ", toId = " << std::to_string(toId) << ", value set to NaN." << std::endl;
                distances.at(i).at(j) = kevintools::NAN_DOUBLE;
            }

        }

    }

    return distances;

}

dataprocessing::Distances::Distances(const std::string &distancesPath, const dataprocessing::Stops &stops) :
    stops(stops), allDistances(readDistanceData(distancesPath))
{

}

double dataprocessing::Distances::getMaxTriangleInequalityViolation(Aint3& indices) const {

    double violationMax = -std::numeric_limits<double>::max();
    indices.at(0) = -1;
    indices.at(1) = -1;
    indices.at(2) = -1;

    int n = allDistances.size();

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                double violation = allDistances.at(i).at(k) - allDistances.at(i).at(j) - allDistances.at(j).at(k);
                if(violation > violationMax){
                    violationMax = violation;
                    indices.at(0) = i;
                    indices.at(1) = j;
                    indices.at(2) = k;
                }
            }
        }
    }

    return violationMax;

}
