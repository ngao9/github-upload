#include "raillines.hpp"

#include "kevintools_functions.hpp"

#include <fstream>
#include <string>

using namespace dataprocessing;

VRailLine RailLines::readLineInfo(const std::string& railLinesPath) const {

    std::ifstream railLinesFile;
    railLinesFile.open(railLinesPath);
    if(!railLinesFile.is_open()){
        throw std::runtime_error("Exception: cannot open file: " + railLinesPath);
    }

    stringRailLineMap uniqueLines = stringRailLineMap();

    while(railLinesFile.peek() != EOF){

        Vstring row = kevintools::getNextLineAndSplitIntoTokens(railLinesFile);

        std::string name = row.at(0);

        if(name == "NAME"){
            continue; //skip header row
        }

        int nb = std::stoi(row.at(1));
        Vint stopIDs;
        for(int k = 2; k < row.size(); ++k){
            if(row.at(k) != "") {
                stopIDs.push_back(std::stoi(row.at(k)));
            }
        }

        Vdouble forwardSchedule;
        row = kevintools::getNextLineAndSplitIntoTokens(railLinesFile);
        if(row.at(0) != "schedule_left_right"){
            throw std::runtime_error("Exception: missing schedule_left_right row.");
        }
        for(int k = 2; k < row.size(); ++k){
            if(row.at(k) != "") {
                forwardSchedule.push_back(std::stod(row.at(k)));
            }
        }

        Vdouble backwardSchedule;
        row = kevintools::getNextLineAndSplitIntoTokens(railLinesFile);
        if(row.at(0) != "schedule_right_left"){
            throw std::runtime_error("Exception: missing schedule_right_left row.");
        }
        for(int k = 2; k < row.size(); ++k){
            if(row.at(k) != "") {
                backwardSchedule.push_back(std::stod(row.at(k)));
            }
        }

        if(uniqueLines.count(name) != 0){
            throw std::runtime_error("Exception: rail line name " + name + " is not unique.");
        }
        uniqueLines.insert(std::make_pair(name, RailLine{name, nb, stopIDs, forwardSchedule, backwardSchedule}));

    }

    railLinesFile.close();

    VRailLine _allLines;
    _allLines.reserve(uniqueLines.size());
    for(auto it = uniqueLines.begin(); it != uniqueLines.end(); ++it ) {
        _allLines.push_back( it->second );
    }

    return _allLines;

}

RailLines::RailLines(const std::string& railLinesPath) :
    allLines(readLineInfo(railLinesPath))
{

}