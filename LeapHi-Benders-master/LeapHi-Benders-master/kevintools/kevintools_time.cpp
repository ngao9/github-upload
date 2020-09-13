#include "kevintools_time.hpp"

using namespace kevintools;

std::clock_t kevintools::startClock(){
    return std::clock();
}

std::string kevintools::timeSinceStartString(std::clock_t start){
    return " (at " + std::to_string((std::clock() - start) / (double) CLOCKS_PER_SEC) + "s)";
}