#pragma once

#include "kevintools_typedefs.hpp"
#include <ctime>
#include <string>

namespace kevintools{

    std::clock_t startClock();
    std::string timeSinceStartString(std::clock_t start);

}