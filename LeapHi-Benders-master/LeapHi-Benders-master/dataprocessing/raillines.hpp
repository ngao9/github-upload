#pragma once

#include "dataprocessing_typedefs.hpp"

namespace dataprocessing {

    class RailLines {

    private:
        VRailLine readLineInfo(const std::string& railLinesPath) const;

    public:
        RailLines(const std::string& railLinesPath);

    //only public for easy access
    public:
        const VRailLine allLines; //[name, nbPerHour, stopIDs, forwardSchedule, backwardSchedule]

    };

}