#pragma once

#include "dataprocessing_typedefs.hpp"

namespace dataprocessing {

    class Stops {

    private:
        VStop readStopInfo(const std::string& stopsPath) const;

    public:
        Stops(const std::string& stopsPath);
        int idToIndex(int id) const; //-1 if not found

    //only public for easy access
    public:
        const VStop allStops; //[id, latlong, route]

    };

}