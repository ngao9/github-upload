#pragma once

#include "dataprocessing_typedefs.hpp"

namespace dataprocessing {

    class Trips {

    private:
        VTrip readTripData(const std::string& tripDataPath) const;

    public:
        Trips(const std::string& tripDataPath);

    //only public for easy access
    public:
        const VTrip allTrips; //[origin, destination, nb, departure timestamps (Vstring)]

    };

}