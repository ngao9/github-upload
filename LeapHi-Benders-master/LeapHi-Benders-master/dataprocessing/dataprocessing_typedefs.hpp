#pragma once

#include "kevintools_functions.hpp"

#include <vector>
#include <array>
#include <map>
#include <unordered_set>
#include <memory>

namespace dataprocessing{

    //Forward declarations
    class RailLines;
    class Stops;
    class Trips;
    class Distances;

    //Basic typedefs
    typedef std::vector<double> Vdouble;
    typedef std::vector<Vdouble> VVdouble;

    typedef std::vector<int> Vint;

    typedef std::array<int, 2> Aint2;
    typedef std::array<int, 3> Aint3;

    typedef std::vector<Aint2> VAint2;
    typedef std::unordered_set<Aint2, kevintools::Aint2_hasher> Aint2USet;

    typedef std::array<double, 2> Adouble2;

    typedef std::vector<std::string> Vstring;


    //Smart pointers
    typedef std::shared_ptr<RailLines> pRailLines;
    typedef std::shared_ptr<Stops> pStops;
    typedef std::shared_ptr<Trips> pTrips;
    typedef std::shared_ptr<Distances> pDistances;


    //Special typedefs
    typedef std::tuple<std::string, int, Vint, Vdouble, Vdouble> RailLine;
    typedef std::vector<RailLine> VRailLine;
    typedef std::map<std::string, RailLine> stringRailLineMap;

    typedef std::tuple<int, Adouble2, Vstring> Stop;
    typedef std::vector<Stop> VStop;
    typedef std::map<int, Stop> intStopMap;

    typedef std::tuple<int, int, int, Vstring> Trip;
    typedef std::vector<Trip> VTrip;

    typedef std::vector<std::string> Vstring;

}