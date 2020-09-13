#pragma once

#include "dataprocessing_typedefs.hpp"

namespace dataprocessing {

    class Distances {

    private:
        VVdouble readDistanceData(const std::string& distancesPath) const;

    public:
        Distances(const std::string& distancesPath, const Stops& stops);
        double getMaxTriangleInequalityViolation(Aint3& indices) const;

    private:
        const Stops& stops;

    //only public for easy access
    public:
        const VVdouble allDistances; //use same index as stops

    };

}