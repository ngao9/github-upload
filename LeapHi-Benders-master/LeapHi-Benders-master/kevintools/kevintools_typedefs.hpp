#pragma once

#include <vector>
#include <memory>
#include <map>

namespace kevintools{

    //Basic typedefs
    typedef std::vector<int> Vint;
    typedef std::vector<double> Vdouble;
    typedef std::vector<Vdouble> VVdouble;

    namespace constraint_sense_wrapper {
        enum constraint_sense {
            LE, GE, EQ
        };
    }

}