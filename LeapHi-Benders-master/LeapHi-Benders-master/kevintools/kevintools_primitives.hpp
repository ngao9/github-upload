#pragma once

#include <limits>

namespace kevintools{

    constexpr int POS_INF_INT = std::numeric_limits<int>::max();
    constexpr int NEG_INF_INT = std::numeric_limits<int>::min();

    inline bool is_inf(int x){
        return (x == POS_INF_INT || x == NEG_INF_INT);
    }

    inline bool is_posinf(int x){
        return x == POS_INF_INT;
    }

    inline bool is_neginf(int x){
        return x == NEG_INF_INT;
    }

    constexpr double POS_INF_DOUBLE = std::numeric_limits<double>::infinity();
    constexpr double NEG_INF_DOUBLE = -std::numeric_limits<double>::infinity();
    constexpr double NAN_DOUBLE = std::numeric_limits<double>::quiet_NaN();

    inline bool is_inf(double x){
        return (x == POS_INF_DOUBLE || x == NEG_INF_DOUBLE);
    }

    inline bool is_posinf(double x){
        return x == POS_INF_DOUBLE;
    }

    inline bool is_neginf(double x){
        return x == NEG_INF_DOUBLE;
    }

}