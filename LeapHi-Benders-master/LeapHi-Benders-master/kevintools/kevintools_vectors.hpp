#pragma once

#include "kevintools_typedefs.hpp"

#include <vector>
#include <cmath>
#include <algorithm>

namespace kevintools{

    template <typename T>
    inline std::vector<T> operator+(const std::vector<T>& v, const T constant){

        std::vector<T> result = std::vector<T>(v.begin(), v.end());
        for(T& elem : result){
            elem += constant;
        }

        return result;

    }

    template <typename T>
    inline std::vector<T> operator*(const std::vector<T>& v, const T constant){

        std::vector<T> result = std::vector<T>(v.begin(), v.end());
        for(T& elem : result){
            elem *= constant;
        }

        return result;

    }

    template <typename T>
    inline Vint castToInt(const std::vector<T>& v){
        return Vint(v.begin(), v.end());
    }

    template <typename T>
    inline Vdouble castToDouble(const std::vector<T>& v){
        return Vdouble(v.begin(), v.end());
    }

    template <typename T>
    inline std::vector<T> removeElementsFromVector(const std::vector<T>& elements, const std::vector<T>& vector){

        std::vector<T> result = vector;

        for (const T& element : elements){
            result.erase(std::find(result.begin(), result.end(), element));
        }

        return result;

    }

    bool almostEqual(const Vdouble& v1, const Vdouble& v2, const double eps);
    bool almostEqual(const VVdouble& v1, const VVdouble& v2, const double eps);
    bool almostLessOrEqual(const Vdouble& v1, const Vdouble& v2, const double eps);

}