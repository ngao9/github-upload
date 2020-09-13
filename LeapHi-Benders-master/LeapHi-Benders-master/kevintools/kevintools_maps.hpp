#pragma once

#include <string>
#include <stdexcept>
#include <cmath>

#include "kevintools_typedefs.hpp"

namespace kevintools{

    template <typename MapType, typename KeyType, typename ValueType>
    inline void safeInsert(MapType& map, const KeyType& key, const ValueType& value){

        auto iterSuccess = map.insert(std::make_pair(key, value));
        const auto& iter = iterSuccess.first;
        const bool& success = iterSuccess.second;

        if(!success) {

            const ValueType& currentValue = iter->second;

            if(currentValue != value){
                std::string message = "Exception: safeInsert tried to overwrite a previously assigned value by a different value.";
                throw std::runtime_error(message);
            }

        }

    }

    inline void safeInsert(std::vector<std::vector<std::vector<double>>>& matrix,
            const int i, const int j, const int k, const double value) {
        const double& oldValue = matrix.at(i).at(j).at(k);
        if (!std::isnan(oldValue) && oldValue != value) {
            std::string message = "Exception: safeInsert tried to overwrite a previously assigned value by a different value.";
            throw std::runtime_error(message);
        }
        matrix.at(i).at(j).at(k) = value;
    };

}