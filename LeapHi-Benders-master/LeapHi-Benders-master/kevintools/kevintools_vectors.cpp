#include "kevintools_vectors.hpp"

bool kevintools::almostEqual(const Vdouble& v1, const Vdouble& v2, const double eps){

    if(v1.size() != v2.size()){
        return false;
    }

    for(int i = 0; i < v1.size(); ++i){
        if(std::fabs(v1.at(i) - v2.at(i)) > eps){
            return false;
        }
    }

    return true;

}

bool kevintools::almostEqual(const VVdouble& v1, const VVdouble& v2, const double eps){

    if(v1.size() != v2.size()){
        return false;
    }

    for(int i = 0; i < v1.size(); ++i){
        if(!almostEqual(v1.at(i), v2.at(i), eps)){
            return false;
        }
    }

    return true;

}

bool kevintools::almostLessOrEqual(const Vdouble& v1, const Vdouble& v2, const double eps){

    if(v1.size() != v2.size()){
        return false;
    }

    for(int i = 0; i < v1.size(); ++i){
        if(v1.at(i) > v2.at(i) + eps){
            return false;
        }
    }

    return true;

}