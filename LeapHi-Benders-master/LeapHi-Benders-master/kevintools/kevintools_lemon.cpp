#include "kevintools_lemon.hpp"

using namespace kevintools;

bool kevintools::contains(const ListDigraph_VNode& vector, const ListDigraph::Node& element){
    return std::find(vector.begin(), vector.end(), element) != vector.end();
}