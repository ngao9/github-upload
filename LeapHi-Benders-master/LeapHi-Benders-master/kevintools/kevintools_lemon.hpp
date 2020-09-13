#pragma once

#include "kevintools_typedefs.hpp"

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/static_graph.h>
#include <lemon/adaptors.h>
#include <lemon/network_simplex.h>
#include <lemon/dijkstra.h>
#include <lemon/bfs.h>

typedef lemon::ListDigraph ListDigraph;
typedef std::vector<ListDigraph::Node> ListDigraph_VNode;

namespace kevintools{

    bool contains(const ListDigraph_VNode& vector, const ListDigraph::Node& element);

}