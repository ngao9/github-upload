#pragma once

#include "networksolver_typedefs.hpp"
#include "kevintools_typedefs.hpp"

#include "network.hpp"
#include "lemonnetwork.hpp"
#include "supergraph.hpp"

#ifdef GUROBI
#include "gurobinetwork.hpp"
#endif

using namespace kevintools::constraint_sense_wrapper;