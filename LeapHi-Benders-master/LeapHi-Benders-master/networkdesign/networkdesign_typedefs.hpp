#pragma once

#include "networksolver_typedefs.hpp"
#include "dataprocessing_typedefs.hpp"
#include "kevintools_typedefs.hpp"
#include "kevintools_cplex.hpp"
#include "kevintools_functions.hpp"

#include <unordered_set>
#include <unordered_map>
#include <random>
#include <set>
#include <queue>
#include <functional>

namespace networkdesign {

    //Constants
    constexpr double mile2km = 1.609344;


    //Forward declarations
    class Callbacks;
    class Instance;
    class NetworkDesignProblem;
    class Runner;


    //Basic typedefs
    typedef std::vector<double> Vdouble;
    typedef std::vector<Vdouble> VVdouble;
    typedef std::vector<VVdouble> VVVdouble;

    typedef std::vector<int> Vint;
    typedef std::vector<Vint> VVint;

    typedef std::set<int> intSet;

    typedef std::array<int, 2> Aint2;
    typedef std::array<int, 3> Aint3;
    typedef std::array<int, 4> Aint4;
    typedef std::array<int, 5> Aint5;

    typedef std::vector<Aint2> VAint2;
    typedef std::vector<Aint3> VAint3;
    typedef std::vector<Aint4> VAint4;
    typedef std::vector<Aint5> VAint5;

    typedef std::unordered_map<Aint3, double, kevintools::Aint3_hasher> Aint3doubleMap;
    typedef std::unordered_map<Aint3, Vint, kevintools::Aint3_hasher> Aint3VintMap;
    typedef std::unordered_set<Aint4, kevintools::Aint4_hasher> Aint4USet;
    typedef std::unordered_set<Aint5, kevintools::Aint5_hasher> Aint5USet;

    typedef std::array<double, 2> Adouble2;
    typedef std::vector<Adouble2> VAdouble2;


    //Smart pointers
    typedef std::shared_ptr<Instance> pInstance;
    typedef std::shared_ptr<NetworkDesignProblem> pNetworkDesignProblem;


    //CPLEX typdefs
    typedef std::vector<IloNumVar> VIloNumVar;
    typedef std::vector<VIloNumVar> VVIloNumVar;
    typedef std::vector<VVIloNumVar> VVVIloNumVar;

    typedef std::vector<IloIntVar> VIloIntVar;

    typedef std::shared_ptr<IloEnv> pIloEnv;
    typedef std::shared_ptr<IloModel> pIloModel;
    typedef std::shared_ptr<IloCplex> pIloCplex;

    typedef std::vector<IloExpr> VIloExpr;
    typedef std::vector<VIloExpr> VVIloExpr;
    typedef std::vector<IloExprArray> VIloExprArray;

    typedef std::vector<IloRange> VIloRange;
    typedef std::vector<VIloRange> VVIloRange;
    typedef std::vector<VVIloRange> VVVIloRange;

#ifdef GUROBI
    //Gurobi typedefs
    typedef std::shared_ptr<GRBEnv> pGRBEnv;
#endif

    //Lemon typedefs
    typedef lemon::ListDigraph ListDigraph;
    typedef std::vector<ListDigraph::Arc> ListDigraph_VArc;
    typedef std::vector<ListDigraph_VArc> ListDigraph_VVArc;
    typedef std::vector<ListDigraph::Node> ListDigraph_VNode;
    typedef std::vector<ListDigraph_VNode> ListDigraph_VVNode;
    typedef std::shared_ptr<ListDigraph::ArcMap<double>> pListDigraph_ArcMap_double;


    //Special typedefs
    typedef Aint3 tripleIndex;
    typedef std::set<tripleIndex> tripleIndexSet;
    typedef std::unordered_map<tripleIndex, IloExpr, kevintools::Aint3_hasher> tripleIndexIloExprMap;

    typedef Aint3doubleMap tripleIndexDoubleMap;
    typedef Aint3VintMap tripleIndexVintMap;

    typedef dataprocessing::Stop stopInfo;
    typedef dataprocessing::VStop VstopInfo;

    typedef dataprocessing::Trip Trip;
    typedef dataprocessing::VTrip VTrip;

    typedef std::pair<double, Vdouble> Solution;

    typedef std::pair<double, Vint> IntSolution;
    typedef std::queue<IntSolution> QIntSolution;

    namespace travel_mode_wrapper {
        enum travel_mode {
            SHUTTLE, BUS, RAIL
        };
    }
    using namespace travel_mode_wrapper;

    enum SubproblemType { GUROBINETWORK, LEMONNETWORK };

    typedef std::function<std::tuple<double, Vdouble, double, Vdouble>(const double, const double, const travel_mode, const int, const int)> ArcFunction;


    ////////////////////////////
    // NetworkSolver typedefs //
    ////////////////////////////

    //Basic typedefs
    typedef networksolver::Network Network;
#ifdef GUROBI
    typedef networksolver::GurobiNetwork GurobiNetwork;
#endif
    typedef networksolver::LemonNetwork LemonNetwork;
    typedef networksolver::SuperGraph SuperGraph;

    //Smart pointers
    typedef networksolver::pNetwork pNetwork;
    typedef networksolver::VpNetwork VpNetwork;
#ifdef GUROBI
    typedef networksolver::pGurobiNetwork pGurobiNetwork;
#endif
    typedef networksolver::pSuperGraph pSuperGraph;

    //Special typedefs
    typedef networksolver::ResourceContribution ResourceContribution;
    typedef networksolver::VResourceContribution VResourceContribution;
    typedef networksolver::ResourceBounds ResourceBounds;
    typedef networksolver::VResourceBounds VResourceBounds;
    using namespace kevintools::constraint_sense_wrapper;


    /////////////////////////////
    // DataProcessing typedefs //
    /////////////////////////////

    typedef dataprocessing::RailLines RailLines;
    typedef dataprocessing::RailLine RailLine;
    typedef dataprocessing::VRailLine VRailLine;
    typedef dataprocessing::Stops Stops;
    typedef dataprocessing::Trips Trips;
    typedef dataprocessing::Distances Distances;
    typedef dataprocessing::pRailLines pRailLines;
    typedef dataprocessing::pStops pStops;
    typedef dataprocessing::pTrips pTrips;
    typedef dataprocessing::pDistances pDistances;

}