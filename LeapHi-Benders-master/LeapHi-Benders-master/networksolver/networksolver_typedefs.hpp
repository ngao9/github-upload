#pragma once

#include "kevintools_gurobi.hpp"
#include "kevintools_lemon.hpp"
#include "kevintools_typedefs.hpp"

#include <memory>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>

namespace networksolver {

    //Forward declarations
    class Network;
#ifdef GUROBI
    class GurobiNetwork;
#endif
    class LemonNetwork;
    class SuperGraph;

    typedef std::shared_ptr<Network> pNetwork;
    typedef std::vector<pNetwork> VpNetwork;
#ifdef GUROBI
    typedef std::shared_ptr<GurobiNetwork> pGurobiNetwork;
    typedef std::vector<pGurobiNetwork> VpGurobiNetwork;
#endif
    typedef std::shared_ptr<LemonNetwork> pLemonNetwork;
    typedef std::vector<pLemonNetwork> VpLemonNetwork;
    typedef std::shared_ptr<SuperGraph> pSuperGraph;


#ifdef GUROBI
    //Gurobi typdefs
    typedef std::vector<GRBVar> VGRBVar;
    typedef std::vector<VGRBVar> VVGRBVar;
    typedef std::vector<VVGRBVar> VVVGRBVar;

    typedef std::shared_ptr<GRBEnv> pGRBEnv;
    typedef std::shared_ptr<GRBModel> pGRBModel;

    typedef std::vector<GRBConstr> VGRBConstr;
    typedef std::vector<VGRBConstr> VVGRBConstr;

    typedef std::vector<GRBLinExpr> VGRBLinExpr;
    typedef std::vector<VGRBLinExpr> VVGRBLinExpr;
#endif


    //Lemon typedefs
    typedef lemon::ListDigraph ListDigraph;
    typedef std::vector<ListDigraph::Arc> ListDigraph_VArc;
    typedef std::vector<ListDigraph_VArc> ListDigraph_VVArc;
    typedef std::vector<ListDigraph::Node> ListDigraph_VNode;
    typedef std::shared_ptr<ListDigraph::ArcMap<double>> pListDigraph_ArcMap_double;

    typedef lemon::StaticDigraph StaticDigraph;
    typedef std::vector<StaticDigraph::Arc> StaticDigraph_VArc;
    typedef std::vector<StaticDigraph::Node> StaticDigraph_VNode;
    typedef std::shared_ptr<StaticDigraph::ArcMap<int>> pStaticDigraph_ArcMap_int;
    typedef std::shared_ptr<StaticDigraph::ArcMap<double>> pStaticDigraph_ArcMap_double;
    typedef std::shared_ptr<StaticDigraph::NodeMap<int>> pStaticDigraph_NodeMap_int;
    typedef std::shared_ptr<StaticDigraph::ArcMap<bool>> pStaticDigraph_ArcMap_bool;

    typedef lemon::NetworkSimplex<StaticDigraph> NetworkSimplex;
    typedef std::shared_ptr<NetworkSimplex> pNetworkSimplex;
    typedef lemon::FilterArcs<const StaticDigraph> FilterArcs;
    typedef std::shared_ptr<FilterArcs> pFilterArcs;
    typedef lemon::Dijkstra<FilterArcs, StaticDigraph::ArcMap<double>> FilterArcsDijkstra;
    typedef std::shared_ptr<FilterArcsDijkstra> pFilterArcsDijkstra;


    //Basic typedefs
    typedef std::vector<double> Vdouble;
    typedef std::vector<Vdouble> VVdouble;

    typedef std::vector<int> Vint;
    typedef std::vector<Vint> VVint;
    typedef std::unordered_map<int,int> intintMap;

    typedef std::vector<bool> Vbool;

    typedef std::array<int, 2> Aint2;
    typedef std::array<int, 3> Aint3;
    typedef std::array<int, 4> Aint4;

    typedef std::array<double, 2> Adouble2;
    typedef std::vector<Adouble2> VAdouble2;

    typedef std::set<int> intSet;


    //Special typedefs
    typedef Aint2 ResourceContribution;
    typedef std::vector<ResourceContribution> VResourceContribution;
    typedef std::vector<VResourceContribution> VVResourceContribution;

    typedef Adouble2 ResourceBounds;
    typedef std::vector<ResourceBounds> VResourceBounds;

    using namespace kevintools::constraint_sense_wrapper;
#ifdef GUROBI
    typedef kevintools::SenseConstraintMap_GRB SenseConstraintMap_GRB;
    typedef kevintools::VSenseConstraintMap_GRB VSenseConstraintMap_GRB;
    typedef kevintools::VVSenseConstraintMap_GRB VVSenseConstraintMap_GRB;
#endif

}