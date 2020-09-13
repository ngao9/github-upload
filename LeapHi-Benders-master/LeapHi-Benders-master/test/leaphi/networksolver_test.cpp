#include "networksolver_cpp_include.hpp"

#include "kevintools_primitives.hpp"
#include "kevintools_gurobi.hpp"
#include "kevintools_vectors.hpp"
#include "kevintools_strings.hpp"

#include <cmath>
#include <algorithm>

#include "gtest/gtest.h"

using namespace networksolver;

typedef std::tuple<pNetwork, pSuperGraph, ListDigraph_VArc, VResourceBounds> NetworkData;
typedef std::tuple<VpNetwork, pSuperGraph, ListDigraph_VArc, VResourceBounds> MultiNetworkData;

ListDigraph::Arc addArc(ListDigraph& graph, int i, int j){
    return graph.addArc(graph.nodeFromId(i), graph.nodeFromId(j));
};

#ifdef GUROBI
NetworkData testGurobiNetwork1(){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(6);
    ListDigraph& graph = superGraph->graph;

    ListDigraph_VArc arcs;
    arcs.push_back(addArc(graph, 0, 3));
    arcs.push_back(addArc(graph, 1, 0));
    arcs.push_back(addArc(graph, 1, 3));
    arcs.push_back(addArc(graph, 1, 4));
    arcs.push_back(addArc(graph, 1, 0));
    arcs.push_back(addArc(graph, 5, 0));

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    ListDigraph_VNode sourceNodes = ListDigraph_VNode();
    ListDigraph_VNode sinkNodes = ListDigraph_VNode();

    pNetwork network = std::make_shared<GurobiNetwork>(superGraph, nodeFilter, arcFilter,
                                                 sourceNodes, Vdouble{}, sinkNodes, Vdouble{},
                                                 0, std::make_shared<GRBEnv>(), VResourceBounds());

    return std::make_tuple(network, superGraph, arcs, VResourceBounds());

}
#endif

#ifdef GUROBI
//https://imada.sdu.dk/~jbj/DM85/mincostnew.pdf
NetworkData testGurobiNetwork2(){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(6);
    ListDigraph& graph = superGraph->graph;

    VVint arcsIndex = {{0,1}, {1,2}, {3,0}, {1,4}, {2,5}, {3,4}, {4,5}};
    Vdouble costs = {2, 3, 1, 1, 1, 5, 4};
    Vint capacities = {10, 8, 15, 15, 15, 15, 15};

    ListDigraph_VArc arcs;
    int nbArcs = arcsIndex.size();
    for(int k = 0; k < nbArcs; ++k) {

        int i = arcsIndex.at(k).at(0);
        int j = arcsIndex.at(k).at(1);
        ListDigraph::Arc arc = addArc(graph, i, j);
        arcs.push_back(arc);

        superGraph->cost[arc] = costs.at(k);

        superGraph->contributions[arc] = VResourceContribution{{k, 1}};

        //Additional constraint resource nbArcs
        if (k == 0) {
            superGraph->contributions[arc].push_back({nbArcs, 3});
        } else if (k == 3) {
            superGraph->contributions[arc].push_back({nbArcs, 5});
        } else if (k == 6) {
            superGraph->contributions[arc].push_back({nbArcs, -1});
        }

        //Lowerbound on arc (3, 4) with resource nbArcs + 1
        if (k == 5) {
            superGraph->contributions[arc].push_back({nbArcs + 1, 1});
        }

    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    ListDigraph_VNode sourceNodes = {graph.nodeFromId(0), graph.nodeFromId(3)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(2), graph.nodeFromId(5)};

    //Resources [0, nbArcs-1] are capacities
    //Resources [nbArcs, nbArcs+1] are additional constraints
    VResourceBounds resources = VResourceBounds();
    resources.reserve(nbArcs + 1);
    for(int k = 0; k < nbArcs; ++k){
        resources.push_back(Adouble2{kevintools::NEG_INF_DOUBLE, (double) capacities.at(k)});
    }
    resources.push_back(Adouble2{33, 33});
    resources.push_back(Adouble2{0, std::numeric_limits<double>::infinity()});

    pNetwork network = std::make_shared<GurobiNetwork>(superGraph, nodeFilter, arcFilter,
                                                 sourceNodes, Vdouble{5, 10}, sinkNodes, Vdouble{-5, -10},
                                                 0, std::make_shared<GRBEnv>(), resources);

    return std::make_tuple(network, superGraph, arcs, resources);

}
#endif

#ifdef GUROBI
NetworkData testGurobiNetwork3(const int maxCardinality){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(4);
    ListDigraph& graph = superGraph->graph;

    VVint arcsIndex = {{0,1}, {0,2}, {3, 0}, {0,3}, {1,2}, {1,3}, {2,3}};
    Vdouble costs = {1, 3, 0, 10, -1, 2, 1};
    Vdouble capacities = {0.9, 1, 0, 1, 1, 1, 1};

    ListDigraph_VArc arcs;
    int nbArcs = arcsIndex.size();
    for(int k = 0; k < nbArcs; ++k){

        int i = arcsIndex.at(k).at(0);
        int j = arcsIndex.at(k).at(1);
        ListDigraph::Arc arc = addArc(graph, i, j);
        arcs.push_back(arc);

        superGraph->cost[arc] = costs.at(k);
        superGraph->contributions[arc] = {{k, 1}};

    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    arcFilter[arcs.at(2)] = false; //filter out arc (3, 0)
    ListDigraph_VNode sourceNodes = {graph.nodeFromId(0)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(3)};

    //Resources [0, nbArcs-1] are capacities
    //Resource [nbArcs] is not used by any arc
    VResourceBounds resources = VResourceBounds();
    for(int k = 0; k < nbArcs; ++k){
        resources.push_back(Adouble2{kevintools::NEG_INF_DOUBLE, capacities.at(k)});
    }

    pNetwork network = std::make_shared<GurobiNetwork>(superGraph, nodeFilter, arcFilter,
                                                 sourceNodes, Vdouble{1}, sinkNodes, Vdouble{-1},
                                                 maxCardinality, std::make_shared<GRBEnv>(), resources);

    return std::make_tuple(network, superGraph, arcs, resources);

}
#endif

#ifdef GUROBI
MultiNetworkData testGurobiMultiNetwork1(const int maxCardinality){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(2);
    ListDigraph& graph = superGraph->graph;

    ListDigraph_VArc arcs;
    arcs.push_back(addArc(graph, 0, 1));
    arcs.push_back(addArc(graph, 0, 1));
    arcs.push_back(addArc(graph, 1, 0));
    arcs.push_back(addArc(graph, 1, 0));

    const ListDigraph::Arc& a0_0 = arcs.at(0);
    const ListDigraph::Arc& a0_1 = arcs.at(1);
    const ListDigraph::Arc& a1_0 = arcs.at(2);
    const ListDigraph::Arc& a1_1 = arcs.at(3);

    superGraph->cost[a0_0] = 10;
    superGraph->cost[a0_1] = 1;
    superGraph->cost[a1_0] = 5;
    superGraph->cost[a1_1] = 1;

    superGraph->contributions[a0_0] = {{0,1}};
    superGraph->contributions[a0_1] = {{1,1}};
    superGraph->contributions[a1_0] = {{2,1}};
    superGraph->contributions[a1_1] = {{3,1}};

    VResourceBounds resourceBounds = VResourceBounds{
            ResourceBounds{-1, 1}, //-1 instead of -inf to allow for changing the bound
            ResourceBounds{-1, 1},
            ResourceBounds{-1, 1},
            ResourceBounds{-std::numeric_limits<double>::infinity(), 1}
    };

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter1(graph, false);
    arcFilter1[a0_0] = true;
    arcFilter1[a0_1] = true;
    ListDigraph::ArcMap<bool> arcFilter2(graph, false);
    arcFilter2[a1_0] = true;
    arcFilter2[a1_1] = true;

    ListDigraph_VNode sourceNodes1 = {graph.nodeFromId(0)};
    ListDigraph_VNode sinkNodes1 = {graph.nodeFromId(1)};
    ListDigraph_VNode sourceNodes2 = {graph.nodeFromId(1)};
    ListDigraph_VNode sinkNodes2 = {graph.nodeFromId(0)};

    pNetwork network1 = std::make_shared<GurobiNetwork>(superGraph, nodeFilter, arcFilter1,
                                                  sourceNodes1, Vdouble{1}, sinkNodes1, Vdouble{-1},
                                                  maxCardinality, std::make_shared<GRBEnv>(), resourceBounds);
    pNetwork network2 = std::make_shared<GurobiNetwork>(superGraph, nodeFilter, arcFilter2,
                                                  sourceNodes2, Vdouble{1}, sinkNodes2, Vdouble{-1},
                                                  maxCardinality, std::make_shared<GRBEnv>(), resourceBounds);

    network1->writeModelToFile("TEST_testGurobiMultiNetwork1_network1.lp");
    network2->writeModelToFile("TEST_testGurobiMultiNetwork1_network2.lp");

    return std::make_tuple(VpNetwork{network1, network2}, superGraph, arcs, resourceBounds);

}
#endif

NetworkData testLemonNetwork1(int costScaleFactor, int capacityScaleFactor){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(4);
    ListDigraph& graph = superGraph->graph;

    VVint arcsIndex = {{0,1}, {1,2}, {1,2}, {1,2}, {2,3}, {0,3}};
    Vdouble costs = {0, 3.71, 2.61, 0.99, 0, 10.08};
    Vdouble capacities = {0, 0.33, 0.6, 0.5, 2, 2};

    ListDigraph_VArc arcs;
    int nbArcs = arcsIndex.size();
    for(int k = 0; k < nbArcs; ++k){

        int i = arcsIndex.at(k).at(0);
        int j = arcsIndex.at(k).at(1);
        ListDigraph::Arc arc = addArc(graph, i, j);
        arcs.push_back(arc);

        superGraph->cost[arc] = costs.at(k);
        superGraph->contributions[arc] = {{k, 1}};

    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    ListDigraph_VNode sourceNodes = {graph.nodeFromId(0)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(3)};

    //Resources [0, nbArcs-1] are capacities
    //Resource [nbArcs] is not used by any arc
    VResourceBounds resources = VResourceBounds();
    for(int k = 0; k < nbArcs; ++k){
        resources.push_back(Adouble2{kevintools::NEG_INF_DOUBLE, capacities.at(k)});
    }

    int maxCardinality = 0;

    pNetwork network = std::make_shared<LemonNetwork>(superGraph, nodeFilter, arcFilter,
                                                       sourceNodes, Vdouble{2}, sinkNodes, Vdouble{-2},
                                                       maxCardinality, costScaleFactor, capacityScaleFactor);

    return std::make_tuple(network, superGraph, arcs, resources);

}

NetworkData testLemonNetwork2(const int maxCardinality){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(4);
    ListDigraph& graph = superGraph->graph;

    VVint arcsIndex = {{0,1}, {0,2}, {3, 0}, {0,3}, {1,2}, {1,3}, {2,3}};
    Vdouble costs = {1, 3, 0, 10, -1, 2, 1};
    Vdouble capacities = {0.9, 1, 0, 1, 1, 1, 1};

    ListDigraph_VArc arcs;
    int nbArcs = arcsIndex.size();
    for(int k = 0; k < nbArcs; ++k){

        int i = arcsIndex.at(k).at(0);
        int j = arcsIndex.at(k).at(1);
        ListDigraph::Arc arc = addArc(graph, i, j);
        arcs.push_back(arc);

        superGraph->cost[arc] = costs.at(k);
        superGraph->contributions[arc] = {{k, 1}};

    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    arcFilter[arcs.at(2)] = false; //filter out arc (3, 0)
    ListDigraph_VNode sourceNodes = {graph.nodeFromId(0)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(3)};

    //Resources [0, nbArcs-1] are capacities
    //Resource [nbArcs] is not used by any arc
    VResourceBounds resources = VResourceBounds();
    for(int k = 0; k < nbArcs; ++k){
        resources.push_back(Adouble2{kevintools::NEG_INF_DOUBLE, capacities.at(k)});
    }

    pNetwork network = std::make_shared<LemonNetwork>(superGraph, nodeFilter, arcFilter,
                                                       sourceNodes, Vdouble{1}, sinkNodes, Vdouble{-1},
                                                       maxCardinality, 2, 10);
    bool requireTightBound = false;
    network->updateResourceBounds(resources, requireTightBound);

    return std::make_tuple(network, superGraph, arcs, resources);

}

NetworkData testLemonNetwork3(){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(5);
    ListDigraph& graph = superGraph->graph;

    VVint arcsIndex = {{0,3}, {0,1}, {1,3}, {0,2}, {2,3}, {0,4}};
    Vdouble costs = {4.7, 2.01, 3.2, 2, 1.3, 10};
    Vdouble capacities = {kevintools::NAN_DOUBLE, 1, 1, 1, 1, 1, 1}; //no capacity constraint first arc

    ListDigraph_VArc arcs;
    int nbArcs = arcsIndex.size();
    for(int k = 0; k < nbArcs; ++k){

        int i = arcsIndex.at(k).at(0);
        int j = arcsIndex.at(k).at(1);
        ListDigraph::Arc arc = addArc(graph, i, j);
        arcs.push_back(arc);

        superGraph->cost[arc] = costs.at(k);
        if(k != 0) { //resource 0 not in superGraph
            superGraph->contributions[arc] = {{k, 1}};
        }

    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    ListDigraph_VNode sourceNodes = {graph.nodeFromId(0)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(3)};

    //Resources [0, nbArcs-1] are capacities
    VResourceBounds resources = VResourceBounds();
    for(int k = 0; k < nbArcs; ++k){
        resources.push_back(Adouble2{kevintools::NEG_INF_DOUBLE, capacities.at(k)});
    }

    int maxCardinality = 0;
    int costScaleFactor = 3;
    int capacityScaleFactor = 2;

    pNetwork network = std::make_shared<LemonNetwork>(superGraph, nodeFilter, arcFilter,
                                                      sourceNodes, Vdouble{1}, sinkNodes, Vdouble{-1},
                                                      maxCardinality, costScaleFactor, capacityScaleFactor);

    return std::make_tuple(network, superGraph, arcs, resources);

}

NetworkData testLemonNetwork4(){

    pSuperGraph superGraph = std::make_shared<SuperGraph>(5);
    ListDigraph& graph = superGraph->graph;

    VVint arcsIndex = {{1,2}, {3,4}, {0,1}, {2,3}};
    Vdouble costs = {4.7, 2.01, 3.2, 2};

    ListDigraph_VArc arcs;
    int nbArcs = arcsIndex.size();
    for(int k = 0; k < nbArcs; ++k){

        int i = arcsIndex.at(k).at(0);
        int j = arcsIndex.at(k).at(1);
        ListDigraph::Arc arc = addArc(graph, i, j);
        arcs.push_back(arc);

        superGraph->cost[arc] = costs.at(k);

    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);
    ListDigraph_VNode sourceNodes = {graph.nodeFromId(0)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(4)};

    VResourceBounds resources = VResourceBounds();

    int maxCardinality = 0;
    int costScaleFactor = 3;
    int capacityScaleFactor = 2;

    pNetwork network = std::make_shared<LemonNetwork>(superGraph, nodeFilter, arcFilter,
                                                      sourceNodes, Vdouble{1}, sinkNodes, Vdouble{-1},
                                                      maxCardinality, costScaleFactor, capacityScaleFactor);

    return std::make_tuple(network, superGraph, arcs, resources);

}

TEST(Network, BasicConstruction){

#ifdef GUROBI
    NetworkData temp = testGurobiNetwork1();
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;

    VVint arcCheck = {{1, 4}, {1,0}, {5, 0}, {0, 3}, {1, 0}, {1, 3}};

    int nbVertices = 0;
    for(ListDigraph::NodeIt vertexIt(graph); vertexIt != lemon::INVALID; ++vertexIt){
        ++nbVertices;
    }

    int nbArcs = 0;
    for(ListDigraph::ArcIt arcIt(graph); arcIt != lemon::INVALID; ++arcIt){
        ++nbArcs;
        int i = graph.id(graph.source(arcIt));
        int j = graph.id(graph.target(arcIt));
        ASSERT_NE(std::find(arcCheck.begin(), arcCheck.end(), Vint{i, j}), arcCheck.end());
    }

    ASSERT_EQ(nbVertices, 6);
    ASSERT_EQ(nbArcs, 6);
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, MultiNetworkConstruction){

#ifdef GUROBI
    double eps = 0.00000001;
    bool requireTightBound = false;

    MultiNetworkData temp = testGurobiMultiNetwork1(0);
    VpNetwork networks = std::get<0>(temp);
    const pNetwork& network1 = networks.at(0);
    const pNetwork& network2 = networks.at(1);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);


    //Original network

    ASSERT_TRUE(network1->solveModel());
    pListDigraph_ArcMap_double arcFlows1 = network1->getArcFlowValues();

    ASSERT_TRUE(network2->solveModel());
    pListDigraph_ArcMap_double arcFlows2 = network2->getArcFlowValues();

    const ListDigraph::Arc& a0_0 = arcs.at(0);
    const ListDigraph::Arc& a0_1 = arcs.at(1);
    const ListDigraph::Arc& a1_0 = arcs.at(2);
    const ListDigraph::Arc& a1_1 = arcs.at(3);

    EXPECT_NEAR((*arcFlows1)[a0_0], 0, eps);
    EXPECT_NEAR((*arcFlows1)[a0_1], 1, eps);
    EXPECT_NEAR((*arcFlows1)[a1_0], 0, eps);
    EXPECT_NEAR((*arcFlows1)[a1_1], 0, eps);
    EXPECT_NEAR(network1->getObjValue(), 1, eps);

    EXPECT_NEAR((*arcFlows2)[a0_0], 0, eps);
    EXPECT_NEAR((*arcFlows2)[a0_1], 0, eps);
    EXPECT_NEAR((*arcFlows2)[a1_0], 0, eps);
    EXPECT_NEAR((*arcFlows2)[a1_1], 1, eps);
    EXPECT_NEAR(network2->getObjValue(), 1, eps);


    //Different bounds

    resourceBounds.at(1) = ResourceBounds{0, 0.5};
    resourceBounds.at(2) = ResourceBounds{0.6, 1};
    network1->updateResourceBounds(resourceBounds, requireTightBound);
    network2->updateResourceBounds(resourceBounds, requireTightBound);

    ASSERT_TRUE(network1->solveModel());
    arcFlows1 = network1->getArcFlowValues();

    ASSERT_TRUE(network2->solveModel());
    arcFlows2 = network2->getArcFlowValues();

    EXPECT_NEAR((*arcFlows1)[a0_0], 0.5, eps);
    EXPECT_NEAR((*arcFlows1)[a0_1], 0.5, eps);
    EXPECT_NEAR((*arcFlows1)[a1_0], 0, eps);
    EXPECT_NEAR((*arcFlows1)[a1_1], 0, eps);
    EXPECT_NEAR(network1->getObjValue(), 5.5, eps);

    EXPECT_NEAR((*arcFlows2)[a0_0], 0, eps);
    EXPECT_NEAR((*arcFlows2)[a0_1], 0, eps);
    EXPECT_NEAR((*arcFlows2)[a1_0], 0.6, eps);
    EXPECT_NEAR((*arcFlows2)[a1_1], 0.4, eps);
    EXPECT_NEAR(network2->getObjValue(), 3.4, eps);


    //Non-supported bound changes

    resourceBounds.at(1) = ResourceBounds{1, std::numeric_limits<double>::infinity()}; //was finite
    resourceBounds.at(3) = ResourceBounds{0, 1}; //was inf
    EXPECT_ANY_THROW(network1->updateResourceBounds(resourceBounds, requireTightBound));
    EXPECT_ANY_THROW(network2->updateResourceBounds(resourceBounds, requireTightBound));
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, ResourceConstraintDuals){

#ifdef GUROBI
    double eps = 0.00000001;
    bool requireTightBound = false;

    MultiNetworkData temp = testGurobiMultiNetwork1(0);
    VpNetwork networks = std::get<0>(temp);
    const pNetwork& network1 = networks.at(0);
    const pNetwork& network2 = networks.at(1);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    //stabilize the dual
    for(int r = 0; r < 3; ++r){
        resourceBounds.at(r) = ResourceBounds{0.2, 0.9};
    }
    resourceBounds.at(3) = ResourceBounds{-std::numeric_limits<double>::infinity(), 0.9};
    network1->updateResourceBounds(resourceBounds, requireTightBound);
    network2->updateResourceBounds(resourceBounds, requireTightBound);

    network1->writeModelToFile("TEST_GurobiNetwork_ResourceConstraintDuals_network1.lp");
    network2->writeModelToFile("TEST_GurobiNetwork_ResourceConstraintDuals_network2.lp");

    ASSERT_TRUE(network1->solveModel());
    pListDigraph_ArcMap_double arcFlows1 = network1->getArcFlowValues();

    ASSERT_TRUE(network2->solveModel());
    pListDigraph_ArcMap_double arcFlows2 = network2->getArcFlowValues();

    typedef std::vector<std::map<constraint_sense , Vdouble> > dualContainer;
    dualContainer resourceDuals = dualContainer(3); //resourceDuals[network][sense] = dual values

    for(int i = 0; i < 2; ++i){
        for(constraint_sense sense : std::vector<constraint_sense>{LE, GE, EQ}){
            resourceDuals.at(i)[sense] = networks.at(i)->getResourceConstraintDuals(sense);
        }
    }

    dualContainer resourceDualsCompare = dualContainer(3);
    resourceDualsCompare.at(0)[LE] = Vdouble{0, 0};
    resourceDualsCompare.at(0)[GE] = Vdouble{9, 0};
    resourceDualsCompare.at(0)[EQ] = Vdouble{0, 0};
    resourceDualsCompare.at(1)[LE] = Vdouble{0, 0};
    resourceDualsCompare.at(1)[GE] = Vdouble{4, 0};
    resourceDualsCompare.at(1)[EQ] = Vdouble{0, 0};

    for(int i = 0; i < 2; ++i){
        for(constraint_sense sense : std::vector<constraint_sense>{LE, GE, EQ}){
            for(int j = 0; j < 2; ++j){

                double val = resourceDuals.at(i).at(sense).at(j);
                double compare = resourceDualsCompare.at(i).at(sense).at(j);

                EXPECT_TRUE( (!std::isnan(val) && !std::isnan(compare)) || (std::isnan(val) && std::isnan(compare)) );

                if(!std::isnan(val)) {
                    EXPECT_NEAR(resourceDuals.at(i).at(sense).at(j), resourceDualsCompare.at(i).at(sense).at(j), eps);
                }

            }
        }
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, MinCostFlowWithSideConstraint){

#ifdef GUROBI
    double eps = 0.00000001;
    bool requireTightBound = false;

    NetworkData temp = testGurobiNetwork2();
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->writeModelToFile("TEST_GurobiNetwork_MinCostFlowWithSideConstraint_1.lp");

    ASSERT_TRUE(network->solveModel());

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{10, 8, 5, 2, 3, 5, 7};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    resourceBounds.at(7) = {32.5, 32.5};
    network->updateResourceBounds(resourceBounds, requireTightBound);
    network->writeModelToFile("TEST_GurobiNetwork_MinCostFlowWithSideConstraint_2.lp");

    ASSERT_TRUE(network->solveModel());

    arcFlowValues = network->getArcFlowValues();
    compareValues = Vdouble{9.9375, 8, 4.9375, 1.9375, 3, 5.0625, 7};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    resourceBounds.at(8).at(0) = 5.5;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    network->writeModelToFile("TEST_GurobiNetwork_MinCostFlowWithSideConstraint_3.lp");

    ASSERT_TRUE(network->solveModel());

    arcFlowValues = network->getArcFlowValues();
    compareValues = Vdouble{9.5, 7.125, 4.5, 2.375, 2.125, 5.5, 7.875};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(Kevintools, addExpressionBoundConstraint){

#ifdef GUROBI
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    GRBLinExpr expression = GRBLinExpr();
    GRBVar x = model.addVar(0, 1, 0, GRB_CONTINUOUS, "");
    expression += x;

    VResourceBounds noThrowBounds = VResourceBounds();
    noThrowBounds.push_back(ResourceBounds{-std::numeric_limits<double>::infinity(), -9});
    noThrowBounds.push_back(ResourceBounds{15, std::numeric_limits<double>::infinity()});
    noThrowBounds.push_back(ResourceBounds{15, 15});
    noThrowBounds.push_back(ResourceBounds{15, 18});
    for(ResourceBounds bounds : noThrowBounds){
        EXPECT_NO_THROW(kevintools::addExpressionBoundConstraint(model, expression, bounds));
    }

    VResourceBounds throwBounds = VResourceBounds();
    throwBounds.push_back(ResourceBounds{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()});
    throwBounds.push_back(ResourceBounds{5, -std::numeric_limits<double>::infinity()});
    throwBounds.push_back(ResourceBounds{std::numeric_limits<double>::infinity(), 10});
    throwBounds.push_back(ResourceBounds{11, 10});
    for(ResourceBounds bounds : throwBounds){
        EXPECT_ANY_THROW(kevintools::addExpressionBoundConstraint(model, expression, bounds));
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, testCardinalityZero){

#ifdef GUROBI
    double eps = 0.00000001;
    int maxCardinality = 0;

    NetworkData temp = testGurobiNetwork3(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->writeModelToFile("TEST_testGurobiNetwork3_cardinality0.lp");
    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0.9, 0, 1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, testCardinalityTwo){

#ifdef GUROBI
    double eps = 0.00000001;
    int maxCardinality = 2;

    NetworkData temp = testGurobiNetwork3(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->writeModelToFile("TEST_testGurobiNetwork3_cardinality2.lp");
    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0, 0.9, 0.1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    Vdouble duals = Vdouble(7);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    EXPECT_NEAR(duals.at(0), -1, eps);
    for(int i = 1; i < duals.size(); ++i){
        EXPECT_NEAR(duals.at(i), 0, eps);
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, testCardinalityThree){

#ifdef GUROBI
    double eps = 0.00000001;
    int maxCardinality = 3;

    NetworkData temp = testGurobiNetwork3(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->writeModelToFile("TEST_testGurobiNetwork3_cardinality3.lp");
    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0.9, 0, 1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    Vdouble duals = Vdouble(7);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    EXPECT_NEAR(duals.at(0), -3, eps);
    for(int i = 1; i < duals.size(); ++i){
        EXPECT_NEAR(duals.at(i), 0, eps);
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

TEST(GurobiNetwork, testCardinalityTen){

#ifdef GUROBI
    double eps = 0.00000001;
    int maxCardinality = 10;

    NetworkData temp = testGurobiNetwork3(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->writeModelToFile("TEST_testGurobiNetwork3_cardinality10.lp");
    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0.9, 0, 1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    Vdouble duals = Vdouble(7);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    EXPECT_NEAR(duals.at(0), -3, eps);
    for(int i = 1; i < duals.size(); ++i){
        EXPECT_NEAR(duals.at(i), 0, eps);
    }
#else
    std::cout << "Test skipped because Gurobi is missing." << std::endl;
    ASSERT_TRUE(false);
#endif

}

Vdouble solveComparePrimalsAndGetLEduals(const ListDigraph& graph,
                                         const Vdouble& compareValues, const double eps,
                                         const pNetwork& network, const ListDigraph_VArc& arcs,
                                         const bool checkArcFlowPath=false) {

    std::cout << "Testing compareValues = " << kevintools::to_string(compareValues) << "." << std::endl;
    EXPECT_TRUE(network->solveModel());

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    if(checkArcFlowPath) {

        const ListDigraph_VArc& arcFlowPath = network->getArcFlowPath();

        int nbPositiveFlows = 0;
        for (int i = 0; i < arcs.size(); ++i) {
            if ((*arcFlowValues)[arcs.at(i)] > 0.5) {
                ++nbPositiveFlows;
            }
        }
        EXPECT_EQ(nbPositiveFlows, arcFlowPath.size());

        for(int i = 0; i < arcFlowPath.size(); ++i){
            EXPECT_TRUE((*arcFlowValues)[arcFlowPath.at(i)] > 0.5);
        }

        for(int i = 0; i + 1 < arcFlowPath.size(); ++i){
            ListDigraph::Node thisTarget = graph.target(arcFlowPath.at(i));
            ListDigraph::Node nextSource = graph.source(arcFlowPath.at(i+1));
            EXPECT_EQ(thisTarget, nextSource);
        }

    }

    Vdouble duals = Vdouble(6);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();
    EXPECT_TRUE(networkDuals.size() == resourceIndices.size());

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    return duals;

}

Vdouble solveAndGetLEduals(const double eps, const pNetwork& network, const ListDigraph_VArc& arcs) {

    network->solveModel();

    Vdouble duals = Vdouble(6);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();
    EXPECT_TRUE(networkDuals.size() == resourceIndices.size());

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    return duals;

}

TEST(LemonNetwork, rounding){

    double eps = 0.00000001;
    bool requireTightBound = false;

    int costScaleFactor = 10;
    int capacityScaleFactor = 2;

    pNetwork _network;
    pLemonNetwork network;
    pSuperGraph superGraph;
    ListDigraph_VArc arcs;
    VResourceBounds resourceBounds;

    std::tie(_network, superGraph, arcs, resourceBounds) = testLemonNetwork1(costScaleFactor, capacityScaleFactor);
    network = std::dynamic_pointer_cast<LemonNetwork>(_network);
    const ListDigraph& graph = superGraph->graph;

    Vdouble zero = Vdouble(resourceBounds.size(), 0);

    //Test for different values of z = resourceBounds.at(0).at(1)
    //First, go over breakpoints
    //Check objectives if no rounding is necessary

    resourceBounds.at(0).at(1) = 0;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z00_duals = solveComparePrimalsAndGetLEduals(graph, {0, 0, 0, 0, 0, 2}, eps, network, arcs);
    EXPECT_NEAR(network->getObjValue(), 20, eps);
    EXPECT_TRUE(kevintools::almostLessOrEqual(z00_duals, zero, eps));
    EXPECT_FALSE(kevintools::almostEqual(z00_duals, zero, eps));

    resourceBounds.at(0).at(1) = 0.5;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z05_duals = solveComparePrimalsAndGetLEduals(graph, {0.5, 0, 0, 0.5, 0.5, 1.5}, eps, network, arcs);
    EXPECT_NEAR(network->getObjValue(), 1.5*10 + 0.5*0.9, eps);
    EXPECT_TRUE(kevintools::almostLessOrEqual(z05_duals, zero, eps));
    EXPECT_FALSE(kevintools::almostEqual(z05_duals, zero, eps));

    resourceBounds.at(0).at(1) = 1;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z10_duals = solveComparePrimalsAndGetLEduals(graph, {1, 0, 0.5, 0.5, 1, 1}, eps, network, arcs);
    //No objective checking, as arcs with rounded capacity are used
    EXPECT_TRUE(kevintools::almostLessOrEqual(z10_duals, zero, eps));
    EXPECT_FALSE(kevintools::almostEqual(z10_duals, zero, eps));

    resourceBounds.at(0).at(1) = 1.5;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z15_duals = solveComparePrimalsAndGetLEduals(graph, {1.5, 0.5, 0.5, 0.5, 1.5, 0.5}, eps, network, arcs);
    //No objective checking, as arcs with rounded capacity are used
    EXPECT_TRUE(kevintools::almostLessOrEqual(z15_duals, zero, eps));
    //Duals may be all zero

    resourceBounds.at(0).at(1) = 2;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z20_duals = solveComparePrimalsAndGetLEduals(graph, {1.5, 0.5, 0.5, 0.5, 1.5, 0.5}, eps, network, arcs);
    //No objective checking, as arcs with rounded capacity are used
    EXPECT_TRUE(kevintools::almostLessOrEqual(z20_duals, zero, eps));
    //Duals may be all zero

    //Check if objective is correctly calculated from the duals

    //Value 0.3 is rounded to 0.5
    resourceBounds.at(0).at(1) = 0.3;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z05_duals_check = solveAndGetLEduals(eps, network, arcs);
    EXPECT_EQ(z05_duals, z05_duals_check); //if fails, could be another optimal dual solution
    double z03_obj_check = 15.45 +
            (0.3 - 0.5) * z05_duals.at(0) +
            (0.33 - 0.5) * z05_duals.at(1) +
            (0.6 - 0.5) * z05_duals.at(2);
    EXPECT_NEAR(network->getObjValue(), z03_obj_check, eps);
    EXPECT_TRUE(z03_obj_check <= 1.7*10 + 0.3*0.9 + eps);
    std::cout << "getObjValue() = " << network->getObjValue() << ", z03_obj_check = " << z03_obj_check << " <= " << 1.7*10 + 0.3*0.9 << std::endl;

    //Value 0.62 is rounded to 0.5
    resourceBounds.at(0).at(1) = 0.62;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    z05_duals_check = solveAndGetLEduals(eps, network, arcs);
    EXPECT_EQ(z05_duals, z05_duals_check); //if fails, could be another optimal dual solution
    double z062_obj_check = 15.45 +
                           (0.62 - 0.5) * z05_duals.at(0) +
                           (0.33 - 0.5) * z05_duals.at(1) +
                           (0.6 - 0.5) * z05_duals.at(2);
    EXPECT_NEAR(network->getObjValue(), z062_obj_check, eps);
    EXPECT_TRUE(z062_obj_check <= 1.38*10 + 0.5*0.9 + 0.12*2.6 + eps);
    std::cout << "getObjValue() = " << network->getObjValue() << ", z062_obj_check = " << z062_obj_check << " <= " << 1.38*10 + 0.5*0.9 + 0.12*2.6 << std::endl;

    std::cout << "z05_duals = " << kevintools::to_string(z05_duals) << "." << std::endl;

    //Value 1.4 is rounded to 1.5
    resourceBounds.at(0).at(1) = 1.4;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z15_duals_check = solveAndGetLEduals(eps, network, arcs);
    EXPECT_EQ(z15_duals, z15_duals_check); //if fails, could be another optimal dual solution
    double z14_obj_check = 8.6 +
                           (1.4 - 1.5) * z15_duals.at(0) +
                           (0.33 - 0.5) * z15_duals.at(1) +
                           (0.6 - 0.5) * z15_duals.at(2);
    EXPECT_NEAR(network->getObjValue(), z14_obj_check, eps);
    EXPECT_TRUE(z14_obj_check <= 0.6*10 + 0.5*0.9 + 0.6*2.6 + 0.3*3.7 + eps);
    std::cout << "getObjValue() = " << network->getObjValue() << ", z14_obj_check = " << z14_obj_check << " <= " << 0.6*10 + 0.5*0.9 + 0.6*2.6 + 0.3*3.7 << std::endl;

    //Value 1.7 is rounded to 1.5
    resourceBounds.at(0).at(1) = 1.7;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    z15_duals_check = solveAndGetLEduals(eps, network, arcs);
    EXPECT_EQ(z15_duals, z15_duals_check); //if fails, could be another optimal dual solution
    double z17_obj_check = 8.6 +
                           (1.7 - 1.5) * z15_duals.at(0) +
                           (0.33 - 0.5) * z15_duals.at(1) +
                           (0.6 - 0.5) * z15_duals.at(2);
    EXPECT_NEAR(network->getObjValue(), z17_obj_check, eps);
    EXPECT_TRUE(z17_obj_check <= 0.57*10 + 0.5*0.9 + 0.6*2.6 + 0.33*3.7 + eps);
    std::cout << "getObjValue() = " << network->getObjValue() << ", z17_obj_check = " << z17_obj_check << " <= " << 0.57*10 + 0.5*0.9 + 0.6*2.6 + 0.33*3.7 << std::endl;

    std::cout << "z15_duals = " << kevintools::to_string(z15_duals) << "." << std::endl;

}

TEST(LemonNetwork, testCardinalityZero){

    double eps = 0.00000001;
    int maxCardinality = 0;

    NetworkData temp = testLemonNetwork2(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0.9, 0, 1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

}

TEST(LemonNetwork, testCardinalityTwo){

    double eps = 0.00000001;
    int maxCardinality = 2;

    NetworkData temp = testLemonNetwork2(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0, 0.9, 0.1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    Vdouble duals = Vdouble(7);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    EXPECT_NEAR(duals.at(0), -1, eps);
    for(int i = 1; i < duals.size(); ++i){
        EXPECT_NEAR(duals.at(i), 0, eps);
    }

}

TEST(LemonNetwork, testCardinalityThree){

    double eps = 0.00000001;
    int maxCardinality = 3;

    NetworkData temp = testLemonNetwork2(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0.9, 0, 1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    Vdouble duals = Vdouble(7);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    EXPECT_NEAR(duals.at(0), -3, eps);
    for(int i = 1; i < duals.size(); ++i){
        EXPECT_NEAR(duals.at(i), 0, eps);
    }

}

TEST(LemonNetwork, testCardinalityTen){

    double eps = 0.00000001;
    int maxCardinality = 10;

    NetworkData temp = testLemonNetwork2(maxCardinality);
    pNetwork network = std::get<0>(temp);
    pSuperGraph superGraph = std::get<1>(temp);
    const ListDigraph& graph = superGraph->graph;
    ListDigraph_VArc arcs = std::get<2>(temp);
    VResourceBounds resourceBounds = std::get<3>(temp);

    network->solveModel();

    pListDigraph_ArcMap_double arcFlowValues = network->getArcFlowValues();
    Vdouble compareValues = Vdouble{0.9, 0.1, 0, 0, 0.9, 0, 1};

    for(int i = 0; i < arcs.size(); ++i){
        EXPECT_NEAR((*arcFlowValues)[arcs.at(i)], compareValues.at(i), eps);
    }

    Vdouble duals = Vdouble(7);
    const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
    const Vint& resourceIndices = network->getResourceIndices();

    for(int i = 0; i < resourceIndices.size(); ++i) {
        const int& resource = resourceIndices.at(i);
        duals.at(resource) += networkDuals.at(i);
    }

    EXPECT_NEAR(duals.at(0), -3, eps);
    for(int i = 1; i < duals.size(); ++i){
        EXPECT_NEAR(duals.at(i), 0, eps);
    }

}

TEST(LemonNetwork, integerData){

    double eps = 0.00000001;

    pNetwork _network;
    pLemonNetwork network;
    pSuperGraph superGraph;
    ListDigraph_VArc arcs;
    VResourceBounds resourceBounds;

    std::tie(_network, superGraph, arcs, resourceBounds) = testLemonNetwork3();
    network = std::dynamic_pointer_cast<LemonNetwork>(_network);
    const ListDigraph& graph = superGraph->graph;

    bool requireTightBound = true;
    bool checkArcFlowPath = true;

    //All arcs open
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble duals = solveComparePrimalsAndGetLEduals(graph, {0, 0, 0, 1, 1, 0}, eps, network, arcs, checkArcFlowPath);
    EXPECT_NEAR(network->getObjValue(), 3.3, eps);
    std::cout << "duals = " << kevintools::to_string(duals) << "." << std::endl;
    EXPECT_TRUE(kevintools::almostEqual(duals, Vdouble{0, 0, 0, 0, 0, 0}, eps));

    //Close arc {0,2}
    resourceBounds.at(3).at(1) = 0;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    duals = solveComparePrimalsAndGetLEduals(graph, {1, 0, 0, 0, 0, 0}, eps, network, arcs, checkArcFlowPath);
    EXPECT_NEAR(network->getObjValue(), 4.7, eps);
    std::cout << "duals = " << kevintools::to_string(duals) << "." << std::endl;
    EXPECT_TRUE(kevintools::almostEqual(duals, Vdouble{0, 0, 0, -2.7, 0, 0}, eps));

    //Additionally close arc {2,3}
    resourceBounds.at(4).at(1) = 0;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    duals = solveComparePrimalsAndGetLEduals(graph, {1, 0, 0, 0, 0, 0}, eps, network, arcs, checkArcFlowPath);
    EXPECT_NEAR(network->getObjValue(), 4.7, eps);
    std::cout << "duals = " << kevintools::to_string(duals) << "." << std::endl;
    EXPECT_TRUE(kevintools::almostEqual(duals, Vdouble{0, 0, 0, -2.7, 0, 0}, eps));

    //Reopen arc {0,2}, keep {2,3} closed
    resourceBounds.at(3).at(1) = 1;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    duals = solveComparePrimalsAndGetLEduals(graph, {1, 0, 0, 0, 0, 0}, eps, network, arcs, checkArcFlowPath);
    EXPECT_NEAR(network->getObjValue(), 4.7, eps);
    std::cout << "duals = " << kevintools::to_string(duals) << "." << std::endl;
    EXPECT_TRUE(kevintools::almostEqual(duals, Vdouble{0, 0, 0, 0, -1.4, 0}, eps));

}

TEST(LemonNetwork, largeScaleFactors) {

    double eps = 0.0000001;
    bool requireTightBound = false;

    int costScaleFactor = std::numeric_limits<int>::max()/20.16;
    int capacityScaleFactor = std::numeric_limits<int>::max()/2;

    pNetwork _network;
    pLemonNetwork network;
    pSuperGraph superGraph;
    ListDigraph_VArc arcs;
    VResourceBounds resourceBounds;

    std::tie(_network, superGraph, arcs, resourceBounds) = testLemonNetwork1(costScaleFactor, capacityScaleFactor);
    network = std::dynamic_pointer_cast<LemonNetwork>(_network);
    const ListDigraph& graph = superGraph->graph;

    Vdouble zero = Vdouble(resourceBounds.size(), 0);

    resourceBounds.at(0).at(1) = 0.5;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    Vdouble z05_duals = solveComparePrimalsAndGetLEduals(graph, {0.5, 0, 0, 0.5, 0.5, 1.5}, eps, network, arcs);
    EXPECT_NEAR(network->getObjValue(), 1.5*10.08 + 0.5*0.99, eps);
    EXPECT_TRUE(kevintools::almostLessOrEqual(z05_duals, zero, eps));
    EXPECT_FALSE(kevintools::almostEqual(z05_duals, zero, eps));

}

TEST(LemonNetwork, tooLargeScaleFactors) {

    double eps = 0.00000001;
    bool requireTightBound = false;

    int costScaleFactor = 1E9;
    int capacityScaleFactor = 1E9;

    pNetwork _network;
    pLemonNetwork network;
    pSuperGraph superGraph;
    ListDigraph_VArc arcs;
    VResourceBounds resourceBounds;

    std::tie(_network, superGraph, arcs, resourceBounds) = testLemonNetwork1(costScaleFactor, capacityScaleFactor);
    network = std::dynamic_pointer_cast<LemonNetwork>(_network);
    const ListDigraph& graph = superGraph->graph;

    Vdouble zero = Vdouble(resourceBounds.size(), 0);

    resourceBounds.at(0).at(1) = 0.5;
    network->updateResourceBounds(resourceBounds, requireTightBound);
    EXPECT_FALSE(network->solveModel());

}

TEST(LemonNetwork, arcFlowPath){

    double eps = 0.00000001;

    pNetwork _network;
    pLemonNetwork network;
    pSuperGraph superGraph;
    ListDigraph_VArc arcs;
    VResourceBounds resourceBounds;

    std::tie(_network, superGraph, arcs, resourceBounds) = testLemonNetwork4();
    network = std::dynamic_pointer_cast<LemonNetwork>(_network);
    const ListDigraph& graph = superGraph->graph;

    bool requireTightBound = true;
    bool checkArcFlowPath = true;

    network->updateResourceBounds(resourceBounds, requireTightBound);
    solveComparePrimalsAndGetLEduals(graph, {1, 1, 1, 1}, eps, network, arcs, checkArcFlowPath);

}