#pragma once

#include "network.hpp"

namespace networksolver {

    class LemonNetwork : public Network {

    //Building the network
    public:
        LemonNetwork(const pSuperGraph& superGraph,
                     const ListDigraph::NodeMap<bool>& nodeFilter,
                     const ListDigraph::ArcMap<bool>& arcFilter,
                     const ListDigraph_VNode& sourceNode,
                     const Vdouble& initialSourceFlows,
                     const ListDigraph_VNode& sinkNodes,
                     const Vdouble& initialSinkFlows,
                     int maxCardinality,
                     int costScaleFactor,
                     int capacityScaleFactor);


    //Update the network
    public:
        void updateResourceBounds(const VResourceBounds& resourceBounds, bool requireTightBound) override; //call before solve
        void updateSourceFlows(const Vdouble& sourceFlows) override;
        void updateSinkFlows(const Vdouble& sinkFlows) override;


    //Solving and retrieving solutions
    public:
        bool solveModel() override;
        int getStatus() const override;
        double getObjValue() const override;
        void writeModelToFile(const std::string& filename) const override;
        pListDigraph_ArcMap_double getArcFlowValues() const override;                       //aggregated values over all associated layered arcs
        ListDigraph_VArc getArcFlowPath(const StaticDigraph::Node& source, const StaticDigraph::Node& sink) const override;
        Vdouble getResourceConstraintDuals(const constraint_sense& sense) const override;   //same size as resourceIndices which indicate the resource numbers
    private:
        double objValue; //set by solveModel()
        void setObjValue(); //requires setResourceConstraintDuals() to be called first
        Vdouble resourceConstraintDuals; //set by solveModel()
        void setResourceConstraintDuals();
        void setResourceConstraintDualsMCF();
        void setResourceConstraintDualsDijkstra();
    //Method getObjValue is guaranteed to give a lower bound on the objective value.
    //The combination of objValue and resourceConstraintDuals is guaranteed to give a valid cut for the value function.
    //Note: this cut need not be tight. If costScaleFactor and capacityScaleFactor increase, the cut will become better, and tight for limit->inf.
    //getArcFlowValues returns the corresponding primal solution, which is optimal for the scaled problem,
        //but not guaranteed to be optimal or even feasible for the original problem.
    //In case of a single source, single sink, flow 1, and binary or unbounded resource upper bounds:
        //Dijkstra is used, and the objective is guaranteed to be tight, and the primal solution is optimal.


    //Lemon
    private:
        const int costScaleFactor;      //scaling factor used in rounding costs to integers
        const int capacityScaleFactor;  //scaling factor used in rounding capacities to integers
        pStaticDigraph_ArcMap_int capacityMap;
        void updateCapacityAndCapacityCorrections(const networksolver::VResourceBounds& externalResourceBounds);
        void updateArcFilter(const networksolver::VResourceBounds& externalResourceBounds);
        pStaticDigraph_ArcMap_int roundedCostMap;
        pStaticDigraph_ArcMap_double unroundedCostMap;
        void initializeRoundedCostMap();
        void initializeUnroundedCostMap();
        pStaticDigraph_NodeMap_int supplyMap;
        bool useSupplyMap() const;
        void updateSupplyMap();
        pNetworkSimplex mcfAlgorithm;
        pFilterArcsDijkstra dijkstraAlgorithm;
        int status; //0 for successful, 1 for error
        Vdouble capacityCorrections; //same index as relevantSuperArcs
        void verifyFlowsIntegral() const;
        void verifySuperGraphAssumptions(); //Verify:
                                                    //1. Every resource is associated with at most one arc.
                                                    //2. Every arc is associated with at most resource.
                                                    //3. All resource contributions have value 1.
                                                    //4. All arc cost >= 0 (warning only).
        bool negativeArcCosts;
        bool useDijkstra;
        bool preventSolve;
        pStaticDigraph_ArcMap_bool arcFilter;
        pFilterArcs filteredGraph;
    public:
        ListDigraph_VVArc getAllShortestPathsShuttleFiltered(const intSet& superHubIDs);

    };

}