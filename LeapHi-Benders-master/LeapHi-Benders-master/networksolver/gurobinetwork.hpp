#pragma once

#include "network.hpp"

namespace networksolver {

    class GurobiNetwork : public Network{


    //Building the network
    public:
        GurobiNetwork(const pSuperGraph& superGraph,
                      const ListDigraph::NodeMap<bool>& nodeFilter,
                      const ListDigraph::ArcMap<bool>& arcFilter,
                      const ListDigraph_VNode& sourceNode,
                      const Vdouble& initialSourceFlows,
                      const ListDigraph_VNode& sinkNodes,
                      const Vdouble& initialSinkFlows,
                      int maxCardinality,
                      const pGRBEnv& env,
                      const VResourceBounds& initialResourceBounds);
        void initializeResourceConstraints(const VResourceBounds& externalResourceBounds);


    //Update the network
    public:
        void updateResourceBounds(const VResourceBounds& resourceBounds, bool requireTightBound) override;
            //requireTightBound is ignored, as the bounds given by GurobiNetwork are always tight
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


    //Gurobi
    private:
        const pGRBEnv env;      //env is shared between subProblems
        const pGRBModel model;  //every network has its own model
        StaticDigraph::ArcMap<GRBVar> arcFlows;               //maps layeredGraph arcs to flow variables
        StaticDigraph::NodeMap<GRBConstr> flowConstraints;    //maps layerGraph nodes to flow constraints
        void initializeVariablesAndObjective();
        void initializeFlowConstraints();
        VGRBLinExpr buildResourceExpressionsLayer(int layer) const;
        VVSenseConstraintMap_GRB resourceConstraints;   //resourceConstraints[layer] := resource constraints that belong to 'layer'
                                                        //resourceConstraints[layer] has the same index as resourceIndices
        void updateResourceConstraints(const VResourceBounds& resourceBounds);

    };

}