#pragma once

#include "networksolver_typedefs.hpp"

namespace networksolver {

    class Network {

        //Class for solving network flow problems with cardinality and resource constraints in continuous variables.
        //Cardinality constraints are modeled through a layered graph.
        //Resource constraints are handled as arc capacities if possible, and side constraints otherwise.
        //Note that each network is solved independently, i.e., there is no interaction between Networks derived from the same superGraph.

        //The parameter maxCardinality is the maximum cardinality of any (fractional) path through the network.
        //If maxCardinality = 0, there is no restriction: only layer 0 is used, nbLayers() = 1 and useLayers() = false.

        //If layers are used (maxCardinality >= 2), then a cardinality constraint is enforced through layers.
        //e.g. if maxCardinality = 3, then at most three arcs are used to get from sink to source
        //In this case, nbLayers() >= 1 and useLayers() = true.
        //The network graph is represented similar to https://arxiv.org/abs/1912.02308v1, i.e., nbLayers = maxCardinality + 1.
        //If the source(s) have in-arcs or the sink(s) have out-arcs, these nodes are duplicated as necessary.

        //To look up information about the layeredGraph arcs, use layeredToSuperArcs to map to superGraph arcs,
        //which stores resource contributions, costs, etc.

        //IMPORTANT ASSUMPTIONS in case layers are used
        //1) It is assumed that it is suboptimal to revisit arcs! (not checked)
        //      If arcs can be revisited, the joint contribution to the original arc can violate the resource bounds.
        //2) It is assumed that resource constraints can be enforced per layer independently! (not checked)
        //      This is the case, for example, when each resource belongs to a single arc, and 1) holds.
        //      Then a constraint can be added in each layer separately.
        //      Note that side constraints may violate this assumption.


    //Building the network
    public:
        Network(const pSuperGraph& superGraph,
                const ListDigraph::NodeMap<bool>& nodeFilter,
                const ListDigraph::ArcMap<bool>& arcFilter,
                const ListDigraph_VNode& sourceNodes,
                const Vdouble& initialSourceFlows,
                const ListDigraph_VNode& sinkNodes,
                const Vdouble& initialSinkFlows,
                int maxCardinality);

    //Update the network
    public:
        virtual void updateResourceBounds(const VResourceBounds& resourceBounds, bool requireTightBound) =0;
        virtual void updateSourceFlows(const Vdouble& sourceFlows);
        virtual void updateSinkFlows(const Vdouble& sinkFlows);


    //Solving and retrieving solutions
    public:
        virtual bool solveModel() =0;
        virtual int getStatus() const =0;
        virtual double getObjValue() const =0;
        virtual void writeModelToFile(const std::string& filename) const =0;
        virtual pListDigraph_ArcMap_double getArcFlowValues() const =0;                     //aggregated values over all associated layered arcs
        ListDigraph_VArc getArcFlowPath() const;
    protected:
        virtual ListDigraph_VArc getArcFlowPath(const StaticDigraph::Node& source, const StaticDigraph::Node& sink) const;
    public:
        virtual Vdouble getResourceConstraintDuals(const constraint_sense& sense) const =0; //same size as resourceIndices which indicate the resource numbers


    //Layered network
    protected:
        StaticDigraph layeredGraph;
        void buildLayeredNetwork(const ListDigraph::NodeMap<bool>& nodeFilter,
                                 const ListDigraph::ArcMap<bool>& arcFilter,
                                 const ListDigraph_VNode& sourceNodes,
                                 const ListDigraph_VNode& sinkNodes);
        StaticDigraph::ArcMap<int> layeredToSuperArcIds;
        const ListDigraph::Arc layeredToSuperArc(const StaticDigraph::Arc& layeredArc) const;
                                                        //map arcs to superGraph
        StaticDigraph::ArcMap<int> originLayer;         //map arcs to origin layer
        StaticDigraph::ArcMap<int> destinationLayer;    //map arcs to destination layer
        Vint relevantSuperArcIds;
        const ListDigraph::Arc relevantSuperArc(int index) const;
        StaticDigraph::ArcMap<int> layeredToRelevantSuperArcIndex;
        void setRelevantSuperArcs();
        bool useLayers() const;
        int nbLayers() const;


    //Layered network - build nodes
    //Add nodes to an empty graph to build a layered copy, taking the filters into account
    private:
        typedef ListDigraph::NodeMap<ListDigraph_VNode> NodeToVNodeMap;
        std::shared_ptr<NodeToVNodeMap> addNodesToFilteredLayeredCopy(const ListDigraph& original, ListDigraph& copy,
                                                                      const ListDigraph::NodeMap<bool>& nodeFilter,
                                                                      const ListDigraph::ArcMap<bool>& arcFilter,
                                                                      const ListDigraph_VNode& sourceNodes,
                                                                      const ListDigraph_VNode& sinkNodes) const;
            //Output maps original nodes to nbLayers copy nodes, corresponding to different layers
            //(*pOriginalToCopyNodes)[node].at(i) = lemon::INVALID if no corresponding node in layer i
        bool hasInArcInFilter(const ListDigraph::Node& node, const ListDigraph& graph, const ListDigraph::ArcMap<bool>& arcFilter) const;
        bool hasOutArcInFilter(const ListDigraph::Node& node, const ListDigraph& graph, const ListDigraph::ArcMap<bool>& arcFilter) const;


    //Layered network - build arcs
    //Add arcs to a graph built by addNodesToFilteredLayeredCopy() to build a layered copy, taking the filters into account
    private:
        typedef ListDigraph::ArcMap<int> ArcToIntMap;
        std::shared_ptr<ArcToIntMap> addArcsToFilteredLayeredCopy(const ListDigraph& original, ListDigraph& copy,
                                                                  const std::shared_ptr<NodeToVNodeMap>& pOriginalToCopyNodes,
                                                                  const ListDigraph::NodeMap<bool>& nodeFilter,
                                                                  const ListDigraph::ArcMap<bool>& arcFilter) const;
            //Output maps copy arcs to original arc ids, different layers map to the same original arc


    //Layered network - get origin layer
    private:
        std::shared_ptr<ArcToIntMap> getCopyOriginLayerMap(const ListDigraph& original, ListDigraph& copy,
                                                           const std::shared_ptr<NodeToVNodeMap>& pOriginalToCopyNodes) const;
            //Output maps copy arc to origin layer
        std::shared_ptr<ArcToIntMap> getCopyDestinationLayerMap(const ListDigraph& original, ListDigraph& copy,
                                                           const std::shared_ptr<NodeToVNodeMap>& pOriginalToCopyNodes) const;
            //Output maps copy arc to destination layer


    //Resources
    protected:
        Vint resourceIndices;   //Indices of relevant resources that are defined externally (e.g., in NetworkDesignProblem).
                                    //Different networks can use the same resources (e.g., when the rhs is the same)
                                    //but the resources are considered to be different on the network level.
        Vint externalToLocalResourceMap;
        void setExternalToLocalResourceMap();
        void clearExternalToLocalResourceMap();
    protected:
        void setResourceIndicesFromNetworkArcs();
        VResourceBounds externalToLocalResourceBounds(const networksolver::VResourceBounds& externalResourceBounds) const;
            //Select subset from externalResourceBounds based on resourceIndices
    public:
        const Vint& getResourceIndices() const;


    //Sinks and sources
    protected:
        bool isLayeredSource(const StaticDigraph::Node& node) const;
        bool isLayeredSink(const StaticDigraph::Node& node) const;
        double getLayeredSourceFlow(const StaticDigraph::Node& node) const;
        double getLayeredSinkFlow(const StaticDigraph::Node& node) const;
        StaticDigraph_VNode layeredSources; //sources in layeredGraph
        StaticDigraph_VNode layeredSinks;   //sinks in layeredGraph
        Vdouble sourceFlows;                //same index as layeredSources
        Vdouble sinkFlows;                  //same index as layeredSinks
    public:
        double getTotalSourceFlow() const;
        double getTotalSinkFlow() const;


    protected:
        const pSuperGraph superGraph;
        const int maxCardinality; //maximum cardinality of any path through the network, -1 if unbounded

    };

}