#pragma once

#include "networksolver_typedefs.hpp"

namespace networksolver {

    class SuperGraph{

    public:
        explicit SuperGraph(int nbVertices);

    public:
        ListDigraph graph;
        ListDigraph::ArcMap<double> cost;       //cost of traversing this (parallel) arc
        ListDigraph::ArcMap<Vint> properties;   //integer arc properties
        ListDigraph::ArcMap<VResourceContribution> contributions; //vector of contributions of (parallel) arc
                                                                  //contributions[l] := (index of resource, increase in value)

    public:
        ListDigraph::Arc getUniqueMatchingArc(
                const ListDigraph::Node& from, const ListDigraph::Node& to, const Vint& matchingProperties) const;
        ListDigraph_VArc getMatchingArcs(
                const ListDigraph::Node& from, const ListDigraph::Node& to, const Vint& matchingProperties, const Vbool& matchOnElements) const;

    //Create an Arc Lookup Table to speed up getUniqueMatchingArc() and getMatchingArcs()
    //When the graph is modified, the lookup table is invalidated without warning
    //It is therefore recommended to delete the lookup table before modifying the graph
    private:
        std::shared_ptr<lemon::AllArcLookUp<ListDigraph>> lookup;
    public:
        void newArcLookupTable();
        void deleteArcLookupTable();

    };

}