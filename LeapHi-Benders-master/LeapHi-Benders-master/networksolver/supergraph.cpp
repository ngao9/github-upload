#include "supergraph.hpp"

using namespace networksolver;

SuperGraph::SuperGraph(const int nbVertices) :
    graph(), cost(graph), properties(graph), contributions(graph), lookup(nullptr){

    std::vector<ListDigraph::Node> allVertices;
    allVertices.reserve(nbVertices);

    for(int i = 0; i < nbVertices; ++i){
        allVertices.push_back(graph.addNode());
    }

    for(int i = 0; i < nbVertices; ++i){
        if(graph.id(allVertices.at(i)) != i){
            throw std::runtime_error("Exception: id does not match index.");
        }
    }

}

ListDigraph::Arc SuperGraph::getUniqueMatchingArc(const ListDigraph::Node& from, const ListDigraph::Node& to, const Vint& matchingProperties) const {

    int nbMatches = 0;
    ListDigraph::Arc result = lemon::INVALID;

    auto processArc = [&](const ListDigraph::Arc& arc){

        if(properties[arc].size() != matchingProperties.size()){
            throw std::runtime_error("Exception: properties[arc].size() != matchingProperties.size().");
        }

        if(properties[arc] == matchingProperties){
            ++nbMatches;
            result = arc;
        }

    };

    if(lookup) {
        for (ListDigraph::Arc arc = lookup->operator()(from, to); arc != lemon::INVALID; arc = lookup->operator()(from, to, arc)) {
            processArc(arc);
        }
    }
    else{
        for(lemon::ConArcIt<ListDigraph> arcIt(graph, from, to); arcIt != lemon::INVALID; ++arcIt){
            processArc(arcIt);
        }
    }

    if(nbMatches != 1){
        throw std::runtime_error("Exception: nbMatches == " + std::to_string(nbMatches) + " != 1.");
    }

    return result;

}

ListDigraph_VArc SuperGraph::getMatchingArcs(const ListDigraph::Node& from, const ListDigraph::Node& to, const Vint& matchingProperties, const Vbool& matchOnElements) const {

    if(matchingProperties.size() != matchOnElements.size()){
        throw std::runtime_error("Exception: matchingProperties.size() != matchOnElements.size().");
    }

    const int& nbProperties = matchingProperties.size();

    ListDigraph_VArc results = ListDigraph_VArc();

    auto processArc = [&](const ListDigraph::Arc& arc){

        if(properties[arc].size() != nbProperties){
            throw std::runtime_error("Exception: properties[arc].size() != nbProperties.");
        }

        bool match = true;

        for(int i = 0; i < nbProperties; ++i){
            if(matchOnElements[i]){
                if(properties[arc].at(i) != matchingProperties.at(i)){
                    match = false;
                }
            }
        }

        if(match){
            results.push_back(arc);
        }

    };

    if(lookup) {
        for (ListDigraph::Arc arc = lookup->operator()(from, to); arc != lemon::INVALID; arc = lookup->operator()(from, to, arc)) {
            processArc(arc);
        }
    }
    else{
        for(lemon::ConArcIt<ListDigraph> arcIt(graph, from, to); arcIt != lemon::INVALID; ++arcIt){
            processArc(arcIt);
        }
    }

    return results;

}

void SuperGraph::newArcLookupTable() {
    lookup = std::make_shared<lemon::AllArcLookUp<ListDigraph>>(graph);
}

void SuperGraph::deleteArcLookupTable() {
    lookup = nullptr;
}