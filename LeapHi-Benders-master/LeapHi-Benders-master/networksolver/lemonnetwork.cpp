#include "lemonnetwork.hpp"

#include "networksolver_cpp_include.hpp"
#include "kevintools_primitives.hpp"

#include <numeric>

using namespace networksolver;

LemonNetwork::LemonNetwork(const pSuperGraph& superGraph,
                           const ListDigraph::NodeMap<bool>& nodeFilter,
                           const ListDigraph::ArcMap<bool>& arcFilter,
                           const ListDigraph_VNode& sourceNode,
                           const Vdouble& initialSourceFlows,
                           const ListDigraph_VNode& sinkNodes,
                           const Vdouble& initialSinkFlows,
                           const int maxCardinality,
                           const int costScaleFactor,
                           const int capacityScaleFactor) :
        Network(superGraph, nodeFilter, arcFilter, sourceNode, initialSourceFlows,
            sinkNodes, initialSinkFlows, maxCardinality),
        objValue(-1), resourceConstraintDuals(resourceIndices.size(), -1),
        costScaleFactor(costScaleFactor), capacityScaleFactor(capacityScaleFactor),
        capacityMap(nullptr), roundedCostMap(nullptr), supplyMap(nullptr),
        status(NetworkSimplex::INFEASIBLE), capacityCorrections(),
        negativeArcCosts(false), useDijkstra(false), preventSolve(true),
        unroundedCostMap(nullptr), filteredGraph(nullptr),
        mcfAlgorithm(nullptr), dijkstraAlgorithm(nullptr)
{

    verifyFlowsIntegral();
    verifySuperGraphAssumptions();

    initializeRoundedCostMap();
    initializeUnroundedCostMap();
    if(useSupplyMap()){
        updateSupplyMap();
    }

}

void LemonNetwork::updateResourceBounds(const VResourceBounds& resourceBounds,  const bool requireTightBound) {

    useDijkstra = requireTightBound;

    if(useDijkstra){
        if(useSupplyMap()){ //more than one source or sink
            useDijkstra = false;
        }
        if(sourceFlows.at(0) != 1 || sinkFlows.at(0) != -1){
            useDijkstra = false;
        }
        if(negativeArcCosts){
            useDijkstra = false;
        }
    }

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt){

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const VResourceContribution& contributions = superGraph->contributions[superArc];

        for (const ResourceContribution& contribution : contributions) {

            const int& resource = contribution.at(0);
            const ResourceBounds& newBounds = resourceBounds.at(resource);

            if (!kevintools::is_neginf(newBounds.at(0))) {
                throw std::runtime_error("Exception: only LE is currently supported by LemonNetwork.");
            } else if (kevintools::is_inf(newBounds.at(1))) {
                throw std::runtime_error("Exception: infinite bounds not currently supported by LemonNetwork.");
            }

            if(useDijkstra) {
                double doubleBound = newBounds.at(1);
                int intBound = ((int) newBounds.at(1));

                if (doubleBound != intBound || (intBound != 0 && intBound != 1)) {
                    useDijkstra = false;
                }
            }

        }

    }

    if(requireTightBound && !useDijkstra){
        throw std::runtime_error("Exception: requireTightBound = true, but Dijkstra cannot be used to generate a tight bound.");
    }

    if(useDijkstra){
        capacityMap = nullptr;
        capacityCorrections.clear();
        updateArcFilter(resourceBounds);
    }
    else{
        updateCapacityAndCapacityCorrections(resourceBounds);
        arcFilter = nullptr;
        filteredGraph = nullptr;
    }

    preventSolve = false;

}

void LemonNetwork::updateSourceFlows(const Vdouble& sourceFlows) {

    Network::updateSourceFlows(sourceFlows);

    if(useSupplyMap()){
        updateSupplyMap();
    }

    preventSolve = true;
    mcfAlgorithm = nullptr;

}

void LemonNetwork::updateSinkFlows(const Vdouble& sinkFlows) {

    Network::updateSinkFlows(sinkFlows);

    if(useSupplyMap()){
        updateSupplyMap();
    }

    preventSolve = true;
    mcfAlgorithm = nullptr;

}

bool LemonNetwork::solveModel() {

    if(preventSolve){
        throw std::runtime_error("Exception: updateResourceBounds() need to be called before solveModel().");
    }

    const StaticDigraph::Node& source = layeredSources.at(0);
    const StaticDigraph::Node& sink = layeredSinks.at(0);

    if(useDijkstra){

        filteredGraph = std::make_shared<FilterArcs>(layeredGraph, *arcFilter);
        dijkstraAlgorithm = std::make_shared<FilterArcsDijkstra>(*filteredGraph, *unroundedCostMap);
        //mcfAlgorithm = nullptr;

        if(useSupplyMap()){
            throw std::runtime_error("Exception: useSupplyMap = true should not be possible at this point.");
        }

        dijkstraAlgorithm->run(source);
        status = (int) !dijkstraAlgorithm->reached(sink);
        setResourceConstraintDuals();
        setObjValue(); //requires setResourceConstraintDuals() to be called first

    }
    else{

        dijkstraAlgorithm = nullptr;

        if(mcfAlgorithm == nullptr) {

            mcfAlgorithm = std::make_shared<NetworkSimplex>(layeredGraph);
            mcfAlgorithm->costMap(*roundedCostMap);

            if (useSupplyMap()) {
                mcfAlgorithm->supplyMap(*supplyMap);
            } else {
                int flow = (int) std::round(sourceFlows.at(0) * capacityScaleFactor);
                mcfAlgorithm->stSupply(source, sink, flow);
            }

        }

        mcfAlgorithm->upperMap(*capacityMap);

        status = (int) (mcfAlgorithm->run() != NetworkSimplex::OPTIMAL);
        setResourceConstraintDuals();
        setObjValue(); //requires setResourceConstraintDuals() to be called first

    }

    preventSolve = true;

    return (status == 0);

}

int LemonNetwork::getStatus() const {
    return status;
}

double LemonNetwork::getObjValue() const {
    return objValue;
}

void LemonNetwork::writeModelToFile(const std::string& filename) const{
    throw std::runtime_error("Exception: writeModelToFile() not implemented for LemonNetwork.");
}

pListDigraph_ArcMap_double LemonNetwork::getArcFlowValues() const {

    //pointer to prevent copying
    pListDigraph_ArcMap_double pSuperArcFlowValues = std::make_shared<ListDigraph::ArcMap<double>>(superGraph->graph);

    StaticDigraph::ArcMap<int> arcFlows(layeredGraph);

    const StaticDigraph::Node& source = layeredSources.at(0);
    const StaticDigraph::Node& sink = layeredSinks.at(0);

    if(useDijkstra){
        for (StaticDigraph::Node node = sink; node != source; node = dijkstraAlgorithm->predNode(node)) {
            const StaticDigraph::Arc& arc = dijkstraAlgorithm->predArc(node);
            arcFlows[arc] = 1;
        }
    }
    else {
        mcfAlgorithm->flowMap(arcFlows);
    }

    for (StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {
        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        if(useDijkstra) {
            (*pSuperArcFlowValues)[superArc] += arcFlows[arcIt];
        }
        else{
            (*pSuperArcFlowValues)[superArc] += ((double) arcFlows[arcIt]) / capacityScaleFactor;
        }
    }

    return pSuperArcFlowValues;

}

ListDigraph_VArc LemonNetwork::getArcFlowPath(const StaticDigraph::Node& source, const StaticDigraph::Node& sink) const{

    if(!useDijkstra){
        throw std::runtime_error("Exception: can only call LemonNetwork::getArcFlowPath after solving with Dijkstra.");
    }

    if(source == sink){
        throw std::runtime_error("Exception: source == sink.");
    }

    StaticDigraph_VArc layeredReverseOrder = StaticDigraph_VArc();
    for (StaticDigraph::Node node = sink; node != source; node = dijkstraAlgorithm->predNode(node)) {
        const StaticDigraph::Arc& arc = dijkstraAlgorithm->predArc(node);
        layeredReverseOrder.push_back(arc);
    }
    int nbArcs = layeredReverseOrder.size();

    ListDigraph_VArc superCorrectOrder;  superCorrectOrder.reserve(nbArcs);
    for(int i = nbArcs - 1; i >= 0; --i){
        const StaticDigraph::Arc& layeredArc = layeredReverseOrder.at(i);
        const ListDigraph::Arc& superArc = layeredToSuperArc(layeredArc);
        superCorrectOrder.push_back(superArc);
    }

    return superCorrectOrder;

}

Vdouble LemonNetwork::getResourceConstraintDuals(const constraint_sense& sense) const{

    if(sense != LE){
        throw std::runtime_error("Exception: only LE is currently supported by LemonNetwork.");
    }

    return resourceConstraintDuals;

}

void LemonNetwork::setObjValue() {

    const StaticDigraph::Node& sink = layeredSinks.at(0);

    if(useDijkstra) {
        objValue = dijkstraAlgorithm->dist(sink);
    }
    else{

        //mcfAlgorithm->totalCost() may result in integer overflow, so calculate manually

        StaticDigraph::ArcMap<int> arcFlows(layeredGraph);
        mcfAlgorithm->flowMap(arcFlows);

        objValue = 0;
        for (StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {
            objValue += (((double) arcFlows[arcIt]) / capacityScaleFactor) * (((double) (*roundedCostMap)[arcIt]) / costScaleFactor);
        }

    }

    //Perform capacity correction, not needed for Dijkstra
    if(!useDijkstra) {

        if(capacityCorrections.size() != relevantSuperArcIds.size()){
            throw std::runtime_error("Exception: capacityCorrections.size() != relevantSuperArcs.size().");
        }

        for (int i = 0; i < capacityCorrections.size(); ++i){

            const ListDigraph::Arc& superArc = relevantSuperArc(i);
            const VResourceContribution& contributions = superGraph->contributions[superArc];
            if (contributions.empty()) {
                continue;
            }
            const int& resource = contributions.at(0).at(0);

            const int& localResource = externalToLocalResourceMap.at(resource);
            if (localResource != -1) {
                double arcDual = resourceConstraintDuals.at(localResource);
                objValue += capacityCorrections.at(i) * arcDual;
            }

        }

    }

}

void LemonNetwork::setResourceConstraintDuals() {

    for(int i = 0; i < resourceConstraintDuals.size(); ++i){
        resourceConstraintDuals.at(i) = 0;
    }

    if(useDijkstra){
        setResourceConstraintDualsDijkstra();
    }
    else{
        setResourceConstraintDualsMCF();
    }

}

void LemonNetwork::setResourceConstraintDualsMCF() {

    if(useDijkstra){
        throw std::runtime_error("Exception: cannot call setResourceConstraintDualsMCF when useDijkstra = true.");
    }

    StaticDigraph::NodeMap<int> potentialMap(layeredGraph);
    mcfAlgorithm->potentialMap(potentialMap);

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt){

        const StaticDigraph::Node& from = layeredGraph.source(arcIt);
        const StaticDigraph::Node& to = layeredGraph.target(arcIt);

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const VResourceContribution& contributions = superGraph->contributions[superArc];

        int arcReducedCost = (*roundedCostMap)[arcIt] + potentialMap[from] - potentialMap[to];
        int scaledArcDual = (arcReducedCost >= 0 ? 0 : arcReducedCost);
        double arcDual = ((double) scaledArcDual) / costScaleFactor;

        if(contributions.empty()){
            continue;
        }
        const int& resource = contributions.at(0).at(0);
        const int& localResource = externalToLocalResourceMap.at(resource);
        if (localResource != -1) {
            resourceConstraintDuals.at(localResource) += arcDual; //aggregate all arcs with same resource
        }
        else{
            throw std::runtime_error("Exception: resource not registered in resourceIndices.");
        }

    }

}

void LemonNetwork::setResourceConstraintDualsDijkstra() {

    if(!useDijkstra){
        throw std::runtime_error("Exception: cannot call setResourceConstraintDualsMCF when useDijkstra = false.");
    }

    const StaticDigraph::Node& source = layeredSources.at(0);
    const StaticDigraph::Node& sink = layeredSinks.at(0);
    double sinkPotential = dijkstraAlgorithm->dist(sink);

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt){

        const StaticDigraph::Node& from = layeredGraph.source(arcIt);
        const StaticDigraph::Node& to = layeredGraph.target(arcIt);

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const VResourceContribution& contributions = superGraph->contributions[superArc];

        double fromPotential = (from == source || dijkstraAlgorithm->reached(from)) ?
                dijkstraAlgorithm->dist(from) : sinkPotential;
        double toPotential = (to == source || dijkstraAlgorithm->reached(to)) ?
                dijkstraAlgorithm->dist(to) : sinkPotential;
        //Note: according to the documentation, reached(source) may be false

        double arcReducedCost = (*unroundedCostMap)[arcIt] + fromPotential - toPotential;
        double arcDual = (arcReducedCost >= 0 ? 0 : arcReducedCost);

        if(contributions.empty()){
            continue;
        }
        const int& resource = contributions.at(0).at(0);
        const int& localResource = externalToLocalResourceMap.at(resource);
        if (localResource != -1) {
            resourceConstraintDuals.at(localResource) += arcDual; //aggregate all arcs with same resource
        }
        else{
            throw std::runtime_error("Exception: resource not registered in resourceIndices.");
        }

    }

}

void LemonNetwork::updateCapacityAndCapacityCorrections(const networksolver::VResourceBounds& externalResourceBounds) {

    if(capacityMap == nullptr){
        capacityMap = std::make_shared<StaticDigraph::ArcMap<int>>(layeredGraph);
    }
    capacityCorrections = Vdouble(relevantSuperArcIds.size(), 0);

    int upperBound = capacityScaleFactor * ((int) std::round(getTotalSourceFlow()));

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt){

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const VResourceContribution& contributions = superGraph->contributions[superArc];

        if(contributions.empty()){
            (*capacityMap)[arcIt] = upperBound;
        }
        else {
            for (const ResourceContribution& contribution : contributions) {

                const int& resource = contribution.at(0);
                const ResourceBounds& newBounds = externalResourceBounds.at(resource);

                if (!kevintools::is_neginf(newBounds.at(0))) {
                    throw std::runtime_error("Exception: only LE is currently supported by LemonNetwork.");
                } else if (kevintools::is_inf(newBounds.at(1))) {
                    throw std::runtime_error("Exception: infinite bounds not currently supported by LemonNetwork.");
                }

                int scaledCapacity = (int) std::round(newBounds.at(1) * capacityScaleFactor);
                (*capacityMap)[arcIt] = scaledCapacity;

                const int& relevantSuperArcsIndex = layeredToRelevantSuperArcIndex[arcIt];
                capacityCorrections.at(relevantSuperArcsIndex) = newBounds.at(1) - ((double) scaledCapacity) / capacityScaleFactor; ///overwritten multiple times with same value

            }
        }

    }

}

void LemonNetwork::updateArcFilter(const networksolver::VResourceBounds& externalResourceBounds){

    arcFilter = std::make_shared<StaticDigraph::ArcMap<bool>>(layeredGraph);

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt){

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const VResourceContribution& contributions = superGraph->contributions[superArc];

        if(contributions.empty()){
            (*arcFilter)[arcIt] = true;
        }
        else {
            for (const ResourceContribution& contribution : contributions) {

                const int& resource = contribution.at(0);
                const ResourceBounds& newBounds = externalResourceBounds.at(resource);

                if(newBounds.at(1) == 0){
                    (*arcFilter)[arcIt] = false;
                }
                else if(newBounds.at(1) == 1){
                    (*arcFilter)[arcIt] = true;
                }
                else{
                    throw std::runtime_error("Exception: newBounds.at(1) expect to be 0 or 1.");
                }

            }
        }

    }

}

void LemonNetwork::initializeRoundedCostMap() {

    roundedCostMap = std::make_shared<StaticDigraph::ArcMap<int>>(layeredGraph);

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const double& cost = superGraph->cost[superArc];

        if(kevintools::is_inf(cost)){
            throw std::runtime_error("Exception: infinite cost not currently supported by LemonNetwork.");
        }

        (*roundedCostMap)[arcIt] = (int) std::floor(cost * costScaleFactor);

    }

}

void LemonNetwork::initializeUnroundedCostMap() {

    unroundedCostMap = std::make_shared<StaticDigraph::ArcMap<double>>(layeredGraph);

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        (*unroundedCostMap)[arcIt] = superGraph->cost[superArc];

    }

}

bool LemonNetwork::useSupplyMap() const {
    return (layeredSources.size() != 1 || layeredSinks.size() != 1);
}

void LemonNetwork::updateSupplyMap() {

    if(!useSupplyMap()){
        std::cerr << "Warning: updateSupplyMap called while useSupplyMap is false." << std::endl;
    }

    supplyMap = std::make_shared<StaticDigraph::NodeMap<int>>(layeredGraph);

    for(int i = 0; i < layeredSources.size(); ++i){
        (*supplyMap)[layeredSources.at(i)] = (int) std::round(sourceFlows.at(i) * capacityScaleFactor);
    }

    for(int i = 0; i < layeredSinks.size(); ++i){
        (*supplyMap)[layeredSinks.at(i)] = (int) std::round(sinkFlows.at(i) * capacityScaleFactor);
    }

}

void LemonNetwork::verifyFlowsIntegral() const {

    for(const double& flow : sourceFlows){
        if(flow != std::round(flow)){
            throw std::runtime_error("Exception: only integer source flows are supported.");
        }
    }

    for(const double& flow : sinkFlows){
        if(flow != std::round(flow)){
            throw std::runtime_error("Exception: only integer sink flows are supported.");
        }
    }

}

void LemonNetwork::verifySuperGraphAssumptions() {

    int count = 0;
    intSet resourcesSeen = intSet();

    for(ListDigraph::ArcIt superArcIt(superGraph->graph); superArcIt != lemon::INVALID; ++superArcIt) {

        const VResourceContribution& contributions = superGraph->contributions[superArcIt];
        if(contributions.empty()){
            continue;
        }
        if(contributions.size() > 1){
            throw std::runtime_error("Exception: superGraph assumption 2 violated.");
        }

        const int& resource = contributions.at(0).at(0);
        const int& value = contributions.at(0).at(1);

        resourcesSeen.insert(resource);
        ++count;
        if(resourcesSeen.size() != count){
            throw std::runtime_error("Exception: superGraph assumption 1 violated.");
        }

        if(value != 1){
            throw std::runtime_error("Exception: superGraph assumption 3 violated.");
        }

        if(superGraph->cost[superArcIt] < 0){
            negativeArcCosts = true;
        }

    }

    if(negativeArcCosts){
        std::cerr << "Warning: superGraph assumption 4 violated, may throw exception on Dijkstra algorithm. If no exception is thrown, all should be fine." << std::endl;
    }

}

ListDigraph_VVArc LemonNetwork::getAllShortestPathsShuttleFiltered(const intSet& superHubIDs){

    ListDigraph_VVArc result = ListDigraph_VVArc();
    result.reserve(layeredSources.size() * layeredSinks.size());

    useDijkstra = true;

    const ListDigraph& graph = superGraph->graph;

    arcFilter = std::make_shared<StaticDigraph::ArcMap<bool>>(layeredGraph, true);
    filteredGraph = std::make_shared<FilterArcs>(layeredGraph, *arcFilter);
    dijkstraAlgorithm = std::make_shared<FilterArcsDijkstra>(*filteredGraph, *unroundedCostMap);

    for(const StaticDigraph::Node& origin : layeredSources){

        StaticDigraph::OutArcIt originArc(layeredGraph, origin);
        int superOriginId = graph.id(graph.source(layeredToSuperArc(originArc)));

        for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt){

            const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
            int superFromId = graph.id(graph.source(superArc));
            int superToId = graph.id(graph.target(superArc));
            int mode = superGraph->properties[superArc].at(0);

            if(originLayer[arcIt] == 0 && superFromId != superOriginId) {
                (*arcFilter)[arcIt] = false; //Arc out of layer 0 node that is not the origin
            }
            else if(mode == 0 && originLayer[arcIt] != 0 && destinationLayer[arcIt] != nbLayers()-1){
                (*arcFilter)[arcIt] = false; //Shuttle arc that does not start from layer 0 or ends in last layer
            }
            else if(mode == 0 && superFromId == superOriginId && destinationLayer[arcIt] != nbLayers()-1 && !superHubIDs.count(superToId)){
                (*arcFilter)[arcIt] = false; //Shuttle arc out of the origin to non-last layer that is not a hub
            }
            else {
                (*arcFilter)[arcIt] = true;
            }

        }

        supplyMap = std::make_shared<StaticDigraph::NodeMap<int>>(layeredGraph);
        for(const StaticDigraph::Node& destination : layeredSinks){
            if(origin != destination){
                (*supplyMap)[destination] = -1;
            }
        }
        (*supplyMap)[origin] = -getTotalSinkFlow();

        dijkstraAlgorithm->run(origin);
        for(const StaticDigraph::Node& destination : layeredSinks){

            StaticDigraph::InArcIt destinationArc(layeredGraph, destination);
            int superDestinationId = graph.id(graph.target(layeredToSuperArc(destinationArc)));

            if(superOriginId != superDestinationId){

                if(!dijkstraAlgorithm->reached(destination)){
                    throw std::runtime_error("Exception: shortest path algorithm failed.");
                }

                result.push_back(getArcFlowPath(origin, destination));

            }

        }

    }

    return result;

}
