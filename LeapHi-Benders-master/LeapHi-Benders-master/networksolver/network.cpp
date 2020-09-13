#include "network.hpp"

#include "networksolver_cpp_include.hpp"

#include <numeric>

using namespace networksolver;

Network::Network(const pSuperGraph& superGraph,
                 const ListDigraph::NodeMap<bool>& nodeFilter,
                 const ListDigraph::ArcMap<bool>& arcFilter,
                 const ListDigraph_VNode& sourceNodes,
                 const Vdouble& initialSourceFlows,
                 const ListDigraph_VNode& sinkNodes,
                 const Vdouble& initialSinkFlows,
                 int maxCardinality) :
        layeredToSuperArcIds(layeredGraph), originLayer(layeredGraph), destinationLayer(layeredGraph),
        relevantSuperArcIds(), layeredToRelevantSuperArcIndex(layeredGraph),
        resourceIndices(), externalToLocalResourceMap(),
        layeredSources(sourceNodes.size(), lemon::INVALID), sourceFlows(initialSourceFlows),
        layeredSinks(sinkNodes.size(), lemon::INVALID), sinkFlows(initialSinkFlows),
        superGraph(superGraph), maxCardinality(maxCardinality)
{

    if(sourceNodes.size() != initialSourceFlows.size()){
        throw std::runtime_error("Exception: sourceNodes.size() != initialSourceFlows.size().");
    }

    if(sinkNodes.size() != initialSinkFlows.size()){
        throw std::runtime_error("Exception: sinkNodes.size() != initialSinkFlows.size().");
    }

    if(!(maxCardinality == 0 || maxCardinality >= 2)){ //TODO: implement maxCardinality == 1 case
        throw std::runtime_error("Exception Network construction: maxCardinality must be 0 or >= 2.");
    }

    double totalSourceFlow = getTotalSourceFlow();
    double totalSinkFlow = getTotalSinkFlow();
    if(std::fabs(totalSourceFlow + totalSinkFlow) > 0.000001){
        throw std::runtime_error("Exception: totalSourceFlow != -totalSinkFlow.");
    }

    buildLayeredNetwork(nodeFilter, arcFilter, sourceNodes, sinkNodes);
    setRelevantSuperArcs();
    setResourceIndicesFromNetworkArcs();
    setExternalToLocalResourceMap();

}

void Network::updateSourceFlows(const Vdouble& sourceFlows) {

    this->sourceFlows = sourceFlows;

    if(sourceFlows.size() != layeredSources.size()){
        throw std::runtime_error("Exception: sourceFlows incorrect size");
    }

}

void Network::updateSinkFlows(const Vdouble& sinkFlows) {

    this->sinkFlows = sinkFlows;

    if(sinkFlows.size() != layeredSinks.size()){
        throw std::runtime_error("Exception: sinkFlows incorrect size");
    }

}

ListDigraph_VArc Network::getArcFlowPath() const{

    if(layeredSources.size() != 1 || layeredSinks.size() != 1){
        throw std::runtime_error("Exception: getArcFlowPath() called on network with multiple sources or sinks.");
    }

    return getArcFlowPath(layeredSources.at(0), layeredSinks.at(0));

}

ListDigraph_VArc Network::getArcFlowPath(const StaticDigraph::Node& source, const StaticDigraph::Node& sink) const {

    //TODO: implement
    throw std::runtime_error("Exception: getArcFlowPath not implemented.");
    return ListDigraph_VArc();

}

void Network::buildLayeredNetwork(const ListDigraph::NodeMap<bool>& nodeFilter,
                                  const ListDigraph::ArcMap<bool>& arcFilter,
                                  const ListDigraph_VNode& sourceNodes,
                                  const ListDigraph_VNode& sinkNodes){

    for(const ListDigraph::Node& sourceNode : sourceNodes){
        if(!nodeFilter[sourceNode]) {
            throw std::runtime_error("Exception Network construction: nodeFilter[sourceNode] == true expected.");
        }
    }
    for(const ListDigraph::Node& sinkNode : sinkNodes){
        if(!nodeFilter[sinkNode]) {
            throw std::runtime_error("Exception Network construction: nodeFilter[sinkNode] == true expected.");
        }
    }
    bool disjoint = true;
    for(const ListDigraph::Node& sourceNode : sourceNodes){
        for(const ListDigraph::Node& sinkNode : sinkNodes){
            if(sourceNode == sinkNode){
                disjoint = false;
            }
        }
    }
    if(!disjoint){
        std::cerr << "Warning Network construction: sinks and sources not disjoint." << std::endl;
    }

    //First build the constructionGraph as a ListDigraph (which is modifiable), then copy to layeredGraph (which is static).
    ListDigraph constructionGraph;

    std::shared_ptr<NodeToVNodeMap> pSuperToConstructionNodes =
            addNodesToFilteredLayeredCopy(superGraph->graph, constructionGraph, nodeFilter, arcFilter, sourceNodes, sinkNodes);

    std::shared_ptr<ArcToIntMap> pConstructionToSuperArcs =
            addArcsToFilteredLayeredCopy(superGraph->graph, constructionGraph, pSuperToConstructionNodes, nodeFilter, arcFilter);

    std::shared_ptr<ArcToIntMap> pConstructionOriginLayer =
            getCopyOriginLayerMap(superGraph->graph, constructionGraph, pSuperToConstructionNodes);

    std::shared_ptr<ArcToIntMap> pConstructionDestinationLayer =
            getCopyDestinationLayerMap(superGraph->graph, constructionGraph, pSuperToConstructionNodes);

    //Copy to layeredGraph
    lemon::DigraphCopy<ListDigraph, StaticDigraph> graphCopier(constructionGraph, layeredGraph);
    graphCopier.arcMap(*pConstructionToSuperArcs, layeredToSuperArcIds);
    graphCopier.arcMap(*pConstructionOriginLayer, originLayer);
    graphCopier.arcMap(*pConstructionDestinationLayer, destinationLayer);
    ListDigraph::NodeMap<StaticDigraph::Node> constructionToLayeredNodes(constructionGraph);
    graphCopier.nodeRef(constructionToLayeredNodes);
    graphCopier.run();

    //Define sources and sinks in layeredGraph
    for(int i = 0; i < sourceNodes.size(); ++i){
        const ListDigraph::Node& superSourceNode = sourceNodes.at(i);
        const ListDigraph::Node& constructionSourceNode = (*pSuperToConstructionNodes)[superSourceNode].at(0);
        layeredSources.at(i) = constructionToLayeredNodes[constructionSourceNode];
    }
    for(int i = 0; i < sinkNodes.size(); ++i){
        const ListDigraph::Node& superSinkNode = sinkNodes.at(i);
        const ListDigraph::Node& constructionSinkNode = (*pSuperToConstructionNodes)[superSinkNode].at(nbLayers()-1);
        layeredSinks.at(i) = constructionToLayeredNodes[constructionSinkNode];
    }

}

const ListDigraph::Arc Network::layeredToSuperArc(const StaticDigraph::Arc& layeredArc) const{
    return superGraph->graph.arcFromId(layeredToSuperArcIds[layeredArc]);
}

const ListDigraph::Arc Network::relevantSuperArc(const int index) const {
    return superGraph->graph.arcFromId(relevantSuperArcIds.at(index));
}

void Network::setRelevantSuperArcs() {

    relevantSuperArcIds = Vint();
    intintMap idToRelevantIndex = intintMap(); //map relevant super arcs ids to new index

    for (StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {
        
        const int superArcId = layeredToSuperArcIds[arcIt];

        auto iter = idToRelevantIndex.find(superArcId);
        if (iter != idToRelevantIndex.end()) {
            layeredToRelevantSuperArcIndex[arcIt] = iter->second;
        }
        else {
            idToRelevantIndex.insert({superArcId, relevantSuperArcIds.size()});
            layeredToRelevantSuperArcIndex[arcIt] = relevantSuperArcIds.size();
            relevantSuperArcIds.push_back(superArcId);
        }

    }

}

bool Network::useLayers() const {
    return maxCardinality > 0;
}

int Network::nbLayers() const {
    return maxCardinality + 1;
}

std::shared_ptr<Network::NodeToVNodeMap> Network::addNodesToFilteredLayeredCopy(const ListDigraph& original, ListDigraph& copy,
                                                                                const ListDigraph::NodeMap<bool>& nodeFilter,
                                                                                const ListDigraph::ArcMap<bool>& arcFilter,
                                                                                const ListDigraph_VNode& sourceNodes,
                                                                                const ListDigraph_VNode& sinkNodes) const
{

    std::shared_ptr<NodeToVNodeMap> pOriginalToCopyNodes = std::make_shared<NodeToVNodeMap>(original);

    auto createNodes = [&copy, pOriginalToCopyNodes](int begin, int end, const ListDigraph::Node& originalNode){
        for(int layer = begin; layer < end; ++layer) {
            ListDigraph::Node copyNode = copy.addNode();
            (*pOriginalToCopyNodes)[originalNode].at(layer) = copyNode;
        }
    };

    for(ListDigraph::NodeIt originalNodeIt(original); originalNodeIt != lemon::INVALID; ++originalNodeIt){
        if(nodeFilter[originalNodeIt]){

            (*pOriginalToCopyNodes)[originalNodeIt] = ListDigraph_VNode(nbLayers(), lemon::INVALID);
            const bool isSource = kevintools::contains(sourceNodes, originalNodeIt);
            const bool isSink = kevintools::contains(sinkNodes, originalNodeIt);

            if(!useLayers()){
                createNodes(0, 1, originalNodeIt);
            }
            else if(isSource && isSink){
                createNodes(0, nbLayers(), originalNodeIt);
            }
            else if (isSource && !isSink){

                if(hasInArcInFilter(originalNodeIt, original, arcFilter)){
                    createNodes(0, nbLayers()-1, originalNodeIt);
                }
                else{
                    createNodes(0, 1, originalNodeIt);
                }

            }
            else if(!isSource && isSink){

                if(hasOutArcInFilter(originalNodeIt, original, arcFilter)){
                    createNodes(1, nbLayers(), originalNodeIt);
                }
                else{
                    createNodes(nbLayers()-1, nbLayers(), originalNodeIt);
                }

            }
            else{
                createNodes(1, nbLayers()-1, originalNodeIt);
            }

        }
    }

    return pOriginalToCopyNodes;

}

bool Network::hasInArcInFilter(const ListDigraph::Node& node, const ListDigraph& graph, const ListDigraph::ArcMap<bool>& arcFilter) const{

    for(ListDigraph::InArcIt it(graph, node); it != lemon::INVALID; ++it){
        if(arcFilter[it]){
            return true;
        }
    }
    return false;

}

bool Network::hasOutArcInFilter(const ListDigraph::Node& node, const ListDigraph& graph, const ListDigraph::ArcMap<bool>& arcFilter) const{

    for(ListDigraph::OutArcIt it(graph, node); it != lemon::INVALID; ++it){
        if(arcFilter[it]){
            return true;
        }
    }
    return false;

}

std::shared_ptr<Network::ArcToIntMap> Network::addArcsToFilteredLayeredCopy(const ListDigraph& original, ListDigraph& copy,
                                                                            const std::shared_ptr<NodeToVNodeMap>& pOriginalToCopyNodes,
                                                                            const ListDigraph::NodeMap<bool>& nodeFilter,
                                                                            const ListDigraph::ArcMap<bool>& arcFilter) const
{

    std::shared_ptr<ArcToIntMap> pCopyToOriginalArcs = std::make_shared<ArcToIntMap>(copy);

    auto addArcBetween = [&copy, pCopyToOriginalArcs]
            (int fromLayer, int toLayer, int originalArcId, const ListDigraph_VNode& copyFroms, const ListDigraph_VNode& copyTos) {

        const ListDigraph::Node& copyFrom = copyFroms.at(fromLayer);
        const ListDigraph::Node& copyTo = copyTos.at(toLayer);

        if (copyFrom != lemon::INVALID && copyTo != lemon::INVALID) {
            ListDigraph::Arc copyArc = copy.addArc(copyFrom, copyTo);
            (*pCopyToOriginalArcs)[copyArc] = originalArcId;
        }

    };

    for(ListDigraph::ArcIt originalArcIt(original); originalArcIt != lemon::INVALID; ++originalArcIt){
        if(arcFilter[originalArcIt]){

            const ListDigraph::Node& originalFrom = original.source(originalArcIt);
            const ListDigraph::Node& originalTo = original.target(originalArcIt);

            const ListDigraph_VNode& copyFroms = (*pOriginalToCopyNodes)[originalFrom];
            const ListDigraph_VNode& copyTos = (*pOriginalToCopyNodes)[originalTo];

            if(copyFroms.empty() || copyTos.empty()){
                throw std::runtime_error("Exception: arc has adjacent node that is not in the node filter.");
            }

            int originalArcId = original.id(originalArcIt);
            if(!useLayers()){
                addArcBetween(0, 0, originalArcId, copyFroms, copyTos);
            }
            else {
                for (int layer = 0; layer < nbLayers() - 2; ++layer) {
                    addArcBetween(layer, layer + 1, originalArcId, copyFroms, copyTos);
                }
                for (int layer = 0; layer < nbLayers() - 1; ++layer) {
                    addArcBetween(layer, nbLayers()-1, originalArcId, copyFroms, copyTos);
                }
            }

        }
    }

    return pCopyToOriginalArcs;

}

std::shared_ptr<Network::ArcToIntMap> Network::getCopyOriginLayerMap(const ListDigraph& original, ListDigraph& copy,
                                                                     const std::shared_ptr<NodeToVNodeMap>& pOriginalToCopyNodes) const
{

    std::shared_ptr<ArcToIntMap> pCopyOriginLayer = std::make_shared<ArcToIntMap>(copy);

    for(ListDigraph::NodeIt originalNodeIt(original); originalNodeIt != lemon::INVALID; ++originalNodeIt){

        const ListDigraph_VNode& copyNodes = (*pOriginalToCopyNodes)[originalNodeIt];

        for (int i = 0; i < copyNodes.size(); ++i) {

            const ListDigraph::Node& copyNode = copyNodes.at(i);

            if (copyNode != lemon::INVALID) {
                for (ListDigraph::OutArcIt copyArcIt(copy, copyNode); copyArcIt != lemon::INVALID; ++copyArcIt) {
                    (*pCopyOriginLayer)[copyArcIt] = i;
                }
            }

        }

    }

    return pCopyOriginLayer;

}

std::shared_ptr<Network::ArcToIntMap> Network::getCopyDestinationLayerMap(const ListDigraph& original, ListDigraph& copy,
                                                                     const std::shared_ptr<NodeToVNodeMap>& pOriginalToCopyNodes) const
{

    std::shared_ptr<ArcToIntMap> pCopyDestinationLayer = std::make_shared<ArcToIntMap>(copy);

    for(ListDigraph::NodeIt originalNodeIt(original); originalNodeIt != lemon::INVALID; ++originalNodeIt){

        const ListDigraph_VNode& copyNodes = (*pOriginalToCopyNodes)[originalNodeIt];

        for (int i = 0; i < copyNodes.size(); ++i) {

            const ListDigraph::Node& copyNode = copyNodes.at(i);

            if (copyNode != lemon::INVALID) {
                for (ListDigraph::InArcIt copyArcIt(copy, copyNode); copyArcIt != lemon::INVALID; ++copyArcIt) {
                    (*pCopyDestinationLayer)[copyArcIt] = i;
                }
            }

        }

    }

    return pCopyDestinationLayer;

}

void Network::setResourceIndicesFromNetworkArcs() {

    intSet uniqueResourceIndices = intSet();

    for (StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);

        const VResourceContribution& resourceContributions = superGraph->contributions[superArc];
        const int& nbContributions = resourceContributions.size();

        for (int l = 0; l < nbContributions; ++l) {
            const ResourceContribution& contribution = resourceContributions.at(l);
            const int& resouceIndex = contribution.at(0);
            uniqueResourceIndices.insert(resouceIndex);
        }

    }

    resourceIndices = Vint(uniqueResourceIndices.begin(), uniqueResourceIndices.end());

}

VResourceBounds Network::externalToLocalResourceBounds(const networksolver::VResourceBounds& externalResourceBounds) const {

    const int& nbIndices = resourceIndices.size();
    VResourceBounds localResourceBounds = VResourceBounds(nbIndices);

    for(int i = 0; i < nbIndices; ++i) {
        const int& index = resourceIndices.at(i);
        localResourceBounds.at(i) = externalResourceBounds.at(index);
    }

    return localResourceBounds;

}

void Network::setExternalToLocalResourceMap(){

    if(resourceIndices.empty()){
        clearExternalToLocalResourceMap();
        return;
    }

    int maxExternalResource = *std::max_element(resourceIndices.begin(), resourceIndices.end());
    externalToLocalResourceMap = Vint(maxExternalResource + 1, -1);

    const int& nbLocalResources = resourceIndices.size();
    for(int r = 0; r < nbLocalResources; ++r){
        const int& localResource = r;
        const int& externalResource = resourceIndices.at(r);
        externalToLocalResourceMap.at(externalResource) = localResource;
    }

}

void Network::clearExternalToLocalResourceMap() {
    externalToLocalResourceMap.clear();
}

const Vint& Network::getResourceIndices() const {
    return resourceIndices;
}

bool Network::isLayeredSource(const StaticDigraph::Node& node) const{
    return std::find(layeredSources.begin(), layeredSources.end(), node) != layeredSources.end();
}

bool Network::isLayeredSink(const StaticDigraph::Node& node) const{
    return std::find(layeredSinks.begin(), layeredSinks.end(), node) != layeredSinks.end();
}

double Network::getLayeredSourceFlow(const StaticDigraph::Node& node) const {

    if(!isLayeredSource(node)){
        throw std::runtime_error("Exception: getLayeredSourceFlow called for non source flow.");
    }

    int index = std::distance(layeredSources.begin(), std::find(layeredSources.begin(), layeredSources.end(), node));
    return sourceFlows.at(index);

}

double Network::getLayeredSinkFlow(const StaticDigraph::Node& node) const {

    if(!isLayeredSink(node)){
        throw std::runtime_error("Exception: getLayeredSinkFlow called for non sink flow.");
    }

    int index = std::distance(layeredSinks.begin(), std::find(layeredSinks.begin(), layeredSinks.end(), node));
    return sinkFlows.at(index);

}

double Network::getTotalSourceFlow() const {
    return std::accumulate(sourceFlows.begin(), sourceFlows.end(), 0);
}

double Network::getTotalSinkFlow() const {
    return std::accumulate(sinkFlows.begin(), sinkFlows.end(), 0);
}
