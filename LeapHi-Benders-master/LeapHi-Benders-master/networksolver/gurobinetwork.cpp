#include "gurobinetwork.hpp"

#include "networksolver_cpp_include.hpp"
#include "kevintools_primitives.hpp"

using namespace networksolver;

GurobiNetwork::GurobiNetwork(const pSuperGraph& superGraph,
                             const ListDigraph::NodeMap<bool>& nodeFilter,
                             const ListDigraph::ArcMap<bool>& arcFilter,
                             const ListDigraph_VNode& sourceNode,
                             const Vdouble& initialSourceFlows,
                             const ListDigraph_VNode& sinkNodes,
                             const Vdouble& initialSinkFlows,
                             const int maxCardinality,
                             const pGRBEnv& env,
                             const VResourceBounds& initialResourceBounds) :
    Network(superGraph, nodeFilter, arcFilter, sourceNode, initialSourceFlows,
            sinkNodes, initialSinkFlows, maxCardinality),
    env(env), model(std::make_shared<GRBModel>(*env)),
    arcFlows(layeredGraph), flowConstraints(layeredGraph), resourceConstraints()
{
    initializeVariablesAndObjective();
    initializeFlowConstraints();

    model->update(); //ensure that the model is updated before modifying flows
    updateSourceFlows(initialSourceFlows);
    updateSinkFlows(initialSinkFlows);

    initializeResourceConstraints(initialResourceBounds);
    model->update(); //ensure that the model is updated after initializing resource constraints

    model->set(GRB_IntParam_OutputFlag, 0);
    model->set(GRB_IntParam_Method, 1); //dual simplex
}

void GurobiNetwork::initializeResourceConstraints(const networksolver::VResourceBounds& externalResourceBounds){

    const int& nbResources = resourceIndices.size();
    if(nbResources == 0){
        return;
    }

    VResourceBounds localResourceBounds = externalToLocalResourceBounds(externalResourceBounds);
    resourceConstraints.reserve(nbLayers());

    for(int q = 0; q < nbLayers(); ++q){
        VGRBLinExpr resourceExpressions = buildResourceExpressionsLayer(q);
        resourceConstraints.push_back(kevintools::addExpressionBoundConstraints(model, resourceExpressions, localResourceBounds));
    }

}

void GurobiNetwork::updateResourceBounds(const VResourceBounds& resourceBounds, const bool requireTightBound){
    updateResourceConstraints(resourceBounds);
}

void GurobiNetwork::updateSourceFlows(const Vdouble& sourceFlows){

    Network::updateSourceFlows(sourceFlows);

    for(int i = 0; i < layeredSources.size(); ++i){
        const StaticDigraph::Node& node = layeredSources.at(i);
        GRBConstr& constraint = flowConstraints[node];
        constraint.set(GRB_DoubleAttr_RHS, sourceFlows.at(i));
    }

    model->update();

}

void GurobiNetwork::updateSinkFlows(const Vdouble& sinkFlows){

    Network::updateSinkFlows(sinkFlows);

    for(int i = 0; i < layeredSinks.size(); ++i){
        const StaticDigraph::Node& node = layeredSinks.at(i);
        GRBConstr& constraint = flowConstraints[node];
        constraint.set(GRB_DoubleAttr_RHS, sinkFlows.at(i));
    }

    model->update();

}

bool GurobiNetwork::solveModel() {
    model->optimize();
    return (getStatus() == 2);
}

int GurobiNetwork::getStatus() const {
    return model->get(GRB_IntAttr_Status);
}

double GurobiNetwork::getObjValue() const {
    return model->get(GRB_DoubleAttr_ObjVal);
}

void GurobiNetwork::writeModelToFile(const std::string& filename) const{
    model->write(filename);
}

pListDigraph_ArcMap_double GurobiNetwork::getArcFlowValues() const {

    //pointer to prevent copying
    std::shared_ptr<ListDigraph::ArcMap<double>> pSuperArcFlowValues = std::make_shared<ListDigraph::ArcMap<double>>(superGraph->graph);

    for (StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {
        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        (*pSuperArcFlowValues)[superArc] += arcFlows[arcIt].get(GRB_DoubleAttr_X);
    }

    return pSuperArcFlowValues;

}

ListDigraph_VArc GurobiNetwork::getArcFlowPath(const StaticDigraph::Node& source, const StaticDigraph::Node& sink) const{
    return Network::getArcFlowPath(source, sink);
}

Vdouble GurobiNetwork::getResourceConstraintDuals(const constraint_sense& sense) const{

    const int& nbResources = resourceIndices.size();
    Vdouble result;
    result.reserve(nbResources);

    for(int i = 0; i < nbResources; ++i) {

        double aggregateDualOverLayers = 0;

        for (int q = 0; q < nbLayers(); ++q) {

            const SenseConstraintMap_GRB& constraints = resourceConstraints.at(q).at(i);

            auto it = constraints.find(sense);
            if (it == constraints.end()) {
                aggregateDualOverLayers += 0; //use zero duals if constraint does not exist
            }
            else {
                aggregateDualOverLayers += it->second.get(GRB_DoubleAttr_Pi);
            }

        }

        result.push_back(aggregateDualOverLayers);

    }

    return result;

}

void GurobiNetwork::initializeVariablesAndObjective() {

    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const int& i = superGraph->graph.id(superGraph->graph.source(superArc));
        const int& j = superGraph->graph.id(superGraph->graph.target(superArc));
        const int& layer = originLayer[arcIt];
        const Vint& properties = superGraph->properties[superArc];

        std::string name = "x_" + std::to_string(layer) + "_" + std::to_string(i) + "_" + std::to_string(j);
        for(int k = 0; k < properties.size(); ++k){
            name += "_" + std::to_string(properties.at(k));
        }
        const double& cost = superGraph->cost[superArc];
        arcFlows[arcIt] = model->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS, name);

    }

}

void GurobiNetwork::initializeFlowConstraints() {

    for(StaticDigraph::NodeIt nodeIt(layeredGraph); nodeIt != lemon::INVALID; ++nodeIt){

        GRBLinExpr flowExpression = GRBLinExpr();

        for(StaticDigraph::OutArcIt arcIt(layeredGraph, nodeIt); arcIt != lemon::INVALID; ++arcIt){
            flowExpression += arcFlows[arcIt];
        }

        for(StaticDigraph::InArcIt arcIt(layeredGraph, nodeIt); arcIt != lemon::INVALID; ++arcIt){
            flowExpression -= arcFlows[arcIt];
        }

        double rhs = 0;
        if(isLayeredSource(nodeIt)){
            rhs = getLayeredSourceFlow(nodeIt);
        }
        else if(isLayeredSink(nodeIt)){
            rhs = getLayeredSinkFlow(nodeIt);
        }

        flowConstraints[nodeIt] = kevintools::addExpressionEqConstraint(model, flowExpression, rhs);

    }

}

VGRBLinExpr GurobiNetwork::buildResourceExpressionsLayer(const int layer) const{

    const int& nbResources = resourceIndices.size();

    VGRBLinExpr resourceExpressions = VGRBLinExpr(nbResources);
    
    for(StaticDigraph::ArcIt arcIt(layeredGraph); arcIt != lemon::INVALID; ++arcIt) {

        const ListDigraph::Arc& superArc = layeredToSuperArc(arcIt);
        const int& arcLayer = originLayer[arcIt];

        if (arcLayer != layer) {
            continue;
        }

        const VResourceContribution& resourceContributions = superGraph->contributions[superArc];
        const int& nbContributions = resourceContributions.size();

        for (int l = 0; l < nbContributions; ++l) {

            const ResourceContribution& contribution = resourceContributions.at(l);
            const int& resource = contribution.at(0);

            const int& localResource = externalToLocalResourceMap.at(resource);
            if (localResource != -1) {
                
                const int& coefficient = contribution.at(1);
                if (coefficient != 0) {
                    resourceExpressions.at(localResource) += (coefficient * arcFlows[arcIt]);
                }

            }

        }

    }

    return resourceExpressions;

}

void GurobiNetwork::updateResourceConstraints(const VResourceBounds& externalResourceBounds){

    const int& nbIndices = resourceIndices.size();
    VResourceBounds localResourceBounds = externalToLocalResourceBounds(externalResourceBounds);

    for(int q = 0; q < nbLayers(); ++q) {
        for(int i = 0; i < nbIndices; ++i) {

            SenseConstraintMap_GRB& constraints = resourceConstraints.at(q).at(i);

            if(constraints.empty()){
                continue;
            }

            ResourceBounds oldBounds = {kevintools::NEG_INF_DOUBLE, kevintools::POS_INF_DOUBLE};
            for (auto& constraint : constraints) {

                const constraint_sense& sense = constraint.first;
                GRBConstr& constr = constraint.second;

                if (sense == LE) {
                    oldBounds.at(1) = constr.get(GRB_DoubleAttr_RHS);
                } else if (sense == GE) {
                    oldBounds.at(0) = constr.get(GRB_DoubleAttr_RHS);
                } else if (sense == EQ) {
                    oldBounds.at(0) = constr.get(GRB_DoubleAttr_RHS);
                    oldBounds.at(1) = oldBounds.at(0);
                }

            }

            const int& index = resourceIndices.at(i);
            ResourceBounds newBounds = externalResourceBounds.at(index);

            if (oldBounds != newBounds) {

                for (int j = 0; j < 2; ++j) {
                    if (kevintools::is_inf(oldBounds.at(j)) && !kevintools::is_inf(newBounds.at(j))) {
                        throw std::runtime_error("Exception: changing finite bound into inf is not yet supported.");
                    }
                    if (!kevintools::is_inf(oldBounds.at(j)) && kevintools::is_inf(newBounds.at(j))) {
                        throw std::runtime_error("Exception: changing inf bound into finite bound is not yet supported.");
                    }
                }


            }

            for (auto& constraint : constraints) {

                const constraint_sense& sense = constraint.first;
                GRBConstr& constr = constraint.second;
                const double& lb = newBounds.at(0);
                const double& ub = newBounds.at(1);

                if (sense == LE) {
                    constr.set(GRB_DoubleAttr_RHS, ub);
                } else if (sense == GE) {
                    constr.set(GRB_DoubleAttr_RHS, lb);
                } else if (sense == EQ) {
                    if (lb != ub) {
                        throw std::runtime_error("Exception: lb != ub");
                    }
                    constr.set(GRB_DoubleAttr_RHS, lb);
                }

            }

        }
    }

}