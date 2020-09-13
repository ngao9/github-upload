#include "networkdesignproblem.hpp"

#include "networksolver_cpp_include.hpp"
#include "kevintools_primitives.hpp"
#include "kevintools_time.hpp"

#include "instance.hpp"
#include "stops.hpp"

#include <algorithm>

using namespace networkdesign;

void NetworkDesignProblem::setCPLEXsettings(const double mipGap, const double timeLimit) {

    cplex->setParam(IloCplex::Param::MIP::Tolerances::MIPGap, mipGap); //fraction (i.e., 0.01 = 1%)
    cplex->setParam(IloCplex::Param::TimeLimit, timeLimit); //seconds
    //cplex->setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);

    cplex->use(MyUserCutCallback(env, callback));
    cplex->use(MyLazyConstraintCallback(env, callback));
    cplex->use(MyHeuristicCallback(env, callback));

}

void NetworkDesignProblem::createResources() {

    const Vint& hubIndices = in.hubIndices;

    //Resources correspond to travel between hubs with BUS
    for(int i : hubIndices){
        for(int j : hubIndices) {
            if(i != j){
                if (in.arcModeExists(i, j, BUS)) {
                    createResourcesArcMode(i, j, BUS);
                }
            }
        }
    }

}

void NetworkDesignProblem::createResourcesArcMode(const int i, const int j, const travel_mode& mode) {

    if(mode != BUS){
        throw std::runtime_error("Exception: mode == BUS expected.");
    }

    bool dummy = false;
    const Vint& possibleNb = in.nbPossible(i, j, mode);

    for(const int& nb : possibleNb) {
        resourceMap.push_back(Aint5{i, j, (int) mode, nb, (int) dummy});
        if(resourceMap.size()-1 != nbResources){
            throw std::runtime_error("Exception: resource not registered correctly.");
        }
        ++nbResources;
    }

}

NetworkDesignProblem::NetworkDesignProblem(const Instance& in,
                                           const double mipGap, const double timeLimit,
                                           const SubproblemType& subProblemType):
        in(in),
        env(), model(std::make_shared<IloModel>(env)), cplex(std::make_shared<IloCplex>(*model)),
        resourceVariables(), callback(std::make_shared<Callbacks>(*this)),
        superGraph(in.buildSuperGraph()),
#ifdef GUROBI
        subProblemEnv(nullptr),
#endif
        subproblemType(subProblemType), subProblems(), nbPeoples(),
        nbResources(0), resourceMap(), resourceBounds()
{

    if(subProblemType == GUROBINETWORK){
#ifdef GUROBI
        subProblemEnv = std::make_shared<GRBEnv>();
#else
        throw std::runtime_error("Exception: GUROBINETWORK cannot be used because Gurobi cannot be found.");
#endif
    }

    setCPLEXsettings(mipGap, timeLimit);
    createResources();
    verifyResourceMapUnique();
    setSuperGraphContributions();
    initializeResourceBounds();
    buildMasterProblem();
    callback->initializeCorePoint();
    //cplex->exportModel("master_problem.lp");
    buildSubProblems();

}

void NetworkDesignProblem::buildMasterProblem(){

    const int& n = in.n;
    const Vint& hubIndices = in.hubIndices;
    const int& nbHubs = hubIndices.size();

    resourceVariables.reserve(nbResources);
    VIloExpr flowExpressionsBus;  flowExpressionsBus.reserve(n); //for flow balance
    VIloExpr flowExpressionsTrain; flowExpressionsTrain.reserve(n);
    for(int i = 0; i < n; ++i){
        flowExpressionsBus.push_back(IloExpr(env));
        flowExpressionsTrain.push_back(IloExpr(env));
    }

    tripleIndexIloExprMap oneFrequencyPerModeExpressions = tripleIndexIloExprMap();
    IloExpr objective(env);

    for(int r = 0; r < nbResources; ++r){

        const Aint5& resourceInfo = resourceMap.at(r);

        const int& i = resourceInfo.at(0);
        const int& j = resourceInfo.at(1);
        const int& nbVehicles = resourceInfo.at(3);
        const travel_mode& mode = (travel_mode) resourceInfo.at(2);

        if(isDummyResource(resourceInfo)){

            //belongs to resource without variable in the master problem
            IloIntVar dummyVariable(env, 0, 0);
            resourceVariables.push_back(dummyVariable);
            kevintools::setName(dummyVariable, "z_" + std::to_string(r) + "_dummy");
            model->add(dummyVariable); //needed to ensure that the variable is in the model

            continue;

        }

        double cost = 0;
        if(mode == BUS || mode == RAIL){
            cost = in.arcObjective(i, j, mode, nbVehicles, Instance::DESIGN);
        } else{
            throw std::runtime_error("Exception: non-bus/train resources not supposed to be in the master problem.");
        }

        IloIntVar resourceVariable(env, 0, 1);
        resourceVariables.push_back(resourceVariable);
        kevintools::setName(resourceVariable, "z_" + std::to_string(r));
        objective += cost * resourceVariable;

        if(mode == BUS){
            flowExpressionsBus.at(i) += resourceVariables.at(r) * nbVehicles;
            flowExpressionsBus.at(j) -= resourceVariables.at(r) * nbVehicles;

        } else if(mode == RAIL){
            flowExpressionsTrain.at(i) += resourceVariables.at(r) * nbVehicles;
            flowExpressionsTrain.at(j) -= resourceVariables.at(r) * nbVehicles;
        }
        tripleIndex key = tripleIndex{ i, j, (int)mode };
        auto iter = oneFrequencyPerModeExpressions.find(key);
        if (iter == oneFrequencyPerModeExpressions.end()) { //not found
            IloExpr newExpressions(env);
            newExpressions += resourceVariables.at(r);
            oneFrequencyPerModeExpressions.insert(std::make_pair(key, newExpressions));
        }
        else {
            iter->second += resourceVariables.at(r);
        }

    }

    kevintools::addExpressionEqConstraints(model, flowExpressionsBus, Vdouble(n, 0));
    kevintools::clearExpressions(flowExpressionsBus);

    kevintools::addExpressionEqConstraints(model, flowExpressionsTrain, Vdouble(n, 0));
    kevintools::clearExpressions(flowExpressionsTrain);

    for (auto& keyVal : oneFrequencyPerModeExpressions) {

        IloExpr& expression = keyVal.second;
        kevintools::addExpressionBoundConstraint(model, expression, -IloInfinity, 1);

    }

    theta = IloNumVar(env, 0, IloInfinity, ILOFLOAT);
    kevintools::setName(theta, "theta");
    objective += theta;

    model->add(IloMinimize(env, objective));
    objective.end();

}

void NetworkDesignProblem::buildSubProblems() {

    const VAint3& odPairs = in.odPairs;
    const int nbODpairs = odPairs.size();

    subProblems = VpNetwork(nbODpairs);
    nbPeoples = Vint(nbODpairs);

    superGraph->newArcLookupTable();

    #pragma omp parallel for default(none) shared(odPairs)
    for(int i = 0; i < nbODpairs; ++i){

        const Aint3& odPair = odPairs.at(i);
        const int& origin = odPair.at(0);
        const int& destination = odPair.at(1);
        const int& nbPeople = odPair.at(2);

        subProblems.at(i) = buildSubProblem(origin, destination);
        nbPeoples.at(i) = nbPeople;

    }

    superGraph->deleteArcLookupTable();

}

pNetwork NetworkDesignProblem::buildSubProblem(const int origin, const int destination){

    if(origin == destination){
        throw std::runtime_error("Exception: OD pair found for which origin == destination.");
    }

    const Vint& hubIndices = in.hubIndices;

    const ListDigraph& graph = superGraph->graph;

    ListDigraph::NodeMap<bool> nodeFilter(graph); //default false
    nodeFilter[graph.nodeFromId(origin)] = true;
    nodeFilter[graph.nodeFromId(destination)] = true;
    for(const int& hubIndex : hubIndices){
        nodeFilter[graph.nodeFromId(hubIndex)] = true;
    }

    ListDigraph::ArcMap<bool> arcFilter(graph); //default false

    //Travel between hubs
    for(int i : hubIndices){
        for(int j : hubIndices) {
            addArcModeToFilter(i, j, BUS, arcFilter);
            addArcModeToFilter(i, j, RAIL, arcFilter);
        }
    }

    //Shuttle arcs
    addArcModeToFilter(origin, destination, SHUTTLE, arcFilter);
    for(int i : hubIndices){
        addArcModeToFilter(origin, i, SHUTTLE, arcFilter);
        addArcModeToFilter(i, destination, SHUTTLE, arcFilter);
    }

    //Remove arcs into origin and out of destination
    for(int i : hubIndices) {
        for(const travel_mode& mode : std::vector<travel_mode>{SHUTTLE, BUS, RAIL}){
            removeArcModeFromFilter(i, origin, mode, arcFilter);
            removeArcModeFromFilter(destination, i, mode, arcFilter);
        }
    }

    ListDigraph_VNode sourceNodes = {graph.nodeFromId(origin)};
    ListDigraph_VNode sinkNodes = {graph.nodeFromId(destination)};

    int maxCardinality = in.maximumNumberOfTransfers + 1;
        //in.maximumNumberOfTransfers = -1 is handled correctly
        //because maxCardinality = 0 corresponds to using a non-layered graph

    pNetwork network = nullptr;
    if(subproblemType == GUROBINETWORK){
#ifdef GUROBI
         network = std::make_shared<GurobiNetwork>(superGraph, nodeFilter, arcFilter,
                                                           sourceNodes, Vdouble{1}, sinkNodes, Vdouble{-1},
                                                           maxCardinality, subProblemEnv, resourceBounds);
#endif
    }
    else if(subproblemType == LEMONNETWORK){

        int costScaleFactor = 1E6;
        int capacityScaleFactor = 1E6;

        network = std::make_shared<LemonNetwork>(superGraph, nodeFilter, arcFilter,
                                                          sourceNodes, Vdouble{1}, sinkNodes, Vdouble{-1},
                                                          maxCardinality, costScaleFactor, capacityScaleFactor);
    }
    else{
        throw std::runtime_error("Exception: no subproblem network created.");
    }

    return network;

}

void NetworkDesignProblem::addArcModeToFilter(const int i, const int j, const travel_mode& mode, ListDigraph::ArcMap<bool>& arcFilter){

    const ListDigraph::Node& iNode = superGraph->graph.nodeFromId(i);
    const ListDigraph::Node& jNode = superGraph->graph.nodeFromId(j);

    //Match all arcs that correspond to different nb of vehicles
    ListDigraph_VArc matchingArcs = superGraph->getMatchingArcs(iNode, jNode, {(int) mode, 0}, {true, false});

    for(const ListDigraph::Arc& arc : matchingArcs){
        arcFilter[arc] = true;
    }

}

void NetworkDesignProblem::removeArcModeFromFilter(const int i, const int j, const travel_mode& mode, ListDigraph::ArcMap<bool>& arcFilter){

    const ListDigraph::Node& iNode = superGraph->graph.nodeFromId(i);
    const ListDigraph::Node& jNode = superGraph->graph.nodeFromId(j);

    //Match all arcs that correspond to different nb of vehicles
    ListDigraph_VArc matchingArcs = superGraph->getMatchingArcs(iNode, jNode, {(int) mode, 0}, {true, false});

    for(const ListDigraph::Arc& arc : matchingArcs){
        arcFilter[arc] = false;
    }

}

int NetworkDesignProblem::solve(const double eps){

    cplex->solve();

    const Vdouble& doubleSolution = getSolution();

    Vint intSolution;
    intSolution.reserve(doubleSolution.size());
    for(int i = 0; i < doubleSolution.size(); ++i){

        double doubleValue = doubleSolution.at(i);
        int intValue = (int) std::round(doubleValue);

        if(std::fabs(doubleValue - intValue) > eps){
            std::cerr << "Warning: solution rounded to integer (" << doubleValue << " != " << intValue << ", eps = " << eps << ")." << std::endl;
        }

        intSolution.push_back(intValue);

    }

    callback->reoptimizeSubproblems(intSolution, true); //reoptimize all subproblems for correct subproblem flows

    return cplex->getStatus();

}

std::pair<Vdouble, VAint5> NetworkDesignProblem::getSolutionAndResourceMap() const {
    return std::make_pair(getSolution(), resourceMap);
}

void NetworkDesignProblem::printSolution(const Stops& stops) const{

    double eps = 0.000001;

    for(int r = 0; r < nbResources; ++r){

        const Aint5& resourceInfo = resourceMap.at(r);
        const int& iIndex = resourceInfo.at(0);
        const int& iID = std::get<0>(stops.allStops.at(iIndex));
        const int& jIndex = resourceInfo.at(1);
        const int& jID = std::get<0>(stops.allStops.at(jIndex));
        const travel_mode& mode = (travel_mode) resourceInfo.at(2);
        const int& nbVehicles = resourceInfo.at(3);
        double val = cplex->getValue(resourceVariables.at(r));

        if(val > eps) {
            std::cout << "z_" << std::to_string(r) << " = " << std::to_string(val) << ": " <<
                         "from " << iID << " (index = " << iIndex << ") " <<
                         "to " << jID << " (index = " << jIndex << "), " <<
                         "mode = " << mode << ", nbVehicles = " << nbVehicles << std::endl;
        }

    }

}

std::tuple<double, Vdouble, double, Vdouble> NetworkDesignProblem::getNetworkScores(const Stops& stops, const RailLines& railLines) const {

    const std::pair<Vdouble, VAint5>& solution = getSolutionAndResourceMap();
    const Vdouble& values = solution.first;

    double designObjective = 0;
    double designCost = 0;
    double designConvenience = 0;

    const Vint& railStationIndices = in.getRailStationIndices(stops, railLines);

    for(int i : railStationIndices){
        for(int j : railStationIndices){
            if(in.arcModeExists(i, j, RAIL)) {

                const Vint& nbPossible = in.nbPossible(i, j, RAIL);
                if (nbPossible.size() != 1) {
                    throw std::runtime_error("Exception: multiple arcs with same origin, destination, and mode.");
                }
                int nb = nbPossible.at(0);

                const ListDigraph::Node from = superGraph->graph.nodeFromId(i);
                const ListDigraph::Node to = superGraph->graph.nodeFromId(j);
                ListDigraph::Arc superArc = superGraph->getUniqueMatchingArc(from, to, {(int) RAIL, nb});

                designObjective += in.arcObjective(i, j, RAIL, nb, Instance::DESIGN);
                designCost += in.arcCost(i, j, RAIL, nb, Instance::DESIGN);
                designConvenience += in.arcConvenience(i, j, RAIL, nb, Instance::DESIGN);

            }
        }
    }

    for(int r = 0; r < nbResources; ++r){

        Aint5 resource = resourceMap.at(r);
        const int& i = resource.at(0);
        const int& j = resource.at(1);
        const int& modeInt = resource.at(2);
        const int& nb = resource.at(3);

        const double& value = values.at(r);

        if(modeInt != (int) BUS){
            throw std::runtime_error("Exception: mode == BUS expected. Verify that rail is not counted twice (here and in railLines loop).");
        }

        if( (!isDummyResource(r) && std::round(value) == 1) || //selected variable
            (isDummyResource(r) && i != -1 && j != -1 && modeInt != -1 && nb != -1) ){ //meaningful dummy (like rail)
            designObjective += in.arcObjective(i, j, (travel_mode) modeInt, nb, Instance::DESIGN);
            designCost += in.arcCost(i, j, (travel_mode) modeInt, nb, Instance::DESIGN);
            designConvenience += in.arcConvenience(i, j, (travel_mode) modeInt, nb, Instance::DESIGN);
        }

    }

    double passengerObjective = 0;
    double passengerCost = 0;
    double passengerConvenience = 0;

    const ListDigraph& graph = superGraph->graph;

    for(int p = 0; p < subProblems.size(); ++p) {

        const pNetwork& subProblem = subProblems.at(p);
        const ListDigraph_VArc& path = subProblem->getArcFlowPath();
        const int& nbPeople = nbPeoples.at(p);

        for (int k = 0; k < path.size(); ++k) {

            const ListDigraph::Arc& superArc = path.at(k);
            const int& i = graph.id(graph.source(superArc));
            const int& j = graph.id(graph.target(superArc));
            const Vint& properties = superGraph->properties[superArc];
            const travel_mode& mode = (travel_mode) properties.at(0);
            const int& nb = properties.at(1);

            passengerObjective += nbPeople * in.arcObjective(i, j, mode, nb, Instance::PASSENGER);
            passengerCost += nbPeople * in.arcCost(i, j, mode, nb, Instance::PASSENGER);
            passengerConvenience += nbPeople * in.arcConvenience(i, j, mode, nb, Instance::PASSENGER);

        }

    }

    return std::make_tuple(designObjective, Vdouble{designCost, designConvenience},
                           passengerObjective, Vdouble{passengerCost, passengerConvenience});

};

Vdouble NetworkDesignProblem::getSolution() const {

    Vdouble solution;
    solution.reserve(nbResources);

    for(int r = 0; r < nbResources; ++r) {
        solution.push_back(cplex->getValue(resourceVariables.at(r)));
    }

    return solution;

}

void NetworkDesignProblem::removeSuboptimalBusesFromSuperGraph(){

    ListDigraph& graph = superGraph->graph;

    const Vdouble& solution = getSolution();

    for(int r = 0; r < nbResources; ++r){

        const Aint5& resourceInfo = resourceMap.at(r);
        const int& i = resourceInfo.at(0);
        const int& j = resourceInfo.at(1);
        const travel_mode& mode = (travel_mode) resourceInfo.at(2);
        const int& nbVehicles = resourceInfo.at(3);
        const bool& dummy = resourceInfo.at(4);

        if(!dummy && mode == BUS){
            if(solution.at(r) < 0.5){

                const ListDigraph::Node from = graph.nodeFromId(i);
                const ListDigraph::Node to = graph.nodeFromId(j);
                Vint properties = {(int) mode, nbVehicles};

                ListDigraph::Arc superArc = superGraph->getUniqueMatchingArc(from, to, properties);
                graph.erase(superArc);

            }
        }

    }

}

ListDigraph_VVArc NetworkDesignProblem::generateAllPaths(){

    if(!in.keepAllShuttleArcs){
        throw std::runtime_error("Exception: can only call generateAllPaths() when instance is generated with keepAllShuttleArcs = true.");
    }

    removeSuboptimalBusesFromSuperGraph();

    const ListDigraph& graph = superGraph->graph;

    ListDigraph_VNode sourceNodes = ListDigraph_VNode();
    ListDigraph_VNode sinkNodes = ListDigraph_VNode();
    for(ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node){
        sourceNodes.push_back(node);
        sinkNodes.push_back(node);
    }

    ListDigraph::NodeMap<bool> nodeFilter(graph, true);
    ListDigraph::ArcMap<bool> arcFilter(graph, true);

    int maxCardinality = in.maximumNumberOfTransfers + 1;
    if(maxCardinality <= 0){
        throw std::runtime_error("Exception: generateAllPaths() not implemented for the case without layers.");
    }

    int costScaleFactor = 1;
    int capacityScaleFactor = 1;

    std::shared_ptr<LemonNetwork> network =
            std::make_shared<LemonNetwork>(superGraph, nodeFilter, arcFilter,
                                           sourceNodes, Vdouble(sourceNodes.size(), 0),
                                           sinkNodes, Vdouble(sinkNodes.size(), 0),
                                           maxCardinality, costScaleFactor, capacityScaleFactor);

    intSet superHubIDs = intSet(in.hubIndices.begin(), in.hubIndices.end());
    return network->getAllShortestPathsShuttleFiltered(superHubIDs);

}

void NetworkDesignProblem::verifyResourceMapUnique() {

    int nbRes = 0;

    Aint5USet resourceSet = Aint5USet();

    for(int r = 0; r < resourceMap.size(); ++r){

        const Aint5& resource = resourceMap.at(r);

        if(!isDummyResource(resource)){
            ++nbRes;
            resourceSet.insert(resource);
            if(resourceSet.size() != nbRes){
                throw std::runtime_error("Exception: resource with index " + std::to_string(r) + " is a duplicate.");
            }
        }

    }

}

void NetworkDesignProblem::setSuperGraphContributions() {

    const ListDigraph& graph = superGraph->graph;
    auto& contributions = superGraph->contributions;

    for(int r = 0; r < nbResources; ++r){

        const Aint5& resourceInfo = resourceMap.at(r);
        const int& i = resourceInfo.at(0);
        const int& j = resourceInfo.at(1);
        const travel_mode& mode = (travel_mode) resourceInfo.at(2);
        const int& nbVehicles = resourceInfo.at(3);
        const bool& dummy = resourceInfo.at(4);

        if(dummy){
            continue;
        }

        const ListDigraph::Node from = graph.nodeFromId(i);
        const ListDigraph::Node to = graph.nodeFromId(j);
        Vint properties = {(int) mode, nbVehicles};

        const ListDigraph::Arc& arc = superGraph->getUniqueMatchingArc(from, to, properties);
        contributions[arc].push_back({r, 1});

    }

}

void NetworkDesignProblem::initializeResourceBounds() {

    //All resources (dummy or not) get initial bound (-inf, 0)
    resourceBounds = VResourceBounds(nbResources, {kevintools::NEG_INF_DOUBLE, 0});

}

bool NetworkDesignProblem::isDummyResource(const int index) const {

    if(index < 0 || index >= nbResources){
        throw std::runtime_error("Exception: invalid resource index.");
    }

    return isDummyResource(resourceMap.at(index));

}

bool NetworkDesignProblem::isDummyResource(const Aint5& resource) const {
    return (bool) resource.at(4);
}


