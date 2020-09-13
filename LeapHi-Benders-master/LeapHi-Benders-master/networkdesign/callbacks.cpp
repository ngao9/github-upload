#include "callbacks.hpp"

#include "networksolver_cpp_include.hpp"

#include "networkdesignproblem.hpp"
#include "instance.hpp"

#include "kevintools_vectors.hpp"
#include "kevintools_primitives.hpp"

using namespace networkdesign;

Callbacks::Callbacks(NetworkDesignProblem& problem) :
        inst(NULL), problem(problem), corePoint(), feasibleSolutions()
{
}

void Callbacks::initializeCorePoint(){

    const int& nbResources = problem.nbResources;
    const VAint5& resourceMap = problem.resourceMap;
    const VIloIntVar& resourceVariables = problem.resourceVariables;
    const pIloModel& model = problem.model;

    corePoint = Vdouble(nbResources, 0);

    for(int r = 0; r < nbResources; ++r){

        const Aint5& resourceInfo = resourceMap.at(r);
        int i = resourceInfo.at(0);
        int j = resourceInfo.at(1);

        if(problem.isDummyResource(resourceInfo)){

            if(problem.resourceBounds.size() != nbResources){
                throw std::runtime_error("Exception: cannot set core point value for dummy resource if problem.resourceBounds is not set.");
            }

            double ub = problem.resourceBounds.at(r).at(1);
            if(kevintools::is_inf(ub)){
                throw std::runtime_error("Exception: cannot set core point value for dummy resource with infinite upperbound.");
            }

            corePoint.at(r) = ub;

        }
        else{

            const travel_mode& mode = (travel_mode) resourceInfo.at(2);

            if(mode == BUS || mode == RAIL){
                corePoint.at(r) = 0.5 / ((double) problem.in.nbPossible(i, j, mode).size());
            }
            else{
                std::runtime_error("Exception: unexpected mode encountered while determining corePoint.");
            }

        }

    }

}

void Callbacks::userCutCallback(IloCplex::CallbackI* instance) {

    this->inst = instance;

    MyUserCutCallbackI* uccInst = dynamic_cast<MyUserCutCallbackI*>(inst);
    if(!uccInst){
        throw std::runtime_error("Exception: downcasting to MyUserCutCallbackI failed.");
    }

    const int& nbVars = problem.nbResources;

    const Solution& values = getCurrentValues();
    const double& thetaVal = values.first;
    const Vdouble& z = values.second;

    const int& nbVariables = problem.nbResources;

    Vdouble point = Vdouble(nbVariables);
    for(int r = 0; r < nbVariables; ++r){
        point.at(r) = (1 - 2 * EPS_STABALIZE) * z.at(r) + 2 * EPS_STABALIZE * corePoint.at(r);
    }

    reoptimizeSubproblems(point);

    const IloNumVar& theta = problem.theta;

    VIloExpr bendersExpressions = getBendersExpressions(point);
    IloExpr aggregatedExpression = kevintools::sum(bendersExpressions, inst->getEnv());
    kevintools::clearExpressions(bendersExpressions);

    double cutVal = uccInst->getValue(aggregatedExpression);

    if(cutVal - thetaVal > EPS_USERCUT_VIOLATION) {
        uccInst->add(theta >= aggregatedExpression);
    }
    aggregatedExpression.end();

}

void Callbacks::lazyConstraintCallback(IloCplex::CallbackI* instance) {

    this->inst = instance;

    MyLazyConstraintCallbackI* lccInst = dynamic_cast<MyLazyConstraintCallbackI*>(inst);
    if(!lccInst){
        throw std::runtime_error("Exception: downcasting to MyLazyConstraintCallbackI failed.");
    }

    const int& nbVars = problem.nbResources;

    const IntSolution& values = getCurrentIntValues();
    const double& thetaVal = values.first;
    const Vint& z = values.second;

    bool requireTightBound = true;
    reoptimizeSubproblems(z, requireTightBound);

    const IloNumVar& theta = problem.theta;

    VIloExpr bendersExpressions = getBendersExpressions(kevintools::castToDouble(z));
    IloExpr aggregatedExpression = kevintools::sum(bendersExpressions, inst->getEnv());
    kevintools::clearExpressions(bendersExpressions);

    double cutVal = lccInst->getValue(aggregatedExpression);

    if(cutVal - thetaVal > EPS_LAZYCONSTRAINT_VIOLATION) {
        lccInst->add(theta >= aggregatedExpression);
    }
    aggregatedExpression.end();

    IntSolution newFeasibleSolution = std::make_pair(cutVal, z);
    feasibleSolutions.push(newFeasibleSolution);

}

void Callbacks::heuristicCallback(IloCplex::CallbackI* instance) {

    this->inst = instance;

    MyHeuristicCallbackI* hcInst = dynamic_cast<MyHeuristicCallbackI*>(inst);
    if(!hcInst){
        throw std::runtime_error("Exception: downcasting to MyHeuristicCallbackI failed.");
    }

    if(feasibleSolutions.size() == 0){
        return;
    }

    //pass first feasible solution in the list to CPLEX
    //cannot set multiple in the same iteration

    const IntSolution& feasibleSolution = feasibleSolutions.front();

    const VIloIntVar& variables = problem.resourceVariables;
    const int& nbResources = problem.nbResources;
    const IloNumVar& theta = problem.theta;

    IloNumVarArray allVariables(problem.env, nbResources + 1);
    IloNumArray allValues(problem.env, nbResources + 1);

    for(int i = 0; i < nbResources; ++i){
        allVariables[i] = variables.at(i);
        allValues[i] = feasibleSolution.second.at(i);
    }
    allVariables[nbResources] = theta;
    allValues[nbResources] = feasibleSolution.first;

    hcInst->setSolution(allVariables, allValues);
    //Note that manually set solutions may violate lazy constraints: https://or.stackexchange.com/q/3327/145

    feasibleSolutions.pop();

}

Solution Callbacks::getCurrentValues(){

    IloCplex::ControlCallbackI* ccInst = dynamic_cast<IloCplex::ControlCallbackI*>(inst);
    if(!ccInst){
        throw std::runtime_error("Exception: downcasting to ControlCallbackI failed.");
    }

    const VIloIntVar& variables = problem.resourceVariables;
    const int& nbVariables = problem.nbResources;

    Vdouble values = Vdouble();
    values.reserve(nbVariables);

    for(int r = 0; r < nbVariables; ++r){
        values.push_back(ccInst->getValue(variables.at(r)));
    }

    double theta = ccInst->getValue(problem.theta);

    return std::make_pair(theta, values);

}

IntSolution Callbacks::getCurrentIntValues() {

    IloCplex::ControlCallbackI* ccInst = dynamic_cast<IloCplex::ControlCallbackI*>(inst);
    if(!ccInst){
        throw std::runtime_error("Exception: downcasting to ControlCallbackI failed.");
    }

    const VIloIntVar& variables = problem.resourceVariables;
    const int& nbVariables = problem.nbResources;

    Vint values = Vint();
    values.reserve(nbVariables);

    for(int r = 0; r < nbVariables; ++r){
        values.push_back((int) (ccInst->getValue(variables.at(r)) + 0.5));
    }

    double theta = ccInst->getValue(problem.theta);

    return std::make_pair(theta, values);

}

void Callbacks::reoptimizeSubproblems(const Vint& upperBounds, const bool requireTightBound) const{
    reoptimizeSubproblems(kevintools::castToDouble(upperBounds), requireTightBound);
}

void Callbacks::reoptimizeSubproblems(const Vdouble& upperBounds, const bool requireTightBound) const{

    const VpNetwork& networks = problem.subProblems;

    VResourceBounds& resourceBounds = problem.resourceBounds;
    const int& nbResources = problem.nbResources;

    for(int r = 0; r < nbResources; ++r){
        if(!problem.isDummyResource(r)){ //only update bounds if resource in the master problem
            resourceBounds.at(r).at(1) = upperBounds.at(r);
        }
    }

    #pragma omp parallel for default(none) shared(networks, resourceBounds)
    for(int i = 0; i < networks.size(); ++i){
        const pNetwork& network = networks.at(i);
        network->updateResourceBounds(resourceBounds, requireTightBound);
        if(!network->solveModel()){
            throw std::runtime_error("Exception: network could not be solved.");
        }
    }

}

VIloExpr Callbacks::getBendersExpressions(const Vdouble& currentValues) {

    //For given network:
    //nbPeople * (opt + sum_j dual_j (z_j - z*_j))

    const VIloIntVar& variables = problem.resourceVariables;
    const VpNetwork& networks = problem.subProblems;
    const int& nbNetworks = networks.size();

    VIloExpr expressions;
    expressions.reserve(nbNetworks);

    for(int i = 0; i < nbNetworks; ++i){

        const pNetwork& network = networks.at(i);
        const int& nbPeople = problem.nbPeoples.at(i);
        IloExpr expression(inst->getEnv());

        const double& obj = network->getObjValue();
        expression += nbPeople * obj;

        const Vdouble& networkDuals = network->getResourceConstraintDuals(LE);
        const Vint& resourceIndices = network->getResourceIndices();

        for(int j = 0; j < resourceIndices.size(); ++j) {
            const int& resource = resourceIndices.at(j);
            const double& dual = networkDuals.at(j);
            if (dual != 0) {
                expression += nbPeople * dual * variables.at(resource);
                expression -= nbPeople * dual * currentValues.at(resource);
            }
        }

        expressions.push_back(expression);

    }

    return expressions;

}


////////////////////////////////
// CPLEX CALLBACK DEFINITIONS //
////////////////////////////////

IloCplex::CallbackI* MyUserCutCallbackI::duplicateCallback() const {
   return (new (getEnv()) MyUserCutCallbackI(*this));
}

MyUserCutCallbackI::MyUserCutCallbackI(IloEnv env, const std::shared_ptr<Callbacks> callbacks) :
    IloCplex::UserCutCallbackI(env), callbacks(callbacks) {}

void MyUserCutCallbackI::main(){
    callbacks->userCutCallback(this);
}

IloCplex::CallbackI* MyLazyConstraintCallbackI::duplicateCallback() const {
    return (new (getEnv()) MyLazyConstraintCallbackI(*this));
}

MyLazyConstraintCallbackI::MyLazyConstraintCallbackI(IloEnv env, const std::shared_ptr<Callbacks> callbacks) :
        IloCplex::LazyConstraintCallbackI(env), callbacks(callbacks) {}

void MyLazyConstraintCallbackI::main(){
    callbacks->lazyConstraintCallback(this);
}

IloCplex::CallbackI* MyHeuristicCallbackI::duplicateCallback() const {
    return (new (getEnv()) MyHeuristicCallbackI(*this));
}

MyHeuristicCallbackI::MyHeuristicCallbackI(IloEnv env, const std::shared_ptr<Callbacks> callbacks) :
        IloCplex::HeuristicCallbackI(env), callbacks(callbacks) {}

void MyHeuristicCallbackI::main(){
    callbacks->heuristicCallback(this);
}