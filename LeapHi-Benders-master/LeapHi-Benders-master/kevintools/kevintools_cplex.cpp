#include "kevintools_cplex.hpp"
#include "kevintools_primitives.hpp"

bool kevintools::isEmpty(const IloExpr& expression){
    return !expression.getLinearIterator().ok();
}

void kevintools::setName(const IloNumVar& var, const std::string& name){
    var.setName(name.c_str());
}

IloRange kevintools::addExpressionEqConstraint(IloModel &model, const IloExpr &expression, const double rhs) {

    if(isEmpty(expression)){
        return {}; //do not add empty constraints
    }

    if(kevintools::is_posinf(rhs) || kevintools::is_neginf(rhs)){
        throw std::runtime_error("Exception: lower bound and/or upper bound invalid.");
    }

    IloRange constraint(model.getEnv(), rhs, expression, rhs);
    model.add(constraint);

    return constraint;

}

IloRange kevintools::addExpressionEqConstraint(const pIloModel& model, const IloExpr &expression, const double rhs) {
    return addExpressionEqConstraint(*model, expression, rhs);
}

IloRangeArray kevintools::addExpressionEqConstraints(IloModel& model, const VIloExpr& expressions, const Vdouble& rhs){

    if(expressions.size() != rhs.size()){
        throw std::runtime_error("Exception: expressions.size() != rhs.size().");
    }

    const int& nbExpressions = expressions.size();
    IloRangeArray constraints = IloRangeArray(model.getEnv(), nbExpressions);

    for(int i = 0; i < nbExpressions; ++i) {
        constraints[i] = kevintools::addExpressionEqConstraint(model, expressions.at(i), rhs.at(i));
    }

    return constraints;

}

IloRangeArray kevintools::addExpressionEqConstraints(const pIloModel& model, const VIloExpr &expressions, const Vdouble &rhs){
    return addExpressionEqConstraints(*model, expressions, rhs);
}

IloRangeArray kevintools::addExpressionEqConstraints(IloModel& model, const IloExprArray &expressions, const Vdouble &rhs){

    if(expressions.getSize() != rhs.size()){
        throw std::runtime_error("Exception: expressions.size() != rhs.size().");
    }

    const int& nbExpressions = expressions.getSize();
    IloRangeArray constraints = IloRangeArray(model.getEnv(), nbExpressions);

    for(int i = 0; i < nbExpressions; ++i) {
        if(expressions[i].isValid()) {
            constraints[i] = kevintools::addExpressionEqConstraint(model, expressions[i], rhs.at(i));
        }
    }

    return constraints;

}

kevintools::SenseConstraintMap_CPX kevintools::addExpressionBoundConstraint(IloModel& model, const IloExpr& expression, const double lb, const double ub){

    //TODO: check for IloInfinity instead of numerical_limits<double>::max()

    if(isEmpty(expression)){
        return SenseConstraintMap_CPX(); //do not add empty constraints
    }

    SenseConstraintMap_CPX result;

    if(lb > ub || kevintools::is_posinf(lb) || kevintools::is_neginf(ub) || (kevintools::is_inf(lb) && kevintools::is_inf(ub))){
        throw std::runtime_error("Exception: lower bound and/or upper bound invalid.");
    }

    if(lb == ub){
        constraint_sense sense = EQ;
        IloRange constraint(model.getEnv(), lb, expression, ub);
        model.add(constraint);
        result.insert(std::make_pair(sense, constraint));
    }
    else{
        if(!kevintools::is_posinf(ub)){
            constraint_sense sense = LE;
            IloRange constraint(model.getEnv(), -IloInfinity, expression, ub);
            model.add(constraint);
            result.insert(std::make_pair(sense, constraint));
        }
        if(!kevintools::is_neginf(lb)){
            constraint_sense sense = GE;
            IloRange constraint(model.getEnv(), lb, expression, IloInfinity);
            model.add(constraint);
            result.insert(std::make_pair(sense, constraint));
        }
    }

    return result;

}

kevintools::SenseConstraintMap_CPX kevintools::addExpressionBoundConstraint(const pIloModel& model, const IloExpr& expression, const double lb, const double ub){
    return addExpressionBoundConstraint(*model, expression, lb, ub);
}

kevintools::SenseConstraintMap_CPX kevintools::addExpressionBoundConstraint(IloModel& model, const IloExpr& expression, const std::array<double, 2>& bounds){
    return kevintools::addExpressionBoundConstraint(model, expression, bounds.at(0), bounds.at(1));
}

kevintools::SenseConstraintMap_CPX kevintools::addExpressionBoundConstraint(const pIloModel& model, const IloExpr &expression, const std::array<double, 2>& bounds){
    return kevintools::addExpressionBoundConstraint(*model, expression, bounds.at(0), bounds.at(1));
}

kevintools::VSenseConstraintMap_CPX kevintools::addExpressionBoundConstraints(IloModel& model, const VIloExpr &expressions, const std::vector<std::array<double, 2>> &bounds){

    if(expressions.size() != bounds.size()){
        throw std::runtime_error("Exception: expressions.size() != bounds.size().");
    }

    const int& nbExpressions = expressions.size();
    VSenseConstraintMap_CPX constraints;
    constraints.reserve(nbExpressions);

    for(int i = 0; i < nbExpressions; ++i) {
        constraints.push_back(addExpressionBoundConstraint(model, expressions.at(i), bounds.at(i)));
    }

    return constraints;

}

kevintools::VSenseConstraintMap_CPX kevintools::addExpressionBoundConstraints(const pIloModel& model, const VIloExpr &expressions, const std::vector<std::array<double, 2>> &bounds){
    return addExpressionBoundConstraints(*model, expressions, bounds);
}

kevintools::VSenseConstraintMap_CPX kevintools::addExpressionBoundConstraints(IloModel& model, const IloExprArray &expressions, const std::vector<std::array<double, 2>> &bounds){

    if(expressions.getSize() != bounds.size()){
        throw std::runtime_error("Exception: expressions.size() != bounds.size().");
    }

    const int& nbExpressions = expressions.getSize();
    VSenseConstraintMap_CPX constraints;
    constraints.reserve(nbExpressions);

    for(int i = 0; i < nbExpressions; ++i) {
        constraints.push_back(addExpressionBoundConstraint(model, expressions[i], bounds.at(i)));
    }

    return constraints;

}

void kevintools::clearExpressions(IloExprArray &expressions) {
    expressions.endElements();
}

void kevintools::clearExpressions(kevintools::VIloExprArray &expressions) {
    for(IloExprArray& element : expressions){
        clearExpressions(element);
    }
    expressions.end();
}

void kevintools::clearExpressions(VIloExpr& expressions) {
    for(IloExpr& expression : expressions){
        expression.end();
    }
}

void kevintools::clearExpressions(VVIloExpr& expressions) {
    for(VIloExpr& element : expressions){
        clearExpressions(element);
    }
}

void kevintools::initializeIfNeeded(IloNumExpr& expression, IloEnv env){

    if(!expression.isValid()){
        expression = IloNumExpr(env);
    }

}

IloExpr kevintools::sum(const VIloExpr& expressions, IloEnv env){

    IloExpr sum(env);
    for(const IloExpr& expression : expressions){
        sum += expression;
    }

    return sum;

}