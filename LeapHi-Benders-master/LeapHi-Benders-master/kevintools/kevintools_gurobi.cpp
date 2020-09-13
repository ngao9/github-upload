#include "kevintools_gurobi.hpp"
#include "kevintools_primitives.hpp"

#ifdef GUROBI

GRBConstr kevintools::addExpressionEqConstraint(GRBModel &model, const GRBLinExpr &expression, const double rhs) {

    if(expression.size() == 0){
        return {}; //do not add empty constraints
    }

    if(kevintools::is_posinf(rhs) || kevintools::is_neginf(rhs)){
        throw std::runtime_error("Exception: lower bound and/or upper bound invalid.");
    }

    return model.addConstr(expression == rhs);

}

GRBConstr kevintools::addExpressionEqConstraint(const pGRBModel& model, const GRBLinExpr &expression, const double rhs) {
    return addExpressionEqConstraint(*model, expression, rhs);
}

kevintools::VGRBConstr kevintools::addExpressionEqConstraints(GRBModel& model, const VGRBLinExpr& expressions, const Vdouble& rhs){

    if(expressions.size() != rhs.size()){
        throw std::runtime_error("Exception: expressions.size() != rhs.size().");
    }

    const int& nbExpressions = expressions.size();
    VGRBConstr constraints;
    constraints.reserve(nbExpressions);

    for(int i = 0; i < nbExpressions; ++i) {
        constraints.push_back(kevintools::addExpressionEqConstraint(model, expressions.at(i), rhs.at(i)));
    }

    return constraints;

}

kevintools::VGRBConstr kevintools::addExpressionEqConstraints(const pGRBModel& model, const VGRBLinExpr &expressions, const Vdouble &rhs){
    return addExpressionEqConstraints(*model, expressions, rhs);
}

kevintools::SenseConstraintMap_GRB kevintools::addExpressionBoundConstraint(GRBModel& model, const GRBLinExpr& expression, const double lb, const double ub){

    if(expression.size() == 0){
        return SenseConstraintMap_GRB(); //do not add empty constraints
    }

    SenseConstraintMap_GRB result;

    if(lb > ub || kevintools::is_posinf(lb) || kevintools::is_neginf(ub) || (kevintools::is_inf(lb) && kevintools::is_inf(ub))){
        throw std::runtime_error("Exception: lower bound and/or upper bound invalid.");
    }

    if(lb == ub){
        constraint_sense sense = EQ;
        GRBConstr constraint = model.addConstr(expression == lb);
        result.insert(std::make_pair(sense, constraint));
    }
    else{
        if(!kevintools::is_posinf(ub)){
            constraint_sense sense = LE;
            GRBConstr constraint = model.addConstr(expression <= ub);
            result.insert(std::make_pair(sense, constraint));
        }
        if(!kevintools::is_neginf(lb)){
            constraint_sense sense = GE;
            GRBConstr constraint = model.addConstr(expression >= lb);
            result.insert(std::make_pair(sense, constraint));
        }
    }

    return result;

}

kevintools::SenseConstraintMap_GRB kevintools::addExpressionBoundConstraint(const pGRBModel& model, const GRBLinExpr& expression, const double lb, const double ub){
    return addExpressionBoundConstraint(*model, expression, lb, ub);
}

kevintools::SenseConstraintMap_GRB kevintools::addExpressionBoundConstraint(GRBModel& model, const GRBLinExpr& expression, const std::array<double, 2>& bounds){
    return kevintools::addExpressionBoundConstraint(model, expression, bounds.at(0), bounds.at(1));
}

kevintools::SenseConstraintMap_GRB kevintools::addExpressionBoundConstraint(const pGRBModel& model, const GRBLinExpr &expression, const std::array<double, 2>& bounds){
    return kevintools::addExpressionBoundConstraint(*model, expression, bounds.at(0), bounds.at(1));
}

kevintools::VSenseConstraintMap_GRB kevintools::addExpressionBoundConstraints(GRBModel& model, const VGRBLinExpr &expressions, const std::vector<std::array<double, 2>> &bounds){

    if(expressions.size() != bounds.size()){
        throw std::runtime_error("Exception: expressions.size() != bounds.size().");
    }

    const int& nbExpressions = expressions.size();
    VSenseConstraintMap_GRB constraints;
    constraints.reserve(nbExpressions);

    for(int i = 0; i < nbExpressions; ++i) {
        constraints.push_back(addExpressionBoundConstraint(model, expressions.at(i), bounds.at(i)));
    }

    return constraints;

}

kevintools::VSenseConstraintMap_GRB kevintools::addExpressionBoundConstraints(const pGRBModel& model, const VGRBLinExpr &expressions, const std::vector<std::array<double, 2>> &bounds){
    return addExpressionBoundConstraints(*model, expressions, bounds);
}

#endif