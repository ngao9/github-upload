#pragma once

#include "kevintools_typedefs.hpp"
#include <array>

#ifdef GUROBI

    #include "gurobi_c++.h"

namespace kevintools{

    //Gurobi typdefs
    typedef std::vector<GRBVar> VGRBVar;
    typedef std::vector<VGRBVar> VVGRBVar;
    typedef std::vector<VVGRBVar> VVVGRBVar;
    typedef std::shared_ptr<GRBEnv> pGRBEnv;
    typedef std::shared_ptr<GRBModel> pGRBModel;
    typedef std::vector<GRBConstr> VGRBConstr;
    typedef std::vector<VGRBConstr> VVGRBConstr;
    typedef std::vector<GRBLinExpr> VGRBLinExpr;

    //SenseConstraintMap_GRB
    using namespace constraint_sense_wrapper;
    typedef std::map<constraint_sense, GRBConstr> SenseConstraintMap_GRB;
    typedef std::vector<SenseConstraintMap_GRB> VSenseConstraintMap_GRB;
    typedef std::vector<VSenseConstraintMap_GRB> VVSenseConstraintMap_GRB;

    //Methods
    GRBConstr addExpressionEqConstraint(GRBModel& model, const GRBLinExpr& expression, double rhs);
    GRBConstr addExpressionEqConstraint(const pGRBModel& model, const GRBLinExpr& expression, double rhs);

    VGRBConstr addExpressionEqConstraints(GRBModel& model, const VGRBLinExpr &expressions, const Vdouble &rhs);
    VGRBConstr addExpressionEqConstraints(const pGRBModel& model, const VGRBLinExpr &expressions, const Vdouble &rhs);

    SenseConstraintMap_GRB addExpressionBoundConstraint(GRBModel& model, const GRBLinExpr &expression, double lb, double ub);
    SenseConstraintMap_GRB addExpressionBoundConstraint(const pGRBModel& model, const GRBLinExpr &expression, double lb, double ub);
    SenseConstraintMap_GRB addExpressionBoundConstraint(GRBModel& model, const GRBLinExpr &expression, const std::array<double, 2>& bounds);
    SenseConstraintMap_GRB addExpressionBoundConstraint(const pGRBModel& model, const GRBLinExpr &expression, const std::array<double, 2>& bounds);

    VSenseConstraintMap_GRB addExpressionBoundConstraints(GRBModel& model, const VGRBLinExpr &expressions, const std::vector<std::array<double, 2>> &bounds);
    VSenseConstraintMap_GRB addExpressionBoundConstraints(const pGRBModel& model, const VGRBLinExpr &expressions, const std::vector<std::array<double, 2>> &bounds);

}

#endif