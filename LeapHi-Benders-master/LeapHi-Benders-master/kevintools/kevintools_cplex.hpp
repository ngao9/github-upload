#pragma once

#include "kevintools_typedefs.hpp"
#include "ilcplex/ilocplex.h"

#include <array>

namespace kevintools{

    //CPLEX typdefs
    typedef std::vector<IloNumVar> VIloNumVar;
    typedef std::vector<VIloNumVar> VVIloNumVar;
    typedef std::vector<VVIloNumVar> VVVIloNumVar;
    typedef std::shared_ptr<IloEnv> pIloEnv;
    typedef std::shared_ptr<IloModel> pIloModel;
    typedef std::shared_ptr<IloCplex> pIloCplex;
    typedef std::vector<IloIntExpr> VIloIntExpr;
    typedef std::vector<VIloIntExpr> VVIloIntExpr;
    typedef std::vector<IloExpr> VIloExpr;
    typedef std::vector<VIloExpr> VVIloExpr;
    typedef std::vector<IloExprArray> VIloExprArray;
    typedef std::vector<IloRange> VIloRange;
    typedef std::vector<VIloRange> VVIloRange;
    typedef std::vector<IloRangeArray> VIloRangeArray;

    //SenseConstraintMap_CPX
    using namespace constraint_sense_wrapper;
    typedef std::map<constraint_sense, IloRange> SenseConstraintMap_CPX;
    typedef std::vector<SenseConstraintMap_CPX> VSenseConstraintMap_CPX;
    typedef std::vector<VSenseConstraintMap_CPX> VVSenseConstraintMap_CPX;

    //Methods
    bool isEmpty(const IloExpr& expression);

    void setName(const IloNumVar& var, const std::string& name);

    IloRange addExpressionEqConstraint(IloModel& model, const IloExpr& expression, double rhs);
    IloRange addExpressionEqConstraint(const pIloModel& model, const IloExpr& expression, double rhs);

    IloRangeArray addExpressionEqConstraints(IloModel& model, const VIloExpr &expressions, const Vdouble &rhs);
    IloRangeArray addExpressionEqConstraints(const pIloModel& model, const VIloExpr &expressions, const Vdouble &rhs);
    IloRangeArray addExpressionEqConstraints(IloModel& model, const IloExprArray &expressions, const Vdouble &rhs);

    SenseConstraintMap_CPX addExpressionBoundConstraint(IloModel& model, const IloExpr &expression, double lb, double ub);
    SenseConstraintMap_CPX addExpressionBoundConstraint(const pIloModel& model, const IloExpr &expression, double lb, double ub);
    SenseConstraintMap_CPX addExpressionBoundConstraint(IloModel& model, const IloExpr &expression, const std::array<double, 2>& bounds);
    SenseConstraintMap_CPX addExpressionBoundConstraint(const pIloModel& model, const IloExpr &expression, const std::array<double, 2>& bounds);

    VSenseConstraintMap_CPX addExpressionBoundConstraints(IloModel& model, const VIloExpr &expressions, const std::vector<std::array<double, 2>> &bounds);
    VSenseConstraintMap_CPX addExpressionBoundConstraints(const pIloModel& model, const VIloExpr &expressions, const std::vector<std::array<double, 2>> &bounds);
    VSenseConstraintMap_CPX addExpressionBoundConstraints(IloModel& model, const IloExprArray &expressions, const std::vector<std::array<double, 2>> &bounds);

    void clearExpressions(IloExprArray& expressions);
    void clearExpressions(VIloExprArray& expressions);

    void clearExpressions(VIloExpr& expressions);
    void clearExpressions(VVIloExpr& expressions);

    void initializeIfNeeded(IloNumExpr& expression, IloEnv env);

    IloExpr sum(const VIloExpr& expressions, IloEnv env);

}