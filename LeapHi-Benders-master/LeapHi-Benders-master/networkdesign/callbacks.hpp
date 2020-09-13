#pragma once

#include "networkdesign_typedefs.hpp"
#include "kevintools_cplex.hpp"

namespace networkdesign {

    class Callbacks : public std::enable_shared_from_this<Callbacks>{


    //Initialization
    public:
        explicit Callbacks(NetworkDesignProblem& problem);
        void initializeCorePoint(); //can only be initialized after the master problem is built


    //References
    private:
        IloCplex::CallbackI* inst;      //pointer to the invoking callback instance
        NetworkDesignProblem& problem;  //reference to the master problem


    //CPLEX legacy callbacks
    public:
        void userCutCallback(IloCplex::CallbackI* instance);
        void lazyConstraintCallback(IloCplex::CallbackI* instance);
        void heuristicCallback(IloCplex::CallbackI* instance);
    private:
        QIntSolution feasibleSolutions; //queue of solutions found by lazyConstraint,
                                        //to be inserted by heuristicCallback


    //Access to solutions
    private:
        Solution getCurrentValues();
        IntSolution getCurrentIntValues();


    //Cut generation
    public:
        void reoptimizeSubproblems(const Vint& upperbounds, bool requireTightBound=false) const;
        void reoptimizeSubproblems(const Vdouble& upperBounds, bool requireTightBound=false) const;
    private:
        VIloExpr getBendersExpressions(const Vdouble& currentValues);


    //Constants
    private:
        static constexpr double EPS_USERCUT_VIOLATION = 1E-6;
        static constexpr double EPS_LAZYCONSTRAINT_VIOLATION = 1E-6;
        static constexpr double EPS_STABALIZE = 0.00001;
        static constexpr double EPS_PRINTING = 1E-6;


    //Core point
    private:
        Vdouble corePoint;

    };


    ////////////////////////////////
    // CPLEX CALLBACK DEFINITIONS //
    ////////////////////////////////

    //Manually define MyUserCutCallbackI based on ILOUSERCUTCALLBACK1 macro.
    class MyUserCutCallbackI : public IloCplex::UserCutCallbackI {
        const std::shared_ptr<Callbacks> callbacks;
        friend class Callbacks; //to allow access to protected methods
    public:
        IloCplex::CallbackI* duplicateCallback() const;
        MyUserCutCallbackI(IloEnv env, std::shared_ptr<Callbacks> callbacks);
        void main();
    };
    inline IloCplex::Callback MyUserCutCallback(IloEnv env, const std::shared_ptr<Callbacks> callbacks) {
        return (IloCplex::Callback(new (env) MyUserCutCallbackI(env, callbacks)));
    }

    //Manually define MyLazyConstraintCallbackI based on ILOLAZYCONSTRAINTCALLBACK1 macro.
    class MyLazyConstraintCallbackI : public IloCplex::LazyConstraintCallbackI {
        std::shared_ptr<Callbacks> callbacks;
        friend class Callbacks;
    public:
        IloCplex::CallbackI* duplicateCallback() const;
        MyLazyConstraintCallbackI(IloEnv env, std::shared_ptr<Callbacks> callbacks);
        void main();
    };
    inline IloCplex::Callback MyLazyConstraintCallback(IloEnv env, const std::shared_ptr<Callbacks> callbacks) {
        return (IloCplex::Callback(new (env) MyLazyConstraintCallbackI(env, callbacks)));
    }

    //Manually define MyHeuristicCallbackI based on ILOHEURISTICCALLBACK1 macro.
    class MyHeuristicCallbackI : public IloCplex::HeuristicCallbackI {
        std::shared_ptr<Callbacks> callbacks;
        friend class Callbacks;
    public:
        IloCplex::CallbackI* duplicateCallback() const;
        MyHeuristicCallbackI(IloEnv env, std::shared_ptr<Callbacks> callbacks);
        void main();
    };
    inline IloCplex::Callback MyHeuristicCallback(IloEnv env, const std::shared_ptr<Callbacks> callbacks) {
        return (IloCplex::Callback(new (env) MyHeuristicCallbackI(env, callbacks)));
    }

}