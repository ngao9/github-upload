#pragma once

#include "networkdesign_typedefs.hpp"

#include "callbacks.hpp"

namespace networkdesign {

    class NetworkDesignProblem {

    private:
        friend class Callbacks;

    private:
        const Instance& in;


    //Build NetworkDesignProblem
    private:
        void setCPLEXsettings(double mipGap, double timeLimit);
        void createResources();
        void createResourcesArcMode(int i, int j, const travel_mode& mode);
    public:
        NetworkDesignProblem(const Instance& in,
                             double mipGap, double timeLimit,
                             const SubproblemType& subProblemType);


    //Master problem
    private:
        void buildMasterProblem();
        IloEnv env;
        const pIloModel model;          //model of the master problem
        const pIloCplex cplex;          //cplex object of the master problem
        IloNumVar theta;
        VIloIntVar resourceVariables;   //same index as resourceMap
        std::shared_ptr<Callbacks> callback;


    //Subproblems
    public: //public for easy access only
        pSuperGraph superGraph; //Supergraph to which all network subProblems are subgraphs
        VpNetwork subProblems;  //Network subProblems
        Vint nbPeoples;         //Number of people associated with each network
    private:
        void buildSubProblems();
        pNetwork buildSubProblem(int origin, int destination);
#ifdef GUROBI
        pGRBEnv subProblemEnv;  //Shared between subproblems
#endif
        void addArcModeToFilter(int i, int j, const travel_mode& mode, ListDigraph::ArcMap<bool>& arcFilter);
        void removeArcModeFromFilter(int i, int j, const travel_mode& mode, ListDigraph::ArcMap<bool>& arcFilter);
        const SubproblemType subproblemType;


    //Solving and solution information
    public:
        int solve(double eps=1E-6);
        std::pair<Vdouble, VAint5> getSolutionAndResourceMap() const; //solution with resourceMap (see below)
        void printSolution(const Stops& stops) const;
        std::tuple<double, Vdouble, double, Vdouble> getNetworkScores(const Stops& stops, const RailLines& railLines) const;
    private:
        Vdouble getSolution() const;


    //Generate all possible paths through the generated design
    private:
        void removeSuboptimalBusesFromSuperGraph();
    public:
        ListDigraph_VVArc generateAllPaths();


    //Resources
    private:
        int nbResources;
        VAint5 resourceMap;             //resourceMap[i] = (i, j, mode, nbVehicles, dummy (0/1))
        VResourceBounds resourceBounds;
        void verifyResourceMapUnique(); //verify that elements resourceMap are unique, except for dummy's
        void setSuperGraphContributions();
        void initializeResourceBounds();
    public:
        bool isDummyResource(int index) const;
        bool isDummyResource(const Aint5& resource) const;

    };

}