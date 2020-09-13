#pragma once

#include <rapidjson/document.h>
#include "networkdesign_typedefs.hpp"
#include "dataprocessing_typedefs.hpp"

namespace networkdesign {

    class Runner {

    //lambdas
    private:
        std::function<int(const double)> min2sec = [](const double min)->int{ return (int)std::ceil(60 * min); };
        std::function<int(const double)> hour2sec = [](const double hour)->int{ return (int)std::ceil(3600 * hour); };
        std::function<int(const double)> km2m = [](const double km)->int{ return (int)std::ceil(1000 * km); };
        std::function<int(const travel_mode, const int)> waitingTimeSec = [this](const travel_mode mode, const int nb)->int{
            if(mode == SHUTTLE){ return 0; }
            else{ return min2sec(fixedTransferTime + (60.0 * timeHorizonInHours)/(2.0 * nb)); }
        };


    //Construct and run instance
    private:
        bool argsProvided(int argc);
        void constructFromHardCoded();
        void constructFromArguments(int argc, char* argv[]);
        void buildInstance();
    public:
        Runner(int argc, char* argv[], const std::clock_t& start);


    //JSON output
    private:
        void prepareParametersJSON(rapidjson::Document& output);
        void prepareDesignJSON(rapidjson::Document& output, ListDigraph::ArcMap<int>& legIndex);
        void prepareRoutesJSON(rapidjson::Document& output);
        void prepareTripsJSON(rapidjson::Document& output, const ListDigraph::ArcMap<int>& legIndex);
        void prepareTripSplittingJSON(rapidjson::Document& output, const ListDigraph::ArcMap<int>& legIndex);
        void prepareScoresJSON(rapidjson::Document& output);
    public:
        void writeJSONtoFile();


    //Helper methods
    private:
        void printTripStatistics();
        void addGeneralInformationToLeg(rapidjson::Value& leg, rapidjson::Document& output, int i, int j, travel_mode mode, int nb);
        void addDesignInformationToLeg(rapidjson::Value& leg, rapidjson::Document& output, int i, int j, travel_mode mode, int nb);
        void addTripInformationToLeg(rapidjson::Value& leg, rapidjson::Document& output, int i, int j, travel_mode mode, int nb, const ListDigraph::Arc& superArc, const ListDigraph::ArcMap<int>& legIndex);
        void addTripSplittingInformationToLeg(rapidjson::Value& leg, rapidjson::Document& output, int i, int j, travel_mode mode, int nb, const ListDigraph::Arc& superArc, const ListDigraph::ArcMap<int>& legIndex);


    //Parameters
    public:

        //Score parameters
        double alpha;
        double shuttleCostPerKm;
        double busCostPerHour;
        double fixedTransferTime;
        double passengerFactor;

        //Other parameters and prepared information
        pStops stops;
        pRailLines railLines;
        bool ignoreFrequencyCorrection;
        pTrips trips;
        int timeHorizonInHours;
        Vint possibleNbBuses;
        pDistances roadTravelTimes;
        pDistances roadTravelDistances;
        double travelTimeFactorShuttle;
        double travelTimeFactorBus;
        VAdouble2 desiredHubLatLons;
        int maximumNumberOfTransfers;

        //Output parameters
        std::string outputFile;
        bool generateTripSplittings;
        bool enableGNUplot;

        //Optimization parameters
        double mipGap;
        double timeLimit;
        SubproblemType subproblemType;


    //Instance and solution
    public:
        std::function<std::tuple<double, Vdouble, double, Vdouble>
                (const double, const double, const travel_mode, const int, const int)> arcScores();
        pInstance instance;
        Vint uniqueNbPossibleForBus;
        pNetworkDesignProblem networkDesignProblem;
        Vdouble solution;
        VAint5 resourceMap;
        int status;
        std::tuple<double, Vdouble, double, Vdouble> networkScores;
        ListDigraph_VVArc passengerPaths;
        ListDigraph_VVArc tripSplittings;

    };

}