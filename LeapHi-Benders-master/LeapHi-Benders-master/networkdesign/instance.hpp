#pragma once

#include "networkdesign_typedefs.hpp"
#include "dataprocessing_typedefs.hpp"

namespace networkdesign {

    class Instance {

    public:
        Instance();


    //Load instance for which the rail is fixed
    private:
        void setLocations(const Stops& stops);
        void setODpairs(const Stops& stops, const Trips& trips);
        void setHubIndices(const VAdouble2& desiredHubLatLons, const Vint& additionalHubIndices, const Stops& stops);
    public:
        void constructFixedRailInstance(const Stops& stops, const RailLines& railLines,
                                        const Trips& trips, int timeHorizonInHours,
                                        const Vint& possibleNbBuses,
                                        const Distances& roadTravelTimes, const Distances& roadTravelDistances,
                                        double travelTimeFactorShuttle, double travelTimeFactorBus,
                                        const VAdouble2& desiredHubLocations,
                                        int maximumNumberOfTransfers,
                                        bool ignoreFrequencyCorrection,
                                        bool keepAllShuttleArcs);
        //Register rail as fixed resources (no dummies needed in the in master)
        //Note: networkDesignProblem should allow rail arcs to be travelled (see NetworkDesignProblem::initializeResourceBounds())
        //
        //SHUTTLE:  travel times:   travelTimeFactorShuttle * roadTravelTimes
        //          distances:      roadTravelTimes
        //BUS:      travel times:   travelTimeFactorBus * roadTravelTimes
        //          distances:      roadTravelTimes
        //RAIL:     travel times:   immediately from rail schedule
        //          distances:      0
        //
        //desiredNonRailHubLocations[i] = {lat, lon} of desired location hub i.
        //The actual location is the stop (rail or other) closest to the location (Earth distance).
        //Rail stations are always hubs. Desired locations overlapping with rail stations are automatically ignored.


    //Adding arcs
    public:
        std::tuple<double, int, Vint> getTimeNumberPerHourAndLineIndicesFromSchedule(const RailLines& railLines, const int fromId, const int toId, const bool ignoreFrequencyCorrection) const;
            //returns {time, nbPerHour}. If fromId and toId not on any same line, returns {NaN, 0, Vint()}.
    private:
        void addRailArcsBasedOnSchedule(const RailLines& railLines, const Vint& railStationIndices, const Stops& stops, bool ignoreFrequencyCorrection);
        void addArcTimeLookupDistanceEstimate(int i, int j, const travel_mode& mode,
                                              const Distances& roadTravelTimes,
                                              double travelTimeFactor, double kmph);
        void addArcTimeLookupDistanceLookup(int i, int j, const travel_mode& mode,
                                            const Distances& roadTravelTimes, const Distances& roadTravelDistances,
                                              double travelTimeFactor);


    //Arc score functions
    private:
        ArcFunction arcScoresFunction; //scores incurred in design level and passenger level (!= what is paid by customer)
        //std::tuple<double designObjective, Vdouble{designCost, designConvenience},
        //double passengerObjective, Vdouble{passengerCost, passengerConvenience}> =
                // arcScoresFunction(double time, double distance, travel_mode mode, int nb, int timeHorizonInHours)
    public:
        enum score_level {DESIGN, PASSENGER, BOTH};
        void setArcScoresFunction(const ArcFunction& function);
        std::tuple<double, Vdouble, double, Vdouble> arcScores(int i, int j, travel_mode mode, int nb) const;
        double arcCost(int i, int j, travel_mode mode, int nb, score_level level) const;
        double arcConvenience(int i, int j, travel_mode mode, int nb, score_level level) const;
        double arcObjective(int i, int j, travel_mode mode, int nb, score_level level) const;


    //Information about times, distances, and frequencies
    private:
        VVVdouble travelTimes;
        VVVdouble travelDistances;
        tripleIndexVintMap nbPossibles; //see nbPossible() for further information
    public:
        double travelTime(int i, int j, travel_mode mode) const;         //travel time from i to j with mode
        double travelDistance(int i, int j, travel_mode mode) const;     //travel distance from i to j with mode
        const Vint& nbPossible(int i, int j, travel_mode mode) const;    //possible number of vehicles over time horizon
        bool arcModeExists(int i, int j, travel_mode mode) const;        //direct travel is not allowed if key does not exists
                                                                                    //true <=> both travel distance and travel time are defined
        tripleIndexSet allExistingArcModes() const;                      //all {i, j, mode} for which arcModeExists = true


    //Build SuperGraph based on instance
    public:
        pSuperGraph buildSuperGraph() const; //based on instance, does not include resource contributions


    //Helper functions
    public:
        Vint uniqueNbPossibleForMode(const travel_mode& mode) const;
        Vint getUniqueRailStationIDs(const Stops &stops, const RailLines &RailLines) const;
        Vint getRailStationIndices(const Stops &stops, const RailLines &RailLines) const;
    private:
        int getClosestLocation(double lat, double lon) const;
        Vint getXnearestFromSubset(int x, int index, const Vint& subset, const Distances& distances) const;


    //Public for easy access, private otherwise
    public:
        int n;                          //number of vertices
        int nbModes;                    //number of modes
        int timeHorizonInHours;         //time horizon in hours
        VVdouble locations;             //locations[i] = {lon, lat} (not a typo), index corresponds to vertices
        VAint3 odPairs;                 //odPairs[i] = (origin, destination, number of people)
        Vint hubIndices;                //indices that correspond to the potential hubs, includes rail stations
        int maximumNumberOfTransfers;   //= maximum arcs per passenger path - 1
                                        //value -1 if no transfer limit imposed
        bool keepAllShuttleArcs;        //keep shuttle arcs for origins and destinations not in odPairs
                                        //may be required by networkdesignproblem

    };

}