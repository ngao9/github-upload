#include "instance.hpp"

#include "dataprocessing_cpp_include.hpp"
#include "networksolver_cpp_include.hpp"

#include "kevintools_primitives.hpp"
#include "kevintools_functions.hpp"
#include "kevintools_strings.hpp"
#include "kevintools_vectors.hpp"
#include "kevintools_maps.hpp"

#include <algorithm>

using namespace networkdesign;

Instance::Instance() :
    arcScoresFunction(nullptr),
    travelTimes(), travelDistances(), nbPossibles(),
    n(-1), nbModes(-1), timeHorizonInHours(-1), locations(), odPairs(),
    hubIndices(), maximumNumberOfTransfers(-1), keepAllShuttleArcs(false)
{
}

void Instance::setLocations(const Stops& stops) {

    locations = VVdouble(n, Vdouble(2));
    for(int i = 0; i < n; ++i){
        Adouble2 latlong = std::get<1>(stops.allStops.at(i));
        locations.at(i).at(1) = latlong.at(0);
        locations.at(i).at(0) = latlong.at(1);
    }

}

void Instance::setODpairs(const Stops& stops, const Trips& trips) {

    const VTrip& allTrips = trips.allTrips;
    int nbODpairs = allTrips.size();

    odPairs = VAint3();
    odPairs.reserve(nbODpairs);

    for(int k = 0; k < nbODpairs; ++k){

        int from = stops.idToIndex(std::get<0>(allTrips.at(k)));
        int to = stops.idToIndex(std::get<1>(allTrips.at(k)));

        if(from == to){
            std::cerr << "Warning: ignored trip from stop " << std::get<0>(allTrips.at(k)) << " to " << std::get<1>(allTrips.at(k)) << "." << std::endl;
            continue;
            //throw std::runtime_error("Exception: from == to.");
        }

        int nbOfPeople = std::get<2>(allTrips.at(k));

        odPairs.push_back(Aint3{from, to, nbOfPeople});

    }

}

void Instance::setHubIndices(const VAdouble2& desiredHubLatLons, const Vint& additionalHubIndices, const Stops& stops) {

    std::set<int> uniqueHubIndices = std::set<int>(additionalHubIndices.begin(), additionalHubIndices.end());

    for(const Adouble2& desiredHubLatLon : desiredHubLatLons){

        const double& lat = desiredHubLatLon.at(0);
        const double& lon = desiredHubLatLon.at(1);

        uniqueHubIndices.insert(getClosestLocation(lat, lon));

    }

    hubIndices = Vint(uniqueHubIndices.begin(), uniqueHubIndices.end());

}

void Instance::constructFixedRailInstance(const Stops& stops, const RailLines& railLines,
                                          const Trips& trips, const int timeHorizonInHours,
                                          const Vint& possibleNbBuses,
                                          const Distances& roadTravelTimes, const Distances& roadTravelDistances,
                                          double travelTimeFactorShuttle, double travelTimeFactorBus,
                                          const VAdouble2& desiredHubLocations,
                                          const int maximumNumberOfTransfers,
                                          const bool ignoreFrequencyCorrection,
                                          const bool keepAllShuttleArcs){

    n = stops.allStops.size();
    nbModes = 3; //TODO: remove hardcoding

    this->travelTimes = VVVdouble(n, VVdouble(n, Vdouble(nbModes, kevintools::NAN_DOUBLE)));
    this->travelDistances = VVVdouble(n, VVdouble(n, Vdouble(nbModes, kevintools::NAN_DOUBLE)));

    this->timeHorizonInHours = timeHorizonInHours;
    this->maximumNumberOfTransfers = maximumNumberOfTransfers;

    this->keepAllShuttleArcs = keepAllShuttleArcs;

    setLocations(stops);
    setODpairs(stops, trips);

    const Vint& railStationHubIndices = getRailStationIndices(stops, railLines);
    setHubIndices(desiredHubLocations, railStationHubIndices, stops);
    Vint nonRailStationHubIndices = kevintools::removeElementsFromVector(railStationHubIndices, hubIndices);

    //Register travelTimes (minutes), travelDistances (km), and nbPossibles (over time horizon).
    //Duplicate keys are automatically ignored by the addArc methods.
    //The SuperGraph constructed by buildSuperGraph() is based on all defined arcs, except self loops.
    //Choosing an appropriate subgraph for each subproblem is the responsibility of the NetworkDesignProblem.

    //Travel between rail stations with rail transit
    addRailArcsBasedOnSchedule(railLines, railStationHubIndices, stops, ignoreFrequencyCorrection);

    //Travel between non-rail stations with buses
    for(int i : nonRailStationHubIndices){
        for(int j : nonRailStationHubIndices){
            addArcTimeLookupDistanceLookup(i, j, BUS, roadTravelTimes, roadTravelDistances, travelTimeFactorBus);
        }
    }

    //Travel from non-rail stations to the nearest rail stations
    for(int i : nonRailStationHubIndices){
        const Vint& nearestRailStationIndices = getXnearestFromSubset(3, i, railStationHubIndices, roadTravelTimes);
        for(int j : nearestRailStationIndices){
            addArcTimeLookupDistanceLookup(i, j, BUS, roadTravelTimes, roadTravelDistances, travelTimeFactorBus);
            addArcTimeLookupDistanceLookup(j, i, BUS, roadTravelTimes, roadTravelDistances, travelTimeFactorBus);
            for(int k : nearestRailStationIndices){
                addArcTimeLookupDistanceLookup(j, k, BUS, roadTravelTimes, roadTravelDistances, travelTimeFactorBus);
            }
        }
    }

    //Arcs potentially needed for OD pairs
    for(const Aint3& odPair : odPairs){

        const int origin = odPair.at(0);
        const int destination = odPair.at(1);

        addArcTimeLookupDistanceLookup(origin, destination, SHUTTLE, roadTravelTimes, roadTravelDistances, travelTimeFactorShuttle);

        for(int i : hubIndices){
            addArcTimeLookupDistanceLookup(origin, i, SHUTTLE, roadTravelTimes, roadTravelDistances, travelTimeFactorShuttle);
            addArcTimeLookupDistanceLookup(i, destination, SHUTTLE, roadTravelTimes, roadTravelDistances, travelTimeFactorShuttle);
        }

    }

    //Keep all shuttle arcs (may be required by networkdesignproblem)
    if(keepAllShuttleArcs) {
        for (int origin = 0; origin < n; ++origin) {
            for (int destination = 0; destination < n; ++destination) {

                const int& originID = std::get<0>(stops.allStops.at(origin));
                const int& destinationID = std::get<0>(stops.allStops.at(destination));

                if (originID == destinationID) {
                    continue;
                }

                addArcTimeLookupDistanceLookup(origin, destination, SHUTTLE, roadTravelTimes, roadTravelDistances, travelTimeFactorShuttle);

            }
        }
    }

    //Set default number of possibles
    kevintools::safeInsert<tripleIndexVintMap, tripleIndex, Vint>(nbPossibles, {-1, -1, SHUTTLE}, {0});
    kevintools::safeInsert<tripleIndexVintMap, tripleIndex, Vint>(nbPossibles, {-1, -1, BUS}, possibleNbBuses);
    //Rail is explicitly defined for each arc in addRailArcsBasedOnSchedule()

}

std::tuple<double, int, Vint> Instance::getTimeNumberPerHourAndLineIndicesFromSchedule(const RailLines &railLines, const int fromId, const int toId, const bool ignoreFrequencyCorrection) const {

    VRailLine allLines = railLines.allLines;
    const int& nbRailLines = allLines.size();

    double minimumTime = kevintools::NAN_DOUBLE;
    double consistentTime = kevintools::NAN_DOUBLE;
    int nbPerHourOfMinimum = 0;
    int totalNbPerHour = 0;
    Vint lineIndices = Vint();

    for(int i = 0; i < nbRailLines; ++i){

        const int& nb = std::get<1>(allLines.at(i));
        const Vint& stopIDs = std::get<2>(allLines.at(i));
        const Vdouble& forwardSchedule = std::get<3>(allLines.at(i));
        const Vdouble& backwardSchedule = std::get<4>(allLines.at(i));

        auto fromIter = std::find(stopIDs.begin(), stopIDs.end(), fromId);
        auto toIter = std::find(stopIDs.begin(), stopIDs.end(), toId);

        if(fromIter == stopIDs.end() || toIter == stopIDs.end()){
            continue; //fromId and toID not both on this line
        }

        if(fromIter == toIter){
            throw std::runtime_error("Exception: fromIter == toIter");
        }

        int fromIndex = std::distance(stopIDs.begin(), fromIter);
        int toIndex = std::distance(stopIDs.begin(), toIter);

        double time;
        if(fromIndex < toIndex){
            time = forwardSchedule.at(toIndex) - forwardSchedule.at(fromIndex);
        }
        else{
            time = backwardSchedule.at(toIndex) - backwardSchedule.at(fromIndex);
        }

        if(time < 0){
            throw std::runtime_error("Exception: schedule reports negative travel time.");
        }

        if(!std::isnan(minimumTime) && time == minimumTime && nb > nbPerHourOfMinimum){
            nbPerHourOfMinimum = nb;
        }

        if(std::isnan(minimumTime) || time < minimumTime){
            minimumTime = time;
            nbPerHourOfMinimum = nb;
            if(ignoreFrequencyCorrection) {
                lineIndices = {i};
            }
        }

        if(std::isnan(consistentTime)){
            consistentTime = time;
        }

        if(!ignoreFrequencyCorrection && std::fabs(consistentTime - time) > 0.000001){
            throw std::runtime_error("Exception: consistentTime - time = " + std::to_string(consistentTime - time) + " != 0.");
            //This error is thrown if two rail lines report different travel times for the same origin and destination.
            //Note that diverging paths are allowed, as long as the total travel time is the same.
            //There may be instances for which the travel time is actually different, but this is not currently supported by this method.
        }

        if(!ignoreFrequencyCorrection){
            lineIndices.push_back(i);
        }

        totalNbPerHour += nb;

    }

    if(!ignoreFrequencyCorrection) {
        return std::make_tuple(consistentTime, totalNbPerHour, lineIndices);
    }
    else{
        return std::make_tuple(minimumTime, nbPerHourOfMinimum, lineIndices);
    }

}

void Instance::addRailArcsBasedOnSchedule(const RailLines& railLines, const Vint& railStationIndices, const Stops& stops, const bool ignoreFrequencyCorrection){

    for(int i : railStationIndices){
        for(int j : railStationIndices){

            if(i == j){ continue; }

            const int& fromId = std::get<0>(stops.allStops.at(i));
            const int& toId = std::get<0>(stops.allStops.at(j));

            std::tuple<double, int, Vint> timeNbPerHour = getTimeNumberPerHourAndLineIndicesFromSchedule(
                    railLines, fromId, toId, ignoreFrequencyCorrection);
            const double& time = std::get<0>(timeNbPerHour);
            const int& nbPerHour = std::get<1>(timeNbPerHour);
            int nb = nbPerHour * timeHorizonInHours;

            if(std::isnan(time)){
                continue; //i and j do not share a line, see getTimeNumberPerHourAndLineIndicesFromSchedule()
            }

            double distance = 0; //ignore fixed costs

            kevintools::safeInsert(travelTimes, i, j, (int) RAIL, time);
            kevintools::safeInsert(travelDistances, i, j, (int) RAIL, distance);

            tripleIndex key = tripleIndex{i, j, RAIL};
            kevintools::safeInsert<tripleIndexVintMap, tripleIndex, Vint>(nbPossibles, key, Vint{nb});

        }
    }

}

void Instance::addArcTimeLookupDistanceEstimate(const int i, const int j, const travel_mode& mode,
                                                const Distances& roadTravelTimes,
                                                const double travelTimeFactor, const double kmph) {

    if(std::isnan(roadTravelTimes.allDistances.at(i).at(j))){
        return;
    }

    double timeLookup = travelTimeFactor * roadTravelTimes.allDistances.at(i).at(j);
    double distanceEstimate = kmph * (roadTravelTimes.allDistances.at(i).at(j) / 60.0);

    kevintools::safeInsert(travelTimes, i, j, (int) mode, timeLookup);
    kevintools::safeInsert(travelDistances, i, j, (int) mode, distanceEstimate);

}

void Instance::addArcTimeLookupDistanceLookup(const int i, const int j, const travel_mode& mode,
                                              const Distances& roadTravelTimes, const Distances& roadTravelDistances,
                                              const double travelTimeFactor) {

    if(std::isnan(roadTravelTimes.allDistances.at(i).at(j)) || std::isnan(roadTravelDistances.allDistances.at(i).at(j))){
        return;
    }

    double timeLookup = travelTimeFactor * roadTravelTimes.allDistances.at(i).at(j);
    double distanceLookup = roadTravelDistances.allDistances.at(i).at(j);

    kevintools::safeInsert(travelTimes, i, j, (int) mode, timeLookup);
    kevintools::safeInsert(travelDistances, i, j, (int) mode, distanceLookup);

}

void Instance::setArcScoresFunction(const ArcFunction& function) {
    arcScoresFunction = function;
}

std::tuple<double, Vdouble, double, Vdouble> Instance::arcScores(int i, int j, travel_mode mode, int nb) const {

    if(arcScoresFunction == nullptr){
        throw std::runtime_error("Exception: arcScoreFunction not set.");
    }

    if(!arcModeExists(i, j, mode)){
        throw std::runtime_error("Exception: arc does not exist.");
    }

    const double& time = travelTime(i, j, mode);
    const double& distance = travelDistance(i, j, mode);

    if(mode != SHUTTLE) {
        const Vint& nbPoss = nbPossible(i, j, mode);
        if (std::find(nbPoss.begin(), nbPoss.end(), nb) == nbPoss.end()) {
            throw std::runtime_error("Exception: nb not valid.");
        }
    }

    return arcScoresFunction(time, distance, mode, nb, timeHorizonInHours);

}

double Instance::arcCost(const int i, const int j, const travel_mode mode, const int nb, const score_level level) const {

    std::tuple<double, Vdouble, double, Vdouble> scores = arcScores(i, j, mode, nb);

    double cost = 0;
    if(level == DESIGN || level == BOTH){
        cost += std::get<1>(scores).at(0);
    }
    if(level == PASSENGER || level == BOTH){
        cost += std::get<3>(scores).at(0);
    }

    return cost;

}

double Instance::arcConvenience(const int i, const int j, const travel_mode mode, const int nb, const score_level level) const {

    std::tuple<double, Vdouble, double, Vdouble> scores = arcScores(i, j, mode, nb);

    double convenience = 0;
    if(level == DESIGN || level == BOTH){
        convenience += std::get<1>(scores).at(1);
    }
    if(level == PASSENGER || level == BOTH){
        convenience += std::get<3>(scores).at(1);
    }

    return convenience;

}

double Instance::arcObjective(const int i, const int j, const travel_mode mode, const int nb, const score_level level) const {

    std::tuple<double, Vdouble, double, Vdouble> scores = arcScores(i, j, mode, nb);

    double objective = 0;
    if(level == DESIGN || level == BOTH){
        objective += std::get<0>(scores);
    }
    if(level == PASSENGER || level == BOTH){
        objective += std::get<2>(scores);
    }

    return objective;

}

double Instance::travelTime(const int i, const int j, const travel_mode mode) const {

    const double& time = travelTimes.at(i).at(j).at((int) mode);
    if(std::isnan(time)) {
        throw std::runtime_error("Exception: travel time not found.");
    }
    else{
        return time;
    }

}

double Instance::travelDistance(const int i, const int j, const travel_mode mode) const{

    const double& distance = travelDistances.at(i).at(j).at((int) mode);
    if(std::isnan(distance)) {
        throw std::runtime_error("Exception: travel distance not found.");
    }
    else{
        return distance;
    }

}

const Vint& Instance::nbPossible(const int i, const int j, const travel_mode mode) const{

    //If not found, use result for [-1, -1, mode]
    //In this way, nbPossibles only has to keep track of the exceptions

    auto iter = nbPossibles.find(tripleIndex{i, j, (int) mode});

    if(iter == nbPossibles.end()){

        iter = nbPossibles.find(tripleIndex{-1, -1, (int) mode});

        if(iter == nbPossibles.end()){
            throw std::runtime_error("Exception: nbPossible not found.");
        }
        else{
            return iter->second;
        }

    }
    else{
        return iter->second;
    }

}

bool Instance::arcModeExists(const int i, const int j, const travel_mode mode) const {

    bool travelTimeExists = !std::isnan(travelTimes.at(i).at(j).at((int) mode));
    bool travelDistanceExists = !std::isnan(travelDistances.at(i).at(j).at((int) mode));

    return (travelTimeExists && travelDistanceExists);

}

tripleIndexSet Instance::allExistingArcModes() const {

    tripleIndexSet result;

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int mode = 0; mode < nbModes; ++mode){

                if(arcModeExists(i, j, (travel_mode) mode)){
                    result.insert(tripleIndex{i, j, mode});
                }

            }
        }
    }

    return result;

}

pSuperGraph Instance::buildSuperGraph() const {

    pSuperGraph superGraph = std::make_shared<SuperGraph>(n);
    ListDigraph& graph = superGraph->graph;
    auto& costs = superGraph->cost;
    auto& intProperties = superGraph->properties;
    auto& contributions = superGraph->contributions;

    const tripleIndexSet& arcModes = allExistingArcModes();

    for(const Aint3& arcMode : arcModes){

        const int& i = arcMode.at(0);
        const int& j = arcMode.at(1);
        const travel_mode& mode = (travel_mode) arcMode.at(2);

        if(i == j){
            continue;
        }

        const ListDigraph::Node& from = graph.nodeFromId(i);
        const ListDigraph::Node& to = graph.nodeFromId(j);

        const Vint& possibleNb = nbPossible(i, j, mode); //for SHUTTLES, possibleNb={0}

        for(const int& nb : possibleNb) {
            const ListDigraph::Arc& arc = graph.addArc(from, to);
            costs[arc] = arcObjective(i, j, mode, nb, PASSENGER);
            intProperties[arc] = Vint{(int) mode, nb};      //Register mode and nb as integer properties of arc
            contributions[arc] = VResourceContribution();   //Initialize only
        }

    }

    return superGraph;

}

Vint Instance::uniqueNbPossibleForMode(const travel_mode& mode) const {

    intSet uniqueNbPossible = intSet();

    for(auto iter = nbPossibles.begin(); iter != nbPossibles.end(); ++iter){

        if(iter->first.at(2) != mode){
            continue;
        }

        const Vint& nbPossible = iter->second;
        uniqueNbPossible.insert(nbPossible.begin(), nbPossible.end());

    }

    Vint result = Vint(uniqueNbPossible.begin(), uniqueNbPossible.end());
    std::sort(result.begin(), result.end());

    return result;

}

Vint Instance::getUniqueRailStationIDs(const Stops& stops, const RailLines& railLines) const{

    intSet uniqueIDs;
    for(const RailLine& line : railLines.allLines){
        for(const int id : std::get<2>(line)){
            uniqueIDs.insert(id);
        }
    }

    return Vint(uniqueIDs.begin(), uniqueIDs.end());

}

Vint Instance::getRailStationIndices(const Stops &stops, const RailLines &railLines) const {

    const Vint& railStationIDs = getUniqueRailStationIDs(stops, railLines);

    Vint railStationIndices;
    railStationIndices.reserve(railStationIDs.size());

    for(const int id : railStationIDs){
        railStationIndices.push_back(stops.idToIndex(id));
    }

    return railStationIndices;

}

int Instance::getClosestLocation(const double lat, const double lon) const{

    int minIndex = -1;
    double minVal = std::numeric_limits<double>::max();

    Vdouble here = Vdouble{lon, lat};

    for(int k = 0; k < locations.size(); ++k){

        const Vdouble& there = locations.at(k);
        double distance = kevintools::distanceEarth(here.at(1), here.at(0), there.at(1), there.at(0));

        if(distance < minVal){
            minIndex = k;
            minVal = distance;
        }

    }

    return minIndex;

}

Vint Instance::getXnearestFromSubset(const int x, const int index, const Vint& subset, const Distances& distances) const{

    Vint result = Vint();
    if(subset.size() == 0){
        return Vint();
    }

    std::map<double, int> distanceIndexMap;

    for(const int& subsetIndex : subset){
        double distance = distances.allDistances.at(index).at(subsetIndex);
        if(std::isnan(distance)){ continue; }
        distanceIndexMap.insert(std::make_pair(distance, subsetIndex));
    }

    for(auto distanceIndex : distanceIndexMap){

        result.push_back(distanceIndex.second);

        if(result.size() == x || result.size() == subset.size()){
            break;
        }

    }

    return result;

}