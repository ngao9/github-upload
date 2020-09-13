#include "runner.hpp"

#include "networkdesign_cpp_include.hpp"
#include "kevintools_vectors.hpp"

#include <ctime>
#include <sstream>
#include <numeric>
#include <forward_list>

#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/pointer.h>
namespace rj = rapidjson;

#ifdef HAS_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

using namespace networkdesign;

#define rjNull rj::Value()
#define rjArray rj::Value(rj::kArrayType)
#define rjObject rj::Value(rj::kObjectType)
#define rjDouble(X) rj::Value().SetDouble(X)
#define rjString(X,Y) rj::Value().SetString(X.c_str(), Y)
#define rjNumString(X,Y) rj::Value().SetString(std::to_string(X).c_str(), Y)

bool Runner::argsProvided(const int argc){

    if(argc == 1){
        std::cerr << "Warning: no command line arguments provided, using hardcoded parameters. Include --help flag for more information." << std::endl;
        return false;
    }
    else{
#ifdef HAS_BOOST_PROGRAM_OPTIONS
        return true;
#else
        std::cerr << "Exception: command line arguments provided, but executable was not compiled with Boost. Either use hardcoded parameters, or recompile with Boost." << std::endl;
        throw std::runtime_error("Exception: command line arguments provided, but executable was not compiled with Boost. Either use hardcoded parameters, or recompile with Boost.");
#endif
    }

}

void Runner::constructFromHardCoded() {

    alpha = 0.05;
    shuttleCostPerKm = 1/mile2km;
    busCostPerHour = 34;
    double fixedTransferTimeInSeconds = 0;
    fixedTransferTime = fixedTransferTimeInSeconds/60.0;
    passengerFactor = 2;

    stops = std::make_shared<Stops>("../data/clustered_stops_1500feet_jan1_oct31.csv");
    railLines = std::make_shared<RailLines>("../data/rail_lines.csv");
    ignoreFrequencyCorrection = false;
    trips = std::make_shared<Trips>("../data/odx_clustered_april3_jan1_oct31.csv");
    timeHorizonInHours = 4;
    Vint busesPerHour = Vint{2, 3, 4};
    possibleNbBuses = kevintools::operator*(busesPerHour, timeHorizonInHours);
    roadTravelTimes = std::make_shared<Distances>("../data/duration_matrix_jan1_oct31.csv", *stops);
    roadTravelDistances = std::make_shared<Distances>("../data/distance_matrix_jan1_oct31.csv", *stops);
    travelTimeFactorShuttle = 1;
    travelTimeFactorBus = 1;
    desiredHubLatLons = VAdouble2{   Adouble2{33.640553, -84.446198},
                                     Adouble2{33.789669, -84.387414},
                                     Adouble2{33.756478, -84.41723},
                                     Adouble2{33.775554, -84.281487},
                                     Adouble2{33.772159, -84.428873},
                                     Adouble2{33.860329, -84.339245},
                                     Adouble2{33.847944, -84.367716},
                                     Adouble2{33.887929, -84.305556},
                                     Adouble2{33.766111, -84.387504},
                                     Adouble2{33.650282, -84.448612},
                                     Adouble2{33.774699, -84.295417},
                                     Adouble2{33.756179, -84.397215},
                                     Adouble2{33.902787, -84.28061},
                                     Adouble2{33.92113, -84.344268},
                                     Adouble2{33.765241, -84.312937},
                                     Adouble2{33.677424, -84.440542},
                                     Adouble2{33.761803, -84.340825},
                                     Adouble2{33.753899, -84.39156},
                                     Adouble2{33.748696, -84.395741},
                                     Adouble2{33.750161, -84.385915},
                                     Adouble2{33.754517, -84.469824},
                                     Adouble2{33.769856, -84.228906},
                                     Adouble2{33.757451, -84.352762},
                                     Adouble2{33.772647, -84.251607},
                                     Adouble2{33.749951, -84.375675},
                                     Adouble2{33.70088, -84.428768},
                                     Adouble2{33.845388, -84.358235},
                                     Adouble2{33.823501, -84.369313},
                                     Adouble2{33.910757, -84.35189},
                                     Adouble2{33.781247, -84.386342},
                                     Adouble2{33.771472, -84.387258},
                                     Adouble2{33.94495, -84.357275},
                                     Adouble2{33.7173, -84.42503},
                                     Adouble2{33.758082, -84.387595},
                                     Adouble2{33.931671, -84.351069},
                                     Adouble2{33.756613, -84.403902},
                                     Adouble2{33.736564, -84.413653},
                                     Adouble2{33.753247, -84.445568},
                                     Adouble2{34.040119, -84.321837},
                                     Adouble2{33.705968, -84.137274},
                                     Adouble2{33.589076, -84.548713},
                                     Adouble2{33.572306, -84.358125},
                                     Adouble2{33.644274, -84.325321},
                                     Adouble2{33.669692, -84.247208},
                                     Adouble2{33.723549, -84.585516},
                                     Adouble2{33.536876, -84.418154},
                                     Adouble2{34.0613, -84.25215},
                                     Adouble2{33.506489, -84.357707},
                                     Adouble2{33.549047, -84.604923},
                                     Adouble2{33.448848, -84.342723},
                                     Adouble2{33.521085, -84.667191},
                                     Adouble2{33.833683, -84.272498},
                                     Adouble2{33.70673, -84.308166},
                                     Adouble2{33.821117, -84.194302},
                                     Adouble2{33.70415, -84.506394},
                                     Adouble2{33.584586, -84.469255},
                                     Adouble2{33.911219, -84.424701} };
    maximumNumberOfTransfers = 3;

    outputFile = "../data/output.json";
    generateTripSplittings = true;
    enableGNUplot = true;

    mipGap = 0;
    timeLimit = 60;
    subproblemType = SubproblemType::LEMONNETWORK;

#ifndef GNUPLOT
    enableGNUplot = false;
#endif

}

void Runner::constructFromArguments(int argc, char** argv) {
#ifdef HAS_BOOST_PROGRAM_OPTIONS

    po::options_description description("LeapHi-Benders Parameters");
    description.add_options()
            ("help", "See this description")
            ("stopFile", po::value<std::string>(), "Filepath for with stops file")
            ("railLineFile", po::value<std::string>(), "Filepath rail lines files")
            ("ignoreFrequencyCorrection", "Ignore rail frequency correction")
            ("tripsFile", po::value<std::string>(), "Filepath for trips file")
            ("timeHorizonInHours", po::value<int>(), "Time horizon in hours")
            ("busesPerHour", po::value<Vint>()->multitoken(), "Number of buses per hour")
            ("roadTravelTimesFile", po::value<std::string>(), "Filepath for road travel times file")
            ("roadTravelDistancesFile", po::value<std::string>(), "Filepath for road travel distances file")
            ("travelTimeFactorShuttle", po::value<double>(), "Multiplicationfactor to obtain shuttle travel times")
            ("travelTimeFactorBus", po::value<double>(), "Multiplicationfactor to obtain bus travel times")
            ("desiredHubLats", po::value<Vdouble>()->multitoken()->zero_tokens(), "Latitudes of desired hub locations")
            ("desiredHubLons", po::value<Vdouble>()->multitoken()->zero_tokens(), "Longitudes of desired hub locations")
            ("maximumNumberOfTransfers", po::value<int>(), "Maximum number of transfers")
            ("outputFile", po::value<std::string>(), "Output file") //all output: if already exists, will be overwritten without warning
            ("generateTripSplittings", "Generate trip splittings (creates overhead in optimization)")
            ("enableGNUplot", "Enable GNUplot plotting")
            ("alpha", po::value<double>(), "Weight for cost and convenience (higher is more convenient)")
            ("shuttleCostPerKm", po::value<double>(), "Shuttle cost per kilometer")
            ("busCostPerHour", po::value<double>(), "Bus cost per hour")
            ("fixedTransferTime", po::value<double>(), "Fixed transfer time in seconds at the start of each non-shuttle arc")
            ("passengerFactor", po::value<double>(), "Multiply the number of passengers with this factor")
            ("mipGap", po::value<double>(), "Mip gap master problem (fraction, 0.01 = 1%)")
            ("timeLimit", po::value<double>(), "Time limit in seconds")
            ("useGurobi", po::value<int>(), "Use Gurobi instead of LEMON to solve the subproblems");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).style(
            po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
    po::notify(vm);

    if(vm.count("help")){
        std::cerr << description << std::endl;
        exit(0);
    }

    alpha = vm["alpha"].as<double>();
    shuttleCostPerKm = vm["shuttleCostPerKm"].as<double>();
    busCostPerHour = vm["busCostPerHour"].as<double>();
    fixedTransferTime = vm["fixedTransferTime"].as<double>()/60.0;
    passengerFactor = vm["passengerFactor"].as<double>();

    stops = std::make_shared<Stops>(vm["stopFile"].as<std::string>());
    railLines = std::make_shared<RailLines>(vm["railLineFile"].as<std::string>());
    ignoreFrequencyCorrection = vm.count("ignoreFrequencyCorrection");
    trips = std::make_shared<Trips>(vm["tripsFile"].as<std::string>());
    timeHorizonInHours = vm["timeHorizonInHours"].as<int>();
    Vint busesPerHour = vm["busesPerHour"].as<Vint>();
    possibleNbBuses = kevintools::operator*(busesPerHour, timeHorizonInHours);
    roadTravelTimes = std::make_shared<Distances>(vm["roadTravelTimesFile"].as<std::string>(), *stops);
    roadTravelDistances = std::make_shared<Distances>(vm["roadTravelDistancesFile"].as<std::string>(), *stops);
    travelTimeFactorShuttle = vm["travelTimeFactorShuttle"].as<double>();
    travelTimeFactorBus = vm["travelTimeFactorBus"].as<double>();
    Vdouble desiredHubLats = Vdouble();
    if(vm.count("desiredHubLats")) {
        desiredHubLats = vm["desiredHubLats"].as<Vdouble>();
    }
    Vdouble desiredHubLons = Vdouble();
    if(vm.count("desiredHubLons")) {
        desiredHubLons = vm["desiredHubLons"].as<Vdouble>();
    }
    desiredHubLatLons = VAdouble2();
    if(desiredHubLats.size() != desiredHubLons.size()){
        throw std::runtime_error("Exception: desiredHubLats.size() != desiredHubLons.size()");
    }
    for(int i = 0; i < desiredHubLats.size(); ++i){
        desiredHubLatLons.push_back({desiredHubLats.at(i), desiredHubLons.at(i)});
    }
    maximumNumberOfTransfers = vm["maximumNumberOfTransfers"].as<int>();

    outputFile = vm["outputFile"].as<std::string>();
    generateTripSplittings = vm.count("generateTripSplittings");
    enableGNUplot = false;
#ifdef GNUPLOT
    enableGNUplot = vm.count("enableGNUplot");
#endif

    mipGap = vm["mipGap"].as<double>();
    timeLimit = vm["timeLimit"].as<double>();
    subproblemType = SubproblemType::LEMONNETWORK;
    if(vm.count("useGurobi")){
        subproblemType = SubproblemType::GUROBINETWORK;
    }

#endif
}

void Runner::buildInstance() {

    instance = std::make_shared<Instance>();
    instance->constructFixedRailInstance(*stops, *railLines, *trips, timeHorizonInHours, possibleNbBuses, *roadTravelTimes, *roadTravelDistances,
                                         travelTimeFactorShuttle, travelTimeFactorBus, desiredHubLatLons, maximumNumberOfTransfers, ignoreFrequencyCorrection, generateTripSplittings);
    instance->setArcScoresFunction(arcScores());
    uniqueNbPossibleForBus = instance->uniqueNbPossibleForMode(BUS);

}

Runner::Runner(const int argc, char* argv[], const std::clock_t& start) {

    auto time = [](const std::clock_t& start){
        return std::to_string((int) (100 * (std::clock() - start) / ((double) CLOCKS_PER_SEC))/100.0) + "s: ";
    };

    auto printHeader = [](const std::string header, std::ostream& stream){
        stream << std::endl;
        for(int i = 0; i < header.length(); ++i){
            stream << "=";
        }
        stream << std::endl << header << std::endl;
        for(int i = 0; i < header.length(); ++i){
            stream << "=";
        }
        stream << std::endl;
    };

    std::ostringstream timeSummary;

    printHeader("TIME SUMMARY", timeSummary);
    timeSummary << time(start) << "start reading parameters" << std::endl;
    argsProvided(argc) ? constructFromArguments(argc, argv) : constructFromHardCoded();
    timeSummary << time(start) << "complete reading parameters" << std::endl;

    printTripStatistics();

    timeSummary << time(start) << "start building instance" << std::endl;
    buildInstance();
    timeSummary << time(start) << "complete building instance" << std::endl;

    timeSummary << time(start) << "start building network design problem" << std::endl;
    networkDesignProblem = std::make_shared<NetworkDesignProblem>(*instance, mipGap, timeLimit, subproblemType);
    timeSummary << time(start) << "complete building network design problem" << std::endl;

    timeSummary << time(start) << "start solving network design problem" << std::endl;
    printHeader("CPLEX LOG", std::cout);
    status = networkDesignProblem->solve();
    std::tie(solution, resourceMap) = networkDesignProblem->getSolutionAndResourceMap();
    timeSummary << time(start) << "solver completed with status = " << std::to_string(status) << std::endl;

    printHeader("SOLUTION", timeSummary);
    networkDesignProblem->printSolution(*stops);

    timeSummary << time(start) << "start calculating network scores" << std::endl;
    networkScores = networkDesignProblem->getNetworkScores(*stops, *railLines);
    timeSummary << time(start) << "complete calculating network scores" << std::endl;

    timeSummary << time(start) << "start retrieving optimal passenger paths" << std::endl;
    const VpNetwork& networks = networkDesignProblem->subProblems;
    passengerPaths.reserve(networks.size());
    for(int p = 0; p < networks.size(); ++p) {
        const pNetwork& network = networks.at(p);
        passengerPaths.push_back(network->getArcFlowPath());
    }
    timeSummary << time(start) << "complete retrieving optimal passenger paths" << std::endl;

    if(generateTripSplittings) {
        timeSummary << time(start) << "start generating all possible trips" << std::endl;
        tripSplittings = networkDesignProblem->generateAllPaths();
        timeSummary << time(start) << "completed generating all possible trips" << std::endl;
    }

    std::cout << timeSummary.str();

}

void Runner::prepareParametersJSON(rj::Document& output){

    auto& allocator = output.GetAllocator();
    output.AddMember("parameters", rjObject, allocator);

    rj::Value& parameters = output["parameters"];
    parameters.AddMember("alpha", alpha, allocator);
    parameters.AddMember("shuttleCostPerKm", shuttleCostPerKm, allocator);
    parameters.AddMember("shuttleCostPerMile", shuttleCostPerKm * mile2km, allocator);
    parameters.AddMember("busCostPerHour", busCostPerHour, allocator);
    parameters.AddMember("fixedTransferTime", 60 * fixedTransferTime, allocator);
    parameters.AddMember("passengerFactor", passengerFactor, allocator);
    parameters.AddMember("frequencies_per_hour", rjArray, allocator);
    parameters.AddMember("maximum_number_of_transfers", maximumNumberOfTransfers, allocator);
    parameters.AddMember("time_horizon", hour2sec(timeHorizonInHours), allocator);
    parameters.AddMember("modes", rjArray, allocator);
    parameters.AddMember("time_units", "seconds", allocator);
    parameters.AddMember("distance_units", "meters", allocator);
    parameters.AddMember("travelTimeFactorShuttle", travelTimeFactorShuttle, allocator);
    parameters.AddMember("travelTimeFactorBus", travelTimeFactorBus, allocator);

    //frequencies_per_hour
    for(const int& nb : uniqueNbPossibleForBus){
        double frequency_per_hour = nb/((double) timeHorizonInHours);
        parameters["frequencies_per_hour"].PushBack(frequency_per_hour, allocator);
    }

    //modes
    for(int i = 0; i < 3; ++i){
        if     (i == (int) SHUTTLE){ parameters["modes"].PushBack("shuttle", allocator); }
        else if(i == (int) RAIL)   { parameters["modes"].PushBack("rail", allocator); }
        else if(i == (int) BUS)    { parameters["modes"].PushBack("bus", allocator); }
        else                       { throw std::runtime_error("Exception: unknown mode."); }
    }

}

void Runner::prepareDesignJSON(rj::Document& output, ListDigraph::ArcMap<int>& legIndex) {

    auto& allocator = output.GetAllocator();

    output.AddMember("design", rjObject, allocator);
    rj::Value& design = output["design"];

    design.AddMember("hubs", rjObject, allocator);
    design.AddMember("legs", rjArray, allocator);

    //hubs
    rj::Value& hubs = design["hubs"];

    const Vint& hubIndices = instance->hubIndices;
    const Vint& railStationIndices = instance->getRailStationIndices(*stops, *railLines);

    for(const int& hubIndex : hubIndices){

        int hubID = std::get<0>(stops->allStops.at(hubIndex));
        hubs.AddMember(rjNumString(hubID, allocator), rjObject, allocator);
        rj::Value& hub = hubs[rjNumString(hubID, allocator)];

        bool isRail = std::find(railStationIndices.begin(), railStationIndices.end(), hubIndex) != railStationIndices.end();

        hub.AddMember("fixed", isRail, allocator);
        hub.AddMember("opened", true, allocator);

    }

    //legs
    rj::Value& legs = design["legs"];

    for(int i : railStationIndices){
        for(int j : railStationIndices){
            if(instance->arcModeExists(i, j, RAIL)) {

                const Vint& nbPossible = instance->nbPossible(i, j, RAIL);
                if (nbPossible.size() != 1) {
                    throw std::runtime_error("Exception: multiple arcs with same origin, destination, and mode.");
                }
                int nb = nbPossible.at(0);

                //Update legIndex
                pSuperGraph superGraph = networkDesignProblem->superGraph;
                const ListDigraph::Node from = superGraph->graph.nodeFromId(i);
                const ListDigraph::Node to = superGraph->graph.nodeFromId(j);
                ListDigraph::Arc superArc = superGraph->getUniqueMatchingArc(from, to, {(int) RAIL, nb});
                legIndex[superArc] = legs.Size();

                //Create leg
                rj::Value leg(rj::kObjectType);
                addGeneralInformationToLeg(leg, output, i, j, RAIL, nb);
                addDesignInformationToLeg(leg, output, i, j, RAIL, nb);
                legs.PushBack(leg, allocator);

            }
        }
    }

    for(int r = 0; r < solution.size(); ++r){

        const Aint5& mapping = resourceMap.at(r);

        if(!networkDesignProblem->isDummyResource(r) && solution.at(r) > 0.5){

            const int& i = mapping.at(0);
            const int& j = mapping.at(1);
            const int& mode = mapping.at(2);
            const int& nb = mapping.at(3);

            if(mode != (int) BUS){
                throw std::runtime_error("Exception: mode == BUS expected.");
            }

            //Update legIndex
            pSuperGraph superGraph = networkDesignProblem->superGraph;
            const ListDigraph::Node from = superGraph->graph.nodeFromId(i);
            const ListDigraph::Node to = superGraph->graph.nodeFromId(j);
            ListDigraph::Arc superArc = superGraph->getUniqueMatchingArc(from, to, {(int) BUS, nb});
            legIndex[superArc] = legs.Size();

            //Create leg
            rj::Value leg(rj::kObjectType);
            addGeneralInformationToLeg(leg, output, i, j, BUS, nb);
            addDesignInformationToLeg(leg, output, i, j, BUS, nb);
            legs.PushBack(leg, allocator);

        }
    }

}

void Runner::prepareRoutesJSON(rapidjson::Document& output){

    if(!output.HasMember("design")){
        throw std::runtime_error("Exception: prepareDesignJSON() must be called before prepareRoutesJSON().");
    }

    auto& allocator = output.GetAllocator();

    output.AddMember("routes", rjObject, allocator);
    rj::Value& routes = output["routes"];

    //Rail routes
    for(const RailLine& line : railLines->allLines){

        const std::string& name = std::get<0>(line);
        const int& nbPerHour = std::get<1>(line);
        const Vint& stopIDs = std::get<2>(line);

        routes.AddMember(rjString(name, allocator), rjObject, allocator);
        rj::Value& route = routes[rjString(name, allocator)];

        route.AddMember("mode", (int) RAIL, allocator);
        route.AddMember("frequency_per_hour", nbPerHour, allocator);
        route.AddMember("stop_ids", rjArray, allocator);

        //stop_ids
        for(int k = 0; k < stopIDs.size(); ++k){
            route["stop_ids"].PushBack(rjNumString(stopIDs.at(k), allocator), allocator);
        }
        for(int k = stopIDs.size() - 2; k >= 1; --k){
            route["stop_ids"].PushBack(rjNumString(stopIDs.at(k), allocator), allocator);
        }

    }

    //Bus routes
    pSuperGraph superGraph = networkDesignProblem->superGraph;
    const ListDigraph& graph = superGraph->graph;

    //Create graph with bus arcs only
    ListDigraph::ArcMap<int> nbMap(graph, 0);
    ListDigraph::ArcMap<bool> arcFilter(graph, false);
    for(int r = 0; r < solution.size(); ++r) {

        const Aint5& mapping = resourceMap.at(r);

        if (!networkDesignProblem->isDummyResource(r) && solution.at(r) > 0.5) {

            const int& i = mapping.at(0);
            const int& j = mapping.at(1);
            const int& mode = mapping.at(2);
            const int& nb = mapping.at(3);

            if (mode != (int) BUS) {
                throw std::runtime_error("Exception: mode == BUS expected.");
            }

            const ListDigraph::Node from = superGraph->graph.nodeFromId(i);
            const ListDigraph::Node to = superGraph->graph.nodeFromId(j);
            ListDigraph::Arc superArc = superGraph->getUniqueMatchingArc(from, to, {(int) BUS, nb});

            nbMap[superArc] += nb;
            arcFilter[superArc] = true;

        }

    }

    //Prepare BFS algorithm
    typedef lemon::FilterArcs<const ListDigraph> FilteredGraph;
    FilteredGraph filteredGraph(graph, arcFilter);
    lemon::Bfs<FilteredGraph> bfs(filteredGraph);

    auto nbArcs = [](const FilteredGraph& filteredGraph)->int{
        int count = 0;
        for(FilteredGraph::ArcIt arcIt(filteredGraph); arcIt != lemon::INVALID; ++arcIt){
            ++count;
        }
        return count;
    };

    //Construct a cycle and remove it from the graph
    int busRoutesFound = 0;
    while(nbArcs(filteredGraph) > 0){

        //Select arbitrary arc
        FilteredGraph::ArcIt arcIt(filteredGraph);
        FilteredGraph::Node from = graph.source(arcIt);
        FilteredGraph::Node to = graph.target(arcIt);

        //Use BFS to find the return path to complete the cycle
        bfs.run(to, from);
        if(!bfs.reached(from)){
            throw std::runtime_error("Exception: no cyclic route could be found.");
        }

        //Determine cycle arcs and nodes
        std::forward_list<FilteredGraph::Arc> cycleArcs = std::forward_list<FilteredGraph::Arc>();
        std::forward_list<FilteredGraph::Node> cycleNodes = std::forward_list<FilteredGraph::Node>();

        FilteredGraph::Node node = from;
        while(node != lemon::INVALID){
            if(bfs.predArc(node) != lemon::INVALID){
                cycleArcs.push_front(bfs.predArc(node));
            }
            cycleNodes.push_front(node);
            node = bfs.predNode(node);
        }
        cycleArcs.push_front(arcIt);

        //Determine maximum possible flow over cycle
        int minNb = std::numeric_limits<int>::max();
        for(const FilteredGraph::Arc& cycleArc : cycleArcs){
            if(nbMap[cycleArc] < minNb){
                minNb = nbMap[cycleArc];
            }
        }

        //Add route information
        std::string name = "BUS_" + std::to_string(busRoutesFound++);
        routes.AddMember(rjString(name, allocator), rjObject, allocator);
        rj::Value& route = routes[rjString(name, allocator)];

        route.AddMember("mode", (int) BUS, allocator);
        route.AddMember("frequency_per_hour", minNb/((double)(timeHorizonInHours)), allocator);
        route.AddMember("stop_ids", rjArray, allocator);

        for(const FilteredGraph::Node& cycleNode : cycleNodes){
            int index = graph.id(cycleNode);
            int id = std::get<0>(stops->allStops.at(index));
            route["stop_ids"].PushBack(rjNumString(id, allocator), allocator);
        }

        //Register route at the legs of the cycle
        for(const FilteredGraph::Arc& cycleArc : cycleArcs){

            int i = graph.id(graph.source(cycleArc));
            int j = graph.id(graph.target(cycleArc));

            int iID = std::get<0>(stops->allStops.at(i));
            int jID = std::get<0>(stops->allStops.at(j));

            for(auto& leg : output["design"]["legs"].GetArray()){
                if(leg["mode"].GetInt() == (int) BUS &&
                   leg["board_stop_id"].GetString() == std::to_string(iID) &&
                   leg["alight_stop_id"].GetString() == std::to_string(jID)){
                    leg["routes"].PushBack(rjString(name, allocator), allocator);
                }
            }

        }

        //Decrease remaining flows and remove at least one arc
        for(const FilteredGraph::Arc& cycleArc : cycleArcs){
            nbMap[cycleArc] -= minNb;
            if(nbMap[cycleArc] == 0){
                arcFilter[cycleArc] = false;
            }
        }

    }

}

void Runner::prepareTripsJSON(rj::Document& output, const ListDigraph::ArcMap<int>& legIndex){

    auto& allocator = output.GetAllocator();

    output.AddMember("trips", rjArray, allocator);

    const VpNetwork& networks = networkDesignProblem->subProblems;
    const pSuperGraph& superGraph = networkDesignProblem->superGraph;
    const ListDigraph& graph = superGraph->graph;

    for(int p = 0; p < networks.size(); ++p) {

        output["trips"].PushBack(rjObject, allocator);
        rj::Value& trip = output["trips"][p];

        trip.AddMember("passengers", networkDesignProblem->nbPeoples.at(p), allocator);
        trip.AddMember("origin_stop_id", rjObject, allocator);
        trip.AddMember("destination_stop_id", rjObject, allocator);
        trip.AddMember("departure_times", rjArray, allocator);
        trip.AddMember("travel_time", 0, allocator);
        trip.AddMember("waiting_time", 0, allocator);
        trip.AddMember("total_time", 0, allocator);
        trip.AddMember("distance", 0, allocator);
        trip.AddMember("objective", 0.0, allocator);
        trip.AddMember("legs", rjArray, allocator);

        const ListDigraph_VArc& path = passengerPaths.at(p);

        //stop_ids
        int originID = std::get<0>(trips->allTrips.at(p));
        int destinationID = std::get<1>(this->trips->allTrips.at(p));

        trip["origin_stop_id"] = rjNumString(originID, allocator);
        trip["destination_stop_id"] = rjNumString(destinationID, allocator);

        //departure times
        for(const std::string& time : std::get<3>(trips->allTrips.at(p))){
            trip["departure_times"].PushBack(rjString(time, allocator), allocator);
        }

        //legs and aggregate information
        rj::Value& legs = trip["legs"];
        for (int k = 0; k < path.size(); ++k) {

            const ListDigraph::Arc& superArc = path.at(k);
            const int& i = graph.id(graph.source(superArc));
            const int& j = graph.id(graph.target(superArc));
            const travel_mode& mode = (travel_mode) superGraph->properties[superArc].at(0);
            const int& nb = superGraph->properties[superArc].at(1);

            rj::Value leg(rj::kObjectType);
            addGeneralInformationToLeg(leg, output, i, j, mode, nb);
            addTripInformationToLeg(leg, output, i, j, mode, nb, superArc, legIndex);

            if(mode == BUS) {
                const auto& routes = output["design"]["legs"][leg["design_leg_index"].GetInt()]["routes"].GetArray();
                for (const auto& route : routes) {
                    leg["routes"].PushBack(rj::Value().SetString(route.GetString(), allocator), allocator);
                }
            }

            legs.PushBack(leg, allocator);

            trip["travel_time"] = trip["travel_time"].GetInt() + leg["travel_time"].GetInt();
            trip["waiting_time"] = trip["waiting_time"].GetInt() + leg["waiting_time"].GetInt();
            trip["total_time"] = trip["total_time"].GetInt() + leg["total_time"].GetInt();
            trip["distance"] = trip["distance"].GetInt() + leg["distance"].GetInt();
            trip["objective"] = trip["objective"].GetDouble() + leg["objective"].GetDouble();

        }

    }

}

void Runner::prepareTripSplittingJSON(rj::Document& output, const ListDigraph::ArcMap<int>& legIndex){

    auto& allocator = output.GetAllocator();

    output.AddMember("trip_splitting", rjObject, allocator);
    rj::Value& trip_splitting = output["trip_splitting"];

    if(!generateTripSplittings){
        trip_splitting = rjNull;
        return;
    }

    const pSuperGraph& superGraph = networkDesignProblem->superGraph;
    const ListDigraph& graph = superGraph->graph;
    const VstopInfo& allStops = stops->allStops;

//    //Legend
//    trip_splitting.AddMember("legend", rjObject, allocator);
//    trip_splitting["legend"].AddMember("b", "board_stop_id", allocator);
//    trip_splitting["legend"].AddMember("a", "alight_stop_id", allocator);
//    trip_splitting["legend"].AddMember("m", "mode", allocator);
//    trip_splitting["legend"].AddMember("i", "design_leg_index", allocator);

    //Initialize matrix
    for(const auto& originStop : allStops){

        int originID = std::get<0>(originStop);
        trip_splitting.AddMember(rjNumString(originID, allocator), rjObject, allocator);
        rj::Value& row = trip_splitting[rjNumString(originID, allocator)];

        for(const auto& destinationStop : allStops) {
            int destinationID = std::get<0>(destinationStop);
            row.AddMember(rjNumString(destinationID, allocator), rjNull, allocator);
        }

    }

    //Fill matrix
    for(const ListDigraph_VArc& path : tripSplittings){

        int origin = graph.id(graph.source(path.at(0)));
        int destination = graph.id(graph.target(path.at(path.size()-1)));

        int originID = std::get<0>(allStops.at(origin));
        int destinationID = std::get<0>(allStops.at(destination));

        rj::Value& cell = trip_splitting[rjNumString(originID, allocator)][rjNumString(destinationID, allocator)];
        cell.SetObject();

        cell.AddMember("legs", rjArray, allocator);
        rj::Value& legs = cell["legs"];

        for (int k = 0; k < path.size(); ++k) {

            const ListDigraph::Arc& superArc = path.at(k);
            const int& i = graph.id(graph.source(superArc));
            const int& j = graph.id(graph.target(superArc));
            const travel_mode& mode = (travel_mode) superGraph->properties[superArc].at(0);
            const int& nb = superGraph->properties[superArc].at(1);

            rj::Value leg(rj::kObjectType);
            addTripSplittingInformationToLeg(leg, output, i, j, mode, nb, superArc, legIndex);
            legs.PushBack(leg, allocator);

        }

    }

}

void Runner::prepareScoresJSON(rapidjson::Document& output) {

    auto& allocator = output.GetAllocator();

    output.AddMember("scores", rjObject, allocator);
    rj::Value& scores = output["scores"];

    const VTrip& allTrips = trips->allTrips;
    auto sumNbPeople = [](const int& sum, const Trip& trip){ return sum + std::get<2>(trip); };
    int totalNbTrips = allTrips.size();
    int totalNbPeople = std::accumulate(allTrips.begin(), allTrips.end(), 0, sumNbPeople);
    double correctedTotalNbPeople = passengerFactor * totalNbPeople;

    double designObjective, passengerObjective;
    Vdouble designCostConv, passengerCostConv;
    std::tie(designObjective, designCostConv, passengerObjective, passengerCostConv) = networkScores;

    double totalObjective = designObjective + passengerObjective;
    double totalCost = designCostConv.at(0) + passengerCostConv.at(0);
    double totalConvenience = designCostConv.at(1) + passengerCostConv.at(1);

    scores.AddMember("total_objective", totalObjective, allocator);
    scores.AddMember("total_cost", totalCost, allocator);
    scores.AddMember("total_convenience", totalConvenience, allocator);
    scores.AddMember("design_objective", designObjective, allocator);
    scores.AddMember("design_cost", designCostConv.at(0), allocator);
    scores.AddMember("design_convenience", designCostConv.at(1), allocator);
    scores.AddMember("passenger_objective", passengerObjective, allocator);
    scores.AddMember("passenger_cost", passengerCostConv.at(0), allocator);
    scores.AddMember("passenger_convenience", passengerCostConv.at(1), allocator);
    scores.AddMember("total_cost_per_passenger", totalCost/correctedTotalNbPeople, allocator);
    scores.AddMember("passenger_cost_per_passenger", passengerCostConv.at(0)/correctedTotalNbPeople, allocator);
    scores.AddMember("passenger_convenience_per_passenger", passengerCostConv.at(1)/correctedTotalNbPeople, allocator);
    scores.AddMember("original_number_of_passengers", totalNbPeople, allocator);
    scores.AddMember("scaled_number_of_passengers", correctedTotalNbPeople, allocator);

}

void Runner::writeJSONtoFile() {

    rj::Document output;
    output.SetObject();

    prepareParametersJSON(output);
    ListDigraph::ArcMap<int> legIndex(networkDesignProblem->superGraph->graph, -1);
    prepareDesignJSON(output, legIndex);
    prepareRoutesJSON(output);
    prepareTripsJSON(output, legIndex);
    prepareTripSplittingJSON(output, legIndex);
    prepareScoresJSON(output);

    rj::StringBuffer buffer;

    if(!generateTripSplittings) {
        rj::PrettyWriter<rj::StringBuffer> writer(buffer);
        output.Accept(writer);
    }
    else{
        rj::Writer<rj::StringBuffer> writer(buffer);
        output.Accept(writer);
    }

    std::ofstream outputStream;
    outputStream.open(outputFile);
    outputStream << buffer.GetString();
    outputStream.close();

}

void Runner::printTripStatistics() {

    std::cout << std::endl;
    std::cout << "===============" << std::endl;
    std::cout << "TRIP STATISTICS" << std::endl;
    std::cout << "===============" << std::endl;

    const VTrip& allTrips = trips->allTrips;
    auto sumNbPeople = [](const int& sum, const Trip& trip){ return sum + std::get<2>(trip); };
    int totalNbTrips = allTrips.size();
    int totalNbPeople = std::accumulate(allTrips.begin(), allTrips.end(), 0, sumNbPeople);

    std::cout << "Total number of trips: " << totalNbTrips << std::endl;
    std::cout << "Total number of people: " << totalNbPeople << std::endl;

}

void Runner::addGeneralInformationToLeg(rj::Value& leg, rj::Document& output,
                                        const int i, const int j, const travel_mode mode, const int nb){

    auto& allocator = output.GetAllocator();

    int iID = std::get<0>(stops->allStops.at(i));
    int jID = std::get<0>(stops->allStops.at(j));

    if(iID != jID && instance->arcModeExists(i, j, mode)){

        leg.AddMember("board_stop_id", rjNumString(iID, allocator), allocator);
        leg.AddMember("alight_stop_id", rjNumString(jID, allocator), allocator);
        leg.AddMember("mode", (int) mode, allocator);
        leg.AddMember("frequency_per_hour", rjObject, allocator);
        leg.AddMember("frequency_per_hour_index", rjObject, allocator);
        leg.AddMember("routes", rjArray, allocator);

        //frequency_per_hour
        if (mode == BUS || mode == RAIL) {
            leg["frequency_per_hour"] = nb/((double) timeHorizonInHours);
        }
        else if(mode == SHUTTLE){
            leg["frequency_per_hour"] = 0;
        }

        //frequency_per_hour_index
        if (mode == BUS) {
            auto iter = std::find(uniqueNbPossibleForBus.begin(), uniqueNbPossibleForBus.end(), nb);
            int index = std::distance(uniqueNbPossibleForBus.begin(), iter);
            leg["frequency_per_hour_index"] = index;
        }
        else{
            leg["frequency_per_hour_index"] = rjNull;
        }

        //routes
        if(mode == RAIL){
            int iID = std::get<0>(stops->allStops.at(i));
            int jID = std::get<0>(stops->allStops.at(j));

            std::tuple<double, int, Vint> lineInfo =
                    instance->getTimeNumberPerHourAndLineIndicesFromSchedule(*railLines, iID, jID, ignoreFrequencyCorrection);

            Vint lineIndices = std::get<2>(lineInfo);

            for(const int& lineIndex : lineIndices){
                std::string name = std::get<0>(railLines->allLines.at(lineIndex));
                leg["routes"].PushBack(rjString(name, allocator), allocator);
            }
        }
        else if(mode == BUS){
            //keep initialized at empty array
        }
        else{
            leg["routes"] = rjNull;
        }

    }

}

void Runner::addDesignInformationToLeg(rj::Value& leg, rj::Document& output,
                                       const int i, const int j, const travel_mode mode, const int nb){

    auto& allocator = output.GetAllocator();

    int iID = std::get<0>(stops->allStops.at(i));
    int jID = std::get<0>(stops->allStops.at(j));

    if(iID != jID && instance->arcModeExists(i, j, mode)){
        leg.AddMember("fixed", mode == RAIL, allocator);
        leg.AddMember("travel_time", min2sec(instance->travelTime(i, j, mode)), allocator);
        leg.AddMember("waiting_time", waitingTimeSec(mode, nb), allocator);
        leg.AddMember("total_time", leg["travel_time"].GetInt() + leg["waiting_time"].GetInt(), allocator);
        leg.AddMember("distance", km2m(instance->travelDistance(i, j, mode)), allocator);
        leg.AddMember("objective", instance->arcObjective(i, j, mode, nb, Instance::DESIGN), allocator);
        leg.AddMember("cost", instance->arcCost(i, j, mode, nb, Instance::DESIGN), allocator);
    }

}

void Runner::addTripInformationToLeg(rj::Value& leg, rj::Document& output,
                                     const int i, const int j, const travel_mode mode, const int nb,
                                     const ListDigraph::Arc& superArc,
                                     const ListDigraph::ArcMap<int>& legIndex){

    auto& allocator = output.GetAllocator();

    int iID = std::get<0>(stops->allStops.at(i));
    int jID = std::get<0>(stops->allStops.at(j));

    if(iID != jID && instance->arcModeExists(i, j, mode)){

        leg.AddMember("travel_time", min2sec(instance->travelTime(i, j, mode)), allocator);
        leg.AddMember("waiting_time", waitingTimeSec(mode, nb), allocator);
        leg.AddMember("total_time", leg["travel_time"].GetInt() + leg["waiting_time"].GetInt(), allocator);
        leg.AddMember("distance", km2m(instance->travelDistance(i, j, mode)), allocator);
        leg.AddMember("objective", instance->arcObjective(i, j, mode, nb, Instance::PASSENGER)/passengerFactor, allocator);
        leg.AddMember("cost", instance->arcCost(i, j, mode, nb, Instance::PASSENGER)/passengerFactor, allocator);
        leg.AddMember("convenience", instance->arcConvenience(i, j, mode, nb, Instance::PASSENGER)/passengerFactor, allocator);
        leg.AddMember("design_leg_index", rjObject, allocator);

        //design_leg_index
        int index = legIndex[superArc];
        if(index == -1){
            leg["design_leg_index"] = rjNull;
        }
        else{
            leg["design_leg_index"] = index;
        }

    }

}

void Runner::addTripSplittingInformationToLeg(rj::Value& leg, rj::Document& output,
                                     const int i, const int j, const travel_mode mode, const int nb,
                                     const ListDigraph::Arc& superArc,
                                     const ListDigraph::ArcMap<int>& legIndex){

    auto& allocator = output.GetAllocator();

    int iID = std::get<0>(stops->allStops.at(i));
    int jID = std::get<0>(stops->allStops.at(j));

    int index = legIndex[superArc];
    if(index == -1){
        leg.AddMember("b", rjNumString(iID, allocator), allocator);
        leg.AddMember("a", rjNumString(jID, allocator), allocator);
        leg.AddMember("m", (int) mode, allocator);
        if(mode != SHUTTLE){
            throw std::runtime_error("Exception: non-design arc found that is not of type SHUTTLE.");
        }
    }
    else{
        leg.AddMember("i", index, allocator);
    }

}

std::function<std::tuple<double, Vdouble, double, Vdouble>
        (const double, const double, const travel_mode, const int, const int)> Runner::arcScores() {

    return [this](const double time, const double distance, const travel_mode mode, const int nb, const int timeHorizonInHours) {

        //time in minutes
        //distance in km

        double designCost;
        double designConvenience;

        double passengerCost;           //for a single passenger
        double passengerConvenience;    //for a single passenger

        if(mode == SHUTTLE){
            designCost = 0;
            designConvenience = 0;
            passengerCost = passengerFactor * shuttleCostPerKm * distance;
            passengerConvenience = passengerFactor * time;
        }
        else if(mode == BUS){
            designCost = (busCostPerHour / 60.0) * nb * time;
            designConvenience = 0;
            passengerCost = passengerFactor * 0;
            passengerConvenience = passengerFactor * (time + fixedTransferTime + (60.0 * timeHorizonInHours)/(2.0 * nb)); //TODO: only transfer time when transfering?
        }
        else if(mode == RAIL){
            designCost = 0; //note: currently all rail distances in Instance are 0
            designConvenience = 0;
            passengerCost = passengerFactor * 0;
            passengerConvenience = passengerFactor * (time + fixedTransferTime + (60.0 * timeHorizonInHours)/(2.0 * nb)); //TODO: only transfer time when transfering?
        }

        double designObjective = (1-alpha) * designCost + alpha * designConvenience;
        double passengerObjective = (1-alpha) * passengerCost + alpha * passengerConvenience;

        return std::make_tuple(designObjective, Vdouble{designCost, designConvenience}, passengerObjective, Vdouble{passengerCost, passengerConvenience});

    };

}