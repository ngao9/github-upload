#include "networkdesign_cpp_include.hpp"

#include "kevintools_vectors.hpp"

#include "gtest/gtest.h"

using namespace networkdesign;

VVdouble railFrequencyCorrectionTestCases(){

    return VVdouble{
        {1001, 1002, 5, 13},
        {1001, 1003, 8.7, 3},
        {1001, 1004, 11.3, 10},
        {1001, 1005, 12.1, 13},
        {1001, 1006, 14.1, 7},
        {1001, 1007, 16.1, 6},
        {1002, 1001, 5, 13},
        {1002, 1003, 3.7, 3},
        {1002, 1004, 6.3, 10},
        {1002, 1005, 7.1, 13},
        {1002, 1006, 9.1, 7},
        {1002, 1007, 11.1, 6},
        {1003, 1001, 7.7, 3},
        {1003, 1002, 2.7, 3},
        {1003, 1005, 3.4, 3},
        {1003, 1006, 5.4, 3},
        {1004, 1001, 10.3, 10},
        {1004, 1002, 5.3, 10},
        {1004, 1005, 0.8, 10},
        {1004, 1006, 2.8, 4},
        {1004, 1007, 4.8, 6},
        {1005, 1001, 11.1, 13},
        {1005, 1002, 6.1, 13},
        {1005, 1003, 3.4, 3},
        {1005, 1004, 0.8, 10},
        {1005, 1006, 2, 7},
        {1005, 1007, 4, 6},
        {1006, 1001, 13.1, 7},
        {1006, 1002, 8.1, 7},
        {1006, 1003, 5.4, 3},
        {1006, 1004, 2.8, 4},
        {1006, 1005, 2, 7},
        {1007, 1001, 15.1, 6},
        {1007, 1002, 10.1, 6},
        {1007, 1004, 4.8, 6},
        {1007, 1005, 4, 6},
        {1009, 1003, 2.6, 5},
        {1009, 1004, 11.6, 5},
        {1009, 1008, 14.6, 5},
        {1003, 1009, 2.6, 5},
        {1003, 1004, 9, 5},
        {1003, 1008, 12, 5},
        {1004, 1009, 11.6, 5},
        {1004, 1003, 9, 5},
        {1004, 1008, 3, 5},
        {1008, 1009, 14.6, 5},
        {1008, 1003, 12, 5},
        {1008, 1004, 3, 5}
    };

}

VVdouble railFrequencyIgnoreCorrectionTestCases(){

    return VVdouble{
            {1001, 1002, 5, 6},
            {1001, 1003, 8.7, 3},
            {1001, 1004, 11.3, 6},
            {1001, 1005, 12.1, 6},
            {1001, 1006, 14.1, 4},
            {1001, 1007, 16.1, 6},
            {1002, 1001, 5, 6},
            {1002, 1003, 3.7, 3},
            {1002, 1004, 6.3, 6},
            {1002, 1005, 7.1, 6},
            {1002, 1006, 9.1, 4},
            {1002, 1007, 11.1, 6},
            {1003, 1001, 7.7, 3},
            {1003, 1002, 2.7, 3},
            {1003, 1005, 3.4, 3},
            {1003, 1006, 5.4, 3},
            {1004, 1001, 10.3, 6},
            {1004, 1002, 5.3, 6},
            {1004, 1005, 0.8, 6},
            {1004, 1006, 2.8, 4},
            {1004, 1007, 4.8, 6},
            {1005, 1001, 11.1, 6},
            {1005, 1002, 6.1, 6},
            {1005, 1003, 3.4, 3},
            {1005, 1004, 0.8, 6},
            {1005, 1006, 2, 4},
            {1005, 1007, 4, 6},
            {1006, 1001, 13.1, 4},
            {1006, 1002, 8.1, 4},
            {1006, 1003, 5.4, 3},
            {1006, 1004, 2.8, 4},
            {1006, 1005, 2, 4},
            {1007, 1001, 15.1, 6},
            {1007, 1002, 10.1, 6},
            {1007, 1004, 4.8, 6},
            {1007, 1005, 4, 6},
            {1009, 1003, 2.6, 5},
            {1009, 1004, 11.6, 5},
            {1009, 1008, 14.6, 5},
            {1003, 1009, 2.6, 5},
            {1003, 1004, 9, 5},
            {1003, 1008, 12, 5},
            {1004, 1009, 11.6, 5},
            {1004, 1003, 9, 5},
            {1004, 1008, 3, 5},
            {1008, 1009, 14.6, 5},
            {1008, 1003, 12, 5},
            {1008, 1004, 3, 5}
    };

}

Stops constructFixedRailInstance(Instance& in, const std::string& railLinesPath, const int timeHorizonInHours, const bool ignoreFrequencyCorrection){

    Stops stops("../test/leaphi/dataprocessing_testdata/stops.csv");
    RailLines railLines(railLinesPath);
    Trips trips("../test/leaphi/dataprocessing_testdata/odx.csv");
    Distances roadTravelTimes("../test/leaphi/dataprocessing_testdata/euclidean_distance_matrix.csv", stops);
    VAdouble2 desiredHubLatLons = VAdouble2{   Adouble2{33.570107, -84.547634},
                                               Adouble2{33.708081, -84.117133},
                                               Adouble2{33.869261, -84.2246},
                                               Adouble2{33.69298, -84.401602},
                                               Adouble2{33.560939, -84.360986},
                                               Adouble2{33.996391, -84.349409},
                                               Adouble2{33.626606, -84.374352},
                                               Adouble2{33.69597, -84.273397},
                                               Adouble2{33.818697, -84.196099},
                                               Adouble2{34.068248, -84.287365}  };

    double mile2km = 1.609344;
    in.constructFixedRailInstance(  stops, railLines,
                                    trips, timeHorizonInHours,
                                    Vint{12, 24},
                                    roadTravelTimes, roadTravelTimes,
                                    1, 1,
                                    desiredHubLatLons,
                                    -1,
                                    ignoreFrequencyCorrection,
                                    false);

    return stops;

}

TEST(Instance, railFrequencyCorrection){

    std::string path = "../test/leaphi/dataprocessing_testdata/rail_lines.csv";

    Instance in = Instance();
    int timeHorizonInHours = 3;
    bool ignoreFrequencyCorrection = false;
    Stops stops = constructFixedRailInstance(in, path, timeHorizonInHours, ignoreFrequencyCorrection);

    VVdouble testCases = railFrequencyCorrectionTestCases();

    for(const Vdouble& testCase : testCases){

        int fromId = (int) std::round(testCase.at(0));
        int toId = (int) std::round(testCase.at(1));
        double time = testCase.at(2);
        int nb = (int) std::round(testCase.at(3));

        int i = stops.idToIndex(fromId);
        int j = stops.idToIndex(toId);

        EXPECT_NEAR(in.travelTime(i, j, RAIL), time, 0.000001);
        EXPECT_EQ(in.nbPossible(i, j, RAIL), Vint{3 * nb});

    }

    int numberRailArcs = 0;
    for(int i = 0; i < in.n; ++i){
        for(int j = 0; j < in.n; ++j){
            numberRailArcs += in.arcModeExists(i, j, RAIL);
        }
    }

    EXPECT_EQ(numberRailArcs, testCases.size());

}

TEST(Instance, railFrequencyCorrectionExceptions){

    std::string path = "../test/leaphi/dataprocessing_testdata/rail_lines_inconsistent_travel.csv";

    Instance in = Instance();
    int timeHorizonInHours = 3;
    bool ignoreFrequencyCorrection = false;
    EXPECT_ANY_THROW(constructFixedRailInstance(in, path, timeHorizonInHours, ignoreFrequencyCorrection));

}

TEST(Instance, railFrequencyIgnoreCorrection){

    std::string path = "../test/leaphi/dataprocessing_testdata/rail_lines.csv";

    Instance in = Instance();
    int timeHorizonInHours = 3;
    bool ignoreFrequencyCorrection = true;
    Stops stops = constructFixedRailInstance(in, path, timeHorizonInHours, ignoreFrequencyCorrection);

    VVdouble testCases = railFrequencyIgnoreCorrectionTestCases();

    for(const Vdouble& testCase : testCases){

        int fromId = (int) std::round(testCase.at(0));
        int toId = (int) std::round(testCase.at(1));
        double time = testCase.at(2);
        int nb = (int) std::round(testCase.at(3));

        int i = stops.idToIndex(fromId);
        int j = stops.idToIndex(toId);

        EXPECT_NEAR(in.travelTime(i, j, RAIL), time, 0.000001);
        EXPECT_EQ(in.nbPossible(i, j, RAIL), Vint{3 * nb});

    }

    int numberRailArcs = 0;
    for(int i = 0; i < in.n; ++i){
        for(int j = 0; j < in.n; ++j){
            numberRailArcs += in.arcModeExists(i, j, RAIL);
        }
    }

    EXPECT_EQ(numberRailArcs, testCases.size());

}

TEST(Instance, uniqueNbPossibleForMode){

    std::string path = "../test/leaphi/dataprocessing_testdata/rail_lines.csv";

    Instance in = Instance();
    int timeHorizonInHours = 3;
    bool ignoreFrequencyCorrection = false;
    Stops stops = constructFixedRailInstance(in, path, timeHorizonInHours, ignoreFrequencyCorrection);

    Vint hardcoded = kevintools::operator*(Vint{3, 4, 5, 6, 7, 10, 13}, in.timeHorizonInHours);
    EXPECT_EQ(in.uniqueNbPossibleForMode(RAIL), hardcoded);

}