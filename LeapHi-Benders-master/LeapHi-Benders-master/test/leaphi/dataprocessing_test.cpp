#include "dataprocessing_cpp_include.hpp"

#include "kevintools_strings.hpp"
#include "kevintools_vectors.hpp"

#include <cmath>

#include "gtest/gtest.h"

using namespace dataprocessing;

void printStop(const Stop& stop){
    std::cout <<    "Stop: id = " << std::get<0>(stop) << ", " <<
                    "latlong = [" << std::get<1>(stop).at(0) << ", " << std::get<1>(stop).at(1) << "], " <<
                    "route = " << kevintools::to_string(std::get<2>(stop)) << "." << std::endl;
}

bool contains(const VStop& vector, const Stop& stop, double eps){

    for(const Stop& element : vector){

        if( std::get<0>(element) == std::get<0>(stop) &&
            std::get<2>(element) == std::get<2>(stop)){

            Adouble2 latlongElement = std::get<1>(element);
            Adouble2 latlongStop = std::get<1>(stop);

            if( std::fabs(latlongElement.at(0) - latlongStop.at(0)) < eps &&
                std::fabs(latlongElement.at(1) - latlongStop.at(1)) < eps){
                return true;
            }

        }

    }

    printStop(stop);
    return false;

}

bool contains(const VRailLine& vector, const RailLine& railLine, double eps){

    for(const RailLine& element : vector){

        if( std::get<0>(element) == std::get<0>(railLine) &&
            std::get<1>(element) == std::get<1>(railLine) &&
            std::get<2>(element) == std::get<2>(railLine)){

            const Vdouble& fwdElement = std::get<3>(element);
            const Vdouble& fwdRailLine = std::get<3>(railLine);

            const Vdouble& bwdElement = std::get<4>(element);
            const Vdouble& bwdRailLine = std::get<4>(railLine);

            if( kevintools::almostEqual(fwdElement, fwdRailLine, eps) &&
                    kevintools::almostEqual(bwdElement, bwdRailLine, eps)){
                return true;
            }

        }

    }

    std::cout << "RailLine " << std::get<0>(railLine) << " could not be found." << std::endl;
    return false;

}

bool contains(const VTrip& vector, const Trip& trip) {

    for(const Trip& element : vector){
        if(element == trip){
            return true;
        }
    }

    std::cout << "Trip [" <<    std::get<0>(trip) << ", " <<
                                std::get<1>(trip) << ", " <<
                                std::get<2>(trip) << ", " <<
                                kevintools::to_string(std::get<3>(trip)) <<
                                " could not be found." << std::endl;
    return false;

}

TEST(Stops, readStops){

    std::string stopPath = "../test/leaphi/dataprocessing_testdata/stops.csv";

    EXPECT_ANY_THROW(Stops(stopPath + "y"));

    Stops stops(stopPath);
    VStop allStops = stops.allStops;

    ASSERT_EQ(allStops.size(), 25 + 9);

    VStop hardcodedStops = {
            Stop{-13, Adouble2{6.17321604767598, 44.5820051843649}, {""}},
            Stop{-14, Adouble2{92.3705757514448, 58.9906143286085}, {""}},
            Stop{-27, Adouble2{83.8480510643171, 65.5400955821967}, {""}},
            Stop{-12, Adouble2{34.4047859264227, 25.2930321321321}, {""}},
            Stop{-6, Adouble2{73.816965108346, 37.4040316059738}, {""}},
            Stop{-31, Adouble2{83.3980958959282, 79.2712298184782}, {""}},
            Stop{-19, Adouble2{93.14878122826, 28.7132182096138}, {""}},
            Stop{-26, Adouble2{39.7930755668276, 94.1092590070672}, {""}},
            Stop{-22, Adouble2{65.6080806850356, 20.4238291233626}, {""}},
            Stop{-7, Adouble2{98.9811831250573, 18.062545897057}, {""}},
            Stop{-1, Adouble2{65.7270289916874, 18.2829135195354}, {""}},
            Stop{-20, Adouble2{21.6083453199187, 80.8314777403796}, {""}},
            Stop{-29, Adouble2{44.2097068893388, 62.1644090507262}, {""}},
            Stop{-4, Adouble2{81.0342923094372, 94.4071844089071}, {""}},
            Stop{-16, Adouble2{49.9337313765204, 31.2079754071943}, {""}},
            Stop{-5, Adouble2{82.9750219559047, 45.4492751621573}, {""}},
            Stop{-21, Adouble2{31.2432818877265, 69.9293755212409}, {""}},
            Stop{-25, Adouble2{45.6398586654972, 35.9112368473649}, {""}},
            Stop{-23, Adouble2{19.4301420881986, 78.2592891170099}, {""}},
            Stop{-9, Adouble2{27.0936147489554, 75.7995133511119}, {""}},
            Stop{-15, Adouble2{73.7454334430732, 31.2274045957946}, {""}},
            Stop{-10, Adouble2{28.4911288432609, 18.3152880048963}, {""}},
            Stop{-11, Adouble2{53.6197542897234, 55.3756203515614}, {""}},
            Stop{-17, Adouble2{72.1806630828485, 95.4565190680976}, {""}},
            Stop{-30, Adouble2{50.2404911246121, 90.7006014924693}, {""}},
            Stop{1008, Adouble2{30.061476266176, 50.0007637028585}, {""}},
            Stop{1003, Adouble2{40.1907420536788, 50.902612966929}, {""}},
            Stop{1007, Adouble2{60.2309554320684, 69.9314927797603}, {""}},
            Stop{1006, Adouble2{40.0101358163779, 69.1642858473511}, {""}},
            Stop{1009, Adouble2{70.0790416736001, 49.3303608573783}, {""}},
            Stop{1002, Adouble2{49.7513639238148, 40.3811034823089}, {""}},
            Stop{1005, Adouble2{50.2581326495969, 60.6360552932413}, {""}},
            Stop{1004, Adouble2{60.0104012715822, 50.6123322746007}, {""}},
            Stop{1001, Adouble2{50.7232821129649, 29.1007648951054}, {""}}
    };

    for(const Stop& hardcodedStop : hardcodedStops){
        EXPECT_TRUE(contains(allStops, hardcodedStop, 0.000001));
    }

}

TEST(RailLines, readRailLines){

    std::string path = "../test/leaphi/dataprocessing_testdata/rail_lines.csv";

    EXPECT_ANY_THROW(RailLines(path + "y"));

    RailLines railLines(path);
    VRailLine allLines = railLines.allLines;

    ASSERT_EQ(allLines.size(), 4);

    VRailLine hardcodedRailLines = {
            RailLine{"GOLD", 3, Vint{1001, 1002, 1003, 1005, 1006}, Vdouble{0, 5, 8.7, 12.1, 14.1}, Vdouble{13.1, 8.1, 5.4, 2, 0}},
            RailLine{"RED", 4, Vint{1006, 1005, 1004, 1002, 1001}, Vdouble{0, 2, 2.8, 8.1, 13.1}, Vdouble{14.1, 12.1, 11.3, 5, 0}},
            RailLine{"GREEN", 5, Vint{1009, 1003, 1004, 1008}, Vdouble{0, 2.6, 11.6, 14.6}, Vdouble{14.6, 12, 3, 0}},
            RailLine{"BLUE", 6, Vint{1001, 1002, 1004, 1005, 1007}, Vdouble{0, 5, 11.3, 12.1, 16.1}, Vdouble{15.1, 10.1, 4.8, 4, 0}}
    };

    for(const RailLine& hardcodedRailLine : hardcodedRailLines){
        EXPECT_TRUE(contains(allLines, hardcodedRailLine, 0.000001));
    }

}

TEST(Distances, readEuclideanDistances){

    std::string stopPath = "../test/leaphi/dataprocessing_testdata/stops.csv";
    Stops stops(stopPath);
    const VStop& allStops = stops.allStops;
    const int& n = allStops.size();

    std::string path = "../test/leaphi/dataprocessing_testdata/euclidean_distance_matrix.csv";
    EXPECT_ANY_THROW(Distances(path + "y", stops));

    Distances euclideanDistances = Distances(path, stops);

    VVdouble calculatedDistances = VVdouble(n, Vdouble(n, 0));

    for(int i = 0; i < n; ++i){

        const double& x1 = std::get<1>(allStops.at(i)).at(0);
        const double& x2 = std::get<1>(allStops.at(i)).at(1);

        for(int j = 0; j < n; ++j){

            const double& y1 = std::get<1>(allStops.at(j)).at(0);
            const double& y2 = std::get<1>(allStops.at(j)).at(1);

            calculatedDistances.at(i).at(j) = std::sqrt((x1-y1)*(x1-y1) + (x2-y2)*(x2-y2));

        }
    }

    EXPECT_TRUE(kevintools::almostEqual(calculatedDistances, euclideanDistances.allDistances, 0.000001));

}

TEST(Distances, missingValues){

    std::string stopPath = "../test/leaphi/dataprocessing_testdata/stops.csv";
    Stops stops(stopPath);
    const VStop& allStops = stops.allStops;
    const int& n = allStops.size();

    std::string path = "../test/leaphi/dataprocessing_testdata/euclidean_distance_matrix_missing.csv";
    EXPECT_ANY_THROW(Distances(path + "y", stops));

    Distances euclideanDistances = Distances(path, stops);

    for(int i = 0; i < n; ++i){

        const int& fromId = std::get<0>(allStops.at(i));

        for(int j = 0; j < n; ++j){

            const int& toId = std::get<0>(allStops.at(j));

            if(fromId == -11 && toId == -29){
                EXPECT_TRUE(std::isnan(euclideanDistances.allDistances.at(i).at(j)));
            }
            else{
                EXPECT_FALSE(std::isnan(euclideanDistances.allDistances.at(i).at(j)));
            }

        }
    }

}

TEST(Trips, readTrips){

    std::string path = "../test/leaphi/dataprocessing_testdata/odx.csv";

    EXPECT_ANY_THROW(Trips(path + "y"));

    Trips trips(path);
    VTrip allTrips = trips.allTrips;

    ASSERT_EQ(allTrips.size(), 3);

    VTrip hardcodedTrips = {
            Trip{-7, -20, 2, {"2018-04-03 08:31:38", "2018-04-03 09:36:41"}},
            Trip{-10, -4, 1, {"2018-04-03 09:25:01"}},
            Trip{-13, 1005, 1, {"2018-04-03 07:30:48"}}
    };

    for(const Trip& hardcodedTrip : hardcodedTrips){
        EXPECT_TRUE(contains(allTrips, hardcodedTrip));
    }

    std::string path2 = "../test/leaphi/dataprocessing_testdata/odx_non_aggregated.csv";
    EXPECT_ANY_THROW(Trips trips2(path2));

    std::string path3 = "../test/leaphi/dataprocessing_testdata/odx_no_destination.csv";
    EXPECT_ANY_THROW(Trips trips3(path3));

}