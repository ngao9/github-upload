#pragma once

#include <vector>
#include <array>
#include <string>

namespace kevintools{

    double euclideanDistance(const std::vector<double>& x, const std::vector<double>& y);

    class Aint2_hasher {
    public:
        std::size_t operator()(const std::array<int,2>& vec) const {
            std::size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    class Aint3_hasher {
    public:
        std::size_t operator()(const std::array<int,3>& vec) const {
            std::size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    class Aint4_hasher {
    public:
        std::size_t operator()(const std::array<int,4>& vec) const {
            std::size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    class Aint5_hasher {
    public:
        std::size_t operator()(const std::array<int,5>& vec) const {
            std::size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str);

    double deg2rad(double deg);
    double rad2deg(double rad);
    double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d);

}