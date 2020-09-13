#pragma once

#include <vector>
#include <sstream>
#include <iterator>

namespace kevintools{

    inline std::string to_string(const std::vector<int>& v){

        std::stringstream ssResult;
        std::copy(v.begin(), v.end(), std::ostream_iterator<int>(ssResult, ", "));
        std::string result = ssResult.str();
        if (v.size() > 0) {
            result.erase(result.end() - 2, result.end());
        }

        return "[ " + result + "]";

    }

    inline std::string to_string(const std::vector<double>& v) {

        std::stringstream ssResult;
        std::copy(v.begin(), v.end(), std::ostream_iterator<double>(ssResult, ", "));
        std::string result = ssResult.str();
        if (v.size() > 0) {
            result.erase(result.end() - 2, result.end());
        }

        return "[ " + result + "]";

    }

    inline std::string to_string(const std::vector<std::string>& v) {

        std::stringstream ssResult;
        std::copy(v.begin(), v.end(), std::ostream_iterator<std::string>(ssResult, ", "));
        std::string result = ssResult.str();
        if (v.size() > 0) {
            result.erase(result.end() - 2, result.end());
        }

        return "[ " + result + "]";

    }

    inline std::vector<int> sizes(const std::vector<std::vector<int> >& v){
        std::vector<int> result(v.size(), 0);
        for (size_t i = 0; i < v.size(); ++i){
            result.at(i) = v.at(i).size();
        }
        return result;
    }

    void removeCharsFromString(std::string &str, const std::string &charsToRemove);
    void trimWhitespace(std::string& str);

}