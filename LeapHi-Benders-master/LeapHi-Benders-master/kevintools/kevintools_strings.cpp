#include "kevintools_strings.hpp"

#include <cstring>
#include <algorithm>

void kevintools::removeCharsFromString(std::string &str, const std::string& charsToRemove) {
    for(int i = 0; i < charsToRemove.length(); ++i ) {
        str.erase(std::remove(str.begin(), str.end(), charsToRemove[i]), str.end());
    }
}

void kevintools::trimWhitespace(std::string& str){
    const std::string& whitespace = " \t";
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
       str = ""; // no content
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    str = str.substr(strBegin, strRange);
}