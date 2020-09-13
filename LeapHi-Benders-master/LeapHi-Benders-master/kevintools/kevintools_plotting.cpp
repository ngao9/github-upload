#include "kevintools_plotting.hpp"

#include <sstream>

std::string kevintools::rgb(const int r, const int g, const int b){
    int intCode = r * 65536 + g * 256 + b;
    std::stringstream stream;
    stream << "#" << std::hex << intCode;
    return stream.str();
}

std::string kevintools::argb(const int alpha, const int r, const int g, const int b){
    int intCode = alpha * 16777216 + r * 65536 + g * 256 + b;
    std::stringstream stream;
    stream << "#" << std::hex << intCode;
    return stream.str();
}

std::string kevintools::argb(const Vint& argb){
    return kevintools::argb(argb.at(0), argb.at(1), argb.at(2), argb.at(3));
}

void kevintools::scatter(GnuplotPipe& gp, const VVdouble& data, const bool replot, const int pointtype, const std::string color){

    const int& n = data.size();

    for(int i = 0; i < n; ++i){

        std::string x = std::to_string(data.at(i).at(0));
        std::string y = std::to_string(data.at(i).at(1));

        gp.sendLine(x + ", " + y, true);

    }

    gp.writeBufferToFile("kevin_plotting_" + std::to_string(gp.counter) + ".dat");
    gp.sendEndOfData(0);

    std::string command = "'kevin_plotting_" + std::to_string(gp.counter) + ".dat'"
                   " with points pointtype " + std::to_string(pointtype) +
                   " lt rgb '" + color + "'";

    if(replot){
        gp.sendLine("replot " + command);
    } else{
        gp.sendLine("plot " + command);
    }

    ++gp.counter;

}

void kevintools::lines(GnuplotPipe& gp, const VVdouble& from, const VVdouble& to, const bool replot, const int pointtype, const std::string& color, const int linewidth){

    const int& n = from.size();

    for(int i = 0; i < n; ++i){

        std::string xFrom = std::to_string(from.at(i).at(0));
        std::string yFrom = std::to_string(from.at(i).at(1));
        std::string xTo = std::to_string(to.at(i).at(0));
        std::string yTo = std::to_string(to.at(i).at(1));

        gp.sendLine(xFrom + ", " + yFrom, true);
        gp.sendLine(xTo + ", " + yTo, true);
        gp.sendNewDataBlock();

    }

    gp.writeBufferToFile("kevin_plotting_" + std::to_string(gp.counter) + ".dat");
    gp.sendEndOfData(0);

    std::string command = "'kevin_plotting_" + std::to_string(gp.counter) + ".dat'"
                          " with lines pointtype " + std::to_string(pointtype) +
                          " lt rgb '" + color + "'"
                          " lw " + std::to_string(linewidth);

    if(replot){
        gp.sendLine("replot " + command);
    } else{
        gp.sendLine("plot " + command);
    }

    ++gp.counter;

}

void kevintools::arrows(GnuplotPipe& gp, const VVdouble& from, const VVdouble& to, const bool replot, const std::string& color, const int linewidth){

    const int& n = from.size();

    for(int i = 0; i < n; ++i){

        std::string xFrom = std::to_string(from.at(i).at(0));
        std::string yFrom = std::to_string(from.at(i).at(1));
        std::string xDiff = std::to_string(to.at(i).at(0) - from.at(i).at(0));
        std::string yDiff = std::to_string(to.at(i).at(1) - from.at(i).at(1));

        gp.sendLine(xFrom + ", " + yFrom + ", " + xDiff + ", " + yDiff + ", ", true);
        gp.sendNewDataBlock();

    }

    gp.writeBufferToFile("kevin_plotting_" + std::to_string(gp.counter) + ".dat");
    gp.sendEndOfData(0);

    std::string command = "'kevin_plotting_" + std::to_string(gp.counter) + ".dat'"
                            " using 1:2:3:4"
                            " with vectors"
                            " front head filled size screen 0.02,15"
                            " lt rgb '" + color + "'"
                            " lw " + std::to_string(linewidth);

    if(replot){
        gp.sendLine("replot " + command);
    } else{
        gp.sendLine("plot " + command);
    }

    ++gp.counter;

}