/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * URL: https://github.com/martinruenz/gnuplot-cpp
 * AUTHOR: Martin Rünz, 2015
 */

#pragma once

#define STRING2(x) #x
#define STRING(x) STRING2(x)

#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

class GnuplotPipe {
public:
    inline GnuplotPipe(bool persist = true, bool noplot = false, std::string gnuplotPath = std::string(STRING(GNUPLOT_PATH)) + "/gnuplot") :
        noplot(noplot) {
        if(noplot){
            std::cerr << "Guplot plotting disabled." << std::endl;
            return;
        }
        std::cout << "Opening gnuplot... ";
        #ifdef _WIN32
                std::string command1 = "\"" + gnuplotPath + "\" -persist";
                std::string command2 = "\"" + gnuplotPath + "\"";
                pipe = _popen(persist ? command1.c_str() : command2.c_str(), "w");
        #else
                pipe = popen(persist ? "gnuplot -persist" : "gnuplot", "w");
        #endif
        if (!pipe)
            std::cout << "failed!" << std::endl;
        else
            std::cout << "succeded." << std::endl;

        counter = 0;
    }
    inline virtual ~GnuplotPipe(){
        if(noplot) return;
        #ifdef _WIN32
                if (pipe) _pclose(pipe);
        #else
                if (pipe) pclose(pipe);
        #endif
    }

    void sendLine(const std::string& text, bool useBuffer = false){
        if(noplot) return;
        if (!pipe) return;
        if (useBuffer)
            buffer.push_back(text + "\n");
        else
            fputs((text + "\n").c_str(), pipe);
    }
    void sendEndOfData(unsigned repeatBuffer = 1){
        if(noplot) return;
        if (!pipe) return;
        for (unsigned i = 0; i < repeatBuffer; i++) {
            for (auto& line : buffer) fputs(line.c_str(), pipe);
            fputs("e\n", pipe);
        }
        fflush(pipe);
        buffer.clear();
    }
    void sendNewDataBlock(){
        if(noplot) return;
        sendLine("\n", !buffer.empty());
    }

    void writeBufferToFile(const std::string& fileName){
        if(noplot) return;
        std::ofstream fileOut(fileName);
        for (auto& line : buffer) fileOut << line;
        fileOut.close();
    }

private:
    GnuplotPipe(GnuplotPipe const&) = delete;
    void operator=(GnuplotPipe const&) = delete;

    FILE* pipe;
    std::vector<std::string> buffer;
    const bool noplot;

public:
    int counter;
};
