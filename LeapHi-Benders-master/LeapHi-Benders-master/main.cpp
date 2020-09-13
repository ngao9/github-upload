#include "networkdesign_cpp_include.hpp"

#include "kevintools_vectors.hpp"

#include <ctime>

using namespace networkdesign;
using namespace dataprocessing;

int main(int argc, char* argv[]){

    std::clock_t start = std::clock();
    Runner runner(argc, argv, start);
    runner.writeJSONtoFile();

}