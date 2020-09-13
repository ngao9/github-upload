 project(leaphi)

 SET(SOURCE_FILES_SRC networksolver_typedefs.hpp networksolver_cpp_include.hpp network.cpp supergraph.cpp lemonnetwork.cpp)
 if(GUROBI_FOUND)
  SET(SOURCE_FILES_SRC ${SOURCE_FILES_SRC} gurobinetwork.cpp)
 endif()
 add_library(networksolver_lib STATIC ${SOURCE_FILES_SRC})
 target_link_libraries(networksolver_lib kevintools_lib)