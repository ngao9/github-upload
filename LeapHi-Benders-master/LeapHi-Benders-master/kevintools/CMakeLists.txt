 project(leaphi)

 SET(SOURCE_FILES_SRC kevintools_strings.cpp kevintools_primitives.cpp kevintools_functions.cpp kevintools_plotting.cpp kevintools_vectors.cpp kevintools_typedefs.hpp kevintools_gurobi.cpp kevintools_cplex.cpp kevintools_time.cpp kevintools_maps.cpp kevintools_lemon.cpp)

 add_library(kevintools_lib STATIC ${SOURCE_FILES_SRC})

 target_link_libraries(kevintools_lib ${GUROBI_LIBRARIES})
 target_link_libraries(kevintools_lib ${LEMON_LINK_LIBRARIES})