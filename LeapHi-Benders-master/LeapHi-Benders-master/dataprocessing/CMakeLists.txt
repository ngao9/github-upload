 project(leaphi)

 SET(SOURCE_FILES_SRC stops.cpp raillines.cpp trips.cpp trips.hpp dataprocessing_typedefs.hpp dataprocessing_cpp_include.hpp distances.cpp)

 add_library(dataprocessing_lib STATIC ${SOURCE_FILES_SRC})
 target_link_libraries(kevintools_lib)
