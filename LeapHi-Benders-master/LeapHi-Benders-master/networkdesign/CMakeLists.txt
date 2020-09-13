 project(leaphi)

 SET(SOURCE_FILES_SRC instance.cpp networkdesignproblem.cpp callbacks.cpp networkdesign_typedefs.hpp networkdesign_cpp_include.hpp runner.cpp)

 add_library(networkdesign_lib STATIC ${SOURCE_FILES_SRC})

 target_link_libraries(networkdesign_lib networksolver_lib kevintools_lib dataprocessing_lib ${requiredlibs} ${OpenMP_LIBRARIES})

 if(NOT WIN32)
    target_link_libraries(networkdesign_lib ${Boost_PROGRAM_OPTIONS_LIBRARY})
 endif()