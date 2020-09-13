project(leaphi-test)

set(SOURCE_FILES networksolver_test.cpp dataprocessing_test.cpp networkdesign_test.cpp)

add_executable(leaphi-test ${SOURCE_FILES})

target_link_libraries(leaphi-test networksolver_lib networkdesign_lib kevintools_lib dataprocessing_lib)
target_link_libraries(leaphi-test gtest gtest_main)

#add_test(NAME all COMMAND state-test)
install(TARGETS leaphi-test DESTINATION bin)