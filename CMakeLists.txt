cmake_minimum_required(VERSION 2.6)

project(BLonD-minimal-cpp)

include_directories(includes)

#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBRARIES OFF)

#release
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fopenmp -O2 -Ofast") #-Ofast
#debug
#set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fopenmp -O0 -g") #-Ofast

set(CMAKE_CXX_COMPILER g++)

set(SOURCE_FILES 
    beams/Beams.cpp 
    beams/Slices.cpp
    input_parameters/GeneralParameters.cpp
    input_parameters/RfParameters.cpp
    trackers/Tracker.cpp
    llrf/PhaseLoop.cpp
    )

add_library(blond STATIC ${SOURCE_FILES})

    
add_executable(TC1_Acceleration tests/TC1_Acceleration.cpp)
add_executable(LHC_restart tests/LHC_restart.cpp)

set(TEST_SOURCE_FILES
    unit-tests/testMyMath.cpp
    unit-tests/testGP.cpp
    unit-tests/testRFP.cpp
    unit-tests/testBigaussian.cpp
    unit-tests/testBeam.cpp
    unit-tests/testSlices.cpp
    unit-tests/testTracker.cpp
    unit-tests/testTC1_Acceleration.cpp
    )

foreach(src ${TEST_SOURCE_FILES})
    string( REPLACE ".cpp" "" testname ${src})
    string( REPLACE "unit-tests/" "" testname ${testname})
    add_executable( ${testname} ${src})
    target_link_libraries( ${testname} blond gtest gtest_main)
    add_test(${testname} ${testname})
endforeach( src ${TEST_SOURCE_FILES})


target_link_libraries(TC1_Acceleration blond)
target_link_libraries(LHC_restart blond)


enable_testing()
include(CTest)
if(UNIX)
  set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
endif()

#add_test(test1 test1)