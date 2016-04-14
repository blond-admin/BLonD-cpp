cmake_minimum_required(VERSION 2.6)

project(BLonD-minimal-cpp)

include_directories(includes)

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

add_library(blond ${SOURCE_FILES})

    
add_executable(TC1_Acceleration tests/TC1_Acceleration.cpp)
add_executable(LHC_restart tests/LHC_restart.cpp)
#add_executable(TC8_Phase_Loop tests/TC8_Phase_Loop.cpp)

target_link_libraries(TC1_Acceleration blond)
target_link_libraries(LHC_restart blond)