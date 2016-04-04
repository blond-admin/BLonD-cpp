cmake_minimum_required(VERSION 2.6)

project(BLonD-minimal-cpp)

include_directories(includes)


set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fopenmp -O2 -Ofast") #-Ofast
set(CMAKE_CXX_COMPILER g++)

set(SOURCE_FILES 
    beams/Beams.cpp 
    beams/Slices.cpp
    input_parameters/GeneralParameters.cpp
    input_parameters/RfParameters.cpp
    trackers/Tracker.cpp
    )

add_library(blond ${SOURCE_FILES})

    
add_executable(TC1_Acceleration tests/TC1_Acceleration.cpp)

target_link_libraries(TC1_Acceleration blond)