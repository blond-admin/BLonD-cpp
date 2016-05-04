
# Requirements list
cmake_minimum_required(VERSION 3.0.2)
project(BLonD-minimal-cpp CXX)
include(CheckCXXCompilerFlag)

# Absolutely required packages
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

# Optional packages
find_package(GTest)

# TODO is 4.8 really the minimum requirement?
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
endif()

check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)
if(!COMPILER_HAS_CXX11)
	message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
endif()

check_cxx_compiler_flag(-ffast-math COMPILER_HAS_FAST_MATH)


include_directories(include)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fopenmp -O3") #-Ofast

if(COMPILER_HAS_FAST_MATH)
	add_compile_options(-ffast-math)
endif()


# Generate simulation executables

set(SOURCE_FILES 
    beams/Beams.cpp 
    beams/Slices.cpp
    input_parameters/GeneralParameters.cpp
    input_parameters/RfParameters.cpp
    trackers/Tracker.cpp
    llrf/PhaseLoop.cpp
    )

# SHARED for dynamic LINKING
# STATIC for static LINKING
add_library(blond SHARED ${SOURCE_FILES})

    
add_executable(TC1_Acceleration tests/TC1_Acceleration.cpp)
add_executable(LHC_restart tests/LHC_restart.cpp)


target_link_libraries(TC1_Acceleration blond)
target_link_libraries(LHC_restart blond)


# Generate unit tests

#if (GTest_FOUND)

	set(TEST_SOURCE_FILES
	    unit-tests/testGP.cpp
	    unit-tests/testRFP.cpp
	    unit-tests/testBigaussian.cpp
	    unit-tests/testBeam.cpp
	    unit-tests/testSlices.cpp
	    unit-tests/testTracker.cpp
	    unit-tests/testTC1_Acceleration.cpp
	    unit-tests/testPL.cpp
	    unit-tests/testLHC_Restart.cpp
	    )


	add_executable(testMyMath unit-tests/testMyMath.cpp)
	foreach(src ${TEST_SOURCE_FILES})
	    string( REPLACE ".cpp" "" testname ${src})
	    string( REPLACE "unit-tests/" "" testname ${testname})
	    add_executable( ${testname} ${src})
	    target_link_libraries( ${testname} blond gtest gtest_main)
	    add_test(${testname} ${testname})
	endforeach( src ${TEST_SOURCE_FILES})

	target_link_libraries(testMyMath gtest gtest_main)


	add_test(testMyMath testMyMath)

	enable_testing()

	if(UNIX)
	  set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
	  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
	endif()

#endif()
