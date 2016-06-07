# Requirements list
cmake_minimum_required(VERSION 3.0.2)
project(BLonD-minimal-cpp CXX)
include(CheckCXXCompilerFlag)

# TODO is 4.8 really the minimum requirement?
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
endif()

check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)
if(!COMPILER_HAS_CXX11)
	message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
endif()



set(GOOGLETEST_ROOT external/googletest/googletest CACHE STRING "Google Test source root")
if(EXISTS ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT})
    set(TEST_ENABLED "True")
    message(STATUS "Found Google Test library source")
else()
    message(STATUS "Testing disabled - please provide GOOGLETEST_ROOT path or use '--recursive' argument in git clone")
endif()

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

#set(CMAKE_BUILD_TYPE Release)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
      set(CMAKE_CXX_FLAGS         "-std=c++11 -fopenmp -Wall -pedantic -fno-rtti -Wunused")
      set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
      set(CMAKE_CXX_FLAGS_RELEASE "-march=native")
      
      #CHECK_CXX_COMPILER_FLAG(-ffast-math COMPILER_HAS_FAST_MATH)
      #if(COMPILER_HAS_FAST_MATH)
      #   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
      #endif()

      #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

      CHECK_CXX_COMPILER_FLAG(-Ofast COMPILER_HAS_OFAST)
      if(COMPILER_HAS_OFAST)
         set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
      else()
         set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      endif()
      
      set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
      if(TEST_ENABLED)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
	     set(CMAKE_EXE_LINKER_FLAGS_DEBUG "--coverage")
      endif()
endif()

# Absolutely required packages
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)


if(TEST_ENABLED)
    include_directories(
        ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}
        ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/include
        )

    set(GOOGLETEST_SOURCES
        ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest-all.cc
        ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest_main.cc
        )

    foreach(_source ${GOOGLETEST_SOURCES})
        set_source_files_properties(${_source} PROPERTIES GENERATED 1)
    endforeach()

    add_library(googletest ${GOOGLETEST_SOURCES})
endif()

include_directories(include)
#link_directories(/usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 /afs/cern.ch/user/k/kiliakis/opt/lib)
link_directories(/usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)

# Generate simulation executables
set(SOURCE_FILES 
    beams/Beams.cpp 
    beams/Slices.cpp
    input_parameters/GeneralParameters.cpp
    input_parameters/RfParameters.cpp
    trackers/Tracker.cpp
    llrf/PhaseLoop.cpp
    impedances/Intensity.cpp
    impedances/InducedVoltage.cpp
    )

# -fopt-info-vec=file
foreach(source ${SOURCE_FILES})
   set_source_files_properties(${source} PROPERTIES COMPILE_FLAGS "-fopt-info-vec=report.txt")
endforeach()

# SHARED for dynamic LINKING
# STATIC for static LINKING
add_library(blond SHARED ${SOURCE_FILES})
    
add_executable(TC1_Acceleration tests/TC1_Acceleration.cpp)
add_executable(LHC_restart tests/LHC_restart.cpp)
add_executable(TC5_Wake_impedance tests/TC5_Wake_impedance.cpp)

SET(LIBRARIES 
   blond
   gsl
   gslcblas
   m)

target_link_libraries(TC1_Acceleration ${LIBRARIES})
target_link_libraries(LHC_restart ${LIBRARIES})
target_link_libraries(TC5_Wake_impedance ${LIBRARIES})

# Generate unit tests
if(TEST_ENABLED)
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
        unit-tests/testResonator.cpp
        unit-tests/testInputTableIntensity.cpp
        unit-tests/testInducedVoltage.cpp
        unit-tests/testTC5_Wake_impedance.cpp
        unit-tests/testInducedVoltageFreq.cpp
    )

    foreach(src ${TEST_SOURCE_FILES})
        string( REPLACE ".cpp" "" testname ${src})
        string( REPLACE "unit-tests/" "" testname ${testname})
        add_executable( ${testname} ${src})
	     target_link_libraries( ${testname} ${LIBRARIES} googletest)
        add_test(${testname} ${testname})
    endforeach( src ${TEST_SOURCE_FILES})

    #TODO: rename testMyMath, merge with TEST_SOURCE_FILES
    add_executable(testMyMath unit-tests/testMyMath.cpp)
    target_link_libraries(testMyMath googletest gsl gslcblas m)
    add_test(testMyMath testMyMath )

    include(CTest)
    enable_testing()

    if(UNIX)
      set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
      set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
    endif()
endif()








