# Requirements list
cmake_minimum_required(VERSION 2.8)
project(BLonD-minimal-cpp CXX)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)
include(CTest)

# Code editors usabilety feature
set_property(GLOBAL PROPERTY USE_FOLDERS On)

# Variables default values list
set(CMAKE_BUILD_TYPE Release)
set(SHARED_LIBRARY_ENABLED "True")
set(TEST_ENABLED "True" CACHE STRING "Cerate unit-test projects")
set(BENCHMARK_ENABLED "True" CACHE STRING "Cerate benchmark projects")
set(ARCHITECTURE "x64" CACHE STRING "Set architecture 'x86' or 'x32'")

# Packages
set(EXTERNAL_INSTALL_DIR external/install CACHE STRING "External libraries installation folder (with subfolders lib and include)")

# Abilety to build external libraries from source or install them from releases
set(WITH_GSL "True" CACHE STRING "Should Build GNU Scientific Library from source") #todo: default to false, always default true on windows
set(GSL_ROOT external/gsl CACHE STRING "GSL - GNU Scientific Library source root")

set(WITH_FFTW "True" CACHE STRING "Should Build FFTW library from source") #todo: default to false, always default true on windows
set(FFTW_ROOT external/gsl CACHE STRING "GSL - GNU Scientific Library source root")

set(WITH_GOOGLETEST "True" CACHE STRING "Should Build GoogleTest library from source")
set(GOOGLETEST_ROOT external/googletest/googletest CACHE STRING "Google Test source root")

set(WITH_BENCHMARK "True" CACHE STRING "Should Build GoogleBenchmark library from source")
set(BENCHMARK_ROOT external/benchmark  CACHE STRING "Google Benchmark Library root")

# Package settings
set(USE_FFTW_OMP "False" CACHE STRING "Should use OpenMP flavour of FFTW library")


# CPU for build count
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CTEST_BUILD_FLAGS -j${N})
  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()

#Compiler, Linker options

#Compiler settings
#MSVC vs Other
if(MSVC)
    if(USE_FFTW_OMP)
         message(WARNINIG "ERROR: No OMP in FFTW on Windows MSVC")
         set(USE_FFTW_OMP "False")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj -Zm500  /MP${N} /openmp /fp:fast")
else()
    check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)

	if(!COMPILER_HAS_CXX11)
	    message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
	endif()

    set(CMAKE_CXX_FLAGS "-fPIC -std=c++11 -fopenmp -Wall -pedantic -Wno-rtti -Wunused")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-march=native")

    if(TEST_ENABLED STREQUAL "True")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "--coverage")
    else()
        message(STATUS "code test coveredge is enabeled only in Debug mode")
    endif()

	check_cxx_compiler_flag(-ffast-math COMPILER_HAS_FAST_MATH)
	if(COMPILER_HAS_FAST_MATH)
	    add_compile_options(-ffast-math)
	endif()


    if(USE_FFTW_OMP)
       set(LIBRARIES ${LIBRARIES} fftw3_omp)
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW_OMP")
    endif()
endif()

#GLANG and GCC specific
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    # TODO is 4.8 really the minimum requirement?
    # TODO do we really want to check for this its 2016 already?
    if((CMAKE_C_COMPILER_ID == "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0"))
        message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
    endif()

    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    else()
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
    endif()
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -fPIC -std=c++11 -fopenmp -Wall -pedantic -Wno-rtti -Wunused")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-march=native")

    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
endif()

# Default packages
# required packages
find_package(OpenMP REQUIRED)

if(!WITH_GSL)
  find_package(GSL REQUIRED)
endif()

if(!WITH_FFTW)
  find_path(FFTW_INCLUDES fftw3.h)
  find_library(FFTW_LIBRARIES NAMES fftw3 libfftw3)
endif()

#not required
if(!WITH_GOOGLETEST AND TEST_ENABLED)
  find_package( GTEST )
endif()

if(!WITH_BENCHMARK AND BENCHMARK_ENABLED)
  find_package( GTEST )
endif()

    if(NOT (EXISTS ${PROJECT_SOURCE_DIR}/external/build))
        execute_process( COMMAND cmake -E make_directory ${PROJECT_SOURCE_DIR}/external/build )
    endif()

    message(STATUS "Building required started...")
    execute_process(
    COMMAND cmake 
        -DCOMPILE_FLAGS=${COMPILE_FLAGS}
        -DGOOGLETEST_ROOT=${GOOGLETEST_ROOT}
    -DBENCHMARK_ROOT=${BENCHMARK_ROOT}
        -DGOOGLETEST_ROOT=${GOOGLETEST_ROOT}
        -DGSL_ROOT=${GSL_ROOT}
        -DORIGIN=${PROJECT_SOURCE_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -G ${CMAKE_GENERATOR}
        ${PROJECT_SOURCE_DIR}/external/
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/external/build
    )

            execute_process(
                COMMAND cmake --build .
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/external/build
            )
            message(STATUS "Building required complete")
    endif()

    # path to external libraries and headers
    include_directories(${PROJECT_SOURCE_DIR}/${EXTERNAL_INSTALL_DIR}/include)
    link_directories(${PROJECT_SOURCE_DIR}/${EXTERNAL_INSTALL_DIR}/lib)
else()
    # Absolutely required packages
    find_package(GSL REQUIRED)
    if(TEST_ENABLED)
        find_package(GTest REQUIRED)
    endif()
endif()


# Code editors usabilety function (structure code as it is placed in folders)
function(assign_source_group)
    foreach(_source IN ITEMS ${ARGN})
        if(IS_ABSOLUTE ${_source})
            file(RELATIVE_PATH _source_rel ${CMAKE_CURRENT_SOURCE_DIR} ${_source})
        else()
            set(source_rel ${_source})
        endif()
        get_filename_component(_source_path ${_source_rel} PATH)
        string(REPLACE "/" "\\" _source_path_msvc ${_source_path})
        source_group("${_source_path_msvc}" FILES ${_source})
    endforeach()
endfunction(assign_source_group)


# main library
include_directories(include)
SET(LIBRARIES 
   gsl
   gslcblas
)

file(GLOB_RECURSE SOURCES
    "src/*"
    "include/*"
)
set(SOURCE_FILES ${SOURCES} src/globals.cpp src/globals.h)
assign_source_group(${SOURCE_FILES})

if(SHARED_LIBRARY_ENABLED)
    add_library(blond SHARED ${SOURCE_FILES})
    target_compile_definitions(blond PRIVATE SHARED=1)
else()
    add_library(blond STATIC ${SOURCE_FILES})
    target_compile_definitions(blond PRIVATE STATIC=1)
endif()
target_link_libraries( blond ${LIBRARIES})

#SET_TARGET_PROPERTIES(blond PROPERTIES LINKER_LANGUAGE CXX)

# demo executables
LIST(APPEND LIBRARIES blond)

file(GLOB EXECUTABLE_SOURCES "demos/*.cpp" )
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES} )
# demo dependencies 

foreach(src ${EXECUTABLE_FILES})
    file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
    string( REPLACE ".cpp" "" testname ${testfile})
    string( REPLACE "demos/" "" testname ${testname})
    add_executable( ${testname} ${src})
    target_link_libraries( ${testname} ${LIBRARIES})
endforeach()


# testing
if(TEST_ENABLED)
    message(STATUS "generating test executables")
    LIST(APPEND LIBRARIES gtest)

    if(UNIX)
      set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
      set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
    endif()

    file(GLOB UNIT_TEST_SOURCES "./unit-tests/*.cpp" )
    set(UNIT_TEST_FILES ${UNIT_TEST_SOURCES} )

    foreach(src ${UNIT_TEST_FILES})
        file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
        string( REPLACE ".cpp" "" testname ${testfile})
        string( REPLACE "unit-tests/" "" testname ${testname})
        add_executable( ${testname} ${src})
        target_link_libraries( ${testname} ${LIBRARIES})
        add_test(${testname} ${testname})
    endforeach()
    
    enable_testing()
endif()


# benchmarks
if(BENCHMARK_ENABLED)
    message(STATUS "generating benchmark executables")
    LIST(APPEND LIBRARIES benchmark)

    if(MSVC)
    LIST(APPEND LIBRARIES Shlwapi)
    else()
    LIST(APPEND LIBRARIES pthread)
    endif()

    file(GLOB BENCHMARK_SOURCES "./benchmarks/*.cpp" )
    set(BENCHMARK_FILES ${BENCHMARK_SOURCES} )

    foreach(src ${BENCHMARK_FILES})
        file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
        string( REPLACE ".cpp" "" testname ${testfile})
        string( REPLACE "benchmarks/" "" testname ${testname})
        add_executable( ${testname} ${src})
        target_link_libraries( ${testname} ${LIBRARIES})
        add_test(${testname} ${testname})
    endforeach()
    
    enable_testing()
endif()
