# Requirements list
cmake_minimum_required(VERSION 3.0.2)
project(BLonD-minimal-cpp CXX)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)
include(CTest)


# Variables default values list
set(CMAKE_BUILD_TYPE Release)
set(TEST_ENABLED "False")
set(GOOGLETEST_ROOT external/googletest/googletest CACHE STRING "Google Test source root")
set(GSL_ROOT external/gsl CACHE STRING "GSL - GNU Scientific Library source root")
set(EXTERNAL_ROOT external/install CACHE STRING "External libraries installation folder (with subfolders lib and include)")


# CPU for build count
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CTEST_BUILD_FLAGS -j${N})
  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()


# Code editors usabilety feature
set_property(GLOBAL PROPERTY USE_FOLDERS On) 

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "-std=c++11 -fopenmp -Wall -pedantic -Wno-rtti -Wunused")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-march=native")
    
    if ((CMAKE_C_COMPILER_ID == "Clang") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8))
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    else()
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
    endif()

    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
    if((TEST_ENABLED STREQUAL "True") AND ( CMAKE_BUILD_TYPE STREQUAL "Debug"))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS "--coverage")
    endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj -Zm500  /MP${N}")
endif()

# TODO is 4.8 really the minimum requirement?
# TODO do we really want to check for this its 2016 already?
if(CMAKE_CXX_COMPILER_ID MATCHES GNU AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
endif()
check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)
if(!COMPILER_HAS_CXX11)
	message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
endif()


# optimization checks
check_cxx_compiler_flag(-ffast-math COMPILER_HAS_FAST_MATH)
if(COMPILER_HAS_FAST_MATH)
	add_compile_options(-ffast-math)
endif()


# Autoenabeling tests if gtest is presented
if(EXISTS ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT})
    set(TEST_ENABLED "True")
    message(STATUS "Found Google Test library source")
else()
    message(STATUS "Testing disabled - please provide GOOGLETEST_ROOT path or use '--recursive' argument in git clone")
endif()


# required packages
find_package(OpenMP REQUIRED)
#GSL (build from source into external instal)
if(EXISTS ${PROJECT_SOURCE_DIR}/${GSL_ROOT})
    message(STATUS "Found GNU Scientific Library ")
    ExternalProject_Add(
        gsl-src
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/${GSL_ROOT}
	PATCH_COMMAND git apply --3way -v --allow-overlap ${PROJECT_SOURCE_DIR}/external/gsl.diff && cmake -E copy ${PROJECT_SOURCE_DIR}/${GSL_ROOT}/contrib/CMakeLists.txt ${PROJECT_SOURCE_DIR}/${GSL_ROOT}/CMakeLists.txt
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DDO_TESTS="False"
                   -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/install
        TMP_DIR ${PROJECT_SOURCE_DIR}/${GSL_ROOT}/build/tmp
        TEST_BEFORE_INSTALL False
        TEST_AFTER_INSTALL False
        BINARY_DIR ${PROJECT_SOURCE_DIR}/${GSL_ROOT}/build
    )
else()
  message(FATAL_ERROR "Not GNU Scientific Library - please provide GSL_ROOT path or use '--recursive' argument in git clone")
endif()


# path to external libraries and headers
include_directories(${EXTERNAL_ROOT}/include)
link_directories(${EXTERNAL_ROOT}/lib)


# Code editors usabilety function (structure code as it is placed in folders)
function(assign_source_group)
    foreach(_source IN ITEMS ${ARGN})
	message(${_source})
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
   gslcblas)

file(GLOB SOURCES
    "beams/*"
    "input_parameters/*"
    "trackers/*"
    "llrf/*"
    "impedances/*"
    "include/*"
)
set(SOURCE_FILES ${SOURCES} )
assign_source_group(${SOURCE_FILES})

# SHARED for dynamic LINKING
# STATIC for static LINKING
add_library(blond SHARED ${SOURCE_FILES})
target_link_libraries( blond ${LIBRARIES})


# demo executables
LIST(APPEND LIBRARIES blond)

file(GLOB EXECUTABLE_SOURCES "tests/*.cpp" )
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES} )
# demo dependencies 

foreach(src ${EXECUTABLE_FILES})
    message(STATUS ${src})
    file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
    string( REPLACE ".cpp" "" testname ${testfile})
    string( REPLACE "tests/" "" testname ${testname})
    add_executable( ${testname} ${src})
    add_dependencies(${testname} blond gsl-src)
    target_link_libraries( ${testname} ${LIBRARIES})
endforeach()


# testing
if(TEST_ENABLED)
    message(STATUS "generating test executables")
    # packages required for testing
    #GoogleTest (build from source into external instal)
    ExternalProject_Add(
        googletest-src
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} 
                   -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
                   -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
                   -Dgtest_force_shared_crt=ON
                   -DBUILD_GTEST=ON
                   -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/install
        TMP_DIR ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/build/tmp
        TEST_BEFORE_INSTALL False
        TEST_AFTER_INSTALL False
        BINARY_DIR ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/build
    )
    LIST(APPEND LIBRARIES gtest)


    if(UNIX)
      set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
      set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
    endif()

    file(GLOB UNIT_TEST_SOURCES "./unit-tests/*.cpp" )
    set(UNIT_TEST_FILES ${UNIT_TEST_SOURCES} )

    foreach(src ${UNIT_TEST_FILES})
        message(STATUS ${src})
        file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
        string( REPLACE ".cpp" "" testname ${testfile})
        string( REPLACE "unit-tests/" "" testname ${testname})
        add_executable( ${testname} ${src})
	add_dependencies(${testname} googletest-src gsl-src)
        target_link_libraries( ${testname} ${LIBRARIES} gtest)
        add_test(${testname} ${testname})
    endforeach()
    
    enable_testing()
endif()
