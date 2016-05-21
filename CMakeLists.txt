# Requirements list
cmake_minimum_required(VERSION 3.0.2)
project(BLonD-minimal-cpp CXX)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(CTest)

# TODO is 4.8 really the minimum requirement?
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
    message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
endif()

check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)
if(!COMPILER_HAS_CXX11)
	message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
endif()

check_cxx_compiler_flag(-ffast-math COMPILER_HAS_FAST_MATH)
if(COMPILER_HAS_FAST_MATH)
	add_compile_options(-ffast-math)
endif()

set(TEST_ENABLED "False")
set(GOOGLETEST_ROOT external/googletest/googletest CACHE STRING "Google Test source root")
if(EXISTS ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT})
    set(TEST_ENABLED "True")
    message(STATUS "Found Google Test library source")
else()
    message(STATUS "Testing disabled - please provide GOOGLETEST_ROOT path or use '--recursive' argument in git clone")
endif()

set(GSL_ROOT external/gsl CACHE STRING "GSL - GNU Scientific Library source root")

    ExternalProject_Add(
        gsl
        GIT_REPOSITORY git://git.savannah.gnu.org/gsl.git
        DOWNLOAD_DIR external/gsl
        SOURCE_DIR external/gsl
        PATCH_COMMAND file(COPY external/gsl/contrib/CMakeLists.txt external/gsl/contrib/ )
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        BINARY_DIR ./
        INSTALL_DIR ./
    )

if(EXISTS ${PROJECT_SOURCE_DIR}/${GSL_ROOT})
    message(STATUS "Found GNU Scientific Library ")
else()
  # todo resolve by downloading external project 
  #  message(FATAL_ERROR "Not GNU Scientific Library - please provide GSL_ROOT path or use '--recursive' argument in git clone")
endif()


set(CMAKE_BUILD_TYPE Release)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "-std=c++11 -fopenmp -Wall -pedantic -Wno-rtti -Wunused")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")
    if(TEST_ENABLED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
	    set(CMAKE_EXE_LINKER_FLAGS "--coverage")
    endif()
endif()

# Absolutely required packages
find_package(OpenMP REQUIRED)

if(TEST_ENABLED)
    ExternalProject_Add(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        DOWNLOAD_DIR external/googletest
        SOURCE_DIR external/googletest
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        BINARY_DIR ./
        INSTALL_DIR ./
    )

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
include_directories(${source_dir}/include)

endif()

include_directories(include)

file(GLOB SOURCES
    "beams/*.cpp"
    "input_parameters/*.cpp"
    "trackers/*.cpp"
    "llrf/*.cpp"
)

set(SOURCE_FILES 
    ${SOURCES}
)

# SHARED for dynamic LINKING
# STATIC for static LINKING
add_library(blond SHARED ${SOURCE_FILES})

file(GLOB EXECUTABLE_SOURCES
    "tests/*.cpp"
)

set(EXECUTABLE_FILES
    ${EXECUTABLE_SOURCES}
)

foreach(src ${EXECUTABLE_FILES})
    file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
    string( REPLACE ".cpp" "" testname ${testfile})
    string( REPLACE "tests/x "" testname ${testname})
    add_executable( ${testname} ${src})
    target_link_libraries( ${testname} blond)
endforeach()

# Generate unit tests
if(TEST_ENABLED)
    file(GLOB UNIT_TEST_SOURCES
        "./unit-tests/*.cpp"
    )

    set(UNIT_TEST_FILES
        $(UNIT_TEST_SOURCES)
    )

    foreach(src ${UNIT_TEST_FILES})
        file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
        string( REPLACE ".cpp" "" testname ${testfile})
        string( REPLACE "unit-tests/" "" testname ${testname})
        add_executable( ${testname} ${src})
        target_link_libraries( ${testname} blond googletest)
        add_test(${testname} ${testname})
    endforeach()
    
    enable_testing()

    if(UNIX)
      set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
      set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
    endif()
endif()








