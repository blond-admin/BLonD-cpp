# usage: cmake -DWITH_FFTW=True -DWITH_GOOGLETEST=True -DWITH_BENCHMARK=True
# cmake -DWITH_FFTW=True -DWITH_GOOGLETEST=True -DWITH_BENCHMARK=True -DCMAKE_BUILD_TYPE=Debug to use asan
# for clion like editors that build in remote folders:
# -DDATAFILES_DIR_DEMOS="path/to/demos/input_files" -DDATAFILES_DIR_TESTS="unit-tests/references"; Default "../demos/input_files" and "../unit-tests/references" 

# Requirements list
cmake_minimum_required(VERSION 2.8)
project(BLonD-minimal-cpp CXX)
include(CheckCXXCompilerFlag)
include(ExternalProject)
include(ProcessorCount)
include(CTest)

# Code editors usabilety feature
set_property(GLOBAL PROPERTY USE_FOLDERS On)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "We have optemized configuration for Debug and Release configurations")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "We have optemized configuration for Debug and Release configurations")
# Variables default values list
set(CMAKE_BUILD_TYPE Release CACHE STRING "default build type for Make-like generators")
set(SHARED_LIBRARY_ENABLED "True")
set(EXTERNAL_INSTALL_DIR ${PROJECT_SOURCE_DIR}/external/install)
set(TEST_ENABLED "True" CACHE STRING "Cerate unit-test projects")
set(BENCHMARK_ENABLED "True" CACHE STRING "Cerate benchmark projects")
set(ARCHITECTURE "x64" CACHE STRING "Set architecture 'x86' or 'x32'") # todo: find usages and  fix comment
if(NOT ABSOLUTE_DATA_PATHS)
   set(DATAFILES_DIR_DEMOS "../demos/input_files" CACHE STRING "Set demos input files path")
   set(DATAFILES_DIR_TESTS "../unit-tests/references" CACHE STRING "Set unit-tests input files path")
else ()
   set(DATAFILES_DIR_DEMOS "${CMAKE_CURRENT_LIST_DIR}/demos/input_files" CACHE STRING "Set demos input files path")
   set(DATAFILES_DIR_TESTS "${CMAKE_CURRENT_LIST_DIR}/unit-tests/references" CACHE STRING "Set unit-tests input files path")
endif()

# Packages
set(EXTERNAL_INSTALL_DIR ${PROJECT_SOURCE_DIR}/external/install CACHE STRING "External libraries installation folder (with subfolders lib and include)")

# Abilety to build external libraries from source or install them from releases
if (NOT DEFINED WITH_FFTW)
    if (MSVC)
        set(WITH_FFTW "True" CACHE STRING "Should Build FFTW library from source")
    else ()
        set(WITH_FFTW "False" CACHE STRING "Should Build FFTW library from source")
    endif ()
endif ()

if (NOT DEFINED WITH_GOOGLETEST)
    if (MSVC)
        set(WITH_GOOGLETEST "True" CACHE STRING "Should Build GoogleTest library from source")
    else ()
        set(WITH_GOOGLETEST "False" CACHE STRING "Should Build GoogleTest library from source")
    endif ()
endif ()

set(GOOGLETEST_ROOT ${PROJECT_SOURCE_DIR}/external/googletest/ CACHE STRING "Google Test source root")

if (NOT DEFINED WITH_BENCHMARK)
    if (MSVC)
        set(WITH_BENCHMARK "True" CACHE STRING "Should Build GoogleBenchmark library from source")
    else ()
        set(WITH_BENCHMARK "False" CACHE STRING "Should Build GoogleBenchmark library from source")
    endif ()
endif ()

set(BENCHMARK_ROOT ${PROJECT_SOURCE_DIR}/external/benchmark CACHE STRING "Google Benchmark Library root")

set(WITH_FORMAT "False" CACHE STRING "Clang-Format all project files on build (requires clang-format)")

# Package settings
set(USE_FFTW_OMP "False" CACHE STRING "Should use OpenMP flavour of FFTW library (linux only)") #todo: test and configure
add_definitions(-DDEMO_FILES="${DATAFILES_DIR_DEMOS}") #demos datafiles

# CPU for build count
ProcessorCount(N)
if (NOT N EQUAL 0)
    set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif ()

# Project architecture
set(ARCHITECTURE "Default" CACHE STRING "Two architectures avaliable 'Dafault' and 'x64', it is impartant for external libraries")
if ((${CMAKE_GENERATOR} MATCHES ".*(64).*") OR (${CMAKE_SIZEOF_VOID_P} MATCHES ".*(64).*"))
    set(ARCHITECTURE "x64")
    message(STATUS "architecture: ${ARCHITECTURE} detected from (${CMAKE_GENERATOR} and ${CMAKE_SIZEOF_VOID_P}) ")
endif ()

#Compiler, Linker options
#Compiler settings
#MSVC vs Other
if (MSVC)
    add_definitions(-D_USE_MATH_DEFINES) #needed for VDT
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj -Zm500  /MP${N} /openmp /fp:precise ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  /GR-")
else ()
    check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)

    if (NOT COMPILER_HAS_CXX11)
        message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
    endif ()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11 -fopenmp -Wall -pedantic -fno-rtti -Wunused -ffast-math")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3") #-fsanitize=address -fno-omit-frame-pointer ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

    if (TEST_ENABLED STREQUAL "True")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} --coverage")
    else ()
        message(STATUS "code test coveredge is enabeled only in Debug mode")
    endif ()

    if (USE_FFTW_OMP)
        set(LIBRARIES ${LIBRARIES} fftw3_omp)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW_OMP")
    endif ()
endif ()

#GLANG and GCC specific
if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    # TODO is 4.8 really the minimum requirement?
    # TODO do we really want to check for this its 2016 already?
    if ((CMAKE_C_COMPILER_ID == "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0"))
        message(FATAL_ERROR "ERROR: Insufficient compiler version, need at lest 4.8 to continue")
    endif ()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")

elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    else ()
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
    endif ()

    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE}-flto")
endif ()

# Default packages
# required packages
find_package(OpenMP REQUIRED)

if (NOT WITH_FFTW)
    find_path(FFTW_INCLUDES fftw3.h)
    find_library(FFTW_LIB NAMES fftw fftw3 libfftw3 libfftw-3.3 REQUIRED)
endif ()
# message(STATUS "---------------- WITH_FFTW ${WITH_FFTW} ${FFTW_LIB}")

#not required
if ((NOT WITH_GOOGLETEST) AND TEST_ENABLED)
    find_package(GTEST REQUIRED)
    find_library(GTEST_LIB
            NAMES gtest libgtest
            )
    find_library(GTEST_MAIN_LIB
            NAMES gtest_main libgtest_main
            )

endif ()

if ((NOT WITH_BENCHMARK) AND BENCHMARK_ENABLED)
    find_path(BENCHMARK_INCLUDES benchmark/benchmark.h)
    find_library(BENCHMARK_LIB NAME benchmark REQUIRED)
endif ()

if (NOT (EXISTS ${PROJECT_SOURCE_DIR}/external/build))
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_SOURCE_DIR}/external/build)
else ()
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_SOURCE_DIR}/external/build/)
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_SOURCE_DIR}/external/build)
endif ()

message(STATUS "Building required started...")

execute_process(
        COMMAND ${CMAKE_COMMAND}
        #options
        -DWITH_FFTW=${WITH_FFTW}
        -DWITH_BENCHMARK=${WITH_BENCHMARK}
        -DWITH_GOOGLETEST=${WITH_GOOGLETEST}
        #destinations
        -DINSTALL_DIR=${EXTERNAL_INSTALL_DIR}
        -DCOMPILE_FLAGS=${COMPILE_FLAGS}
        -DGOOGLETEST_ROOT=${GOOGLETEST_ROOT}
        -DBENCHMARK_ROOT=${BENCHMARK_ROOT}
        -DGOOGLETEST_ROOT=${GOOGLETEST_ROOT}
        #forward all collected compilation arguments
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} # relevant only for Make-like build files
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
        -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
        -DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}
        -DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}
        -v
        -G "${CMAKE_GENERATOR}"
        ${PROJECT_SOURCE_DIR}/external/
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/external/build
)

message(STATUS "building external projects begins")
execute_process(
        COMMAND cmake --build . --config Release
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/external/build
)
if (WIN32 AND NOT MINGW)
    execute_process(
            COMMAND cmake --build . --config Debug
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/external/build
    )
endif ()
message(STATUS "building external projects ended")

message(STATUS "Building required libraries complete")
# path to external libraries and headers
include_directories(${EXTERNAL_INSTALL_DIR}/include)
link_directories(${EXTERNAL_INSTALL_DIR}/lib)

if (WITH_FFTW)
    if (WIN32 AND NOT MINGW)
        link_directories(${EXTERNAL_INSTALL_DIR}/lib/fftw/$(Configuration)/)
        find_library(FFTW_LIB
                NAMES fftw fftw3 libfftw3 libfftw-3.3
                PATHS "${EXTERNAL_INSTALL_DIR}/lib/fftw/Debug"
                NO_DEFAULT_PATH
                )
        get_filename_component(FFTW_LIB ${FFTW_LIB} NAME)
    else ()
        find_library(FFTW_LIB
                NAMES fftw fftw3 libfftw3 libfftw-3.3
                PATHS "${EXTERNAL_INSTALL_DIR}/lib/"
                NO_DEFAULT_PATH
                )
    endif ()

    add_custom_command(TARGET ALL PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${EXTERNAL_INSTALL_DIR}/bin/" $<TARGET_FILE_DIR>)
endif ()


# Code editors usabilety function (structure code as it is placed in folders)
function(assign_source_group)
    foreach (_source IN ITEMS ${ARGN})
        if (IS_ABSOLUTE ${_source})
            file(RELATIVE_PATH _source_rel ${CMAKE_CURRENT_SOURCE_DIR} ${_source})
        else ()
            set(source_rel ${_source})
        endif ()
        get_filename_component(_source_path ${_source_rel} PATH)
        string(REPLACE "/" "\\" _source_path_msvc ${_source_path})
        source_group("${_source_path_msvc}" FILES ${_source})
    endforeach ()
endfunction(assign_source_group)

include_directories(include
   $ENV{HOME}/opt/include
   $ENV{HOME}/Downloads/Python-2.7.12/Include
   $ENV{HOME}/Downloads/Python-2.7.12
   )

link_directories(/usr/lib
   $ENV{HOME}/opt/lib/python2.7/config
   )

# main library
include_directories(include)

SET(LIBRARIES
        ${FFTW_LIB}
        python2.7
        util
        dl
        )

file(GLOB_RECURSE SOURCES
        "src/*"
        "include/*"
        )

set(SOURCE_FILES ${SOURCES})
assign_source_group(${SOURCE_FILES})

# foreach(source ${SOURCE_FILES})
#   set_source_files_properties(${source} PROPERTIES COMPILE_FLAGS "-fopt-info-vec=report.txt")
# endforeach()

if (SHARED_LIBRARY_ENABLED)
    add_library(blond SHARED ${SOURCE_FILES})
    target_compile_definitions(blond PRIVATE SHARED=1)
else ()
    add_library(blond STATIC ${SOURCE_FILES})
    target_compile_definitions(blond PRIVATE STATIC=1)
endif ()
target_link_libraries(blond ${LIBRARIES})

#SET_TARGET_PROPERTIES(blond PROPERTIES LINKER_LANGUAGE CXX)

# demo executables
LIST(APPEND LIBRARIES blond)

file(GLOB EXECUTABLE_SOURCES "demos/*.cpp")
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES})



# demo dependencies

foreach (src ${EXECUTABLE_FILES})
    file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
    string(REPLACE ".cpp" "" testname ${testfile})
    string(REPLACE "demos/" "" testname ${testname})
    add_executable(${testname} ${src})
    target_link_libraries(${testname} ${LIBRARIES})
endforeach ()


# testing
if (TEST_ENABLED)
    add_definitions(-DTEST_FILES="${DATAFILES_DIR_TESTS}") #demos datafiles

    if (WITH_GOOGLETEST)
        if (WIN32 AND NOT MINGW)
            link_directories(${EXTERNAL_INSTALL_DIR}/lib/gtest/$(Configuration)/)
            find_library(GTEST_LIB
                    NAMES gtest libgtest
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/gtest/Debug"
                    NO_DEFAULT_PATH
                    )

            find_library(GTEST_MAIN_LIB
                    NAMES gtest_main libgtest_main
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/gtest/Debug"
                    NO_DEFAULT_PATH
                    )
            get_filename_component(GTEST_MAIN_LIB ${GTEST_MAIN_LIB} NAME)
            get_filename_component(GTEST_LIB ${GTEST_LIB} NAME)
        else ()
            find_library(GTEST_LIB
                    NAMES gtest libgtest
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/"
                    NO_DEFAULT_PATH
                    )
            find_library(GTEST_MAIN_LIB
                    NAMES gtest_main libgtest_main
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/"
                    NO_DEFAULT_PATH
                    )
        endif ()
    endif ()


    message(STATUS "generating test executables")
    LIST(APPEND LIBRARIES ${GTEST_LIB})
    LIST(APPEND LIBRARIES ${GTEST_MAIN_LIB})

    if (UNIX)
        set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
        set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
    endif ()

    file(GLOB UNIT_TEST_SOURCES "./unit-tests/*.cpp")
    set(UNIT_TEST_FILES ${UNIT_TEST_SOURCES})

    foreach (src ${UNIT_TEST_FILES})
        file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
        string(REPLACE ".cpp" "" testname ${testfile})
        string(REPLACE "unit-tests/" "" testname ${testname})
        add_executable(${testname} ${src})
        target_link_libraries(${testname} ${LIBRARIES})
        add_test(${testname} ${testname})
    endforeach ()

    enable_testing()
endif ()


# benchmarks
if (BENCHMARK_ENABLED)
    if (WITH_BENCHMARK)
        if (WIN32 AND NOT MINGW)
            link_directories(${EXTERNAL_INSTALL_DIR}/lib/fftw/$(Configuration)/)
            find_library(BENCHMARK_LIB
                    NAMES benchmark libbenchmark
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/benchmark/Debug"
                    NO_DEFAULT_PATH
                    )
            get_filename_component(BENCHMARK_LIB ${BENCHMARK_LIB} NAME)
        else ()
            find_library(BENCHMARK_LIB
                    NAMES benchmark libbenchmark
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/"
                    NO_DEFAULT_PATH
                    )
        endif ()
    endif ()

    message(STATUS "generating benchmark executables")
    LIST(APPEND LIBRARIES ${BENCHMARK_LIB})

    if (MSVC)
        LIST(APPEND LIBRARIES Shlwapi)
    else ()
        LIST(APPEND LIBRARIES pthread)
    endif ()

    file(GLOB BENCHMARK_SOURCES "./benchmarks/*.cpp")
    set(BENCHMARK_FILES ${BENCHMARK_SOURCES})

    foreach (src ${BENCHMARK_FILES})
        file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
        string(REPLACE ".cpp" "" testname ${testfile})
        string(REPLACE "benchmarks/" "" testname ${testname})
        add_executable(${testname} ${src})
        target_link_libraries(${testname} ${LIBRARIES})
        add_test(${testname} ${testname})
    endforeach ()

    enable_testing()
endif ()

install(TARGETS blond
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        INCLUDES DESTINATION include
        CONFIGURATIONS Debug|Release
        )

include(clang-format.cmake)