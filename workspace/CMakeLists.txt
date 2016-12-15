# Requirements list
cmake_minimum_required(VERSION 2.8.4)
project(MyBLonD++ CXX)
include(CheckCXXCompilerFlag)

set(USAGE_MESSAGE "usage: cmake -DCMAKE_BUILD_TYPE=Debug/Release \
                                -DBLOND_DIR=path/to/blond_cpp \
                                -DMAIN_FILES_DIR=path/to/mymain/files \
                                ..")

# Variables default values list
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(NOT BLOND_DIR)
    message(FATAL_ERROR "ERROR: You must specify the path to BLonD directory\
                        \n${USAGE_MESSAGE}")
endif()

if(NOT MAIN_FILES_DIR)
    message(FATAL_ERROR "ERROR: You must specify the path to your main files\
                        \n${USAGE_MESSAGE}")
endif()


set(SHARED_LIBRARY_ENABLED "False")
set(EXTERNAL_INSTALL_DIR ${BLOND_DIR}/external/install)
set(USE_OMP "True" CACHE STRING
    "OpenMP is used by default to parallelize the code")
set(USE_FFTW_OMP "False" CACHE STRING 
    "Should use OpenMP flavour of FFTW library (linux only)")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BLOND_BUILD_DIR "build")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BLOND_BUILD_DIR "debug")
endif()

add_definitions(-DPROJECT_HOME="${PROJECT_SOURCE_DIR}"
                -DPYHOME="${EXTERNAL_INSTALL_DIR}"
                -DPYPATH="${BLOND_DIR}/python/")


#Compiler settings
check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CXX11)
check_cxx_compiler_flag(-Ofast HAS_OFAST)
check_cxx_compiler_flag(-ffast-math HAS_FFAST_MATH)

if (NOT COMPILER_HAS_CXX11)
    message(FATAL_ERROR "ERROR: Compiler is lacking support for c++11")
endif ()

set(WARNINGS "-Wall -pedantic")
set(CXXFLAGS "-std=gnu++11")
set(DEBUG "-O0 -g -fno-omit-frame-pointer")
set(RELEASE "-march=native")

if(CYGWIN)
    set(CXXFLAGS "${CXXFLAGS} -mcmodel=large")
else ()
    set(CXXFLAGS "${CXXFLAGS} -INTERFACE_POSITION_INDEPENDENT_CODE")
endif ()

if(HAS_OFAST)
   set(RELEASE "${RELEASE} -Ofast")
elseif(HAS_FFAST_MATH)
   set(RELEASE "${RELEASE} -O3 -ffast-math")
else()
   set(RELEASE "${RELEASE} -O3")
endif()

#gcc, clang and icc specific
if (CMAKE_CXX_COMPILER_ID MATCHES GNU)

    set(CXXFLAGS "${CXXFLAGS} -funsafe-math-optimizations")
    set(RELEASE "${RELEASE} -Ofast")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
    # set(DEBUG "${DEBUG} -D_GLIBCXX_DEBUG")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

    set(CXXFLAGS "${CXXFLAGS} -funsafe-math-optimizations")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")

elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    # set(VEC "-qopt-report=5 -qopt-report-file=vec_report.txt")
else ()
    message(FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} is not supported")
endif()

# Libraries list
# openmp
find_package(OpenMP QUIET)
if ((OPENMP_FOUND) AND (USE_OMP))
    message(STATUS "OpenMP found")
    set(OPENMP "-DUSE_OMP -fopenmp")
else ()
    message(STATUS "OpenMP will not be used")
    set(USE_FFTW_OMP "False")
endif ()

# fftw
find_path(FFTW_INCLUDES fftw3.h
    PATHS ${EXTERNAL_INSTALL_DIR}/include
    NO_DEFAULT_PATH)

if(FFTW_INCLUDES STREQUAL "FFTW_INCLUDES-NOTFOUND")
    message(FATAL_ERROR "ERROR: fftw3.h not found")
endif()

# omp fftw
if (USE_FFTW_OMP)
    find_library(FFTW_OMP_LIB 
        NAMES fftw_omp fftw3_omp fftw3_omp-3 fftw3-3_omp
        PATHS ${EXTERNAL_INSTALL_DIR}/lib
        NO_DEFAULT_PATH
        REQUIRED)
    set(OPENMP "${OPENMP} -DUSE_FFTW_OMP")
    set(FFTW_LIB ${FFTW_LIB} ${FFTW_OMP_LIB})
else ()
    find_library(FFTW_LIB NAMES fftw fftw3 fftw3-3 fftw-3.3
        PATHS ${EXTERNAL_INSTALL_DIR}/lib
        NO_DEFAULT_PATH
        REQUIRED)
endif ()

# blond library
find_path(BLOND_HEADER NAMES blond.h
    PATHS "${BLOND_DIR}/include/blond/"
    NO_DEFAULT_PATH)

if(BLOND_HEADER STREQUAL "BLOND_HEADER-NOTFOUND")
    message(FATAL_ERROR "ERROR: blond.h not found")
endif()


find_library(BLOND_LIB NAMES blond libblond
    PATHS "${BLOND_DIR}/${BLOND_BUILD_DIR}"
    NO_DEFAULT_PATH
    REQUIRED)


include_directories(
    ${PROJECT_SOURCE_DIR}/include/
    ${BLOND_DIR}/include/
    ${BLOND_DIR}/external/install/include/
    ${BLOND_DIR}/external/install/include/python2.7/
    ${BLOND_DIR}/external/install/lib/python2.7/site-packages
    )

link_directories(
    ${EXTERNAL_INSTALL_DIR}/lib
   )


set(CMAKE_CXX_FLAGS "${WARNINGS} ${CXXFLAGS} ${OPENMP}")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE} ${VEC}")
set(CMAKE_CXX_FLAGS_DEBUG "${DEBUG}")

list(APPEND PYTHON_LIB python2.7)
list(APPEND HDF5_LIB hdf5 hdf5_cpp hdf5_hl hdf5_hl_cpp)
list(APPEND REST_LIBS util dl)
list(APPEND LIBRARIES
        ${FFTW_LIB}
        ${PYTHON_LIB}
        ${HDF5_LIB}
        ${REST_LIBS}
        ${BLOND_LIB}
        )

# build the executables
file(GLOB_RECURSE EXECUTABLE_SOURCES "${MAIN_FILES_DIR}*.cpp")
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES})
foreach (src ${EXECUTABLE_FILES})
    file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
    string(REPLACE ".cpp" "" testname ${testfile})
    string(REPLACE "${MAIN_FILES_DIR}" "" testname ${testname})
    add_executable(${testname} ${src})
    target_link_libraries(${testname} ${LIBRARIES})
endforeach ()
