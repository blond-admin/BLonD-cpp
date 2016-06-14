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
      set(CMAKE_CXX_FLAGS_DEBUG   "-O3 -g3")# -fsanitize=address")
      set(CMAKE_CXX_FLAGS_RELEASE "-march=native")
      
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
#find_package(GSL REQUIRED)
find_path(FFTW_INCLUDES fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3 libfftw3)


#set(LIBRARIES blond gsl gslcblas fftw3 m)
set(LIBRARIES blond fftw3 m)

if(USE_FFTW_OMP)
   set(LIBRARIES ${LIBRARIES} fftw3_omp)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW_OMP")
endif()

#if(USE_FFTW)
#message(STATUS "Using FFTW")
#set(LIBRARIES "${LIBRARIES} -lfftw3 -lfftw3_omp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW")
#endif()

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
link_directories(/usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
link_directories("$ENV{HOME}/opt/lib")

file(GLOB SOURCES 
      beams/*.cpp
      input_parameters/*.cpp
      trackers/*.cpp
      llrf/*.cpp
      impedances/*.cpp
      )

set(SOURCE_FILES ${SOURCES})

foreach(source ${SOURCE_FILES})
   set_source_files_properties(${source} PROPERTIES COMPILE_FLAGS "-fopt-info-vec=report.txt")
endforeach()

# SHARED for dynamic LINKING
# STATIC for static LINKING
add_library(blond SHARED ${SOURCE_FILES})




file(GLOB EXECUTABLE_SOURCES "tests/*.cpp")
set(EXECUTABLE_FILES ${EXECUTABLE_SOURCES})

foreach(src ${EXECUTABLE_FILES})
   file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
   string( REPLACE ".cpp" "" testname ${testfile})
   string( REPLACE "tests/" "" testname ${testname})
   add_executable( ${testname} ${src})
   target_link_libraries( ${testname} ${LIBRARIES})
endforeach()
    
# Generate unit tests
if(TEST_ENABLED)
   file(GLOB TEST_SOURCES unit-tests/*.cpp)
   set(TEST_SOURCE_FILES ${TEST_SOURCES})

   foreach(src ${TEST_SOURCE_FILES})
      file(RELATIVE_PATH testfile ${CMAKE_CURRENT_SOURCE_DIR} ${src})
      string( REPLACE ".cpp" "" testname ${testfile})
      string( REPLACE "unit-tests/" "" testname ${testname})
      add_executable( ${testname} ${src})
      target_link_libraries( ${testname} ${LIBRARIES} googletest)
      add_test(${testname} ${testname})
   endforeach( src ${TEST_SOURCE_FILES})

   include(CTest)
   enable_testing()

   if(UNIX)
      set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
      set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full")
   endif()
endif()

