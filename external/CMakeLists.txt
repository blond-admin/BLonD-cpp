cmake_minimum_required(VERSION 2.8)
project(externals)
include(ExternalProject)
add_custom_target(external) 

#Log required variables
message(STATUS "Building external packages")
message( STATUS "into ${INSTALL_DIR}" )
message( STATUS "WITH_GSL ${WITH_GSL}, GSL_PATCH_PATH ${GSL_PATCH_PATH} " )
message( STATUS "WITH_GOOGLETEST ${WITH_GOOGLETEST}" )
message( STATUS "WITH_BENCHMARK ${WITH_BENCHMARK}" )
message( STATUS "WITH_FFTW ${WITH_FFTW}" )
message( STATUS "Generator ${CMAKE_GENERATOR}" )



#Check if required libs and includes are already presented
set(INSTALL_LIB_DIR "${INSTALL_DIR}/lib")
set(INSTALL_BIN_DIR "${INSTALL_DIR}/bin")
set(INSTALL_INC_DIR "${INSTALL_DIR}/include")
message(STATUS "starting FFTW")
include("fftw.cmake")
include("gtest.cmake")
include("benchmark.cmake")
message(STATUS "ended FFTW")
