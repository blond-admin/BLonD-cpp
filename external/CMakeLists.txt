cmake_minimum_required(VERSION 2.8)
project(externals)
include(ExternalProject)
add_custom_target(external) 

#Log required variables
mseeage(STATUS "Building external packages")
mseeage( STATUS "into ${INSTALL_DIR}" )
mseeage( STATUS "WITH_GSL ${WITH_GSL}, GSL_PATCH_PATH ${GSL_PATCH_PATH} " )
mseeage( STATUS "WITH_GOOGLETEST ${WITH_GOOGLETEST}" )
mseeage( STATUS "WITH_BENCHMARK ${WITH_BENCHMARK}" )
mseeage( STATUS "WITH_FFTW ${WITH_FFTW}" )


function(GetnerateLibName arg)
    set(arg "${CMAKE_STATIC_LIBRARY_PREFIX}${arg}${CMAKE_STATIC_LIBRARY_SUFFIX}" PARENT_SCOPE)
endfunction()


#Check if required libs and includes are already presented
set(INSTALL_LIB_DIR "${INSTALL_DIR}/lib/")

if(WITH_GSL)
    set(gslLib "gsl")
    set(gslcblasLib "gslcblas")
    GetnerateLibName(${gslLib})
    GetnerateLibName(${gslcblasLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${gslLib}) AND (NOT EXISTS ${INSTALL_LIB_DIR}${gslcblasLib}) )
        set(BUILD_GSL "True")
    endif()
endif()
if(WITH_GOOGLETEST)
    set(gtestLib "gtest")
    GetnerateLibName(${gtestLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${gtestLib}))
        set(BUILD_GOOGLETEST "True")
    endif()
endif()
if(WITH_BENCHMARK)
    set(benchmarkLib "benchmark")
    GetnerateLibName(${benchmarkLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${gtestLib}))
        set(BUILD_BENCHMARK "True")
    endif()
endif()
if(WITH_FFTW)
    set(fftwLib "fftw")
    GetnerateLibName(${fftwLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${fftwLib}))
        set(BUILD_FFTW "True")
    endif()
endif()


if(BUILD_GSL)
        message(STATUS "generating GNU Scientific Library libraries")
    set(gslLib "${CMAKE_STATIC_LIBRARY_PREFIX}gsl${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    set(gslcblasLib "${CMAKE_STATIC_LIBRARY_PREFIX}gslcblas${CMAKE_STATIC_LIBRARY_SUFFIX}" )
    ExternalProject_Add(
        gsl-src
        GIT_REPOSITORY git://git.savannah.gnu.org/gsl.git
        DOWNLOAD_DIR ${GSL_ROOT}
        SOURCE_DIR ${GSL_ROOT}
        PATCH_COMMAND cmake
            -DPROJECT_SOURCE_DIR=${ORIGIN}
            -DGSL_ROOT=${GSL_ROOT}
            -P ${ORIGIN}/external/gsl.patch.cmake
        CMAKE_ARGS  -DCMAKE_BUILD_TYPE=Release
        -DDO_TESTS="False"
        -DCMAKE_INSTALL_PREFIX=${ORIGIN}/external/install
        -DCMAKE_RULE_MESSAGES=OFF
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
        -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
        -DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}
        TMP_DIR ${ORIGIN}/${GSL_ROOT}/build/tmp
        TEST_BEFORE_INSTALL False
        TEST_AFTER_INSTALL False
        BINARY_DIR ${ORIGIN}/${GSL_ROOT}/build
    )
    add_dependencies(external gsl-src)
endif()

if(BUILD_GOOGLETEST)
    message(STATUS "generating test libraries")
    # packages required for testing
    #GoogleTest (build from source into external instal)
    ExternalProject_Add(
        googletest-src
            GIT_REPOSITORY https://github.com/google/googletest
            DOWNLOAD_DIR ${GOOGLETEST_ROOT}
            SOURCE_DIR ${GOOGLETEST_ROOT}/
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                   -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
                   -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
                   -Dgtest_force_shared_crt=ON
                   -DBUILD_GTEST=ON
                   -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                   -DCMAKE_RULE_MESSAGES=OFF
        TMP_DIR ${GOOGLETEST_ROOT}/build/tmp
        TEST_BEFORE_INSTALL False
        TEST_AFTER_INSTALL False
        BINARY_DIR ${GOOGLETEST_ROOT}/build
    )
	add_dependencies(external googletest-src)
endif()

if(BUILD_BENCHMARK)
    message(STATUS "generating benchmark libraries")
    # packages required for benchmarking
    #Google Benchmark (build from source into external instal)
    message(STATUS "generating test executables")
    ExternalProject_Add(
        benchmark-src
            GIT_REPOSITORY https://github.com/google/benchmark
            DOWNLOAD_DIR ${BENCHMARK_ROOT}
            SOURCE_DIR ${BENCHMARK_ROOT}/
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
		 -DBENCHMARK_ENABLE_LTO=true
		-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
        TMP_DIR ${BENCHMARK_ROOT}/build/tmp
        TEST_BEFORE_INSTALL False
        TEST_AFTER_INSTALL False
        BINARY_DIR ${BENCHMARK_ROOT}/build
    )
	add_dependencies(external benchmark-src)
endif()

if(BUILD_FFTW)

endif()
