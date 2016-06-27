if(WITH_GSL)
    set(gslLib "gsl")
    set(gslcblasLib "gslcblas")
    GetnerateLibName(${gslLib})
    GetnerateLibName(${gslcblasLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${gslLib}) AND (NOT EXISTS ${INSTALL_LIB_DIR}${gslcblasLib}) )
        set(BUILD_GSL "True")
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