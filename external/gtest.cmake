if(WITH_GOOGLETEST)
    set(gtestLib "gtest")
    GetnerateLibName(${gtestLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${gtestLib}))
        set(BUILD_GOOGLETEST "True")
    endif()
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