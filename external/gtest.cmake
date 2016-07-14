if(WITH_GOOGLETEST)
    message(STATUS "searching for  ${INSTALL_INC_DIR}/gtest/gtest.h")
    if( (NOT EXISTS ${INSTALL_INC_DIR}/gtest/gtest.h))
        set(BUILD_GOOGLETEST "True")
    else()
        set(BUILD_GOOGLETEST "False")
    endif()
endif()

if(BUILD_GOOGLETEST)
    message(STATUS "not found, building gtest")
    set(GOOGLETEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/gtest/")
    ExternalProject_Add(
            googletest-src
            GIT_REPOSITORY "https://github.com/google/googletest.git"
            DOWNLOAD_DIR ${GOOGLETEST_ROOT}
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
            -Dgtest_force_shared_crt=ON
            -DBUILD_GTEST=ON
            -DBUILD_GMOCK=OFF
            -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
            -DCMAKE_RULE_MESSAGES=OFF
            TMP_DIR ${GOOGLETEST_ROOT}/build/tmp
            TEST_BEFORE_INSTALL False
            TEST_AFTER_INSTALL False
            BINARY_DIR ${GOOGLETEST_ROOT}/build
            INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/install
    )

    if(WIN32 AND NOT MINGW)
        #on windows projects require debug+release libraries
        execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/gtest.lib ${INSTALL_LIB_DIR}/gtest/Release/gtest.lib)
        execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/gtest_main.lib ${INSTALL_LIB_DIR}/gtest/Release/gtest_main.lib)
        set(GOOGLETEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/gtest-d/")
        ExternalProject_Add(
                googletest-src-dbg
                GIT_REPOSITORY "https://github.com/google/googletest.git"
                DOWNLOAD_DIR ${GOOGLETEST_ROOT}
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=Debug
                -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
                -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
                -Dgtest_force_shared_crt=ON
                -DBUILD_GTEST=ON
                -DBUILD_GMOCK=OFF
                -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                -DCMAKE_RULE_MESSAGES=OFF
                TMP_DIR ${GOOGLETEST_ROOT}/build/tmp
                TEST_BEFORE_INSTALL False
                TEST_AFTER_INSTALL False
                BINARY_DIR ${GOOGLETEST_ROOT}/build
                INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/install/libs/Debug
        )

        execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/gtest.lib ${INSTALL_LIB_DIR}/gtest/Debug/gtest.lib)
        execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/gtest_main.lib ${INSTALL_LIB_DIR}/gtest/Debug/gtest_main.lib)
    endif()
endif()