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
            GIT_REPOSITORY https://github.com/google/googletest
            DOWNLOAD_DIR ${GOOGLETEST_ROOT}
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
            -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
            -Dgtest_force_shared_crt=ON
            -DBUILD_GTEST=ON
            -DBUILD_GMOCK=OFF
            -Dgtest_disable_pthreads=ON
            -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
            -DCMAKE_RULE_MESSAGES=OFF
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}
            -DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}
            TMP_DIR ${GOOGLETEST_ROOT}/build/tmp
            TEST_BEFORE_INSTALL False
            TEST_AFTER_INSTALL False
            BINARY_DIR ${GOOGLETEST_ROOT}/build
            INSTALL_DIR ${INSTALL_DIR}
    )

    if(WIN32 AND NOT MINGW)
        #on windows projects require debug+release libraries
        #adding additional step to move release build libraries
        ExternalProject_Add_Step(googletest-src AFTER_INSTALL
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/gtest.lib ${INSTALL_LIB_DIR}/gtest/Release/gtest.lib
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/gtest_main.lib ${INSTALL_LIB_DIR}/gtest/Release/gtest_main.lib
                COMMENT "Copy release files into release directory"
                DEPENDEES Install
                )

        #building debug version
        set(GOOGLETEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/gtest-d/")
        ExternalProject_Add(
                googletest-src-dbg
                GIT_REPOSITORY https://github.com/google/googletest
                DOWNLOAD_DIR ${GOOGLETEST_ROOT}
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=Debug
                -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
                -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
                -Dgtest_force_shared_crt=ON
                -DBUILD_GTEST=ON
                -DBUILD_GMOCK=OFF
                -Dgtest_disable_pthreads=ON
                -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                -DCMAKE_RULE_MESSAGES=OFF
                -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
                -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
                -DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}
                -DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}
                TMP_DIR ${GOOGLETEST_ROOT}/build/tmp
                TEST_BEFORE_INSTALL False
                TEST_AFTER_INSTALL False
                BINARY_DIR ${GOOGLETEST_ROOT}/build
                INSTALL_DIR ${INSTALL_DIR}
        )

        ExternalProject_Add_Step(googletest-src-dbg AFTER_INSTALL
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/gtest.lib ${INSTALL_LIB_DIR}/gtest/Debug/gtest.lib
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/gtest_main.lib ${INSTALL_LIB_DIR}/gtest/Debug/gtest_main.lib
                COMMAND ${CMAKE_COMMAND} -E remove -f ${INSTALL_LIB_DIR}/gtest.lib
                COMMAND ${CMAKE_COMMAND} -E remove -f ${INSTALL_LIB_DIR}/gtest_main.lib
                COMMENT "Copy release files into release directory"
                DEPENDEES Install
                )

        message(STATUS "note: on windows lib builds are seprate for debug and release")
    endif()
endif()