if(WITH_BENCHMARK)
    message(STATUS "searching for ${INSTALL_INC_DIR}/benchmark/benchmark.h")
    if( (NOT EXISTS ${INSTALL_INC_DIR}/benchmark/benchmark.h))
        set(BUILD_BENCHMARK "True")
    endif()
endif()


if(BUILD_BENCHMARK)
    message(STATUS "not found, building benchmark")
    set(BENCHMARK_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/gbench/")

    ExternalProject_Add(
            benchmark-src
            GIT_REPOSITORY https://github.com/google/benchmark
            DOWNLOAD_DIR ${BENCHMARK_ROOT}
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
            -DBENCHMARK_ENABLE_LTO=true
            -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}
            -DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}
            TMP_DIR ${BENCHMARK_ROOT}/build/tmp
            TEST_BEFORE_INSTALL False
            TEST_AFTER_INSTALL False
            BINARY_DIR ${BENCHMARK_ROOT}/build
            INSTALL_DIR ${INSTALL_DIR}
    )

    if(WIN32 AND NOT MINGW)
        #on windows projects require debug+release libraries

        #adding additional step to move release build libraries
        ExternalProject_Add_Step(benchmark-src AFTER_INSTALL
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/benchmark.lib ${INSTALL_LIB_DIR}/benchmark/Release/benchmark.lib
                COMMENT "Copy release files into release directory"
                DEPENDEES Install
                )

        #building debug version
        set(BENCHMARK_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/gbench-d/")
        ExternalProject_Add(
                benchmark-src-d
                GIT_REPOSITORY https://github.com/google/benchmark
                DOWNLOAD_DIR ${BENCHMARK_ROOT}
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=Debug
                -DBENCHMARK_ENABLE_LTO=true
                -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
                -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
                -DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}
                -DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}
                TMP_DIR ${BENCHMARK_ROOT}/build/tmp
                TEST_BEFORE_INSTALL False
                TEST_AFTER_INSTALL False
                BINARY_DIR ${BENCHMARK_ROOT}/build
                INSTALL_DIR ${INSTALL_DIR}
        )

        ExternalProject_Add_Step(benchmark-src-d AFTER_INSTALL
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/benchmark.lib ${INSTALL_LIB_DIR}/benchmark/Debug/benchmark.lib
                COMMAND ${CMAKE_COMMAND} -E remove -f ${INSTALL_LIB_DIR}/benchmark.lib
                COMMENT "Copy release files into release directory"
                DEPENDEES Install
                )

        message(STATUS "note: on windows lib builds are seprate for debug and release")
    endif()
endif()