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
            TMP_DIR ${BENCHMARK_ROOT}/build/tmp
            TEST_BEFORE_INSTALL False
            TEST_AFTER_INSTALL False
            BINARY_DIR ${BENCHMARK_ROOT}/build
            INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/install/
    )

    if(WIN32 AND NOT MINGW)
        #on windows projects require debug+release libraries
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/benchmark.lib ${INSTALL_LIB_DIR}/benchmark/Release/benchmark.lib)
        set(BENCHMARK_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/gbench-d/")

        ExternalProject_Add(
                benchmark-src-d
                GIT_REPOSITORY https://github.com/google/benchmark
                DOWNLOAD_DIR ${BENCHMARK_ROOT}
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=Debug
                -DBENCHMARK_ENABLE_LTO=true
                -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                TMP_DIR ${BENCHMARK_ROOT}/build/tmp
                TEST_BEFORE_INSTALL False
                TEST_AFTER_INSTALL False
                BINARY_DIR ${BENCHMARK_ROOT}/build
                INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/install/
        )

        execute_process(
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/benchmark.lib ${INSTALL_LIB_DIR}/benchmark/Debug/benchmark.lib
                COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_SOURCE_DIR}/install/lib/benchmark.lib
                )

    endif()
endif()