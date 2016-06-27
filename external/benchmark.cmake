if(WITH_BENCHMARK)
    set(benchmarkLib "benchmark")
    GetnerateLibName(${benchmarkLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${gtestLib}))
        set(BUILD_BENCHMARK "True")
    endif()
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