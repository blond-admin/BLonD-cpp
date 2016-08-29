#Cmake compilation options
You can force the compiler (on the first run in clean build folder) using `-DCMAKE_C_COMPILER` and `-DCMAKE_CXX_COMPILER`

## Main build parameters
 * Build test projects `TEST_ENABLED`
 * Build benchmark projects `BENCHMARK_ENABLED`
 * Build external libraries `BUILD_EXTERNALS`
   * Build Google Benchmark `WITH_BENCHMARK`
   * Build Google Test `WITH_GOOGLETEST`	
   * Build boost `WITH_COMPUTE`	
   * build FFTW `WITH_FFTW`	
 * Format on solution build ``WITH_FORMAT`	
 * `OPENCL_SDK_PATH` Path of OpenCL SDK (or CUDA)	
 * To set data reference files path's one can set `DATAFILES_DIR_DEMOS` and `DATAFILES_DIR_TESTS`.
 * Set bool `ABSOLUTE_DATA_PATHS` to embed absolute test data path's instead of default relative to build folder


