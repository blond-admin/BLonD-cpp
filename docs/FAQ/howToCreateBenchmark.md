#How to create a benchmark
##General algorithm
 1. Take a unit test
 2. Keep unit-test structure
 3. Keep unit-test data access
 4. Select a part that you want to speed up
 5. Isolate data inputs and test results e.g:
    * In: RfP->omega_RF, RfP->voltage; 
    * Out: Beam->dE;
 6. Make a benchmark that will use test case in one of the runs to always keep the code work-proofed
 7. Benchmark old code
 8. Create a prototype that can be a class or function accepting isolated input data
 9. Explore benchmark on different devices

##Code
###Project files
Add a `benchUnitTestNmae.cmake` file to find packages with find Package + set define configuration code alike
```cmake
find_package(MPI REQUIRED)
if(MPI)
	add_definitions("-DWITH_MPI=1")
endif()
```
###Source code
Add a `benchUnitTestNmae.cpp` and follow [Google Benchmark](https://github.com/google/benchmark) instructions to 
create a benchmark. Use `#ifdef` directives to filter package dependent code for example:
 ```cpp
 //...
 #ifdef WITH_MPI
 //MPI DEPENDENT TEST CODE
 #endif // WITH_MPI
 ```
 
 For more complex scenarios it is useful to create files that use prototype approach wrapping 
 input and generating output based on same interface showcasing different frameworks. 
##Run
Execute benchmarks as described [here](./howToRunBenchmarks.md) 