C++ Version of [CERN's BLonD code][1]

## Status 

NOT Stable - Under heavy development

[![Build Status](https://travis-ci.org/kiliakis/BLonD-minimal-cpp.svg?branch=master)](https://travis-ci.org/kiliakis/BLonD-minimal-cpp)
[![Build status](https://ci.appveyor.com/api/projects/status/onoxb8504u8kwfkk?svg=true)](https://ci.appveyor.com/project/kiliakis/blond-minimal-cpp)
[![Coverage Status](https://coveralls.io/repos/github/kiliakis/BLonD-minimal-cpp/badge.svg?branch=master)](https://coveralls.io/github/kiliakis/BLonD-minimal-cpp?branch=master)

## Requirements
* cmake version >= 2.8 [install](https://cmake.org/install/)

####Linux
* gcc version >= 4.8.0 [install](https://gcc.gnu.org/wiki/InstallingGCC)

####Windows
* Visual Studio version >= 2015 with C++ compiler [install](https://www.visualstudio.com/)

## Installation (Windows and Linux)
1. clone the code into a directory (lets call it BLonD-minimal-cpp/):
```bash  
git clone --branch=master https://github.com/kiliakis/BLonD-minimal-cpp
```
2. To compile all dependencies and build blond library run the commands:
```bash
cd BLonD-minimal-cpp
mkdir build  
cd build 
cmake -DWITH_FFTW=True -DWITH_GOOGLETEST=True -DWITH_BENCHMARK=True .. # Configuration
cmake --build . # Compilation
ctest -VV # Testing
```
What was happening here:
   1. we opened folder with downloaded Blond source files
   2. created a folder to hold solution and project files
   3. On configuration step:
     1. Downloaded build and installed external libraries (FFTW, GoogleTest, GoogleBenchmark) into `BLonD-minimal-cpp\external\install`
     2. Generated solution and project files
   4. Compiled and linked default build configuration
   5. Executed unit tests
3. The executables should be ready!
4. Developer's Notes:
  * On Linux, by default, the Release version of the code is compiled. You can build a debug version by adding `-DCMAKE_BUILD_TYPE=Debug` argument to configuration command, before `..`
  * On Windows by default, the Debug version of the code is compiled. You can build a debug version by adding `--target ALL_BUILD --config Release` argument to Compilation command after `.`
  * On Windows one shall copy contents of `external/install/bin/fftw/$(configuration)` into corresponding to given build configuration folder to be able to test and run executables.
  * To commit properly formatted code, reformatted by clang-format on each build please add `-DWITH_FORMAT=True` argument to configuration command, before `..`, note clang-format shall be [installed](http://llvm.org/releases/download.html)


## Using system Libraries (advanced)
If FFTW, GoogleTest or GoogleBenchmark are already installed in your system you can set `-DWITH_*` to `False` or skip this arguments when calling configuration commands

## Configuration
The following definitions, found in file include/blond/configuration.h, can be commented / uncommented to alter simulation's configuration:
```c
#define TIMING
#define PRINT_RESULTS
``
*Note that a re-compile is needed every time a change is made.* 

## Usage
The following optional command line arguments can be specified in order to specify some basic simulation parameters:

* -m <num>, --threads=\<num\> : Number of OpenMP threads that will be used in the simulation (default: 1)
* -p <num>, --particles=\<num\> : Number of macro particles that will be simulated (default: 10k)
* -t <num>, --turns=\<num\>     : Number of simulation turns (default: 10k)
* -s <num>, --slices=\<num\>    : Number of slices that will be used to generate the beam's histogram (default: 100)

Example: `./testcase -t 1000 -p2000`  
Or type: `./testcase -h` for more

## Running the Unit-Tests
Once you have successfully compiled the code you can run the tests:
```bash
cd BLonD-minimal-cpp/build
ctest -VV
```
Then you can generate unit-test documentation:

## Building Documentation
To generate html documentation with search and graphical class hierarchy's please [install Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) and [Graphviz](http://www.graphviz.org/Download..php) and run:
```bash
cd BLonD-minimal-cpp/build
doxygen Doxyfile
```

#### Building Unit-Test Documentation (Linux only)
To generate html documentation on unit tests coverage please install [genhtml](http://linux.die.net/man/1/genhtml) and [lcov](http://ltp.sourceforge.net/coverage/lcov.php)
```bash
cd BLonD-minimal-cpp/build
lcov --capture --directory .. --output-file coverage.info
genhtml coverage.info --output-directory html
```

## Original BLonD Links

* Repository: https://gitlab.cern.ch/dquartul/BLonD
* Documentation: http://blond-documentation.web.cern.ch/
* Project website: http://blond.web.cern.ch

## Developers

- Alexandre Lasheen (alexandre.lasheen (at) cern.ch)
- Juan Esteban Muller (juan.fem (at) cern.ch)
- Danilo Quartullo (danilo.quartullo (at) cern.ch)
- Helga Timko (Helga.Timko (at) cern.ch)
- Konstantinos Iliakis (konstantinos.iliakis (at) cern.ch)

## Contributors Notice

Dear all contributors, you are kindly requested to format your code using astyle format options found [here] (https://root.cern.ch/coding-conventions#Astyle).

[1]: http://blond.web.cern.ch
