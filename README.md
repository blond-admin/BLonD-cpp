C++ Version of [CERN's BLonD code][1]

## Status

NOT Stable - Under heavy development

[![Build Status](https://travis-ci.org/kiliakis/BLonD-minimal-cpp.svg?branch=master)](https://travis-ci.org/kiliakis/BLonD-minimal-cpp)
[![Coverage Status](https://coveralls.io/repos/github/kiliakis/BLonD-minimal-cpp/badge.svg?branch=master)](https://coveralls.io/github/kiliakis/BLonD-minimal-cpp?branch=master)

## Requirements
* cmake version >= 3.0.2 [install](https://cmake.org/install/)
* gcc version >= 4.8.0 [install](https://gcc.gnu.org/wiki/InstallingGCC)
* GNU Scientific Library [install](http://www.gnu.org/software/gsl/)


## Installation
1. clone the code into a directory (lets call it home/)  
    ```bash  
    git clone --branch=master --recursive https://github.com/kiliakis/BLonD-minimal-cpp.git home  
    ```

2. run the commands 
    ```bash
    cd home
    mkdir build  
    cd build  
    cmake ..  
    make  
    ```

3. The executables should be ready!

## Configuration

The following definitions, found in file include/configuration.h, can be commented / uncommented to alter simulation's configuration:

```c
#define FIXED_PARTICLES
#define TIMING
#define PRINT_RESULTS
```

Note that a re-compile is needed every time a change is made. 

## Usage

The following optional command line arguments can be specified in order to specify some basic simulation parameters:

* -m <num>, --threads=\<num\> : Number of openmp threads that will be used in the simulation (default: 1)
* -p <num>, --particles=\<num\> : Number of macroparticles that will be simulated (default: 10k)
* -t <num>, --turns=\<num\>     : Number of simulation turns (default: 10k)
* -s <num>, --slices=\<num\>    : Number of slices that will be used to generate the beam's histogram (default: 100)

Example: `./testcase -t 1000 -p2000`  
Or type: `./testcase -h` for more

## Buiding and runniung Unit-Tests (googletest)
```bash
cd home  
mkdir build 
cd build   
cmake ..  
make
ctest -VV
```

## Buiding Unit-Test Documentation
```bash
cd home  
mkdir build 
cd build   
cmake ..  
make
ctest -VV
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

