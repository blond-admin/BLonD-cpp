C++ Version of [CERN's BLonD code][1]

## Status 

NOT Stable - Under heavy development

[![Build Status](https://travis-ci.org/kiliakis/BLonD-minimal-cpp.svg?branch=master)](https://travis-ci.org/kiliakis/BLonD-minimal-cpp)
[![Coverage Status](https://coveralls.io/repos/github/kiliakis/BLonD-minimal-cpp/badge.svg?branch=master)](https://coveralls.io/github/kiliakis/BLonD-minimal-cpp?branch=master)

## Requirements
* cmake version >= 3.0.2 [install](https://cmake.org/install/)
* gcc version >= 4.8.0 [install](https://gcc.gnu.org/wiki/InstallingGCC)
* GNU Scientific Library [install](http://www.gnu.org/software/gsl/)
* FFTW3 Library [install](http://www.fftw.org/download.html)


## Installation

#### Libraries Installation

* GNU Scientific Library
If you are using a Linux distribution, in most of the cases, the following set of commands should be enough:
    ```bash
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
    tar -xzvf gsl-2.1.tar.gz gsl-home
    cd gsl-2.1
    ./configure
    make
    make install #you may need sudo access for this step, if you don't have it then try ./configure --prefix=/path/to/install
    ```

* FFTW3 Library
If you are using a Linux distribution, in most of the cases, the following set of commands should be enough:
    ```bash
    wget http://www.fftw.org/fftw-3.3.4.tar.gz
    tar -xzvf fftw-3.3.4.tar.gz
    cd fftw-3.3.4
    ./configure  
    make 
    make install #you may need sudo access for this step, if you don't have it then try ./configure --prefix=/path/to/install
    ```
  *Note that if you want to use the multi-threaded version of this library you must configure as `./configure --enable-openmp`.*

#### BLonD++ installation

1. clone the code into a directory (lets call it BLonD++/)  
    ```bash  
    git clone --branch=master --recursive https://github.com/kiliakis/BLonD-minimal-cpp.git BLonD++    
    ```

2. run the commands 
    ```bash
    cd BLonD++
    mkdir build  
    cd build  
    cmake .. # or cmake -DUSE_FFTW_OMP=True .. for the multithreaded version
    make  
    ```

3. The executables should be ready!

4. Developer's Notes:
  * By default, the Release version of the code is compiled. You can build a debug version with `cmake -DCMAKE_BUILD_TYPE=Debug ..`     

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

* -m <num>, --threads=\<num\> : Number of OpenMP threads that will be used in the simulation (default: 1)
* -p <num>, --particles=\<num\> : Number of macro particles that will be simulated (default: 10k)
* -t <num>, --turns=\<num\>     : Number of simulation turns (default: 10k)
* -s <num>, --slices=\<num\>    : Number of slices that will be used to generate the beam's histogram (default: 100)

Example: `./testcase -t 1000 -p2000`  
Or type: `./testcase -h` for more

## Running the Unit-Tests (googletest)
Once you have successfully compiled the code you can run the tests:
```bash
cd BLonD++/build
ctest -VV
```
Then you can generate unit-test documentation:
## Building Unit-Test Documentation
```bash
cd BLonD++/build
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

