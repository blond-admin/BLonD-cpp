C++ Version of [CERN's BLonD code][1]

## Status 

NOT Stable - Under heavy development

[![Build Status](https://travis-ci.org/kiliakis/BLonD-cpp.svg?branch=master)](https://travis-ci.org/kiliakis/BLonD-cpp)
[![Coverage Status](https://coveralls.io/repos/github/kiliakis/BLonD-cpp/badge.svg?branch=master)](https://coveralls.io/github/kiliakis/BLonD-cpp?branch=master)

## Requirements

#### Linux
* cmake version >= 2.8 [install](https://cmake.org/install/)
* gcc version >= 4.8.0 [install](https://gcc.gnu.org/wiki/InstallingGCC)
* System packages:
  * Debian based distributions: 
    * python-dev
    * libfreetype6
    * libpng-dev
    * python-tk
  * Redhat based distributions:
    * python-devel
    * freetype-devel
    * libpng-devel
    * tcl
    * tkinter

#### Windows
* Windows version >= 7
* cygwin x86_64 [install] (https://cygwin.com/install.html)

## Installation

#### Linux
* clone the code into a directory (lets call it BLonD++/):
```bash  
git clone --branch=master https://github.com/kiliakis/BLonD-cpp BLonD++
```
* To compile all dependencies and build blond library run the commands:
```bash
cd BLonD++  
source install-linux.sh  
mkdir build  
cd build 
cmake .. # Configuration
cmake --build . # Compilation
ctest -VV # Testing
```

What was happening here:
   1. We opened folder with downloaded BLonD source files
   2. We installed external dependencies and exported some needed environmental variables
   3. Created a folder to hold solution and project files
   3. On configuration step we generated all executables (demos + unit-tests):
   4. Compiled and linked default build configuration
   5. Executed unit tests

* The executables should be ready!
* Developer's Notes:
  * On Linux, by default, the Release version of the code is compiled. You can build a debug version by adding `-DCMAKE_BUILD_TYPE=Debug` argument to configuration command, before `..`
  * To commit properly formatted code, reformatted by clang-format on each build please add `-DWITH_FORMAT=True` argument to configuration command, before `..`, note clang-format shall be [installed](http://llvm.org/releases/download.html)
  * To set data refrence files path's one can set `-DDATAFILES_DIR_DEMOS` and `-DDATAFILES_DIR_TESTS`.


#### Windows
* Download setup-x86_64.exe from cygwin installation page.
* Run setup-x86_64.exe, follow the on-screen instructions.
* In the "Select Packages" screen, the packages needed for a basic cygwin installation are already checked. Except from these, add:
  * git (devel category, bin only)
  * wget (web category, bin only)
* launch a cygwin terminal
* clone the code into a directory (lets call it BLonD++/):
```bash  
git clone --branch=master https://github.com/kiliakis/BLonD-cpp BLonD++
```
* To compile all dependencies and build blond library run the commands:
```bash
cd BLonD++  
source install-cygwin.sh  
mkdir build  
cd build 
cmake .. # Configuration
cmake --build . # Compilation
ctest -VV # Testing
```

## Using system Libraries (advanced)
**TODO**

## Configuration
The following definitions, found in file include/blond/configuration.h, can be commented / uncommented to alter simulation's configuration:
```c
#define TIMING
#define PRINT_RESULTS
```
*Note that a re-compile is needed every time a change is made.* 

## Usage
The following optional command line arguments can be specified in order to specify some basic simulation parameters:
* -m <num>, --threads=\<num\>   : Number of OpenMP threads that will be used in the simulation (default: 1)
* -p <num>, --particles=\<num\> : Number of macro particles that will be simulated (default: 10k)
* -t <num>, --turns=\<num\>     : Number of simulation turns (default: 10k)
* -s <num>, --slices=\<num\>    : Number of slices that will be used to generate the beam's histogram (default: 100)  

Example: `./testcase -t 1000 -p2000`  
Or type: `./testcase -h` for more  

## Running the Unit-Tests
Once you have successfully compiled the code you can run the tests:
```bash
cd BLonD++/build
ctest -VV
```
Then you can generate unit-test documentation:

## Building Documentation
To generate html documentation with search and graphical class hierarchy's please [install Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) and [Graphviz](http://www.graphviz.org/Download.php) and run:
```bash
cd BLonD++/build
doxygen Doxyfile
```

#### Building Unit-Test Documentation (Linux only)
To generate html documentation on unit tests coverage please install [genhtml](http://linux.die.net/man/1/genhtml) and [lcov](http://ltp.sourceforge.net/coverage/lcov.php)
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
- Oleg Jakushkin (oleg.jakushkin (at) gmail.com)

## Contributors Notice

Dear all contributors, you are kindly requested to format your code using astyle format options found [here] (https://root.cern.ch/coding-conventions#Astyle).

[1]: http://blond.web.cern.ch
