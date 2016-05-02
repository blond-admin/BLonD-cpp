C++ Version of [CERN's BLonD code][1]

## STATUS

NOT Stable - Under heavy development

## REQUIREMENTS

* gcc version >= 4.8.0 [install](https://gcc.gnu.org/wiki/InstallingGCC)
* cmake version >= 3.0.2 [install](https://cmake.org/install/)
* GNU Scientific Library [install](http://www.gnu.org/software/gsl/)

#### Optional

* Google Test [install](https://github.com/google/googletest)

## INSTALATION

1. clone the code into a directory (lets call it home/)
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

The following environmental variables can be specified in order to dynamically alter some basic simulation parameters:

* N_THREADS : number of openmp threads that will be used for the simulation (default: 1)
* N_PARTICLES : number of macroparticles that will be simulated (default: 10k)
* N_TURNS : number of simulation turns (default: 10k)
* N_SLICES : number of slices that will be used to generate the beam's histogram (default: 100)

Example: `N_THREADS=4 N_TURNS=10000 ./testcase`

## LINKS

* Repository: https://gitlab.cern.ch/dquartul/BLonD
* Documentation: http://blond-documentation.web.cern.ch/
* Project website: http://blond.web.cern.ch

## DEVELOPERS

- Alexandre Lasheen (alexandre.lasheen (at) cern.ch)
- Juan Esteban Muller (juan.fem (at) cern.ch)
- Danilo Quartullo (danilo.quartullo (at) cern.ch)
- Helga Timko (Helga.Timko (at) cern.ch)
- Konstantinos Iliakis (konstantinos.iliakis (at) cern.ch)


[1]: http://blond.web.cern.ch
