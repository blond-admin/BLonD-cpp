#How to add a library (advanced)
It is hard to add an external library. One shall take into consideration all supported platforms and CI scripts.

Know project structure: there is main library, unit-tests, benchmarks. Follow this guide if you want to add library for 
the main library. In case you would prefer to add a benchmark proving it is effective please read [how to add 
benchmark](./howToCreateBenchmark.md). 

If your library is optional and can be used by main BLonD library do not forget to add defines to the 
main configuration file `BLonD-minimal-cpp/CMakeLists.txt` in case when the library is found using 
`add_definitions("-DWITH_YOUR_LIB=1")`.

##Look around
Add to main project configuration `BLonD-minimal-cpp/CMakeLists.txt` a snippet that would search for `YOUR_LIB`  
```cmake
if (NOT WITH_YOUR_LIB)
    find_package(YOUR_LIB_PATH REQUIRED)
    find_library(YOUR_LIB
            NAMES YOUR_LIB_name1
            )
endif ()
```

##Add library compilation script
Pass `-DWITH_YOUR_LIB=${WITH_YOUR_LIB}` to `execute_process->COMMAND` in the `BUILD_EXTERNALS` section of 
the main configuration file `BLonD-minimal-cpp/CMakeLists.txt`. We build external libraries only if user wants us to.
Also we try to build libraries only once it is needed.
###Add boost library
The headers are already added if `BUILD_EXTERNALS` option is passed to the CMake configuration script.
To add boost library that requires compilation edit `BLonD-minimal-cpp/external/boost.cmake` file extending 
`ExternalProject_Add->BUILD_COMMAND` by adding to it `--with-desired_library`. 

###CMake based library
Create a file named `YOUR_LIB.cmake` in `BLonD-minimal-cpp/external/` folder. It shall be based on these templeate replacing `YOUR_LIB` and cmake configuration with relevant data:
```cmake
if(WITH_YOUR_LIB)
    message(STATUS "searching for ${INSTALL_INC_DIR}/YOUR_LIB/YOUR_LIB_important.h")
    if( (NOT EXISTS ${INSTALL_INC_DIR}/YOUR_LIB/YOUR_LIB_important.h))
        set(BUILD_YOUR_LIB "True")
    endif()
endif()


if(BUILD_YOUR_LIB)
    message(STATUS "not found, building YOUR_LIB")
    set(YOUR_LIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/YOUR_LIB/")

    ExternalProject_Add(
            YOUR_LIB-src
            GIT_REPOSITORY https://github.com/YOUR_LIB
            DOWNLOAD_DIR ${YOUR_LIB_ROOT}
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_EXE_LINKER_FLAGS_RELEASE=${CMAKE_EXE_LINKER_FLAGS_RELEASE}
            -DCMAKE_EXE_LINKER_FLAGS_DEBUG=${CMAKE_EXE_LINKER_FLAGS_DEBUG}
            TMP_DIR ${YOUR_LIB_ROOT}/build/tmp
            TEST_BEFORE_INSTALL False
            TEST_AFTER_INSTALL False
            BINARY_DIR ${YOUR_LIB_ROOT}/build
            INSTALL_DIR ${INSTALL_DIR}
    )

    if(WIN32 AND NOT MINGW)
        #on windows projects require debug+release libraries

        #adding additional step to move release build libraries
        ExternalProject_Add_Step(YOUR_LIB-src AFTER_INSTALL
                COMMAND ${CMAKE_COMMAND} -E copy ${INSTALL_LIB_DIR}/YOUR_LIB.lib ${INSTALL_LIB_DIR}/YOUR_LIB/$(Configuration)/YOUR_LIB.lib
                COMMAND ${CMAKE_COMMAND} -E remove -f ${INSTALL_LIB_DIR}/YOUR_LIB.lib
                COMMENT "Copy release files into release directory"
                DEPENDEES Install
                )

        message(STATUS "note: on windows lib builds are seprate for debug and release")
    endif()
endif()
```
As you can see we 
 1. checked if library is already installed
 2. Only if library is not istalled:
    1. Added external project
    2. specifide a way to get required source files (see [ExternalProject documentation](https://cmake.org/cmake/help/v3.0/module/ExternalProject.html) for details)
    3. passed to it same compiler arguments as are passed to main blond project
    4. Added logic for library artifacts moving on windows
    
    An example of such build 
    script can be found here `BLonD-minimal-cpp/external/benchmark.cmake`.
###Other library
To build  [ExternalProjects](https://cmake.org/cmake/help/v3.0/module/ExternalProject.html) that are specific for each 
platform we use expand `BUILD_COMMAND` using array of CMake `COMMAND` lines for each platform. An example of such build 
script can be found here `BLonD-minimal-cpp/external/fftw.cmake`.
###Search for externaly build library files
 1. Add `include("YOUR_LIB.cmake")` to the end of `BLonD-minimal-cpp/external/CMakeLists.txt` file to enable its build..
 2. After the `BUILD_EXTERNALS` section of the main configuration file `BLonD-minimal-cpp/CMakeLists.txt` add code that would
detect library files following this example:
```cmake
if (WITH_YOUR_LIB)
        if (WIN32 AND NOT MINGW)
            link_directories(${EXTERNAL_INSTALL_DIR}/lib/YOUR_LIB/$(Configuration)/)
            find_library(YOUR_LIB_LIB_PATH
                    NAMES YOUR_LIB_name
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/YOUR_LIB/Debug"
                    NO_DEFAULT_PATH
                    )
            get_filename_component(YOUR_LIB_PATH ${YOUR_LIB_PATH} NAME)
        else ()
            find_library(YOUR_LIB_PATH
                    NAMES YOUR_LIB_name
                    PATHS "${EXTERNAL_INSTALL_DIR}/lib/"
                    NO_DEFAULT_PATH
                    )
        endif ()
    endif ()
```

##Link project to library
In the main project configuration `BLonD-minimal-cpp/CMakeLists.txt` add library path to the list of `LIBRARIES` using 
```cmake
 LIST(APPEND LIBRARIES
         ${YOUR_LIB_PATH}
         )
```
##Test CI
Edit `BLonD-minimal-cpp/appveyor.yml` and `.travis.yml` scripts to build with your library.



