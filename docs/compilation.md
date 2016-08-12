#Advanced build recipes
Having source code downloaded into `BLonD-minimal-cpp`

##Windows
Currently we build and test only Visual Studio solutions on windows. Yet one may desire to use only opensource
toolchain such as MinGW or another development environment such as Clion.

###MSYS2 (MinGW, Bash shell, x64)
Download and install [MSYS2](https://sourceforge.net/projects/msys2/)
Open MSYS2 console and install packages
```bash
update-core # requires session restart
# start session with msys2_shell.cmd -mingw64
pacman -Suu # updates old packages
#pacman -Ss # to list avaliable packages
#install MinGW with C, C++, OpenMP, OpenCL, CMake packages:
pacman -S base-devel cmake mingw-w64-x86_64-toolchain mingw-w64-x86_64-opencl-headers mingw-w64-x86_64-cmake
```
change directory into BLonD-minimal-cpp.
Make sure `BLonD-minimal-cpp/external/install` and `BLonD-minimal-cpp/external/build` folders are empty.
Create a build directory and compile project as usual. Note: build will be slower than when MSVC is used. Also note that
Windows Defender and other anty-virus software slows compilation down.
```
mkdir build
cd build
cmake -DWITH_FFTW=True -DWITH_GOOGLETEST=True -DWITH_BENCHMARK=True -DBUILD_EXTERNALS=True  -G"MSYS Makefiles"  ..
```

### Clion
To use IDE you need to
1. compile project using MSYS2 toolchain
2. set CLion "`Preferences/Settings -> Build, Execution, Deployment -> Toolchains`" to used in MSYS2 MinGW folder for gcc
 and set cmake from mingw bin folder.
3. set CLion "`Preferences/Settings -> Build, Execution, Deployment -> Cmake -> options`" to
`-DWITH_FFTW=True -DWITH_GOOGLETEST=True -DWITH_BENCHMARK=True`
Notes:
1. Default build configuration is "Debug".
2. Set "`Preferences/Settings -> Build, Execution, Deployment -> Cmake -> options`" -j1 to ensure cmake log and build
order in cases of repeated failure

### Mixing CLion and MSVC
1. Build project under CLion
2. Delete `BLonD-minimal-cpp/external/install/include` and `BLonD-minimal-cpp/external/build` folders
3. CLion is not using same directory for build as we used under MSYS so remove all from  `BLonD-minimal-cpp/build`
folder
4. `Shift+RMB -> Open Command Prompt here` to open command prompt in build folder
5. Run cmake as usual `cmake -DWITH_FFTW=True -DWITH_GOOGLETEST=True -DWITH_BENCHMARK=True -DBUILD_EXTERNALS=True` it
will take time
6. Open generated .sln file in MSVC
