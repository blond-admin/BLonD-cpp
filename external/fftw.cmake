#based on https://github.com/hinerm/ITK/blob/master/CMake/itkExternal_FFTW.cmake
# Encapsulates building FFTW as an External Project.
function(GetnerateLibName arg)
    set(arg "${CMAKE_STATIC_LIBRARY_PREFIX}${arg}${CMAKE_STATIC_LIBRARY_SUFFIX}" PARENT_SCOPE)
endfunction()

if (WITH_FFTW)
    set(fftwLib "fftw")
    GetnerateLibName(${fftwLib})
    message(STATUS "searching for ${INSTALL_INC_DIR}/fftw3.h")
    if ((NOT EXISTS {INSTALL_INC_DIR}/fftw3.h))
        set(BUILD_FFTW "True")
    endif ()
endif ()

if (BUILD_FFTW)
    message(STATUS "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    if (WIN32 AND NOT MINGW)
        message(STATUS "Can't build fftw as external project on Windows, downloading official release")
        message(STATUS "Can't build fftw as external project on Windows")
        get_filename_component(VS_BIN_DIR "${CMAKE_LINKER}" DIRECTORY)

        file(DOWNLOAD ftp://ftp.fftw.org/pub/fftw/fftw-3.3-libs-visual-studio-2010.zip ./fftw-prefix/src/project.zip)

        if (ARCHITECTURE MATCHES "x64")
            set(FFTW_ARCHITECTURE "64")
        else ()
            set(FFTW_ARCHITECTURE "32")
            set(ARCHITECTURE "x86")
        ENDIF ()
        message("detected architecture set for ARCHITECTURE  ${ARCHITECTURE} for FFTW_ARCHITECTURE  ${FFTW_ARCHITECTURE}")

        set(FFTW_BUILD_RELATIVE_PATH "../fftw/fftw-3.3-libs/libfftw-3.3")
        ExternalProject_add(fftw
                URL "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.3.zip"
                CONFIGURE_COMMAND ""
                BUILD_COMMAND
                    COMMAND ${CMAKE_COMMAND} -E tar xfv ../project.zip --format=zip
                    COMMAND ${CMAKE_COMMAND} -E rename ./fftw-3.3-libs ../fftw/fftw-3.3-libs
                    COMMAND ${CMAKE_VS_MSBUILD_COMMAND} ${FFTW_BUILD_RELATIVE_PATH}/libfftw-3.3.vcxproj /p:PlatformToolset=v140 /p:configuration=Debug /p:platform=${ARCHITECTURE}
                    COMMAND ${CMAKE_VS_MSBUILD_COMMAND} ${FFTW_BUILD_RELATIVE_PATH}/libfftw-3.3.vcxproj /p:PlatformToolset=v140 /p:configuration=Release /p:platform=${ARCHITECTURE}
                INSTALL_COMMAND
                    COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_BUILD_RELATIVE_PATH}/Release/libfftw-3.3.lib ${INSTALL_LIB_DIR}/fftw/Release/libfftw-3.3.lib
                    COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_BUILD_RELATIVE_PATH}/Release/libfftw-3.3.dll ${INSTALL_BIN_DIR}/fftw/Release/libfftw-3.3.dll
                    COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_BUILD_RELATIVE_PATH}/Debug/libfftw-3.3.lib ${INSTALL_LIB_DIR}/fftw/Debug/libfftw-3.3.lib
                    COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_BUILD_RELATIVE_PATH}/Debug/libfftw-3.3.dll ${INSTALL_BIN_DIR}/fftw/Debug/libfftw-3.3.dll
                    COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_BUILD_RELATIVE_PATH}/../api/fftw3.h ${INSTALL_INC_DIR}/fftw3.h
                )

    else ()

    set(FFTW_OPTIMIZATION_CONFIGURATION "--enable-sse2 --enable-avx --enable-openmp" CACHE INTERNAL "architecture flags: --enable-sse --enable-sse2 --enable-altivec --enable-mips-ps --enable-cell")

    if (${CMAKE_CURRENT_BINARY_DIR} MATCHES ".*[ \t].*")
        message(FATAL_ERROR
                "Can't build fftw in a directory with whitespace in its name")
    endif ()
    #
    # build fftw as an external project
    if (BUILD_SHARED_LIBS)
        set(FFTW_SHARED_FLAG --enable-shared)
    endif ()

    ExternalProject_add(fftwf
            PREFIX fftwf
            URL "http://www.fftw.org/fftw-3.3.4.tar.gz"
            CONFIGURE_COMMAND
            env
            "CC=${CMAKE_C_COMPILER} ${CMAKE_C_COMPILER_ARG1}"
            "CFLAGS=${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}"
            "LDFLAGS=$ENV{LDFLAGS}"
            "LIBS=$ENV{LIBS}"
            "CPP=$ENV{CPP}"
            "CPPFLAGS=$ENV{CPPFLAGS}"
            ${ITK_BINARY_DIR}/fftwf/src/fftwf/configure
            ${FFTW_SHARED_FLAG}
            ${FFTW_OPTIMIZATION_CONFIGURATION}
            ${FFTW_THREADS_CONFIGURATION}
            --disable-fortran
            --enable-float
            --prefix=${INSTALL_DIR}
            )
    endif ()

else ()
    message(STATUS "fftw3 already installed!")
endif ()
