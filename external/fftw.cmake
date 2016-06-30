#based on https://github.com/hinerm/ITK/blob/master/CMake/itkExternal_FFTW.cmake
# Encapsulates building FFTW as an External Project.
if (WITH_FFTW)
    set(fftwLib "fftw")
    GetnerateLibName(${fftwLib})
    message(STATUS "searching for ${INSTALL_LIB_DIR}${fftwLib}")
    if ((NOT EXISTS ${INSTALL_LIB_DIR}${fftwLib}))
        set(BUILD_FFTW "True")
    endif ()
endif ()

if (BUILD_FFTW)

    if (WIN32 AND NOT MINGW)
        message(STATUS "Can't build fftw as external project on Windows, downloading official release")
        message(STATUS "Can't build fftw as external project on Windows")
        get_filename_component(VS_BIN_DIR "${CMAKE_LINKER}" DIRECTORY)

        file(DOWNLOAD ftp://ftp.fftw.org/pub/fftw/fftw-3.3-libs-visual-studio-2010.zip ./fftw-prefix/src/project.zip)

        if (ARCHITECTURE MATCHES "x64")
            set(FFTW_ARCHITECTURE "64")
        else()
            set(FFTW_ARCHITECTURE "32")
            set(ARCHITECTURE  "x86")
        ENDIF()
        message("detected architecture set for ARCHITECTURE  ${ARCHITECTURE} for FFTW_ARCHITECTURE  ${FFTW_ARCHITECTURE}")


            ExternalProject_add(fftw
                    URL "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.3.zip"
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND
                    echo %cd%
                    ${CMAKE_COMMAND} -E tar "xfv" "./project.zip" --format=zip
                    ${CMAKE_COMMAND} -E rename ./fftw-3.3-libs ./fftw/fftw-3.3-libs
                    #${VS_BIN_DIR}/lib /machine:x${FFTW_MACHINE_ARCHITECTURE} /def:libfftw3-3.def
                    ${CMAKE_VS_MSBUILD_COMMAND} ./fftw/fftw-3.3-libs/fftw-3.3-libs.sln /p:PlatformToolset=v140 /p:configuration=Debug /p:platform "${ARCHITECTURE}" /project "libfftw-3.3"
                    ${CMAKE_VS_MSBUILD_COMMAND} ./fftw/fftw-3.3-libs/fftw-3.3-libs.sln /p:PlatformToolset=v140 /p:configuration=Release /p:platform "${ARCHITECTURE}" /project "libfftw-3.3"
                        #cmake -E fftw3.h ${INSTALL_DIR}/bin/
                        #cmake -E libfftw3-3.lib ${INSTALL_DIR}/lib/
                        #cmake -E fftw3.h ${INSTALL_DIR}/include/
                    INSTALL_COMMAND ""
                    )
        return() # todo remove!

    else ()

        endif ()
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
endif ()