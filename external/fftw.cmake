#based on https://github.com/hinerm/ITK/blob/master/CMake/itkExternal_FFTW.cmake
# Encapsulates building FFTW as an External Project.

if (WITH_FFTW)
    set(fftwLib "fftw")
    GetnerateLibName(${fftwLib})
    if ((NOT EXISTS ${INSTALL_LIB_DIR}${fftwLib}))
        set(BUILD_FFTW "True")
    endif ()
endif ()

if (BUILD_FFTW)

    if (WIN32 AND NOT MINGW)
        message(STATUS "Can't build fftw as external project on Windows, downloading official release")
        message("Can't build fftw as external project on Windows")
        if (ARCHITECTURE MATCHES "x86")
            ExternalProject_add(fftw
                    URL "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4-dll32.zip"
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND ""
                    INSTALL_COMMAND ""
                    )
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