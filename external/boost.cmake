if (WITH_COMPUTE)
    message(STATUS "searching for ${INSTALL_INC_DIR}/boost/compute.hpp")
    if ((NOT EXISTS ${INSTALL_INC_DIR}/boost/compute.hpp))
        set(BUILD_COMPUTE "True")
    endif ()
endif ()


if (BUILD_COMPUTE)
    message(STATUS "not found, building boost.compute")
    set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/build/boost/")

    message("External project - Boost")

    set(Boost_Bootstrap_Command)
    if (WIN32 AND NOT MINGW)
        set(Boost_Bootstrap_Command bootstrap.bat)
        set(Boost_b2_Command b2.exe)
    else ()
        set(Boost_Bootstrap_Command ./bootstrap.sh)
        set(Boost_b2_Command ./b2)
    endif ()

    if (ARCHITECTURE MATCHES "x64")
        set(BOOST_ARCHITECTURE "64")
    else ()
        set(BOOST_ARCHITECTURE "32")
        set(ARCHITECTURE "x86")
    ENDIF ()

    set(Boost_Version 1_61_0)
    set(Boost_Version_s 1_61)
    ExternalProject_Add(boost
            URL http://kent.dl.sourceforge.net/project/boost/boost/1.61.0/boost_${Boost_Version}.zip
            DOWNLOAD_DIR ${BOOST_ROOT}
            BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND
            COMMAND ${Boost_Bootstrap_Command}
            BUILD_COMMAND
            COMMAND ${Boost_b2_Command}
            --with-thread
            --prefix=${INSTALL_DIR} -d0
            address-model=32_64 --threading=multi --link=shared,static --variant=debug,release install -j8
            INSTALL_COMMAND
                COMMAND ${CMAKE_COMMAND} -E echo "Buist compiled and build OK"
            )


    if (NOT WIN32)
        set(Boost_LIBRARY_DIR ${INSTALL_DIR}/lib/boost/)
        set(Boost_INCLUDE_DIR ${INSTALL_DIR}/include/)
    else ()
        execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${INSTALL_INC_DIR}/boost-${Boost_Version_s}/boost ${INSTALL_INC_DIR}/boost )
        set(Boost_LIBRARY_DIR ${INSTALL_DIR}/lib/)
        set(Boost_INCLUDE_DIR ${INSTALL_DIR}/include/boost-${Boost_Version}/)
    endif ()

endif ()