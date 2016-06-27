cmake_minimum_required(VERSION 2.8)
project(externals)
include(ExternalProject)
add_custom_target(external) 

#Log required variables
mseeage(STATUS "Building external packages")
mseeage( STATUS "into ${INSTALL_DIR}" )
mseeage( STATUS "WITH_GSL ${WITH_GSL}, GSL_PATCH_PATH ${GSL_PATCH_PATH} " )
mseeage( STATUS "WITH_GOOGLETEST ${WITH_GOOGLETEST}" )
mseeage( STATUS "WITH_BENCHMARK ${WITH_BENCHMARK}" )
mseeage( STATUS "WITH_FFTW ${WITH_FFTW}" )



#Check if required libs and includes are already presented
set(INSTALL_LIB_DIR "${INSTALL_DIR}/lib/")


if(WITH_FFTW)
    set(fftwLib "fftw")
    GetnerateLibName(${fftwLib})
    if( (NOT EXISTS ${INSTALL_LIB_DIR}${fftwLib}))
        set(BUILD_FFTW "True")
    endif()
endif()







if(BUILD_FFTW)

endif()
