#OpenCL
if (WITH_COMPUTE)
    add_definitions(-DWITH_OPENCL=1)
    #OpenCL is required for GPU computations
    find_package(OpenCL)
    if (OpenCL)
        include_directories(${OpenCL_INCLUDE_DIRS})
        link_directories(${OpenCL_LIBRARY})
        LIST(APPEND LIBRARIES
                ${OpenCL_LIBRARY}
                )
    elseif (OPENCL_SDK_PATH)
        if (WIN32 AND NOT MINGW)
            if (${ARCHITECTURE} STREQUAL "x64")
                message(status "ARCHITECTURE ${ARCHITECTURE}")
                link_directories("${OPENCL_SDK_PATH}/lib/x64")
                link_directories("${OPENCL_SDK_PATH}/lib/x86_64")
            else ()
                link_directories("${OPENCL_SDK_PATH}/lib/Win32")
                link_directories("${OPENCL_SDK_PATH}/lib/x86")
                link_directories("${OPENCL_SDK_PATH}/lib/")
            endif ()
            set(OPENCL_LIB "OpenCL.lib")
        elseif (MINGW)
            LIST(APPEND LIBRARIES
                    "C:/Windows/System32/OpenCL.dll"
                    )
        else ()
            message(STATUS "ARCHITECTURE ${ARCHITECTURE}")
            if (ARCHITECTURE STREQUAL "x64")
                find_library(OPENCL_LIB
                        NAMES cl CL OpenCL OpenCL.lib libOpenCL libOpenCL.so
                        PATHS "${OPENCL_SDK_PATH}" "${OPENCL_SDK_PATH}/lib/x64" "${OPENCL_SDK_PATH}/lib/x86_64"
                        NO_DEFAULT_PATH
                        )
            else ()
                find_library(OPENCL_LIB
                        NAMES cl CL OpenCL OpenCL.lib libOpenCL libOpenCL.so
                        PATHS "${OPENCL_SDK_PATH}" "${OPENCL_SDK_PATH}/lib/Win32" "${OPENCL_SDK_PATH}/lib/x86" "${OPENCL_SDK_PATH}/lib/"
                        NO_DEFAULT_PATH
                        )
            endif ()
        endif ()
        message("Note: OPENCL_LIB ${OPENCL_LIB}")
        message("Note: OPENCL_SDK_PATH ${OPENCL_SDK_PATH}")

        include_directories("${OPENCL_SDK_PATH}/include")

        get_filename_component(OPENCL_DIR ${OPENCL_LIB} DIRECTORY)
        link_directories(${OPENCL_DIR})
        LIST(APPEND LIBRARIES
                ${OPENCL_LIB}
                )
    endif ()
endif ()