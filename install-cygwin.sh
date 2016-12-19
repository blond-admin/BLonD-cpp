#!/bin/bash

verlte() {
    [  "$1" = "`echo -e "$1\n$2" | sort -V | head -n1`" ]
}

verlt() {
    [ "$1" = "$2" ] && return 1 || verlte $1 $2
}

LIB="/usr/lib"
BIN="/usr/bin"
PYTHON="/usr/bin/python"
INCLUDE="/usr/include"
BLOND_HOME="$(pwd)"
EXTERNAL="${BLOND_HOME}/external"
INSTALL="${EXTERNAL}/install"
log="${EXTERNAL}/log.out"

mkdir -p ${EXTERNAL}/tmp
mkdir -p ${INSTALL}/lib
mkdir -p ${INSTALL}/include
echo > ${log}


echo -e "---- BLonD++ on Windows auto-configuration script\n"
echo -e "---- Start of Part1: Cygwin packages installation\n\n"

echo -e "---- Installing apt-cyg package manager for Cygwin\n"
wget rawgit.com/transcode-open/apt-cyg/master/apt-cyg -O${EXTERNAL}/tmp/apt-cyg
install ${EXTERNAL}/tmp/apt-cyg ${BIN}
echo -e "---- apt-cyg has been installed successfully\n\n"

echo -e "---- Installing all the necessary packages\n"

echo -e "---- This is going to take some time (15'-30') and some space (3-4GB)\n"
apt-cyg install $(cat ${EXTERNAL}/cygwin-packages.txt)
echo -e "---- Cygwin's packages has been installed successfully\n\n"

echo -e "---- End of Part1: Cygwin packages installation\n\n"

echo -e "---- Start of Part2: External depedencies installation\n\n"


INSTALL_FFTW=true
INSTALL_GTEST=true
INSTALL_PYTHON=true
INSTALL_HDF5=true
INSTALL_GBENCH=false
SKIP_QUESTIONS=false

while getopts ":f:t:p:h:b:s" opt; do
    case $opt in 
        f)
            INSTALL_FFTW=$OPTARG
            ;;
        t)
            INSTALL_GTEST=$OPTARG
            ;;
        p)
            INSTALL_PYTHON=$OPTARG
            ;;
        h)
            INSTALL_HDF5=$OPTARG
            ;;
        b)
            INSTALL_GBENCH=$OPTARG
            # echo "-b specified with arg $OPTARG"
            ;;
        s) SKIP_QUESTIONS=true
            ;;
        :)
            echo -e "Option -$OPTARG requires and argument"
            ;;
        \?)
            echo -e "Invalid option -$OPTARG"
            ;;
    esac
done


# -----------------
# fftw installation
# -----------------


if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.dll.a ]; then
    echo -e "\n\n---- Looks like fftw3 is already installed,"
    if [ "${SKIP_QUESTIONS}" = "true" ] ; then
        echo -e "---- fftw3 installation will be skipped"
        INSTALL_FFTW=false
    else
        echo -e "---- are you sure you want to reinstall it?"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) INSTALL_FFTW=true; break;;
                No ) INSTALL_FFTW=false; break;;
            esac
        done
    fi
fi


if [ "${INSTALL_FFTW}" = "true" ] ; then
    echo -e "\n\n---- Installing fftw3"
    cp -l ${LIB}/libfftw3.dll.a ${INSTALL}/lib/
    cp -l ${INCLUDE}/fftw3.h ${INSTALL}/include/

    if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.dll.a ]; then
        echo -e "---- fftw3 has been installed successfully\n\n"
    else
        echo -e "---- fftw3 has failed to install successfully"
        echo -e "---- You will have to manually install this library"
        echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"

    fi
    # echo -e "---- Installation of fftw3 is completed\n\n"
fi

# -----------------------
# end of fftw intallation
# -----------------------


cd ${BLOND_HOME}


# -----------------
# HDF5 installation
# -----------------

if [ -e ${INSTALL}/include/hdf5.h ] && \
    [ -e ${INSTALL}/lib/libhdf5.dll.a ] && \
    [ -e ${INSTALL}/lib/libhdf5_cpp.dll.a ] && \
    [ -e ${INSTALL}/lib/libhdf5_hl.dll.a ] && \
    [ -e ${INSTALL}/lib/libhdf5_hl_cpp.dll.a ]
then
    echo -e "\n\n---- Looks like HDF5 is already installed,"
    if [ "${SKIP_QUESTIONS}" = "true" ] ; then
        echo -e "---- HDF5 installation will be skipped"
        INSTALL_HDF5=false
    else
        echo -e "---- are you sure you want to reinstall it?"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) INSTALL_HDF5=true; break;;
                No ) INSTALL_HDF5=false; break;;
            esac
        done
    fi
fi


if [ "${INSTALL_HDF5}" = "true" ]; then
    echo -e "\n\n---- Installing HDF5"
    cp -l ${LIB}/libhdf5.dll.a ${INSTALL}/lib/
    cp -l ${LIB}/libhdf5_cpp.dll.a ${INSTALL}/lib/
    cp -l ${LIB}/libhdf5_hl.dll.a ${INSTALL}/lib/
    cp -l ${LIB}/libhdf5_hl_cpp.dll.a ${INSTALL}/lib/
    cp -l ${INCLUDE}/H5*.h ${INSTALL}/include/
    cp -l ${INCLUDE}/hdf5.h ${INSTALL}/include/
    cp -l ${INCLUDE}/hdf5_hl.h ${INSTALL}/include/

    if [ -e ${INSTALL}/include/hdf5.h ] && \
        [ -e ${INSTALL}/lib/libhdf5.dll.a ] && \
        [ -e ${INSTALL}/lib/libhdf5_cpp.dll.a ] && \
        [ -e ${INSTALL}/lib/libhdf5_hl.dll.a ] && \
        [ -e ${INSTALL}/lib/libhdf5_hl_cpp.dll.a ]
    then
        echo -e "---- HDF5 has been installed successfully\n\n"
    else
        echo -e "---- HDF5 has failed to install successfully"
        echo -e "---- You will have to manually install this library"
        echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
    fi
    # echo -e "---- Installation of HDF5 is completed\n\n"
fi

# -----------------------
# end of HDF5 intallation
# -----------------------


cd ${BLOND_HOME}

# -----------------------
# googletest installation
# -----------------------

if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.dll.a ] \
    && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
    echo -e "\n\n---- Looks like googletest is already installed,"
    if [ "${SKIP_QUESTIONS}" = "true" ] ; then
        echo -e "---- googletest installation will be skipped"
        INSTALL_GTEST=false
    else
        echo -e "---- are you sure you want to reinstall it?"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) INSTALL_GTEST=true; break;;
                No ) INSTALL_GTEST=false; break;;
            esac
        done
    fi
fi

if [ "${INSTALL_GTEST}" = "true" ] ; then

    echo -e "\n\n---- Installing googletest"

    git clone --branch=master https://github.com/google/googletest.git ${EXTERNAL}/tmp/googletest
    cd ${EXTERNAL}/tmp/googletest/googletest
    cp -rl include/* "${INSTALL}/include/"
    mkdir -p build &>> $log
    cd build
    cmake .. &>> $log
    make &>> $log
    cp -l libgtest_main.a "${INSTALL}/lib/libgtest_main.a"
    cp -l libgtest.a "${INSTALL}/lib/libgtest.dll.a"

    cd ${BLOND_HOME}

    if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.dll.a ] \
        && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
        echo -e "---- Googletest has been installed successfully\n\n"
    else
        echo -e "---- Googletest has failed to install successfully"
        echo -e "---- You will have to manually install this library"
        echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
    fi
    # echo -e "---- Installation of googletest is completed\n\n"
fi

# ------------------------------
# end of googletest installation
# ------------------------------

cd ${BLOND_HOME}



# -----------------------
# googlebench installation
# -----------------------


# -------------------
# Python installation
# -------------------

if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/libpython2.7.dll.a ]; then
    echo -e "\n\n---- Looks like Python2.7 is already installed,"
    if [ "${SKIP_QUESTIONS}" = "true" ] ; then
        echo -e "---- Python2.7 installation will be skipped"
        INSTALL_PYTHON=false
    else
        echo -e "---- are you sure you want to reinstall it?"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) INSTALL_PYTHON=true; break;;
                No ) INSTALL_PYTHON=false; break;;
            esac
        done
    fi
fi

if [ "${INSTALL_PYTHON}" = "true" ] ; then
    echo -e "\n\n---- Installing Python2.7\n\n"
    cp -l ${LIB}/libpython2.7.dll.a ${INSTALL}/lib/
    cp -rl ${LIB}/python2.7 ${INSTALL}/lib
    cp -rl ${INCLUDE}/python2.7 ${INSTALL}/include

    cd ${BLOND_HOME}

    if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/libpython2.7.dll.a ]; then
        echo -e "---- Python has been installed successfully\n"
        ${PYTHON} -m ensurepip
        echo -e "---- Pip has been installed successfully\n\n"
    else
        echo -e "\n\n---- Python has failed to install successfully"
        echo -e "---- You will have to manually install this library"
        echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
    fi

fi



# --------------------------
# end of Python installation
# --------------------------


# -----------------------
# Python external modules installation
# -----------------------
PIP_INSTALLED=$(${PYTHON} -m pip -V | echo $?)
if [ "$PIP_INSTALLED" == "1" ]; then
    echo -e "\n\n---- PIP is needed in order to install the required python modules\n\n"
else
    echo -e "\n\n---- Installing Python's external modules..."
    ${PYTHON} -m pip install -r ${EXTERNAL}/python-packages.txt
    echo -e "\n\n---- Python's external modules have been installed successfully\n\n"
fi
# ----------------------------------
# end of Python Modules installation
# ----------------------------------




echo -e "\n\n---- Exteranl depedencies installation procedure completed"
echo -e "---- You may now try to build the project :)"
echo -e "---- You can consult the ${log} file for any errors\n"

echo -e "---- End of Part2: External depedencies installation\n\n"

rm -rf ${EXTERNAL}/tmp &> /dev/null
