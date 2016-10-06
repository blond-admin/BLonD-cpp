#!/bin/bash

echo -e "Installing necessary libraries...\n"

BLOND_HOME="$(pwd)"
EXTERNAL="${BLOND_HOME}/external"
INSTALL="${EXTERNAL}/install"
log="${EXTERNAL}/log.out"

mkdir -p ${EXTERNAL}/tmp
mkdir -p ${INSTALL}/lib
mkdir -p ${INSTALL}/include
echo > ${log}

INSTALL_FFTW=true
INSTALL_GTEST=true
INSTALL_PYTHON=true
INSTALL_HDF5=true
INSTALL_GBENCH=false
SKIP_QUESTIONS=false

while getopts ":f:t:p:h:b:s:" opt; do
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


if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.la ]; then
   echo -e "\n\n---- Looks like fftw3 is already installed,"
   if [ "${SKIP_QUESTIONS}" = "true" ] ; then
      echo -e "---- fftw3 installation will be skipped"
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
   wget www.fftw.org/fftw-3.3.4.tar.gz -O${EXTERNAL}/tmp/fftw3.tar.gz
   tar -xzvf ${EXTERNAL}/tmp/fftw3.tar.gz -C${EXTERNAL}/tmp &>> $log
   cd ${EXTERNAL}/tmp/fftw-3.3.4
   ./configure --disable-alloca \
               --disable-fortran \
               --disable-static \
               --enable-shared \
               --enable-openmp \
               --enable-sse2 \
               --enable-avx \
               --with-our-malloc \
               --with-incoming-stack-boundary=2 \
               --prefix="${INSTALL}" &>> $log
               # --enable-threads \
               # --with-combined-threads \
   # ./configure --enable-openmp --prefix="${INSTALL}" &>> $log
   make &>> $log
   make install &>> $log

   if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.a ]; then
      echo -e "---- fftw3 has been installed successfully\n\n"
   else
      echo -e "---- fftw3 has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"

   fi
   echo -e "---- Installation of fftw3 is completed\n\n"
fi

# -----------------------
# end of fftw intallation
# -----------------------


cd ${BLOND_HOME}


# -----------------
# HDF5 installation
# -----------------

if [ -e ${INSTALL}/include/hdf5.h ] && [ -e ${INSTALL}/lib/libhdf5.so ]; then
   echo -e "\n\n---- Looks like HDF5 is already installed,"
   if [ "${SKIP_QUESTIONS}" = "true" ] ; then
      echo -e "---- HDF5 installation will be skipped"
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
   wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.17.tar.gz -O${EXTERNAL}/tmp/hdf5.tar.gz
   tar -xzvf ${EXTERNAL}/tmp/hdf5.tar.gz -C${EXTERNAL}/tmp &>> $log
   cd ${EXTERNAL}/tmp/hdf5-1.8.17
   ./configure --enable-cxx \
               --enable-hl \
               --enable-static  \
               --enable-shared \
               --enable-fast-install \
               --prefix="${INSTALL}" &>> $log
   make &>> $log
   make install &>> $log

if [ -e ${INSTALL}/include/hdf5.h ] && [ -e ${INSTALL}/lib/libhdf5.so ]; then
      echo -e "---- HDF5 has been installed successfully\n\n"
   else
      echo -e "---- HDF5 has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"

   fi
   echo -e "---- Installation of HDF5 is completed\n\n"
fi

# -----------------------
# end of HDF5 intallation
# -----------------------


cd ${BLOND_HOME}

# -----------------------
# googletest installation
# -----------------------

if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.a ] \
   && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
   echo -e "\n\n---- Looks like googletest is already installed,"
   if [ "${SKIP_QUESTIONS}" = "true" ] ; then
      echo -e "---- googletest installation will be skipped"
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

   git clone https://github.com/google/googletest.git ${EXTERNAL}/tmp/googletest
   cd ${EXTERNAL}/tmp/googletest/googletest
   cp -r include/* "${INSTALL}/include/"
   mkdir -p build &>> $log
   cd build
   cmake .. &>> $log
   make &>> $log
   cp *.a "${INSTALL}/lib"

   cd ${BLOND_HOME}

   if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.a ] \
      && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
      echo -e "---- Googletest has been installed successfully\n\n"
   else
      echo -e "---- Googletest has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"

   fi
   echo -e "---- Installation of googletest is completed\n\n"
fi

# ------------------------------
# end of googletest installation
# ------------------------------

cd ${BLOND_HOME}



# -----------------------
# googlebench installation
# -----------------------


if [ "${INSTALL_GBENCH}" = "true" ] ; then

   echo -e "\n\n---- Installing googlebench"

   git clone https://github.com/google/benchmark.git ${EXTERNAL}/tmp/googlebench
   cd ${EXTERNAL}/tmp/googlebench/
   cp -r include/* "${INSTALL}/include/"
   mkdir -p build &>> $log
   cd build
   CC=gcc CXX=g++ cmake -DCMAKE_BUILD_TYPE=Release \
                        -DBENCHMARK_CXX_LINKER_FLAGS="-lrt" \
                        .. &>> $log
   make &>> $log
   cp src/*.a "${INSTALL}/lib"

   cd ${BLOND_HOME}

   if [ -e ${INSTALL}/include/benchmark/benchmark.h ] && [ -e ${INSTALL}/lib/libbenchmark.a ]; then
      echo -e "---- googlebench has been installed successfully\n\n"
   else
      echo -e "---- googlebench has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
   fi
   echo -e "---- Installation of googlebench is completed\n\n"
fi

# ------------------------------
# end of googlebench installation
# ------------------------------


cd ${BLOND_HOME}

# -------------------
# Python installation
# -------------------
# PY_VERSION=$(which python &>/dev/null && python --version 2>&1 || false)
PY_VERSION=$(which python 2>/dev/null)
if [ -z ${PY_VERSION} ]; then

   if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/python2.7/config/libpython2.7.a ]; then
      echo -e "\n\n---- Looks like Python2.7 is already installed,"
      if [ "${SKIP_QUESTIONS}" = "true" ] ; then
         echo -e "---- Python2.7 installation will be skipped"
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
      wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz -O${EXTERNAL}/tmp/Python-2.7.12.tgz
      tar -xzvf ${EXTERNAL}/tmp/Python-2.7.12.tgz -C${EXTERNAL}/ &>> $log
      cd ${EXTERNAL}/Python-2.7.12
      ./configure --enable-unicode=ucs4 \
                  --with-threads \
                  --enable-shared \
                  --prefix="$(pwd)" &>> $log
      make &>> $log
      make install &>> $log

      # export PATH="${INSTALL}/bin:$PATH"
      cd ${BLOND_HOME}

      if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/python2.7/config/libpython2.7.a ]; then
         echo -e "---- Python has been installed successfully\n\n"
         PYTHON="${INSTALL}/Python-2.7.12/bin/python"
      else
         echo -e "\n\n---- Python has failed to install successfully"
         echo -e "---- You will have to manually install this library"
         echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
      fi

      #echo -e "---- Installing of Python2.7 is completed\n\n"
   fi

else
   PYTHON="python"
fi


# --------------------------
# end of Python installation
# --------------------------


# -----------------------
# Python external modules installation
# -----------------------
PIP_INSTALLED=`echo $?`
if [ "$PIP_INSTALLED" == "1" ]; then
   echo -e "\n\n---- PIP is needed in order to install required python modules"
   echo -e "---- If you are on Fedora/CentOS/RHEL try: yum install python-pip"
   echo -e "---- If you are on Debian/Ubuntu try: apt-get install python-pip"
   echo -e "---- and then re-run this script."
   echo -e "---- For more information, please visit this site: https://packaging.python.org/install_requirements_linux/\n\n"
else
   echo -e "\n\n---- Installing Python's external modules..."
   ${PYTHON} -m pip install --user virtualenv 2>> $log
   # pip install virtualenv &>> $log
   ${PYTHON} -m virtualenv --python=${PYTHON} ${INSTALL} &>> $log
   source ${INSTALL}/bin/activate &>> $log
   ${PYTHON} -m pip install -r ${EXTERNAL}/python-packages.txt 2>> $log
   export PYTHONPATH="${BLOND_HOME}/python:$PYTHONPATH"
   echo -e "\n\n---- Python's external modules have been installed successfully\n\n"
fi
# ----------------------------------
# end of Python Modules installation
# ----------------------------------



echo -e "\n\n---- Exteranl depedencies installation procedure completed"
echo -e "---- You may now try to build the project :)"
echo -e "---- You can consult the ${log} file for any errors\n\n"

rm -rf ${EXTERNAL}/tmp &> /dev/null
