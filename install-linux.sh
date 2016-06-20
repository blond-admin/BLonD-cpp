#!/bin/bash

echo -e "Installing necessary libraries...\n"

BLOND_HOME=$(pwd)

mkdir -p tmp
mkdir -p external/opt/lib
mkdir -p external/opt/include
OPT="${BLOND_HOME}/external/opt"

INSTALL_FFTW=true
INSTALL_GTEST=true

if [ -e ${OPT}/include/fftw3.h ] && [ -e ${OPT}/lib/libfftw3.a ]; then
   echo -e "---- Looks like fftw3 is already installed,"
   echo -e "----  are you sure you want to reinstall it?"
   select yn in "Yes" "No"; do
      case $yn in
         Yes ) INSTALL_FFTW=true; break;;
         No ) INSTALL_FFTW=false; break;;
      esac
   done
fi


if [ "${INSTALL_FFTW}" = "true" ] ; then
   echo -e "\n\n---- Installing fftw3\n\n"
   wget www.fftw.org/fftw-3.3.4.tar.gz -Otmp/fftw3.tar.gz
   tar -xzvf tmp/fftw3.tar.gz -Cexternal
   cd external/fftw-3.3.4
   ./configure --enable-openmp --prefix="${BLOND_HOME}/external/opt"
   make 
   make install

   if [ -e ${OPT}/include/fftw3.h ] && [ -e ${OPT}/lib/libfftw3.a ]; then
      echo -e "\n\n---- fftw3 is successfully installed\n\n"
   else
      echo -e "\n\n---- fftw3 has failed to install successfully"
      echo -e "---- You will have to manually install this library\n\n"
   fi
fi

cd ${BLOND_HOME}

if [ -e ${OPT}/include/gtest/gtest.h ] && [ -e ${OPT}/lib/libgtest.a ] \
   && [ -e ${OPT}/lib/libgtest_main.a ]; then
   echo -e "---- Looks like googletest is already installed,"
   echo -e "---- are you sure you want to reinstall it?"
   select yn in "Yes" "No"; do
      case $yn in
         Yes ) INSTALL_GTEST=true; break;;
         No ) INSTALL_GTEST=false; break;;
      esac
   done
fi

if [ "${INSTALL_GTEST}" = "true" ] ; then

   echo -e "\n\n---- Installing googletest\n\n"

   git clone https://github.com/google/googletest.git external/googletest
   cd external/googletest/googletest
   cp -r include/* "${OPT}/include/"
   mkdir -p build
   cd build && cmake .. && make
   cp *.a "${OPT}/lib"

   cd ${BLOND_HOME}

   if [ -e ${OPT}/include/gtest/gtest.h ] && [ -e ${OPT}/lib/libgtest.a ] \
      && [ -e ${OPT}/lib/libgtest_main.a ]; then
      echo -e "\n\n---- Googletest is successfully installed\n\n"
   else
      echo -e "\n\n---- Googletest has failed to install successfully"
      echo -e "---- You will have to manually install this library\n\n"
   fi
fi

rm -r tmp &> /dev/null
