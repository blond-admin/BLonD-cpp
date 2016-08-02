#!/bin/bash

echo -e "Installing necessary libraries...\n"

BLOND_HOME=$(pwd)

mkdir -p external/tmp
mkdir -p external/install/lib
mkdir -p external/install/include
INSTALL="${BLOND_HOME}/external/install"
log="${BLOND_HOME}/external/log.out"
touch ${log}

INSTALL_FFTW=true
INSTALL_GTEST=true
INSTALL_PYTHON=true

# -----------------
# fftw installation
# -----------------

if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.a ]; then
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
   wget www.fftw.org/fftw-3.3.4.tar.gz -Otmp/fftw3.tar.gz 2>> $log
   tar -xzvf external/tmp/fftw3.tar.gz -Cexternal 2>> $log
   cd external/fftw-3.3.4 2>> &log
   ./configure --enable-openmp --prefix="${BLOND_HOME}/external/install" 2>> $log
   make 2>> $log
   make install 2>> $log

   if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.a ]; then
      echo -e "\n\n---- fftw3 is successfully installed\n\n"
   else
      echo -e "\n\n---- fftw3 has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"

   fi
fi

# -----------------------
# end of fftw intallation
# -----------------------


cd ${BLOND_HOME}

# -----------------------
# googletest installation
# -----------------------

if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.a ] \
   && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
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

   git clone https://github.com/google/googletest.git external/googletest 2>> $log
   cd external/googletest/googletest 2>> $log
   cp -r include/* "${INSTALL}/include/" 2>> $log
   mkdir -p build 2>> $log
   cd build && cmake .. && make 2>> $log
   cp *.a "${INSTALL}/lib" 2>> $log

   cd ${BLOND_HOME}

   if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.a ] \
      && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
      echo -e "\n\n---- Googletest is successfully installed\n\n"
   else
      echo -e "\n\n---- Googletest has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"

   fi
fi

# ------------------------------
# end of googletest installation
# ------------------------------

cd ${BLOND_HOME}

# -------------------
# Python installation
# -------------------

if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/python2.7/config/libpython2.7.a ]; then
   echo -e "---- Looks like python is already installed,"
   echo -e "---- are you sure you want to reinstall it?"
   select yn in "Yes" "No"; do
      case $yn in
         Yes ) INSTALL_PYTHON=true; break;;
         No ) INSTALL_PYTHON=false; break;;
      esac
   done
fi

if [ "${INSTALL_PYTHON}" = "true" ] ; then
   echo -e "\n\n---- Installing python\n\n"
   wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz -O"${BLOND_HOME}/externall/tmp" 2>> $log
   tar -xzvf external/tmp/Python-2.7.12.tgz -C"${BLOND_HOME}/external" 2>> $log
   cd external/Python-2.7.12 2>> $log
   ./configure --enable-unicode=ucs4 --prefix="${BLOND_HOME}/external/install" 2>> $log
   make 2>> $log
   make install 2>> $log

   cd ${BLOND_HOME}

   if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/python2.7/config/libpython2.7.a ]; then
      echo -e "\n\n---- Python is successfully installed\n\n"
   else
      echo -e "\n\n---- Python has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
   fi
fi

# --------------------------
# end of Python installation
# --------------------------


# ---------------------------
# Python Modules installation
# ---------------------------
PYTHON_MODULES=( "numpy" "scipy" "matplotlib")

cd ${BLOND_HOME}
PIP_INSTALLED=`which pip`

if [ -z "$PIP_INSTALLED" ]; then
   echo -e "\n\n---- PIP is needed in order to install required python modules"
   echo -e "---- If you are on Fedora/CentOS/RHEL try: yum install python-pip"
   echo -e "---- If you are on Debian/Ubuntu try: apt-get install python-pip"
   echo -e "---- and then re-run this script."
   echo -e "---- For more information, please visit this site: https://packaging.python.org/install_requirements_linux/ \n\n"
else
   for module in "${PYTHON_MODULES[@]}"; do
      echo -e "\n\n---- Installing ${module}\n\n"
      pip install --upgrade --target="${BLOND_HOME}/external/lib/python2.7/site-packages" ${module} 2>> $log
   done
fi

# ----------------------------------
# end of Python Modules installation
# ----------------------------------


echo -e "\n\n---- Exteranl depedencies installation procedure completed"
echo -e "---- Now try to build the project"
echo -e "---- You can consult the ${log} file for any errors\n\n"

rm -rf external/tmp &> /dev/null
