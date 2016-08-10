#!/bin/bash

echo -e "Installing necessary libraries...\n"

BLOND_HOME=$(pwd)
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

# -----------------
# fftw installation
# -----------------

if [ -e ${INSTALL}/include/fftw3.h ] && [ -e ${INSTALL}/lib/libfftw3.la ]; then
   echo -e "\n\n---- Looks like fftw3 is already installed,"
   echo -e "---- are you sure you want to reinstall it?\n\n"
   select yn in "Yes" "No"; do
      case $yn in
         Yes ) INSTALL_FFTW=true; break;;
         No ) INSTALL_FFTW=false; break;;
      esac
   done
fi


if [ "${INSTALL_FFTW}" = "true" ] ; then
   echo -e "\n\n---- Installing fftw3"
   wget www.fftw.org/fftw-3.3.4.tar.gz -O${EXTERNAL}/tmp/fftw3.tar.gz
   tar -xzvf ${EXTERNAL}/tmp/fftw3.tar.gz -C${EXTERNAL} &>> $log
   cd ${EXTERNAL}/fftw-3.3.4
   ./configure --disable-alloca \
               --disable-fortran \
               --disable-static \
               --enable-shared \
               --enable-threads \
               --with-combined-threads \
               --enable-sse2 \
               --enable-avx \
               --with-our-malloc \
               --with-incoming-stack-boundary=2 \
               --prefix="${INSTALL}" &>> $log
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

# -----------------------
# googletest installation
# -----------------------

if [ -e ${INSTALL}/include/gtest/gtest.h ] && [ -e ${INSTALL}/lib/libgtest.a ] \
   && [ -e ${INSTALL}/lib/libgtest_main.a ]; then
   echo -e "\n\n---- Looks like googletest is already installed,"
   echo -e "---- are you sure you want to reinstall it?\n\n"
   select yn in "Yes" "No"; do
      case $yn in
         Yes ) INSTALL_GTEST=true; break;;
         No ) INSTALL_GTEST=false; break;;
      esac
   done
fi

if [ "${INSTALL_GTEST}" = "true" ] ; then

   echo -e "\n\n---- Installing googletest"

   git clone https://github.com/google/googletest.git ${EXTERNAL}/googletest
   cd ${EXTERNAL}/googletest/googletest
   cp -r include/* "${INSTALL}/include/"
   mkdir -p build &>> $log
   cd build && cmake .. && make &>> $log
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

# -------------------
# Python installation
# -------------------

if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/python2.7/config/libpython2.7.a ]; then
   echo -e "\n\n---- Looks like Python2.7 is already installed,"
   echo -e "---- are you sure you want to reinstall it?\n\n"
   select yn in "Yes" "No"; do
      case $yn in
         Yes ) INSTALL_PYTHON=true; break;;
         No ) INSTALL_PYTHON=false; break;;
      esac
   done
fi






if [ "${INSTALL_PYTHON}" = "true" ] ; then
   echo -e "\n\n---- Installing Python2.7\n\n"
   wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz -O${EXTERNAL}/tmp/Python-2.7.12.tgz
   tar -xzvf ${EXTERNAL}/tmp/Python-2.7.12.tgz -C"${EXTERNAL}" &>> $log
   cd ${EXTERNAL}/Python-2.7.12
   ./configure --enable-unicode=ucs4 \
               --with-threads \
               --enable-shared \
               --prefix="${INSTALL}" &>> $log
   make &>> $log
   make install &>> $log

   cd ${BLOND_HOME}

   if [ -e ${INSTALL}/include/python2.7/Python.h ] && [ -e ${INSTALL}/lib/python2.7/config/libpython2.7.a ]; then
      echo -e "---- Python has been installed successfully\n\n"
   else
      echo -e "\n\n---- Python has failed to install successfully"
      echo -e "---- You will have to manually install this library"
      echo -e "---- into directory ${BLOND_HOME}/external/install\n\n"
   fi

   #echo -e "---- Installing of Python2.7 is completed\n\n"
fi

PYTHON=${INSTALL}/bin/python2.7
export PATH="${INSTALL}/bin:$PATH"
# export PYTHONPATH="${BLOND_HOME}/python"
export LD_LIBRARY_PATH="${INSTALL}/lib:$LD_LIBRARY_PATH"
# --------------------------
# end of Python installation
# --------------------------



# ------------------------------
# Python setuptools installation
# ------------------------------
echo -e "\n\n---- Installing setuptools.."

$PYTHON -c "import setuptools" &> /dev/null
SETUPTOOLS_INSTALLED=`echo $?`
if [ "$SETUPTOOLS_INSTALLED" == "1" ]; then
    wget --no-check-certificate https://bootstrap.pypa.io/ez_setup.py -O${EXTERNAL}/tmp/ez_setup.py
    $PYTHON ${EXTERNAL}/tmp/ez_setup.py --insecure &>> $log
fi

echo -e "---- Installation of setuptools is completed\n\n"

# ------------------------------
# End of setuptools installation
# ------------------------------


# -----------------------
# Python pip installation
# -----------------------

echo -e "\n\n---- Installing pip.."

# PIP_INSTALLED=`which pip`
$PYTHON -c "import pip" &> /dev/null
PIP_INSTALLED=`echo $?`
if [ "$PIP_INSTALLED" == "1" ]; then
   wget https://pypi.python.org/packages/e7/a8/7556133689add8d1a54c0b14aeff0acb03c64707ce100ecd53934da1aa13/pip-8.1.2.tar.gz -O${EXTERNAL}/tmp/pip-8.1.2.tar.gz
   tar -xzvf ${EXTERNAL}/tmp/pip-8.1.2.tar.gz -C${EXTERNAL} &>> $log
   cd ${EXTERNAL}/pip-8.1.2
   $PYTHON setup.py install --prefix="${INSTALL}" &>> ${log}
fi

echo -e "---- Installation of pip is completed\n\n"

# -----------------------
# End of pip installation
# -----------------------


# -----------------------
# Python external modules installation
# -----------------------

#PYTHON_MODULES=( "numpy" )
PYTHON_MODULES=( "numpy" "scipy" "matplotlib" "h5py" )
$PYTHON -c "import pip" &> /dev/null
PIP_INSTALLED=`echo $?`

if [ "$PIP_INSTALLED" == "1" ]; then
   echo -e "\n\n---- PIP is needed in order to install required python modules"
   echo -e "---- If you are on Fedora/CentOS/RHEL try: yum install python-pip"
   echo -e "---- If you are on Debian/Ubuntu try: apt-get install python-pip"
   echo -e "---- and then re-run this script."
   echo -e "---- For more information, please visit this site: https://packaging.python.org/install_requirements_linux/ \n\n"
else
   for module in "${PYTHON_MODULES[@]}"; do
      echo -e "\n\n---- Installing ${module}"
      $PYTHON -c "import $module" &> /dev/null
      IS_INSTALLED=`echo $?`
      if [ "$IS_INSTALLED" == "1" ]; then
          $PYTHON -m pip install \
          --target="${INSTALL}/lib/python2.7/site-packages" \
          --global-option=build_ext \
          --global-option="-L/usr/lib" \
          --global-option="-L/usr/lib64" \
          --global-option="-L${INSTALL}/lib" \
          ${module}
      fi
      echo -e "---- Installation of ${module} is completed\n\n"
   done
fi


# ----------------------------------
# end of Python Modules installation
# ----------------------------------


echo -e "\n\n---- Exteranl depedencies installation procedure completed"
echo -e "---- Now try to build the project"
echo -e "---- You can consult the ${log} file for any errors\n\n"

rm -rf ${EXTERNAL}/tmp &> /dev/null
rm -rf ${BLOND_HOME}/setuptools-*.zip &> /dev/null
