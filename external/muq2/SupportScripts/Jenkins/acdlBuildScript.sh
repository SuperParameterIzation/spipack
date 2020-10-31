#!/bin/bash

shopt -s nocasematch

#######################################
##### Use the directory name to determine build type
#######################################
dir=`pwd`
bin_dir=$(echo "$dir" | tr '[:lower:]' '[:upper:]')

if echo $bin_dir| grep -q "PYTHON"; then
  with_python=1
  temp="${bin_dir##*PYTHON}"
  PYTHON_VERSION="${temp:0:1}"
else
  with_python=0
  PYTHON_VERSION="NA"
fi

####################################
##### SET MACHINE SPECIFIC PATHS
####################################
if [[ `hostname` == "reynolds" ]]; then

  BOOST_SOURCE=/home/mparno/util/boost_1_63_0.tar.gz
  HDF5_SOURCE=/home/mparno/util/CMake-hdf5-1.8.19.tar.gz
  NANOFLANN_SOURCE=/home/mparno/util/nanoflann-src.zip
  NLOPT_SOURCE=/home/mparno/util/nlopt-2.4.2.tar.gz
  STANMATH_SOURCE=/home/mparno/util/stanmath-v2.18.0.zip

  GTEST_DIR=/home/mparno/util/gtest_install

elif [[ `hostname` == "macys.mit.edu" ]]; then

  BOOST_SOURCE=/Users/mparno/util/boost_1_63_0.tar.gz

  # export the path so we can get cmake
  export PATH=/usr/local/bin:/usr/local/sbin:$PATH

  GTEST_DIR=/Users/jenkins/util/gtest/

fi

#######################################
##### EXTRACT COMPILER FROM WORKSPACE
#######################################

my_cc_compiler="gcc"
my_cxx_compiler="g++"

echo "C Compiler = $my_cc_compiler"
echo "CXX Compiler = $my_cxx_compiler"


#######################################
##### CHANGE INTO BUILD DIR
#######################################
BUILD_DIR="$dir/build"
INSTALL_DIR="$dir/install"
echo "BUILD_DIR = $BUILD_DIR"
echo "SOURCE = $dir"

if [ -d "$BUILD_DIR" ]; then
  echo "Build directory, $BUILD_DIR already exists."
else
  echo "Making build directory, $BUILD_DIR"
  mkdir "$BUILD_DIR"
fi

# cd into build directory and remove all previous files
cd "$BUILD_DIR"
if [ -d "CMakeFiles" ]; then
  rm CMakeCache.txt
  rm -rf CMakeFiles
  rm -rf modules
fi

#######################################
##### RUN CMAKE
#######################################
cmake \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
-DCMAKE_CXX_COMPILER=$my_cxx_compiler \
-DCMAKE_C_COMPILER=$my_cc_compiler \
-DMUQ_USE_GTEST=ON \
-DMUQ_GTEST_DIR=$GTEST_DIR \
-DBOOST_EXTERNAL_SOURCE=$BOOST_SOURCE \
-DHDF5_EXTERNAL_SOURCE=$HDF5_SOURCE \
-DNANOFLANN_EXTERNAL_SOURCE=$NANOFLANN_SOURCE \
-DNLOPT_EXTERNAL_SOURCE=$NLOPT_SOURCE \
-DSTANMATH_EXTERNAL_SOURCE=$STANMATH_SOURCE \
-DMUQ_USE_PYTHON=$with_python \
-DPYBIND11_PYTHON_VERSION=$PYTHON_VERSION \
../

#######################################
##### BUILD MUQ
#######################################
make install > OutputFromMake.txt
make examples >> OutputFromMake.txt
tail -200 OutputFromMake.txt

cd "$dir"

exit 0
