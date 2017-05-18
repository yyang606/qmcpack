#!/bin/bash
# load intel17.0.0, qmcpack and espress
module load env/dev
export CC=mpiicc
export CXX=mpiicpc

# build complex executable
if [ ! -d complex_build ]; then
  mkdir complex_build
  cd complex_build
  cmake -DCMAKE_BUILD_TYPE=release -DQMC_COMPLEX=1 ..
  make -j24
  cd ..
else
  cd complex_build
  make -j4 qmcpack
  cd ..
fi

# build real executable
if [ ! -d build ]; then
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=release -DQMC_COMPLEX=0 ..
  make -j24
  cd ..
else
  cd build
  make -j4 qmcpack
  cd ..
fi

# build debug executable
if [ ! -d debug_build ]; then
  mkdir debug_build
  cd debug_build
  cmake -DCMAKE_BUILD_TYPE=debug -DCMAKE_C_FLAGS=-O0 -DCMAKE_CXX_FLAGS=-O0 ..
  make -j24
  cd ..
else
  cd debug_build
  make -j4 qmcpack
  cd ..
fi

# make symbolic links
cd build/bin
ln -s ../../complex_build/bin/qmcpack qmcpack_comp
ln -s ../../debug_build/bin/qmcpack debug_qmcpack
cd ..
