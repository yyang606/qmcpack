#!/bin/bash

module load modules/2.1-20230222
module load gcc/10.4.0 
module load openmpi/4.0.7
module load intel-oneapi-mkl/2023.0.0
module load hdf5/mpi-1.10.9
module load boost/1.80.0

module load git cmake

export OMP_NUM_THREADS=1

export PATH=/mnt/home/yyang/software/qmcpack/build/bin:$PATH

export CC=mpicc
export CXX=mpicxx
cmake -D QMC_COMPLEX=1 ..
echo Configuring and building froom inside the build directory.
echo Check the results of the CMake configuration to ensure that the preferred
echo compilers and libraries have been selected. See README and documentation 
echo for guidance.
#cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
make -j 128
#cd build; cmake ..; make -j 128; cd ..
