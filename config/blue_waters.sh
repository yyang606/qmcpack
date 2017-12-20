#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on Blue Waters, University of Illinois Urbnana-Champaign ##
##                                                            ##
## Last modified: Dec 20, 2017                                ##
################################################################

# Load required modules (assuming default settings have not been modified)
source $MODULESHOME/init/bash
module swap PrgEnv-cray PrgEnv-gnu

# use parallel HDF5 to speed up I/O of large jobs
if (echo $LOADEDMODULES | grep -q hdf5)
then
module unload cray-hdf5
fi
module load cray-hdf5-parallel

# FFT library is important for orbital splining
module load fftw

# miscellaneous
module load boost
module load libxml2
module load cmake
module load bwpy # numpy, h5py are libraries used in ctest

# use AMD optimized math libraries (performance critical!)
XT_FLAGS="-DHAVE_AMDLIBM=1"
# Set cmake variables, shared for cpu builds
CMAKE_FLAGS="-D QMC_INCLUDE=/u/staff/rmokos/libs/amdlibm/include \
             -D QMC_EXTRA_LIBS=/u/staff/rmokos/libs/amdlibm/lib/static/libamdlibm.a"

# always dynamic linking
export CRAYPE_LINK_TYPE=dynamic

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..

export CC=cc
export CXX=CC

# Configure and build cpu real
echo ""
echo ""
echo "building qmcpack for cpu real"
mkdir build_cpu_real
cd build_cpu_real
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      $CMAKE_FLAGS ..
make -j 32
cd ..
ln -s ./build_cpu_real/bin/qmcpack ./qmcpack_cpu_real

# Configure and build cpu complex
echo ""
echo ""
echo "building qmcpack for cpu complex"
mkdir build_cpu_comp
cd build_cpu_comp
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      -D QMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 32
cd ..
ln -s ./build_cpu_comp/bin/qmcpack ./qmcpack_cpu_comp
