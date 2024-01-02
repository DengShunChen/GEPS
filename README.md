# GFS TCO #

## Requirement ##
* CMake 3.21.4
* NVIDIA HPC SDK 23.11
* OpenMPI 4.0.1
* DMS v4
* NetCDF 4.9.0
* FFTW 3.5.5
* Operlib
* W3 2.0.2
* Grib2 libraries
   * g2-1.4.0
   * png-1.6.37
   * jasper-1.900.1

## Quick start ##

### Setup environment ###

```sh
MACHINE="a100" # set MACHINE="fx1000" for arm machine
. /usr/share/Modules/init/bash
module purge
module use modulefiles
module load modulefile.tcogfs.a100
module unuse modulefiles
```

### Build ###

The basic instructions (baseline, CPU only):

```sh
cmake -Bbuild -S. \
	-DCMAKE_BUILD_TYPE=Release \
	-DUSE_RSM=OFF
cd build
make -j`nproc`
```

To enable NVIDIA GPU:

```sh
cmake -Bbuild -S. \
	-DCMAKE_BUILD_TYPE=Release \
	-DUSE_RSM=OFF \
	-DUSE_CUDA=ON \
	-DUSE_ACC=ON
cd build
make -j`nproc`
```

### Run ###