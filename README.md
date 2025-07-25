
[![pipeline status](http://git.rdc.cwb/tco/tco639l72/badges/development/pipeline.svg)](http://git.rdc.cwb/tco/tco639l72/-/commits/development)

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

## Quick start (x86_64 / GPU) ##

### Setup environment ###

```sh
MACHINE="a100"
. /usr/share/Modules/init/bash
module purge
module use modulefiles
module load modulefile.tcogfs.a100
module unuse modulefiles
```

### Build ###

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

```sh
cd job
pjsub TCo383L72_IC_sample_a100
```

## Quick start (ARM) ##

### Setup environment and build ###

```sh
./build.sh fx1000
```

### Run ###

```sh
cd job
pjsub TCo383L72_IC_sample_fx1000
```

### CI ###

Please request at least two GPUs for CI. Some unit tests will validate the correctness of CUDA-aware MPI.
The `pjsub` command in the following block is an example to request whole resources of a GPU node (32 core + 8 GPU). You can customize the options according to your needs.

```sh
pjsub -L "vnode=1,vnode-core=32,ru=rscunit_pg01,rg=gpu-rd-large,gpu=8" --sparam wait-time=100 -g sum --interact
bash qa/L0_test_subroutine/test.sh
```

### Guidance for GPU porting ###

```sh
https://www.notion.so/OpenACC-porting-TCo-GPU-1826f2851fb180e894aec45be810557a?pvs=4
```
