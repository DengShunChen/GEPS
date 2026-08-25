
[![pipeline status](http://git.rdc.cwb/tco/tco639l72/badges/development/pipeline.svg)](http://git.rdc.cwb/tco/tco639l72/-/commits/development)

# GFS TCO #

## Requirement ##
* [GEPS_LIB](../GEPS_LIB) dependency stack (per-compiler install tree)
* CMake 3.21.4
* Toolchains (via `GEPS_COMPILER`):
  * `nvidia` — NVIDIA HPC SDK + CUDA/OpenACC
  * `gnu` — GCC/GFortran + OpenMPI (CPU)
  * `intel` — Intel oneAPI (icx/ifx + MPI) (CPU)
  * `fujitsu` — Fujitsu tcsds (`fccpx`/`frtpx`; alias `tcsds`) (CPU)
  * `rocm` — AMD ROCm (`amdclang`/`amdflang` + OpenMP offload / HIP; alias `amd`)

## Quick start (multi-compiler / GEPS_LIB) ##

Build matching GEPS_LIB first (`./build_all.sh --compiler <name>`), then:

```sh
export GEPS_LIB_ROOT=/path/to/GEPS_LIB   # default: ../GEPS_LIB

./compile nvidia    # NVIDIA GPU (CUDA/OpenACC)
./compile gnu       # CPU
./compile intel     # CPU
./compile fujitsu   # or: ./compile tcsds
./compile rocm      # AMD Instinct MI300X (HIP + OpenMP offload, gfx942)
./compile rocm-cpu  # AMD toolchain, host OpenMP only (no GPU offload)
```

Artifacts land in `build_<compiler>/` (`rocm-cpu` → `build_rocm_cpu`).

`./compile rocm` does **not** use NVIDIA OpenACC: amdflang cannot offload OpenACC to AMDGPU. GPU kernels are compile-time translated (`cmake/acc2omp.py`) to OpenMP target + HIP libraries. Original NVIDIA sources are unchanged.

### Math libraries (BLAS/LAPACK)

| `GEPS_COMPILER` | Math library | Device |
|-----------------|--------------|--------|
| `gnu` | OpenBLAS (`GEPS_LIB` `openblas/*`, else system) | CPU |
| `rocm` | OpenBLAS + hipBLAS/hipFFT/hipSOLVER/RCCL | AMD GPU (`AMD_GPU_ARCHS`, default `gfx942`) |
| `intel` | Intel MKL (`module load mkl/...`, `-qmkl=sequential`) | CPU |
| `nvidia` | NVHPC BLAS/LAPACK + cuBLAS/cuFFT/NCCL | NVIDIA GPU |
| `fujitsu` | Fujitsu SSL2 (`-SSL2BLAMP`) | CPU |

Do **not** use OpenBLAS for intel/nvidia/fujitsu.

## Quick start (legacy x86_64 / GPU modulefile) ##

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

## Quick start (legacy ARM / Makefile) ##

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
