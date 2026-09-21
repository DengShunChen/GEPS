
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

Artifacts land in `build_<compiler>/` (`rocm-cpu` → `build_rocm_cpu`). The ROCm job scripts
currently run `build_rocm_hip/bin/tcogfs.x` (see the ROCm section below); aligning `./compile rocm`
with that tree is an open item.

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

## AMD ROCm / Instinct MI300X ##

Status (2026-09-19): the full TCo383L72 tau=1 case runs on 2× MI300X (gfx942, ROCm 7.2.2,
amdflang 22.0.0) and all six `surf pres tend rms` values match the CPU run to 1e-11..1e-14.
2.43 s per integration step on 2 GPUs versus 10.5 s on 32 CPU cores of the same machine
(NVIDIA reference: 0.34 s). Documents:

* [`ROCM_PORT_HANDOFF.md`](ROCM_PORT_HANDOFF.md) — layer-by-layer diagnosis log, environment facts,
  working rules and the open-items list (§9). Read §0, §4 and §8 first.
* [`ROCM_OPTIMIZATION_REPORT.md`](ROCM_OPTIMIZATION_REPORT.md) — the 740 s → 2.43 s/step work,
  one measured root cause per change.
* [`TICKET_amd_rocm_support.md`](TICKET_amd_rocm_support.md), [`TICKET_internal_infra.md`](TICKET_internal_infra.md)
  — outward-facing issue reports; `qa/rocm_toolchain_poc/` holds stand-alone reproducers with results.

### How the port works

`src/` and `src/nvidia/` are not modified. At build time `cmake/acc2omp.py` translates the
OpenACC / CUDA-Fortran GPU sources to OpenMP target offload, and `src/rocm/` provides the
AMD-side shims and replacement kernels (`hip_compat.cc` for cuBLAS/cuFFT/NCCL → hipBLAS/hipFFT/RCCL,
`geps_acc_wait.f90` for stream waits, tiled NDSL advection kernels, `geps_hsa_pool.cc`).
Fixes that nvfortran tolerates but OpenMP does not (uninitialised `private` scalars, missing
waits after asynchronous collectives, graph-capture regions) live in the translator, so the
NVIDIA path is unaffected.

### Build

```sh
export GEPS_LIB_ROOT=/path/to/GEPS_LIB GEPS_COMPILER=rocm
source $GEPS_LIB_ROOT/env.sh && load_geps_compiler
export COMP_MP=mpifort C_COMP_MP=mpicc CXX_COMP_MP=mpicxx   # required whenever cmake reconfigures
cmake -Bbuild_rocm_hip -S. -DCMAKE_BUILD_TYPE=Release -DUSE_RSM=OFF -DCWBSUM=OFF -DTIMCOMCPL=OFF \
      -DUSE_aeroclx=OFF -DGEPS_COMPILER=rocm -DUSE_CUDA=OFF -DUSE_HIP=ON -DUSE_ACC=OFF -DUSE_PAR=OFF \
      -DUSE_OMP=ON -DUSE_NDMS=OFF -DAMD_GPU_ARCHS=gfx942
cmake --build build_rocm_hip --target tcogfs.x -j 32
src/rocm/build_hsa_pool.sh          # libgeps_hsa_pool.so (LD_PRELOAD allocation cache), kept out of CMake
```

`-DUSE_HIP=OFF` with the same toolchain gives the host-only reference build (`build_rocm_cpu`).

### Run

```sh
cd job
./TCo383L72_IC_sample_rocm_2gpu 2>&1 | tee TCo383L72_IC_sample_rocm.log.$(date +%Y%m%d)_2gpu   # 2 MPI ranks, one per GPU
./TCo383L72_IC_sample_rocm_cpu                                                                  # CPU reference (32 ranks)
./TCo383L72_IC_sample_rocm_rocprof                                                              # rocprofv3 kernel trace
```

A successful run ends with `tau= 1.000` and `PROGRAM CWBGFS HAS ENDED`; compare the six
`surf pres tend rms` lines against the CPU log. `job/cmp_ckpt.py` / `job/cmp_mf.py` compare
optional checkpoints. A 2-GPU tau=1 run takes ~6.5 min wall (model start to end, dominated by
initialisation and NNMI; 2.4 s per step); the CPU smoke ~4 min. Do not run both at once — they share
the DMS key.

### Runtime knobs (set by `job/regression.ksh`, all overridable)

| Variable | Default | Meaning |
|---|---|---|
| `GEPS_ROCM_DEVICES` | `0` | `HIP_VISIBLE_DEVICES` / `ROCR_VISIBLE_DEVICES` for the rank (the 2-GPU script binds one card per rank) |
| `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD` | `0` | libomptarget's pool only reuses identical sizes and grows to OOM on the radiation blocks; keep at 0 |
| `LIBOMPTARGET_STACK_SIZE` | `65536` | runtime-sized `private` arrays corrupt across lanes below 49152 |
| `GEPS_HSA_POOL` | `1` | preload `libgeps_hsa_pool.so`; `0` disables (costs ~4 s/step in HSA alloc/free) |
| `GEPS_FFT_PLAN_SETS` | `8` | rocFFT plan-set LRU budget; the model needs 7 |
| `GEPS_FFT_THREADS` | `4` | host threads issuing per-latitude FFTs, one stream each |
| `GEPS_PROF_WRAP` | — | command prefix for `tcogfs.x`, e.g. `rocprofv3 --kernel-trace` (needs `OMP_TOOL=disabled`) |

Translator knobs are read by `cmake/acc2omp.py` at build time: `GEPS_ACC2OMP_SPMD` (default `1`,
fold physics longitude loops into `collapse(3)`), `GEPS_ACC2OMP_PROBES` (default `0`, debug probes),
`GEPS_ACC2OMP_TIMERS=<files>` (per-routine `system_clock` timers).

### Known open items

Single-GPU (NPEY=1) run not yet exercised at this resolution; launch latency and kernel
efficiency remain the gap to the NVIDIA reference; upstream reports for the compiler/runtime
issues in the tickets. Full list in `ROCM_PORT_HANDOFF.md` §9.

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
