# GEPS/TCoGFS on MI300X with afar 24.3.0: build and run guide

For the AMD team, in reply to the request of 2026-10-16. It covers the full application, not the
standalone PoCs in `qa/rocm_toolchain_poc/`. Following §1–§5 reproduces our runs; §6 lists how to
trigger each problem we reported.

GEPS/TCoGFS is the CWA global NWP model (Fortran + MPI + OpenMP target offload, spectral
TCo383 grid, 72 levels). The NVIDIA code (OpenACC + CUDA Fortran, `src/nvidia/`) is translated at
build time into OpenMP target code by `cmake/acc2omp.py`; AMD-specific shims and kernels are in
`src/rocm/`.

| Step | Script | Time |
|---|---|---|
| 0. Configure paths | `cp env.example.sh env.sh` and edit | — |
| 1. Install afar 24.3.0 and system packages | `01_install_afar.sh` | ~20 min (9.1 GB download) |
| 2. Build dependency libraries | `02_build_deps.sh` | ~30 min |
| 3. Build the model | `03_build_app.sh` | ~20 min |
| 4. Run 1 h / 24 h and check | `04_run.sh`, `04_run.sh 24` | ~3 min / ~10 min |

## 0. Test system

Our results come from this machine. The run needs **two MI300X with 192 GB each**: one MPI rank
per GPU. The single-rank working set does not fit on one card.

| | |
|---|---|
| GPU | 2× MI300X (gfx942), SPX / NPS1, PCIe between the two cards |
| Driver | amdgpu / KFD 6.16.13 |
| Host | AMD EPYC 9684X (192 threads), 2.2 TB RAM |
| OS | Ubuntu 24.04.5 LTS, kernel 6.8.0-138-generic, running in a container as root |
| Compiler | afar 24.3.0: `therock-afar-24.3.0-multiarch-10.1.0-592954c` (AMD flang 24.0.0git, llvm-project 3ba19712e9fb), the SDK's own HIP 7.16 / HSA 1.21 / libomptarget / rocFFT / rocBLAS / RCCL |
| Disk | SDK 23 GB unpacked, GEPS_LIB 3 GB, GEPS 2 GB, input data 10 GB, run output ~5 GB |

The scripts call `apt-get` and run the model as root (`OMPI_ALLOW_RUN_AS_ROOT=1`). On a system
without root, install the packages in `01_install_afar.sh` another way.

## 1. Source code and input data

| Item | Location | Version |
|---|---|---|
| Model | `https://github.com/DengShunChen/GEPS`, branch `feature/multi-compiler` | the commit that adds this directory (see the email) |
| Dependency libraries | `https://github.com/DengShunChen/GEPS_LIB`, branch `main` | `fe66918d557d0997ee0f317e9c4d35e9a49cebc6` |
| Input data | `geps_TCo383L72_20220815_input.tar.bz2` + `.sha256` | initial state 2022-08-15 00Z |

Both repositories are private; CWA provides access separately (see the email). Library source
tarballs (OpenMPI, netCDF, ...) are included in GEPS_LIB, so no further downloads are needed.

```bash
. ./env.sh            # after step 0
mkdir -p $WORK && cd $WORK
git clone -b feature/multi-compiler https://github.com/DengShunChen/GEPS.git
git clone https://github.com/DengShunChen/GEPS_LIB.git && git -C GEPS_LIB checkout fe66918d557d0997ee0f317e9c4d35e9a49cebc6
mkdir -p dmsdb && sha256sum -c geps_TCo383L72_20220815_input.tar.bz2.sha256 \
  && tar -I lbzip2 -xf geps_TCo383L72_20220815_input.tar.bz2 -C dmsdb
cp GEPS/qa/amd_app_repro/env.example.sh GEPS/qa/amd_app_repro/env.sh   # edit WORK if needed
```

`dmsdb/` then holds `ncep_ana_n1.ufs/TCo383l72_22081500/` (analysis and 6-h first guess) and
`bckdms.ufs/BCK_TCo383_GI30S_xnew/` (static fields). The data are in CWA's DMS format; the DMS
library and tools are built from GEPS_LIB in step 2. Static climatology files are in `GEPS/fix/`
(in the repository).

## 2. Dependency libraries

`02_build_deps.sh` builds everything into `GEPS_LIB/install/rocm/` with the afar compilers. All of
it must be rebuilt with afar, because Fortran `.mod` files are tied to the compiler version.

| Library | Version | Notes |
|---|---|---|
| OpenMPI | 4.1.6 | `src/openmpi/compile.sh`, `--enable-mpi-fortran=all`, wrappers point at amdflang/amdclang |
| netCDF-C / netCDF-Fortran | 4.9.0 / 4.6.0 | with HDF5 1.13.2, szip 2.1.1 |
| OpenBLAS | 0.3.27 | host BLAS/LAPACK |
| FFTW | 3.3.5 | host FFT |
| zlib, libpng, jasper | 1.2.11, 1.6.37, bundled | GRIB2 packing |
| g2lib, w3lib | 1.4.0, 2.0.2 | GRIB |
| DMS (dms-4.0.0), dms38key, operlib | CWA | data I/O |
| MCT, geps_coupler | bundled | coupler (built but unused in this case) |
| CMake | 3.21.4 | bootstrapped from source first |

Known build traps, already handled by the script: `gtar` is called by the libpng build (a symlink
to `tar` is put on `PATH`); operlib needs `ksh`; netCDF needs `m4`.

Check: every package prints `[OK]` and `install/rocm/` contains `openmpi-4.1.6`, `netcdf-*`,
`dms-4.0.0`, among others.

## 3. Building the model

`03_build_app.sh` runs `GEPS/qa/afar_app/build_afar_hip.sh`, which:

1. loads the afar toolchain and GEPS_LIB modules (`GEPS_LIB/env.sh`, `module load geps/1.0`);
2. configures CMake:
   ```
   cmake -Bbuild_afar_hip -S. -DCMAKE_BUILD_TYPE=Release -DGEPS_COMPILER=rocm -DUSE_HIP=ON -DUSE_OMP=ON \
     -DUSE_CUDA=OFF -DUSE_ACC=OFF -DUSE_PAR=OFF -DUSE_NDMS=OFF -DUSE_RSM=OFF -DCWBSUM=OFF -DTIMCOMCPL=OFF \
     -DUSE_aeroclx=OFF -DAMD_GPU_ARCHS=gfx942 \
     -DCMAKE_Fortran_COMPILER=$ROCM_PATH/bin/amdflang -DCMAKE_C_COMPILER=$ROCM_PATH/bin/amdclang \
     -DCMAKE_CXX_COMPILER=$ROCM_PATH/bin/amdclang++
   ```
   Compiler flags are in `cmake/compilers/ROCM.cmake` (`-O3 -fopenmp --offload-arch=gfx942`, hipBLAS,
   rocFFT, RCCL);
3. pre-builds `rank.mod` (a module-order cycle in our CMake on a clean tree);
4. builds `tcogfs.x`. During this, `acc2omp.py` translates every file in `src/nvidia/` into
   `build_afar_hip/src/acc2omp/`; these generated files are what amdflang compiles;
5. builds the HSA allocation cache `src/rocm/geps_hsa_pool.cc` (§4, item 3) into
   `build_afar_hip/lib/libgeps_hsa_pool.so` and `build_afar_hip/lib/hsa_front/libhsa-runtime64.so.1`.

Output: `build_afar_hip/bin/tcogfs.x`. The ROCm 7.2.2 build is the same with
`ROCM_PATH=/opt/rocm` and a GEPS_LIB built with 7.2.2.

## 4. Running

`04_run.sh [hours]` calls `GEPS/job/TCo383L72_IC_sample_rocm_2gpu` (2 MPI ranks, `NPEX=1 NPEY=2`),
which checks the inputs, links the DMS data into a working database, and sources
`job/regression.ksh`. That file sets the runtime environment below and launches
`mpiexec -n 2 build_afar_hip/bin/tcogfs.x`. Each setting can be overridden from the command line.

| Variable | Value | Why |
|---|---|---|
| `HIP_VISIBLE_DEVICES`, `ROCR_VISIBLE_DEVICES` | `0,1` (`GEPS_ROCM_DEVICES`) | one rank per GPU, rank *r* uses device *r* |
| `OMP_TARGET_OFFLOAD` | `MANDATORY` | |
| `OMP_NUM_THREADS` | `1` | |
| `HSA_XNACK` | `1` | |
| `HSA_ENABLE_SDMA` | `0` | `omp_target_memcpy` D2D 58 → 1240 GB/s, −3% step time (PoC 04) |
| `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD` | `0` | libomptarget's pool reuses only exact sizes and grows until OUT_OF_RESOURCES (PoC 03) |
| `LIBOMPTARGET_STACK_SIZE` | `65536` | needed on 7.2.2 (PoC 06); afar is correct at 1024 too |
| `LD_LIBRARY_PATH` | `build_afar_hip/lib/hsa_front` first | **HSA allocation cache.** With the threshold at 0 every `data create` is a raw HSA alloc/free (3.6k pairs per step). The cache keeps freed blocks ≥ 256 KB by size class. afar's libomptarget `dlopen`s HSA, so the cache is a `libhsa-runtime64.so.1` that sits in front of the real one (on 7.2.2 it is an `LD_PRELOAD`). `GEPS_HSA_POOL=0` disables it |
| `GEPS_HSA_POOL_CAP_GB` | `48` | cache limit per rank. 96 drives VRAM to 191.5/192 GiB and triggers KFD eviction |
| `ROCFFT_RTC_CACHE_PATH` | `GEPS/qa/afar_app/rocfft_cache/afar_gfx942.db` | the SDK's rocFFT has no AOT kernel database. The first run JIT-compiles ~4 min and fills a 38 MB file; later runs start normally |
| `GEPS_FFT_PLAN_SETS` | `8` | rocFFT plan cache in our FFT wrapper |
| `GEPS_TAUE` | hours | forecast length (default 1) |
| `GEPS_PROF_WRAP` | e.g. `rocprofv3 --kernel-trace --stats -d DIR -o r_%pid% --` | prefix for the binary; also set `OMP_TOOL=disabled` |

Run time on our system with afar: 1 h forecast about 2 min 40 s wall time in total (start-up and
normal-mode initialisation about 2.3 min, then 6 steps in about 15 s); 24 h forecast about 10 min.

## 5. Test case and verification

| | 1 h | 24 h |
|---|---|---|
| Case | TCo383L72 (JCAP 383, 72 levels), initial state 2022-08-15 00Z | same |
| Steps | 6 (10-min step) | 144 |
| Radiation | steps 1 (and every 6 steps) | every 6 steps |
| History output | — | every 6 h (the following step takes ~22 s) |
| Command | `./04_run.sh` | `./04_run.sh 24` |

A successful run ends with `PROGRAM CWBGFS HAS ENDED`. Correctness is checked with
`surf pres tend rms(GPU)`, printed once per step. `check_run.py <log>` prints this, the step times,
and the differences from `reference.json`:

| Step | afar 24.3.0 | ROCm 7.2.2 | CPU (32 cores) |
|---|---|---|---|
| 1 | 0.43807396825058964 | 0.43807396825059497 | 0.4380739682505907 |
| 6 | 0.3925504241467753 | 0.39255042414676805 | |
| 144 (24 h) | 0.35008264843432285 | 0.35047487674616895 | |

- The same binary on the same machine gives bit-identical values run to run.
- afar and 7.2.2 differ by 1.2e-14 at step 1. Over 144 steps the difference grows to 1e-3, as
  expected for a chaotic model, so the 24 h value is only comparable within one toolchain.
- We have only run on this system. A step-1 difference much larger than 1e-14 relative to these
  values would point to a problem rather than rounding.

Reference timings (2× MI300X, our system):

| | 1 h: mean of steps 3–6 | 24 h: sum of steps 3–144 | 24 h median step | 24 h peak VRAM |
|---|---|---|---|---|
| afar 24.3.0 | 2.31 s | 431.8 s | 2.35 s | 167.2 GiB |
| ROCm 7.2.2 | 2.21 s | 415.7 s | 2.28 s | 176.9 GiB |

Step 2 includes warm-up (3.8–4.9 s). Radiation steps take about 3.6 s (7.2.2) and 3.9 s (afar).

## 6. Reproducing the reported problems

Each item changes one setting relative to the default run in §4.

| # | Problem | How to trigger | What to look for |
|---|---|---|---|
| 1 | libomptarget pool (PoC 03) | `GEPS_HSA_POOL=0 LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=default ./04_run.sh` | VRAM ~187 GB after step 1, then `HSA_STATUS_ERROR_OUT_OF_RESOURCES` |
| 2 | Raw HSA alloc/free cost | `GEPS_HSA_POOL=0 ./04_run.sh` | ~12.6 s/step instead of ~2.3 (measured on 7.2.2; afar 13.5) |
| 3 | HSA loaded through `dlopen` | `nm -D $ROCM_PATH/lib/libomptarget.so* \| grep hsa_` is empty on afar. To see the effect, rename `build_afar_hip/lib/hsa_front`; `regression.ksh` then falls back to `LD_PRELOAD` of `libgeps_hsa_pool.so`, as on 7.2.2 | the preload is never called: 13.5 s/step, as in item 2 |
| 4 | rocFFT without AOT kernels | delete `qa/afar_app/rocfft_cache/afar_gfx942.db` and run | ~4 min extra at start-up. With no cache file at all, start-up does not finish in 30 min |
| 5 | Device malloc slows under VRAM pressure | revert our fix: `git -C $GEPS_SRC show abac239f -- src/nvidia/sflx_gpu.f90 \| git -C $GEPS_SRC apply -R`, rebuild, `GEPS_HSA_POOL_CAP_GB=96 ./04_run.sh 24` | from step ~34 every step costs +0.8 s (24 h total ~537 s instead of 432). Per-routine timers: build with `GEPS_ACC2OMP_TIMERS=intgrt_gpu,diabat_gpu ./03_build_app.sh build_afar_timers` and run with `GEPS_ROCM_BUILD_DIR=build_afar_timers`; `pbl_noah_gpu` and `tranrs1_gpu` grow |
| 6 | KFD eviction near the VRAM limit | `GEPS_HSA_POOL_CAP_GB=96 ./04_run.sh 24` (with the sflx fix) | one radiation step of 11–14 s; `qa/tools/kfdmon.sh LOG OUT` shows `evicted_ms` and `page_out` rising in that step (counters under `/sys/class/kfd/kfd/proc/<host pid>/`; set the two GPU ids in the script from `/sys/class/kfd/kfd/topology/nodes/*/gpu_id`) |
| 7 | Nested reduction (PoC 01, fixed in afar) | `GEPS_ACC2OMP_VEC_REDMAX=1 ./03_build_app.sh build_afar_redmax`, run with `GEPS_ROCM_BUILD_DIR=build_afar_redmax` | afar: bit-identical to the default build. 7.2.2: step-1 value off by 3e-5 |
| 8 | Private arrays on the device stack (PoC 06, fixed in afar) | `LIBOMPTARGET_STACK_SIZE=1024 ./04_run.sh` | afar: bit-identical. 7.2.2: Inf skin temperatures / wrong values below 49152 |
| 9 | RPC client in the image (PoC 05) | `$ROCM_PATH/llvm/bin/llvm-objcopy --dump-section=.llvm.offloading=img.bin $GEPS_SRC/build_afar_hip/bin/tcogfs.x`, then look for `__llvm_rpc_client` in the device ELF | present on afar even without device I/O (heap-allocated private arrays). Absent on 7.2.2 |
| 10 | rocprofv3 overhead | `OMP_TOOL=disabled GEPS_PROF_WRAP="rocprofv3 --kernel-trace --stats -d /tmp/prof -o r_%pid% --" ./04_run.sh` | afar's rocprofv3: 8.5 s/step; 7.2.2's: 3.2 s/step |

A run that stops early leaves `HSA_STATUS_ERROR_*` or `memory access fault` in the log;
`check_run.py` prints these.

## 7. Files in this directory

| File | Purpose |
|---|---|
| `env.example.sh` | paths; copy to `env.sh` |
| `01_install_afar.sh` … `04_run.sh` | steps 1–4 |
| `check_run.py`, `reference.json` | verification against our results |
| `pack_data.sh` | CWA side: builds the input-data tarball |
| `REPLY_zh.md` | CWA side: reply email draft |

Background: `ROCM_OPTIMIZATION_REPORT.md` (how the port reached 2.2 s/step),
`qa/rocm_toolchain_poc/AFAR_24.3.0_RETEST.md` (afar re-test report),
`qa/afar_app/INSTALL_AFAR.md` (our own afar installation log).
