# Re-test with AMD afar 24.3.0 (2026-09-29)

Follow-up to `README.md`. AMD asked us to install afar 24.3.0, re-run the seven PoCs, and check
whether the fixes for PoC 01, 02, 06 and 07 carry over to the application (GEPS/TCoGFS). They
also asked how PoC 03, 04 and 05 now affect the application.

## Environment

| | |
|---|---|
| Package | `therock-afar-24.3.0-multiarch-10.1.0-592954c.tar.bz2`, unpacked at `/opt/afar/current` |
| Compiler | AMD flang 24.0.0git (ROCm/llvm-project 3ba19712e9fb), "AMD AFAR drop #24.3.0 09/17/26" |
| Runtime | the SDK's own stack: `.info/version` 10.2.0, HIP 7.16, HSA 1.21, libomptarget 24.0git, rocFFT/rocBLAS/RCCL from the SDK |
| GPU / driver | 2x MI300X (gfx942, SPX/NPS1), amdgpu KFD 6.16.13 (unchanged) |
| Baseline | ROCm 7.2.2 + AMD flang 22.0.0git (roc-7.2.2), same machine, same day |

Standalone logs: `results_afar/` (afar) and `results/` (7.2.2). Run with
`ROCM_PATH=/opt/afar/current RESULTS=results_afar ./run_all.sh`.

Application: TCo383L72, 2 MPI ranks on 2 MI300X, six time steps. "s/step" below is the mean of
steps 3–6 (step 2 includes warm-up). The application was rebuilt entirely with afar, including all
its dependency libraries (OpenMPI, netCDF, OpenBLAS and others), because Fortran `.mod` files are
tied to the compiler version.

## Summary

| PoC | Standalone, afar 24.3.0 | Application, afar 24.3.0 |
|---|---|---|
| 01 nested reduction | **fixed**: all 5 variants correct at -O2 and -O0 | **fixed**: with the natural nested `parallel do reduction(max)` in hdiffu/rayleifr (our two-kernel workaround removed), results are **bit-identical** to the workaround build. On 7.2.2 the same change moved step-1 `sptend` by 3e-5 relative. 2.33 vs 2.35 s/step |
| 02 descriptor launch cost | **fixed**: 17 assumed-shape 710.7 → 37.1 µs/launch, 17 allocatable 173 → 39.7 µs. The per-array copies are now one batched 1456-byte transfer | whole-toolchain step time 2.38 → 2.35 s (−1%). Not separable from other changes: rocprofv3 does not see these small copies. 44% of kernel launch sites are in routines with allocatable or assumed-shape arrays (upper bound on exposure), so the pattern is no longer a significant cost for GEPS |
| 03 pool, varying sizes | **still fails**: 30.96 GB not returned after 16 iterations; threshold 0 → 1 GB alloc+free 21.3 ms | **still fails**: default memory manager and no cache of ours → VRAM 187 GB, `HSA_STATUS_ERROR_OUT_OF_RESOURCES` after step 1 (7.2.2 behaved the same). Our workaround still has to be used, **but needed rework** (see new item 1) |
| 04 hipMemcpy on OMP pointers | **still fails**: 0.051–0.056 GB/s. `HSA_ENABLE_SDMA=0` does **not** change hipMemcpy; it raises `omp_target_memcpy` D2D 58 → **1240 GB/s**, also on 7.2.2 (1232 GB/s) | `HSA_ENABLE_SDMA=0`: 7.2.2 2.38 → **2.25 s/step (−5.3%)**, afar 2.35 → **2.27 (−3.4%)**, bit-identical. We will adopt it. hipMemcpy is not used in the app since 09-18 |
| 05 RPC client | still linked, 4 threads; per-launch timing unchanged (did not reproduce standalone before either) | The afar image of GEPS **contains `__llvm_rpc_client` even with no device I/O**: the PoC 06 fix brings it in (new item 3). The RPC server thread is running. Step time equals the RPC-free 7.2.2 build (2.35 vs 2.38), whereas on 7.2.2 the RPC client cost +0.7 s/step. No visible penalty now. A clean A/B is not possible without rewriting kernels: the heap-allocated arrays are sized by runtime grid dimensions (`dxfact(lons)`, `var(latf)`, `wss(jtrun,2)`, `dqmono(levs,nvars,nxy)`, …) |
| 06 device stack, private arrays | **fixed**: nbad = 0 at every stack size including 1024. The compiler warns "using device heap allocation instead" | **fixed**: `LIBOMPTARGET_STACK_SIZE=1024` gives results **bit-identical** to 65536 (on 7.2.2 anything < 49152 gave Inf skin temperatures), same speed. 23 arrays in 4 files take the heap path |
| 07 nowait + firstprivate | **fixed**: all four variants compile and run correctly | compile test on real app code (`get_phi_gpu.f90`, two kernels with `nowait` + `taskwait`): 7.2.2 still gives the NYI error, afar compiles. Not run: our translator does not emit `nowait` yet |

Numerics, afar against 7.2.2: step-1 `sptend` 0.43807396825058964 vs 0.43807396825059497
(1.2e-14 relative), step 6 1.8e-15. Every afar run is bit-identical to every other afar run.

## New issues found while moving the application to afar

1. **libomptarget no longer links HSA.** It `dlopen`s `libhsa-runtime64.so.1` and resolves every
   entry point with `dlsym(handle, …)`. `LD_PRELOAD` interposition of `hsa_amd_memory_pool_allocate/free`
   (our PoC 03 workaround) therefore never runs. The model silently fell back to raw HSA
   alloc/free: **13.5 s/step instead of 2.4**. We rebuilt the cache as a front
   `libhsa-runtime64.so.1` that owns the soname and depends on the real runtime. With it the step time is
   back to normal, with the same 48 GB cap as on 7.2.2. (We first thought afar needed 96 GB, because
   48 GB gave alternating 2.5/4.3 s steps. That turned out to be the device-heap problem described
   below, not the cache.) A size-class pool inside
   libomptarget (PoC 03) would make this whole layer unnecessary.
2. **rocFFT in the SDK has no precompiled kernel database.** 7.2.2 ships
   `lib/rocfft/rocfft_kernel_cache.db`. Without it every FFT length is JIT-compiled one at a time at plan
   creation. Without a cache file they are compiled again every time a plan is re-created. GEPS
   needs hundreds of lengths, and start-up (normal mode initialisation) was still running after
   30 minutes when we stopped it. Workaround: `ROCFFT_RTC_CACHE_PATH=<file>`. The first run then
   compiles for about 4 minutes (a 38 MB database), and later runs start normally. Please ship the
   AOT cache for gfx942, or document this.
3. **The PoC 06 fix links the RPC client.** Runtime-sized private arrays now go to the device heap,
   and that is enough to put `__llvm_rpc_client` into the image and start the RPC server thread
   (`privarr_auto.x`: 1 RPC symbol; `privarr_fixed.x`: 0). So every program with such arrays gets
   the PoC 05 mechanism, even without device I/O. GEPS shows no step-time penalty with afar, but the
   compiler also warns the heap path "can severely degrade performance". A stack path sized from the
   known array extents would avoid both.
4. rocprofv3 from the SDK slows GEPS to 8.5 s/step (7.2.2's rocprofv3: 3.2 vs 2.4 unprofiled).
   Not investigated; it made a trace-based PoC 05 comparison unusable.

Porting notes (not AMD defects): flang 24 implements the Fortran 2023 `SPLIT` intrinsic, so GEPS's
own external `split` needed `external split` (4 routines). A clean build tree exposed a
module-order cycle in our CMake (`rank.mod`).

## 24-hour forecast (144 steps, 2026-09-30)

Same set-up as above (2x MI300X, `HSA_ENABLE_SDMA=0`); both toolchains complete the run.

| (first 24 h runs, before the fixes) | ROCm 7.2.2 (48 GB allocation cache) | afar 24.3.0 (96 GB allocation cache) |
|---|---|---|
| Peak VRAM | 171.9 / 192 GiB | **191.4 / 192 GiB** |
| Regular steps | ~2.5 s throughout | 2.6 s in hours 0–2, **~3.5 s from hour 6 on** |
| Radiation steps | ~3.6 s | 5.5–6.7 s |
| History output (6/12/18 h) | 22 s | 24.5 s |

In the 1-hour runs the two toolchains are equally fast. Over 24 hours afar was about 30% slower in
total (537 vs 415 s) and left essentially no VRAM headroom. The sections below give the cause
(device heap and VRAM oversubscription, not the cache) and the results after the fix.

### Root cause of the long-run slowdown

It is **not** gradual, and it is not our allocation cache: from one time step onward every step
costs +0.8 s. With a 96 GB cache this happens at step 34, reproducibly. The extra time sits almost
entirely in the two routines whose runtime-sized private arrays afar puts on the **device heap**
(the PoC 06 fix): `tranrs1` goes from 0.04 to 0.72 s per step, and `pbl_noah`/`sflx` from 0.09 to
0.22 s, with 1.0–2.6 s on radiation steps. On ROCm 7.2.2, which keeps these arrays on the stack,
both stay constant (0.05 s).

The onset follows free VRAM:
- `LIBOMPTARGET_HEAP_SIZE=4G` moves the slowdown from tranrs1 to sflx, at the same step.
- `16G` makes it start earlier (step ~20).
- A smaller cache (64 GB) postpones it past 6 hours (the end of that run).

**Device `malloc` gets 10–40× slower once free VRAM is low.** Confirmation: giving the 12 sflx
arrays a fixed size (`dimension(nsold)`, where `nsold = 4` is a parameter) keeps them off the heap.
pbl_noah then runs at 0.05 s, regular steps stay at 2.18–2.24 s for all 6 hours of the run with no jump, and
results are bit-identical.

We would like to ask for:
- a stack or pre-sized path for private arrays whose extents are known at kernel entry;
- device `malloc` that keeps its speed under memory pressure, or at least a diagnostic when it
  degrades.

### After moving sflx's private arrays off the heap (24 hours)

We made the sflx change permanent in the source; it is bit-identical on 7.2.2 for 24 hours. Then
we repeated the 24-hour runs:

| | Median step | Integration total | Peak VRAM |
|---|---|---|---|
| ROCm 7.2.2 | 2.28 s | 416 s | 176.9 GiB |
| afar 24.3.0, 48 GB cache | 2.35 s | 432 s | 167.2 GiB |
| afar 24.3.0, before the change (96 GB cache) | 3.07 s | 537 s | 191.4 GiB |

With the heap path avoided, afar is within 4% of 7.2.2 over 24 hours and uses less VRAM. The
remaining gap is radiation steps (~3.9 vs 3.6 s) and one 13.7 s radiation call per run (not yet
understood). The remaining afar heap users (ndslfv, gwdps) are not on the critical path. Our view
of the device-heap request above is unchanged: any application with runtime-sized private arrays
will hit this.

### The one-time 13.7 s radiation step

This happens only when VRAM reaches ~191.5 of 192 GiB (96 GB cache). It did not appear with a 48 GB
cache, and `LIBOMPTARGET_STACK_SIZE=1024` does not remove it, so scratch memory is not the cause.
KFD's per-process counters show the cause directly. In that step, `evicted_ms` rises by 19.6 s
and `page_out` by 363,320 pages: every page that XNACK migrated into VRAM at start-up
(`HSA_XNACK=1`) is evicted back to system memory. During the eviction the process's queues are
suspended. With the 48 GB cache both counters stay at 0. So it is VRAM oversubscription, which we
avoid by keeping the cache at 48 GB.

Related observation: 43 kernels in the afar image call device `malloc`, with no compile-time
warning. They are mostly `teams distribute` kernels with an inner `parallel do`, so we assume
variables shared with the inner parallel region are globalized. These kernels are the other place
where device-heap speed under memory pressure matters.
