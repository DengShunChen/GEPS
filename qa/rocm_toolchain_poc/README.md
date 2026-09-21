# ROCm 7.2.2 / amdflang toolchain issues — self-contained reproducers

Seven independent issues met while porting a Fortran + OpenMP-offload
numerical weather prediction model to MI300X. Every directory below is a
standalone program (no application code, no MPI, no BLAS) with a `run.sh`
that builds it with a stock ROCm install and prints a verdict or a
measurement. `run_all.sh` runs them all (~6 min) and keeps the logs in
`results/`; `clean.sh` removes build products.

```bash
./run_all.sh          # everything, logs in results/
./01_nested_reduction/run.sh   # or any single one
```

## Environment the numbers below were taken on

| | |
|---|---|
| ROCm | 7.2.2 (`/opt/rocm` -> `/opt/rocm-7.2.2`) |
| Compiler | AMD flang 22.0.0git (roc-7.2.2 26084 f58b06dce1f9c15707c5f808fd002e18c2accf7e); `amdclang++` from the same package for the one C++ case |
| OpenMP runtime | `/opt/rocm-7.2.2/lib/llvm/lib/libomptarget.so.22.0git` |
| GPU | AMD Instinct MI300X, gfx942, SPX / NPS1, 192 GiB |
| Kernel | 6.8.0-138-generic, amdgpu (KFD) 6.16.13 |
| Flags | `-O2 -fopenmp --offload-arch=gfx942`, `OMP_TARGET_OFFLOAD=MANDATORY` |

## Summary

| # | Directory | Component | Kind | One line |
|---|---|---|---|---|
| 1 | `01_nested_reduction` | amdflang codegen | **wrong results, silent** | `parallel do reduction(+/max:x)` nested under `target teams distribute` returns 0 for every team; at `-O0` even the flat combined form returns 0 |
| 2 | `02_descriptor_launch_cost` | amdflang + libomptarget | performance | every launch re-copies the descriptor of every allocatable / assumed-shape array it touches: 17 assumed-shape dummies cost 710 µs per launch (68 tiny synchronous copies) vs 22 µs for explicit-shape |
| 3 | `03_pool_varying_size` | libomptarget memory manager | OOM / performance | the pool only reuses blocks of identical size, so varying-size maps grow VRAM monotonically (31 GB held after 16 iterations); with the pool disabled a 1 GB alloc/free pair costs ~21 ms |
| 4 | `04_hipmemcpy_omp_ptr` | HIP runtime (CLR) | performance, silent | `hipMemcpy*` between device buffers allocated by libomptarget runs at **0.05 GB/s** (CPU large-BAR path) vs 1300 GB/s for `hipMalloc` buffers; data is correct so nothing warns |
| 5 | `05_rpc_client_latency` | flang runtime + libomptarget | performance (app-level) | a `print` on a never-taken branch links `__llvm_rpc_client` into the image and starts an RPC server thread; in the application this delayed kernel-completion observation by ~one kernel duration (see caveat) |
| 6 | `06_device_stack_private_arrays` | libomptarget / device codegen | **wrong results, silent** | 8 runtime-sized private arrays of 4 doubles per lane (256 B) are overwritten by other lanes with the default `LIBOMPTARGET_STACK_SIZE`; fixed-size arrays of the same size are fine |
| 7 | `07_nowait_firstprivate_nyi` | amdflang | compile error | `target ... nowait` fails with "not yet implemented: Unhandled clause privatization for deferred target tasks" as soon as the region reads any scalar not explicitly `map`ped |

Issues 1 and 6 produce wrong numbers with no diagnostic; in the application
they showed up as a 3e-5 drift in a global diagnostic and as `Inf` surface
temperatures on a varying ~10 % of land points respectively, and each cost
about a day to trace back to the toolchain.

---

## 1. Nested `parallel do reduction` under `target teams distribute` returns 0

`01_nested_reduction/nested_reduction.f90` computes the same per-level sum
five ways. Host OpenMP gets all five right; on the GPU:

```
(a) nested  parallel do reduction(+)   WRONG got(1)=  0.000000E+00 expected   1.001000E+03
(b) nested  parallel do reduction(max) WRONG got(1)=  0.000000E+00 expected   1.001000E+00
(c) teams/distribute split, nested (+) WRONG got(1)=  0.000000E+00 expected   1.001000E+03
(d) flat teams distribute parallel do  OK    got(1)=  1.001000E+03
(e) inner loop serial, no reduction    OK    got(1)=  1.001000E+03
```

Same at `-O1`/`-O3`. At **`-O0` variant (d) — the plain combined
`target teams distribute parallel do reduction(+:s)` — is wrong as well.**
The pattern that fails is the textbook one:

```fortran
!$omp target teams distribute private(wt)
do k = 1, lev
   wt = 0.0d0
   !$omp parallel do reduction(+:wt)
   do i = 1, n
      wt = wt + a(i, k)
   end do
   wsum(k) = wt            ! <- 0 for every k
end do
```

Workaround we use: never nest a reduction under `teams distribute`; write
team-level reductions as a flat kernel or as two kernels.

## 2. Per-launch descriptor re-mapping for allocatable / assumed-shape arrays

`02_descriptor_launch_cost/desc_launch.f90` launches the same trivial
64-iteration kernel through four interfaces, all data already resident
(`enter data` once, 2000 launches):

```
   2 explicit-shape dummies                11.1 us/launch
  17 explicit-shape dummies                21.9
  17 assumed-shape dummies (descriptor)   710.7
  17 allocatables in scope                173.0
  17 allocatables + map(present,alloc:)   170.6
```

`LIBOMPTARGET_INFO=-1` shows why. Per launch:

```
   8 arguments (2 explicit):     0 host->device copies,              0 device->host
  62 arguments (17 explicit):    0 host->device copies,              0 device->host
  87 arguments (17 allocatable): 17 host->device copies (1280 B),    0 device->host
  88 arguments (17 assumed):     51 host->device copies (1552 B),   17 device->host
```

Every allocatable costs one synchronous 72–88 byte descriptor copy per
launch (~10 µs); every assumed-shape dummy costs three copies up plus one
8-byte copy back (~40 µs). `-fno-defer-desc-map`, `map(present,alloc:)`
and `has_device_addr` do not change it. The descriptor of an array that is
already present and unchanged should not need to be re-copied on every
launch. A model with ~40k launches per time step spends seconds here.

## 3. Memory manager: no reuse across sizes, and the cost of turning it off

`03_pool_varying_size/pool_growth.f90` maps, uses and releases one array per
iteration, each a slightly different size (1.5 → 2.3 GB). Free VRAM (via
`hipMemGetInfo`) after each `exit data map(release:)`:

```
default pool:                       iter 16: not returned 31.45 GB   RESULT: FAIL
LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0: not returned  0.50 GB   RESULT: PASS
same size every iteration:          not returned  2.37 GB (one block) RESULT: PASS
```

Released blocks are kept but only ever handed out for an identical size,
so a program whose work arrays follow the current block size (the normal
Fortran pattern) exhausts VRAM. The only fix is to disable the pool — and
`alloc_cost.f90` shows what that costs: every `enter/exit data` pair becomes
a raw `hsa_amd_memory_pool_allocate/_free`:

```
   size        enter data (alloc)   exit data (free)     [us]
    4 MB             150.9             168.5     (pool: 28.8 / 1.6)
   64 MB             106.5             164.2
 1024 MB           21010.2             696.4     (~21 ms per alloc+free pair)
```

In the application (3.6k create/delete pairs per step) that was 4.1 s of an
11.8 s step; we ended up writing our own `LD_PRELOAD` size-class cache
around the HSA pool calls. A size-class (best-fit / bucketed) pool inside
libomptarget would remove both problems.

## 4. `hipMemcpy*` on libomptarget-allocated device memory runs on the CPU

`04_hipmemcpy_omp_ptr/memcpy_omp_ptr.cc` (C++, `amdclang++ -fopenmp`) copies
64 MB device-to-device between `hipMalloc` buffers and between
`omp_target_alloc` buffers:

```
  copy API                     hipMalloc bufs       omp bufs   [GB/s]
  hipMemcpy(Default)                   1325.6          0.051   <-- CPU path
  hipMemcpy(DeviceToDevice)            1347.1          0.051   <-- CPU path
  hipMemcpyDtoD                        1336.6          0.051   <-- CPU path
  hipMemcpyAsync(D2D)+sync             1217.4          0.051   <-- CPU path
  omp_target_memcpy                      54.1         57.850
  hipMemcpy on omp buffers: 0 wrong bytes of 67108864 (data is correct, only slow)
```

`hipPointerGetAttributes` reports both pointers as `type=2` (device) on the
same device, yet every HIP copy API falls back to a host memcpy through the
large-BAR mapping, 25 000× slower, with correct results and no warning.
`HSA_XNACK=1` makes no difference. This also bites every HIP library handed
an OpenMP-mapped buffer: we found it inside RCCL's out-of-place
`ncclAllGather` (local-chunk copy) and in our own D2D state copies; it was
the single largest cost in the port twice (a 343 MB array took 13 s).

Secondary: `omp_target_memcpy` D2D on the same buffers tops out at ~57 GB/s,
25× below `hipMemcpy` on `hipMalloc` memory.

## 5. Device-side Fortran I/O links the RPC client and starts a server thread

`05_rpc_client_latency/rpc_latency.F90` built with `-DDEVIO` differs from
the plain build only by one `print` on a branch that can never be taken.
`check_image.sh` extracts the embedded gfx942 ELF and lists its symbols:

```
rpc_plain.x: 0 rpc/Fortran-I/O symbols in device image
rpc_devio.x: 1 rpc/Fortran-I/O symbols in device image
    __llvm_rpc_client
rpc_plain: 3 threads      rpc_devio: 4 threads     (the extra one is libomptarget's RPC server)
```

**Caveat — the timing effect does not reproduce in this microbenchmark.**
Per-launch wall time for 0.25–8 ms kernels is the same with and without the
RPC client here (synchronous and `nowait` variants both included). In the
full application (2 MPI ranks, many streams, ~40k kernels/step) the RPC
build showed, in `rocprofv3 --kernel-trace`, the GPU idle for 2–6 ms after
every kernel longer than ~1 ms with the host in `hsa_signal_wait_scacquire`
past the kernel's end timestamp; removing the three device-side
`write`/`stop` statements (image no longer contains `__llvm_rpc_client`)
took the step from 3.49 s to 2.80 s with bit-identical results. We include
this PoC for the mechanism (any unreachable I/O statement changes the
runtime's kernel-wait path for the whole program) and for the image check,
and would welcome guidance on what the RPC-enabled wait path does
differently.

## 6. Runtime-sized private arrays are corrupted with the default device stack

`06_device_stack_private_arrays/privarr.F90`: each lane owns 8 private
automatic arrays `a1(nsoil)..a8(nsoil)` (nsoil = 4, i.e. 256 bytes per
lane), fills them with a lane-unique pattern, does some work with a
data-dependent trip count, then checks its own arrays:

```
automatic, default stack:  nbad = 16938192 of 19070976   (varies run to run)
automatic, STACK_SIZE=1024: nbad = 16916064
automatic, STACK_SIZE=4096: nbad = 14397899
automatic, STACK_SIZE=16384: nbad = 0
-DFIXED a1..a8(4), default: nbad = 0
```

89 % of the private elements were overwritten by other lanes. The same
data in one 2-D automatic array `a(nsoil,20)` is fine, so this is per-array,
not per-byte, and 256 B/lane is far below any plausible default stack. There
is no fault, no message — just wrong numbers on a varying subset of lanes.
Either the default should be large enough for this, or an overflow should
be detected. (In the application the fix was `LIBOMPTARGET_STACK_SIZE=65536`,
which in turn makes some dispatches request 64 KB of private segment per
lane and hit scratch limits — see the ticket for that side effect.)

## 7. `target nowait` + implicitly firstprivate scalar: not yet implemented

`07_nowait_firstprivate_nyi/nowait_scalar.F90`:

```fortran
!$omp target teams distribute parallel do nowait
do i = 1, n
   a(i) = a(i) + c          ! c is a plain host scalar
end do
```

```
error: not yet implemented: Unhandled clause privatization for deferred target tasks in omp.target operation
```

`firstprivate(c)` fails the same way; `map(to:c)` compiles and runs
correctly; a region that reads no scalar compiles. Since nearly every
Fortran kernel reads a dummy argument or a loop bound, `nowait` is
effectively unusable without rewriting every directive with explicit
`map(to:)` for each scalar.

---

## Not included on purpose

Two more items from our port notes are application bugs that ROCm merely
exposed, not toolchain defects, so they have no PoC here: an out-of-bounds
read `plold(i,0,jj)` that is harmless on CPU/NVIDIA and faults on ROCm, and
the OpenACC→OpenMP semantic difference for `private` scalars that are
initialised outside the kernel (OpenACC keeps the outer value, OpenMP
`private` does not — `firstprivate` is the correct translation).

Earlier reports from the same site — `rocprofv3` SIGSEGV in `ompt_post_init`
with OpenBLAS, and `HSA_STATUS_ERROR_OUT_OF_RESOURCES` on a DPX-partitioned
card — are in `../../TICKET_amd_rocm_support.md` and `../rocm_repro/`.
