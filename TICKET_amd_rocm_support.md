# AMD ROCm Support — bug reports from a gfx942 / MI300X site

Date: 2026-09-10
Reporter: numerical weather prediction model port (Fortran + OpenMP target offload + HIP)

Two independent issues below. **Issue 2 has a 15-line, self-contained, non-confidential
reproducer** and is the one we would most like looked at.

---

## Common environment

| | |
|---|---|
| ROCm | **7.2.2** (`/opt/rocm` -> `/opt/rocm-7.2.2`) |
| GPU | AMD Instinct MI300X, `gfx942`, Device ID `0x74a1` |
| GPU config | 2x physical cards, NPS1 / **SPX** (unpartitioned), 304 CU each, 192 GiB each |
| KFD nodes | node 4 (PCI `0000:46:00.0`), node 5 (PCI `0000:66:00.0`); `num_xcc=8`, `simd_count=1216` |
| Kernel | `6.8.0-138-generic` |
| amdgpu (KFD) | `6.16.13` |
| Fortran compiler | `AMD flang 22.0.0git` (`roc-7.2.2 26084 f58b06dce1f9c15707c5f808fd002e18c2accf7e`) |
| OpenMP runtime | `/opt/rocm-7.2.2/lib/llvm/lib/libomp.so`, `libomptarget.so.22.0git` |
| `rocprofv3` | version `1.1.0`, git_revision `671d39a71e33c49fba50b12e30b1aea45c5ed366` |
| OpenBLAS | 0.3.27 (built with the same amdflang toolchain) |
| MPI | self-built OpenMPI 4.1.6 |

> Note on source code: the application is a proprietary operational NWP model and its
> source cannot be shared. Everything below is reproduced with **standalone, non-confidential
> test cases we wrote for this report**, which are included inline.

---

## Issue 1 — `hsa_amd_memory_pool_allocate` returns `HSA_STATUS_ERROR_OUT_OF_RESOURCES` on a **partitioned** MI300X while VRAM is abundant

### Status: resolved for us by reconfiguring the GPU, reported for your awareness

We hit this for two days on a machine whose single MI300X was **DPX-partitioned**.
The same binary on **unpartitioned SPX** cards does not hit it at all. We can no longer
reproduce on demand (the machine has been reconfigured), so we are filing this as an
observation with a clean A/B rather than as an actionable repro.

### Symptom (on the DPX-partitioned configuration)

```
"PluginInterface" error: Failure to allocate device memory: "unknown or internal error"
error in hsa_amd_memory_pool_allocate: HSA_STATUS_ERROR_OUT_OF_RESOURCES:
The runtime failed to allocate the necessary resources. This error may also occur when
the core runtime library needs to spawn threads or create internal OS-specific events.
omptarget error: Call to getTargetPointer returned null pointer (device failure or illegal mapping).
omptarget fatal error 1: failure of target construct while offloading is mandatory
```

Triggered from an `!$omp target enter data create(...)` of ~17 plain scratch arrays,
on its **second** invocation. Offloading is `OMP_TARGET_OFFLOAD=MANDATORY`, `HSA_XNACK=1`.

### The A/B that we think is the interesting part

Identical binary, identical inputs, identical environment variables. **Only the GPU
partitioning mode differs:**

| | DPX-partitioned (1 card) | SPX unpartitioned |
|---|---|---|
| Result | `HSA_STATUS_ERROR_OUT_OF_RESOURCES` | **no such error at all** |
| VRAM in use over the whole run (`rocm-smi`, polled every 2s) | **flat at ~46 GB**, then allocation fails | 36 GB -> 77 GB -> 83 GB -> peak **135 GB** |
| VRAM reported total | 206,141,652,992 B | 206,141,652,992 B |

So on the partitioned device the runtime stopped being able to allocate at roughly
46 GB, while `rocm-smi` reported ~192 GiB total and ~150 GB free.

### What we verified was *not* the cause

- **Not actual VRAM exhaustion.** Polled `rocm-smi --showmeminfo vram` every 0.5s for the
  whole run; usage sat flat at ~46 GB out of 206 GB.
- **Not the number of arrays in the `enter data` clause.** We consolidated 15 of the 17
  arrays into a single buffer via F2008 pointer bounds-remapping (3 tracked entities
  instead of 17). **Identical failure, identical location** — so the exhausted resource
  does not scale with the mapped-entity count of the failing construct.
- **Not the libomptarget memory pool threshold.** `LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0`
  made it *worse* (failed earlier, in unrelated code, as a memory access fault);
  `=4294967296` behaved exactly like the default.
- **Not a GPU/driver fault.** A standalone HIP program doing `hipHostRegister` +
  `hipMalloc` + `hipMemcpy` of the same 510,935,040-byte buffer succeeded every time.

### The KFD detail that surprised us

We assumed the partition exposed smaller fixed queue quotas. **It does not** — the
per-node properties are byte-identical between the two configurations:

| `/sys/class/kfd/kfd/topology/nodes/N/properties` | DPX | SPX |
|---|---|---|
| `num_sdma_engines` | 2 | 2 |
| `num_sdma_queues_per_engine` | 8 | 8 |
| `num_gws` | 64 | 64 |
| VRAM total | 206,141,652,992 | 206,141,652,992 |

### Questions we would value an answer to

1. On a **partitioned** MI300X, is there a documented limit on how much a *single process*
   can allocate through `hsa_amd_memory_pool_allocate`, that is **lower** than the pool size
   reported by `hsa_amd_memory_pool_get_info` / `rocm-smi`? If so, where is it exposed?
2. `HSA_STATUS_ERROR_OUT_OF_RESOURCES` here clearly is not about VRAM bytes. **Which**
   resource class should we suspect, and is there a supported way to observe its usage?
3. Is the interaction of `HSA_XNACK=1` with a partitioned device known to reduce the
   allocatable ceiling?

---

## Issue 2 — `rocprofv3` SIGSEGVs during OMPT initialization, before the program starts

### Status: **fully reproducible, 15-line self-contained test case below**

Any binary built with `amdflang -fopenmp --offload-arch=gfx942` **that also links OpenBLAS**
crashes immediately when launched under `rocprofv3`. The program itself is fine — it only
crashes when the profiler is attached — and the crash happens during library initialization,
long before any GPU activity.

This makes `rocprofv3` unusable for us, which is why Issue 1 above has no profiler data.

### Reproducer

`rocprofv3_ompt_segv.f90`:

```fortran
program ompt_repro2
  implicit none
  integer, parameter :: n = 4
  real(kind=8) :: a(n,n), b(n,n), c(n,n)
  integer :: i
  real(kind=8) :: s
  a = 1.0d0; b = 2.0d0; c = 0.0d0
  ! OpenBLAS init happens on this call (gotoblas_init -> omp_get_num_places)
  call dgemm('N','N', n, n, n, 1.0d0, a, n, b, n, 0.0d0, c, n)
  s = 0.0d0
  !$omp target teams distribute parallel do reduction(+:s)
  do i = 1, 1024
     s = s + real(i, 8)
  end do
  print *, 'dgemm ok c(1,1)=', c(1,1), ' target sum=', s
end program
```

Build and run:

```bash
OB=/path/to/openblas-0.3.27/lib
/opt/rocm/bin/amdflang -fopenmp --offload-arch=gfx942 -o rocprofv3_ompt_segv \
    rocprofv3_ompt_segv.f90 -L$OB -lopenblas -Wl,-rpath,$OB

./rocprofv3_ompt_segv                                     # works:  "dgemm ok ... target sum= 524800."
rocprofv3 --kernel-trace -d out -- ./rocprofv3_ompt_segv  # SIGSEGV
```

### Crash

```
W[rocprofv3] tool initialization ::     0.001168 sec
*** Aborted at 1789024707 (unix time) ***
PC: @     0x720136ad90d4 pthread_mutex_lock
*** SIGSEGV (@0x2a8) received by PID 121711 (TID 0x72012c4e8f00) from PID 680; stack trace: ***
    @     0x720136adafb3 (unknown)
    @     0x72013d8d1a20 (unknown)
    @     0x720136a7e330 (unknown)
    @     0x720136ad90d4 pthread_mutex_lock
    @     0x720138089ab4 omp_get_num_devices
    @     0x72013a9d2413 ompt_post_init
    @     0x72013a950558 __kmp_do_middle_initialize()
    @     0x72013a95052c __kmp_middle_initialize
    @     0x72013a9cd228 __kmp_api_omp_get_num_places
    @     0x72013b02ec97 blas_get_cpu_number
    @     0x72013b02f6c7 gotoblas_init
    @     0x72013da2371f (unknown)
    @     0x72013da23824 (unknown)
    @     0x72013da3d5a0 (unknown)
```

Faulting address `0x2a8` looks like a null-pointer dereference at a fixed struct offset:
`ompt_post_init` calls `omp_get_num_devices`, which locks a mutex inside a structure that
is not yet initialized at that point in library startup.

### Bisecting results (each verified independently)

| Configuration | Result |
|---|---|
| offload (`--offload-arch=gfx942`) **+ OpenBLAS** | **SIGSEGV** |
| offload, `omp_get_num_places()` called from `main()`, no OpenBLAS | works |
| offload, `omp_get_num_places()` called from an ELF constructor in the executable, no OpenBLAS | works |
| **no** offload (no `libomptarget` linked) + OpenBLAS | works |
| any of the above without `rocprofv3` | works |

So the trigger needs **both** `libomptarget` linked **and** the first OpenMP API call to
come from a **shared library's** constructor (`gotoblas_init` inside `libopenblas.so`),
which runs before `libomp`/`libomptarget` are fully initialized. A constructor compiled
into the executable itself does **not** trigger it — the shared-library init ordering
appears to matter.

Trace flags make no difference (`--kernel-trace`, `--hsa-trace`, `--sys-trace` all crash
the same way).

### Side note

`rocprofv3` will not even load without `libdw.so.1` present (Ubuntu package `libdw1t64`);
it fails with `error while loading shared libraries: libdw.so.1`. Worth listing as an
explicit runtime dependency if it is not already.

---

## Retracted: earlier claim that `rocprof` v1 HSA API interception is broken

An earlier draft of this report contained a third issue claiming that the legacy `rocprof`
records **zero** HSA API calls. **That was our error and we are retracting it.**

Our "minimal case" ran `rocprof --hsa-trace` against `rocminfo` with a filter listing
`hsa_queue_create` / `hsa_amd_memory_pool_allocate` / `hsa_signal_create`. `rocminfo` never
calls any of those — it only enumerates agents — so zero records was the correct result.

Re-tested against a small HIP program that actually creates a queue and allocates, the
counters are all recorded correctly:

| HSA API | Calls |
|---|---|
| `hsa_queue_create` | 1 |
| `hsa_amd_memory_pool_allocate` | 4 |
| `hsa_amd_memory_pool_free` | 1 |
| `hsa_signal_create` | 9 |
| `hsa_amd_signal_create` | 84 |

Legacy `rocprof` HSA tracing works as documented on this system. No action needed.
