# ROCm debugging tools for the GEPS GPU port

## Leak hunting: `LIBOMPTARGET_INFO=8`

Logs only OpenMP map-entry create/remove, each line carrying `Size=` and
`Name=<source file>`. Cheap enough to run on a full smoke test.

```bash
cd GEPS/job
LIBOMPTARGET_INFO=8 ./TCo383L72_IC_sample_rocm > /tmp/omptinfo.log 2>&1
```

Write the log to the container overlay (`/tmp`, ~1.7 TB), NOT to
`/mlsteam/workspace` (NFS, only ~57 GB free) - a full run produces ~1.3 GB.

Then:

```bash
python3 analyze_omptinfo.py  /tmp/omptinfo.log   # what is never removed, by source file
python3 omptinfo_timeline.py /tmp/omptinfo.log   # outstanding bytes per model phase
```

`omptinfo_timeline.py` is the one that distinguishes the two failure modes:
if outstanding bytes stay flat while free VRAM falls, the memory is being held
below the map layer (libomptarget pool / HSA allocator), not leaked by GEPS.

### Reading the output: two traps

**1. A run that ends in `abort()` inflates "never removed".** What is still
mapped at death is (legitimate long-lived maps) + (the working set of the call
that was executing). Neither is a leak. GEPS always aborts, so this always
applies.

**2. `create >> remove` only means a leak if that code ran many complete
calls.** We got this wrong once: `diabat_gpu.f90` showed 145 creates / 3
removes and was reported as a 49 GB leak. It was called exactly once and died
midway - its `exit data` (src/diabat_gpu.f90:4129) never ran. Nothing was
wrong with it.

So before calling anything a leak, check *when* its events happen:

```bash
# tag each create/remove with its log line number and compare against
# the phase markers the model itself prints
grep -an "Name=.*<file>" /tmp/omptinfo.log | head
```

A real leak looks like: the same code called N times, outstanding bytes rising
monotonically with N. That is what `omptinfo_timeline.py` shows per phase.

## Do NOT use `rocprof --hsa-trace` for this

Two independent reasons, both verified 2026-09-10:

1. Its HSA API trace is buffered and only written on normal exit. `SIGABRT`
   loses the whole thing (kernel/copy traces survive - they stream to
   `rpl_data_*`). GEPS always ends in abort, so the HSA table is always empty.
   `SIGTERM`/`SIGINT` do preserve it.
2. The profiler's own allocations dwarf what we are measuring: without it the
   tau=0 output stage uses ~36 GB of VRAM, under `rocprof --hsa-trace` the same
   stage exhausts all 192 GiB and aborts before reaching `iteration=1`.

## Tool bug reproducers

`run_repro.sh` covers two independent findings - see `../../TICKET_amd_rocm_support.md`:

1. `rocprofv3_ompt_segv.f90` - rocprofv3 SIGSEGVs in `ompt_post_init` for any
   amdflang + `--offload-arch` + OpenBLAS binary. Runs fine without the profiler.
2. `hsa_trace_mini.hip` - shows legacy `rocprof` HSA API interception works
   normally (so the interceptor is not broken; see point 1 above for why GEPS
   runs record nothing).

## `hipsolver_dsyevj_check.cc` - is the eigen solver to blame?

Mirrors the call sequence of `geps_rocm_dsyevj_dev`
(`src/rocm/hip_compat.cc`) on a known 4x4 symmetric matrix, inlined so it does
not drag in the whole compat layer. Seconds to build and run:

```bash
hipcc -O1 --offload-arch=gfx942 -o dsyevj_check hipsolver_dsyevj_check.cc \
      -lhipsolver -lamdhip64 && HIP_VISIBLE_DEVICES=0 ./dsyevj_check
```

Expected (verified 2026-09-11): `rc=0 devinfo=0`, eigenvalues
`0.254719 1.822717 3.177283 4.745281`, `max|eigenvector entry| = 0.777951`,
every column norm^2 `= 1.000000`.

So hipSOLVER itself is fine. If GEPS still sees `max|mx| > 1e12` coming out of
`eigen_mx_gpu`, the solve is not writing where we read - suspect
`!$omp target data use_device_addr(...)` and `device_ptr()` silently passing a
host address through. See the layer-8 section of ROCM_PORT_HANDOFF.md.
