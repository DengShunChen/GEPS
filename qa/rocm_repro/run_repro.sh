#!/bin/bash
# Two independent ROCm tool checks. Run from this directory.
#   1) rocprofv3_ompt_segv.f90  -> reproduces a rocprofv3 SIGSEGV (tool bug)
#   2) hsa_trace_mini.hip       -> shows rocprof v1 HSA API tracing DOES work
set -u
OB=${OPENBLAS_LIB:-/mlsteam/workspace/data/geps/GEPS_LIB/install/rocm/openblas-0.3.27/lib}
export HIP_VISIBLE_DEVICES=0 ROCR_VISIBLE_DEVICES=0

echo "=== 1. rocprofv3 OMPT SIGSEGV ==="
/opt/rocm/bin/amdflang -fopenmp --offload-arch=gfx942 -o rocprofv3_ompt_segv \
    rocprofv3_ompt_segv.f90 -L"$OB" -lopenblas -Wl,-rpath,"$OB" || exit 1
echo "--- bare run (expected: works) ---";      ./rocprofv3_ompt_segv; echo "rc=$?"
echo "--- under rocprofv3 (expected: SIGSEGV in ompt_post_init) ---"
rocprofv3 --kernel-trace -d rpv3_out -- ./rocprofv3_ompt_segv 2>&1 | tail -18; echo "rc=${PIPESTATUS[0]}"

echo; echo "=== 2. rocprof v1 HSA API trace (expected: counters ARE recorded) ==="
hipcc --offload-arch=gfx942 -o hsa_trace_mini hsa_trace_mini.hip 2>/dev/null || exit 1
rocprof --hsa-trace -i rocprof_hsa_filter.txt -o hsa_out.csv ./hsa_trace_mini >/dev/null 2>&1
echo "--- hsa_stats.csv ---"; cat hsa_out.hsa_stats.csv 2>/dev/null | head -12
