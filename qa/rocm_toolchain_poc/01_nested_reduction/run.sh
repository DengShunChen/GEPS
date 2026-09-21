#!/bin/bash -x
# PoC 01: nested `parallel do reduction` under `target teams distribute` returns 0.
cd "$(dirname "$0")" && source ../common.sh
for opt in -O2 -O0; do
  banner "amdflang $opt"
  $FC $opt -fopenmp --offload-arch=$ARCH nested_reduction.f90 -o nested_reduction$opt.x || exit 1
  ./nested_reduction$opt.x
done
banner "host OpenMP only (no offload) - reference behaviour"
$FC -O2 -fopenmp nested_reduction.f90 -o nested_reduction_host.x && OMP_TARGET_OFFLOAD=DISABLED ./nested_reduction_host.x
