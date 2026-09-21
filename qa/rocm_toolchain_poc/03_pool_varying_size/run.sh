#!/bin/bash
# PoC 03: libomptarget memory manager grows monotonically for varying-size maps.
cd "$(dirname "$0")" && source ../common.sh
banner "build"
$HIPCC -O2 --offload-arch=$ARCH -c freegb.cc -o freegb.o || exit 1
$FC $FFLAGS pool_growth.f90 freegb.o -L$ROCM/lib -lamdhip64 -o pool_growth.x || exit 1
banner "default memory manager (pool enabled)"
./pool_growth.x 16
banner "LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 (pool disabled)"
LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 ./pool_growth.x 16
banner "03b: cost of one enter/exit data pair, pool enabled vs disabled"
$FC $FFLAGS alloc_cost.f90 -o alloc_cost.x && { echo "-- default (pool)"; ./alloc_cost.x; echo "-- LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0"; LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 ./alloc_cost.x; }
banner "same program, every iteration the SAME size (pool reuse works)"
sed 's/n = 200000000_8 + it\*7000000_8/n = 250000000_8/' pool_growth.f90 > pool_same.f90
$FC $FFLAGS pool_same.f90 freegb.o -L$ROCM/lib -lamdhip64 -o pool_same.x && ./pool_same.x 6
