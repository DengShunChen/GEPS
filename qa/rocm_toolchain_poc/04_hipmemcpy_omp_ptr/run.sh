#!/bin/bash
# PoC 04: hipMemcpy* between libomptarget-allocated device buffers runs on the CPU.
cd "$(dirname "$0")" && source ../common.sh
CXX=${CXX:-$ROCM/bin/amdclang++}
banner "build"
$CXX -O2 -fopenmp --offload-arch=$ARCH -D__HIP_PLATFORM_AMD__ -I$ROCM/include memcpy_omp_ptr.cc -o memcpy_omp_ptr.x -L$ROCM/lib -lamdhip64 || exit 1
banner "64 MB copies"
./memcpy_omp_ptr.x 64
banner "HSA_XNACK=1 (as used by the application) - same"
HSA_XNACK=1 ./memcpy_omp_ptr.x 64
