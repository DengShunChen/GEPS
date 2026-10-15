#!/bin/bash
# Build tcogfs.x with the AMD afar 24.3.0 flang (TheRock 10.x SDK) into a separate tree,
# leaving build_rocm_hip (ROCm 7.2.2) untouched. Same configure line as ROCM_PORT_HANDOFF.md §8.
#   usage: build_afar_hip.sh [build dir name, default build_afar_hip]
#   acc2omp switches pass through the environment, e.g.
#     GEPS_ACC2OMP_VEC_REDMAX=1 build_afar_hip.sh build_afar_hip_redmax    (PoC 01 on the app)
#     GEPS_ACC2OMP_KEEP_DEVIO=sflx build_afar_hip.sh build_afar_hip_devio  (PoC 05 on the app)
set -e
GEPS=$(cd "$(dirname "$0")/../.." && pwd)
BLD=${1:-build_afar_hip}
export ROCM_PATH=${ROCM_PATH:-/opt/afar/current}
export GEPS_LIB_ROOT=${GEPS_LIB_ROOT:-$(cd "$GEPS/../GEPS_LIB_afar" && pwd)}
export GEPS_OPENMPI_MODULEFILES=$GEPS_LIB_ROOT/install/modulefiles
export GEPS_COMPILER=rocm
. /usr/share/modules/init/bash && module purge
source "$GEPS_LIB_ROOT/env.sh" && load_geps_compiler
module use "$GEPS_LIB_ROOT/install/modulefiles" && module load geps/1.0
export ROCM_PATH=/opt/afar/current GEPS_LIB_ROOT
export COMP_MP=mpifort C_COMP_MP=mpicc CXX_COMP_MP=mpicxx
cd "$GEPS"
echo "[afar] $(amdflang --version | head -1)  mpifort=$(command -v mpifort)  -> $BLD"
if [ ! -f "$BLD/CMakeCache.txt" ]; then
  cmake -B"$BLD" -S. -DCMAKE_BUILD_TYPE=Release -DUSE_RSM=OFF -DCWBSUM=OFF -DTIMCOMCPL=OFF -DUSE_aeroclx=OFF \
        -DGEPS_COMPILER=rocm -DUSE_CUDA=OFF -DUSE_HIP=ON -DUSE_ACC=OFF -DUSE_PAR=OFF -DUSE_OMP=ON -DUSE_NDMS=OFF \
        -DAMD_GPU_ARCHS=gfx942 -DCMAKE_Fortran_COMPILER="$ROCM_PATH/bin/amdflang" \
        -DCMAKE_C_COMPILER="$ROCM_PATH/bin/amdclang" -DCMAKE_CXX_COMPILER="$ROCM_PATH/bin/amdclang++"
fi
# Clean-tree ordering bug: rocm/geps_acc_wait.f90 (target geps_rocm_gpu) does `use rank`, but
# mod_rank.f90 belongs to the main library, which depends on geps_rocm_gpu. Incremental trees
# never notice (rank.mod is left over); a fresh tree needs rank.mod (and nccl.mod) first.
if [ ! -f "$BLD/src/rank.mod" ]; then
  fl=$(grep "^Fortran_FLAGS" "$BLD/src/CMakeFiles/tcogfs.dir/flags.make" | cut -d= -f2-)
  in=$(grep "^Fortran_INCLUDES" "$BLD/src/CMakeFiles/geps_rocm_gpu.dir/flags.make" | cut -d= -f2-)
  df=$(grep "^Fortran_DEFINES" "$BLD/src/CMakeFiles/tcogfs.dir/flags.make" | cut -d= -f2-)
  tmpo=$(mktemp -d)
  for f in src/rocm/cudafor.f90 src/rocm/nccl.f90 src/mod_rank.f90; do
    (cd "$BLD/src" && "$ROCM_PATH/bin/amdflang" $df $in $fl -c "$GEPS/$f" -o "$tmpo/$(basename $f).o") || exit 1
  done
fi
cmake --build "$BLD" --target tcogfs.x -j ${NJ:-32}
CXX=$ROCM_PATH/llvm/bin/clang++ src/rocm/build_hsa_pool.sh "$GEPS/$BLD/lib"
ls -la "$BLD/bin/tcogfs.x"
