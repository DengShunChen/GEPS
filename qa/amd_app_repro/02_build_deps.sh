#!/bin/bash
# Build GEPS_LIB (OpenMPI, netCDF, OpenBLAS, ...) with the afar compilers into $GEPS_LIB_ROOT/install/rocm.
# About 30 minutes. Logs: $GEPS_LIB_ROOT/logs/.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
. "$HERE/env.sh"
cd "$GEPS_LIB_ROOT"
mkdir -p logs .shim
ln -sf "$(command -v tar)" .shim/gtar        # the libpng build calls gtar
export PATH=$PWD/.shim:$PATH

echo "[deps] OpenMPI 4.1.6 with $ROCM_PATH ($(date +%T))"
env -u GEPS_LIB_ROOT ROCM_PATH="$ROCM_PATH" GEPS_OPENMPI_DISTCLEAN=1 GEPS_OPENMPI_JOBS="$NJ" \
  bash src/openmpi/compile.sh > logs/openmpi.out 2>&1

echo "[deps] other libraries ($(date +%T))"
env -u GEPS_LIB_ROOT ROCM_PATH="$ROCM_PATH" GEPS_OPENMPI_MODULEFILES="$PWD/install/modulefiles" \
  bash -c '. /usr/share/modules/init/bash; module purge; ./build_all.sh --compiler rocm \
    libpng w3lib operlib zlib netcdf jasper fftw openblas g2lib dms dms38key mct geps_coupler' \
  > logs/build_all.out 2>&1
grep -E "^\[OK\]|FATAL" logs/build_all.out
ls install/rocm
