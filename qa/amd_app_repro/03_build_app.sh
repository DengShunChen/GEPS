#!/bin/bash
# Build tcogfs.x with afar into $GEPS_SRC/$GEPS_ROCM_BUILD_DIR, plus the HSA allocation cache
# (lib/hsa_front/libhsa-runtime64.so.1). A clean build takes about 20 minutes; most of it is the link.
# acc2omp switches pass through, e.g. GEPS_ACC2OMP_VEC_REDMAX=1 (README §6, item 6).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
. "$HERE/env.sh"
"$GEPS_SRC/qa/afar_app/build_afar_hip.sh" "${1:-$GEPS_ROCM_BUILD_DIR}" 2>&1 | tee "$HERE/build_app.log"
ls -la "$GEPS_SRC/${1:-$GEPS_ROCM_BUILD_DIR}/bin/tcogfs.x" "$GEPS_SRC/${1:-$GEPS_ROCM_BUILD_DIR}/lib/hsa_front/libhsa-runtime64.so.1"
