#!/bin/bash
# Build and run every PoC in order; output goes to the terminal and to $RESULTS/<n>.log.
# Total run time on one MI300X: ~6 minutes (PoC 03 and 04 dominate).
# Another toolchain:  ROCM_PATH=/opt/afar/therock-... RESULTS=results_afar ./run_all.sh
cd "$(dirname "$0")"
RESULTS=${RESULTS:-results}
ROCM=${ROCM_PATH:-/opt/rocm}
mkdir -p "$RESULTS"
for d in 0*_*/; do
  d=${d%/}
  echo; echo "==================== $d ===================="
  ./$d/run.sh 2>&1 | tee "$RESULTS/$d.log"
done
echo; echo "environment:" | tee "$RESULTS/env.log"
{ echo "ROCM_PATH=$ROCM"; $ROCM/bin/amdflang --version | head -1; cat $ROCM/.info/version 2>/dev/null | sed 's/^/ROCm /'; \
  rocminfo 2>/dev/null | grep -m1 'Name:.*gfx'; uname -r; } | tee -a "$RESULTS/env.log"
