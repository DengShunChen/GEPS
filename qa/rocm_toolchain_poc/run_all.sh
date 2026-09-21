#!/bin/bash
# Build and run every PoC in order; output goes to the terminal and to results/<n>.log.
# Total run time on one MI300X: ~6 minutes (PoC 03 and 04 dominate).
cd "$(dirname "$0")"
mkdir -p results
for d in 0*_*/; do
  d=${d%/}
  echo; echo "==================== $d ===================="
  ./$d/run.sh 2>&1 | tee results/$d.log
done
echo; echo "environment:" | tee results/env.log
{ /opt/rocm/bin/amdflang --version | head -1; cat /opt/rocm/.info/version 2>/dev/null | sed 's/^/ROCm /'; \
  rocminfo 2>/dev/null | grep -m1 'Name:.*gfx'; uname -r; } | tee -a results/env.log
