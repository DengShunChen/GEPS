#!/bin/bash
# Run TCo383L72 on 2 MI300X and check the result.
#   04_run.sh            1-hour forecast (6 steps), about 3 minutes once the rocFFT cache exists
#   04_run.sh 24         24-hour forecast (144 steps), about 10 minutes
#   VAR=value 04_run.sh  any runtime setting from README §4 overrides the default, e.g.
#                        GEPS_HSA_POOL_CAP_GB=96 ./04_run.sh 24
# The log lands in $GEPS_SRC/job/TCo383L72_IC_sample_rocm.log.<tag>.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
. "$HERE/env.sh"
H=${1:-1}
TAG=${TAG:-afar_${H}h_$(date +%Y%m%d_%H%M%S)}
cd "$GEPS_SRC/job"
LOG=TCo383L72_IC_sample_rocm.log.$TAG
[ "$H" != 1 ] && export GEPS_TAUE=$H
echo "[run] ${H} h forecast -> job/$LOG"
GEPS_FFT_PLAN_SETS=8 ./TCo383L72_IC_sample_rocm_2gpu > "$LOG" 2>&1
python3 "$HERE/check_run.py" "$LOG"
