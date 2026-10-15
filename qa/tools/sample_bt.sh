#!/bin/bash
# poor-man's host profiler: wait for the forecast phase, then take gdb backtraces of rank 0
LOG=$1; OUT=$2
until grep -q "forcast begin" "$LOG" 2>/dev/null; do sleep 5; done
sleep 45   # skip step 1 (physics init)
PID=$(pgrep -f "build_rocm_hip/bin/tcogfs.x" | head -1)
echo "sampling pid $PID" > "$OUT"
for i in $(seq 1 40); do
  timeout 20 rocgdb -p $PID -batch -ex "bt 14" 2>/dev/null | grep -E "^#" >> "$OUT"
  echo "---" >> "$OUT"
done
