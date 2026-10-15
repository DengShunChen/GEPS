#!/bin/bash
LOG=$1; OUT=$2; N=${3:-40}
until grep -q "forcast begin" "$LOG" 2>/dev/null; do sleep 3; done
sleep 40
PID=$(for p in $(pgrep -x tcogfs.x); do grep -q "State:.*[RS]" /proc/$p/status 2>/dev/null && [ "$(awk '{print $3}' /proc/$p/stat)" != Z ] && echo $p; done | head -1)
echo "pid=$PID" > "$OUT"
for i in $(seq 1 $N); do
  timeout 30 /opt/rocm/bin/rocgdb -p $PID -batch -ex "thread 1" -ex "bt 20" 2>/dev/null | grep -E "^#" >> "$OUT"
  echo "---" >> "$OUT"
done
