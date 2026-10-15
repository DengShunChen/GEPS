#!/bin/bash
LOG=$1; OUT=$2
until grep -q "forcast begin" "$LOG" 2>/dev/null; do sleep 3; done
sleep 2
PID=$(pgrep -x tcogfs.x | while read p; do [ -d /proc/$p ] && grep -q State.*R /proc/$p/status 2>/dev/null && echo $p; done | head -1)
[ -z "$PID" ] && PID=$(pgrep -x tcogfs.x | tail -1)
echo "pid=$PID" > "$OUT"
for i in $(seq 1 30); do
  timeout 20 rocgdb -p $PID -batch -ex "bt 16" 2>/dev/null | grep -E "^#" >> "$OUT"
  echo "---" >> "$OUT"
done
