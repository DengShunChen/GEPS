#!/bin/bash
# kfdmon.sh LOG OUT: once per second, for every KFD process dir (named by HOST pid), the
# current tau and evicted_ms / faults / page_in / page_out per GPU. Stops when no live tcogfs.x.
LOG=$1; OUT=$2; : > $OUT
live() { for p in $(pgrep -x tcogfs.x); do [ "$(awk '{print $3}' /proc/$p/stat 2>/dev/null)" != Z ] && return 0; done; return 1; }
until live; do sleep 1; done
while live; do
  tau=$(grep -o "tau= *[0-9.]*, *Timing" "$LOG" 2>/dev/null | tail -1 | awk '{print $2}' | tr -d ,)
  for d in /sys/class/kfd/kfd/proc/*; do
    s="$(date +%s) tau=$tau hpid=$(basename $d)"
    for g in 29122 43483; do
      [ -f $d/stats_$g/evicted_ms ] || continue
      s="$s g$g ev=$(cat $d/stats_$g/evicted_ms) f=$(cat $d/counters_$g/faults) in=$(cat $d/counters_$g/page_in) out=$(cat $d/counters_$g/page_out)"
    done
    echo "$s" >> $OUT
  done
  sleep 1
done
