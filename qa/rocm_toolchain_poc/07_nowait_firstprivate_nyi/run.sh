#!/bin/bash
# PoC 07: `target nowait` + implicitly firstprivate scalar does not compile.
cd "$(dirname "$0")" && source ../common.sh
for v in "" -DEXPLICIT -DMAPTO -DNOSCALAR; do
  banner "variant: ${v:-(default, implicit firstprivate)}"
  if $FC $FFLAGS $v nowait_scalar.F90 -o nowait_scalar$v.x 2> build$v.log; then
    ./nowait_scalar$v.x
  else
    grep -m1 'error:' build$v.log
  fi
done
