#!/bin/bash
# PoC 06: runtime-sized private arrays corrupted unless LIBOMPTARGET_STACK_SIZE is raised.
cd "$(dirname "$0")" && source ../common.sh
banner "build"
$FC $FFLAGS privarr.F90 -o privarr_auto.x || exit 1
$FC $FFLAGS -DFIXED privarr.F90 -o privarr_fixed.x || exit 1
banner "automatic arrays, default stack (run 3x: count varies)"
for r in 1 2 3; do ./privarr_auto.x | head -1; done
banner "automatic arrays, LIBOMPTARGET_STACK_SIZE sweep"
for s in 1024 4096 16384 65536; do LIBOMPTARGET_STACK_SIZE=$s ./privarr_auto.x | head -1; done
banner "fixed-size arrays, default stack (control)"
./privarr_fixed.x
