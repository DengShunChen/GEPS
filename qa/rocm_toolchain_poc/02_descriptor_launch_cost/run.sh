#!/bin/bash
# PoC 02: per-launch cost of descriptor re-mapping for allocatable/assumed-shape arrays.
cd "$(dirname "$0")" && source ../common.sh
banner "build"
$FC $FFLAGS desc_launch.f90 -o desc_launch.x || exit 1
banner "timing (2000 launches each)"
./desc_launch.x
banner "data movement per launch (LIBOMPTARGET_INFO=-1, 1 timed launch each)"
echo "arguments  : number of kernel arguments = which kernel"
echo "   8  = k_explicit2 (2 explicit-shape),  62 = k_explicit (17 explicit-shape)"
echo "  88  = k_assumed  (17 assumed-shape),   87 = 17 allocatables referenced in scope"
LIBOMPTARGET_INFO=-1 ./desc_launch.x 1 2>&1 | grep -E 'Copying data|Entering OpenMP kernel' | awk '
/Entering OpenMP kernel/ { flush(); n=$0; sub(/.*with /,"",n); sub(/ arguments.*/,"",n); h2d=0; d2h=0; b=0; next }
/Copying data from host to device/ { h2d++; match($0,/Size=[0-9]+/); b+=substr($0,RSTART+5,RLENGTH-5) }
/Copying data from device to host/ { d2h++ }
function flush() { if (n!="") printf "  %3s arguments: %2d host->device copies (%4d bytes), %2d device->host copies\n", n, h2d, b, d2h }
END { flush() }' | sort | uniq -c | sed 's/^ *\([0-9]*\) /  x\1 /'
banner "-fno-defer-desc-map does not change it"
$FC $FFLAGS -fno-defer-desc-map desc_launch.f90 -o desc_launch_nodefer.x 2>/dev/null && ./desc_launch_nodefer.x | sed -n '4,5p'
