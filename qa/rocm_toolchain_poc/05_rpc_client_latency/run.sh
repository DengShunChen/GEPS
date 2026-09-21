#!/bin/bash
# PoC 05: device-side Fortran I/O in the image starts an RPC server thread.
cd "$(dirname "$0")" && source ../common.sh
banner "build"
$FC $FFLAGS rpc_latency.F90 -o rpc_plain.x || exit 1
$FC $FFLAGS -DDEVIO rpc_latency.F90 -o rpc_devio.x || exit 1
$FC $FFLAGS -DNOWAIT rpc_latency.F90 -o rpc_plain_nowait.x || exit 1
$FC $FFLAGS -DNOWAIT -DDEVIO rpc_latency.F90 -o rpc_devio_nowait.x || exit 1
chmod +x check_image.sh
banner "device image contents"
./check_image.sh rpc_plain.x
./check_image.sh rpc_devio.x
banner "threads at runtime (the extra thread is libomptarget's RPC server)"
for x in rpc_plain rpc_devio; do
  ./$x.x >/dev/null & pid=$!; sleep 1.5; echo "$x: $(ls /proc/$pid/task | wc -l) threads"; wait $pid
done
for v in "" _nowait; do
  banner "plain$v (no device I/O)";  ./rpc_plain$v.x
  banner "devio$v (never-taken print in the kernel)"; ./rpc_devio$v.x
done
