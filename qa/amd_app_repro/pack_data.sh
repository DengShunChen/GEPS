#!/bin/bash
# CWA side: package the input data AMD needs (about 10 GB uncompressed) into one tarball
# with a checksum. Not needed by AMD.
#   pack_data.sh [dmsdb dir] [output dir]
set -euo pipefail
SRC=${1:-/mlsteam/workspace/data/geps/dmsdb}
OUT=${2:-$PWD}
NAME=geps_TCo383L72_20220815_input
cd "$SRC"
tar -cf - ncep_ana_n1.ufs/TCo383l72_22081500 bckdms.ufs/BCK_TCo383_GI30S_xnew | lbzip2 > "$OUT/$NAME.tar.bz2"
(cd "$OUT" && sha256sum "$NAME.tar.bz2" > "$NAME.tar.bz2.sha256" && cat "$NAME.tar.bz2.sha256")
ls -la "$OUT/$NAME.tar.bz2"
