#!/bin/bash
# Install the afar 24.3.0 SDK at /opt/afar/current and the system packages the build needs.
# Skip the SDK part if afar 24.3.0 is already installed; set ROCM_PATH in env.sh to it.
set -euo pipefail
PKG=therock-afar-24.3.0-multiarch-10.1.0-592954c

# ksh: operlib's Makefile uses SHELL=ksh; m4: netCDF configure; lbzip2: parallel unpack;
# environment-modules: GEPS_LIB/GEPS load their toolchain through `module`.
apt-get install -y build-essential git curl rsync environment-modules ksh m4 lbzip2 libtirpc-dev libcurl4-openssl-dev

if [ ! -x /opt/afar/$PKG/bin/amdflang ]; then
  mkdir -p /opt/afar && cd /opt/afar
  [ -f afar.tar.bz2 ] || curl -sSfL -o afar.tar.bz2 https://repo.radeon.com/rocm/misc/flang/$PKG.tar.bz2   # 9.1 GB
  tar -I lbzip2 -xf afar.tar.bz2
fi
ln -sfn /opt/afar/$PKG /opt/afar/current
/opt/afar/current/bin/amdflang --version | head -1
