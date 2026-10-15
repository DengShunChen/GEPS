#!/bin/bash
# Update package lists
apt update
# Upgrade installed packages
apt upgrade -y
# Install development tools
apt install -y build-essential git git-lfs curl wget unzip openssh-server vim htop \
       	rsync environment-modules \
	libtirpc3t64 libhwloc15 libevent-core-2.1-7t64 \
	libevent-pthreads-2.1-7t64 libdw1t64 \
	libtirpc-dev libcurl4-openssl-dev
# 2026-10-15: packages installed by hand since 2026-09-10 (from /var/log/apt/history.log)
#   afar 24.3.0 toolchain (GEPS/qa/afar_app/INSTALL_AFAR.md): lbzip2 ksh m4
#   pptx/pdf reports: nodejs npm libreoffice-impress fonts-noto-cjk poppler-utils
#   bandwidth checks: rocm-bandwidth-test
apt install -y lbzip2 ksh m4 \
	nodejs npm libreoffice-impress fonts-noto-cjk poppler-utils \
	rocm-bandwidth-test
# Claude Code memory and session transcripts live on the overlay rootfs; restore the NFS copy
# (made by /mlsteam/workspace/claude_backup/backup.sh) if this is a fresh container
if [ ! -d /root/.claude/projects ] && [ -d /mlsteam/workspace/claude_backup/projects ]; then
	mkdir -p /root/.claude && rsync -a /mlsteam/workspace/claude_backup/projects/ /root/.claude/projects/
fi
# start ssh server
/etc/init.d/ssh start
# init ssh key
cp /mlsteam/workspace/.ssh/authorized_keys /root/.ssh/authorized_keys
chmod 600 /root/.ssh/authorized_keys

# for VS Code remote development ENV
export PATH=/opt/venv/bin:/opt/rocm/bin:$PATH

# Install Jupyter Notebook
#pip install jupyter

# 
export PATH="/mlsteam/workspace/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin:$PATH"

cd /mlsteam/workspace/init/amdgpu_top
cargo --version
rustc --version

cargo build --release --locked &&
./target/release/amdgpu_top --smi
