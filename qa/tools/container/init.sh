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
