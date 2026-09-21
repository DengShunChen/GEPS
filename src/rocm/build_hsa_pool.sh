#!/bin/bash
# Builds libgeps_hsa_pool.so (LD_PRELOAD allocation cache, see geps_hsa_pool.cc).
# Kept out of CMake on purpose: touching src/CMakeLists.txt reconfigures the
# GPU tree (ROCM_PORT_HANDOFF.md §8). Output: build_rocm_hip/lib/libgeps_hsa_pool.so
set -e
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/../../build_rocm_hip/lib}
mkdir -p "$out"
CXX=${CXX:-/opt/rocm/llvm/bin/clang++}
"$CXX" -O2 -fPIC -shared -std=c++17 -I/opt/rocm/include "$here/geps_hsa_pool.cc" -o "$out/libgeps_hsa_pool.so" -ldl
echo "built $out/libgeps_hsa_pool.so"
