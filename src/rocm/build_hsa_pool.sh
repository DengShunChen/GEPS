#!/bin/bash
# Builds libgeps_hsa_pool.so (LD_PRELOAD allocation cache, see geps_hsa_pool.cc).
# Kept out of CMake on purpose: touching src/CMakeLists.txt reconfigures the
# GPU tree (ROCM_PORT_HANDOFF.md §8). Output: build_rocm_hip/lib/libgeps_hsa_pool.so
# Also builds <out>/hsa_front/libhsa-runtime64.so, the same cache for runtimes whose
# libomptarget dlopen()s HSA instead of linking it (afar 24.x; see GEPS_HSA_FRONT).
set -e
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/../../build_rocm_hip/lib}
mkdir -p "$out"
ROCM=${ROCM_PATH:-/opt/rocm}
CXX=${CXX:-$ROCM/llvm/bin/clang++}
"$CXX" -O2 -fPIC -shared -std=c++17 -I$ROCM/include "$here/geps_hsa_pool.cc" -o "$out/libgeps_hsa_pool.so" -ldl
echo "built $out/libgeps_hsa_pool.so"
# only for runtimes whose libomptarget does not link HSA (ROCm 7.2.2 links it: preload suffices)
omptarget=$(ls "$ROCM"/lib/llvm/lib/libomptarget.so.* 2>/dev/null | head -1)
if readelf -d "$omptarget" 2>/dev/null | grep -q "libhsa-runtime64"; then
  echo "skip hsa_front: $omptarget links HSA directly"; exit 0
fi
mkdir -p "$out/hsa_front"
# the front owns the soname libhsa-runtime64.so.1; the real runtime is its dependency
# under the name libhsa-real.so (link-time stub fixes the DT_NEEDED name, runtime symlink)
real=$(readlink -f "$ROCM/lib/libhsa-runtime64.so.1")
stub=$(mktemp -d)
echo "" | "$CXX" -x c -shared -fPIC -o "$stub/libhsa-real.so" -Wl,-soname,libhsa-real.so -
printf 'ROCR_1 { global: hsa_amd_memory_pool_allocate; hsa_amd_memory_pool_free; };\n' > "$stub/ver.map"
ln -sfn "$real" "$out/hsa_front/libhsa-real.so"
# the real runtime finds its own deps via $ORIGIN/rocm_sysdeps/lib; $ORIGIN is now hsa_front/
[ -d "$(dirname "$real")/rocm_sysdeps" ] && ln -sfn "$(dirname "$real")/rocm_sysdeps" "$out/hsa_front/rocm_sysdeps"
"$CXX" -O2 -fPIC -shared -std=c++17 -DGEPS_HSA_FRONT -I$ROCM/include "$here/geps_hsa_pool.cc" \
  -o "$out/hsa_front/libhsa-runtime64.so.1" -Wl,-soname,libhsa-runtime64.so.1 \
  -L"$stub" -Wl,--no-as-needed -lhsa-real -Wl,-rpath,'$ORIGIN' -ldl \
  -Wl,--version-script,"$stub/ver.map"
rm -rf "$stub"
echo "built $out/hsa_front/libhsa-runtime64.so.1 (front for $real)"
