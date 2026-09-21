# Shared settings for every PoC. Source it, do not execute it.
# Requires only a stock ROCm install; no application code, no MPI, no BLAS.
ROCM=${ROCM_PATH:-/opt/rocm}
FC=${FC:-$ROCM/bin/amdflang}
HIPCC=${HIPCC:-$ROCM/bin/hipcc}
ARCH=${OFFLOAD_ARCH:-gfx942}
FFLAGS="-O2 -fopenmp --offload-arch=$ARCH"
export HIP_VISIBLE_DEVICES=${HIP_VISIBLE_DEVICES:-0}
export ROCR_VISIBLE_DEVICES=${ROCR_VISIBLE_DEVICES:-0}
export OMP_TARGET_OFFLOAD=MANDATORY
export PATH=$ROCM/bin:$ROCM/lib/llvm/bin:$PATH
export LD_LIBRARY_PATH=$ROCM/lib:$ROCM/lib/llvm/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
banner() { echo; echo "################ $* ################"; }
