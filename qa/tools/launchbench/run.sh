#!/bin/bash
export GEPS_LIB_ROOT=/mlsteam/workspace/data/geps/GEPS_LIB
export GEPS_COMPILER=rocm
. /usr/share/modules/init/bash && module purge
source $GEPS_LIB_ROOT/env.sh >/dev/null && load_geps_compiler >/dev/null
cd $(dirname $0)
[ lb.x -nt lb.f90 ] || amdflang -O2 -fopenmp --offload-arch=gfx942 lb.f90 -o lb.x || exit 1
export HIP_VISIBLE_DEVICES=${DEV:-0} ROCR_VISIBLE_DEVICES=${DEV:-0} OMP_TARGET_OFFLOAD=MANDATORY
echo "== default env (threshold 0, stack 65536, xnack 1) =="
HSA_XNACK=1 LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 LIBOMPTARGET_STACK_SIZE=65536 ./lb.x
echo "== no env =="
./lb.x
echo "== xnack 0 only =="
HSA_XNACK=0 LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 LIBOMPTARGET_STACK_SIZE=65536 ./lb.x
echo "== threshold default, others as job =="
HSA_XNACK=1 LIBOMPTARGET_STACK_SIZE=65536 ./lb.x
