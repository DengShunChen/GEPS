#!/bin/bash

MACHINE="a100"

MDIR=${PWD}
source ${MDIR}/qa/utils/setup.sh
#cd ${MDIR}/build_${machine}/test
#mpiexec -n 8 test_trngra3 #flase -> test_nccl test_tranrs, ok test_initial test_tendget(mini intgrt)
cd ${MDIR}/build_${machine}/test && ctest . --output-on-failure

