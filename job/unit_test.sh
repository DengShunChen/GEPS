#!/bin/bash

set -x

MDIR=$(pwd)
echo "MDIR=${MDIR}"
cd ${MDIR}

group=sum

if [ $(id -gn) != $group ]; then
  exec sg $group "$0 $*"
fi

pjsub -N unit-test -L "vnode=1,vnode-core=32,ru=rscunit_pg01,rg=gpu-rd-large,gpu=8" --sparam wait-time=100 --interact -g ${group} qa/L0_test_subroutine/test.sh

