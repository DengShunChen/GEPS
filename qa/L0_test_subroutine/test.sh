#!/bin/bash

MACHINE="a100"

MDIR=${PWD}
source ${MDIR}/qa/utils/setup.sh
cd ${MDIR}/build_${machine}/test && ctest . --output-on-failure