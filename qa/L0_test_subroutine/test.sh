#!/bin/bash

MDIR=${PWD}
source ${MDIR}/qa/utils/setup.sh
cd ${MDIR}/build/test && ctest . --output-on-failure