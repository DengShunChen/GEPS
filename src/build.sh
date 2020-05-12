#!/bin/bash
#------------------------------------------------------------------------------
#
#  build GFS forecast model
# 
# 1. select targe machine(fx10, fx100, or pcc) 
# 2. run this building up script
#   ./build.sh
#                                                                Deng-Shun Chen
#                                                                   2020-04-28
#------------------------------------------------------------------------------
#
# MACHINE : fx10, fx100, pcc
#
if [ $# == 1 ] ; then
  export MACHINE=${1}
else
  echo "usage: $0 [MACHINE]"
  exit
fi
machines='fx100 fx10 pcc'
[[ $machines =~ (^|[[:space:]])$MACHINE($|[[:space:]]) ]] && known='True' || known='False'
if [ "${known}" == 'True' ] ; then
  echo "${HOSTNAME} : Build ${MACHINE} executable"
else
  echo "Fatal Error : $0: Unknown machine --> ${MACHINE}" ; exit
fi

if [ "${MACHINE}" == 'fx10' ] &&  [ "${HOSTNAME}" != 'login07' ] && [ "${HOSTNAME}" != 'login08' ] && [ "${HOSTNAME}" != 'login05' ] && [ "${HOSTNAME}" != 'login06' ] ; then
   echo "Fatal Error : Build ${MACHINE} executable, please move to login07/08 for inside HPC, login05/06 for outside HPC !" 
   exit
fi
if [ "${MACHINE}" == 'fx100' ] &&  [ "${HOSTNAME}" != 'login11' ] && [ "${HOSTNAME}" != 'login12' ] && [ "${HOSTNAME}" != 'login15' ] && [ "${HOSTNAME}" != 'login16' ]; then
   echo "Fatal Error : Build ${MACHINE} executable, please move to login11/12 for inside HPC, login15/16 for outside HPC !" 
   exit
fi

set -x

# load libs
export MDIR=$(cd ..;pwd)
. /usr/share/Modules/init/bash
module use  ${MDIR}/modulefiles
module load modulefile.tcogfs.${MACHINE}
module av
module list



# compile
make clean
make -j12

