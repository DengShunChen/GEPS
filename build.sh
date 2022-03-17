#!/bin/bash
#------------------------------------------------------------------------------
#
#  build GFS forecast model
# 
# 1. select targe machine(fx10, fx100, or pcc) 
# 2. run this building up script
#   ./build.sh [MACHINE]
#                                                                Deng-Shun Chen
#                                                                   2020-04-28
#-------------
# add feature:                                        
#   1. for RSM-IO                                                  CHEN,YING-JU
#      use "modulefile.tcogfs.${MACHINE}.rsm" files                  2021-07-30
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
machines='fx1000 fx100 fx10 pcc'
[[ $machines =~ (^|[[:space:]])$MACHINE($|[[:space:]]) ]] && known='True' || known='False'
if [ "${known}" == 'True' ] ; then
  echo "${HOSTNAME} : Build ${MACHINE} executable"
else
  echo "Fatal Error : $0: Unknown machine --> ${MACHINE}" ; exit
fi

if [ "${MACHINE}" == 'fx10' ] ; then 
  hostnames='login07 login08 login05 login06'
  [[ $hostnames =~ (^|[[:space:]])$HOSTNAME($|[[:space:]]) ]] && known='True' || known='False' 
  if [ "${known}" == 'False' ] ; then
    echo "Fatal Error : Build ${MACHINE} executable, please move to login07/08 for inside HPC, login05/06 for outside HPC !" 
    exit
  fi
fi
if [ "${MACHINE}" == 'fx100' ] ; then
  hostnames='login11 login12 login15 login16'
  [[ $hostnames =~ (^|[[:space:]])$HOSTNAME($|[[:space:]]) ]] && known='True' || known='False' 
  if [ "${known}" == 'False' ] ; then
   echo "Fatal Error : Build ${MACHINE} executable, please move to login11/12 for inside HPC, login15/16 for outside HPC !" 
   exit
  fi
fi
if [ "${MACHINE}" == 'fx1000' ] ; then
  hostnames='h6ln12 h6ln13 h6ln15 h6ln16 h6ln17 h6ln18 h6ln19 h6ln23'
  [[ $hostnames =~ (^|[[:space:]])$HOSTNAME($|[[:space:]]) ]] && known='True' || known='False'
  if [ "${known}" == 'False' ] ; then
   echo "Fatal Error : Build ${MACHINE} executable, please move to login12/13/15/16/17/18/19 for inside HPC, login23 for outside HPC !"
   exit
  fi
fi

set -x

# load libs
export MDIR=$(pwd)
. /usr/share/Modules/init/bash
module purge
module use  ${MDIR}/modulefiles
module av
module show modulefile.tcogfs.${MACHINE}
module load modulefile.tcogfs.${MACHINE}
module list
module unuse ${MDIR}/modulefiles

# compile
cd src/
make clean
make -j24

