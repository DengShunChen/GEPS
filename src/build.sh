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

# load libs
export MDIR=$(cd ..;pwd)
. /usr/share/Modules/init/bash
module use ${MDIR}/modulefiles/${MACHINE}
module add modulefile.tcogfs.${MACHINE}

case $MACHINE in
fx10)
  export FFLAGS_COM="-c -X9 -Free -CcdRR8 -Cpp -Kfast,ocl,autoobjstack -x- -fw"
  export DEFINE_FLAGS="-DTIMING -DUSE_FFTW -DW3TAG"
  export MSG="-Ec -Qa,d,i,p,t,x -Koptmsg=2"
  export LIBS_COM="-SSL2BLAMP"
  export OMP="-Kopenmp,parallel"
;;
fx100)
  export FFLAGS_COM="-c -X9 -Free -CcdRR8 -Cpp -Kfast,ocl,autoobjstack -x- -fw"
  export DEFINE_FLAGS="-DTIMING -DUSE_FFTW -DW3TAG"
  export MSG="-Ec -Qa,d,i,p,t,x -Koptmsg=2"
  export LIBS_COM="-SSL2BLAMP"
  export OMP="-Kopenmp,parallel"
;;
pcc)
  export FFLAGS_COM="-c -O3 -Free -r8 -fpp -convert big_endian -w -ftz -ip -fp-model precise -align all -fno-alias -FR -xHost -fp-model fast=2 -no-heap-arrays -no-prec-div -no-prec-sqrt -fno-common -xCORE-AVX2 -auto"
  export DEFINE_FLAGS="-DTIMING -DUSE_FFTW "
  export MSG=""
  export LIBS_COM="-mkl=cluster -qopt-matmul"
  export LDFLAGS_COM='-static-intel -qopt-report'
  export OMP="-parallel"
;;
*)
  echo "${0} : Unknown machine --> ${MACHINE}" ; exit
;;

esac

make clean
make -j24

