#!/bin/bash
#------------------------------------------------------------------------------
#
#  build GFS with w3tag 
#
# 1. Add -DW3TAG into Makfile
#   FFLAGS  = -c -X9 -Free -CcdRR8 -Cpp -DTIMING -DMULTIPLE -DUSE_FFTW  -DW3TAG \
#             -Kfast,parallel,openmp,ocl,autoobjstack \
#             -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)
# 2. run this building up script
#   ./build.sh
#                                                                Deng-Shun Chen
#                                                                   2020-04-28
#------------------------------------------------------------------------------
MACHINE=fx100

. /usr/share/Modules/init/bash
module use /nwpr/gfs/xb80/modulefiles/nceplibs/${MACHINE}
module add w3nco

make clean
make -j24

