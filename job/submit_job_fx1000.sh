#!/bin/ksh
#PJM -L "node=9:noncont"
#PJM -L rscgrp=small
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=10:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=384"
#PJM -j 
#PJM -g sum
#PJM -N TCo383-MASOP
#PJM -o %j.log
#PJM -e %j.err

set -x

OMP=1

export I_MPI_PIN_CELL=core
export JCAP=383
export MDIR=`cd ../ ; pwd`
if [ $JCAP = 639 ] ; then
  export NPEX=4
  export NPEY=384
elif [ $JCAP = 383 ] ; then
  export NPEX=4
  export NPEY=96
fi

export MPI=$((${NPEX}*${NPEY}))
export machine=fx1000
mdir=`cd ../ ; pwd`
export dtg=`cat ${mdir}/job/run_date |cut -c3-10`
#==============================================================================================================#
source /users/xa09/sample/setup_mpi+omp.fx1000 $OMP

./regression.ksh

echo "Ending at: " `date`

