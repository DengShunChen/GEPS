#!/bin/ksh
#PJM -L "node=24:noncont"
##PJM -L rscgrp=large-x
#PJM -L elapse=8:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=192"
#PJM -j
#PJM -N TCo383L72_fx10
#PJM -s

set -x

MPI=192
OMP=2
machine=fx10

source /users/xa09/sample/setup_mpi+omp.${machine} $OMP

#==============================================================================================================#
LANG_HOME=${LANG_HOME-"/opt/FJSVmxlang"}
TOFU_HOME=${TOFU_HOME-"/opt/FJSVpxtof/sparc64fx"}
PATH=${LANG_HOME}/bin:$PATH; export PATH
LD_LIBRARY_PATH=$LANG_HOME/lib64:$TOFU_HOME/lib64; export LD_LIBRARY_PATH
LD_LIBRARY_PATH=/opt/FJSVpxtof/sparc64fx/lib64:${LD_LIBRARY_PATH}; export LD_LIBRARY_PATH
LD_LIBRARY_PATH=/users/xa09/operlib/lib:/package/fx10/dms/dms.v4/lib:/opt/FJSVXosSclib/lib64:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
LANG=C; export LANG

OMPI_NOT_USE_PLE_TMPDIR=1; export OMPI_NOT_USE_PLE_TMPDIR

OMP_NUM_THREADS=$PARALLEL; export OMP_NUM_THREADS
THREAD_STACK_SIZE=1048576; export THREAD_STACK_SIZE
FLIB_FASTOMP=TRUE; export FLIB_FASTOMP
PSZ=${PSZ-"32MB"}
LPG="/opt/FJSVxosmmm/sbin/lpgparm -t $PSZ -d $PSZ -h $PSZ -s $PSZ -p $PSZ"
XOS_MMM_L_ARENA_FREE=2; export XOS_MMM_L_ARENA_FREE

#-- enviornment
 user=`whoami`
 datamv='login07'
 dmsdb_home=/nwpr/gfs/${user}/data2/dmsdb
 MDIR=/nwpr/gfs/${user}/tco639l72
 DMSPATH=/package/${machine}/dms/dms.v4/bin
 GFSDIR=$MDIR
 GFSFIX=$MDIR/fix
 GFSWRK=${GFSDIR}/wrk
 rm -rf $GFSWRK
 mkdir -p $GFSWRK

#-- dms data
 dd='180901'
 dtg=${dd}'00'

 #EXP='TCo383L72_fpvs'
 EXP='TCo383L72_omix'
 indmsowner='xb118'
 indmshead='TCo383L72_'
 indmstail='GIMGM'
 indmsdb='TCo383L72'

 outdmshead=${EXP}_git
 outdmstail=${indmstail}
 outdmsdb=${indmsdb}

#-- executable
 EXEC='MTCo383L72_'${machine}

#---------------------------------------------------------#
 indmsfile=${indmsdb}/20${dd}@${dmsdb_home}/ana
 outdmsfile=${outdmshead}${dtg}${outdmstail}@${dmsdb_home}/archive

 if [ $indmsowner != $user ] ; then
  ${DMSPATH}/rdmsdbcrt -p ufs  $indmsdb
  ssh $datamv -n ${DMSPATH}/rdmscpy -l34 -k "*" ${indmsfile}@${indmsowner}@${datamv} ${indmsfile}
 fi

#-- write out running date tag
 echo $dtg > ${GFSWRK}/crdate

#-- write out file list
cat > ${GFSWRK}/filist << EOF
 &filst
 ifilin='ANADMS',
 ifilout='FCSTDMS',
 cwbout='${GFSWRK}/tmpdir/gfs_cwbout',
 bckfile='BCKOPS',
 phyout='${GFSWRK}/tmpdir/gfs_phyout',
 namlsts='${GFSWRK}/namlsts_tco383',
 crdate='${GFSWRK}/crdate',
 ocards='${GFSWRK}/ocards',
 cntrl='${GFSWRK}/gfsctl',
 &end
EOF


 export GFSDIR DMSPATH
 export NWPETC=${GFSDIR}/etc
 export NWPETCGLB=${GFSWRK}
 export GLB_TYPHINI="/nwpr/gfs/xb118/data2/obstyphoon"
 export GLB_WTYPHINI="/nwpr/gfs/xb118/data2/obstyphoon"
 export FIXDIR=${GFSFIX}

 export ANADMS=${indmsfile}
 export FCSTDMS=${outdmsfile}
 export BCKOPS=BCK_TCo383_GI30S@${dmsdb_home}/bckdms

${DMSPATH}/rdmspurge -f FCSTDMS
${DMSPATH}/rdmscrt -l34 FCSTDMS

 export FLIB_CNTL_BARRIER_ERR=FALSE



export O3FORC=${O3FORC:-${FIXDIR}/global_o3prdlos.f77}
export O3CLIM=${O3CLIM:-${FIXDIR}/global_o3clim.txt}
export AEROSOL_FILE=${AEROSOL_FILE:-${FIXDIR}/global_climaeropac_global.txt}
export EMMISSIVITY_FILE=${EMMISSIVITY_FILE:-${IXDIR}/global_sfc_emissivity_idx.txt}


cd $GFSWRK
#====================================================================
ln -fs $O3FORC fort.28
ln -fs $O3CLIM fort.48
ln -fs $AEROSOL_FILE  aerosol.dat
ln -fs $EMMISSIVITY_FILE sfc_emissivity_idx.txt

ln -fs $FIXDIR/* .
cp $NWPETC/gfsctl $GFSWRK/gfsctl
cp $NWPETC/ocards $GFSWRK/ocards
cp $NWPETC/namlsts_tco383 $GFSWRK/namlsts_tco383

#  cd $GFSWRK

# cd ${GFSDIR}

 #FCT_MODEL=$MDIR/src_fpvs/$EXEC
 FCT_MODEL=$MDIR/src/$EXEC
 /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} 

 if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
 fi

