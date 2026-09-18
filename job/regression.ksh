#!/bin/ksh

#-- enviornment
 user=`whoami`
# datamv='login11'
 rocm_lab=0
 if [ ${machine} = a100 ]; then
         mach='x86_64'
 elif [ ${machine} = fx1000 ]; then
         mach=${machine}
 elif [ ${machine} = rocm -o ${machine} = rocm_cpu ]; then
         mach='rocm'
         rocm_lab=1
 fi
 if [ ${rocm_lab} = 1 ]; then
   dmsdb_home=${GEPS_DMSDB:-/mlsteam/workspace/data/geps/dmsdb}
   if [ ! -f ${HOME}/.dmsrc ]; then
     printf '[DMSPATH]\n%s\n' "${dmsdb_home}" > ${HOME}/.dmsrc
   fi
   : ${GEPS_LIB_ROOT:=${MDIR}/../GEPS_LIB}
   DMSPATH=${GEPS_LIB_ROOT}/install/rocm/dms-4.0.0/bin
 else
   dmsdb_home=$(cat ~/.dmsrc |xargs | cut -d' ' -f 2)
   DMSPATH=/package/${mach}/dms/dms.v4/bin
 fi
 GFSDIR=$MDIR
 GFSFIX=$MDIR/fix
 export GFSWRK=${GFSDIR}/work_${machine}
 rm -rf $GFSWRK
 mkdir -p $GFSWRK

 #about grib2 output
 # grib2 output folder
 GRBDIR="."
 #output grib2 0:off 1:on
 OUTGRB="1"
 OUTDMS=0
 if [ ${rocm_lab} = 1 ]; then
   # wrt_grb2/addgrid SIGSEGV at tau=0 (imax=jmax=-1). Keep DMS output for smoke.
   OUTGRB="0"
   OUTDMS=1
 fi
 
#-- dms data
 JCAP=${JCAP:-639}

 if [ $JCAP = 639  ] ; then
   DMSFLAG=GJ
 elif [ $JCAP = 383  ] ; then
   DMSFLAG=GI
 elif [ $JCAP = 199  ] ; then
   DMSFLAG=GK
 fi

 dtg='22081500'
 if [ -x /users/xb80/bin/Caldtg.ksh ]; then
   fgdtg=$(/users/xb80/bin/Caldtg.ksh ${dtg} -6)
 else
   fgdtg='22081418'
 fi

 idmshead='MASOPS'
 idmsbody=''
 idmstail=''
 idmsdb="TCo${JCAP}L72"

 odmshead='STOC'
 odmsbody=${dtg}
 odmstail="${DMSFLAG}MG"
 odmsdb=${idmsdb}

# bckhead="BCK_TCo${JCAP}_${DMSFLAG}30S_dyclm"
 bckhead="BCK_TCo${JCAP}_${DMSFLAG}30S_xnew"

 ksgeo=1
#-- executable
 EXEC='MTCo639L72_'${machine}

#---------------------------------------------------------#
 idmsfile=${idmshead}${idmsbody}${idmstail}@${idmsdb}
 odmsfile=${odmshead}${odmsbody}${odmstail}@${odmsdb}

 if [ ${rocm_lab} = 1 ]; then
   if [ ! -d ${dmsdb_home}/${idmsdb}.ufs ]; then
     ${DMSPATH}/rdmsdbcrt -p ufs $idmsdb
   fi
   if [ ! -d ${dmsdb_home}/${idmsdb}.ufs/${idmshead}${idmsbody}${idmstail} ]; then
     ${DMSPATH}/rdmscrt $idmsfile
   fi
 else
   ${DMSPATH}/rdmsdbcrt -p ufs $idmsdb
   ${DMSPATH}/rdmscrt $idmsfile
 fi

  export LNCP='ln -fs'

  # TCo IC data path
  if [ ${rocm_lab} = 1 ]; then
    export source="${dmsdb_home}/ncep_ana_n1.ufs/TCo${JCAP}l72_${dtg}"
  else
    export source="/data/common/gfs/dms_data/ncep_ana_n1.ufs/TCo${JCAP}l72_${dtg}"
  fi
  #export source="/data/common/gfs/dms_data/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"
  #export source="/data/common/gfs/dms_data/ERA5_n1.ufs/TCo${JCAP}l72_${dtg}"
  #export source="/data/common/gfs/dms_data/NCEPCFSR_n1.ufs/TCo${JCAP}l72_${dtg}"

  if [ ! -d "${source}" ]; then
    echo "ERROR: IC source missing: ${source}"
    exit 1
  fi

  # link/copy DMS files
  export target="${dmsdb_home}/${idmsdb}.ufs"

  # analysis
  echo ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
  echo ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}

 if [ ${rocm_lab} = 1 ]; then
   if [ ! -d ${dmsdb_home}/bckdms.ufs ]; then
     ${DMSPATH}/rdmsdbcrt -p ufs bckdms
   fi
   export source="${dmsdb_home}/bckdms.ufs"
 else
   ${DMSPATH}/rdmsdbcrt -p ufs bckdms
   export source="/data/common/gfs/dms_data/bckdms.ufs"
 fi
 export target="${dmsdb_home}/bckdms.ufs"

 if [ ! -e ${target}/${bckhead} ] ; then
   ${DMSPATH}/rdmscrt ${bckhead}@bckdms
   ${LNCP} ${source}/${bckhead}/* ${target}/${bckhead}
 fi
 if [ ${rocm_lab} = 1 -a ! -e ${target}/${bckhead} ]; then
   echo "ERROR: BCK missing: ${target}/${bckhead}"
   exit 1
 fi

#-- MERRA2 aerosol climatological data :
#     for AERO_TCo383L72  : ksgeo=1 , TER=30S_xnew
#     for AERO_TCo199L128 : ksgeo=2 , TER=30S_xnew
 aerohead="AERO_TCo${JCAP}L72"
 export AERODMS=${aerohead}@aeroclx
 export source="/data/common/gfs/dms_data/aeroclx.ufs"
 export target="${dmsdb_home}/aeroclx.ufs"
# ${DMSPATH}/rdmsdbcrt -p ufs aeroclx
# if [ ! -e ${target}/${aerohead} ] ; then
#   ${DMSPATH}/rdmscrt ${AERODMS}
#   ${LNCP} ${source}/${aerohead}/* ${target}/${aerohead}
# fi

#----------------------------------------------------------------#
RSMMPMD='f'
MODLST_RSM=''
if [[ ${RSMMPMD} = 't' ]] ;then
  MODLST_RSM='outrsm=t,rsmoutinv=6,'
fi
#----------------------------------------------------------------#

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
 namlsts='${GFSWRK}/namlsts',
 crdate='${GFSWRK}/crdate',
 ocards='${GFSWRK}/ocards',
 cntrl='${GFSWRK}/gfsctl',
! ifilin_aero='AERODMS',
 ifilin_nc='${GFSWRK}',
 ifilout_grb='${GRBDIR}', 
 &end
EOF

export GFSDIR DMSPATH
export NWPETC=${GFSDIR}/etc
export NWPETCGLB=${GFSWRK}
#export GLB_TYPHINI="/ncs/ncsatyp/TYP/M00/dtg/tdty" #inner
export GLB_TYPHINI="/nwpr/gfs/a361/MODEL/typhoon"
export FIXDIR=${GFSFIX}

export ANADMS=${idmsfile}
export FCSTDMS=${odmsfile}
export BCKOPS=${bckhead}@bckdms

if [ ${rocm_lab} = 1 ]; then
  if [ -d ${dmsdb_home}/${odmsdb}.ufs/${odmshead}${odmsbody}${odmstail} ]; then
    ${DMSPATH}/rdmspurge -f FCSTDMS || true
  fi
else
  ${DMSPATH}/rdmspurge -f FCSTDMS
fi
${DMSPATH}/rdmscrt -l34 FCSTDMS

export FLIB_CNTL_BARRIER_ERR=FALSE
#export O3FORC=${O3FORC:-${FIXDIR}/global_o3prdlos.f77}
if [ -f ${FIXDIR}/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77 ]; then
  export O3FORC=${O3FORC:-${FIXDIR}/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77}
else
  export O3FORC=${O3FORC:-${FIXDIR}/global_o3prdlos.f77}
fi
export O3CLIM=${O3CLIM:-${FIXDIR}/global_o3clim.txt}
export AEROSOL_FILE=${AEROSOL_FILE:-${FIXDIR}/global_climaeropac_global.txt}
export EMMISSIVITY_FILE=${EMMISSIVITY_FILE:-${FIXDIR}/global_sfc_emissivity_idx.txt}

cd $GFSWRK
#====================================================================
#ln -fs $O3FORC fort.28
ln -fs $O3FORC global_o3prdlos
ln -fs $O3CLIM fort.48
ln -fs $AEROSOL_FILE  aerosol.dat
ln -fs $EMMISSIVITY_FILE sfc_emissivity_idx.txt

ln -fs $FIXDIR/* .
cp $NWPETC/gfsctl $GFSWRK/gfsctl
if [ "${GITLAB_CICD:-0}" = 1 ] ; then
  echo -e "00\n06" > $GFSWRK/gfsctl
fi
if [ "${SMOKE:-0}" = 1 ] ; then
  echo -e "00\n01" > $GFSWRK/gfsctl
fi
cp $NWPETC/ocards $GFSWRK/ocards
cp $NWPETC/namlsts $GFSWRK/namlsts

if [ $JCAP = 639  ] ; then
  MODLST_RES='dt=450., hfilt=1., cgw=4.2e-5, cgwd=1.20, cmbk=1.00,spl2=50.,itter=2,'
  MODLST_PHY='isot=2,ivegsrc=2,'
  MODEL_BASIC='nco=640,'
  if [ ${machine} = a100 ]; then
    MODLST_PHY="nmgwcv=2, nmmiph=15"
  elif [ ${rocm_lab} = 1 ]; then
    MODLST_PHY="nmgwcv=2, nmmiph=15"
  fi
elif [ $JCAP = 383  ] ; then
  MODLST_RES='dt=600., hfilt=1., cgw=2.6e-5, cgwd=2.40, cmbk=0.60,'
  MODLST_PHY="nmgwcv=1,"
  MODEL_BASIC='nco=384,'
  if [ ${machine} = a100 ]; then
    MODLST_PHY="nmgwcv=2, nmmiph=15"
  elif [ ${rocm_lab} = 1 ]; then
    MODLST_PHY="nmgwcv=2, nmmiph=15"
  fi
elif [ $JCAP = 199  ] ; then
  MODLST_RES='dt=1200., hfilt=1., cgwd=2.40, cmbk=0.60, mom4ice=f, '
  MODEL_BASIC='nco=200,'
  if [ ${machine} = a100 ]; then
    MODLST_PHY="nmgwcv=2, nmmiph=15"
  elif [ ${rocm_lab} = 1 ]; then
    MODLST_PHY="nmgwcv=2, nmmiph=15"
  fi
fi

TAUE=120.0
TAUO=1.0
TAUP=6.0
if [ "${SMOKE:-0}" = 1 ]; then
  # keep taup=6 so getrdy can read the 6h FG (22081418). taup=1 looks for 23Z keys we don't have.
  TAUE=1.0
  TAUO=1.0
  TAUP=6.0
fi

cat > ${GFSWRK}/namlsts << EOF
 &model_param
  nco=640,
  lev=72,
  ncld=7,
  octahedral=true,
  nout=50000,
  io_quilting=false,
  npex=${NPEX},
  npey=${NPEY},
  ${MODEL_BASIC}
 &end

 &modlst
  taui=0.0, taue=${TAUE}, tauo=${TAUO}, taup=${TAUP}, taureg=0.,
  dt=450.0,
  cstar=f, update=t, lsimpl=t,
  hfilt=1.,
  ksgeo=${ksgeo}, yesdia=t,
  dopbl=t, docup=t, dorad=t, dolsp=t, doshl=t, dodry=f,
  dograv=true, docgrav=true,
  donnmi=true,
  dosppt=false, doskeb=false, doshum=false,
  cutfreq=3, nnmivm=3,
  doincr=f,
  hdiff=t, frad=1.0, ldiag=0,
  idg=40, jdg=108,
  itypbl=0, numreduce=5, ptmeans=800.,
  irad=2, nmland=2,
  nmcup=7, nmshl=4, nmpbl=4, nmmiph=12,
  nmgwor=2, nmgwcv=2,
  ktcup=20, cgw=4.2e-5,
  mtnvar=14, doo3l=t,
  ioutsigr=1,
  ggdef='${DMSFLAG}0G', gmdef='${DMSFLAG}MG',
  domfc=0., out_green=t, otgreen=3., out_hp=false,
  outdms   =${OUTDMS:-0}, outgrb2  =${OUTGRB}, 
  ndsladvh2=false,
  isot=1, ivegsrc=1, cgwd=1.20, cmbk=1.00, tofd=t,
  spl1=5., spl2=100.,
  two_loop=t, mass_dp=t, doclx=f,
  itter=1, factop=80.,alpha=0.7,
!  naero=1,
  ${MODLST_RES}
  ${MODLST_PHY}
  ${MODLST_RSM}
 &end

 &typ
  write_mem=0,
  trk_intv=3,
  write_tau=6,
 &end

 &stochy_physics
  ncep_seeds = true,
  use_zmtnblck = true,
  sppt_logit = true,
  sppt_sigtop1 = 0.1,
  sppt_sigtop2 = 0.025,
  sppt_sfclimit = true,
  sppt_sigbot1 = 0.975,
  sppt_sigbot2 = 0.9,
  sppt = 0.80,0.4,0.10,0.08,0.04
  sppt_seed = -999,-999,-999,-999,-999
  sppt_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
  sppt_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
  shum = 0.04,-999,-999,-999,-999
  shum_seed = -999,-999,-999,-999,-999
  shum_decort = 2.16E4,1.728E5,2.592E6,7.776E6,3.1536E7
  shum_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
  shum_sigefold = 0.2,
  skeb_sigtop1 = 0.1,
  skeb_sigtop2 = 0.025, 
  skeb_sigbot1 = 0.975,
  skeb_sigbot2 = 0.9,
  skeb_vdof = 5,
  skebnorm = 2,
  skebfilt = 12,
  skeb = 20.0,-999,-999,-999,-999
  skeb_seed = -999,-999,-999,-999,-999
  skeb_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
  skeb_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
  ssst = 0.80,-999,-999,-999,-999
  ssst_seed = -999,-999,-999,-999,-999
  ssst_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
  ssst_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 /

 &gce_3ice
  SL_sedi=false, sat_predict=true, new_saturation=true,
  use_cpm=false, use_declination=true
 /

EOF

 if [ $CMAKE_BUILD = 1 ] ; then
   if [ ${machine} = rocm ]; then
	# GEPS_PROF_WRAP="rocprofv3 --kernel-trace --stats -d /path -o r_%pid% --" prefixes the binary (2026-09-18)
	FCT_MODEL="${GEPS_PROF_WRAP:-} $MDIR/build_rocm_hip/bin/tcogfs.x"
   elif [ ${machine} = rocm_cpu ]; then
	FCT_MODEL=$MDIR/build_rocm_cpu/bin/tcogfs.x
   else
	FCT_MODEL=$MDIR/build_${machine}/bin/tcogfs.x
   fi
 else
	FCT_MODEL=$MDIR/src/$EXEC
 fi
 if [[ ${RSMMPMD} = 't' ]] ;then
      FCT_RSM="../rsm/run/exe/rsm.x"
      if [[ ! -f ${FCT_RSM}  ]];then echo 'rsm exe not found'; exit -1 ;fi
      export GMPI=${MPI}
      export RMPI=216  #RSM use cpu core
 fi

 if [ ${machine} = a100 ]; then
  export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC="1"
  export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC_POOLSIZE="60G"
  export NVCOMPILER_ACC_USE_GRAPH="1"
  export NVCOMPILER_ACC_CUDA_NOCOPY="1"
 fi
 if [ ${machine} = rocm ]; then
  # 2026-09-10: lab now has TWO PHYSICAL MI300X, both NPS1/SPX (unpartitioned):
  #   GPU[0] = KFD node 4, PCI 0000:46:00.0   GPU[1] = KFD node 5, PCI 0000:66:00.0
  # GPU[1] is a real second card, not a DPX partition (pre-09-10 layout, when the
  # only card was DPX-split and .1 had to be avoided). Each card now exposes the
  # full 192GiB to HIP; the old "~96GB visible" note no longer applies.
  # Still binding GPU0 only: the model runs 1 MPI rank, multi-GPU needs RCCL work.
  # GEPS_ROCM_PCI is informational only - nothing reads it.
  export GEPS_ROCM_PCI="${GEPS_ROCM_PCI:-0000:46:00.0}"
  # 2026-09-16: GEPS_ROCM_DEVICES="0,1" runs on both cards (one MPI rank per
  # card, device_init picks mod(rank, ndevices)); default stays GPU0 only.
  export HIP_VISIBLE_DEVICES="${GEPS_ROCM_DEVICES:-0}"
  export ROCR_VISIBLE_DEVICES="${GEPS_ROCM_DEVICES:-0}"
  export GPU_DEVICE_ORDINAL="${GEPS_ROCM_DEVICES:-0}"
  export OMP_DEFAULT_DEVICE=0
  export OMP_NUM_THREADS=1
  export OMP_TARGET_OFFLOAD=MANDATORY
  export HSA_XNACK=1
  # 2026-09-16: libomptarget's memory manager keeps large freed blocks in a
  # pool and only reuses exact sizes; the radiation blocks all differ in
  # size, so device memory grew monotonically until OUT_OF_RESOURCES
  # (reproduced standalone: exit data of 1.5-1.9 GB never returned VRAM,
  # returned fully with the threshold at 0). Safe now that the layer-6/7/10
  # races are fixed - the earlier "=0 crashes early" was the masked
  # use-after-free.
  export LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD="${LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD:-0}"
  # layer 18 (2026-09-17): runtime-sized private arrays in target regions
  # (sflx_gpu's dimension(nsoil) work arrays) live on the device stack; with
  # the default size neighbouring lanes overwrite each other -> Inf skin
  # temperature on ~10% of land columns, count varying run to run.  A
  # standalone test needs >= 49152 for 20 such arrays; 65536 leaves margin.
  export LIBOMPTARGET_STACK_SIZE="${LIBOMPTARGET_STACK_SIZE:-65536}"
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
 elif [ ${machine} = rocm_cpu ]; then
  unset HIP_VISIBLE_DEVICES ROCR_VISIBLE_DEVICES
  export OMP_TARGET_OFFLOAD=disabled
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
 fi
 if [ ${rocm_lab} = 1 ]; then
  ulimit -s unlimited
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} -Wl,-T
  else
    mpiexec -n $MPI ${FCT_MODEL} -Wl,-T
  fi
 else
  /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} -Wl,-T
 fi

 ###--- NVIDIA Nsight Systems tool for A100 ---###
 #nsys profile -t cuda,nvtx,openacc --cuda-memory-usage=true mpiexec -n ${MPI} ${FCT_MODEL} -Wl,-T

 ###--- run rsm ---###
 #/usr/bin/time -p mpiexec -stdin rfcstparm.all \
 # -n ${GMPI} ${FCT_MODEL} -Wl,-T : -n ${RMPI} ${FCT_RSM}

 if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!" ; exit 9
 fi
