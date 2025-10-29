#!/bin/ksh

#-- enviornment
 user=`whoami`
# datamv='login11'
 if [ ${machine} = a100 ]; then
         mach='x86_64'
 elif [ ${machine} = fx1000 ]; then
         mach=${machine}
 fi
 dmsdb_home=$(cat ~/.dmsrc |xargs | cut -d' ' -f 2)
 DMSPATH=/package/${mach}/dms/dms.v4/bin
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
 fgdtg=$(/users/xb80/bin/Caldtg.ksh ${dtg} -6)

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

 ${DMSPATH}/rdmsdbcrt -p ufs $idmsdb
 ${DMSPATH}/rdmscrt $idmsfile

  export LNCP='ln -fs'

  # TCo IC data path
  export source="/data/common/gfs/dms_data/ncep_ana_n1.ufs/TCo${JCAP}l72_${dtg}"
  #export source="/data/common/gfs/dms_data/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"
  #export source="/data/common/gfs/dms_data/ERA5_n1.ufs/TCo${JCAP}l72_${dtg}"
  #export source="/data/common/gfs/dms_data/NCEPCFSR_n1.ufs/TCo${JCAP}l72_${dtg}"

  # link/copy DMS files
  export target="${dmsdb_home}/${idmsdb}.ufs"

  # analysis
  echo ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
  echo ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}

 ${DMSPATH}/rdmsdbcrt -p ufs bckdms

 export source="/data/common/gfs/dms_data/bckdms.ufs"
 export target="${dmsdb_home}/bckdms.ufs"

 if [ ! -e ${target}/${bckhead} ] ; then
   ${DMSPATH}/rdmscrt ${bckhead}@bckdms
   ${LNCP} ${source}/${bckhead}/* ${target}/${bckhead}
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
export GLB_TYPHINI="/nwpr/gfs/a361/MODEL/typhoon"
export FIXDIR=${GFSFIX}

export ANADMS=${idmsfile}
export FCSTDMS=${odmsfile}
export BCKOPS=${bckhead}@bckdms

${DMSPATH}/rdmspurge -f FCSTDMS
${DMSPATH}/rdmscrt -l34 FCSTDMS

export FLIB_CNTL_BARRIER_ERR=FALSE
#export O3FORC=${O3FORC:-${FIXDIR}/global_o3prdlos.f77}
export O3FORC=${O3FORC:-${FIXDIR}/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77}
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
if [ $GITLAB_CICD = 1 ] ; then
  echo -e "00\n06" > $GFSWRK/gfsctl
fi
cp $NWPETC/ocards $GFSWRK/ocards
cp $NWPETC/namlsts $GFSWRK/namlsts

if [ $JCAP = 639  ] ; then
  MODLST_RES='dt=450., hfilt=1., cgw=4.2e-5, cgwd=1.20, cmbk=1.00,spl2=50.,itter=2,'
  MODLST_PHY='isot=2,ivegsrc=2,'
  MODEL_BASIC='nco=640,'
elif [ $JCAP = 383  ] ; then
  MODLST_RES='dt=600., hfilt=1., cgw=2.6e-5, cgwd=2.40, cmbk=0.60,'
  MODLST_PHY="nmgwcv=1,"
  MODEL_BASIC='nco=384,'
elif [ $JCAP = 199  ] ; then
  MODLST_RES='dt=1200., hfilt=1., cgwd=2.40, cmbk=0.60, mom4ice=f, '
  MODEL_BASIC='nco=200,'
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
  taui=0.0, taue=120.0, tauo=1.0, taup=6.0, taureg=0.,
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
  outdms   =0, outgrb2  =${OUTGRB}, 
  ndsladvh2=false,
  isot=1, ivegsrc=1, cgwd=1.20, cmbk=1.00, tofd=t,
  spl1=5., spl2=100.,
  two_loop=t, mass_dp=t, doclx=f,
  itter=1, factop=80.,alpha=0.7,
!  naero=1,
  ${MODLST_RES}
  ${MODLST_PHY}
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
  use_cpm=false, use_declination=false
 /

EOF

 if [ $CMAKE_BUILD = 1 ] ; then
	FCT_MODEL=$MDIR/build_${machine}/bin/tcogfs.x
 else
	FCT_MODEL=$MDIR/src/$EXEC
 fi

 if [ ${machine} = a100 ]; then
  export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC="1"
  export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC_POOLSIZE="40G"
  export NVCOMPILER_ACC_USE_GRAPH="1"
  export NVCOMPILER_ACC_CUDA_NOCOPY="1"
 fi
 /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} -Wl,-T

 if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!" ; exit 9
 fi
