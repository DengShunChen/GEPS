MDIR=$PWD

machine=a100
TERR=GMTED30
JCAP=383

NPEX=1
NPEY=8

MPI=8
OMP=1
export OMP_NUM_THREADS=$OMP

# << Enviroment variables of OpenACC with NVCOMPILER >>
export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC="1"
export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC_POOLSIZE="40G"
export NVCOMPILER_ACC_USE_GRAPH="1"
export NVCOMPILER_ACC_CUDA_NOCOPY="1"

. /usr/share/Modules/init/bash

module purge
module use ${MDIR}/modulefiles
module load modulefile.tcogfs.a100
module unuse ${MDIR}/modulefiles

#-- enviornment
 user=`whoami`
# datamv='login11'
 dmsdb_home=$(cat ~/.dmsrc |xargs | cut -d' ' -f 2)
 DMSPATH=/package/${machine}/dms/dms.v4/bin
 GFSDIR=$MDIR
 GFSFIX=$MDIR/fix
 GFSWRK=$MDIR/build_${machine}/test

#-- dms data
 JCAP=${JCAP:-639}

 if [ $JCAP = 639  ] ; then
   DMSFLAG=GJ
 elif [ $JCAP = 383  ] ; then
   DMSFLAG=GI
 fi

 dtg='18090800'
 fgdtg=$(/users/xb80/bin/Caldtg.ksh ${dtg} -6)

 idmshead='MASOPS'
 idmsbody=''
 idmstail=''
 idmsdb="TCo${JCAP}L72"

 odmshead='STOC'
 odmsbody=${dtg}
 odmstail="${DMSFLAG}MG"
 odmsdb=${idmsdb}

#-- executable
 EXEC='MTCo639L72_'${machine}

#---------------------------------------------------------#
 idmsfile=${idmshead}${idmsbody}${idmstail}@${idmsdb}
 odmsfile=${odmshead}${odmsbody}${odmstail}@${odmsdb}

  export LNCP='ln -fs'

  # TCo IC data path
  export source="/data/common/gfs/dms_data/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"

  # link/copy DMS files
  export target="${dmsdb_home}/${idmsdb}.ufs"

  # analysis
  echo ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
  echo ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}

 export source="/data/common/gfs/dms_data/bckdms.ufs"
 export target="${dmsdb_home}/bckdms.ufs"

 if [ ! -e ${target}/BCK_TCo${JCAP}_${DMSFLAG}30S ] ; then
   ${DMSPATH}/rdmscrt BCK_TCo${JCAP}_${DMSFLAG}30S@bckdms
   ${LNCP} ${source}/BCK_TCo${JCAP}_${DMSFLAG}30S/* ${target}/BCK_TCo${JCAP}_${DMSFLAG}30S
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
 &end
EOF

export GFSDIR DMSPATH
export NWPETC=${GFSDIR}/etc
export NWPETCGLB=${GFSWRK}
export GLB_TYPHINI="/nwpr/gfs/a361/MODEL/typhoon"
export FIXDIR=${GFSFIX}

export ANADMS=${idmsfile}
export FCSTDMS=${odmsfile}
export BCKOPS=BCK_TCo${JCAP}_${DMSFLAG}30S@bckdms

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
cp $NWPETC/ocards $GFSWRK/ocards
cp $NWPETC/namlsts $GFSWRK/namlsts

if [ $JCAP = 639  ] ; then
  MODLST_RES='dt=450., hfilt=1., cgw=4.2e-5, cgwd=1.20, cmbk=1.00,'
  MODEL_BASIC='nco=640,'
elif [ $JCAP = 383  ] ; then
  MODLST_RES='dt=720., hfilt=1., cgw=2.6e-5, cgwd=1.60, cmbk=0.30,'
  MODEL_BASIC='nco=384,'
fi

cat > ${GFSWRK}/namlsts << EOF
 &model_param
  nco=640,
  lev=72,
  ncld=7,
  octahedral=true,
  nout=9000,
  io_quilting=false,
  npex=${NPEX},
  npey=${NPEY},
  ${MODEL_BASIC}
 &end

 &modlst
  taui=0.0, taue=120.0, tauo=1.0, taup=6.0, taureg=6.,
  dt=450.0,
  cstar=f, update=t, lsimpl=t,
  hfilt=1.,
  ksgeo=1, yesdia=f,
  dopbl=t, docup=t, dorad=t, dolsp=t, doshl=t, dodry=f,
  dograv=true, docgrav=true,
  donnmi=true,
  dosppt=false, dospptout=false,
  doshum=false,
  cutfreq=3, nnmivm=3,
  doincr=f,
  hdiff=t, frad=1.0, ldiag=0,
  idg=40, jdg=108,
  itypbl=0, numreduce=5, ptmeans=800.,
  irad=2, nmland=2,
  nmcup=6, nmshl=3, nmpbl=4, nmmiph=2,
  nmgwor=2, nmgwcv=2,
  ktcup=20, cgw=4.2e-5,
  mtnvar=14, doo3l=t,
  ioutsigr=1,
  ggdef='${DMSFLAG}0G', gmdef='${DMSFLAG}MG',
  domfc=384., out_green=t, otgreen=3., out_hp=false,
  ndsladvh2=false,
  isot=1, ivegsrc=1, cgwd=1.20, cmbk=1.00,
  spl1=10., spl2=100.,
  two_loop=t,
  itter=2, vd=0.002, factop=60.,
  ${MODLST_RES}
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
  ssst = 0.80,-999,-999,-999,-999
  ssst_seed = -999,-999,-999,-999,-999
  ssst_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
  ssst_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 /

EOF
