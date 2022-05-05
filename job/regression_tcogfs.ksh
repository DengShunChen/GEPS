#!/bin/ksh

PJM_JOBID=${PJM_JOBID:-???}

#-- enviornment
 user=`whoami`
 datamv='h6dm13'
 dmsdb_home=$(cat ~/.dmsrc |xargs | cut -d' ' -f 2)
 DMSPATH=/package/${machine}/dms/bin
 GFSDIR=$MDIR
 GFSFIX=$MDIR/fix
 GFSWRK=${GFSDIR}/work.${PJM_JOBID}
 rm -rf $GFSWRK
 mkdir -p $GFSWRK

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

 odmshead="J${PJM_JOBID}_"
 odmsbody=${dtg}
 odmstail="${DMSFLAG}MG"
 odmsdb=${idmsdb}

#---------------------------------------------------------#
  idmsfile=${idmshead}${idmsbody}${idmstail}@${idmsdb}
  odmsfile=${odmshead}${odmsbody}${odmstail}@${odmsdb}

  ${DMSPATH}/rdmsdbcrt -p ufs $idmsdb
  ${DMSPATH}/rdmscrt $idmsfile
  if [ $? -ne 0 ] ; then 
    mkdir /users/xb80/data/dmsdb/$(echo ${idmsfile} | cut -d'@' -f2).ufs/$(echo ${idmsfile} | cut -d'@' -f1)
  fi

  export LNCP='ln -fs'

  # maybe no need to change
  # export source="/users/xb80/data/dmsdb/TCo${JCAP}L72_S2TY.ufs/T_exp20${dtg}"           # TCo IC data path
   export source="/users/xb80/data/dmsdb/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"           # TCo IC data path
  #export source="/nwpr/gfs/xb126/data2/Tool/Nemsio2Dms_v2/OUTPUT/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"           # TCo IC data path

  # link/copy DMS files
  export target="${dmsdb_home}/${idmsdb}.ufs"

  # analysis
  echo ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
  echo ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}

# ${DMSPATH}/rdmsdbcrt -p ufs bckdms
# ${DMSPATH}/rdmscrt BCK_TCo${JCAP}_${DMSFLAG}30S@bckdms
# if [ $? -ne 0 ] ; then 
#   mkdir -p /users/xb80/data/dmsdb/bckdms.ufs/BCK_TCo${JCAP}_${DMSFLAG}30S
# fi

# export source="/users/xb80/data/common/gfs/dms_data/bckdms.ufs"
# export target="${dmsdb_home}/bckdms.ufs"
#
# if [ ! -e ${target}/BCK_TCo${JCAP}_${DMSFLAG}30S ] ; then
#   ${LNCP} ${source}/BCK_TCo${JCAP}_${DMSFLAG}30S/* ${target}/BCK_TCo${JCAP}_${DMSFLAG}30S
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
 &end
EOF

export GFSDIR DMSPATH
export NWPETC=${GFSDIR}/etc
export NWPETCGLB=${GFSWRK}
export GLB_TYPHINI="/nwp/npcagfs/TYP/M00/dtg/ty"
export FIXDIR=${GFSFIX}

export ANADMS=${idmsfile}
export FCSTDMS=${odmsfile}
export BCKOPS=BCK_TCo${JCAP}_${DMSFLAG}30S@bckdms
#export BCKOPS=BCK_TCo${JCAP}_${DMSFLAG}30S_xnew@bckdms

${DMSPATH}/rdmspurge -f FCSTDMS
${DMSPATH}/rdmscrt -l34 FCSTDMS
# For fx1000 temp fix
if [ $? -ne 0 ] ; then 
  mkdir /users/xb80/data/dmsdb/$(echo ${odmsfile} | cut -d'@' -f2).ufs/$(echo ${odmsfile} | cut -d'@' -f1)
fi

export FLIB_CNTL_BARRIER_ERR=FALSE
export O3FORC=${O3FORC:-${FIXDIR}/global_o3prdlos.f77}
export O3CLIM=${O3CLIM:-${FIXDIR}/global_o3clim.txt}
export AEROSOL_FILE=${AEROSOL_FILE:-${FIXDIR}/global_climaeropac_global.txt}
export EMMISSIVITY_FILE=${EMMISSIVITY_FILE:-${FIXDIR}/global_sfc_emissivity_idx.txt}


cd $GFSWRK
#====================================================================
ln -fs $O3FORC fort.28
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
  ncld=3,
  octahedral=t,
  nout=9000,
  io_quilting=f,
  npex=${NPEX},
  npey=${NPEY},
  ${MODEL_BASIC}
 &end

 &modlst
  taui=0.0, taue=120.0, tauo=1.0, taup=6.0, taureg=6.,
  dt=450.0,
  cstar=f, update=t, lsimpl=t,
  hfilt=1.,
  ksgeo=2, yesdia=t,
  dopbl=t, docup=t, dorad=t, dolsp=t, doshl=t, dodry=f, dograv=t, docgrav=t,
  donnmi=t, 
  dosppt=false, dospptout=false, doshum=false,
  cutfreq=3, nnmivm=3,
  doincr=f,
  hdiff=t, frad=1.0,
  ldiag=0,
  idg=40, jdg=108,
  itypbl=0, numreduce=5, ptmeans=800.,
  nmcup=6, nmpbl=4, nmland=2, nmshl=3,
  nmgwor=2, nmgwcv=2,
  ktcup=20, cgw=4.2e-5,
  mtnvar=14, doo3l=t,
  irad=2, ioutsigr=1,
  ggdef='${DMSFLAG}0G', gmdef='${DMSFLAG}MG',
  domfc=384., out_green=t, otgreen=3., out_hp=f,
  ndsladvh2=f,
  isot=1, ivegsrc=1, cgwd=1.20, cmbk=1.00,
  spl1=5., spl2=50., af=0.1,
  ${MODLST_RES}
 &end

 &typ
  write_mem=0,
  trk_intv=3,
  write_tau=6,
 &end
 
 &stochy_physics
  ncep_seeds = false,
  sppt = 0.8,0.4,0.2,0.08,0.04
  sppt_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7 
  sppt_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
  shum = 0.8,-999,-999,-999,-999
  shum_decort = 2.16E4,1.728E5,2.592E6,7.776E6,3.1536E7
  shum_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 /
EOF


 FCT_MODEL="$MDIR/build/bin/tcogfs.x "
 #FCT_MODEL="$MDIR/src/MTCo639L72_fx1000 "
 # Version not compatiable -->FCT_MODEL="/users/fjapl10/gfs/GFS/MNH/GJ_gfs/src_tco639l72/src/MTCo639L72_fx1000"

 /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} 

 if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
 fi

