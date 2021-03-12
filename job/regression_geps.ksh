#!/bin/ksh

#-- enviornment
 user=`whoami`
 datamv='login11'
 dmsdb_home=$(cat ~/.dmsrc |xargs | cut -d' ' -f 2)
 DMSPATH="/package/${machine}/dms/dms.v4/bin"
 GFSDIR="$MDIR"
 GFSFIX="$MDIR/fix"
 GFSWRK="${GFSDIR}/work_${machine}"
 Caldtg="/nwpr/gfs/xb80/bin/Caldtg.ksh"

 # GODAS data source
 dataDir=/nwpr/gfs/xb99/data2
 # GODAS data
 GODASDIR=${dataDir}/obs_data/godas
 # GODAS data (Input)
 GODASpenDIR=${dataDir}/obs_data/godas/pentad
 # SIT initial data
 GODASgfsDIR="/nwpr/gfs/xb80/data2/IC_SIT"

 export OCNDMS=${OCNDMS:-"OCNDMS@ncepec"}

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
 dtg10="20${dtg}"

 fgdtg=$(${Caldtg} ${dtg} -6)

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

 ${DMSPATH}/rdmsdbcrt -p ufs $idmsdb
 ${DMSPATH}/rdmscrt $idmsfile

  export LNCP='ln -fs'

  # maybe no need to change
#  export source="/data/common/gfs/dms_data/TCo639L72_S2TY.ufs/T_exp20${dtg}"           # TCo IC data path
  export source="/nwpr/gfs/xb126/data2/Tool/Nemsio2Dms_v2/OUTPUT/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"           # TCo IC data path

  # link/copy DMS files
  export target="${dmsdb_home}/${idmsdb}.ufs"

  # analysis
  echo ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${dtg}* ${target}/${idmshead}${idmsbody}${idmstail}
  echo ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}
       ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}

  ${DMSPATH}/rdmsdbcrt -p ufs bckdms
  ${DMSPATH}/rdmscrt BCK_TCo${JCAP}_${DMSFLAG}30S@bckdms
 
  # Bundary condition
  export source="/data/common/gfs/dms_data/bckdms.ufs"
  export target="${dmsdb_home}/bckdms.ufs"
  
  if [ ! -e ${target}/BCK_TCo${JCAP}_${DMSFLAG}30S ] ; then
    ${LNCP} ${source}/BCK_TCo${JCAP}_${DMSFLAG}30S/* ${target}/BCK_TCo${JCAP}_${DMSFLAG}30S
  fi
 
  # Oceanic initial and bundary condition 
  dtg_yy=`echo $dtg10 | cut -c1-4`
  dtg_yymmdd=`echo $dtg10 | cut -c1-8`
 
  dtg_yynow=${dtg_yy}
  dtg_yybefore1=$(( $dtg_yynow -1 ))
  dtg_yyafter1=$(( $dtg_yynow +1 ))
  dtg_yyafter2=$(( $dtg_yynow +2 ))
 
  godas_mmdd=( 0105 0110 0115 0120 0125 0130
               0204 0209 0214 0219 0224
               0301 0306 0311 0316 0321 0326 0331
               0405 0410 0415 0420 0425 0430
               0505 0510 0515 0520 0525 0530
               0604 0609 0614 0619 0624 0629
               0704 0709 0714 0719 0724 0729
               0803 0808 0813 0818 0823 0823
               0902 0907 0912 0917 0922 0927
               1002 1007 1012 1017 1022 1027
               1101 1106 1111 1116 1121 1126
               1201 1206 1211 1216 1221 1226 1231 )
 
  pre480=`${Caldtg} ${dtg10} -480`
  pre480_yy=`echo $pre480 | cut -c1-4`
  pre480_mmdd=`echo $pre480 | cut -c5-8`
 
  if [ ${pre480_mmdd} -le ${godas_mmdd[0]} ]; then
     godas_base_yy=$(( ${pre480_yy}-1 ))
     godas_base_mmdd=${godas_mmdd[${ngodasmmdd}-1]}
     godas_base=${godas_base_yy}${godas_base_mmdd}
     echo godas_base=${godas_base}
  else
     godas_base_yy=${pre480_yy}
     godas_base_mmdd=${godas_mmdd[0]}
     godas_base=${godas_base_yy}${godas_base_mmdd}
     godas_base_pre1=${godas_base_yy}${godas_base_mmdd}
     godas_base_pre2=${godas_base_yy}${godas_base_mmdd}
 
    for iday in ${godas_mmdd[@]} ; do
      if [ $iday -le ${pre480_mmdd} ]; then
        godas_base_pre2=${godas_base_pre1}
        godas_base_pre1=${godas_base}
        godas_base=${pre480_yy}${iday}
      else
        break
      fi
    done
  fi
 
 if [ -f ${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base}.lev21mixed.a.nc ]; then
   echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base}.lev21mixed.a.nc found."
 else
   echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base}.lev21mixed.a.nc doesn't exist, check pre1 data."
   if [ -f ${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre1}.lev21mixed.a.nc ]; then
     echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre1}.lev21mixed.a.nc found."
     godas_base=${godas_base_pre1}
   else
     echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre1}.lev21mixed.a.nc doesn't exist, check pre2 data."
     if [ -f ${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre2}.lev21mixed.a.nc ]; then
       echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre2}.lev21mixed.a.nc found."
       godas_base=${godas_base_pre2}
     else
       echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre2}.lev21mixed.a.nc doesn't exist, exit"
       echo "01_check_godas failed." > $stafl
       exit 1
     fi
   fi
 fi

  yyyy_now=`echo $dtg10 | cut -c1-4`
  yyyymmdd=`echo $dtg10 | cut -c1-8`
  yyyy_m1=$(($yyyy_now-1))
  yyyy_p1=$(($yyyy_now+1))
  yyyy_p2=$(($yyyy_now+2))
 
  echo godas_base=${godas_base}

  LNCP="ln -fs" 
  rm -f ${GFSWRK}/dailygodas${yyyy_m1}
  ${LNCP} ${GODASgfsDIR}/godas.P.${yyyy_m1}.TCo383.a.nc           ${GFSWRK}/dailygodas${yyyy_m1}
  rm -f ${GFSWRK}/dailygodas${yyyy_now}
  ${LNCP} ${GODASgfsDIR}/godas.P.${godas_base}.pano.TCo383.a.nc   ${GFSWRK}/dailygodas${yyyy_now}
  rm -f ${GFSWRK}/dailygodas${yyyy_p1}
  ${LNCP} ${GODASgfsDIR}/godas.P.${godas_base}.pano.TCo383.a.nc   ${GFSWRK}/dailygodas${yyyy_p1}
  rm -f ${GFSWRK}/dailygodas${yyyy_p2}
  ${LNCP} ${GODASgfsDIR}/godas.P.${yyyy_p2}.TCo383.a.nc           ${GFSWRK}/dailygodas${yyyy_p2}
  rm -f ${GFSWRK}/unit.97
  ${LNCP} ${GODASgfsDIR}/TCo383_${godas_base}pano_${yyyymmdd}.nc  ${GFSWRK}/unit.97

  # get oceanic climatology
  OCNCLMFN="OCNCLM"
  ocnclm_source="/nwpr/gfs/xb99/data2/dmsdb/OISST_clim.ufs"
  if [ ! -e ${dmsdb_home}/${idmsdb}.ufs/${OCNCLMFN}/*  ] ; then
    ${DMSPATH}/rdmscrt -l34  ${OCNCLMFN}@${idmsdb}
    ln -fs ${ocnclm_source}/CLMdaily_TCo383_2010_2019/* ${dmsdb_home}/${idmsdb}.ufs/${OCNCLMFN}/
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
 ifilin_sst='OCNDMS',
 ifilin_nc='${GFSWRK}',
 ifilin_ClmANA='OCNCLM',
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

${DMSPATH}/rdmspurge -f FCSTDMS
${DMSPATH}/rdmscrt -l34 FCSTDMS

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
  MODLST_RES='dt=225., tfilt=0.040, hfilt=1., cgw=4.2e-5,'
  MODEL_BASIC='nco=640,'
elif [ $JCAP = 383  ] ; then
  MODLST_RES='dt=360., tfilt=0.050, hfilt=1, cgw=2.6e-5,'
  MODEL_BASIC='nco=384,'
fi

 export MODLST_SIT="do_sit=false, fsit=-99., ldailyFCTsst=true, ldailyFCTicesndpt=f, dailyClm_option=1, dSITdt_intv=12., weightSIT=1.0,"

 export SIT_NML=" lpre6hr_sit=f, loutsit24=t, outsitmean=6.,
 sit_domain_w= 0., sit_domain_e= 360., sit_domain_s= -30., sit_domain_n= 30., sit_domain_extgrd= 10.,
 ltimeblending= T,
 timebl_option=0, timebl_start= 5., timebl_allsit= 15.,
 lobs_ocn_rerun= F,
 TRIGSIT%COUNTER= 1, TRIGSIT%UNIT= steps, TRIGSIT%ADJUSTMENT= exact, TRIGSIT%OFFSET= 0,
 lsit_ice= T, lsit_salt= T,
 zocn_option= 99,
 ocn_tlz= 200., ocn_k1= 20,
 lssst= T, sit_ice_option= 0,
 maskid= 1, lgodas= T, ldailysst=T, lmixedlayer=F, lamip= T, lwoa0= T, lsftobswt= T, lwarning_msg= 2, lsice_nudg= F, lsit_lw= F,
 ssit_restore_time= -99.,
 usit_restore_time= 604800.,
 dsit_restore_time= 0.,
 ssits_restore_time= -99.,
 usits_restore_time= 604800.,
 dsits_restore_time= 0.,
 ssituv_restore_time= -99.,
 usituv_restore_time= 604800.,
 dsituv_restore_time= 0.,
 nsit_nudg= 0,
 sitbox_nudg_w(1)=  40., sitbox_nudg_e(1)= 180., sitbox_nudg_s(1)= -15., sitbox_nudg_n(1)=  15.,
 sitbox_st_restore_time(1)= -99.,
 sitbox_ut_restore_time(1)= -99.,
 sitbox_dt_restore_time(1)= 0.,
 sitbox_ss_restore_time(1)= -99.,
 sitbox_us_restore_time(1)= -99.,
 sitbox_ds_restore_time(1)= 0.,
 sitbox_suv_restore_time(1)= -99.,
 sitbox_uuv_restore_time(1)= 0.,
 sitbox_duv_restore_time(1)= 0.,
 locaf= F,locn= F,lopen_bound= F,lall_straits= T,lstrict_channel= T,etopo_nres= 1,
 ocn_domain_w= 0.,ocn_domain_e= 360.,ocn_domain_s= -80.,ocn_domain_n= 80.,
 ratio_dt_o2a= 1.,ocn_couple_option= 0,high_current_killer= 4,
 locn_msg= F,ocn_lon_factor= 1,ocn_lat_factor= 1,
 TRIGOCN%COUNTER= 1, TRIGOCN%UNIT= steps, TRIGOCN%ADJUSTMENT= exact, TRIGOCN%OFFSET= 0,
 socn_restore_time= -9.000000000000000E+033,
 uocn_restore_time= -9.000000000000000E+033,
 docn_restore_time= -9.000000000000000E+033,
 nobox_nudg= 0,
 obox_restore_time= -9.000000000000000E+033, obox_nudg_flag= 0, obox_nudg_w= 6*-999.000000000000, obox_nudg_e= 6*-999.000000000000, obox_nudg_s= 6*-999.000000000000, obox_nudg_n= 6*-999.000000000000,
 kocn_dm0z= 1., ncarpet= 1, kcsmag= 1.,kalbw= 1.,ck=0.1,ce=0.7,Prw= 1.,d0= 0.03,csl= -27.,
 por_min= 0.1,csiced= 0.,lasia= F,lsteady_TKE=F,
"


cat > ${GFSWRK}/namlsts << EOF
 &model_param
  nco=640,
  lev=72,
  ncld=3,
  octahedral=true,
  nout=9000,
  io_quilting=false,
  npex=${NPEX},
  npey=${NPEY},
  ${MODEL_BASIC}
 &end

 &modlst
  taui=0.0, taue=120.0, tauo=1.0, taup=6.0, taureg=6.,
  dt=225.0,
  cstar=f, update=t, lsimpl=t,
  tfilt=0.04, hfilt=1.,
  ksgeo=2, yesdia=t,
  dopbl=t, docup=t, dorad=t, dolsp=t, doshl=t, dodry=f, 
  dograv=true, docgrav=true,
  donnmi=true, 
  dosppt=true,  dospptout=false,
  doshum=false, 
  dossst=true,
  cutfreq=3, nnmivm=3,doincr=f,
  hdiff=t, frad=1.0, ldiag=0,
  idg=40, jdg=108,
  itypbl=0, numreduce=5, ptmeans=800., ptop=0.1,
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
  spl1=5.,
  ${MODLST_RES}
  ${MODLST_SIT}
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

 &sit_nml
  ${SIT_NML}
 /


EOF


 FCT_MODEL=$MDIR/src/$EXEC
 /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} 

 if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
 fi

