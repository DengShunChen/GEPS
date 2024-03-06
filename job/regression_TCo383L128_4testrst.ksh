#!/bin/ksh

#-- enviornment
 user=`whoami`
# datamv='login15'
 dmsdb_home=$(cat ~/.dmsrc |xargs | cut -d' ' -f 2)
 DMSPATH=/users/xa09/pkg/fx1000/dms38key/bin
 GFSDIR=$MDIR
 GFSFIX=$MDIR/fix
 GFSWRK=${GFSDIR}/work_${machine}_rst
 GFSWRKR=${GFSDIR}/work_${machine}_rst/GFSRST
 levs=128    # for model layer
 tau=48      # end of integal hours
 rstau=24    # restart hours
 rstauo=24   # interval(hours) for output restart file
 intvh=3     # for restart gfsctl, set 3 or 6 inteval hours 
 
# Caldtg="/nwpr/gfs/xb80/bin/Caldtg.ksh"

 # GODAS data source
# dataDir=/nwpr/gfs/xb99/data2
 # GODAS data
# GODASDIR=${dataDir}/obs_data/godas
 # GODAS data (Input)
# GODASpenDIR=${dataDir}/obs_data/godas/pentad
 # SIT initial data
# GODASgfsDIR="/nwpr/gfs/xb80/data2/IC_SIT"

 rm -rf $GFSWRK
 mkdir -p $GFSWRK
 mkdir -p $GFSWRKR

#-- dms data
 JCAP=${JCAP:-383}

 if [ $JCAP = 639  ] ; then
   DMSFLAG=GJ
 elif [ $JCAP = 383  ] ; then 
   DMSFLAG=GI
 elif [ $JCAP = 199  ] ; then 
   DMSFLAG=GK
 fi

 dtg='23080100'
 dtg10="20${dtg}"

# fgdtg=$(${Caldtg} ${dtg} -6)

 idmshead=TCo${JCAP}l${levs}
 idmsbody=${dtg}
 idmstail=''
 idmsdb="C384T_anal_n1"

 odmshead=O${dtg10}
 odmsbody='test'
 odmstail='rst'
 odmsdb=$idmshead

#-- executable
 EXEC='MTCo639L72_'${machine}

#---------------------------------------------------------#
 idmsfile=${idmshead}_${dtg}@${idmsdb}
 odmsfile=${odmshead}_${odmsbody}_${odmstail}${rstau}@${odmsdb}

 ${DMSPATH}/rdmsdbcrt -p ufs $idmsdb
 ${DMSPATH}/rdmscrt $idmsfile
 ${DMSPATH}/rdmsdbcrt -p ufs $odmsdb
 ${DMSPATH}/rdmscrt $odmsfile

  export LNCP='ln -fs'

  # maybe no need to change
  export source="/data/common/gfs/dms_data/${idmsdb}.ufs/${idmshead}_${idmsbody}"           # TCo IC data path
#  export source="/nwpr/gfs/xb126/data2/Tool/Nemsio2Dms_v2/OUTPUT/ncep_ana.ufs/TCo${JCAP}l72_${dtg}"           # TCo IC data path

  # link/copy DMS files
  export target="${dmsdb_home}/${idmsdb}.ufs/${idmshead}_${idmsbody}"

  # analysis
  echo ${LNCP} ${source}/* ${target}/.
       ${LNCP} ${source}/* ${target}/.
#  echo ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}
#       ${LNCP} ${source}/*${fgdtg}* ${target}/${idmshead}${idmsbody}${idmstail}

 # Bundary conditio
 BCKOPSFN="BCK_TCo${JCAP}_${DMSFLAG}30S_xnew_38"
 export BCKOPS=${BCKOPSFN}@${dmsdb_home}/bckdms
 source="/data/common/gfs/dms_data/bckdms.ufs"
 target="${dmsdb_home}/bckdms.ufs"
 if [ ! -e ${target}/${BCKOPSFN} ] ; then
   ${DMSPATH}/rdmscrt ${BCKOPS}
   ${LNCP} ${source}/${BCKOPSFN}/* ${target}/${BCKOPSFN}/.
 fi

 # get oceanic climatology
 #OCNCLMFN="OCNCLM"
 #export OCNCLM=${OCNCLMFN}@${idmsdb}
 #${DMSPATH}/rdmscrt -l34  ${OCNCLM}
 #source="/nwpr/gfs/xb99/data2/dmsdb/OISST_clim.ufs"
 #target="${dmsdb_home}/${idmsdb}.ufs"
 #if [ ! -e ${target}/${OCNCLMFN}/*  ] ; then
 #  ${LNCP} ${source}/CLMdaily_TCo383_2010_2019/* ${target}/${OCNCLMFN}/
 #fi

 # get MOM SST
 #OCNDMSFN="OCNDMS"
 #export OCNDMS=${OCNDMSFN}@${idmsdb}
 #${DMSPATH}/rdmscrt -l34 ${OCNDMS}
 #export source="/nwpr/gfs/xb80/data2/dmsdb/ncepec.ufs"
 #export target="${dmsdb_home}/${idmsdb}.ufs"
# if [ ! -e ${target}/${OCNDMSFN}/* ] ; then
#   ${LNCP} ${source}/${OCNDMSFN}/* ${target}/${OCNDMSFN}/
# fi

  # Oceanic initial and bundary condition
#  dtg_yy=`echo $dtg10 | cut -c1-4`
#  dtg_yymmdd=`echo $dtg10 | cut -c1-8`

#  dtg_yynow=${dtg_yy}
#  dtg_yybefore1=$(( $dtg_yynow -1 ))
#  dtg_yyafter1=$(( $dtg_yynow +1 ))
#  dtg_yyafter2=$(( $dtg_yynow +2 ))

#  godas_mmdd=( 0105 0110 0115 0120 0125 0130
#               0204 0209 0214 0219 0224
#               0301 0306 0311 0316 0321 0326 0331
#               0405 0410 0415 0420 0425 0430
#               0505 0510 0515 0520 0525 0530
#               0604 0609 0614 0619 0624 0629
#               0704 0709 0714 0719 0724 0729
#               0803 0808 0813 0818 0823 0823
#               0902 0907 0912 0917 0922 0927
#               1002 1007 1012 1017 1022 1027
#               1101 1106 1111 1116 1121 1126
#               1201 1206 1211 1216 1221 1226 1231 )

#  pre480=`${Caldtg} ${dtg10} -480`
#  pre480_yy=`echo $pre480 | cut -c1-4`
#  pre480_mmdd=`echo $pre480 | cut -c5-8`

#  if [ ${pre480_mmdd} -le ${godas_mmdd[0]} ]; then
#     godas_base_yy=$(( ${pre480_yy}-1 ))
#     godas_base_mmdd=${godas_mmdd[${ngodasmmdd}-1]}
#     godas_base=${godas_base_yy}${godas_base_mmdd}
#     echo godas_base=${godas_base}
#  else
#     godas_base_yy=${pre480_yy}
#     godas_base_mmdd=${godas_mmdd[0]}
#     godas_base=${godas_base_yy}${godas_base_mmdd}
#     godas_base_pre1=${godas_base_yy}${godas_base_mmdd}
#     godas_base_pre2=${godas_base_yy}${godas_base_mmdd}

#    for iday in ${godas_mmdd[@]} ; do
#      if [ $iday -le ${pre480_mmdd} ]; then
#        godas_base_pre2=${godas_base_pre1}
#        godas_base_pre1=${godas_base}
#        godas_base=${pre480_yy}${iday}
#      else
#        break
#      fi
#    done
#  fi

# if [ -f ${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base}.lev21mixed.a.nc ]; then
#   echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base}.lev21mixed.a.nc found."
# else
#   echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base}.lev21mixed.a.nc doesn't exist, check pre1 data."
#   if [ -f ${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre1}.lev21mixed.a.nc ]; then
#     echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre1}.lev21mixed.a.nc found."
#     godas_base=${godas_base_pre1}
#   else
#     echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre1}.lev21mixed.a.nc doesn't exist, check pre2 data."
#     if [ -f ${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre2}.lev21mixed.a.nc ]; then
#       echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre2}.lev21mixed.a.nc found."
#       godas_base=${godas_base_pre2}
#     else
#       echo "${GODASpenDIR}/${godas_base_yy}/split/godas.P.${godas_base_pre2}.lev21mixed.a.nc doesn't exist, exit"
#       echo "01_check_godas failed." > $stafl
#       exit 1
#     fi
#   fi
# fi

# yyyy_now=`echo $dtg10 | cut -c1-4`
# yyyymmdd=`echo $dtg10 | cut -c1-8`
# yyyy_m1=$(($yyyy_now-1))
# yyyy_p1=$(($yyyy_now+1))
# yyyy_p2=$(($yyyy_now+2))

# echo godas_base=${godas_base}

# rm -f ${GFSWRK}/dailygodas${yyyy_m1}
# ${LNCP} ${GODASgfsDIR}/godas.P.${yyyy_m1}.TCo383.a.nc           ${GFSWRK}/dailygodas${yyyy_m1}
# rm -f ${GFSWRK}/dailygodas${yyyy_now}
# ${LNCP} ${GODASgfsDIR}/godas.P.${godas_base}.pano.TCo383.a.nc   ${GFSWRK}/dailygodas${yyyy_now}
# rm -f ${GFSWRK}/dailygodas${yyyy_p1}
# ${LNCP} ${GODASgfsDIR}/godas.P.${godas_base}.pano.TCo383.a.nc   ${GFSWRK}/dailygodas${yyyy_p1}
# rm -f ${GFSWRK}/dailygodas${yyyy_p2}
# ${LNCP} ${GODASgfsDIR}/godas.P.${yyyy_p2}.TCo383.a.nc           ${GFSWRK}/dailygodas${yyyy_p2}
# rm -f ${GFSWRK}/unit.97
# ${LNCP} ${GODASgfsDIR}/TCo383_${godas_base}pano_${yyyymmdd}.nc  ${GFSWRK}/unit.97

#----------------------------------------------------------------#

#-- write out running date tag
 echo $dtg > ${GFSWRK}/crdate

#create timectl
cat> ${GFSWRK}/timectl <<eofa
$rstauo
eofa

#-- write out file list
cat > ${GFSWRK}/filist << EOF
 &filst
 ifilin='ANADMS',
 ifilout='FCSTDMS',
 cwbout='${GFSWRKR}/',
 bckfile='BCKOPS',
 phyout='${GFSWRKR}/',
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
export GLB_TYPHINI="/ncs/ncsatyp/TYP/M00/dtg/ty"
export FIXDIR=${GFSFIX}

export ANADMS=${idmsfile}
export FCSTDMS=${odmsfile}
#export BCKOPS=BCK_TCo${JCAP}_${DMSFLAG}30S_xnew@bckdms

#${DMSPATH}/rdmspurge -f FCSTDMS
#${DMSPATH}/rdmscrt -l38 FCSTDMS

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
# copy restrt file to GFSWRK
if [ $rstau -lt 100 ]; then
   rstt=00000${rstau}
elif [ $rstau -lt 1000 ]; then
   rstt=0000${rstau}
elif [ $rstau -lt 10000 ]; then
   rstt=000${rstau}
elif [ $rstau -lt 100000 ]; then
   rstt=00${rstau}
elif [ $rstau -lt 1000000 ]; then
   rstt=0${rstau}
elif [ $rstau -lt 10000000 ]; then
   rstt=${rstau}
fi
mkdir -p $GFSWRKR/phyout_${rstt}
cp ${GFSDIR}/work_${machine}/GFSRST/cwbout_${rstt} $GFSWRKR/
cp ${GFSDIR}/work_${machine}/GFSRST/phyout_${rstt}/* $GFSWRKR/phyout_${rstt}/

if [ $JCAP = 639  ] ; then
  MODLST_RES='dt=450., hfilt=1., cgw=4.2e-5,'
  MODEL_BASIC='nco=640,'
elif [ $JCAP = 383  ] ; then
  MODLST_RES='dt=600., hfilt=1., cgw=4.2e-5, cgwd=2.40, cmbk=0.60, nmgwcv=1, tofd=t, factop=60., itter=1,'
  MODEL_BASIC='nco=384,'
elif [ $JCAP = 199  ] ; then
  MODLST_RES='dt=900., hfilt=1., cgw=4.2e-5, cgwd=1.20, cmbk=0.40,'
  MODEL_BASIC='nco=200,'
fi

export MODLST_SIT="do_sit=false, updatetg=12, fsit=-99., ldailyFCTsst=false,
                   ldailyFCTicesndpt=false, dailyClm_option=0, dSITdt_intv=12., weightSIT=1.0,"

cat > ${GFSWRK}/namlsts << EOF
 &model_param
  nco=640,
  lev=${levs},
  ncld=3,
  octahedral=t,
  nout=900000000,
  io_quilting=f,
  npex=${NPEX},
  npey=${NPEY},
  ${MODEL_BASIC}
 &end

 &modlst
  taui=0.0, taue=720.0, tauo=1.0, taup=6.0, taureg=6.,
  dt=450.0,
  cstar=f, update=t, lsimpl=t,
  hfilt=1.,
  ksgeo=1, yesdia=t,
  dopbl=t, docup=t, dorad=t, dolsp=t, doshl=t, dodry=f,
  dograv=t, docgrav=t,
  donnmi=t, 
  dosppt=f, dospptout=f,
  doshum=f,
  dossst=f, doclx=t,
  cutfreq=3, nnmivm=3, doincr=f,
  hdiff=t, frad=1.0, ldiag=0,
  idg=40, jdg=108,
  itypbl=0, numreduce=8, ptmeans=800., 
  irad=2, nmland=2,
  nmcup=6, nmshl=3, nmpbl=4, nmmiph=2,
  nmgwor=2, nmgwcv=1,
  ktcup=20, cgw=4.2e-5,
  mtnvar=14, doo3l=t, pdfcloud=f,
  ioutsigr=1,
  ggdef='${DMSFLAG}0G', gmdef='${DMSFLAG}MG',
  domfc= 24., out_green=f, otgreen=3., out_hp=f,
  ndsladvh2=f,
  isot=1, ivegsrc=1, cgwd=2.40, cmbk=0.60,
  spl1= 1., spl2=100., mwhd=1.,
  ${MODLST_RES}
  ${MODLST_SIT}
 &end

 &typ
  write_mem=0,
  trk_intv=3,
  write_tau=3,
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
  ssst = 0.80,0.4,0.10,0.08,0.04
  ssst_seed = -999,-999,-999,-999,-999
  ssst_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
  ssst_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 /

 &sit_nml
   lamip= false, lwoa0= false, ldailysst= false,
 /


EOF

# cp $NWPETC/namlsts $GFSWRK/namlsts
 $NWPETC/mkgfs_rst.sh $rstau $tau $NWPETC $intvh
 cp $NWPETC/gfsctl_rst${rstau} $GFSWRK/gfsctl
 $NWPETC/mkocards_rst.sh $rstau $tau $NWPETC
 cp $NWPETC/ocards_rst${rstau} $GFSWRK/ocards

 FCT_MODEL=$MDIR/src/$EXEC
# /usr/bin/time -p mpiexec --of-proc TCo383gfs_rst.out -n $MPI ${FCT_MODEL} 
 /usr/bin/time -p mpiexec -n $MPI ${FCT_MODEL} -Wl,-T
 if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
 fi

