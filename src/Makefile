SHELL = /bin/sh

RANLIB = ranlib
AR = ar cq
RM = rm -f
TIMCOMCPL ?= FALSE

SRCS_f90  = \
  ancrt.f90 \
  bl_ysu.f90 \
  ccoszen.f90 \
  cem1.f90 \
  check.f90 \
  check_global_values.f90 \
  chlen.f90 \
  coftrix.f90 \
  cons.f90 \
  correct.f90 \
  coszenpm.f90 \
  couvtq.f90 \
  cplxmtl.f90 \
  crayfun.f90 \
  creeve.f90 \
  ctemp2.f90 \
  cuadjtq.f90 \
  cuadjtq_n.f90 \
  cuasc.f90 \
  cuasc_n.f90 \
  cubase.f90 \
  cubasmc.f90 \
  cubasmc_n.f90 \
  cuddraf.f90 \
  cuddraf_n.f90 \
  cudlfs.f90 \
  cudlfs_n.f90 \
  cudtdq.f90 \
  cudtdq_n.f90 \
  cududv.f90 \
  cududv_n.f90 \
  cuentr.f90 \
  cuentr_n.f90 \
  cuflx.f90 \
  cuflx_n.f90 \
  cuini.f90 \
  cuini_n.f90 \
  cumastr_driv.f90 \
  cumastr_driv_n.f90 \
  cumastr.f90 \
  cumastr_n.f90 \
  cup92.f90 \
  cupcwb.f90 \
  cutype_n.f90 \
  cwb.f90 \
  dayfil.f90 \
  days.f90 \
  dcbase.f90 \
  dcyc2.f90 \
  diabat.f90 \
  diagrd.f90 \
  dmsread.f90 \
  dmswrit.f90 \
  drychk.f90 \
  drycon.f90 \
  eeps.f90 \
  eigen.f90 \
  eigenlib.f90 \
  eigrs.f90 \
  eispack.f90 \
  evpcup.f90 \
  fftx.f90 \
  fpvs.f90 \
  funcphys.f90 \
  gather_spec.f90 \
  gathv.f90 \
  gausl3.f90 \
  gauslv.f90 \
  geostd.f90 \
  get_model_param.f90 \
  get_phi.f90 \
  getphi.f90 \
  getrdy.f90 \
  gfcst.f90 \
  glbmean.f90 \
  gpvs.f90 \
  grdcon.f90 \
  grdtpw.f90 \
  gridnl_hybrid.f90 \
  gridnl_hybrid_ndsl.f90 \
  grrad.f90 \
  gscond.f90 \
  gscondp.f90 \
  gwdc.f90 \
  gwdp.f90 \
  gwdps.f90 \
  hdiffu.f90 \
  incrini.f90 \
  indexx.f90 \
  inicons.f90 \
  init_block.f90 \
  initial.f90 \
  intgrt.f90 \
  intgrt_3tl.f90 \
  invmtx.f90 \
  ioserver.f90 \
  iounitdef.f90 \
  iw3jdn.f90 \
  joinrs.f90 \
  joinsr.f90 \
  landpack.f90 \
  landsfc.f90 \
  landsfc_n.f90 \
  landsmc.f90 \
  lgndr.f90 \
  lightning_ec.f90 \
  lsp.f90 \
  machine.f90 \
  make_list.f90 \
  matrix_hybrid_cwb.f90 \
  maxp.f90 \
  mersenne_twister.f90 \
  mfpbl.f90 \
  minp.f90 \
  mixpbl.f90 \
  mixpbl_new.f90 \
  mixpbl_n.f90 \
  mo_constants.f90 \
  mo_convect_tables.f90 \
  mo_cumulus_flux.f90 \
  mod_albn.f90 \
  mod_const.f90 \
  mod_eos_ocean.f90 \
  mod_fftcom.f90 \
  mod_fj_pad.f90 \
  mod_grid.f90 \
  mod_grb2_param.f90 \
  mod_index.f90 \
  mod_leap.f90 \
  mod_mpe.f90 \
  mod_netcdf.f90 \
  mod_noah.f90 \
  mod_outflds.f90 \
  mod_ozne_def.f90 \
  mod_param.f90 \
  mod_paramt.f90 \
  mod_phygrid.f90 \
  mod_raddiag.f90 \
  mod_radn.f90 \
  mod_rank.f90 \
  mod_sit_control.f90 \
  mod_sitgrid.f90 \
  mod_sit_vdiff.f90 \
  mod_spec.f90 \
  mod_sst.f90 \
  mod_stochastic.f90 \
  mod_stochastic_physics.f90 \
  mod_typhoon.f90 \
  module_bfmicrophysics.f90 \
  module_mp_radar.f90 \
  module_mp_thompson_gfs.f90 \
  module_mp_thompson_gfs_new.f90 \
  module_mp_thompson_make_NC.f90 \
  module_mp_wsm6_fv3.f90 \
  module_mp_gfdl_fv3_v16.f90 \
  module_mp_gfdl_fv3_v2.f90 \
  module_mp_gfdl_fv3_v3.f90 \
  module_mp_gsfcgce.f90 \
  module_mp_gsfcgce_3ice_nuwrf.f90 \
  module_mp_gsfcgce_4ice_nuwrf.f90 \
  module_gocart_coupling.f90 \
  cld_eff_rad.f90 \
  moninedmf.f90 \
  mpe2d.f90 \
  mpe_barrier.f90 \
  mpe_finalize.f90 \
  mpe_gather_io.f90 \
  mpe_init.f90 \
  mpe_recv_data.f90 \
  mpe_recv_key.f90 \
  mpe_scatter_io.f90 \
  mpe_send_data.f90 \
  mpe_send_int.f90 \
  mpe_send_key.f90 \
  mpe_transpose_rs1.f90 \
  mpe_transpose_rs.f90 \
  mpe_transpose_sr1.f90 \
  mpe_transpose_sr.f90 \
  mpe_unify.f90 \
  mp_scheme.f90 \
  mtxmlp.f90 \
  mtxprt.f90 \
  mxslvc.f90 \
  namelist_soilveg.f90 \
  ncoling.f90 \
  ndslfv_monoadvh.f90 \
  ndslfv_monoadvv.f90 \
  ndslfv_pack.f90 \
  ndslfv_para.f90 \
  nnmi.f90 \
  nor_gwdp.f90 \
  normalize.f90 \
  ntke_func.f90 \
  out24.f90 \
  out2d.f90 \
  out2d_mfc.f90 \
  outflds.f90 \
  outflds_fv3.f90 \
  outflds_green.f90 \
  outflds_hp.f90 \
  outsigr.f90 \
  outsigs.f90 \
  outsit.f90 \
  pbl_noah.f90 \
  pbltke.f90 \
  pbltke_n.f90 \
  physcons.f90 \
  physpara.f90 \
  pmatrx_hybrid_cwb.f90 \
  pnmy.f90 \
  postq.f90 \
  precpd.f90 \
  precpd_n.f90 \
  precpdp.f90 \
  prerrtmg.f90 \
  prexp_hybrid_cwb.f90 \
  qmax2d.f90 \
  qmaxn3_r.f90 \
  qmaxn3_w.f90 \
  qmaxn3p.f90 \
  qprint.f90 \
  qprnt3.f90 \
  qsatq.f90 \
  radiation_aerosols.f90 \
  radiation_astronomy.f90 \
  radiation_clouds.f90 \
  radiation_gases.f90 \
  radiation_surface.f90 \
  rad_initialize.f90 \
  radir99.f90 \
  radlw_datatb.f90 \
  radlw_main.f90 \
  radlw_param.f90 \
  radsl99.f90 \
  radsw_datatb.f90 \
  radsw_main.f90 \
  radsw_param.f90 \
  radtn99.f90 \
  rayleifr.f90 \
  rcloud.f90 \
  rcmpir24c.f90 \
  rcmpir24.f90 \
  rcmpsl2c.f90 \
  rcmpsl2.f90 \
  readalb.f90 \
  readclx.f90 \
  readaeroclx.f90 \
  read_mtnvar.f90 \
  readozpl.f90 \
  reducegrid.f90 \
  reduceintp.f90 \
  reducepick.f90 \
  rkdisir.f90 \
  roptir.f90 \
  rozone.f90 \
  rozphys.f90 \
  ozphys_2015.f90 \
  rqkir.f90 \
  rrtmg.f90 \
  rsmout.f90 \
  rsmout_sigp.f90 \
  rstrandz.f90 \
  rstrantq.f90 \
  samfdeepcnv.f90 \
  samfshalcnv.f90 \
  samfdeepcnv_kh.f90 \
  samfshalcnv_kh.f90 \
  sascnv.f90 \
  sascnv_n.f90 \
  sascnv_sa.f90 \
  scatter_spec.f90 \
  send_data.f90 \
  send_idate.f90 \
  set_soilveg.f90 \
  setupv.f90 \
  sfc_diag.f90 \
  sfc_diff.f90 \
  sfc_drv.f90 \
  sfc_ocean.f90 \
  sfcpbl.f90 \
  sfc_sice.f90 \
  sfcuvt.f90 \
  sflx.f90 \
  shalcnv_new.f90 \
  shalcnv_sa.f90 \
  shlcon.f90 \
  shr_kind_mod.f90 \
  sigful.f90 \
  siimpl.f90 \
  sitocn.f90 \
  smth9.f90 \
  sortml.f90 \
  splin1.f90 \
  splinc.f90 \
  srms.f90 \
  stemp.f90 \
  structrn.f90 \
  stupij.f90 \
  stupiy.f90 \
  syslbl_r.f90 \
  syslbl_w.f90 \
  tendget.f90 \
  tracking.f90 \
  trandv.f90 \
  tranrs1.f90 \
  tranrs.f90 \
  transr1.f90 \
  transr.f90 \
  tranuv.f90 \
  tranuv1.f90 \
  trdih.f90 \
  trdiv2.f90 \
  trdivv.f90 \
  triadj.f90 \
  trigau.f90 \
  trngra3.f90 \
  trngra.f90 \
  ttstd.f90 \
  ujoinrs.f90 \
  ujoinsr.f90 \
  unify_reduceintp.f90 \
  unify_reducepick.f90 \
  uzmean.f90 \
  varpht.f90 \
  vartran.f90 \
  vartrix.f90 \
  voterp.f90 \
  vstruc_hybrid_cwb.f90 \
  vstruc_hybrid.f90 \
  vterp.f90 \
  vterpj.f90 \
  w3fs26.f90 \
  whtrec.f90 \
  whttau.f90 \
  wrte_data.f90 \
  wrte_idate.f90 \
  xyintpo.f90 \
  zilch.f90 \
  zimadv.f90 \
  zx.f90 

CPL = ${clpath}
ifeq ($(strip $(TIMCOMCPL)),TRUE)
  SRCS_f90 := $(SRCS_f90) mod_gfscpl.f90
endif

SRCS_c = cwb_fftw.c

DEPEND_FILES = $(SRCS_f) $(SRCS_f90) $(SRCS_F) $(SRCS_F90)

OBJS_f   = $(SRCS_f:.f=.o)
OBJS_f90 = $(SRCS_f90:.f90=.o)
OBJS_F   = $(SRCS_F:.F=.o)
OBJS_F90 = $(SRCS_F90:.F90=.o)
OBJS_c   = $(SRCS_c:.c=.o)

OBJS = $(OBJS_f) $(OBJS_f90) $(OBJS_F) $(OBJS_F90) $(OBJS_c)

CC   =	${C_COMP}
FC90 =  ${COMP_MPI}

INCLUDE = -I. -I${FFTW_INC} -I${NETCDF_INC} -I$(G2_INC)


CMD  =  MTCo639L72_${MACHINE}

ifeq ($(strip $(TIMCOMCPL)), TRUE)
  INCLUDE := $(INCLUDE) -I$(CPL) -I/data/common/gfs/GEPSv3_lib/$(MACHINE)/MCT/include
  CMD = libtco.a
endif

INSTALL = /usr/bin/install
# select one of the followings :
# 1. mpi single precision for transpose and unify
# 2. mpi single precision for unify only
# 3. mpi double precision 
#
 
# 1. mpi single precision for transpose and unify
#FFLAGS  = -c -X9 -Free -CcdRR8 -DMPISP -Cpp \
#          -Kfast,parallel,openmp,ocl \
#          -Ksimd=1 -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)
 
# 2. mpi single precision for unify only
#FFLAGS  = -c -X9 -Free -CcdRR8 -DMPISP -DUNIFY_ONLY -Cpp \
#          -Kfast,parallel,openmp,ocl \
#          -Ksimd=1 -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)
 
# 3. mpi double precision 
#FFLAGS  = -c -X9 -Free -CcdRR8 -Cpp \
#          -Kfast,parallel,openmp,ocl \
#          -Ksimd=1 -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)
 
# use CWB current fft
# FFLAGS  = -c -X9 -Free -CcdRR8 -Cpp -DTIMING -DMULTIPLE \
#           -Kfast,parallel,openmp,ocl,autoobjstack \
#           -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)

# use FFTW fft
# FFLAGS  = -c -X9 -Free -CcdRR8 -Cpp -DTIMING -DMULTIPLE -DUSE_FFTW -DW3TAG \
           -Kfast,parallel,openmp,ocl,autoobjstack \
	   -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)

 FFLAGS       = ${FFLAGS_COM}       ${LIBS_COM} ${OMP_COM} ${DEFINE_FLAGS} ${MSG_COM} $(INCLUDE) 
 FFLAGS_NOPAR = ${FFLAGS_COM}       ${LIBS_COM} ${DEFINE_FLAGS} ${MSG_COM} $(INCLUDE) 
 FFLAGS_GWDPS = ${FFLAGS_COM_GWDPS} ${LIBS_COM} ${OMP_COM} ${DEFINE_FLAGS} ${MSG_COM} $(INCLUDE)
 CFLAGS       = ${CFLAGS_COM} $(INCLUDE)
#           -Kfast,parallel,openmp,ocl,autoobjstack -Nquickdbg
#---------------------------------------------------------------------------------------------
# For fx10 compiler, some subroutine, such as rcloud, can not use auto-parallel option.
# This problem  may caused by compiler itself.
# If it's using fx-100 compiler(new version), you probably can get -Kparallel flag back.
#--------------------------------------------------------------------------------------------
LDFLAGS = ${OMP_COM} ${LDFLAGS_COM}
G2LIB_COM= $(G2_LIB) $(PNG_LIB) $(JASPER_LIB) $(W3LIB_LIB) $(Z_LIB)
LIBS	= ${NETCDF_LIB} ${HDF5_LIB} ${W3NCO_LIBd} ${DMS_LIB} ${OPERLIB_LIB} ${FFTW_LIB} ${LIBS_COM} ${MPMD_LIB} $(W3LIB_LIB) $(G2LIB_COM)

ifeq ($(strip $(TIMCOMCPL)),TRUE)
  FFLAGS := $(FFLAGS) -DTIMCOMCPL
  LIBS   := $(LIBS) -L$(CPL) -lcpl -L/data/common/gfs/GEPSv3_lib/$(MACHINE)/MCT/lib -lmct -lmpeu
endif

.SUFFIXES: .o .mod .f90 .c
all defalut:	depend $(CMD)

ifeq ($(strip $(TIMCOMCPL)), TRUE)
$(CMD): $(OBJS)
	@echo
	@echo '==== LINK ===================================================='
	$(RM) $(@)
	$(AR) $(@) $(OBJS)
	$(RANLIB) $(@)
else
$(CMD):	$(OBJS)
	@echo
	@echo '==== LINK ===================================================='
	$(FC90) -o $(@)  $(OBJS)  $(LDFLAGS) $(LIBS)
endif

.f90.o:
	@echo
	@echo 'compile --> '$<
	$(FC90) $(FFLAGS) $(LIBS) $<

.c.o:
	@echo
	@echo 'compile -->'$<
	$(CC) $(CFLAGS) $<

clean :
	@echo
	@echo '==== CLEAN ===================================================='
	rm -f *.o *.lst *.mod *.a depend MTCo639L72_${MACHINE}

install :
	$(INSTALL) -d ../bin
	$(INSTALL) -D $(CMD) ../bin

depend: $(DEPEND_FILES) Makefile
	@echo "==== Building dependencies for TCoGFS ====="
	@ls -1 $(DEPEND_FILES) > Srcfiles
	@echo "." > Filepath
	@$(MKDEPENDS) -m Filepath Srcfiles > depend
	@$(RM) -f Filepath Srcfiles

# automaticaly generate make dependencies
MKDEPENDS = ../mkDepends.pl
include ./depend
# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif

# old method for make dependecies
#include ./make.dependency

%.o: %.mod

include ./make.rules.$(MACHINE)

