SHELL=/bin/ksh

MACHINE=fx100

SRCS  = mod_param.f90 mod_paramt.f90 mod_rank.f90 machine.f90 physcons.f90 \
	mod_const.f90 mod_index.f90 mod_fftcom.f90 mod_fj_pad.f90 mod_grid.f90 \
	mod_mpe.f90 mod_noah.f90 mod_phygrid.f90 \
	mod_spec.f90 mod_stochastic.f90 mod_typhoon.f90 get_model_param.f90 \
	mod_albn.f90 mod_ozne_def.f90 mod_radn.f90 mod_raddiag.f90 \
	mod_outflds.f90 \
	mo_constants.f90 mo_convect_tables.f90 mo_cumulus_flux.f90 shr_kind_mod.f90 \
	namelist_soilveg.f90 set_soilveg.f90 sfc_diag.f90 \
	sfc_diff.f90 sfc_drv.f90 sfc_ocean.f90 sfc_sice.f90 sflx.f90 init_block.f90 \
	ancrt.f90 cem1.f90 coftrix.f90 \
	cons.f90 correct.f90 ccoszen.f90 couvtq.f90 cplxmtl.f90 creeve.f90 \
	readozpl.f90 rozphys.f90 readalb.f90 \
	prerrtmg.f90 rrtmg.f90\
	cup92.f90 cupcwb.f90 dayfil.f90 days.f90 dcbase.f90 diabat.f90 \
        gensppt3dv5.f90  spptout.f90 \
	diagrd.f90 dmsread.f90 dmswrit.f90 drychk.f90 drycon.f90 \
	eeps.f90 eigen.f90 eigrs.f90 evpcup.f90 gwdp.f90 \
        gwdc.f90 read_mtnvar.f90 gwdps.f90 \
	gathv.f90 gausl3.f90 gauslv.f90 geostd.f90 \
	getphi.f90 getrdy.f90 gfcst.f90 glbmean.f90 \
	grdcon.f90 grdtpw.f90 hdiffu.f90 \
	incrini.f90 inicons.f90 initial.f90 intgrt.f90 invmtx.f90 lgndr.f90 lsp.f90 \
	maxp.f90 minp.f90 mixpbl.f90 mtxmlp.f90 mtxprt.f90 \
	nnmi.f90 out24.f90 out2d.f90 out2d_mfc.f90 outflds.f90 outflds_green.f90 outsigs.f90 outsigr.f90 pbltke.f90 \
	postq.f90 qmaxn3.f90 qprint.f90 stemp.f90 vterpj.f90 \
	qprnt3.f90 qsatq.f90 radir99.f90 radsl99.f90 radtn99.f90 \
	rcloud.f90 rcmpir24.f90 rcmpir24c.f90 rcmpsl2.f90 rcmpsl2c.f90 \
	rkdisir.f90 roptir.f90 rqkir.f90 mxslvc.f90 \
	readclx.f90 rozone.f90 setupv.f90 \
	sfcpbl.f90 shlcon.f90 sigful.f90 siimpl.f90 sortml.f90 \
	splin1.f90 splinc.f90 srms.f90 structrn.f90 stupij.f90 stupiy.f90 \
	syslbl.f90 tendget.f90 \
	trandv.f90 tranrs.f90 transr.f90 tranuv.f90 trdih.f90 \
	trdiv2.f90 trdivv.f90 triadj.f90 trigau.f90 trngra.f90 trngra3.f90 \
	ttstd.f90 uzmean.f90 varpht.f90 vartran.f90 vartrix.f90 \
	voterp.f90 vterp.f90 whtrec.f90 whttau.f90 \
	zilch.f90 zimadv.f90 zx.f90 check.f90 chlen.f90 \
	ctemp2.f90 ncoling.f90 \
	fftx.f90 eispack.f90 indexx.f90 eigenlib.f90 crayfun.f90 \
	qmaxn3p.f90 transr1.f90 tranrs1.f90 make_list.f90 \
	gather_spec.f90 scatter_spec.f90 \
	joinrs.f90 joinsr.f90 ujoinrs.f90 ujoinsr.f90 \
	mpe_unify.f90 \
	mpe_transpose_sr.f90 mpe_transpose_rs.f90 \
	mpe_transpose_sr1.f90 mpe_transpose_rs1.f90  \
	mpe_barrier.f90  mpe_init.f90 mpe_finalize.f90 \
	mpe_gather_io.f90 mpe_scatter_io.f90 sfcuvt.f90 \
	landpack.f90 landsfc.f90 landsmc.f90 \
	rstrandz.f90 rstrantq.f90 fpvs.f90 gpvs.f90 gscond.f90 precpd.f90 precpd_n.f90\
	sascnv.f90 sascnv_n.f90 sascnv_sa.f90 landsfc_n.f90 mixpbl_n.f90 pbltke_n.f90 \
	reducepick.f90 reducegrid.f90 pnmy.f90 reduceintp.f90 tracking.f90 smth9.f90 \
	gridnl_hybrid_ndsl.f90 vstruc_hybrid.f90 \
	matrix_hybrid_cwb.f90 pmatrx_hybrid_cwb.f90 prexp_hybrid_cwb.f90 \
	vstruc_hybrid_cwb.f90 rayleifr.f90 pbl_noah.f90 cwb.f90 nor_gwdp.f90 normalize.f90 \
	coszenpm.f90 shalcnv_new.f90 mixpbl_new.f90 shalcnv_sa.f90\
	qmax2d.f90 iw3jdn.f90 w3fs26.f90\
	physpara.f90 funcphys.f90 iounitdef.f90 \
	module_bfmicrophysics.f90 \
	random_number.f90\
	radlw_param.f90\
	radsw_param.f90\
	radlw_datatb.f90 radlw_main.f90 \
	radsw_datatb.f90 radsw_main.f90 \
	radiation_astronomy.f90 radiation_gases.f90\
	radiation_aerosols.f90\
	radiation_surface.f90 radiation_clouds.f90\
	grrad.f90 rad_initialize.f90 \
	mpe_send_key.f90 mpe_send_data.f90 \
	ioserver.f90 mpe_recv_data.f90 mpe_recv_key.f90 \
	cuadjtq.f90 cuasc.f90 cubase.f90 cubasmc.f90 cuddraf.f90 cudlfs.f90  \
	cudtdq.f90 cududv.f90 cuentr.f90 cuflx.f90 cuini.f90 cumastr.f90     \
	cumastr_driv.f90 \
        ndslfv_pack.f90 ndslfv_para.f90 ndslfv_monoadvv.f90 ndslfv_monoadvh.f90 \
	mpe2d.f90 \
	cuadjtq_n.f90 cuasc_n.f90 cubasmc_n.f90 cuddraf_n.f90\
	cudlfs_n.f90 cudtdq_n.f90 cududv_n.f90 cuentr_n.f90 cuflx_n.f90 \
	cuini_n.f90 cumastr_driv_n.f90 cumastr_n.f90 cutype_n.f90 ntke_func.f90 \
        get_phi.f90 module_mp_radar.f90 module_mp_wsm6_fv3.f90 lightning_ec.f90 \
        mod_sitgrid.f90 mod_eos_ocean.f90 mod_sit_vdiff.f90 outsit.f90 \
        mod_sst.f90 mod_netcdf.f90 mod_sit_control.f90

OBJS=	${SRCS:.f90=.o}

CC   =	fccpx100
FC90 =  mpifrtpx100
#CC   =	fccpx
#FC90 =  mpifrtpx
#FC   =  /package/$(MACHINE)/dms/bin/dmsmpiompf90
FC   =  /package/$(MACHINE)/dms/bin/dmsmpif90


INCLUDE = -I. -I/users/xa09/$(MACHINE)/fftw-3.3.7/include -I/package/$(MACHINE)/netcdf-3.6.3/include
#INCLUDE = -I. -I/nwpr/cfsres/b827/fx10/fftw-3.3.7/include

CMD  =  MTCo639L72_$(MACHINE)

# MSG  =    -Ec -Qa,d,i,p,t,x -Koptmsg=2 -Nquickdbg
 MSG  =    -Ec -Qa,d,i,p,t,x -Koptmsg=2

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
 FFLAGS  = -c -X9 -Free -CcdRR8 -Cpp -DTIMING -DMULTIPLE -DUSE_FFTW \
           -Kfast,parallel,openmp,ocl,autoobjstack \
	   -SSL2BLAMP -x- -fw $(MSG) $(INCLUDE)

#---------------------------------------------------------------------------------------------
# For fx10 compiler, some subroutine, such as rcloud, can not use auto-parallel option.
# This problem  may caused by compiler itself.
# If it's using fx-100 compiler(new version), you probably can get -Kparallel flag back.
#--------------------------------------------------------------------------------------------
LDFLAGS = -Kopenmp,parallel -L/package/$(MACHINE)/dms/dms.v4/lib -lrdms -lgdbm -L/package/$(MACHINE)/operlib/lib -lnwp \
          -L/users/xa09/$(MACHINE)/fftw-3.3.7/lib -lfftw3_threads -lfftw3 -SSL2BLAMP

LIBS	= -L/package/$(MACHINE)/netcdf-3.6.3/lib -lnetcdf

all :	$(CMD)

$(CMD):	$(OBJS)
	@echo
	@echo '==== LINK ===================================================='
	$(CC) -c -Xg cwb_fftw.c
	$(FC) -o $(@)  $(OBJS) cwb_fftw.o $(LDFLAGS) $(LIBS)

$(OBJS):
	$(FC90) $(FFLAGS) $*.f90

clean :
	@echo
	@echo '==== CLEAN ===================================================='
	rm -f *.o *.lst *.mod MTCo639L72_$(MACHINE)

install :
	$(INSTALL) -d ../bin
	$(INSTALL) -D $(CMD) ../bin

ancrt.o : ancrt.f90 
cem1.o : cem1.f90 mod_paramt.o 
check.o : check.f90 mod_fftcom.o mod_index.o mod_rank.o mod_mpe.o 
chlen.o : chlen.f90 
coftrix.o : coftrix.f90 
cons.o : cons.f90 mod_typhoon.o mod_fftcom.o mod_const.o mod_index.o mod_rank.o mod_mpe.o mod_param.o mo_convect_tables.o mod_spec.o module_mp_wsm6_fv3.o mod_sit_control.o
correct.o : correct.f90 mod_index.o 
ccoszen.o : ccoszen.f90 mod_index.o
coszenpm.o : coszenpm.f90 mod_index.o 
couvtq.o : couvtq.f90 
cplxmtl.o : cplxmtl.f90 
crayfun.o : crayfun.f90 
creeve.o : creeve.f90 
ctemp2.o : ctemp2.f90 mod_rank.o 
cup92.o : cup92.f90 mod_paramt.o funcphys.o physcons.o
cupcwb.o : cupcwb.f90 mod_paramt.o 
cwb.o : cwb.f90 
dayfil.o : dayfil.f90 mod_rank.o 
days.o : days.f90 
dcbase.o : dcbase.f90 
diabat.o : diabat.f90 mod_index.o mod_rank.o mod_mpe.o mod_radn.o physpara.o module_mp_wsm6_fv3.o mod_sitgrid.o mod_sit_vdiff.o mod_eos_ocean.o mod_netcdf.o mod_sit_control.o physcons.o mod_sst.o
diagrd.o : diagrd.f90 mod_index.o mod_rank.o mod_mpe.o 
dmsread.o : dmsread.f90 mod_index.o mod_rank.o mod_mpe.o 
dmswrit.o : dmswrit.f90 mod_index.o mod_rank.o mod_mpe.o 
drychk.o : drychk.f90 
drycon.o : drycon.f90 
eeps.o : eeps.f90 mod_paramt.o 
eigen.o : eigen.f90 
eigenlib.o : eigenlib.f90 
eigrs.o : eigrs.f90 
eispack.o : eispack.f90 
evpcup.o : evpcup.f90 mod_paramt.o 
fftx.o : fftx.f90 
	 $(FC90) $(FFLAGS) -Knoparallel $*.f90
fpvs.o : fpvs.f90 
gather_spec.o : gather_spec.f90 
gathv.o : gathv.f90 
gausl3.o : gausl3.f90 mod_rank.o 
gauslv.o : gauslv.f90 mod_paramt.o
gensppt3dv5.o : gensppt3dv5.f90 mod_rank.o mod_fftcom.o mod_index.o mod_mpe.o
spptout.o : spptout.f90 mod_mpe.o mod_rank.o mod_index.o
geostd.o : geostd.f90 
get_model_param.o : get_model_param.f90 mod_rank.o mod_paramt.o mod_const.o mod_param.o 
getphi.o : getphi.f90 
getrdy.o : getrdy.f90 mod_ozne_def.o mod_albn.o mod_radn.o mod_noah.o mod_typhoon.o mod_phygrid.o mod_fftcom.o mod_spec.o mod_grid.o mod_const.o mod_index.o mod_rank.o mod_mpe.o mod_param.o
gfcst.o : gfcst.f90 mod_const.o mod_param.o 
glbmean.o : glbmean.f90 mod_index.o mod_mpe.o 
gpvs.o : gpvs.f90 
grdcon.o : grdcon.f90 
grdtpw.o : grdtpw.f90 
gscond.o : gscond.f90 
gwdp.o : gwdp.f90 mod_paramt.o
gwdc.o : gwdc.f90 mod_paramt.o mod_rank.o mod_mpe.o mod_index.o
read_mtnvar.o : read_mtnvar.f90 mod_index.o mod_param.o
gwdps.o : gwdps.f90 mod_rank.o mod_mpe.o mod_index.o
hdiffu.o : hdiffu.f90 mod_rank.o mod_mpe.o mod_index.o mod_const.o
incrini.o : incrini.f90 mod_fftcom.o mod_spec.o mod_grid.o mod_const.o mod_index.o mod_rank.o mod_param.o 
indexx.o : indexx.f90 
inicons.o : inicons.f90 mod_index.o 
init_block.o : init_block.f90 mod_noah.o mod_const.o mod_param.o
initial.o : initial.f90 mod_fftcom.o mod_grid.o mod_spec.o mod_const.o mod_index.o mod_rank.o mod_mpe.o mod_param.o 
intgrt.o : intgrt.f90 mod_raddiag.o mod_radn.o mod_albn.o mod_noah.o mod_typhoon.o mod_phygrid.o mod_fftcom.o mod_grid.o mod_spec.o mod_const.o mod_index.o mod_rank.o mod_mpe.o mod_param.o  mod_sitgrid.o mod_sit_vdiff.o mod_netcdf.o mod_sit_control.o mod_sst.o
invmtx.o : invmtx.f90 mod_rank.o 
joinrs.o : joinrs.f90 mod_index.o 
joinsr.o : joinsr.f90 
landpack.o : landpack.f90 
landsfc.o : landsfc.f90 
landsfc_n.o : landsfc_n.f90 
landsmc.o : landsmc.f90 
lgndr.o : lgndr.f90 mod_index.o 
lsp.o : lsp.f90 mod_paramt.o 
machine.o : machine.f90 
make_list.o : make_list.f90 mod_index.o mod_rank.o mod_mpe.o mod_param.o mod_spec.o
matrix_hybrid_cwb.o : matrix_hybrid_cwb.f90 mod_param.o 
maxp.o : maxp.f90 mod_mpe.o mod_index.o 
minp.o : minp.f90 mod_mpe.o mod_index.o 
mixpbl.o : mixpbl.f90 mod_paramt.o 
mixpbl_new.o : mixpbl_new.f90 
mixpbl_n.o : mixpbl_n.f90 
mod_const.o : mod_const.f90 mod_param.o 
mod_fftcom.o : mod_fftcom.f90 
mod_fj_pad.o : mod_fj_pad.f90 
mod_grid.o : mod_grid.f90 mod_param.o mod_index.o
mod_index.o : mod_index.f90  mod_param.o mod_rank.o
mod_mpe.o : mod_mpe.f90 mod_rank.o
mod_noah.o : mod_noah.f90 mod_param.o mod_index.o
mod_param.o : mod_param.f90 
mod_paramt.o : mod_paramt.f90 
mod_phygrid.o : mod_phygrid.f90 mod_param.o mod_index.o
mod_rank.o : mod_rank.f90 
mod_spec.o : mod_spec.f90 mod_param.o mod_index.o
mod_stochastic.o : mod_stochastic.f90 
mod_typhoon.o : mod_typhoon.f90 mod_param.o 
mod_outflds.o : mod_outflds.f90 mod_index.o mod_mpe.o mod_rank.o
mpe_barrier.o : mpe_barrier.f90 
mpe_finalize.o : mpe_finalize.f90 mod_typhoon.o mod_noah.o mod_phygrid.o mod_spec.o mod_grid.o mod_const.o mod_param.o 
mpe_gather_io.o : mpe_gather_io.f90 
mpe_init.o : mpe_init.f90 mod_typhoon.o mod_noah.o mod_phygrid.o mod_spec.o mod_grid.o mod_const.o mod_rank.o mod_param.o 
mpe_scatter_io.o : mpe_scatter_io.f90 
mpe_transpose_rs1.o : mpe_transpose_rs1.f90 
mpe_transpose_rs.o : mpe_transpose_rs.f90 
mpe_transpose_sr1.o : mpe_transpose_sr1.f90 
mpe_transpose_sr.o : mpe_transpose_sr.f90 
mpe_unify.o : mpe_unify.f90 mod_index.o mod_rank.o mod_mpe.o 
#---------------------------------------------------------------------------------------
#  for rrtmg radiation - new programs
#---------------------------------------------------------------------------------------
mod_albn.o : mod_albn.f90 mod_param.o mod_index.o
mod_ozne_def.o : mod_ozne_def.f90 machine.o
mod_raddiag.o : mod_raddiag.f90 mod_param.o mod_index.o
mod_radn.o : mod_radn.f90
module_bfmicrophysics.o : module_bfmicrophysics.f90 machine.o funcphys.o physcons.o 
qmax2d.o : qmax2d.f90 mod_index.o
iw3jdn.o : iw3jdn.f90
w3fs26.o : w3fs26.f90
physcons.o : physcons.f90 machine.o
funcphys.o : funcphys.f90 physcons.o machine.o
physpara.o : physpara.f90 machine.o
iounitdef.o : iounitdef.f90
random_number.o : random_number.f90
radsw_param.o : radsw_param.f90 physpara.o
radlw_param.o : radlw_param.f90 physpara.o
radsw_datatb.o : radsw_datatb.f90 physpara.o radsw_param.o funcphys.o
radlw_datatb.o : radsw_datatb.f90 physpara.o radlw_param.o
radiation_astronomy.o : radiation_astronomy.f90 physpara.o physcons.o iounitdef.o
radiation_aerosols.o : radiation_aerosols.f90  physpara.o physcons.o iounitdef.o radsw_param.o radlw_param.o
radiation_gases.o : radiation_gases.f90 physpara.o funcphys.o physcons.o mod_ozne_def.o
radiation_surface.o : radiation_surface.f90 physpara.o physcons.o iounitdef.o
radiation_clouds.o : radiation_clouds.f90 physpara.o physcons.o module_bfmicrophysics.o iounitdef.o
radlw_main.o : radlw_main.f90 physpara.o physcons.o random_number.o radlw_param.o radlw_datatb.o
radsw_main.o : radsw_main.f90 physpara.o physcons.o random_number.o radsw_param.o radsw_datatb.o
grrad.o : grrad.f90 physpara.o physcons.o funcphys.o radiation_astronomy.o radiation_gases.o radiation_aerosols.o radiation_clouds.o radiation_surface.o radsw_param.o radsw_main.o radlw_param.o radlw_main.o
rad_initialize.o : rad_initialize.f90 physpara.o grrad.o
readozpl.o : readozpl.f90 machine.o mod_ozne_def.o
rozphys.o : rozphys.f90 physcons.o mod_ozne_def.o
prerrtmg.o : prerrtmg.f90 mod_index.o random_number.o grrad.o
rrtmg.o : rrtmg.f90 mod_mpe.o mod_rank.o mod_index.o mod_radn.o physpara.o grrad.o
outsigr.o : outsigr.f90 mod_mpe.o mod_index.o
readalb.o : readalb.f90 mod_index.o mod_rank.o mod_mpe.o
#----------------------------------------------------------------------------
ccoszen.o : ccoszen.f90 mod_index.o

mtxmlp.o : mtxmlp.f90 
mtxprt.o : mtxprt.f90 mod_rank.o 
mxslvc.o : mxslvc.f90 mod_fj_pad.o 
namelist_soilveg.o : namelist_soilveg.f90 
ncoling.o : ncoling.f90 mod_mpe.o mod_index.o 
nnmi.o : nnmi.f90 mod_rank.o 
nor_gwdp.o : nor_gwdp.f90 mod_index.o
normalize.o : normalize.f90 
out24.o : out24.f90 mod_index.o mod_mpe.o mod_rank.o mod_mpe.o
out2d.o : out2d.f90 mod_index.o mod_mpe.o mod_rank.o 
out2d_mfc.o : out2d_mfc.f90 mod_index.o mod_mpe.o mod_rank.o 
outflds.o : outflds.f90 mod_index.o mod_rank.o mod_mpe.o mod_outflds.o
outflds_green.o : outflds_green.f90 mod_index.o mod_rank.o mod_mpe.o
outsigs.o : outsigs.f90 mod_mpe.o mod_index.o 
pbl_noah.o : pbl_noah.f90 mod_paramt.o mod_index.o mod_rank.o mod_mpe.o 
	$(FC90) $(FFLAGS) -Knosimd pbl_noah.f90
pbltke.o : pbltke.f90 mod_paramt.o 
pbltke_n.o : pbltke_n.f90 mod_paramt.o 
physcons.o : physcons.f90 machine.o 
pmatrx_hybrid_cwb.o : pmatrx_hybrid_cwb.f90 mod_param.o 
pnmy.o : pnmy.f90 
postq.o : postq.f90 
precpd.o : precpd.f90
#xb110> for updated Zhao-Carr MP scheme
precpd_n.o : precpd_n.f90 
#xb110<
prexp_hybrid_cwb.o : prexp_hybrid_cwb.f90 
qmaxn3.o : qmaxn3.f90 mod_index.o mod_rank.o 
qmaxn3p.o : qmaxn3p.f90 mod_index.o mod_rank.o mod_mpe.o 
qprint.o : qprint.f90 
qprnt3.o : qprnt3.f90 mod_index.o mod_rank.o mod_mpe.o 
qsatq.o : qsatq.f90 
radir99.o : radir99.f90 
radsl99.o : radsl99.f90 
radtn99.o : radtn99.f90 
rayleifr.o : rayleifr.f90 mod_index.o mod_rank.o mod_mpe.o 
rcloud.o : rcloud.f90
	$(FC90) $(FFLAGS) -Knosimd rcloud.f90 
rcmpir24c.o : rcmpir24c.f90 mod_rank.o mod_fj_pad.o mod_paramt.o 
rcmpir24.o : rcmpir24.f90 mod_paramt.o 
rcmpsl2c.o : rcmpsl2c.f90 mod_fj_pad.o 
rcmpsl2.o : rcmpsl2.f90 
readclx.o : readclx.f90 mod_index.o 
reducegrid.o : reducegrid.f90 
reduceintp.o : reduceintp.f90 mod_index.o 
reducepick.o : reducepick.f90 
rkdisir.o : rkdisir.f90 
roptir.o : roptir.f90 
rozone.o : rozone.f90 
rqkir.o : rqkir.f90 
rstrandz.o : rstrandz.f90 mod_fftcom.o mod_paramt.o mod_index.o 
rstrantq.o : rstrantq.f90 mod_fftcom.o mod_paramt.o mod_index.o 
sascnv.o : sascnv.f90 
sascnv_n.o : sascnv_n.f90 
sascnv_sa.o : sascnv_sa.f90 
scatter_spec.o : scatter_spec.f90 
set_soilveg.o : set_soilveg.f90 namelist_soilveg.o 
setupv.o : setupv.f90 
sfc_diag.o : sfc_diag.f90 physcons.o machine.o 
sfc_diff.o : sfc_diff.f90 physcons.o machine.o 
sfc_drv.o : sfc_drv.f90 physcons.o machine.o 
sfc_ocean.o : sfc_ocean.f90 physcons.o machine.o 
sfcpbl.o : sfcpbl.f90 mod_paramt.o 
sfc_sice.o : sfc_sice.f90 physcons.o machine.o 
sfcuvt.o : sfcuvt.f90 
sflx.o : sflx.f90 namelist_soilveg.o 
shalcnv_new.o : shalcnv_new.f90 
shalcnv_sa.o : shalcnv_sa.f90 
shlcon.o : shlcon.f90 mod_paramt.o 
sigful.o : sigful.f90 mod_index.o mod_rank.o mod_mpe.o 
siimpl.o : siimpl.f90 mod_paramt.o mod_index.o mod_spec.o
smth9.o : smth9.f90 
sortml.o : sortml.f90 
splin1.o : splin1.f90 
splinc.o : splinc.f90 
srms.o : srms.f90 mod_index.o mod_mpe.o 
stemp.o : stemp.f90 
structrn.o : structrn.f90 mod_paramt.o 
stupij.o : stupij.f90 
stupiy.o : stupiy.f90 
syslbl.o : syslbl.f90 
tendget.o : tendget.f90 mod_grid.o mod_spec.o mod_const.o mod_index.o mod_rank.o mod_mpe.o mod_param.o 
tracking.o : tracking.f90 mod_index.o mod_rank.o mod_mpe.o mod_typhoon.o
trandv.o : trandv.f90 mod_fftcom.o mod_paramt.o mod_index.o 
tranrs1.o : tranrs1.f90 mod_fftcom.o mod_paramt.o mod_index.o 
tranrs.o : tranrs.f90 mod_fj_pad.o mod_fftcom.o mod_paramt.o mod_index.o 
transr1.o : transr1.f90 mod_fftcom.o mod_paramt.o mod_index.o 
transr.o : transr.f90 mod_fftcom.o mod_paramt.o mod_index.o 
tranuv.o : tranuv.f90 mod_fftcom.o mod_paramt.o mod_index.o 
trdih.o : trdih.f90 
trdiv2.o : trdiv2.f90 
trdivv.o : trdivv.f90 
triadj.o : triadj.f90 
trigau.o : trigau.f90 
trngra.o : trngra.f90 mod_fftcom.o mod_paramt.o mod_index.o 
ttstd.o : ttstd.f90 
ujoinrs.o : ujoinrs.f90 
ujoinsr.o : ujoinsr.f90 mod_index.o 
uzmean.o : uzmean.f90 mod_mpe.o mod_index.o 
varpht.o : varpht.f90 mod_index.o mod_rank.o mod_mpe.o 
vartran.o : vartran.f90 mod_index.o 
vartrix.o : vartrix.f90 mod_index.o 
voterp.o : voterp.f90 mod_index.o mod_mpe.o 
vstruc_hybrid_cwb.o : vstruc_hybrid_cwb.f90 
vstruc_hybrid.o : vstruc_hybrid.f90 
vterp.o : vterp.f90 
vterpj.o : vterpj.f90 
whtrec.o : whtrec.f90 
whttau.o : whttau.f90 
zilch.o : zilch.f90 
zimadv.o : zimadv.f90 mod_paramt.o mod_index.o 
zx.o : zx.f90 mod_index.o 
mpe_send_key.o : mpe_send_key.f90 mod_rank.o
mpe_send_data.o : mpe_send_data.f90 mod_rank.o
ioserver.o : ioserver.f90 mod_param.o mod_const.o 
mpe_recv_data.o : mpe_recv_data.f90
mpe_recv_key.o : mpe_recv_key.f90
#---------------------------------------------------------------------------------------
#  for cumastr
#---------------------------------------------------------------------------------------
cumastr_driv.o : cumastr_driv.f90 mod_index.o mod_mpe.o mod_rank.o mo_constants.o
cuadjtq.o : cuadjtq.f90 mo_convect_tables.o mo_constants.o shr_kind_mod.o 
cuasc.o : cuasc.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o 
cubase.o : cubase.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o 
cubasmc.o : cubasmc.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o 
cuddraf.o : cuddraf.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o 
cudlfs.o : cudlfs.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o 
cudtdq.o : cudtdq.f90 mo_constants.o shr_kind_mod.o 
cududv.o : cududv.f90 mo_constants.o shr_kind_mod.o 
cuentr.o : cuentr.f90 mo_constants.o shr_kind_mod.o 
cuflx.o : cuflx.f90 mo_constants.o shr_kind_mod.o 
cuini.o : cuini.f90 mo_constants.o shr_kind_mod.o 
cumastr.o : cumastr.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o 
mo_constants.o : mo_constants.f90 shr_kind_mod.o 
mo_convect_tables.o : mo_convect_tables.f90 mo_constants.o shr_kind_mod.o 
mo_cumulus_flux.o : mo_cumulus_flux.f90 shr_kind_mod.o mod_param.f90
shr_kind_mod.o : shr_kind_mod.f90 
#---------------------------------------------------------------------------------------
#  for NDSL
#---------------------------------------------------------------------------------------
trngra3.o : trngra3.f90 mod_fftcom.o mod_paramt.o mod_index.o
gridnl_hybrid_ndsl.o : gridnl_hybrid_ndsl.f90
ndslfv_pack.o : ndslfv_pack.f90 mod_rank.o mod_index.o mod_grid.o mod_param.o
ndslfv_para.o : ndslfv_para.f90 mod_rank.o mod_index.o mod_grid.o
ndslfv_monoadvh.o : ndslfv_monoadvh.f90 mod_rank.o mod_index.o mod_grid.o mod_param.o
ndslfv_monoadvv.o : ndslfv_monoadvv.f90 mod_rank.o mod_index.o mod_grid.o mod_param.o
mpe2d.o : mpe2d.f90 mod_param.o mod_index.o
#xb110, for modified new Tiedtke scheme>
cumastr_driv_n.o : cumastr_driv_n.f90 mod_index.o mod_mpe.o mod_rank.o mod_grid.o mo_constants.o
cuadjtq_n.o:cuadjtq_n.f90 mo_convect_tables.o mo_constants.o shr_kind_mod.o ntke_func.o
cuasc_n.o : cuasc_n.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o ntke_func.o
cutype_n.o : cutype_n.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o ntke_func.o
cubasmc_n.o : cubasmc_n.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o
cuddraf_n.o : cuddraf_n.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o
cudlfs_n.o : cudlfs_n.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o ntke_func.o
cudtdq_n.o : cudtdq_n.f90 mo_constants.o shr_kind_mod.o ntke_func.o
cududv_n.o : cududv_n.f90 mo_constants.o shr_kind_mod.o
cuentr_n.o : cuentr_n.f90 mo_constants.o shr_kind_mod.o
cuflx_n.o : cuflx_n.f90 mo_constants.o shr_kind_mod.o ntke_func.o
cuini_n.o : cuini_n.f90 mo_constants.o shr_kind_mod.o
cumastr_n.o : cumastr_n.f90 mo_cumulus_flux.o mo_constants.o shr_kind_mod.o
ntke_func.o : ntke_func.f90 mo_constants.o
#xb110<
get_phi.o : get_phi.f90
module_mp_radar.o : module_mp_radar.f90
module_mp_wsm6_fv3.o : module_mp_wsm6_fv3.f90 physcons.o module_mp_radar.o
#xb110, for lightning scheme>
lightning_ec.o: lightning_ec.f90 mo_constants.o
#----------SIT---------------
mod_sitgrid.o : mod_sitgrid.f90 mod_sit_control.o mod_netcdf.o
mod_sit_vdiff.o : mod_sit_vdiff.f90 mod_eos_ocean.o mod_netcdf.o mod_sst.o mod_sit_control.o
mod_sst.o: mod_sst.f90 mod_param.o mod_mpe.o mod_index.o mod_const.o mod_sit_control.o mod_netcdf.o mod_eos_ocean.o
mod_netcdf.o: mod_netcdf.f90 mod_mpe.o mod_const.o mod_sit_control.o
mod_eos_ocean.o: mod_eos_ocean.f90 
mod_sit_control.o: mod_sit_control.f90
outsit.o : outsit.f90 mod_index.o mod_rank.o mod_mpe.o mod_sitgrid.o
