################################################################################
# Copyright (c) 2020-2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
################################################################################

# do NOT compile unify_grid.f90 unify_gridr.f90 ozne_def.f90 physcons_v.f90

list(APPEND TCO_SRC
ancrt.f90
bl_ysu.f90
ccoszen.f90
cem1.f90
check.f90
check_global_values.f90
chlen.f90
coftrix.f90
cons.f90
correct.f90
coszenpm.f90
couvtq.f90
cplxmtl.f90
crayfun.f90
creeve.f90
ctemp2.f90
cuadjtq.f90
cuadjtq_n.f90
cuasc.f90
cuasc_n.f90
cubase.f90
cubasmc.f90
cubasmc_n.f90
cuddraf.f90
cuddraf_n.f90
cudlfs.f90
cudlfs_n.f90
cudtdq.f90
cudtdq_n.f90
cududv.f90
cududv_n.f90
cuentr.f90
cuentr_n.f90
cuflx.f90
cuflx_n.f90
cuini.f90
cuini_n.f90
cumastr_driv.f90
cumastr_driv_n.f90
cumastr.f90
cumastr_n.f90
cup92.f90
cupcwb.f90
cutype_n.f90
cwb.f90
dayfil.f90
days.f90
dcbase.f90
dcyc2.f90
diabat.f90
diagrd.f90
dmsread.f90
dmswrit.f90
drychk.f90
drycon.f90
eeps.f90
eigen.f90
eigenlib.f90
eigrs.f90
eispack.f90
evpcup.f90
fftx.f90
fpvs.f90
funcphys.f90
gather_spec.f90
gathv.f90
gausl3.f90
gauslv.f90
geostd.f90
get_model_param.f90
get_phi.f90
getphi.f90
getrdy.f90
glbmean.f90
gpvs.f90
grdcon.f90
grdtpw.f90
gridnl_hybrid.f90
gridnl_hybrid_ndsl.f90
grrad.f90
gscond.f90
gscondp.f90
gwdc.f90
gwdp.f90
gwdps.f90
hdiffu.f90
incrini.f90
indexx.f90
inicons.f90
init_block.f90
initial.f90
intgrt.f90
intgrt_3tl.f90
invmtx.f90
ioserver.f90
iounitdef.f90
iw3jdn.f90
joinrs.f90
joinsr.f90
landpack.f90
landsfc.f90
landsfc_n.f90
landsmc.f90
lgndr.f90
lightning_ec.f90
lsp.f90
machine.f90
make_list.f90
matrix_hybrid_cwb.f90
maxp.f90
mersenne_twister.f90
mfpbl.f90
minp.f90
mixpbl.f90
mixpbl_new.f90
mixpbl_n.f90
mo_constants.f90
mo_convect_tables.f90
mo_cumulus_flux.f90
mod_albn.f90
mod_const.f90
mod_eos_ocean.f90
mod_fftcom.f90
mod_fj_pad.f90
mod_grb2_param.f90
mod_grid.f90
mod_leap.f90
mod_index.f90
mod_mpe.f90
mod_netcdf.f90
mod_noah.f90
mod_outflds.f90
mod_ozne_def.f90
mod_param.f90
mod_paramt.f90
mod_phygrid.f90
mod_raddiag.f90
mod_radn.f90
mod_rank.f90
mod_sit_control.f90
mod_sitgrid.f90
mod_sit_vdiff.f90
mod_spec.f90
mod_sst.f90
mod_stochastic.f90
mod_stochastic_physics.f90
mod_typhoon.f90
module_bfmicrophysics.f90
module_mp_radar.f90
module_mp_thompson_gfs.f90
module_mp_thompson_gfs_new.f90
module_mp_thompson_make_NC.f90
module_mp_wsm6_fv3.f90
module_mp_gfdl_fv3_v16.f90
module_mp_gfdl_fv3_v2.f90
module_mp_gfdl_fv3_v3.f90
module_mp_gsfcgce.f90
module_mp_gsfcgce_3ice_nuwrf.f90
module_mp_gsfcgce_4ice_nuwrf.f90
cld_eff_rad.f90
moninedmf.f90
mpe2d.f90
mpe_barrier.f90
mpe_finalize.f90
mpe_gather_io.f90
mpe_init.f90
mpe_recv_data.f90
mpe_recv_key.f90
mpe_scatter_io.f90
mpe_send_data.f90
mpe_send_int.f90
mpe_send_key.f90
mpe_transpose_rs1.f90
mpe_transpose_rs.f90
mpe_transpose_sr1.f90
mpe_transpose_sr.f90
mpe_unify.f90
mp_scheme.f90
mtxmlp.f90
mtxprt.f90
mxslvc.f90
namelist_soilveg.f90
ncoling.f90
ndslfv_monoadvh.f90
ndslfv_monoadvv.f90
ndslfv_pack.f90
ndslfv_para.f90
nnmi.f90
nor_gwdp.f90
normalize.f90
ntke_func.f90
out24.f90
out2d.f90
out2d_mfc.f90
outflds.f90
outflds_fv3.f90
outflds_green.f90
outflds_hp.f90
outsigr.f90
outsigs.f90
outsit.f90
pbl_noah.f90
pbltke.f90
pbltke_n.f90
physcons.f90
physpara.f90
pmatrx_hybrid_cwb.f90
pnmy.f90
postq.f90
precpd.f90
precpd_n.f90
precpdp.f90
prerrtmg.f90
prexp_hybrid_cwb.f90
qmax2d.f90
qmaxn3_r.f90
qmaxn3_w.f90
qmaxn3p.f90
qprint.f90
qprnt3.f90
qsatq.f90
radiation_aerosols.f90
radiation_astronomy.f90
radiation_clouds.f90
radiation_gases.f90
radiation_surface.f90
rad_initialize.f90
radir99.f90
radlw_datatb.f90
radlw_main.f90
radlw_param.f90
radsl99.f90
radsw_datatb.f90
radsw_main.f90
radsw_param.f90
radtn99.f90
rayleifr.f90
rcloud.f90
rcmpir24c.f90
rcmpir24.f90
rcmpsl2c.f90
rcmpsl2.f90
readalb.f90
readclx.f90
read_mtnvar.f90
readozpl.f90
reducegrid.f90
reduceintp.f90
reducepick.f90
rkdisir.f90
roptir.f90
rozone.f90
rozphys.f90
rqkir.f90
rrtmg.f90
rsmout_sigp.f90
rstrandz.f90
rstrantq.f90
samfdeepcnv.f90
samfdeepcnv_kh.f90
samfshalcnv.f90
samfshalcnv_kh.f90
sascnv.f90
sascnv_n.f90
sascnv_sa.f90
scatter_spec.f90
set_soilveg.f90
setupv.f90
sfc_diag.f90
sfc_diff.f90
sfc_drv.f90
sfc_ocean.f90
sfcpbl.f90
sfc_sice.f90
sfcuvt.f90
sflx.f90
shalcnv_new.f90
shalcnv_sa.f90
shlcon.f90
shr_kind_mod.f90
sigful.f90
siimpl.f90
smth9.f90
sortml.f90
splin1.f90
splinc.f90
srms.f90
stemp.f90
structrn.f90
stupij.f90
stupiy.f90
syslbl_r.f90
syslbl_w.f90
tendget.f90
tracking.f90
trandv.f90
tranrs1.f90
tranrs.f90
transr1.f90
transr.f90
tranuv.f90
tranuv1.f90
trdih.f90
trdiv2.f90
trdivv.f90
triadj.f90
trigau.f90
trngra3.f90
trngra.f90
ttstd.f90
ujoinrs.f90
ujoinsr.f90
unify_reduceintp.f90
unify_reducepick.f90
uzmean.f90
varpht.f90
vartran.f90
vartrix.f90
voterp.f90
vstruc_hybrid_cwb.f90
vstruc_hybrid.f90
vterp.f90
vterpj.f90
w3fs26.f90
whtrec.f90
whttau.f90
zilch.f90
zimadv.f90
zx.f90
cwb_fftw.c
)

list(APPEND RSM_SRC
rsmout.f90
send_data.f90
send_idate.f90
wrte_data.f90
wrte_idate.f90
xyintpo.f90
)

if(${USE_CUDA})
    list(APPEND TCO_SRC
        nvidia/cuda_init.f90
        nvidia/cublas_wrapper.f90
        nvidia/cufft_plan_cache.cc
        nvidia/cufft_wrapper.f90
        nvidia/trngra3_gpu.f90
        nvidia/tranuv_gpu.f90
        nvidia/ujoinsr_gpu.f90
        nvidia/mpe2d_gpu.f90
        nvidia/mpe_transpose_rs1_gpu.f90
        nvidia/mpe_transpose_sr_gpu.f90
        nvidia/transr_gpu.f90
        nvidia/trandv_gpu.f90
        nvidia/joinrs_gpu.f90
        nvidia/mpe_transpose_rs_gpu.f90
        nvidia/tranrs_gpu.f90
        nvidia/rstrandz_gpu.f90
        nvidia/transr1_gpu.f90
        nvidia/prexp_hybrid_cwb_gpu.f90
        nvidia/trngra_gpu.f90
        nvidia/siimpl_gpu.f90
        nvidia/tranrs1_gpu.f90
        nvidia/gridnl_hybrid_ndsl_gpu.f90
        nvidia/hdiffu_gpu.f90
        nvidia/intgrt_gpu.f90
    )
endif()

if(NOT ${USE_RSM})
  set(EXEC "tcogfs.x")
  set(LIB "tcogfs")
else()
  set(EXEC "tcogfs_rsm.x")
  set(LIB "tcogfs_rsm")
endif()

if(NOT ${USE_RSM})
  #add_executable(${EXEC} gfcst.f90 ${TCO_SRC})
  add_library(${LIB} STATIC ${TCO_SRC})
else()
  #add_executable(${EXEC} gfcst.f90 ${TCO_SRC} ${RSM_SRC})
  add_library(${LIB} STATIC ${TCO_SRC} ${RSM_SRC})
endif()

add_executable(${EXEC} gfcst.f90)
target_link_libraries(${EXEC} ${LIB})

set_target_properties(${EXEC}
    PROPERTIES
    LINKER_LANGUAGE Fortran
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

if(NOT ${USE_CUDA})

	# For special source files
	set_source_files_properties(fftx.f90 ndslfv_pack.f90
		PROPERTIES COMPILE_FLAGS "-Knoparallel"
	)
	set_source_files_properties(pbl_noah.f90 rcloud.f90 gscond.f90 gscondp.f90
		PROPERTIES COMPILE_FLAGS "-Knosimd"
	)
	target_link_options(${LIB} PUBLIC "-Kopenmp,parallel")
        target_link_libraries( ${LIB} ${NetCDF_Fortran_LIBRARIES}
                ${DMS_LIBRARIES}
                ${W3_LIBRARIES}
        )
else()
	set_source_files_properties(pbl_noah.f90 rcloud.f90 gscond.f90 gscondp.f90
		PROPERTIES COMPILE_FLAGS "-Mvect=nosimd"
	)
        add_library(mod_ndslfv_monoadv_gpu
          mod_const.f90
          mod_param.f90
          mod_index.f90
          mod_grid.f90
          mod_rank.f90
          mpe_ndslfv_gpu.f90
          mod_mpe.f90
          mod_ndslfv_monoadv_gpu.f90)
        add_dependencies(${LIB} mod_ndslfv_monoadv_gpu)
        target_link_libraries( ${LIB} PUBLIC mod_ndslfv_monoadv_gpu
          ${NetCDF_Fortran_LIBRARIES}
          ${DMS_LIBRARIES}
          ${W3_LIBRARIES}
        )
        target_link_options(${LIB} PUBLIC "-Mconcur")
endif()

install(TARGETS ${EXEC} DESTINATION ${CMAKE_INSTALL_PREFIX})
target_include_directories( ${LIB} PUBLIC "$<BUILD_INTERFACE:${TCOGFS_INCLUDE_DIRS}>")
