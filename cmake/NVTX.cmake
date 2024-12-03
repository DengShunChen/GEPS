# ##############################################################################
# Copyright (c) 2020-2023, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

# See NVTX Fortran Wrapper:
# https://docs.nvidia.com/hpc-sdk/compilers/fortran-cuda-interfaces/index.html#cfnvtx-runtime
add_compile_options(-fPIC -Wl,-export-dynamic -Minstrument=functions)
list(
  APPEND
  EX_FUNC
  fpvs
  reducepick
  cyclic_cell_ppm_intp
  vlog
  fpvsx
  vexp
  reducepickr
  syslbl
  reduceintp
  cwb_fftw_init_threads
  cwb_fftw_plan_with_nthreads
  cwb_fftw_make_planner_thread_safe
  reducepicki
  reducegrid
  zilch
  sigmap
  gathv
  geostd
  trdivv
  qsatq_2d
  def_cfl_step
  qsatq
  pythag
  vertical_cell_ppm_intp
  mpe_mpe_bcast_col_i_scalar
  cplxmtl
  fftfax
  module_radsw_main_swflux
  cyclic_cell_plm_intp
  module_radlw_main_setcoef
  mpe_mpe_bcast_col_i_scalar
  sfc_diff
  tridin
  tridi2
  sfc_sice
  module_radsw_main_setcoef
  module_radsw_main_cldprop
  sfc_diag
  module_radiation_astronomy_coszmn
  module_radiation_surface_setemis
  module_radlw_main_rtrnmr
  module_radsw_main_swflux
  module_radlw_main_setcoef
  module_radlw_main_cldprop
  sfc_ocean
  module_radsw_main_setcoef)
# Replace ";" with "," for Minstrument-exclude-func-list
set(EX_STRING "${EX_FUNC}")
string(REPLACE ";" "," EX_FUNC "${EX_STRING}")

add_compile_options(-Minstrument-exclude-func-list=${EX_FUNC})

link_libraries(-traceback -lnvhpcwrapnvtx)

add_compile_options(-DUSE_NVTX)
