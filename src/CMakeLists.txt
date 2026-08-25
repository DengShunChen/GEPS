# ##############################################################################
# Copyright (c) 2020-2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

# do NOT compile unify_grid.f90 unify_gridr.f90 ozne_def.f90 physcons_v.f90

file(
  GLOB TCO_SRC
  RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "*.f90" "*.c")

list(
  APPEND
  RSM_SRC
  rsmout.f90
  send_data.f90
  send_idate.f90
  wrte_data.f90
  wrte_idate.f90
  xyintpo.f90)

# Main program is not included in the TCO sources.
if(NOT ${TIMCOMCPL})
list(REMOVE_ITEM TCO_SRC "gfcst.f90")
endif()

# RSM sources are not incluced in the TCO sources.
foreach(RSM_ITEM ${RSM_SRC})
  list(REMOVE_ITEM TCO_SRC ${RSM_ITEM})
endforeach()

# Remove files ending with "_gpu" in the TCO sources.
list(FILTER TCO_SRC EXCLUDE REGEX ".*_gpu\\.f90$")
# Remove unnecessary sources in the TCO sources.
list(REMOVE_ITEM TCO_SRC "unify_grid.f90")
list(REMOVE_ITEM TCO_SRC "unify_gridr.f90")
list(REMOVE_ITEM TCO_SRC "module_mp_gsfcgce_3ice_cwb.f90")
list(REMOVE_ITEM TCO_SRC "module_mp_gsfcgce_3ice_nuwrf_ori.f90")
list(REMOVE_ITEM TCO_SRC "module_mp_thompson_gfs_new_ori.f90")

if(GEPS_COMPILER STREQUAL "rocm")
  list(APPEND TCO_SRC "${CMAKE_CURRENT_SOURCE_DIR}/rocm/flang_compat.f90")
endif()

set(GEPS_USE_GPU OFF)
if(${USE_CUDA} OR ${USE_HIP})
  set(GEPS_USE_GPU ON)
endif()

function(geps_acc2omp in_file out_var)
  get_filename_component(_geps_abs "${in_file}" ABSOLUTE)
  file(RELATIVE_PATH _geps_rel "${CMAKE_CURRENT_SOURCE_DIR}" "${_geps_abs}")
  string(REPLACE "/" "_" _geps_flat "${_geps_rel}")
  set(_geps_out "${CMAKE_CURRENT_BINARY_DIR}/hip_src/${_geps_flat}")
  get_filename_component(_geps_out_dir "${_geps_out}" DIRECTORY)
  file(MAKE_DIRECTORY "${_geps_out_dir}")
  add_custom_command(
    OUTPUT "${_geps_out}"
    COMMAND python3 "${PROJECT_SOURCE_DIR}/cmake/acc2omp.py" "${_geps_abs}" "${_geps_out}"
    DEPENDS "${_geps_abs}" "${PROJECT_SOURCE_DIR}/cmake/acc2omp.py"
    COMMENT "acc2omp ${_geps_rel}")
  set(${out_var}
      "${_geps_out}"
      PARENT_SCOPE)
endfunction()

if(${USE_CUDA})
  file(GLOB TCO_GPU_SRC "${CMAKE_CURRENT_SOURCE_DIR}" "nvidia/*.f90"
       "nvidia/*.cc")
  list(APPEND TCO_SRC ${TCO_GPU_SRC})

  file(GLOB TCO_GPU_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*_gpu.f90")
  list(APPEND TCO_SRC ${TCO_GPU_SRC})
  list(REMOVE_ITEM TCO_SRC "mod_ndslfv_monoadv_gpu.f90")
  list(REMOVE_ITEM TCO_SRC "mpe_ndslfv_gpu.f90")
  file(GLOB TCO_GPU_SRC "${CMAKE_CURRENT_SOURCE_DIR}" "cwa/*.f90")
  list(APPEND TCO_SRC ${TCO_GPU_SRC})
elseif(${USE_HIP})
  add_library(
    geps_rocm_gpu STATIC
    rocm/hip_compat.cc
    rocm/flang_device_rt.c
    rocm/cudafor.f90
    rocm/openacc.f90
    rocm/cublas.f90
    rocm/cufft.f90
    rocm/cusolverDn.f90
    rocm/cusparse.f90
    rocm/nccl.f90
    rocm/geps_acc_wait.f90
    rocm/nvtx.f90)
  set_target_properties(
    geps_rocm_gpu PROPERTIES LINKER_LANGUAGE Fortran Fortran_MODULE_DIRECTORY
                             "${CMAKE_Fortran_MODULE_DIRECTORY}")
  target_include_directories(
    geps_rocm_gpu PUBLIC "${CMAKE_Fortran_MODULE_DIRECTORY}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/rocm")
  target_compile_options(
    geps_rocm_gpu PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:--offload-arch=${AMD_GPU_ARCHS}>"
    "$<$<COMPILE_LANGUAGE:C>:SHELL:--offload-arch=${AMD_GPU_ARCHS}>")

  file(GLOB TCO_GPU_NV "${CMAKE_CURRENT_SOURCE_DIR}/nvidia/*.f90")
  file(GLOB TCO_GPU_CC "${CMAKE_CURRENT_SOURCE_DIR}/nvidia/*.cc")
  file(GLOB TCO_GPU_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/*_gpu.f90")
  file(GLOB TCO_GPU_CWA "${CMAKE_CURRENT_SOURCE_DIR}/cwa/*.f90")
  list(REMOVE_ITEM TCO_GPU_ROOT
       "${CMAKE_CURRENT_SOURCE_DIR}/mod_ndslfv_monoadv_gpu.f90"
       "${CMAKE_CURRENT_SOURCE_DIR}/mpe_ndslfv_gpu.f90")
  set(TCO_GPU_ACC ${TCO_GPU_NV} ${TCO_GPU_ROOT} ${TCO_GPU_CWA})
  foreach(_geps_f ${TCO_GPU_ACC})
    geps_acc2omp("${_geps_f}" _geps_t)
    list(APPEND TCO_SRC "${_geps_t}")
  endforeach()
  list(APPEND TCO_SRC ${TCO_GPU_CC})

  # These live only in libmod_ndslfv_monoadv_gpu. Compiling them again into
  # libtcogfs.a duplicates OpenMP offload device symbols (lld error).
  list(REMOVE_ITEM TCO_SRC
       mod_const.f90
       mod_param.f90
       mod_index.f90
       mod_grid.f90
       mod_rank.f90
       mod_mpe.f90
       mpe_ndslfv_gpu.f90
       mod_ndslfv_monoadv_gpu.f90)

  set(_geps_host_acc
      qmaxn3_w.f90
      gpvs.f90
      mod_noah.f90
      mod_phygrid.f90
      mod_grb2_param.f90
      mod_typhoon.f90
      outflds_fv3.f90
      mod_spec.f90
      smth9.f90
      getrdy.f90
      cons.f90
      mpe_init.f90
      mpe_finalize.f90)
  foreach(_geps_f ${_geps_host_acc})
    if("${_geps_f}" IN_LIST TCO_SRC)
      list(REMOVE_ITEM TCO_SRC "${_geps_f}")
      geps_acc2omp("${CMAKE_CURRENT_SOURCE_DIR}/${_geps_f}" _geps_t)
      list(APPEND TCO_SRC "${_geps_t}")
    endif()
  endforeach()
endif()

if(NOT ${USE_RSM} AND NOT ${TIMCOMCPL})
  set(EXEC "tcogfs.x")
  set(LIB "tcogfs")
else()
  set(LIB "tcogfs")
endif()

if(NOT ${USE_RSM} AND NOT ${TIMCOMCPL})
  # add_executable(${EXEC} gfcst.f90 ${TCO_SRC})
  add_library(${LIB} STATIC ${TCO_SRC})
else()
  # add_executable(${EXEC} gfcst.f90 ${TCO_SRC} ${RSM_SRC})
  add_library(${LIB} STATIC ${TCO_SRC})
endif()

if(NOT ${TIMCOMCPL})
  if(${USE_HIP})
    geps_acc2omp("${CMAKE_CURRENT_SOURCE_DIR}/gfcst.f90" _geps_gfcst)
    add_executable(${EXEC} "${_geps_gfcst}")
    add_dependencies(${EXEC} geps_rocm_gpu)
  else()
    add_executable(${EXEC} gfcst.f90)
  endif()
  target_link_libraries(${EXEC} ${LIB})
endif()

if(${USE_HIP})
  # Only translated GPU / OpenACC sources get --offload-arch. Host .f90 stay
  # CPU OpenMP so the device LTO link does not ingest every TU.
  get_target_property(_geps_lib_srcs ${LIB} SOURCES)
  foreach(_geps_s ${_geps_lib_srcs})
    if(_geps_s MATCHES "hip_src/")
      set_property(SOURCE "${_geps_s}" APPEND PROPERTY
        COMPILE_OPTIONS "--offload-arch=${AMD_GPU_ARCHS}")
    endif()
  endforeach()
  if(NOT ${TIMCOMCPL})
    set_property(SOURCE "${_geps_gfcst}" APPEND PROPERTY
      COMPILE_OPTIONS "--offload-arch=${AMD_GPU_ARCHS}")
  endif()
endif()

if(NOT ${TIMCOMCPL})
set_target_properties(
  ${EXEC}
  PROPERTIES LINKER_LANGUAGE Fortran
             ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
             LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
             RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
else()
set_target_properties(
  ${LIB}
  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
             LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

if(GEPS_COMPILER STREQUAL "fujitsu")
  # Fujitsu-only per-file / link flags
  set_source_files_properties(fftx.f90 ndslfv_pack.f90
                              PROPERTIES COMPILE_FLAGS "-Knoparallel")
  set_source_files_properties(pbl_noah.f90 rcloud.f90 gscond.f90 gscondp.f90
                              PROPERTIES COMPILE_FLAGS "-Knosimd")
  target_link_options(${LIB} PUBLIC "-Kopenmp,parallel")
  target_link_libraries(${LIB} ${NetCDF_Fortran_LIBRARIES} ${DMS_LIBRARIES}
                        ${W3_LIBRARIES})
elseif(GEPS_COMPILER STREQUAL "nvidia")
  set_source_files_properties(pbl_noah.f90 rcloud.f90 gscond.f90 gscondp.f90
                              PROPERTIES COMPILE_FLAGS "-Mvect=nosimd")
  add_library(
    mod_ndslfv_monoadv_gpu
    mod_const.f90
    mod_param.f90
    mod_index.f90
    mod_grid.f90
    mod_rank.f90
    mpe_ndslfv_gpu.f90
    mod_mpe.f90
    mod_ndslfv_monoadv_gpu.f90)
  add_dependencies(${LIB} mod_ndslfv_monoadv_gpu)
  target_link_libraries(
    ${LIB} PUBLIC mod_ndslfv_monoadv_gpu ${NetCDF_Fortran_LIBRARIES}
                  ${DMS_LIBRARIES} ${W3_LIBRARIES})
  if(${USE_PAR})
    target_link_options(${LIB} PUBLIC "-Mconcur")
  endif()
elseif(GEPS_COMPILER STREQUAL "rocm" AND ${USE_HIP})
  set(_geps_ndsl
      mod_const.f90
      mod_param.f90
      mod_index.f90
      mod_grid.f90
      mod_rank.f90
      mpe_ndslfv_gpu.f90
      mod_mpe.f90
      mod_ndslfv_monoadv_gpu.f90)
  set(_geps_ndsl_t "")
  foreach(_geps_f ${_geps_ndsl})
    geps_acc2omp("${CMAKE_CURRENT_SOURCE_DIR}/${_geps_f}" _geps_t)
    list(APPEND _geps_ndsl_t "${_geps_t}")
  endforeach()
  add_library(mod_ndslfv_monoadv_gpu ${_geps_ndsl_t})
  target_compile_options(
    mod_ndslfv_monoadv_gpu PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:--offload-arch=${AMD_GPU_ARCHS}>")
  add_dependencies(mod_ndslfv_monoadv_gpu geps_rocm_gpu)
  add_dependencies(${LIB} mod_ndslfv_monoadv_gpu geps_rocm_gpu)
  target_link_libraries(
    ${LIB} PUBLIC mod_ndslfv_monoadv_gpu geps_rocm_gpu
                  ${NetCDF_Fortran_LIBRARIES} ${DMS_LIBRARIES} ${W3_LIBRARIES})
else()
  # gnu / intel / rocm-CPU
  target_link_libraries(${LIB} ${NetCDF_Fortran_LIBRARIES} ${DMS_LIBRARIES}
                        ${W3_LIBRARIES})
endif()

if(NOT ${TIMCOMCPL})
install(TARGETS ${EXEC} DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
target_include_directories(${LIB}
                           PUBLIC "$<BUILD_INTERFACE:${TCOGFS_INCLUDE_DIRS}>")
