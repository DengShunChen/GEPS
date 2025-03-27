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
add_executable(${EXEC} gfcst.f90)
target_link_libraries(${EXEC} ${LIB})
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

if(NOT ${USE_CUDA})

  # For special source files
  set_source_files_properties(fftx.f90 ndslfv_pack.f90
                              PROPERTIES COMPILE_FLAGS "-Knoparallel")
  set_source_files_properties(pbl_noah.f90 rcloud.f90 gscond.f90 gscondp.f90
                              PROPERTIES COMPILE_FLAGS "-Knosimd")
  target_link_options(${LIB} PUBLIC "-Kopenmp,parallel")
  target_link_libraries(${LIB} ${NetCDF_Fortran_LIBRARIES} ${DMS_LIBRARIES}
                        ${W3_LIBRARIES})
else()
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
endif()

if(NOT ${TIMCOMCPL})
install(TARGETS ${EXEC} DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
target_include_directories(${LIB}
                           PUBLIC "$<BUILD_INTERFACE:${TCOGFS_INCLUDE_DIRS}>")
