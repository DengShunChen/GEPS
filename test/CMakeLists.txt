#
# Copyright (c) 2024, NVIDIA CORPORATION. Licensed under the Apache License,
# Version 2.0 (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.
#

add_custom_target(
  run
  COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/job pjsub
          TCo383L72_IC_sample_fx100
  DEPENDS ${EXEC})

if(${USE_CUDA})
  if(${USE_RSM})
    set(TCO_LIB "tcogfs_rsm")
  else()
    set(TCO_LIB "tcogfs")
  endif()

  add_executable(test_dgemm_wrapper "test_dgemm_wrapper.f90" "helper.f90")
  add_test(NAME test_dgemm_wrapper COMMAND test_dgemm_wrapper)
  set_target_properties(test_dgemm_wrapper PROPERTIES LINKER_LANGUAGE Fortran)
  target_include_directories(test_dgemm_wrapper PUBLIC ${TCO_LIB})
  target_link_libraries(test_dgemm_wrapper PUBLIC ${TCO_LIB})

  add_executable(test_fft_wrapper "test_fft_wrapper.f90" "helper.f90")
  add_test(NAME test_fft_wrapper COMMAND test_fft_wrapper)
  set_target_properties(test_fft_wrapper PROPERTIES LINKER_LANGUAGE Fortran)
  target_include_directories(test_fft_wrapper PUBLIC ${TCO_LIB})
  target_link_libraries(test_fft_wrapper PUBLIC ${TCO_LIB})

  add_executable(test_cufft_plan_cache "test_cufft_plan_cache.f90" "helper.f90")
  add_test(NAME test_cufft_plan_cache COMMAND test_cufft_plan_cache)
  set_target_properties(test_cufft_plan_cache PROPERTIES LINKER_LANGUAGE Fortran)
  target_include_directories(test_cufft_plan_cache PUBLIC ${TCO_LIB})
  target_link_libraries(test_cufft_plan_cache PUBLIC ${TCO_LIB})

  add_executable(test_cufft_loop "test_cufft_loop.f90" "helper.f90")
  add_test(NAME test_cufft_loop COMMAND test_cufft_loop)
  set_target_properties(test_cufft_loop PROPERTIES LINKER_LANGUAGE Fortran)
  target_include_directories(test_cufft_loop PUBLIC ${TCO_LIB})
  target_link_libraries(test_cufft_loop PUBLIC ${TCO_LIB})

  function(subroutine_mpi_test target)
    add_executable(${target} ${target}.f90 helper.f90)
    add_test(
      NAME ${target}
      COMMAND mpiexec -n 8 ${target}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test)
    set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                               ${PROJECT_BINARY_DIR}/test)
    target_include_directories(${target} PUBLIC ${TCO_LIB})
    target_link_libraries(${target} PUBLIC ${TCO_LIB})
  endfunction()

  subroutine_mpi_test(test_tranuv)
  subroutine_mpi_test(test_ujoinsr)
  subroutine_mpi_test(test_trngra3)
  subroutine_mpi_test(test_transr)
  subroutine_mpi_test(test_trandv)
  subroutine_mpi_test(test_joinrs)
  subroutine_mpi_test(test_tranrs)
  subroutine_mpi_test(test_rstrandz)
  subroutine_mpi_test(test_transr1)
  subroutine_mpi_test(test_prexp_hybrid_cwb)
  subroutine_mpi_test(test_trngra)
  subroutine_mpi_test(test_siimpl)
  subroutine_mpi_test(test_tranrs1)
  subroutine_mpi_test(test_gridnl_hybrid_ndsl)
  subroutine_mpi_test(test_hdiffu)
  subroutine_mpi_test(test_whdiffu)
  # subroutine_mpi_test(test_cyclic_cell_massadvx)
  # subroutine_mpi_test(test_ndslfv_monoadvh)
  # subroutine_mpi_test(test_ndslfv_monoadvh_fgnl)
  subroutine_mpi_test(test_ndslfv_monoadvv)
  subroutine_mpi_test(test_ndslfv_monoadvv_fgnl)
  subroutine_mpi_test(test_nccl)
  subroutine_mpi_test(test_mpe2d_unify_nx)
  subroutine_mpi_test(test_mpe2d_unify_my1d)
  subroutine_mpi_test(test_ndslfv_monoadvh_refactor)
  subroutine_mpi_test(test_ndslfv_monoadvh_fgnl_refactor)
  subroutine_mpi_test(test_ptotc)
  subroutine_mpi_test(test_tridi2)
  subroutine_mpi_test(test_mfpbl)
  subroutine_mpi_test(test_moninedmf)
  subroutine_mpi_test(test_gwdps)
  #subroutine_mpi_test(test_gwdps_cpus)
  subroutine_mpi_test(test_rozphys)

endif()
