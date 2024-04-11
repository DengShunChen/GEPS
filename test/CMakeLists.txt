# 
# Copyright (c) 2024, NVIDIA CORPORATION.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#



add_custom_target(run
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/job  pjsub TCo383L72_IC_sample_fx100
    DEPENDS ${EXEC}
)

if(${USE_CUDA})
list(APPEND TEST_DGEMM_WRAPPER_SRC
"../src/nvidia/cublas_wrapper.f90"
"test_dgemm_wrapper.f90"
)

add_executable(test_dgemm_wrapper ${TEST_DGEMM_WRAPPER_SRC})
add_test(NAME test_dgemm_wrapper COMMAND test_dgemm_wrapper)
set_target_properties(test_dgemm_wrapper PROPERTIES LINKER_LANGUAGE Fortran)

list(APPEND TEST_FFT_WRAPPER_SRC
"../src/mod_rank.f90"
"../src/cwb_fftw.c"
"../src/fftx.f90"
"../src/nvidia/cufft_plan_cache.cc"
"../src/nvidia/cufft_wrapper.f90"
"test_fft_wrapper.f90"
)

add_executable(test_fft_wrapper ${TEST_FFT_WRAPPER_SRC})
add_test(NAME test_fft_wrapper COMMAND test_fft_wrapper)
set_target_properties(test_fft_wrapper PROPERTIES LINKER_LANGUAGE Fortran)

list(APPEND TEST_CUFFT_PLAN_CACHE_SRC
"../src/mod_rank.f90"
"../src/cwb_fftw.c"
"../src/fftx.f90"
"../src/nvidia/cufft_plan_cache.cc"
"../src/nvidia/cufft_wrapper.f90"
"test_cufft_plan_cache.f90"
)

add_executable(test_cufft_plan_cache ${TEST_CUFFT_PLAN_CACHE_SRC})
add_test(NAME test_cufft_plan_cache COMMAND test_cufft_plan_cache)
set_target_properties(test_cufft_plan_cache PROPERTIES LINKER_LANGUAGE Fortran)

list(APPEND TEST_CUFFT_LOOP_SRC
"../src/mod_rank.f90"
"../src/cwb_fftw.c"
"../src/fftx.f90"
"../src/nvidia/cufft_plan_cache.cc"
"../src/nvidia/cufft_wrapper.f90"
"test_cufft_loop.f90"
)
add_executable(test_cufft_loop ${TEST_CUFFT_LOOP_SRC})
add_test(NAME test_cufft_loop COMMAND test_cufft_loop)
set_target_properties(test_cufft_loop PROPERTIES LINKER_LANGUAGE Fortran)

if (${USE_RSM})
   set(TCO_LIB "tcogfs_rsm")
else()
   set(TCO_LIB "tcogfs")
endif()

function(subroutine_mpi_test target)
   add_executable(${target} ${target}.f90)
   add_test(NAME ${target}
      COMMAND mpiexec -n 8 ${target}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test
   )
   set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
   set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test)
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
endif()

