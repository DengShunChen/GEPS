# Copyright (c) 2024, NVIDIA CORPORATION. Licensed under the Apache License,
# Version 2.0 (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#set(CPL_PATH "/package/x86_64/nvidia/geps_coupler")
#set(TIMCOM_INC_PATH
#    ${CPL_PATH}/include
#)
#set(TIMCOM_LIB_PATH
#    ${CPL_PATH}/lib
#)
#set(TIMCOM_LINK_LIB
#    -lcpl
#)


add_custom_target(
    run
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/job pjsub
            TCo383L72_IC_sample_fx100
    DEPENDS ${EXEC}
)

if(${USE_CUDA})
    if(${USE_RSM})
        set(TCO_LIB "tcogfs_rsm")
    else()
        set(TCO_LIB "tcogfs")
    endif()

    set(TESTS
        test_adjptqintp
        test_ozphys_2015
        test_gwdc
        test_samfdeepcnv_kh
        test_samfshalcnv_kh
        test_gsfcgce_3ice_nuwrf
        test_pbl_noah
        test_sfc_sice
        test_sfc_diff
        test_sfc_drv
        test_sflx
        test_sfc_ocean
        test_fpvs
        test_sfc_diag
        test_dgemm_wrapper
        test_fft_wrapper
        test_cufft_plan_cache
        test_cufft_loop
        test_tranuv
        test_ujoinsr
        test_trngra3
        test_transr
        test_trandv
        test_joinrs
        test_tranrs
        test_rstrandz
        test_transr1
        test_prexp_hybrid_cwb
        test_trngra
        test_siimpl
        test_tranrs1
        test_gridnl_hybrid_ndsl
        test_hdiffu
        test_whdiffu
        # test_cyclic_cell_massadvx
        # test_ndslfv_monoadvh
        # test_ndslfv_monoadvh_fgnl
        test_ndslfv_monoadvv
        test_ndslfv_monoadvv_fgnl
        test_nccl
        test_mpe2d_unify_nx
        test_mpe2d_unify_my1d
        test_mpe2d_unify_lev
        test_ndslfv_monoadvh_refactor
        test_ndslfv_monoadvh_fgnl_refactor
        test_ptotc
        test_tridi2
        test_tridin
        test_mfpbl
        test_moninedmf
        test_gwdps
        # test_gwdps_cpus
        test_rozphys
        test_get_phi
        test_rayleifr
        test_tendget
        test_zx
        test_correct
        test_vartrix
        test_mpe2d_row_broadcast
        test_nnmi
        test_eigen
        #test_initial
        #test_sflx_file
        #test_pbl_noah_file
        #test_gsfcgce_3ice_nuwrf_file
        #test_fall_flux_file
        #test_saticel_s_file
        #test_samfdeepcnv_kh_file
        #test_samfshalcnv_kh_file
        #test_nor_gwdp_file
        #test_ozphys_2015_file
        #test_gwdc_file

    )

    foreach(TEST ${TESTS})
        add_executable(
            ${TEST}
            ${TEST}.f90
            helper.f90
            mod_vcopy.f90
	    #cpl_driver.f90
    	)
        add_test(
            NAME ${TEST}
            COMMAND mpiexec -n 8 ${TEST}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test
        )
        set_target_properties(
            ${TEST}
            PROPERTIES
            LINKER_LANGUAGE Fortran
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test
        )
        target_include_directories(
            ${TEST}
            PUBLIC ${TCO_LIB} #${TIMCOM_INC_PATH}
        )
#        target_link_directories(
#            ${TEST}
#            PUBLIC ${TIMCOM_LIB_PATH}
#        )
        target_link_libraries(
            ${TEST}
	    PUBLIC ${TCO_LIB} #${TIMCOM_LINK_LIB}
        )

    endforeach()
endif()
