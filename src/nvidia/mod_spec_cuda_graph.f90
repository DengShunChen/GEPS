module spec_cuda_graph
   use cudafor
   use const, only: RTYPE
   use param
   use index

   implicit none

   type spec_cg
      logical :: created = .false.
      type(cudaGraph) :: graph
      type(cudaGraphExec) :: graph_exec
      type(cudaGraphNode) :: error_node
      character(len=1) :: buffer
      integer:: buffer_len
   end type spec_cg
   type(spec_cg) trngra3_fft_cg, trngra3_lt_cg, trngra_fft_cg, &
      transr_fft_cg, transr_lt_cg, transr1_fft_cg, transr1_lt_cg, &
      tranrs_fft_cg, tranrs_lt_cg, tranrs1_fft_cg, &
      tranuv_fft_cg, tranuv_lt_cg, rstrandz_fft_cg, rstrandz_lt_cg, &
      trandv_fft_cg, trandv_lt_cg

   !! buffer
   real(kind=RTYPE), dimension(:), allocatable :: cc_cg, gwk1_cg, wcc_fk_cg
   real, dimension(:, :), allocatable :: wc_cg, ws_cg, fj_weight_cg

contains
   subroutine allocate_spec_cg_buffer
      implicit none
      integer ierr, async_id
      async_id = 1

      allocate (cc_cg((nx + 2)*levp*2*my_max), gwk1_cg((nx + 2)*levp*2*my_max), &
                wcc_fk_cg(levp*2*jtmax*my_max*nsizey*2), &
                stat=ierr)
      if (ierr /= 0) then
         write (6, *) 'spec_cuda_graph : allocate fail 1 '
         stop
      end if
      allocate (wc_cg(levp*2*my*jtmax*2, 2), &
                ws_cg(levp*2*jtrun*jtmax*2, 2), &
                fj_weight_cg((jtrun + nsizey)*(my/2)*jtmax, 2), &
                stat=ierr)
      if (ierr /= 0) then
         write (6, *) 'spec_cuda_graph : allocate fail 2 '
         stop
      end if

      !$acc enter data create(cc_cg, gwk1_cg, wcc_fk_cg, &
      !$acc                   wc_cg, ws_cg, fj_weight_cg) async(async_id)
   end subroutine allocate_spec_cg_buffer

   subroutine deallocate_spec_cg_buffer
      implicit none
      integer ierr

      !$acc exit data delete(cc_cg, gwk1_cg, wcc_fk_cg, &
      !$acc                   wc_cg, ws_cg, fj_weight_cg)
      deallocate (cc_cg, gwk1_cg, wcc_fk_cg)
      deallocate (wc_cg, ws_cg, fj_weight_cg)

   end subroutine deallocate_spec_cg_buffer

end module spec_cuda_graph
