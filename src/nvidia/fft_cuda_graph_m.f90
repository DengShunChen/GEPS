module spec_cuda_graph
   use cudafor

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
      tranuv_fft_cg, tranuv_lt_cg, rstrandz_fft_cg, rstrandz_lt_cg

end module spec_cuda_graph
