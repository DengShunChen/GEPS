module fft_cuda_graph
    use cudafor

    implicit none

    logical :: trngra3_graph_created = .false.
    type(cudaGraph) :: trngra3_graph
    type(cudaGraphExec) :: trngra3_graph_exec
    type(cudaGraphNode) :: trngra3_error_node
    character(len=1) :: trngra3_buffer
    integer:: trngra3_buffer_len
    logical :: trngra_graph_created = .false.
    type(cudaGraph) :: trngra_graph
    type(cudaGraphExec) :: trngra_graph_exec
    type(cudaGraphNode) :: trngra_error_node
    character(len=1) :: trngra_buffer
    integer:: trngra_buffer_len
    logical :: tranuv_graph_created = .false.
    type(cudaGraph) :: tranuv_graph
    type(cudaGraphExec) :: tranuv_graph_exec
    type(cudaGraphNode) :: tranuv_error_node
    character(len=1) :: tranuv_buffer
    integer:: tranuv_buffer_len
    logical :: transr_graph_created = .false.
    type(cudaGraph) :: transr_graph
    type(cudaGraphExec) :: transr_graph_exec
    type(cudaGraphNode) :: transr_error_node
    character(len=1) :: transr_buffer
    integer:: transr_buffer_len
    logical :: transr1_graph_created = .false.
    type(cudaGraph) :: transr1_graph
    type(cudaGraphExec) :: transr1_graph_exec
    type(cudaGraphNode) :: transr1_error_node
    character(len=1) :: transr1_buffer
    integer:: transr1_buffer_len
    logical :: tranrs_graph_created = .false.
    type(cudaGraph) :: tranrs_graph
    type(cudaGraphExec) :: tranrs_graph_exec
    type(cudaGraphNode) :: tranrs_error_node
    character(len=1) :: tranrs_buffer
    integer:: tranrs_buffer_len
    logical :: tranrs1_graph_created = .false.
    type(cudaGraph) :: tranrs1_graph
    type(cudaGraphExec) :: tranrs1_graph_exec
    type(cudaGraphNode) :: tranrs1_error_node
    character(len=1) :: tranrs1_buffer
    integer:: tranrs1_buffer_len
    logical :: rstrandz_graph_created = .false.
    type(cudaGraph) :: rstrandz_graph
    type(cudaGraphExec) :: rstrandz_graph_exec
    type(cudaGraphNode) :: rstrandz_error_node
    character(len=1) :: rstrandz_buffer
    integer:: rstrandz_buffer_len

end module fft_cuda_graph
