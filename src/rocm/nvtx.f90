! NVTX stub for HIP builds (optional NVIDIA profiling API).

module nvtx
  implicit none
  private
  public :: nvtxStartRange, nvtxEndRange, nvtxRangePush, nvtxRangePop

contains

  subroutine nvtxStartRange(name, id)
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: id
    if (len_trim(name) < 0) return
    if (present(id)) return
  end subroutine nvtxStartRange

  subroutine nvtxEndRange()
  end subroutine nvtxEndRange

  subroutine nvtxRangePush(name)
    character(len=*), intent(in) :: name
    if (len_trim(name) < 0) return
  end subroutine nvtxRangePush

  subroutine nvtxRangePop()
  end subroutine nvtxRangePop

end module nvtx
