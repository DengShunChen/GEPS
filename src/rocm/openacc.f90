! OpenACC runtime subset on HIP (amdflang has no AMDGPU OpenACC).

module openacc
  use iso_c_binding
  use cudafor, only: cuda_stream_kind, geps_acc_wait
  implicit none
  private
  integer, parameter, public :: acc_device_nvidia = 1
  integer, parameter, public :: acc_device_radeon = 2
  integer, parameter, public :: acc_device_default = 0
  integer, parameter, public :: acc_device_host = 3

  public :: acc_get_num_devices, acc_set_device_num, acc_get_device_num
  public :: acc_get_cuda_stream, acc_set_cuda_stream
  public :: geps_acc_wait
  public :: c_null_ptr
  public :: accx_async_begin_capture, accx_async_end_capture, accx_graph_launch

  type, bind(C), public :: acc_graph_t
    type(c_ptr) :: g = c_null_ptr
    type(c_ptr) :: e = c_null_ptr
  end type

  interface
    function geps_hip_device_count() bind(C, name="geps_hip_device_count")
      import c_int
      integer(c_int) :: geps_hip_device_count
    end function
    function geps_hip_set_device(id) bind(C, name="geps_hip_set_device")
      import c_int
      integer(c_int), value :: id
      integer(c_int) :: geps_hip_set_device
    end function
    function geps_hip_get_device() bind(C, name="geps_hip_get_device")
      import c_int
      integer(c_int) :: geps_hip_get_device
    end function
    function geps_acc_get_stream(async_id) bind(C, name="geps_acc_get_stream")
      import c_int, c_int64_t
      integer(c_int), value :: async_id
      integer(c_int64_t) :: geps_acc_get_stream
    end function
    subroutine geps_acc_set_stream(async_id, stream) bind(C, name="geps_acc_set_stream")
      import c_int, c_int64_t
      integer(c_int), value :: async_id
      integer(c_int64_t), value :: stream
    end subroutine
    function geps_hip_begin_capture(stream, mode) bind(C, name="geps_hip_begin_capture")
      import c_int64_t, c_int
      integer(c_int64_t), value :: stream
      integer(c_int), value :: mode
      integer(c_int) :: geps_hip_begin_capture
    end function
    function geps_hip_end_capture(stream, graph) bind(C, name="geps_hip_end_capture")
      import c_int64_t, c_ptr, c_int
      integer(c_int64_t), value :: stream
      type(c_ptr) :: graph
      integer(c_int) :: geps_hip_end_capture
    end function
    function geps_hip_graph_instantiate(exec, graph, unused) bind(C, name="geps_hip_graph_instantiate")
      import c_ptr, c_int
      type(c_ptr) :: exec
      type(c_ptr), value :: graph
      integer(c_int), value :: unused
      integer(c_int) :: geps_hip_graph_instantiate
    end function
    function geps_hip_graph_launch(exec, stream) bind(C, name="geps_hip_graph_launch")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: exec
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_hip_graph_launch
    end function
  end interface

contains

  integer function acc_get_num_devices(devtype)
    integer, intent(in) :: devtype
    acc_get_num_devices = geps_hip_device_count()
    if (devtype == acc_device_host) acc_get_num_devices = 1
  end function

  subroutine acc_set_device_num(devnum, devtype)
    integer, intent(in) :: devnum, devtype
    integer :: istat
    istat = geps_hip_set_device(devnum)
    if (devtype == acc_device_host) return
  end subroutine

  integer function acc_get_device_num(devtype)
    integer, intent(in) :: devtype
    acc_get_device_num = geps_hip_get_device()
    if (devtype == acc_device_host) acc_get_device_num = 0
  end function

  function acc_get_cuda_stream(async_id) result(stream)
    integer, intent(in) :: async_id
    integer(cuda_stream_kind) :: stream
    stream = int(geps_acc_get_stream(async_id), cuda_stream_kind)
  end function

  subroutine acc_set_cuda_stream(async_id, stream)
    integer, intent(in) :: async_id
    integer(cuda_stream_kind), intent(in) :: stream
    call geps_acc_set_stream(async_id, int(stream, c_int64_t))
  end subroutine

  subroutine accx_async_begin_capture(async_id)
    integer, intent(in) :: async_id
    integer :: istat
    istat = geps_hip_begin_capture(geps_acc_get_stream(async_id), 0)
  end subroutine

  subroutine accx_async_end_capture(async_id, graph)
    integer, intent(in) :: async_id
    type(acc_graph_t), intent(inout) :: graph
    integer :: istat
    istat = geps_hip_end_capture(geps_acc_get_stream(async_id), graph%g)
    istat = geps_hip_graph_instantiate(graph%e, graph%g, 0)
  end subroutine

  subroutine accx_graph_launch(graph, async_id)
    type(acc_graph_t), intent(in) :: graph
    integer, intent(in) :: async_id
    integer :: istat
    istat = geps_hip_graph_launch(graph%e, geps_acc_get_stream(async_id))
  end subroutine

end module openacc
