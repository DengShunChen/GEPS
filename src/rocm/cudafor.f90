! CUDA Fortran subset implemented on HIP (ROCm). Only compiled for USE_HIP.

module cudafor
  use iso_c_binding
  implicit none
  private
  integer, parameter, public :: cuda_stream_kind = c_intptr_t
  integer, parameter, public :: cudaSuccess = 0
  integer, parameter, public :: cudaStreamCaptureModeGlobal = 0
  integer, parameter, public :: cudastreamnonblocking = 1
  integer, parameter, public :: cudaEventDefault = 0

  type, bind(C), public :: cudaGraph
    type(c_ptr) :: g = c_null_ptr
  end type
  type, bind(C), public :: cudaGraphExec
    type(c_ptr) :: e = c_null_ptr
  end type
  type, bind(C), public :: cudaGraphNode
    type(c_ptr) :: n = c_null_ptr
  end type
  type, bind(C), public :: cudaEvent
    type(c_ptr) :: ev = c_null_ptr
  end type

  public :: int_ptr_kind
  public :: cudaGetErrorString
  public :: cudaStreamBeginCapture, cudaStreamEndCapture
  public :: cudaGraphInstantiate, cudaGraphLaunch
  public :: cudaEventCreate, cudaEventDestroy, cudaEventRecord
  public :: cudaStreamWaitEvent, cudaStreamSynchronize
  public :: cudaStreamCreatewithFlags
  public :: cudaMallocAsync, cudaFreeAsync
  public :: cudaMemcpy, cudaMemcpyAsync, cudaMemsetAsync
  integer, parameter, public :: cudaMemcpyHostToHost = 0
  integer, parameter, public :: cudaMemcpyHostToDevice = 1
  integer, parameter, public :: cudaMemcpyDeviceToHost = 2
  integer, parameter, public :: cudaMemcpyDeviceToDevice = 3
  public :: cudaProfilerStart, cudaProfilerStop
  public :: geps_acc_wait

  interface
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
    function geps_hip_event_create(ev) bind(C, name="geps_hip_event_create")
      import c_ptr, c_int
      type(c_ptr) :: ev
      integer(c_int) :: geps_hip_event_create
    end function
    function geps_hip_event_destroy(ev) bind(C, name="geps_hip_event_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: ev
      integer(c_int) :: geps_hip_event_destroy
    end function
    function geps_hip_event_record(ev, stream) bind(C, name="geps_hip_event_record")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: ev
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_hip_event_record
    end function
    function geps_hip_stream_wait_event(stream, ev, flags) bind(C, name="geps_hip_stream_wait_event")
      import c_ptr, c_int64_t, c_int
      integer(c_int64_t), value :: stream
      type(c_ptr), value :: ev
      integer(c_int), value :: flags
      integer(c_int) :: geps_hip_stream_wait_event
    end function
    function geps_hip_stream_create_flags(stream, flags) bind(C, name="geps_hip_stream_create_flags")
      import c_int64_t, c_int
      integer(c_int64_t) :: stream
      integer(c_int), value :: flags
      integer(c_int) :: geps_hip_stream_create_flags
    end function
    function geps_hip_stream_sync(stream) bind(C, name="geps_hip_stream_sync")
      import c_int64_t, c_int
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_hip_stream_sync
    end function
    function geps_hip_memcpy(dst, src, bytes, knd) bind(C, name="geps_hip_memcpy")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: dst, src
      integer(c_int64_t), value :: bytes
      integer(c_int), value :: knd
      integer(c_int) :: geps_hip_memcpy
    end function
    function geps_hip_memcpy_async(dst, src, bytes, knd, stream) bind(C, name="geps_hip_memcpy_async")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: dst, src
      integer(c_int64_t), value :: bytes, stream
      integer(c_int), value :: knd
      integer(c_int) :: geps_hip_memcpy_async
    end function
    function geps_hip_memset_async(dst, bytes, stream) bind(C, name="geps_hip_memset_async")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: dst
      integer(c_int64_t), value :: bytes, stream
      integer(c_int) :: geps_hip_memset_async
    end function
    function geps_hip_memset_i32_async(dst, val, count, stream) bind(C, name="geps_hip_memset_i32_async")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: dst
      integer(c_int), value :: val
      integer(c_int64_t), value :: count, stream
      integer(c_int) :: geps_hip_memset_i32_async
    end function
    function geps_hip_malloc_async(ptr, bytes, stream) bind(C, name="geps_hip_malloc_async")
      import c_ptr, c_int64_t, c_int
      type(c_ptr) :: ptr
      integer(c_int64_t), value :: bytes
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_hip_malloc_async
    end function
    function geps_hip_free_async(ptr, stream) bind(C, name="geps_hip_free_async")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: ptr
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_hip_free_async
    end function
    subroutine geps_hip_error_string(err, buf, n) bind(C, name="geps_hip_error_string")
      import c_int, c_char
      integer(c_int), value :: err, n
      character(c_char), intent(out) :: buf(*)
    end subroutine
    subroutine geps_acc_wait_c(async_id) bind(C, name="geps_hip_wait")
      import c_int
      integer(c_int), value :: async_id
    end subroutine
  end interface

  interface cudaMallocAsync
    module procedure cudaMallocAsync_r8
  end interface
  interface cudaFreeAsync
    module procedure cudaFreeAsync_r8
  end interface
  interface cudaMemsetAsync
    module procedure cudaMemsetAsync_r8, cudaMemsetAsync_r4, cudaMemsetAsync_i4
  end interface

contains

  integer function int_ptr_kind()
    int_ptr_kind = 8
  end function

  subroutine geps_acc_wait(async_id)
    integer, intent(in) :: async_id
    call geps_acc_wait_c(async_id)
  end subroutine

  function cudaGetErrorString(ierr) result(msg)
    integer, intent(in) :: ierr
    character(len=128) :: msg
    character(c_char) :: buf(128)
    integer :: i
    buf = c_null_char
    call geps_hip_error_string(ierr, buf, 128)
    msg = ' '
    do i = 1, 128
      if (buf(i) == c_null_char) exit
      msg(i:i) = buf(i)
    end do
  end function

  integer function cudaStreamBeginCapture(stream, mode)
    integer(cuda_stream_kind), intent(in) :: stream
    integer, intent(in) :: mode
    cudaStreamBeginCapture = geps_hip_begin_capture(int(stream, c_int64_t), mode)
  end function

  integer function cudaStreamEndCapture(stream, graph)
    integer(cuda_stream_kind), intent(in) :: stream
    type(cudaGraph), intent(out) :: graph
    cudaStreamEndCapture = geps_hip_end_capture(int(stream, c_int64_t), graph%g)
  end function

  integer function cudaGraphInstantiate(exec, graph, unused)
    type(cudaGraphExec), intent(out) :: exec
    type(cudaGraph), intent(in) :: graph
    integer, intent(in) :: unused
    cudaGraphInstantiate = geps_hip_graph_instantiate(exec%e, graph%g, unused)
  end function

  integer function cudaGraphLaunch(exec, stream)
    type(cudaGraphExec), intent(in) :: exec
    integer(cuda_stream_kind), intent(in) :: stream
    cudaGraphLaunch = geps_hip_graph_launch(exec%e, int(stream, c_int64_t))
  end function

  integer function cudaEventCreate(ev)
    type(cudaEvent), intent(out) :: ev
    cudaEventCreate = geps_hip_event_create(ev%ev)
  end function

  integer function cudaEventDestroy(ev)
    type(cudaEvent), intent(in) :: ev
    cudaEventDestroy = geps_hip_event_destroy(ev%ev)
  end function

  integer function cudaEventRecord(ev, stream)
    type(cudaEvent), intent(in) :: ev
    integer(cuda_stream_kind), intent(in) :: stream
    cudaEventRecord = geps_hip_event_record(ev%ev, int(stream, c_int64_t))
  end function

  integer function cudaStreamWaitEvent(stream, ev, flags)
    integer(cuda_stream_kind), intent(in) :: stream
    type(cudaEvent), intent(in) :: ev
    integer, intent(in) :: flags
    cudaStreamWaitEvent = geps_hip_stream_wait_event(int(stream, c_int64_t), ev%ev, flags)
  end function

  integer function cudaStreamSynchronize(stream)
    integer(cuda_stream_kind), intent(in) :: stream
    cudaStreamSynchronize = geps_hip_stream_sync(int(stream, c_int64_t))
  end function

  integer function cudaStreamCreatewithFlags(stream, flags)
    integer(cuda_stream_kind), intent(out) :: stream
    integer, intent(in) :: flags
    integer(c_int64_t) :: s
    s = 0
    cudaStreamCreatewithFlags = geps_hip_stream_create_flags(s, flags)
    stream = int(s, cuda_stream_kind)
  end function

  integer function cudaMallocAsync_r8(arr, n, stream)
    real(c_double), pointer, intent(inout) :: arr(:)
    integer, intent(in) :: n
    integer(cuda_stream_kind), intent(in) :: stream
    type(c_ptr) :: p
    p = c_null_ptr
    cudaMallocAsync_r8 = geps_hip_malloc_async(p, int(n, c_int64_t) * 8_c_int64_t, &
                                               int(stream, c_int64_t))
    if (cudaMallocAsync_r8 == 0) call c_f_pointer(p, arr, [n])
  end function

  integer function cudaFreeAsync_r8(arr, stream)
    real(c_double), pointer, intent(inout) :: arr(:)
    integer(cuda_stream_kind), intent(in) :: stream
    if (.not. associated(arr)) then
      cudaFreeAsync_r8 = 0
      return
    end if
    cudaFreeAsync_r8 = geps_hip_free_async(c_loc(arr(lbound(arr, 1))), int(stream, c_int64_t))
    nullify (arr)
  end function

  integer function cudaMemcpy(dst, src, n, knd)
    real(c_double), intent(inout), target :: dst(*)
    real(c_double), intent(in), target :: src(*)
    integer, intent(in) :: n, knd
    cudaMemcpy = geps_hip_memcpy(c_loc(dst), c_loc(src), int(n, c_int64_t) * 8_c_int64_t, knd)
  end function

  integer function cudaMemcpyAsync(dst, src, n, knd, stream)
    real(c_double), intent(inout), target :: dst(*)
    real(c_double), intent(in), target :: src(*)
    integer, intent(in) :: n, knd
    integer(cuda_stream_kind), intent(in) :: stream
    cudaMemcpyAsync = geps_hip_memcpy_async(c_loc(dst), c_loc(src), int(n, c_int64_t) * 8_c_int64_t, &
                                            knd, int(stream, c_int64_t))
  end function

  integer function cudaMemsetAsync_r8(dst, val, n, stream)
    real(c_double), intent(inout), target, contiguous :: dst(..)
    real(c_double), intent(in) :: val
    integer, intent(in) :: n
    integer(cuda_stream_kind), intent(in) :: stream
    if (val /= 0.0_c_double) print *, 'cudaMemsetAsync: zero-fill only, val=', val
    cudaMemsetAsync_r8 = geps_hip_memset_async(geps_cptr_r8(dst), int(n, c_int64_t) * 8_c_int64_t, &
                                               int(stream, c_int64_t))
  end function

  integer function cudaMemsetAsync_r4(dst, val, n, stream)
    real(c_double), intent(inout), target, contiguous :: dst(..)
    real(c_float), intent(in) :: val
    integer, intent(in) :: n
    integer(cuda_stream_kind), intent(in) :: stream
    if (val /= 0.0_c_float) print *, 'cudaMemsetAsync: zero-fill only, val=', val
    cudaMemsetAsync_r4 = geps_hip_memset_async(geps_cptr_r8(dst), int(n, c_int64_t) * 8_c_int64_t, &
                                               int(stream, c_int64_t))
  end function

  integer function cudaMemsetAsync_i4(dst, val, n, stream)
    integer(c_int), intent(inout), target, contiguous :: dst(..)
    integer(c_int), intent(in) :: val
    integer, intent(in) :: n
    integer(cuda_stream_kind), intent(in) :: stream
    cudaMemsetAsync_i4 = geps_hip_memset_i32_async(geps_cptr_i4(dst), val, int(n, c_int64_t), &
                                                   int(stream, c_int64_t))
  end function

  function geps_cptr_r8(dst) result(p)
    real(c_double), intent(inout), target, contiguous :: dst(..)
    type(c_ptr) :: p
    p = c_null_ptr
    select rank (dst)
    rank (1)
      p = c_loc(dst(lbound(dst, 1)))
    rank (2)
      p = c_loc(dst(lbound(dst, 1), lbound(dst, 2)))
    rank (3)
      p = c_loc(dst(lbound(dst, 1), lbound(dst, 2), lbound(dst, 3)))
    rank (4)
      p = c_loc(dst(lbound(dst, 1), lbound(dst, 2), lbound(dst, 3), lbound(dst, 4)))
    rank default
      error stop 'cudaMemsetAsync: unsupported rank'
    end select
  end function

  function geps_cptr_i4(dst) result(p)
    integer(c_int), intent(inout), target, contiguous :: dst(..)
    type(c_ptr) :: p
    p = c_null_ptr
    select rank (dst)
    rank (1)
      p = c_loc(dst(lbound(dst, 1)))
    rank (2)
      p = c_loc(dst(lbound(dst, 1), lbound(dst, 2)))
    rank (3)
      p = c_loc(dst(lbound(dst, 1), lbound(dst, 2), lbound(dst, 3)))
    rank default
      error stop 'cudaMemsetAsync: unsupported integer rank'
    end select
  end function

  subroutine cudaProfilerStart()
  end subroutine
  subroutine cudaProfilerStop()
  end subroutine

end module cudafor
