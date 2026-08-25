! NCCL Fortran subset on RCCL. Types match NVIDIA CUDA Fortran usage.

module nccl
  use iso_c_binding
  use cudafor, only: cuda_stream_kind
  implicit none
  private

  integer, parameter :: nccl_ok = 0
  integer, parameter, public :: ncclSum = 0
  integer, parameter, public :: ncclMax = 2

  type, bind(C), public :: ncclResult
    integer(c_int) :: val = 0
  end type
  type(ncclResult), parameter, public :: ncclSuccess = ncclResult(nccl_ok)

  type, bind(C), public :: ncclUniqueId
    integer(c_int8_t) :: internal(128) = 0
  end type
  type, bind(C), public :: ncclComm
    type(c_ptr) :: c = c_null_ptr
  end type
  type, public :: ncclDataType
    integer :: v = 8
  end type
  type(ncclDataType), parameter, public :: ncclFloat32 = ncclDataType(7)
  type(ncclDataType), parameter, public :: ncclFloat64 = ncclDataType(8)

  public :: operator(/=), operator(==)
  public :: c_null_ptr
  public :: ncclGetUniqueId, ncclCommInitRank, ncclCommSplit, ncclCommDestroy
  public :: ncclGroupStart, ncclGroupEnd
  public :: ncclSend, ncclRecv, ncclAllReduce, ncclAllGather, ncclBroadcast
  public :: ncclGetErrorString

  interface operator(/=)
    module procedure nccl_ne
  end interface
  interface operator(==)
    module procedure nccl_eq
  end interface

  interface
    function geps_nccl_unique_id(id) bind(C, name="geps_nccl_unique_id")
      import c_ptr, c_int
      type(c_ptr), value :: id
      integer(c_int) :: geps_nccl_unique_id
    end function
    function geps_nccl_comm_init_rank(comm, nranks, id, rank) bind(C, name="geps_nccl_comm_init_rank")
      import c_ptr, c_int
      type(c_ptr) :: comm
      integer(c_int), value :: nranks, rank
      type(c_ptr), value :: id
      integer(c_int) :: geps_nccl_comm_init_rank
    end function
    function geps_nccl_comm_split(comm, color, key, newcomm) bind(C, name="geps_nccl_comm_split")
      import c_ptr, c_int
      type(c_ptr), value :: comm
      integer(c_int), value :: color, key
      type(c_ptr) :: newcomm
      integer(c_int) :: geps_nccl_comm_split
    end function
    function geps_nccl_comm_destroy(comm) bind(C, name="geps_nccl_comm_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: comm
      integer(c_int) :: geps_nccl_comm_destroy
    end function
    function geps_nccl_group_start() bind(C, name="geps_nccl_group_start")
      import c_int
      integer(c_int) :: geps_nccl_group_start
    end function
    function geps_nccl_group_end() bind(C, name="geps_nccl_group_end")
      import c_int
      integer(c_int) :: geps_nccl_group_end
    end function
    function geps_nccl_send(s, count, dtype, peer, comm, stream) bind(C, name="geps_nccl_send")
      import c_ptr, c_int, c_int64_t
      type(c_ptr), value :: s, comm
      integer(c_int64_t), value :: count, stream
      integer(c_int), value :: dtype, peer
      integer(c_int) :: geps_nccl_send
    end function
    function geps_nccl_recv(r, count, dtype, peer, comm, stream) bind(C, name="geps_nccl_recv")
      import c_ptr, c_int, c_int64_t
      type(c_ptr), value :: r, comm
      integer(c_int64_t), value :: count, stream
      integer(c_int), value :: dtype, peer
      integer(c_int) :: geps_nccl_recv
    end function
    function geps_nccl_allreduce(s, r, count, dtype, op, comm, stream) bind(C, name="geps_nccl_allreduce")
      import c_ptr, c_int, c_int64_t
      type(c_ptr), value :: s, r, comm
      integer(c_int64_t), value :: count, stream
      integer(c_int), value :: dtype, op
      integer(c_int) :: geps_nccl_allreduce
    end function
    function geps_nccl_allgather(s, r, count, dtype, comm, stream) bind(C, name="geps_nccl_allgather")
      import c_ptr, c_int, c_int64_t
      type(c_ptr), value :: s, r, comm
      integer(c_int64_t), value :: count, stream
      integer(c_int), value :: dtype
      integer(c_int) :: geps_nccl_allgather
    end function
    function geps_nccl_broadcast(s, r, count, dtype, root, comm, stream) bind(C, name="geps_nccl_broadcast")
      import c_ptr, c_int, c_int64_t
      type(c_ptr), value :: s, r, comm
      integer(c_int64_t), value :: count, stream
      integer(c_int), value :: dtype, root
      integer(c_int) :: geps_nccl_broadcast
    end function
    subroutine geps_nccl_error_string(err, buf, n) bind(C, name="geps_nccl_error_string")
      import c_int, c_char
      integer(c_int), value :: err, n
      character(c_char), intent(out) :: buf(*)
    end subroutine
  end interface

contains

  logical function nccl_ne(a, b)
    type(ncclResult), intent(in) :: a, b
    nccl_ne = a%val /= b%val
  end function
  logical function nccl_eq(a, b)
    type(ncclResult), intent(in) :: a, b
    nccl_eq = a%val == b%val
  end function

  function wrap(v) result(r)
    integer(c_int), intent(in) :: v
    type(ncclResult) :: r
    r%val = v
  end function

  integer function dtype_val(dt)
    type(ncclDataType), intent(in) :: dt
    dtype_val = dt%v
  end function

  function ncclGetUniqueId(id) result(r)
    type(ncclUniqueId), intent(out), target :: id
    type(ncclResult) :: r
    r = wrap(geps_nccl_unique_id(c_loc(id)))
  end function

  function ncclCommInitRank(comm, nranks, id, rank) result(r)
    type(ncclComm), intent(out) :: comm
    integer, intent(in) :: nranks, rank
    type(ncclUniqueId), intent(in), target :: id
    type(ncclResult) :: r
    r = wrap(geps_nccl_comm_init_rank(comm%c, nranks, c_loc(id), rank))
  end function

  function ncclCommSplit(comm, color, key, newcomm, config) result(r)
    type(ncclComm), intent(in) :: comm
    integer, intent(in) :: color, key
    type(ncclComm), intent(out) :: newcomm
    type(c_ptr), intent(in) :: config
    type(ncclResult) :: r
    r = wrap(geps_nccl_comm_split(comm%c, color, key, newcomm%c))
    if (.not. c_associated(config)) continue
  end function

  function ncclCommDestroy(comm) result(r)
    type(ncclComm), intent(in) :: comm
    type(ncclResult) :: r
    r = wrap(geps_nccl_comm_destroy(comm%c))
  end function

  function ncclGroupStart() result(r)
    type(ncclResult) :: r
    r = wrap(geps_nccl_group_start())
  end function

  function ncclGroupEnd() result(r)
    type(ncclResult) :: r
    r = wrap(geps_nccl_group_end())
  end function

  function ncclSend(buf, count, dt, peer, comm, stream) result(r)
    type(*), dimension(*), target :: buf
    integer, intent(in) :: count, peer
    type(ncclDataType), intent(in) :: dt
    type(ncclComm), intent(in) :: comm
    integer(cuda_stream_kind), intent(in) :: stream
    type(ncclResult) :: r
    r = wrap(geps_nccl_send(c_loc(buf), int(count, c_int64_t), dtype_val(dt), peer, comm%c, &
                            int(stream, c_int64_t)))
  end function

  function ncclRecv(buf, count, dt, peer, comm, stream) result(r)
    type(*), dimension(*), target :: buf
    integer, intent(in) :: count, peer
    type(ncclDataType), intent(in) :: dt
    type(ncclComm), intent(in) :: comm
    integer(cuda_stream_kind), intent(in) :: stream
    type(ncclResult) :: r
    r = wrap(geps_nccl_recv(c_loc(buf), int(count, c_int64_t), dtype_val(dt), peer, comm%c, &
                            int(stream, c_int64_t)))
  end function

  function ncclAllReduce(sbuf, rbuf, count, dt, op, comm, stream) result(r)
    type(*), dimension(*), target :: sbuf
    type(*), dimension(*), target :: rbuf
    integer, intent(in) :: count, op
    type(ncclDataType), intent(in) :: dt
    type(ncclComm), intent(in) :: comm
    integer(cuda_stream_kind), intent(in) :: stream
    type(ncclResult) :: r
    r = wrap(geps_nccl_allreduce(c_loc(sbuf), c_loc(rbuf), int(count, c_int64_t), dtype_val(dt), &
                                 op, comm%c, int(stream, c_int64_t)))
  end function

  function ncclAllGather(sbuf, rbuf, count, dt, comm, stream) result(r)
    type(*), dimension(*), target :: sbuf
    type(*), dimension(*), target :: rbuf
    integer, intent(in) :: count
    type(ncclDataType), intent(in) :: dt
    type(ncclComm), intent(in) :: comm
    integer(cuda_stream_kind), intent(in) :: stream
    type(ncclResult) :: r
    r = wrap(geps_nccl_allgather(c_loc(sbuf), c_loc(rbuf), int(count, c_int64_t), dtype_val(dt), &
                                 comm%c, int(stream, c_int64_t)))
  end function

  function ncclBroadcast(sbuf, rbuf, count, dt, root, comm, stream) result(r)
    type(*), dimension(*), target :: sbuf
    type(*), dimension(*), target :: rbuf
    integer, intent(in) :: count, root
    type(ncclDataType), intent(in) :: dt
    type(ncclComm), intent(in) :: comm
    integer(cuda_stream_kind), intent(in) :: stream
    type(ncclResult) :: r
    r = wrap(geps_nccl_broadcast(c_loc(sbuf), c_loc(rbuf), int(count, c_int64_t), dtype_val(dt), &
                                 root, comm%c, int(stream, c_int64_t)))
  end function

  function ncclGetErrorString(ierr) result(msg)
    type(ncclResult), intent(in) :: ierr
    character(len=128) :: msg
    character(c_char) :: buf(128)
    integer :: i
    buf = c_null_char
    call geps_nccl_error_string(ierr%val, buf, 128)
    msg = ' '
    do i = 1, 128
      if (buf(i) == c_null_char) exit
      msg(i:i) = buf(i)
    end do
  end function

end module nccl
