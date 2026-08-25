! cuSPARSE subset on hipSPARSE.

module cusparse
  use iso_c_binding
  use cudafor, only: cuda_stream_kind
  implicit none
  private
  integer, parameter, public :: CUSPARSE_STATUS_SUCCESS = 0

  type, bind(C), public :: cusparseHandle
    type(c_ptr) :: h = c_null_ptr
  end type

  public :: cusparseCreate, cusparseSetStream
  public :: cusparseDgtsvInterleavedBatch_bufferSizeExt
  public :: cusparseDgtsvInterleavedBatch

  interface
    function geps_sparse_create(handle) bind(C, name="geps_sparse_create")
      import c_ptr, c_int
      type(c_ptr) :: handle
      integer(c_int) :: geps_sparse_create
    end function
    function geps_sparse_set_stream(handle, stream) bind(C, name="geps_sparse_set_stream")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: handle
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_sparse_set_stream
    end function
    function geps_sparse_dgtsv_interleaved_buf(handle, algo, m, dl, d, du, x, batch, bytes) &
      bind(C, name="geps_sparse_dgtsv_interleaved_buf")
      import c_ptr, c_int, c_int64_t
      type(c_ptr), value :: handle, dl, d, du, x
      integer(c_int), value :: algo, m, batch
      integer(c_int64_t) :: bytes
      integer(c_int) :: geps_sparse_dgtsv_interleaved_buf
    end function
    function geps_sparse_dgtsv_interleaved(handle, algo, m, dl, d, du, x, batch, buf) &
      bind(C, name="geps_sparse_dgtsv_interleaved")
      import c_ptr, c_int
      type(c_ptr), value :: handle, dl, d, du, x, buf
      integer(c_int), value :: algo, m, batch
      integer(c_int) :: geps_sparse_dgtsv_interleaved
    end function
  end interface

contains

  integer function cusparseCreate(handle)
    type(cusparseHandle), intent(out) :: handle
    cusparseCreate = geps_sparse_create(handle%h)
  end function

  integer function cusparseSetStream(handle, stream)
    type(cusparseHandle), intent(in) :: handle
    integer(cuda_stream_kind), intent(in) :: stream
    cusparseSetStream = geps_sparse_set_stream(handle%h, int(stream, c_int64_t))
  end function

  integer function cusparseDgtsvInterleavedBatch_bufferSizeExt(handle, algo, m, dl, d, du, x, &
                                                               batchCount, buffer_size)
    type(cusparseHandle), intent(in) :: handle
    integer, intent(in) :: algo, m, batchCount
    real(c_double), intent(in), target :: dl(*), d(*), du(*), x(*)
    integer(c_int64_t), intent(out) :: buffer_size
    cusparseDgtsvInterleavedBatch_bufferSizeExt = geps_sparse_dgtsv_interleaved_buf( &
      handle%h, algo, m, c_loc(dl), c_loc(d), c_loc(du), c_loc(x), batchCount, buffer_size)
  end function

  integer function cusparseDgtsvInterleavedBatch(handle, algo, m, dl, d, du, x, batchCount, pbuffer)
    type(cusparseHandle), intent(in) :: handle
    integer, intent(in) :: algo, m, batchCount
    real(c_double), intent(inout), target :: dl(*), d(*), du(*), x(*)
    character, intent(inout), target :: pbuffer(*)
    cusparseDgtsvInterleavedBatch = geps_sparse_dgtsv_interleaved( &
      handle%h, algo, m, c_loc(dl), c_loc(d), c_loc(du), c_loc(x), batchCount, c_loc(pbuffer))
  end function

end module cusparse
