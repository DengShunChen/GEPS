! cuBLAS subset on hipBLAS.

module cublas
  use iso_c_binding
  use cudafor, only: cuda_stream_kind
  implicit none
  private
  integer, parameter, public :: CUBLAS_FILL_MODE_LOWER = 122
  integer, parameter, public :: CUBLAS_FILL_MODE_UPPER = 121
  integer, parameter, public :: HIPBLAS_OP_N = 111
  integer, parameter, public :: HIPBLAS_OP_T = 112
  integer, parameter, public :: HIPBLAS_OP_C = 113

  type, bind(C), public :: cublashandle
    type(c_ptr) :: h = c_null_ptr
  end type

  public :: cublasCreate, cublasSetStream, cublasGetHandle, dgemm

  interface
    function geps_blas_create(handle) bind(C, name="geps_blas_create")
      import c_ptr, c_int
      type(c_ptr) :: handle
      integer(c_int) :: geps_blas_create
    end function
    function geps_blas_default() bind(C, name="geps_blas_default")
      import c_ptr
      type(c_ptr) :: geps_blas_default
    end function
    function geps_blas_set_stream(handle, stream) bind(C, name="geps_blas_set_stream")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: handle
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_blas_set_stream
    end function
    function geps_blas_dgemm(handle, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) &
      bind(C, name="geps_blas_dgemm")
      import c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      integer(c_int), value :: ta, tb, m, n, k, lda, ldb, ldc
      real(c_double), value :: alpha, beta
      type(c_ptr), value :: a, b, c
      integer(c_int) :: geps_blas_dgemm
    end function
  end interface

contains

  integer function cublasCreate(handle)
    type(cublashandle), intent(out) :: handle
    cublasCreate = geps_blas_create(handle%h)
  end function

  integer function cublasSetStream(handle, stream)
    type(cublashandle), intent(in) :: handle
    integer(cuda_stream_kind), intent(in) :: stream
    cublasSetStream = geps_blas_set_stream(handle%h, int(stream, c_int64_t))
  end function

  function cublasGetHandle() result(handle)
    type(cublashandle) :: handle
    handle%h = geps_blas_default()
  end function

  integer function trans_op(ch)
    character(len=*), intent(in) :: ch
    select case (ch(1:1))
    case ('t', 'T')
      trans_op = HIPBLAS_OP_T
    case ('c', 'C')
      trans_op = HIPBLAS_OP_C
    case default
      trans_op = HIPBLAS_OP_N
    end select
  end function

  subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    character(len=*), intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(in), target :: a(lda, *), b(ldb, *)
    real(c_double), intent(inout), target :: c(ldc, *)
    type(cublashandle) :: handle
    integer :: istat
    handle = cublasGetHandle()
    istat = geps_blas_dgemm(handle%h, trans_op(transa), trans_op(transb), m, n, k, &
                            alpha, c_loc(a), lda, c_loc(b), ldb, beta, c_loc(c), ldc)
    if (istat /= 0) then
      print *, 'hipblasDgemm failed', istat
    end if
  end subroutine

end module cublas
