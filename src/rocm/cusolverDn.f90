! cuSOLVER dense subset on hipSOLVER.

module cusolverDn
  use iso_c_binding
  use cudafor, only: cuda_stream_kind
  use cublas, only: CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER
  implicit none
  private
  integer, parameter, public :: CUSOLVER_STATUS_SUCCESS = 0
  integer, parameter, public :: CUSOLVER_EIG_MODE_NOVECTOR = 201
  integer, parameter, public :: CUSOLVER_EIG_MODE_VECTOR = 202
  public :: CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER

  type, bind(C), public :: cusolverDnHandle
    type(c_ptr) :: h = c_null_ptr
  end type
  type, bind(C), public :: cusolverDnSyevjInfo
    type(c_ptr) :: p = c_null_ptr
  end type

  public :: cusolverDnCreate, cusolverDnDestroy, cusolverDnSetStream
  public :: cusolverDnCreateSyevjInfo
  public :: cusolverDnXsyevjSetTolerance, cusolverDnXsyevjSetSortEig
  public :: cusolverDnDsyevd_bufferSize, cusolverDnDsyevd
  public :: cusolverDnDsyevj_bufferSize, cusolverDnDsyevj

  interface
    function geps_solver_create(handle) bind(C, name="geps_solver_create")
      import c_ptr, c_int
      type(c_ptr) :: handle
      integer(c_int) :: geps_solver_create
    end function
    function geps_solver_destroy(handle) bind(C, name="geps_solver_destroy")
      import c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: geps_solver_destroy
    end function
    function geps_solver_set_stream(handle, stream) bind(C, name="geps_solver_set_stream")
      import c_ptr, c_int64_t, c_int
      type(c_ptr), value :: handle
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_solver_set_stream
    end function
    function geps_solver_syevj_create(info) bind(C, name="geps_solver_syevj_create")
      import c_ptr, c_int
      type(c_ptr) :: info
      integer(c_int) :: geps_solver_syevj_create
    end function
    function geps_solver_syevj_set_tol(info, tol) bind(C, name="geps_solver_syevj_set_tol")
      import c_ptr, c_double, c_int
      type(c_ptr), value :: info
      real(c_double), value :: tol
      integer(c_int) :: geps_solver_syevj_set_tol
    end function
    function geps_solver_syevj_set_sort(info, sort) bind(C, name="geps_solver_syevj_set_sort")
      import c_ptr, c_int
      type(c_ptr), value :: info
      integer(c_int), value :: sort
      integer(c_int) :: geps_solver_syevj_set_sort
    end function
    function geps_solver_dsyevd_buf(handle, jobz, uplo, n, a, lda, w, lwork) &
      bind(C, name="geps_solver_dsyevd_buf")
      import c_ptr, c_int
      type(c_ptr), value :: handle, a, w
      integer(c_int), value :: jobz, uplo, n, lda
      integer(c_int) :: lwork
      integer(c_int) :: geps_solver_dsyevd_buf
    end function
    function geps_solver_dsyevd(handle, jobz, uplo, n, a, lda, w, work, lwork, devinfo) &
      bind(C, name="geps_solver_dsyevd")
      import c_ptr, c_int
      type(c_ptr), value :: handle, a, w, work, devinfo
      integer(c_int), value :: jobz, uplo, n, lda, lwork
      integer(c_int) :: geps_solver_dsyevd
    end function
    function geps_solver_dsyevj_buf(handle, jobz, uplo, n, a, lda, w, lwork, params) &
      bind(C, name="geps_solver_dsyevj_buf")
      import c_ptr, c_int
      type(c_ptr), value :: handle, a, w, params
      integer(c_int), value :: jobz, uplo, n, lda
      integer(c_int) :: lwork
      integer(c_int) :: geps_solver_dsyevj_buf
    end function
    function geps_solver_dsyevj(handle, jobz, uplo, n, a, lda, w, work, lwork, devinfo, params) &
      bind(C, name="geps_solver_dsyevj")
      import c_ptr, c_int
      type(c_ptr), value :: handle, a, w, work, devinfo, params
      integer(c_int), value :: jobz, uplo, n, lda, lwork
      integer(c_int) :: geps_solver_dsyevj
    end function
  end interface

contains

  integer function cusolverDnCreate(handle)
    type(cusolverDnHandle), intent(out) :: handle
    cusolverDnCreate = geps_solver_create(handle%h)
  end function

  integer function cusolverDnDestroy(handle)
    type(cusolverDnHandle), intent(in) :: handle
    cusolverDnDestroy = geps_solver_destroy(handle%h)
  end function

  integer function cusolverDnSetStream(handle, stream)
    type(cusolverDnHandle), intent(in) :: handle
    integer(cuda_stream_kind), intent(in) :: stream
    cusolverDnSetStream = geps_solver_set_stream(handle%h, int(stream, c_int64_t))
  end function

  integer function cusolverDnCreateSyevjInfo(params)
    type(cusolverDnSyevjInfo), intent(out) :: params
    cusolverDnCreateSyevjInfo = geps_solver_syevj_create(params%p)
  end function

  integer function cusolverDnXsyevjSetTolerance(params, tol)
    type(cusolverDnSyevjInfo), intent(in) :: params
    real(c_double), intent(in) :: tol
    cusolverDnXsyevjSetTolerance = geps_solver_syevj_set_tol(params%p, tol)
  end function

  integer function cusolverDnXsyevjSetSortEig(params, sort)
    type(cusolverDnSyevjInfo), intent(in) :: params
    integer, intent(in) :: sort
    cusolverDnXsyevjSetSortEig = geps_solver_syevj_set_sort(params%p, sort)
  end function

  integer function cusolverDnDsyevd_bufferSize(handle, jobz, uplo, n, a, lda, w, lwork)
    type(cusolverDnHandle), intent(in) :: handle
    integer, intent(in) :: jobz, uplo, n, lda
    real(c_double), intent(in), target :: a(lda, *), w(*)
    integer, intent(out) :: lwork
    integer(c_int) :: lw
    lw = 0
    cusolverDnDsyevd_bufferSize = geps_solver_dsyevd_buf(handle%h, jobz, uplo, n, c_loc(a), lda, c_loc(w), lw)
    lwork = lw
  end function

  integer function cusolverDnDsyevd(handle, jobz, uplo, n, a, lda, w, work, lwork, devinfo)
    type(cusolverDnHandle), intent(in) :: handle
    integer, intent(in) :: jobz, uplo, n, lda, lwork
    real(c_double), intent(inout), target :: a(lda, *), w(*), work(*)
    integer, intent(out), target :: devinfo
    cusolverDnDsyevd = geps_solver_dsyevd(handle%h, jobz, uplo, n, c_loc(a), lda, c_loc(w), &
                                          c_loc(work), lwork, c_loc(devinfo))
  end function

  integer function cusolverDnDsyevj_bufferSize(handle, jobz, uplo, n, a, lda, w, lwork, params)
    type(cusolverDnHandle), intent(in) :: handle
    integer, intent(in) :: jobz, uplo, n, lda
    real(c_double), intent(in), target :: a(lda, *), w(*)
    integer, intent(out) :: lwork
    type(cusolverDnSyevjInfo), intent(in) :: params
    integer(c_int) :: lw
    lw = 0
    cusolverDnDsyevj_bufferSize = geps_solver_dsyevj_buf(handle%h, jobz, uplo, n, c_loc(a), lda, &
                                                         c_loc(w), lw, params%p)
    lwork = lw
  end function

  integer function cusolverDnDsyevj(handle, jobz, uplo, n, a, lda, w, work, lwork, devinfo, params)
    type(cusolverDnHandle), intent(in) :: handle
    integer, intent(in) :: jobz, uplo, n, lda, lwork
    real(c_double), intent(inout), target :: a(lda, *), w(*), work(*)
    integer, intent(out), target :: devinfo
    type(cusolverDnSyevjInfo), intent(in) :: params
    cusolverDnDsyevj = geps_solver_dsyevj(handle%h, jobz, uplo, n, c_loc(a), lda, c_loc(w), &
                                          c_loc(work), lwork, c_loc(devinfo), params%p)
  end function

end module cusolverDn
