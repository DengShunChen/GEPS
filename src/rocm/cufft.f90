! cuFFT subset on hipFFT. Plans are integer(4) ids (hipfftHandle is a pointer).

module cufft
  use iso_c_binding
  use cudafor, only: cuda_stream_kind
  implicit none
  private
  integer, parameter, public :: CUFFT_SUCCESS = 0
  integer, parameter, public :: CUFFT_R2C = int(z'2a')
  integer, parameter, public :: CUFFT_C2R = int(z'2c')
  integer, parameter, public :: CUFFT_D2Z = int(z'6a')
  integer, parameter, public :: CUFFT_Z2D = int(z'6c')

  public :: cufftCreate, cufftDestroy, cufftSetAutoAllocation
  public :: cufftMakePlanMany, cufftSetStream, cufftSetWorkArea
  public :: cufftExecZ2D, cufftExecD2Z, cufftExecC2R, cufftExecR2C

  interface
    function geps_fft_create(plan) bind(C, name="geps_fft_create")
      import c_int
      integer(c_int) :: plan
      integer(c_int) :: geps_fft_create
    end function
    function geps_fft_destroy(plan) bind(C, name="geps_fft_destroy")
      import c_int
      integer(c_int), value :: plan
      integer(c_int) :: geps_fft_destroy
    end function
    function geps_fft_set_auto_alloc(plan, auto_alloc) bind(C, name="geps_fft_set_auto_alloc")
      import c_int
      integer(c_int), value :: plan, auto_alloc
      integer(c_int) :: geps_fft_set_auto_alloc
    end function
    function geps_fft_make_plan_many(plan, rank, n, inembed, istride, idist, onembed, ostride, &
                                     odist, ffttype, batch, work_size) bind(C, name="geps_fft_make_plan_many")
      import c_int, c_int64_t
      integer(c_int), value :: plan, rank, n, inembed, istride, idist, onembed, ostride, odist, ffttype, batch
      integer(c_int64_t) :: work_size
      integer(c_int) :: geps_fft_make_plan_many
    end function
    function geps_fft_set_stream(plan, stream) bind(C, name="geps_fft_set_stream")
      import c_int, c_int64_t
      integer(c_int), value :: plan
      integer(c_int64_t), value :: stream
      integer(c_int) :: geps_fft_set_stream
    end function
    function geps_fft_set_work_area(plan, work) bind(C, name="geps_fft_set_work_area")
      import c_int, c_ptr
      integer(c_int), value :: plan
      type(c_ptr), value :: work
      integer(c_int) :: geps_fft_set_work_area
    end function
    function geps_fft_exec_z2d(plan, inn, out) bind(C, name="geps_fft_exec_z2d")
      import c_int, c_ptr
      integer(c_int), value :: plan
      type(c_ptr), value :: inn, out
      integer(c_int) :: geps_fft_exec_z2d
    end function
    function geps_fft_exec_d2z(plan, inn, out) bind(C, name="geps_fft_exec_d2z")
      import c_int, c_ptr
      integer(c_int), value :: plan
      type(c_ptr), value :: inn, out
      integer(c_int) :: geps_fft_exec_d2z
    end function
    function geps_fft_exec_c2r(plan, inn, out) bind(C, name="geps_fft_exec_c2r")
      import c_int, c_ptr
      integer(c_int), value :: plan
      type(c_ptr), value :: inn, out
      integer(c_int) :: geps_fft_exec_c2r
    end function
    function geps_fft_exec_r2c(plan, inn, out) bind(C, name="geps_fft_exec_r2c")
      import c_int, c_ptr
      integer(c_int), value :: plan
      type(c_ptr), value :: inn, out
      integer(c_int) :: geps_fft_exec_r2c
    end function
  end interface

contains

  integer function cufftCreate(plan)
    integer(c_int), intent(out) :: plan
    cufftCreate = geps_fft_create(plan)
  end function

  integer function cufftDestroy(plan)
    integer(c_int), intent(in) :: plan
    cufftDestroy = geps_fft_destroy(plan)
  end function

  integer function cufftSetAutoAllocation(plan, auto)
    integer(c_int), intent(in) :: plan, auto
    cufftSetAutoAllocation = geps_fft_set_auto_alloc(plan, auto)
  end function

  integer function cufftMakePlanMany(plan, rank, n, inembed, istride, idist, onembed, ostride, &
                                     odist, ffttype, batch, work_size)
    integer(c_int), intent(in) :: plan, rank, n, inembed, istride, idist, onembed, ostride, odist, ffttype, batch
    integer(c_int64_t), intent(out) :: work_size
    cufftMakePlanMany = geps_fft_make_plan_many(plan, rank, n, inembed, istride, idist, onembed, &
                                                ostride, odist, ffttype, batch, work_size)
  end function

  integer function cufftSetStream(plan, stream)
    integer(c_int), intent(in) :: plan
    integer(cuda_stream_kind), intent(in) :: stream
    cufftSetStream = geps_fft_set_stream(plan, int(stream, c_int64_t))
  end function

  integer function cufftSetWorkArea(plan, workspace)
    integer(c_int), intent(in) :: plan
    integer, intent(inout), target :: workspace(*)
    cufftSetWorkArea = geps_fft_set_work_area(plan, c_loc(workspace))
  end function

  integer function cufftExecZ2D(plan, inn, out)
    integer(c_int), intent(in) :: plan
    real(c_double), intent(in), target :: inn(*)
    real(c_double), intent(out), target :: out(*)
    cufftExecZ2D = geps_fft_exec_z2d(plan, c_loc(inn), c_loc(out))
  end function

  integer function cufftExecD2Z(plan, inn, out)
    integer(c_int), intent(in) :: plan
    real(c_double), intent(in), target :: inn(*)
    real(c_double), intent(out), target :: out(*)
    cufftExecD2Z = geps_fft_exec_d2z(plan, c_loc(inn), c_loc(out))
  end function

  integer function cufftExecC2R(plan, inn, out)
    integer(c_int), intent(in) :: plan
    real(c_float), intent(in), target :: inn(*)
    real(c_float), intent(out), target :: out(*)
    cufftExecC2R = geps_fft_exec_c2r(plan, c_loc(inn), c_loc(out))
  end function

  integer function cufftExecR2C(plan, inn, out)
    integer(c_int), intent(in) :: plan
    real(c_float), intent(in), target :: inn(*)
    real(c_float), intent(out), target :: out(*)
    cufftExecR2C = geps_fft_exec_r2c(plan, c_loc(inn), c_loc(out))
  end function

end module cufft
