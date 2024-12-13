program test_rayleifr
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call rayleifr_unit
   call mpe_finalize

end program

subroutine rayleifr_unit
   use param
   use const, only: RTYPE, cosl, rad, dt
   use index

   implicit none

   integer, parameter :: warmup = 5
   integer, parameter :: steps = 5
   real(kind=RTYPE) ut(nxp, lev, my_max), vt(nxp, lev, my_max)
   real(kind=RTYPE) ut_c(nxp, lev, my_max), vt_c(nxp, lev, my_max)
   real(kind=RTYPE) ut_g(nxp, lev, my_max), vt_g(nxp, lev, my_max)
   integer i
   integer async_id

   async_id = 1

   call random_seed()
   call random_number(ut)
   call random_number(vt)
   ut = ut*5.5e-8
   vt = vt*5.5e-8

   ut_c = ut
   vt_c = vt
   do i = 1, steps
      call rayleifr(nx, my, my_max, lev, rad, cosl, dt, ut_c, vt_c)
   end do

   ut_g = ut
   vt_g = vt
   !$acc enter data copyin(nx, nxp, my, my_max, lev, rad, cosl, dt, &
   !$acc&                  jlistnum, jlist1, nxdef_2d) async(async_id)
   !$acc enter data copyin(ut_g, vt_g) async(async_id)
   !$acc wait(async_id)
   do i = 1, steps
      call rayleifr_gpu(nx, my, my_max, lev, rad, cosl, dt, ut_g, vt_g)
   end do
   !$acc exit data delete(nx, nxp, my, my_max, lev, rad, cosl, dt, &
   !$acc& jlistnum, jlist1, nxdef_2d) async(async_id)
   !$acc exit data copyout(ut_g, vt_g) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(ut_g, size(ut_g), ut_c, size(ut_c), 1e-4_4, "Array ut_c")
   call assert_rmse(vt_g, size(vt_g), vt_c, size(vt_c), 1e-4_4, "Array vt_c")
#else
   call assert_allclose(ut_g, size(ut_g), ut_c, size(ut_c), 1e-10, 1e-10, "Array ut_c")
   call assert_allclose(vt_g, size(vt_g), vt_c, size(vt_c), 1e-10, 1e-10, "Array vt_c")
#endif

end subroutine rayleifr_unit
