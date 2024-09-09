!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_rstrandz
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call rstrandz_unit
   call mpe_finalize

end program

subroutine rstrandz_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: warmup = 8
   integer, parameter :: steps = 24
   real(kind=RTYPE) vdmer(nxp, levf, my_max), vdzon(nxp, levf, my_max)
   real(kind=RTYPE) w(my)
   real(kind=RTYPE) cim(jtmax)
   real(kind=RTYPE) onocos(my)
   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax)
   real(kind=RTYPE) hldten(lev, 2, jtrun, jtmax), vorten(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) hldten_gpu(lev, 2, jtrun, jtmax), vorten_gpu(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max), gwk1(nx + 2, lev, 2, my_max)
   integer i
   integer async_id

   async_id = 1

   call random_seed()
   call random_number(vdmer)
   call random_number(vdzon)
   call random_number(w)
   call random_number(cim)
   call random_number(onocos)
   call random_number(poly)
   call random_number(dpoly)

   hldten = 0.
   vorten = 0.
   hldten_gpu = 0.
   vorten_gpu = 0.

   do i = 1, steps
      call rstrandz(jtrun, jtmax, nx, my, my_max, lev, vdmer, vdzon, w, cim, onocos, poly, dpoly, hldten, vorten, nsizey)
   end do

   !$acc enter data copyin(vdmer, vdzon, w, cim, onocos, poly, dpoly, hldten_gpu, vorten_gpu, &
   !$acc& nlist, jlist2, mtrundef, jlist1, nxjlen, nxjlen_all, nxdef, cc, gwk1) async(async_id)
   do i = 1, steps
     call rstrandz_gpu(jtrun, jtmax, nx, my, my_max, lev, vdmer, vdzon, w, cim, onocos, poly, dpoly, hldten_gpu, vorten_gpu, nsizey, cc, gwk1)
      if (i .eq. warmup .OR. i .eq. steps) then
         !$acc wait(async_id)
      end if
   end do
   !$acc exit data copyout(hldten_gpu, vorten_gpu) delete(vdmer, vdzon, w, cim, onocos, poly, dpoly, &
   !$acc& nlist, jlist2, mtrundef, jlist1, nxjlen, nxjlen_all, nxdef, cc, gwk1) async(async_id)
   !$acc wait(async_id)

   call assert_allclose(hldten_gpu, size(hldten_gpu), hldten, size(hldten), 1e-10, 1e-1, "Array hldten")
   call assert_allclose(vorten_gpu, size(vorten_gpu), vorten, size(vorten), 1e-10, 1e-1, "Array vorten")

end subroutine rstrandz_unit
