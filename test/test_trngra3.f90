!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_trngra3
   use param
   use const
   
   implicit none
   
   call mpe_init
   call cons
   call trngra3_unit
   call mpe_finalize

end program

subroutine trngra3_unit
   use param
   use const, only: RTYPE
   use index
   
   implicit none
   
   real(kind=RTYPE) cim(jtmax)
   real(kind=RTYPE) poly(jtrun, my / 2, jtmax), dpoly(jtrun, my / 2, jtmax)
   real(kind=RTYPE) s(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) dlpl(nxp, levF, my_max), dtpl(nxp, levF, my_max)
   real(kind=RTYPE) dlpl_gpu(nxp, levF, my_max), dtpl_gpu(nxp, levF, my_max)
   
   call random_seed()
   call random_number(cim)
   call random_number(poly)
   call random_number(dpoly)
   call random_number(s)
   
   dlpl = 0.
   dtpl = 0.
   dlpl_gpu = 0.
   dtpl_gpu = 0.
   
   call trngra3(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsizey)
   call trngra3_gpu(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu, nsizey)

   if (all(abs(dlpl - dlpl_gpu) <= 1e-10) .AND. all(abs(dtpl - dtpl_gpu) <= 1e-10)) then
      PRINT *, "test_trngra3 passed."
   else
      PRINT *, "test_trngra3 failed."
      call exit(1)
   end if
   
end subroutine trngra3_unit
