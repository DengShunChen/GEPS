subroutine correct_unit
   use rank, only: myrank
   use const, only: RTYPE, tmcor, pmcor
   use param, only: lev, jtrun, jtmax
   use index, only: levp, mlist, Llist, mlistnum
   use cudafor
   use openacc
   use mod_vcopy

   implicit none

   integer, parameter :: steps = 3
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax):: vorten, divten, phiten

   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax):: vornow, divnow, temnow, &
                                                        vornow_c, divnow_c, temnow_c
   real(kind=RTYPE), dimension(jtrun, jtmax, 2) :: plnow, plnow_c

   real(kind=RTYPE) :: err(3)

   real(kind=8) tm_1, tm_2, tm_use, mpi_wtime
   logical equal, assert_specspace_close, assert_specspace_2d_close
   integer async_id, s

   async_id = 1
   equal = .true.

   call random_seed()
   call random_number(vorten)
   call random_number(divten)
   call random_number(phiten)

   call random_number(vornow)
   call random_number(divnow)
   call random_number(temnow)
   call random_number(plnow)

   call vcopy_spec(vornow_c, vornow)
   call vcopy_spec(divnow_c, divnow)
   call vcopy_spec(temnow_c, temnow)
   call vcopy_spec(plnow_c, plnow)

   ! << CPU >>
   tm_1 = mpi_wtime()
   call correct(vorten, divten, phiten, tmcor, pmcor, vornow_c, &
                divnow_c, temnow_c, plnow_c, jtrun, jtmax, lev)
   tm_use = mpi_wtime() - tm_1
   if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
      "Elapsed time of CPU: ", tm_use, " seconds"

   ! << GPU >>
   !$acc enter data copyin(lev, jtrun, jtmax, levp, mlistnum, mlist, Llist, &
   !$acc&                  tmcor, pmcor, vorten, divten, phiten) async(async_id)
   !$acc enter data create(vornow, divnow, temnow, plnow) async(async_id)
   do s = 1, steps
      !$acc update device(vornow, divnow, temnow, plnow) async(async_id)
      !$acc wait(async_id)
      tm_1 = mpi_wtime()
      call correct_gpu(vorten, divten, phiten, tmcor, pmcor, vornow, &
                       divnow, temnow, plnow, jtrun, jtmax, lev)

      !$acc wait(async_id)
      tm_use = mpi_wtime() - tm_1
      if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
         "Elapsed time of GPU: ", tm_use, " seconds"
   end do
   !$acc exit data delete(lev, jtrun, jtmax, levp, mlistnum, mlist, Llist, &
   !$acc&                 tmcor, pmcor, vorten, divten, phiten) async(async_id)
   !$acc exit data copyout(vornow, divnow, temnow, plnow) async(async_id)
   !$acc wait(async_id)

#ifdef VERBOSE
   if (myrank .eq. 0) then
      write (*, '(1X, A6, A23, 3A15)') &
         "name", "max|Err|", "L2-Err", "RMSE", "max|Var|"
   end if
   call varerr_spec3D(err, vornow, vornow_c, "vornow")
   call varerr_spec3D(err, divnow, divnow_c, "divnow")
   call varerr_spec3D(err, temnow, temnow_c, "temnow")
#endif

#ifdef SP
   equal = assert_specspace_close(vornow, vornow_c, 1e-4_4, 1e-4_4, "vornow") &
           .and. equal
   equal = assert_specspace_close(divnow, divnow_c, 1e-4_4, 1e-4_4, "divnow") &
           .and. equal
   equal = assert_specspace_close(temnow, temnow_c, 1e-4_4, 1e-4_4, "temnow") &
           .and. equal
   equal = assert_specspace_2d_close(plnow, plnow_c, 1e-4_4, 1e-4_4, "plnow") &
           .and. equal
#else
   equal = assert_specspace_close(vornow, vornow_c, 1e-10, 1e-10, "vornow") &
           .and. equal
   equal = assert_specspace_close(divnow, divnow_c, 1e-10, 1e-10, "divnow") &
           .and. equal
   equal = assert_specspace_close(temnow, temnow_c, 1e-10, 1e-10, "temnow") &
           .and. equal
   equal = assert_specspace_2d_close(plnow, plnow_c, 1e-10, 1e-10, "plnow") &
           .and. equal
#endif
   if (.not. equal) call exit(1)

end subroutine correct_unit

subroutine VarErr_spec3D(Err, a_spec, b_spec, lab)
   use const, only: RTYPE, poly
   use param, only: nx, my_max, my, lev, jtrun, jtmax
   use index, only: nxp, levp, nxdef, jlist1, jlistnum, &
                    nsizey
   use rank, only: myrank
   use mpe
   implicit none
   character(len=6), intent(in) :: lab
   real(kind=RTYPE), intent(out):: Err(3)
   real(kind=RTYPE), intent(in):: A_spec(levp, 2, jtrun, jtmax), &
                                  B_spec(levp, 2, jtrun, jtmax)

   real(kind=RTYPE), dimension(nxp, lev, my_max):: A_grid, B_grid
   real(kind=RTYPE):: cc(nx + 2, levp, my_max), dummy

   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, A_spec, cc, 1, nsizey)
   call ujoinsr(cc, A_grid, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)

   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, B_spec, cc, 1, nsizey)
   call ujoinsr(cc, B_grid, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)

   call VarErr(Err, a_grid, b_grid, lev, 1, lab)

   return
end subroutine VarErr_Spec3D
! ============================================================
subroutine VarErr(Err, a, b, lev, nvar, lab)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxp, nxdef, nxdef_2d, jlist1, jlistnum
   use rank, only: myrank
   use mpe
   implicit none
   ! Measure the difference of between two arrays.
   ! There are the three measurement methods.
   ! 1. Maximum norm
   ! 2. L2-norm (The vertical integral weights need to be corrected)
   ! 3. RMSE
   character(len=*), intent(in):: lab
   real(kind=RTYPE), intent(out):: Err(3, nvar)
   integer, intent(in)::lev, nvar
   real(kind=RTYPE), intent(in):: A(nxp, lev, nvar, my_max), &
                                  B(nxp, lev, nvar, my_max)

   integer i, j, k, n, nxj, jj, pts, lons
   real(kind=8) vamax
   real(kind=RTYPE) Err_loc(nvar, 3), avar_amax(nvar), sum_local(nvar), pi, tmp

   if (octahedral) then
      pts = (20 + nx)*my*lev
   elseif (numreduce == -99) then
      pts = nx*my*lev
   end if

   pi = 4.*atan(1.)
   Err = 0.
   Err_loc = 0.
   avar_amax = 0.

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      lons = nxdef(j)
      sum_local = 0.
      do n = 1, nvar
         do k = 1, lev
            do i = 1, nxj
               avar_amax(n) = max(avar_amax(n), abs(A(i, k, n, jj)))

               tmp = A(i, k, n, jj) - B(i, k, n, jj)

               Err_loc(n, 1) = max(Err_loc(n, 1), abs(tmp))
               sum_local(n) = sum_local(n) + tmp**2
            end do
         end do
         Err_loc(n, 2) = Err_loc(n, 2) + sum_local(n)*(2.*pi/lons)*weight(j)/lev
         Err_loc(n, 3) = Err_loc(n, 3) + sum_local(n)
      end do
   end do

   call mpe_global_max(Err_loc(1, 1), nvar, RTYPE)
   call mpe_global_sum(Err_loc(1, 2), nvar*2, RTYPE)
   call mpe_global_max(avar_amax, nvar, RTYPE)

   do n = 1, nvar
      Err_loc(n, 2) = sqrt(Err_loc(n, 2))
      Err_loc(n, 3) = sqrt(Err_loc(n, 3)/pts/nvar)
   end do

   if (myrank .eq. 0) then
      do n = 1, nvar
         write (*, '(1X, A6, i0.2, 1pe23.15, 3(1pe15.7))') &
            trim(lab), n, Err_loc(n, 1:3), avar_amax(n)
      end do
   end if

   return
end subroutine VarErr

program test_correct
   implicit none

   call mpe_init
   call cons
   call correct_unit
   call mpe_finalize

end program test_correct
