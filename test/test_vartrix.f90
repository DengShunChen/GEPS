subroutine vartrix_unit(no, iflag)
   use rank, only: myrank
   use const, only: RTYPE, nnmivm, rad, omega, eigval
   use param, only: lev, jtrun, jtmax
   use index, only: levp, mlist, mlistnum, Llist, Lstart, Lend
   use cudafor
   use openacc
   use mod_vcopy

   implicit none

   integer, parameter :: steps = 3
   integer, no, iflag
   real, dimension(no, 2, 2*jtmax*nnmivm) :: x, x_c
   real x_wrk(no*2, 2)
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax):: vorten, divten, phiten, &
                                                        vorten_c, divten_c, phiten_c

   real, dimension(jtrun, jtrun, nnmivm):: a, b, c
   real, dimension(jtrun, jtmax, nnmivm):: h
   integer nw(jtrun, jtmax)
   integer L, k, mf, m, ns, na, nbig, j, chk_i, chk_j
   logical nnmical

   real(kind=RTYPE), dimension(levp*6*jtrun*jtmax):: cg_buf
   logical cg_created
   type(acc_graph_t) :: cg_graph
   type(cudaGraphExec) :: cg_exec

   real(kind=RTYPE) :: err(3)

   real(kind=8) tm_1, tm_2, tm_use, mpi_wtime
   logical equal, assert_specspace_close
   integer async_id, s

   async_id = 1
   equal = .true.
   cg_created = .false.

   call random_seed()
   call random_number(vorten)
   call random_number(divten)
   call random_number(phiten)
   call random_number(x)
   x_c = x                      !
   vorten_c = vorten
   divten_c = divten
   phiten_c = phiten
   call inicons(rad, omega, eigval, lev, jtrun, jtmax, &
                nw, a, b, c, h, nnmivm)

   ! << CPU >>
   tm_1 = mpi_wtime()
   do L = 1, nnmivm
      if ((L .ge. Lstart) .and. (L .le. Lend)) then
         k = L - Lstart + 1
         if (iflag .eq. +2) &
            call vartran(vorten_c, divten_c, phiten_c, nw, jtrun, jtmax, levp, k, &
                         rad, omega, h(1, 1, L), iflag)

         do m = 1, mlistnum
            mf = mlist(m)
            nbig = jtrun - mf + 1
            !  calculate the size of symmetric and antisymmetric matrix which
            !  include the gravity and rossby wave
            ns = 2*int((nbig + 1)/2) + int(nbig/2)
            na = 2*int(nbig/2) + int((nbig + 1)/2)
            j = 1 + (m - 1)*2 + (L - 1)*2*jtmax
            ! symmetric matrix
            call vartrix(vorten_c, divten_c, phiten_c, levp, &
                         k, x_c(1:ns, 1:2, j), ns, m, nbig, jtrun, jtmax, +1, iflag)
            ! antisymmetric matrix
            call vartrix(vorten_c, divten_c, phiten_c, levp, &
                         k, x_c(1:na, 1:2, j + 1), na, m, nbig, jtrun, jtmax, -1, iflag)

         end do

         if (iflag .eq. -2) &
            call vartran(vorten_c, divten_c, phiten_c, nw, jtrun, jtmax, levp, k, &
                         rad, omega, h(1, 1, L), iflag)
      end if
   end do
   tm_use = mpi_wtime() - tm_1
   if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
      "Elapsed time of CPU: ", tm_use*1e+6, " us"

   ! << GPU >>
   !$acc enter data copyin(lev, jtrun, jtmax, levp, mlistnum, mlist, Llist, &
   !$acc&                  Lstart, Lend, rad, omega, h) async(async_id)
   !$acc enter data copyin(x, vorten, divten, phiten) async(async_id)
   !$acc wait(async_id)
   do s = 1, steps
      if (s .gt. 1) then
         if (iflag .eq. 2) then
            !$acc update device(vorten, divten, phiten) async(async_id)
         else if (iflag .eq. -2) then
            !$acc update device(x) async(async_id)
         end if
         !$acc wait(async_id)
      end if
      tm_1 = mpi_wtime()
      call vartrix_gpu(vorten, divten, phiten, levp, x, no, &
                       rad, omega, h, iflag)

      !$acc wait(async_id)
      tm_use = mpi_wtime() - tm_1
      if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
         "Elapsed time of GPU: ", tm_use*1e+6, " us"
   end do
   !$acc exit data delete(lev, jtrun, jtmax, levp, mlistnum, mlist, Llist, &
   !$acc&                 Lstart, Lend, rad, omega, h) async(async_id)
   !$acc exit data copyout(x, vorten, divten, phiten) async(async_id)
   !$acc wait(async_id)

   if (iflag .eq. +2) then
      call assert_allclose_r8(x, size(x), x_c, size(x_c), 1e-10, 1e-10, "x")
   else if (iflag .eq. -2) then
#ifdef VERBOSE
      if (myrank .eq. 0) then
         write (*, '(1X, A6, A23, 3A15)') &
            "name", "max|Err|", "L2-Err", "RMSE", "max|Var|"
      end if
      call varerr_spec(err, vorten, vorten_c, "vorten")
      call varerr_spec(err, divten, divten_c, "divten")
      call varerr_spec(err, phiten, phiten_c, "phiten")
#endif
#ifdef SP
      equal = assert_specspace_close(vorten, vorten_c, 1e-4_4, 1e-4_4, "vorten") &
              .and. equal
      equal = assert_specspace_close(divten, divten_c, 1e-4_4, 1e-4_4, "divten") &
              .and. equal
      equal = assert_specspace_close(phiten, phiten_c, 1e-4_4, 1e-4_4, "phiten") &
              .and. equal
#else
      equal = assert_specspace_close(vorten, vorten_c, 1e-10, 1e-10, "vorten") &
              .and. equal
      equal = assert_specspace_close(divten, divten_c, 1e-10, 1e-10, "divten") &
              .and. equal
      equal = assert_specspace_close(phiten, phiten_c, 1e-10, 1e-10, "phiten") &
              .and. equal
#endif
      if (.not. equal) call exit(1)
   end if

end subroutine vartrix_unit

subroutine VarErr_spec(Err, a_spec, b_spec, lab)
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
end subroutine VarErr_Spec
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

program test_vartrix
   use rank, only: myrank
   use param, only: jtrun
   implicit none
   integer no

   call mpe_init
   call cons
   no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
   if (myrank .eq. 0) print *, "forward"
   call vartrix_unit(no, +2)
   if (myrank .eq. 0) print *, "backward"
   call vartrix_unit(no, -2)
   call mpe_finalize

end program test_vartrix
