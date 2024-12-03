program test_ndslfv_monoadvv
   implicit none
   logical fgnl

   call mpe_init
   call cons
   fgnl = .false.
   call ndslfv_monoadvv_unit(fgnl, .false.)
   fgnl = .true.
   call ndslfv_monoadvv_unit(fgnl, .true.)
   fgnl = .true.
   call ndslfv_monoadvv_unit(fgnl, .false.)

   call mpe_finalize
end program test_ndslfv_monoadvv

subroutine ndslfv_monoadvv_unit(fgnl, forward)
   use const, only: RTYPE, dt, dsigma
   use rank, only: myrank
   use param, only: nx, my_max, lev, ncld
   use index, only: nxp, levp, levf, myf, jlistnum, jlen, nsizex, row_comm, &
                    nxdef, jlist1, nxptot, nxjp, nxdef_2d, nxjp_acc
   use grid, only: latpart, ndslvvar
   use mod_ndslfv_monoadv_gpu, only: ndslfv_monoadvv_gpu, ndslfv_monoadvv_fgnl_gpu
   use mpe
   implicit none
   logical forward, fgnl
   integer, parameter:: steps = 5
   ! integer, parameter:: steps = 120
   integer, parameter:: async_id = 1
   integer, parameter:: max_nvar = 4
   integer:: nvar

   real(kind=RTYPE):: ut(nxp, lev, my_max), &
                      vt(nxp, lev, my_max), &
                      tt(nxp, lev, my_max), &
                      qt(nxp, lev*ncld, my_max)

   real(kind=RTYPE):: ut_cpu(nxp, lev, my_max), &
                      vt_cpu(nxp, lev, my_max), &
                      tt_cpu(nxp, lev, my_max), &
                      qt_cpu(nxp, lev*ncld, my_max)

   real(kind=RTYPE), allocatable, pinned:: ut_gpu(:, :, :), &
                                           vt_gpu(:, :, :), &
                                           tt_gpu(:, :, :), &
                                           qt_gpu(:, :, :)
   real(kind=RTYPE), allocatable, pinned:: pdot(:, :, :), &
                                           ptm(:, :)

   ! real(kind=RTYPE):: ut_gpu(nxp, lev, my_max), &
   !                    vt_gpu(nxp, lev, my_max), &
   !                    tt_gpu(nxp, lev, my_max), &
   !                    qt_gpu(nxp, lev*ncld, my_max)
   ! real(kind=RTYPE):: pdot(nxp, lev + 1, latpart), &
   !                    ptm(nxp, latpart)

   real(kind=RTYPE) dtah
   integer i, k, j, jj, nxj, n

   real(kind=RTYPE):: err_arr(4, max_nvar), vamax, ummax, vmmax
   character(len=6):: name(max_nvar)

   name = (/'t', 'u', 'v', 'q'/)
   dtah = 0.5*dt

   if (fgnl) then
      nvar = 3
   else
      nvar = 4
   end if

   allocate (ut_gpu(nxp, lev, my_max), &
             vt_gpu(nxp, lev, my_max), &
             tt_gpu(nxp, lev, my_max), &
             qt_gpu(nxp, lev*ncld, my_max), &
             pdot(nxp, lev + 1, latpart), &
             ptm(nxp, latpart))

   if (myrank .eq. 0) then
      print *, "========================================"
      print *, "    start test ndslfv_monoadvv with "
      print *, "    fgnl=", fgnl
      print *, "    forward=", forward
      print *, "========================================"
      write (*, '(1X, A5, 1X, f15.2)') "dt=", dt
      write (*, '(1X, A5, 1X, f15.2)') "dtah=", dtah
   end if

   call random_seed()
   call random_number(ptm)
   call random_number(pdot)
   call random_number(ut)
   call random_number(vt)
   call random_number(tt)
   if (not(fgnl)) call random_number(qt)

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         !! map [0,1] to [500, 1050]
         ptm(i, jj) = 500.+550.*ptm(i, jj)
      end do

      do k = 1, lev + 1
         do i = 1, nxj
            pdot(i, k, jj) = pdot(i, k, jj)/dtah/4.
         end do
      end do
   end do

   !$acc data present_or_copyin(jlistnum,jlist1,nxdef_2d, nxjp, &
   !$acc& nxjp_acc, ncld, lev, nxptot, ndslvvar, dsigma,  &
   !$acc& dtah) &
   !$acc& create(ut_gpu, vt_gpu, tt_gpu, qt_gpu, &
   !$acc& pdot, ptm) &
   !$acc& async(async_id)

   !$acc wait(async_id)
   do n = 1, steps
      if (myrank .eq. 0) write (*, '(" << ", i3, " >>")') n

      ! << CPU >>
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k = 1, lev
            do i = 1, nxj
               tt_cpu(i, k, jj) = tt(i, k, jj)
               ut_cpu(i, k, jj) = ut(i, k, jj)
               vt_cpu(i, k, jj) = vt(i, k, jj)
            end do
         end do
      end do
      if (fgnl) then
         call ndslfv_monoadvv_fgnl(ut_cpu, vt_cpu, tt_cpu, pdot, ptm, &
                                   nxjp, dtah, 3, forward)
      else
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do k = 1, ncld*lev
               do i = 1, nxj
                  qt_cpu(i, k, jj) = qt(i, k, jj)
               end do
            end do
         end do
         call ndslfv_monoadvv(tt_cpu, qt_cpu, ut_cpu, vt_cpu, pdot, ptm, &
                              nxjp, dtah, forward)
      end if

      ! << GPU >>
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k = 1, lev
            do i = 1, nxj
               tt_gpu(i, k, jj) = tt(i, k, jj)
               ut_gpu(i, k, jj) = ut(i, k, jj)
               vt_gpu(i, k, jj) = vt(i, k, jj)
            end do
         end do
      end do
      if (fgnl) then
         !$acc update device(ut_gpu,vt_gpu,tt_gpu)
         call ndslfv_monoadvv_fgnl_gpu(tt_gpu, ut_gpu, vt_gpu, pdot, ptm, &
                                       nxjp, dtah, 3, forward)
         !$acc wait(async_id)
         !$acc update self(tt_gpu,ut_gpu,vt_gpu)
      else
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do k = 1, ncld*lev
               do i = 1, nxj
                  qt_gpu(i, k, jj) = qt(i, k, jj)
               end do
            end do
         end do
         !$acc update device(ut_gpu,vt_gpu,tt_gpu,qt_gpu)
         !$acc update device(ptm, pdot)
         call ndslfv_monoadvv_gpu(tt_gpu, qt_gpu, ut_gpu, vt_gpu, pdot, ptm, &
                                  nxjp, dtah, forward)
         !$acc wait(async_id)
         !$acc update self(tt_gpu,qt_gpu,ut_gpu,vt_gpu)
      end if
   end do
   !$acc end data

   call Varerr(err_arr(1, 1), tt_gpu, tt_cpu, lev, 1)
   call Varerr(err_arr(1, 2), ut_gpu, ut_cpu, lev, 1)
   call Varerr(err_arr(1, 3), vt_gpu, vt_cpu, lev, 1)
   if (not(fgnl)) &
      call Varerr(err_arr(1, 4), qt_gpu, qt_cpu, lev, ncld)

   if (myrank .eq. 0) then
      do k = 1, nvar
         write (*, '(1X, A6, 1pe23.15, 3(1pe15.7))') &
            name(k), err_arr(1:4, k)
      end do
   end if

#ifdef SP
   if (all(err_arr(1, 1:nvar) < 1e-5)) then
#else
   if (all(err_arr(1, 1:nvar) < 1e-10)) then
#endif
      if (myrank .eq. 0) &
         write (*, '(4X, A, 1X, 2(1X,A,L)), A'), &
         "test_ndslfv_monoadvv", &
         "fgnl=", fgnl, &
         "forward=", forward, " passed."
   else
      if (myrank .eq. 0) &
         write (*, '(4X, A, 1X, 2(1X,A,L)), A'), &
         "test_ndslfv_monoadvv", &
         "fgnl=", fgnl, &
         "forward=", forward, " failed."
      call exit(1)
   end if
end subroutine ndslfv_monoadvv_unit

subroutine VarErr(Err, a, b, lev, nvar)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxp, nxdef, jlist1, jlistnum, nxdef_2d
   use rank, only: myrank
   use mpe
   implicit none
   real(kind=RTYPE), intent(out):: Err(4)
   integer, intent(in)::lev, nvar
   real(kind=RTYPE), intent(in):: A(nxp, lev*nvar, my_max), &
                                  B(nxp, lev*nvar, my_max)

   integer i, j, k, nxj, jj, pts
   real(kind=RTYPE) tmp, vamax, sum_local
   if (octahedral) then
      pts = .5*(20 + nx)*my*lev
   elseif (numreduce == -99) then
      pts = nx*my*lev
   end if

   Err = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      sum_local = 0.
      do k = 1, lev*nvar
         do i = 1, nxj
            tmp = A(i, k, jj) - B(i, k, jj)
            Err(1) = max(Err(1), abs(tmp))
            sum_local = sum_local + tmp**2
         end do
      end do
      Err(2) = Err(2) + sum_local/nxj/(lev*nvar)*weight(j)
      Err(3) = Err(3) + sum_local
   end do

   call mpe_global_max(Err(1), 1, RTYPE)
   call mpe_global_sum(Err(2), 2, RTYPE)
   Err(2) = sqrt(Err(2))
   Err(3) = sqrt(Err(3)/pts/nvar)
   Err(4) = vamax(a, lev)

   return
end subroutine VarErr
! ------------------------------------------------------------
real(kind=8) function Vamax(a, lev)
   use const, only: RTYPE
   use param, only: nx, my_max, my
   use index, only: nxp, nxdef_2d, jlist1, jlistnum, levp
   use rank, only: myrank
   use mpe
   implicit none
   integer, intent(in)::lev
   real(kind=RTYPE), intent(in):: A(nxp, lev, my_max)

   integer i, j, k, nxj, jj

   vamax = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do k = 1, lev
         do i = 1, nxj
            vamax = max(abs(A(i, k, jj)), vamax)
         end do
      end do
   end do

   call mpe_global_max(vamax, 1, RTYPE)

   return
end function Vamax
