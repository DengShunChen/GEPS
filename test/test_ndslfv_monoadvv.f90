program test_ndslfv_monoadvv
   implicit none

   call mpe_init
   call cons
   call ndslfv_monoadvv_unit(.false.)
   call mpe_finalize
end program test_ndslfv_monoadvv

subroutine ndslfv_monoadvv_unit(forward)
   use const, only: RTYPE, dt, dsigma
   use rank, only: myrank
   use param, only: nx, my_max, lev, ncld
   use index, only: nxp, levp, levf, myf, jlistnum, jlen, nsizex, row_comm, &
                    nxdef, jlist1, nxptot, nxjp, nxdef_2d, nxjp_acc
   use grid, only: latpart, ndslvvar
   use mod_ndslfv_monoadv_gpu, only: ndslfv_monoadvv_gpu
   use mpe
   implicit none
   logical forward
   integer, parameter:: steps = 5
   ! integer, parameter:: steps = 120
   integer, parameter:: async_id = 1
   integer, parameter:: nvar = 4

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

   real(kind=RTYPE):: err_arr(4, nvar), vamax, ummax, vmmax
   character(len=6):: name(nvar)

   name = (/'t', 'u', 'v', 'q'/)
   dtah = 0.5*dt

   allocate (ut_gpu(nxp, lev, my_max), &
             vt_gpu(nxp, lev, my_max), &
             tt_gpu(nxp, lev, my_max), &
             qt_gpu(nxp, lev*ncld, my_max), &
             pdot(nxp, lev + 1, latpart), &
             ptm(nxp, latpart))

   if (myrank .eq. 0) then
      print *, "========================================"
      print *, "    start test ndslfv_monoadvv with "
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
   call random_number(qt)

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
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
         nxj = nxdef(j)
         do k = 1, lev
            do i = 1, nxj
               tt_cpu(i, k, jj) = tt(i, k, jj)
               ut_cpu(i, k, jj) = ut(i, k, jj)
               vt_cpu(i, k, jj) = vt(i, k, jj)
            end do
         end do
         do k = 1, ncld*lev
            do i = 1, nxj
               qt_cpu(i, k, jj) = qt(i, k, jj)
            end do
         end do
      end do
      call ndslfv_monoadvv(tt_cpu, qt_cpu, ut_cpu, vt_cpu, pdot, ptm, &
                           nxjp, dtah, forward)
      ! << GPU >>
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         do k = 1, lev
            do i = 1, nxj
               tt_gpu(i, k, jj) = tt(i, k, jj)
               ut_gpu(i, k, jj) = ut(i, k, jj)
               vt_gpu(i, k, jj) = vt(i, k, jj)
            end do
         end do
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
   end do
   !$acc end data

   call Varerr(err_arr(1, 1), tt_gpu, nxp, tt_cpu, nxp, lev, 1)
   call Varerr(err_arr(1, 2), ut_gpu, nxp, ut_cpu, nxp, lev, 1)
   call Varerr(err_arr(1, 3), vt_gpu, nxp, vt_cpu, nxp, lev, 1)
   call Varerr(err_arr(1, 4), qt_gpu, nxp, qt_cpu, nxp, lev, ncld)

   if (myrank .eq. 0) then
      do k = 1, nvar
         write (*, '(1X, A6, 1pe23.15, 3(1pe15.7))') &
            name(k), err_arr(1:4, k)
      end do
   end if

   if (all(err_arr(1, 1:nvar) < 1e-10)) then
      if (myrank .eq. 0) write (*, '(A,i3,A)') &
         "test_ndslfv_monoadvv (forward=", forward, ") passed."
   else
      if (myrank .eq. 0) write (*, '(A,i3,A)') &
         "test_ndslfv_monoadvv (forward=", forward, ") failed."
      call exit(1)
   end if
end subroutine ndslfv_monoadvv_unit

subroutine VarErr(Err, a, lda, b, ldb, lev, nvar)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxdef, jlist1, jlistnum
   use rank, only: myrank
   use mpe
   implicit none
   real(kind=RTYPE), intent(out):: Err(4)
   integer, intent(in)::lda, ldb, lev, nvar
   real(kind=RTYPE), intent(in):: A(lda, lev*nvar, my_max), &
                                  B(ldb, lev*nvar, my_max)

   integer i, j, k, nxj, jj, pts
   real(kind=RTYPE) tmp, vamax, sum
   if (octahedral) then
      pts = (20 + nx)*my*lev
   elseif (numreduce == -99) then
      pts = nx*my*lev
   end if

   Err = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      sum = 0.
      do k = 1, lev*nvar
         do i = 1, nxj
            tmp = A(i, k, jj) - B(i, k, jj)
            Err(1) = max(Err(1), abs(tmp))
            sum = sum + tmp**2
         end do
      end do
      Err(2) = Err(2) + sum/nxj/(lev*nvar)*weight(j)
      Err(3) = Err(3) + sum
   end do

   call mpe_global_max(Err(1), 1, RTYPE)
   call mpe_global_sum_r8(Err(2), 2, RTYPE)
   Err(2) = sqrt(Err(2))
   Err(3) = sqrt(Err(3)/pts/nvar)
   Err(4) = vamax(a, lda, lev)

   return
end subroutine VarErr
! ------------------------------------------------------------
real(kind=8) function Vamax(a, lda, lev)
   use const, only: RTYPE
   use param, only: nx, my_max, my
   use index, only: nxdef, jlist1, jlistnum, levp
   use rank, only: myrank
   use mpe
   implicit none
   integer, intent(in)::lda, lev
   real(kind=RTYPE), intent(in):: A(lda, lev, my_max)

   integer i, j, k, nxj, jj

   vamax = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do k = 1, lev
         do i = 1, nxj
            vamax = max(abs(A(i, k, jj)), vamax)
         end do
      end do
   end do

   call mpe_global_max(vamax, 1, RTYPE)

   return
end function Vamax
