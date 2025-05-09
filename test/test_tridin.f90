program test_moninedmf
   implicit none
   call mpe_init
   call cons
   call tridin_unit
   call mpe_finalize

end program

subroutine assert_real(actual, n_actual, desired, n_desired, rtol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff
   logical :: equal
   integer :: i

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   rel_diff = 0.0
   abs_diff = 0.0

   rel_diff = maxval(abs((actual - desired)/(desired + eps)))
   abs_diff = maxval(abs(actual - desired))
   !if (abs_diff > atol) equal = .false.
   if (rel_diff > rtol) equal = .false.

   if (.not. equal) then
      print *, "Arrays are not close within tolerance rtol =", rtol
      print *, "Max relative difference = ", rel_diff
      !print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_real

subroutine tridin_unit
   use index, only: nxjp, nxp, jlistnum, jlist1, nxptot, nxjp_acc
   use param, only: lev, my_max, my
   use machine, only: kind_phys
   use const, only: RTYPE
   use rank, only: myrank
   use cusparse
   use cudafor
   use openacc
   !USE, INTRINSIC :: IEEE_ARITHMETIC
   implicit none

   integer, parameter :: ncld = 7
   real(kind=kind_phys) au(nxp, lev - 1, my_max), au_gpu(nxp, lev - 1, my_max)
   real(kind=kind_phys) a1(nxp, lev, my_max), a1_gpu(nxp, lev, my_max)
   real(kind=kind_phys) a2(nxp, lev, my_max, ncld), a2_gpu(nxp, lev, my_max, ncld)
   real(kind=kind_phys) al(nxp, lev - 1, my_max), ad(nxp, lev, my_max)
   integer :: ntrac, myim(my_max), ix, km
   integer :: i, ii, j, jj, async_id, n, k, kk, im, j1, istat, accui, m
   integer, dimension(34) :: seed
   real(kind=kind_phys) dt2, num
   real(kind=kind_phys), allocatable, dimension(:, :) :: alc, adc, auc, a1c
   real(kind=kind_phys), allocatable, dimension(:, :, :) :: a2c
   real(kind=kind_phys) :: alg(nxptot, lev), adg(nxptot, lev), &
                           aug(nxptot, lev), a1g(nxptot, lev), &
                           a2g(nxptot, lev, ncld), a3g(nxptot, lev, ncld+1)
   type(cusparseHandle) :: handle
   integer(8), value :: buffer_size
   character(1), allocatable :: pbuffer(:)
   integer(kind=cuda_stream_kind) stream
   real :: time1, time2, ct, gt, diff, tt1, tt2
   ix = nxp
   km = lev

   async_id = -1
   stream = acc_get_cuda_stream(async_id)
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   call random_seed()
   !call random_seed(put=seed)
   call random_number(a1)
   call random_number(a2)

   if (.true.) then
      al = al*0.1 + 0.25
      ad = ad*0.1 + 1.2
      au = au*0.1 + 0.4
      !al = 0.25
      !ad = 1.2
      !au = 0.4
      a1 = a1*100.-30.
      a2 = a2*100.-50.
   end if

   au_gpu = au
   a1_gpu = a1
   a2_gpu = a2

   ct = 0.
   gt = 0.
   do ii = 1, 2
      do jj = 1, jlistnum
         j = jlist1(jj)
         im = nxjp(j)
         allocate (alc(im, lev - 1))
         allocate (adc(im, lev))
         allocate (auc(im, lev - 1))
         allocate (a1c(im, lev))
         allocate (a2c(im, lev, ncld))

         do i = 1, im
            do k = 1, lev
               adc(i, k) = ad(i, k, jj)
               a1c(i, k) = a1(i, k, jj)
               do m = 1, ncld
                  a2c(i, k, m) = a2(i, k, jj, m)
               end do
            end do
            do k = 1, lev - 1
               alc(i, k) = al(i, k, jj)
               auc(i, k) = au(i, k, jj)
            end do

         end do

         call cpu_time(time1)
         call tridin(im, lev, ncld, alc, adc, auc, a1c, a2c, auc, a1c, a2c)
         call cpu_time(time2)

         do k = 1, lev
            do i = 1, im
               a1(i, k, jj) = a1c(i, k)
               do m = 1, ncld
                  a2(i, k, jj, m) = a2c(i, k, m)
               end do
            end do
         end do

         deallocate (alc)
         deallocate (adc)
         deallocate (auc)
         deallocate (a1c)
         deallocate (a2c)

         if (ii .ne. 1) ct = ct + time2 - time1
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      istat = cusparseCreate(handle)
      do jj = 1, jlistnum
         j1 = jlist1(jj)
         myim(jj) = nxjp(j1)
      end do
      
      !$acc enter data copyin(ad,al,au_gpu,a1_gpu,a2_gpu,jlist1,nxjp,nxjp_acc,myim)
      !$acc enter data create(alg,adg,aug,a3g)
      call cpu_time(time1)
      
!      call set_tridiagonal(al, ad, au_gpu, a1_gpu, alg, adg, aug, a1g)
!      call cusparseDgtsvInterleavedBatch_async(handle, 0, lev, alg, adg, aug, a1g, nxptot, -1)
!      !$acc wait(1)
!      call set_tridiagonal(al, ad, au_gpu, a2_gpu, alg, adg, aug, a2g)
!      call cusparseDgtsvInterleavedBatch_async(handle, 0, lev, alg, adg, aug, a2g, nxptot, -1)
!      !$acc wait(1)
!
!      !$acc parallel loop gang collapse(3) private(jj,k,i,accui)
!      do jj = 1, jlistnum
!         do k = 1, lev
!            do i = 1, nxp
!               if (i .le. myim(jj)) then
!                  accui = nxjp_acc(jj) - 1
!                  a1_gpu(i, k, jj) = a1g(accui + i, k)
!                  a2_gpu(i, k, jj) = a2g(accui + i, k)
!               end if
!            end do
!         end do
!      end do
      
      !$acc parallel loop gang vector collapse(3) private(jj,i,k,accui)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  accui = nxjp_acc(jj) - 1
                  if (k .eq. 1) then
                     alg(accui + i, k) = 0.
                  else
                     alg(accui + i, k) = al(i, k - 1, jj)
                  end if
                  if (k .eq. lev) then
                     aug(accui + i, k) = 0.
                  else
                     aug(accui + i, k) = au_gpu(i, k, jj)
                  end if
                  
                  adg(accui + i, k) = ad(i, k, jj)
                  a3g(accui + i, k, 1) = a1_gpu(i, k, jj)
                  do m = 1, ncld
                     a3g(accui + i, k, m+1) = a2_gpu(i, k, jj, m)
                  end do
               end if
            end do
         end do
      end do
      call tridin_gpu(nxptot, km, ncld, alg, adg, aug, a3g, aug, a3g, -1)
      !$acc parallel loop gang collapse(3) private(jj,k,i,accui)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  accui = nxjp_acc(jj) - 1
                  a1_gpu(i, k, jj) = a3g(accui + i, k, 1)
                  do m = 1, ncld
                     a2_gpu(i, k, jj, m) = a3g(accui + i, k, m+1)
                  end do
               end if
            end do
         end do
      end do


      call cpu_time(time2)
      !$acc exit data delete(alg,adg,aug,a3g)
      !$acc exit data copyout(a1_gpu,a2_gpu) &
      !$acc&          delete(au_gpu,al,ad,myim,jlist1,nxjp,nxjp_acc)
!      istat = cusparseDestroy(handle)
      if (ii .ne. 1) gt = gt + time2 - time1
   end do

   call assert_real(a1_gpu, size(a1_gpu), a1, size(a1), &
                    1e-8, "Array a1")
   call assert_real(a2_gpu, size(a2_gpu), a2, size(a2), &
                    1e-8, "Array a2")
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

contains
   subroutine set_tridiagonal(dli, di, dui, bi, dl, d, du, b)
      real(kind=kind_phys) :: dl(nxptot, km), d(nxptot, km), &
                              du(nxptot, km), b(nxptot, km)
      real(kind=kind_phys) :: dli(ix, km - 1, my_max), di(ix, km, my_max), &
                              dui(ix, km - 1, my_max), bi(ix, km, my_max)
      integer :: accui, jj, k, i

      !$acc parallel loop gang vector collapse(3) private(jj,i,k,accui)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  accui = nxjp_acc(jj) - 1
                  d(accui + i, k) = di(i, k, jj)
                  b(accui + i, k) = bi(i, k, jj)
                  if (k .eq. 1) then
                     dl(accui + i, k) = 0.
                  else
                     dl(accui + i, k) = dli(i, k - 1, jj)
                  end if
                  if (k .eq. lev) then
                     du(accui + i, k) = 0.
                  else
                     du(accui + i, k) = dui(i, k, jj)
                  end if
               end if
            end do
         end do
      end do
   end subroutine set_tridiagonal

end subroutine tridin_unit

