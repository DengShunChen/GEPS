subroutine ptotc_gpu(pdrym, lprint)
   ! Present on device: jlist1, nxdef_2d, dsigma, pt, qt, nxjlen_all, cosl, jlist2, nxdef
   use index
   use rank
   use const
   use param
   use grid
   use openacc
   use cudafor

   implicit none

   integer i, j, k, n, jj, kk, nxj, nxjf, kn
   logical lprint
   real sumtot, sumwat, dsigp &
      , sumtotm, sumwatm, pdrym &
      , qtot
   real(kind=RTYPE), dimension(nxp, my_max) :: sumtotp, sumwatp
   real(kind=RTYPE) :: sumtotp_t, sumwatp_t
   real(kind=RTYPE), dimension(nx, my_max) :: sumtottx, sumwattx
   real(kind=RTYPE), dimension(my_max) :: sumtotpy, sumwatpy
   real(kind=RTYPE) :: sumtotpy_t, sumwatpy_t
   real(kind=RTYPE), dimension(my) :: sumtotty, sumwatty
   integer :: async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(sumtotp, sumwatp, sumtottx, sumwattx, sumtotpy, sumwatpy, sumtotty, sumwatty) async(async_id)
   !$acc parallel loop collapse(2) private(j, nxj, dsigp, qtot, kk, sumtotp_t, sumwatp_t) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nxp
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         if (i .le. nxj) then
            sumtotp_t = 0.
            sumwatp_t = 0.
            do k = 1, lev
               dsigp = dsigma(k, 1)*pt(i, jj) + dsigma(k, 2)
               qtot = 0.
               do n = 1, ncld
                  kk = k + (n - 1)*lev
                  qtot = qtot + qt(i, kk, jj)
               end do
               sumtotp_t = sumtotp_t + dsigp
               sumwatp_t = sumwatp_t + dsigp*qtot
            end do
            sumtotp(i, jj) = sumtotp_t
            sumwatp(i, jj) = sumwatp_t
         end if
      end do
   end do

   call mpe2d_unify_nx_gpu(sumtottx, sumtotp)
   call mpe2d_unify_nx_gpu(sumwattx, sumwatp)
   !$acc parallel loop gang private(sumtotpy_t, sumwatpy_t) async(async_id)
   do jj = 1, jlistnum
      sumtotpy_t = 0.
      sumwatpy_t = 0.
      !$acc loop worker reduction(+:sumtotpy_t) reduction(+:sumwatpy_t)
      do i = 1, nx
         j = jlist1(jj)
         nxjf = nxdef(j)
         if (i .le. nxjf) then
            sumtotpy_t = sumtotpy_t + sumtottx(i, jj)*cosl(j)*nx/nxjf
            sumwatpy_t = sumwatpy_t + sumwattx(i, jj)*cosl(j)*nx/nxjf
         end if
      end do
      sumtotpy(jj) = sumtotpy_t
      sumwatpy(jj) = sumwatpy_t
   end do

   call mpe2d_unify_my1d_gpu(sumtotty, sumtotpy)
   call mpe2d_unify_my1d_gpu(sumwatty, sumwatpy)
   !$acc exit data delete(sumtotp, sumwatp, sumtottx, sumwattx, sumtotpy, sumwatpy) copyout(sumtotty, sumwatty) async(async_id)
   !$acc wait(async_id)

   sumtot = 0.
   sumwat = 0.
   do j = 1, my
      sumtot = sumtot + sumtotty(j)
      sumwat = sumwat + sumwatty(j)
   end do

   kn = nx*my
   sumtotm = sumtot/float(kn)
   sumwatm = sumwat/float(kn)
   pdrym = sumtotm - sumwatm

   if (myrank .eq. 0 .and. lprint) then
      open (35, file='pdry.txt', form='formatted', status='unknown', &
            position='append')
      write (35, *) pdrym, sumwatm, sumtotm
      close (35)
      print *, 'dry air mass = ', pdrym, ' hPa'
      print *, 'water mass = ', sumwatm, ' hPa'
      print *, 'total air mass = ', sumtotm, ' hPa'
   end if

   return
end
