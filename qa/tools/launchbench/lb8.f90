! completion latency vs kernel duration
program lb8
   use omp_lib
   implicit none
   integer, parameter :: n = 27648, nrep = 200
   real(8), allocatable :: a(:)
   integer :: i, r, k, iters, kk
   real(8) :: t0, t1, x, tk
   character(len=32) :: arg
   allocate(a(n)); a = 1.0d0
   !$omp target enter data map(to:a)
   do k = 1, 6
      iters = 2**(k+5) * 100
      ! warm
      !$omp target teams distribute parallel do private(x, i, kk)
      do i = 1, n
         x = a(i)
         do kk = 1, iters
            x = x*0.999999d0 + 1.0d-7
         end do
         a(i) = x
      end do
      t0 = omp_get_wtime()
      do r = 1, nrep
         !$omp target teams distribute parallel do private(x, i, kk)
         do i = 1, n
            x = a(i)
            do kk = 1, iters
               x = x*0.999999d0 + 1.0d-7
            end do
            a(i) = x
         end do
      end do
      t1 = omp_get_wtime()
      print '(a,i8,a,f9.1,a)', 'iters ', iters, ' wall per launch ', (t1-t0)/nrep*1e6, ' us'
   end do
end program lb8
