! traced kernel duration vs wall: uneven work across workgroups
program lb9
   use omp_lib
   implicit none
   integer, parameter :: n = 311296, nrep = 50
   real(8), allocatable :: a(:)
   integer :: i, r, kk, iters, mode
   real(8) :: t0, t1, x
   allocate(a(n)); a = 1.0d0
   !$omp target enter data map(to:a)
   do mode = 1, 3
      do r = 1, nrep + 5
         if (r == 6) t0 = omp_get_wtime()
         !$omp target teams distribute parallel do private(x, kk, iters)
         do i = 1, n
            select case (mode)
            case (1); iters = 2000                          ! uniform
            case (2); iters = merge(16000, 0, i <= n/8)      ! first eighth heavy
            case (3); iters = merge(16000, 0, mod((i-1)/256, 8) == 0)  ! every 8th workgroup heavy
            end select
            x = a(i)
            do kk = 1, iters
               x = x*0.999999d0 + 1.0d-7
            end do
            a(i) = x
         end do
      end do
      t1 = omp_get_wtime()
      print '(a,i2,a,f9.1,a)', 'mode ', mode, ' wall per launch ', (t1-t0)/nrep*1e6, ' us'
   end do
end program lb9
