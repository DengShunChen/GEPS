program s
   implicit none
   integer, parameter :: lev = 8, n = 1000
   real(8) :: a(n, lev), wsum(lev), wt, ref
   integer :: k, i
   a = 1.0d0
   !$omp target enter data map(to:a) map(alloc:wsum)
   !$omp target teams distribute private(wt)
   do k = 1, lev
      wt = 0.0
      !$omp parallel do reduction(+:wt)
      do i = 1, n
         wt = wt + a(i, k)
      end do
      wsum(k) = wt
   end do
   !$omp target update from(wsum)
   print *, 'sum reduction nested under teams distribute: got', wsum(1), ' expected', real(n,8)
   !$omp target teams distribute parallel do reduction(+:wt)
   do i = 1, n
      wt = wt + a(i, 1)
   end do
   print *, 'flat target teams distribute parallel do reduction: got', wt
end program s
