program lb6
   implicit none
   integer, parameter :: n1 = 8, n2 = 8, n3 = 8
   real(8), allocatable :: a1(:,:,:), b1(:,:)
   integer :: ot
   allocate(a1(n1,n2,n3), b1(n2,n3)); a1=1; b1=1
   !$omp target enter data map(to:a1,b1)
   !$omp target teams distribute parallel do
   do ot = 1, n3
      b1(1,ot) = b1(1,ot) + a1(1,1,ot)
   end do
   !$omp target teams distribute parallel do
   do ot = 1, n3
      b1(1,ot) = b1(1,ot) + a1(1,1,ot)
   end do
end program lb6
