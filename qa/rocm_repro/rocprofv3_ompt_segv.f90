! Minimal repro: amdflang + OpenMP TARGET OFFLOAD + OpenBLAS, under rocprofv3.
program ompt_repro2
  implicit none
  integer, parameter :: n = 4
  real(kind=8) :: a(n,n), b(n,n), c(n,n)
  integer :: i
  real(kind=8) :: s
  a = 1.0d0; b = 2.0d0; c = 0.0d0
  ! OpenBLAS init happens on this call (gotoblas_init -> omp_get_num_places)
  call dgemm('N','N', n, n, n, 1.0d0, a, n, b, n, 0.0d0, c, n)
  s = 0.0d0
  !$omp target teams distribute parallel do reduction(+:s)
  do i = 1, 1024
     s = s + real(i, 8)
  end do
  print *, 'dgemm ok c(1,1)=', c(1,1), ' target sum=', s
end program
