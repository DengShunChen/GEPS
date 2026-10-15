program t
  use iso_c_binding
  implicit none
  interface
    function geps_free() bind(C, name="freegb")
      import c_double
      real(c_double) :: geps_free
    end function
  end interface
  real(8), allocatable :: a(:)
  integer :: it, i
  integer(8) :: n
  print *, 'start free GB =', geps_free()
  do it = 1, 8
     n = 200000000_8 + it * 7000000_8        ! 1.6 .. 2.0 GB, all different sizes
     allocate(a(n)); a = 0.0d0
     !$omp target enter data map(to:a)
     !$omp target teams distribute parallel do
     do i = 1, n
        a(i) = a(i) + 1.0d0
     end do
     print *, '   mapped: free GB =', geps_free()
     !$omp target exit data map(release:a)
     deallocate(a)
     print *, 'it', it, ' size GB=', real(n*8)/1073741824.0, ' free GB =', geps_free()
  end do
end program
