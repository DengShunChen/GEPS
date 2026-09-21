! PoC 07 - amdflang: `target ... nowait` fails to compile as soon as the
! region reads any scalar that is not explicitly mapped (implicit
! firstprivate), which is the normal case for a Fortran kernel that uses a
! dummy argument or loop bound:
!
!   error: not yet implemented: Unhandled clause privatization for deferred
!          target tasks in omp.target operation
!
! Variants (all from this file):
!   default      : reads scalar `c` implicitly           -> compile error
!   -DEXPLICIT   : same with firstprivate(c)             -> compile error
!   -DMAPTO      : same with map(to:c)                   -> compiles, correct
!   -DNOSCALAR   : region reads no scalar at all         -> compiles, correct
program nowait_scalar
   implicit none
   integer, parameter :: n = 1000
   real(8) :: a(n), c
   integer :: i
   a = 1.0d0; c = 1.0d0
   !$omp target enter data map(to:a)
#if defined(MAPTO)
   !$omp target teams distribute parallel do nowait map(to:c)
#elif defined(EXPLICIT)
   !$omp target teams distribute parallel do nowait firstprivate(c)
#else
   !$omp target teams distribute parallel do nowait
#endif
   do i = 1, n
#ifdef NOSCALAR
      a(i) = a(i) + 1.0d0
#else
      a(i) = a(i) + c
#endif
   end do
   !$omp taskwait
   !$omp target update from(a)
   print '(a,f4.1,a)', 'a(1) = ', a(1), '  (expected 2.0)'
end program
