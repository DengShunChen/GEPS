! PoC 06 - runtime-sized private (automatic) arrays in a target region are
! silently corrupted across lanes unless LIBOMPTARGET_STACK_SIZE is raised.
!
! Each lane fills 8 private arrays with a lane-unique pattern, does some
! data-dependent work, then verifies its own arrays.  On a correct
! implementation nbad is always 0.  With the default device stack a large,
! run-to-run varying fraction of the elements has been overwritten by other
! lanes, and nothing reports an error.
!
!   -DFIXED : same 8 arrays with compile-time size 4 -> always correct, so the
!             bug is specific to runtime-sized (automatic) private arrays.
!
! Note the arrays only occupy 8*4*8 = 256 bytes per lane, far below any
! plausible default stack, and one 2-D automatic array a(nsoil,20) of larger
! total size does NOT show the problem - so this is per-array, not per-byte.
subroutine work(n, m, nsoil, out, nbad)
   implicit none
   integer, intent(in) :: n, m, nsoil
   real(8), intent(out) :: out(n, m)
   integer, intent(out) :: nbad
#ifdef FIXED
   real(8) :: a1(4), a2(4), a3(4), a4(4), a5(4), a6(4), a7(4), a8(4)
#else
   real(8) :: a1(nsoil), a2(nsoil), a3(nsoil), a4(nsoil), a5(nsoil), a6(nsoil), a7(nsoil), a8(nsoil)
#endif
   real(8) :: s
   integer :: i, j, k, it
   nbad = 0
   !$omp target teams distribute parallel do collapse(2) private(i, j, k, it, s, a1, a2, a3, a4, a5, a6, a7, a8) &
   !$omp&   map(from:out) map(tofrom:nbad)
   do j = 1, m
      do i = 1, n
         do k = 1, nsoil
            a1(k) = dble(i) + 1000.d0*dble(j) + 1.d0 + 0.001d0*k
            a2(k) = dble(i) + 1000.d0*dble(j) + 2.d0 + 0.001d0*k
            a3(k) = dble(i) + 1000.d0*dble(j) + 3.d0 + 0.001d0*k
            a4(k) = dble(i) + 1000.d0*dble(j) + 4.d0 + 0.001d0*k
            a5(k) = dble(i) + 1000.d0*dble(j) + 5.d0 + 0.001d0*k
            a6(k) = dble(i) + 1000.d0*dble(j) + 6.d0 + 0.001d0*k
            a7(k) = dble(i) + 1000.d0*dble(j) + 7.d0 + 0.001d0*k
            a8(k) = dble(i) + 1000.d0*dble(j) + 8.d0 + 0.001d0*k
         end do
         s = 0.d0
         do it = 1, 50 + mod(i*7 + j*13, 97)          ! data-dependent trip count
            s = s + sqrt(dble(it) + a1(1 + mod(it, nsoil)) + a2(1 + mod(it, nsoil)) + a3(1 + mod(it, nsoil)) + a4(1 + mod(it, nsoil)) + a5(1 + mod(it, nsoil)) + a6(1 + mod(it, nsoil)) + a7(1 + mod(it, nsoil)) + a8(1 + mod(it, nsoil)))
         end do
         do k = 1, nsoil
            if (a1(k) /= dble(i) + 1000.d0*dble(j) + 1.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a2(k) /= dble(i) + 1000.d0*dble(j) + 2.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a3(k) /= dble(i) + 1000.d0*dble(j) + 3.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a4(k) /= dble(i) + 1000.d0*dble(j) + 4.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a5(k) /= dble(i) + 1000.d0*dble(j) + 5.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a6(k) /= dble(i) + 1000.d0*dble(j) + 6.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a7(k) /= dble(i) + 1000.d0*dble(j) + 7.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
            if (a8(k) /= dble(i) + 1000.d0*dble(j) + 8.d0 + 0.001d0*k) then
               !$omp atomic
               nbad = nbad + 1
            end if
         end do
         out(i, j) = s + a1(2) + a8(3)
      end do
   end do
end subroutine

program privarr
   implicit none
   integer, parameter :: n = 1552, m = 384, nsoil = 4
   real(8), allocatable :: out(:, :)
   integer :: nbad
   character(len=64) :: env
   allocate(out(n, m))
   call work(n, m, nsoil, out, nbad)
   call get_environment_variable('LIBOMPTARGET_STACK_SIZE', env)
   if (len_trim(env) == 0) env = '(default)'
#ifdef FIXED
   print '(a,a,a,i0,a,i0,a,es22.15)', 'fixed-size a1..a8(4)     LIBOMPTARGET_STACK_SIZE=', trim(env), &
      ': nbad = ', nbad, ' of ', n*m*nsoil*8, '  out(497,210) = ', out(497, 210)
#else
   print '(a,a,a,i0,a,i0,a,es22.15)', 'automatic a1..a8(nsoil)  LIBOMPTARGET_STACK_SIZE=', trim(env), &
      ': nbad = ', nbad, ' of ', n*m*nsoil*8, '  out(497,210) = ', out(497, 210)
#endif
   if (nbad == 0) then
      print '(a)', 'RESULT: PASS'
   else
      print '(a)', 'RESULT: FAIL (private arrays overwritten by other lanes)'
   end if
end program
