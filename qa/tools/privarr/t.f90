      subroutine work(n, m, nsoil, out, nbad)
      implicit none
      integer, intent(in) :: n, m, nsoil
      real(8), intent(out) :: out(n, m)
      integer, intent(out) :: nbad
      integer, parameter :: nsold = 4
      real(8) :: fixed(nsold), auto(nsoil)
      real(8) :: s
      integer :: i, j, k, it
      nbad = 0
      !$omp target teams distribute parallel do collapse(2) private(i, j, k, it, s, fixed, auto) map(from:out) map(tofrom:nbad)
      do j = 1, m
         do i = 1, n
            do k = 1, nsoil
               auto(k) = dble(i) + 1000.d0*dble(j) + 0.001d0*k
            end do
            do k = 1, nsold
               fixed(k) = -(dble(i) + 1000.d0*dble(j)) + 0.001d0*k
            end do
            ! spin a bit so neighbouring lanes get a chance to interleave
            s = 0.d0
            do it = 1, 50 + mod(i*7 + j*13, 97)
               s = s + sqrt(dble(it) + auto(1 + mod(it, nsoil)))
            end do
            ! check that the private arrays still hold this column's values
            do k = 1, nsoil
               if (auto(k) /= dble(i) + 1000.d0*dble(j) + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
            end do
            do k = 1, nsold
               if (fixed(k) /= -(dble(i) + 1000.d0*dble(j)) + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
            end do
            out(i, j) = s + auto(2) + fixed(3)
         end do
      end do
      end subroutine

      program p
      implicit none
      integer, parameter :: n = 1552, m = 384
      real(8), allocatable :: out(:, :)
      integer :: nbad
      allocate(out(n, m))
      call work(n, m, 4, out, nbad)
      print *, 'private-array isolation test: nbad =', nbad, ' sample', out(497, 210)
      end program
