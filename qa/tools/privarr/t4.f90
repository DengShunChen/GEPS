      subroutine work(n, m, nsoil, out, nbad)
      implicit none
      integer, intent(in) :: n, m, nsoil
      real(8), intent(out) :: out(n, m)
      integer, intent(out) :: nbad
      real(8) :: a1(nsoil), a2(nsoil), a3(nsoil), a4(nsoil), a5(nsoil), a6(nsoil), a7(nsoil), a8(nsoil)
      real(8) :: s
      integer :: i, j, k, it
      nbad = 0
      !$omp target teams distribute parallel do collapse(2) private(i, j, k, it, s, a1, a2, a3, a4, a5, a6, a7, a8) map(from:out) map(tofrom:nbad)
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
            do it = 1, 50 + mod(i*7 + j*13, 97)
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

      program p
      implicit none
      integer, parameter :: n = 1552, m = 384
      real(8), allocatable :: out(:, :)
      integer :: nbad
      allocate(out(n, m))
      call work(n, m, 4, out, nbad)
      print *, '20 private automatic arrays: nbad =', nbad, ' sample', out(497, 210)
      end program
