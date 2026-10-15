      subroutine work(n, m, nsoil, out, nbad)
      implicit none
      integer, intent(in) :: n, m, nsoil
      real(8), intent(out) :: out(n, m)
      integer, intent(out) :: nbad
      real(8) :: a1(nsoil), a2(nsoil), a3(nsoil), a4(nsoil), a5(nsoil), a6(nsoil), a7(nsoil), a8(nsoil), a9(nsoil), a10(nsoil), a11(nsoil), a12(nsoil), a13(nsoil), a14(nsoil), a15(nsoil), a16(nsoil), a17(nsoil), a18(nsoil), a19(nsoil), a20(nsoil)
      real(8) :: s
      integer :: i, j, k, it
      nbad = 0
      !$omp target teams distribute parallel do collapse(2) private(i, j, k, it, s, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20) map(from:out) map(tofrom:nbad)
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
               a9(k) = dble(i) + 1000.d0*dble(j) + 9.d0 + 0.001d0*k
               a10(k) = dble(i) + 1000.d0*dble(j) + 10.d0 + 0.001d0*k
               a11(k) = dble(i) + 1000.d0*dble(j) + 11.d0 + 0.001d0*k
               a12(k) = dble(i) + 1000.d0*dble(j) + 12.d0 + 0.001d0*k
               a13(k) = dble(i) + 1000.d0*dble(j) + 13.d0 + 0.001d0*k
               a14(k) = dble(i) + 1000.d0*dble(j) + 14.d0 + 0.001d0*k
               a15(k) = dble(i) + 1000.d0*dble(j) + 15.d0 + 0.001d0*k
               a16(k) = dble(i) + 1000.d0*dble(j) + 16.d0 + 0.001d0*k
               a17(k) = dble(i) + 1000.d0*dble(j) + 17.d0 + 0.001d0*k
               a18(k) = dble(i) + 1000.d0*dble(j) + 18.d0 + 0.001d0*k
               a19(k) = dble(i) + 1000.d0*dble(j) + 19.d0 + 0.001d0*k
               a20(k) = dble(i) + 1000.d0*dble(j) + 20.d0 + 0.001d0*k
            end do
            s = 0.d0
            do it = 1, 50 + mod(i*7 + j*13, 97)
               s = s + sqrt(dble(it) + a1(1 + mod(it, nsoil)) + a2(1 + mod(it, nsoil)) + a3(1 + mod(it, nsoil)) + a4(1 + mod(it, nsoil)) + a5(1 + mod(it, nsoil)) + a6(1 + mod(it, nsoil)) + a7(1 + mod(it, nsoil)) + a8(1 + mod(it, nsoil)) + a9(1 + mod(it, nsoil)) + a10(1 + mod(it, nsoil)) + a11(1 + mod(it, nsoil)) + a12(1 + mod(it, nsoil)) + a13(1 + mod(it, nsoil)) + a14(1 + mod(it, nsoil)) + a15(1 + mod(it, nsoil)) + a16(1 + mod(it, nsoil)) + a17(1 + mod(it, nsoil)) + a18(1 + mod(it, nsoil)) + a19(1 + mod(it, nsoil)) + a20(1 + mod(it, nsoil)))
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
               if (a9(k) /= dble(i) + 1000.d0*dble(j) + 9.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a10(k) /= dble(i) + 1000.d0*dble(j) + 10.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a11(k) /= dble(i) + 1000.d0*dble(j) + 11.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a12(k) /= dble(i) + 1000.d0*dble(j) + 12.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a13(k) /= dble(i) + 1000.d0*dble(j) + 13.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a14(k) /= dble(i) + 1000.d0*dble(j) + 14.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a15(k) /= dble(i) + 1000.d0*dble(j) + 15.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a16(k) /= dble(i) + 1000.d0*dble(j) + 16.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a17(k) /= dble(i) + 1000.d0*dble(j) + 17.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a18(k) /= dble(i) + 1000.d0*dble(j) + 18.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a19(k) /= dble(i) + 1000.d0*dble(j) + 19.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
               if (a20(k) /= dble(i) + 1000.d0*dble(j) + 20.d0 + 0.001d0*k) then
                  !$omp atomic
                  nbad = nbad + 1
               end if
            end do
            out(i, j) = s + a1(2) + a20(3)
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
