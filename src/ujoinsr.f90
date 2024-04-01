subroutine ujoinsr(cc, r1, r2, r3, r4, nx, my_max, lev, jlistnum, num, ncld)

   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(*), r1(*), r2(*), r3(*), r4(*)
   integer nx, my_max, lev, jlistnum, num, ncld

   if (num .eq. 1) call ujoin1sr(cc, r1, nx, my_max, lev &
                                 , jlistnum, ncld)
   if (num .eq. 2) call ujoin2sr(cc, r1, r2, nx, my_max, lev &
                                 , jlistnum, ncld)
   if (num .eq. 3) call ujoin3sr(cc, r1, r2, r3, nx, my_max, lev &
                                 , jlistnum, ncld)
   if (num .eq. 4) call ujoin4sr(cc, r1, r2, r3, r4, nx, my_max, lev &
                                 , jlistnum, ncld)

   return
end

subroutine ujoin1sr(cc, r1, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nx, levp, ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n

   bufA = 0.
   bufB = 0.
   r1 = 0.
   do jj = 1, jlistnum
   do n = 1, ncld
   do k = 1, levp
   do i = 1, nx
      bufA(i, k, n, jj) = cc(i, k, n, jj)
   end do
   end do
   end do
   end do

   call mpe2d_transpose_nx_levp(bufA, bufB, nxp, nx, lev, levp, ncld, myf, my_max, jlistnum, jlen, nsizex, row_comm)

   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      do k = 1, lev
         kk = nk + k
         do i = 1, nxp
            r1(i, kk, jj) = bufB(i, k, n, jj)
         end do
      end do
   end do
   end do

   return
end

subroutine ujoin2sr(cc, r1, r2, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 1 + ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r2(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nx, levp, 1 + ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, 1 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n

   bufA = 0.
   bufB = 0.
   r1 = 0.
   do jj = 1, jlistnum
   do n = 1, 1 + ncld
   do k = 1, levp
   do i = 1, nx
      bufA(i, k, n, jj) = cc(i, k, n, jj)
   end do
   end do
   end do
   end do

   call mpe2d_transpose_nx_levp(bufA, bufB, nxp, nx, lev, levp, 1 + ncld, myf, my_max, jlistnum, jlen, nsizex, row_comm)

   do jj = 1, jlistnum
   do k = 1, lev
   do i = 1, nxp
      r1(i, k, jj) = bufB(i, k, 1, jj)
   end do
   end do
   end do

   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      do k = 1, lev
         kk = nk + k
         do i = 1, nxp
            r2(i, kk, jj) = bufB(i, k, 1 + n, jj)
         end do
      end do
   end do
   end do

   return
end

subroutine ujoin3sr(cc, r1, r2, r3, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 2 + ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r2(nxp, lev, my_max)
   real(kind=RTYPE) r3(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nx, levp, 2 + ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, 2 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n

   bufA = 0.
   bufB = 0.
   r1 = 0.

   do jj = 1, jlistnum
   do n = 1, 2 + ncld
   do k = 1, levp
   do i = 1, nx
      bufA(i, k, n, jj) = cc(i, k, n, jj)
   end do
   end do
   end do
   end do

   call mpe2d_transpose_nx_levp(bufA, bufB, nxp, nx, lev, levp, 2 + ncld, myf, my_max, jlistnum, jlen, nsizex, row_comm)

   do jj = 1, jlistnum
   do k = 1, lev
   do i = 1, nxp
      r1(i, k, jj) = bufB(i, k, 1, jj)
      r2(i, k, jj) = bufB(i, k, 2, jj)
   end do
   end do
   end do

   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      do k = 1, lev
         kk = nk + k
         do i = 1, nxp
            r3(i, kk, jj) = bufB(i, k, 2 + n, jj)
         end do
      end do
   end do
   end do

   return
end

subroutine ujoin4sr(cc, r1, r2, r3, r4, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 3 + ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r2(nxp, lev, my_max)
   real(kind=RTYPE) r3(nxp, lev, my_max)
   real(kind=RTYPE) r4(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nx, levp, 3 + ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, 3 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n

   do jj = 1, jlistnum
   do n = 1, 3 + ncld
   do k = 1, levp
   do i = 1, nx
      bufA(i, k, n, jj) = cc(i, k, n, jj)
   end do
   end do
   end do
   end do

   call mpe2d_transpose_nx_levp(bufA, bufB, nxp, nx, lev, levp, 3 + ncld, myf, my_max, jlistnum, jlen, nsizex, row_comm)

   do jj = 1, jlistnum
   do k = 1, lev
   do i = 1, nxp
      r1(i, k, jj) = bufB(i, k, 1, jj)
      r2(i, k, jj) = bufB(i, k, 2, jj)
      r3(i, k, jj) = bufB(i, k, 3, jj)
   end do
   end do
   end do

   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      do k = 1, lev
         kk = nk + k
         do i = 1, nxp
            r4(i, kk, jj) = bufB(i, k, 3 + n, jj)
         end do
      end do
   end do
   end do

   return
end
