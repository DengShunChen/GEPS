program t
   implicit none
   integer, parameter :: lev = 72, jlistnum = 192, nxp = 1536
   real(8), allocatable :: ut(:,:,:), vt(:,:,:), cosl(:)
   integer, allocatable :: jlist1(:), nxdef_2d(:)
   real(8) :: wmax(lev), wref(lev), wt, xx, rad
   integer :: k, jj, i, j, nxj, nbad
   allocate(ut(nxp,lev,jlistnum), vt(nxp,lev,jlistnum), cosl(384), jlist1(jlistnum), nxdef_2d(384))
   call random_number(ut); call random_number(vt)
   do j = 1, 384
      cosl(j) = 0.1d0 + 0.9d0*abs(cos((j-0.5d0)*3.14159d0/384))
      nxdef_2d(j) = 64 + (nxp-64)*j/384
   end do
   do jj = 1, jlistnum
      jlist1(jj) = 2*jj - 1
   end do
   rad = 6371000.d0
   ! reference (host)
   do k = 1, lev
      wt = 0.0
      do jj = 1, jlistnum
         do i = 1, nxp
            j = jlist1(jj); nxj = nxdef_2d(j)
            if (i .le. nxj) then
               xx = rad/cosl(j)
               wt = max(wt, xx*sqrt(ut(i,k,jj)**2 + vt(i,k,jj)**2))
            end if
         end do
      end do
      wref(k) = wt
   end do
   !$omp target enter data map(to:ut, vt, cosl, jlist1, nxdef_2d) map(alloc:wmax)
   !$omp target teams distribute private(wt) firstprivate(j, nxj, xx)
   do k = 1, lev
      wt = 0.0
      !$omp parallel do private(xx, j, nxj) collapse(2) reduction(max:wt)
      do jj = 1, jlistnum
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            if (i .le. nxj) then
               xx = rad/cosl(j)
               wt = max(wt, xx*sqrt(ut(i, k, jj)**2 + vt(i, k, jj)**2))
            end if
         end do
      end do
      wmax(k) = wt
   end do
   !$omp target update from(wmax)
   nbad = 0
   do k = 1, lev
      if (wmax(k) .ne. wref(k)) then
         nbad = nbad + 1
         if (nbad .le. 5) print *, 'k', k, 'gpu', wmax(k), 'ref', wref(k)
      end if
   end do
   print *, 'nested parallel-do max reduction: bad levels =', nbad, ' of', lev
   ! variant: single combined directive with reduction on teams distribute parallel do (per k? no) -> serial-inner baseline
   !$omp target teams distribute parallel do private(wt, j, nxj, xx, jj, i)
   do k = 1, lev
      wt = 0.0
      do jj = 1, jlistnum
         do i = 1, nxp
            j = jlist1(jj); nxj = nxdef_2d(j)
            if (i .le. nxj) then
               xx = rad/cosl(j)
               wt = max(wt, xx*sqrt(ut(i,k,jj)**2 + vt(i,k,jj)**2))
            end if
         end do
      end do
      wmax(k) = wt
   end do
   !$omp target update from(wmax)
   nbad = 0
   do k = 1, lev
      if (wmax(k) .ne. wref(k)) nbad = nbad + 1
   end do
   print *, 'serial-inner baseline: bad levels =', nbad
end program t
