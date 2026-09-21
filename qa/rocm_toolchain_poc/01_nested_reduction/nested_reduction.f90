! PoC 01 - amdflang: `parallel do reduction(...)` nested inside
! `target teams distribute` produces 0 for every team.
!
! Five variants of the same computation. Only (a) and (b) are wrong.
!   (a) target teams distribute  /  parallel do reduction(+:wt)   -> wrong (0)
!   (b) target teams distribute  /  parallel do reduction(max:wt) -> wrong (0)
!   (c) target teams + distribute split, nested reduction(+:wt)   -> wrong (0)
!   (d) flat: target teams distribute parallel do reduction(+:s)  -> correct
!   (e) same loop nest, inner loop serial (no reduction clause)   -> correct
program nested_reduction
   implicit none
   integer, parameter :: lev = 8, n = 1000
   real(8) :: a(n, lev), wsum(lev), wmax(lev), wt, s
   integer :: k, i, nfail
   nfail = 0
   do k = 1, lev
      do i = 1, n
         a(i, k) = 1.0d0 + 0.001d0*k          ! sum = n*(1+0.001k), max = 1+0.001k
      end do
   end do
   !$omp target enter data map(to:a) map(alloc:wsum, wmax)

   ! (a) nested + reduction
   !$omp target teams distribute private(wt)
   do k = 1, lev
      wt = 0.0d0
      !$omp parallel do reduction(+:wt)
      do i = 1, n
         wt = wt + a(i, k)
      end do
      wsum(k) = wt
   end do
   !$omp target update from(wsum)
   call check('(a) nested  parallel do reduction(+)  ', wsum, [(n*(1.0d0+0.001d0*k), k=1,lev)], nfail)

   ! (b) nested max reduction
   !$omp target teams distribute private(wt)
   do k = 1, lev
      wt = 0.0d0
      !$omp parallel do reduction(max:wt)
      do i = 1, n
         wt = max(wt, a(i, k))
      end do
      wmax(k) = wt
   end do
   !$omp target update from(wmax)
   call check('(b) nested  parallel do reduction(max)', wmax, [(1.0d0+0.001d0*k, k=1,lev)], nfail)

   ! (c) teams / distribute on separate directives
   !$omp target teams private(wt)
   !$omp distribute
   do k = 1, lev
      wt = 0.0d0
      !$omp parallel do reduction(+:wt)
      do i = 1, n
         wt = wt + a(i, k)
      end do
      wsum(k) = wt
   end do
   !$omp end distribute
   !$omp end target teams
   !$omp target update from(wsum)
   call check('(c) teams/distribute split, nested (+)', wsum, [(n*(1.0d0+0.001d0*k), k=1,lev)], nfail)

   ! (d) flat combined directive - correct
   s = 0.0d0
   !$omp target teams distribute parallel do reduction(+:s)
   do i = 1, n
      s = s + a(i, 1)
   end do
   call check('(d) flat teams distribute parallel do ', [s], [n*(1.0d0+0.001d0)], nfail)

   ! (e) nested loop, inner loop serial - correct (so the loop body is fine)
   !$omp target teams distribute parallel do private(wt, i)
   do k = 1, lev
      wt = 0.0d0
      do i = 1, n
         wt = wt + a(i, k)
      end do
      wsum(k) = wt
   end do
   !$omp target update from(wsum)
   call check('(e) inner loop serial, no reduction   ', wsum, [(n*(1.0d0+0.001d0*k), k=1,lev)], nfail)

   if (nfail == 0) then
      print *, 'RESULT: PASS'
   else
      print '(a,i0,a)', 'RESULT: FAIL (', nfail, ' variants wrong)'
   end if
contains
   subroutine check(label, got, want, nfail)
      character(*), intent(in) :: label
      real(8), intent(in) :: got(:), want(:)
      integer, intent(inout) :: nfail
      if (all(abs(got - want) <= 1d-9*abs(want))) then
         print '(2a,es14.6)', label, ' OK    got(1)=', got(1)
      else
         print '(2a,es14.6,a,es14.6)', label, ' WRONG got(1)=', got(1), ' expected ', want(1)
         nfail = nfail + 1
      end if
   end subroutine
end program
