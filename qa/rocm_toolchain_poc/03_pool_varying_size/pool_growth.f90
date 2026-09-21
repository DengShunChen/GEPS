! PoC 03 - libomptarget's memory manager only reuses blocks of the *same*
! size.  A program that maps a sequence of large arrays whose sizes differ
! slightly (a very common pattern: work arrays sized by the current block)
! never gets memory back: every exit data leaves the block in the pool,
! every new size allocates a new block, and free VRAM decreases
! monotonically until allocation fails.
!
! With LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0 (pool disabled) free VRAM
! returns to its starting value after every exit data.
program pool_growth
   use iso_c_binding
   implicit none
   interface
      function freegb() bind(C, name="freegb")
         import c_double
         real(c_double) :: freegb
      end function
   end interface
   real(8), allocatable :: a(:)
   integer :: it, niter
   integer(8) :: i, n
   real(8) :: f0, f1, fmin
   character(len=16) :: arg
   niter = 16
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg); read(arg, *) niter
   end if
   f0 = freegb()
   print '(a,f8.2,a)', 'start: free VRAM ', f0, ' GB'
   do it = 1, niter
      n = 200000000_8 + it*7000000_8            ! 1.5 .. 2.5 GB, every size distinct
      allocate(a(n)); a = 0.0d0
      !$omp target enter data map(to:a)
      !$omp target teams distribute parallel do
      do i = 1, n
         a(i) = a(i) + 1.0d0
      end do
      fmin = freegb()
      !$omp target exit data map(release:a)
      deallocate(a)
      f1 = freegb()
      print '(a,i3,a,f6.2,a,f8.2,a,f8.2,a,f8.2,a)', 'iter ', it, ': mapped ', real(n*8)/1073741824.0, &
         ' GB | free while mapped ', fmin, ' GB | free after exit data ', f1, ' GB | not returned ', f0 - f1, ' GB'
   end do
   ! Holding one block (the last size) in a pool is legitimate; holding
   ! more than two of the largest array means the pool never reuses.
   if (f0 - f1 > 2.0d0*real(n*8)/1073741824.0) then
      print '(a,f8.2,a)', 'RESULT: FAIL - ', f0 - f1, ' GB still held after every array was released (pool never reuses)'
   else
      print '(a,f8.2,a)', 'RESULT: PASS - at most one block held (', f0 - f1, ' GB)'
   end if
end program
