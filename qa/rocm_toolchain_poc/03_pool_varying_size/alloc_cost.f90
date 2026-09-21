! PoC 03b - the other side of the trade-off.  With the pool disabled
! (LIBOMPTARGET_MEMORY_MANAGER_THRESHOLD=0, the only way to stop PoC 03a
! from running out of VRAM) every enter/exit data pair is a raw
! hsa_amd_memory_pool_allocate / _free, and the free alone costs hundreds
! of microseconds whatever the size.
program alloc_cost
   use omp_lib
   implicit none
   integer, parameter :: nrep = 500
   integer(8), parameter :: sizes(5) = [4096_8, 262144_8, 4194304_8, 67108864_8, 1073741824_8]  ! bytes
   real(8), allocatable :: a(:)
   integer :: k, r
   integer(8) :: n
   real(8) :: t0, tin, tout
   print '(a)', '   size        enter data (alloc)   exit data (free)     [us per call, data uninitialised]'
   do k = 1, size(sizes)
      n = sizes(k)/8
      allocate(a(n))
      do r = 1, 5   ! warm
         !$omp target enter data map(alloc:a)
         !$omp target exit data map(delete:a)
      end do
      tin = 0; tout = 0
      do r = 1, nrep
         t0 = omp_get_wtime()
         !$omp target enter data map(alloc:a)
         tin = tin + (omp_get_wtime() - t0)
         t0 = omp_get_wtime()
         !$omp target exit data map(delete:a)
         tout = tout + (omp_get_wtime() - t0)
      end do
      deallocate(a)
      print '(f9.3,a,2f18.1)', real(sizes(k))/1048576.0, ' MB', tin/nrep*1e6, tout/nrep*1e6
   end do
end program
