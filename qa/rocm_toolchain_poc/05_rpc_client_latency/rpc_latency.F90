! PoC 05 - a device-side Fortran `print`/`write`/`stop` anywhere in the
! offload image (even on a branch that is never taken) links
! __llvm_rpc_client into the image.  libomptarget then starts an RPC
! server thread for the whole run.  In the application this made the host
! observe kernel completion ~one kernel-duration late for every kernel
! longer than ~1 ms (0.7 s of a 3.5 s time step).
!
! Build variants from this one file:
!   plain     : no device I/O
!   -DDEVIO   : one `print` on a never-taken branch inside the target region
!   -DNOWAIT  : launch with `nowait` + taskwait instead of a synchronous target
program rpc_latency
   use omp_lib
   implicit none
   integer, parameter :: n = 27648, nrep = 100
   real(8), allocatable :: a(:)
   integer :: i, r, k, iters, kk
   real(8) :: t0, x
   allocate(a(n)); a = 1.0d0
   !$omp target enter data map(to:a)
   print '(a)', '   iters   wall/launch [us]'
   do k = 1, 6
      iters = 2**(k+5)*100
      do r = 1, 3 + nrep
         if (r == 4) t0 = omp_get_wtime()      ! 3 warm-up launches
#ifdef NOWAIT
         ! map(to:iters) is a workaround for PoC 07 (nowait + implicit firstprivate scalar = NYI)
         !$omp target teams distribute parallel do private(x, kk) map(to:iters) nowait
#else
         !$omp target teams distribute parallel do private(x, kk)
#endif
         do i = 1, n
            x = a(i)
            do kk = 1, iters
               x = x*0.999999d0 + 1.0d-7
            end do
#ifdef DEVIO
            if (x < -1.0d300) print *, 'never printed', i   ! never true, but links the RPC client
#endif
            a(i) = x
         end do
#ifdef NOWAIT
         !$omp taskwait
#endif
      end do
      print '(i8,f14.1)', iters, (omp_get_wtime() - t0)/nrep*1e6
   end do
   !$omp target update from(a)
   print '(a,es20.12)', 'checksum a(1) = ', a(1)
end program
