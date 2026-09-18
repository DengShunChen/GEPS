#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine zx_gpu(evec, vorten, divten, phiten, jtrun, jtmax, lev, &
                  wrk, cg_created, cg_graph)
   !
   !  purpose : do vertical transform and inverse transform
   !------------------------------------------------------------------------
   !  **** input ****
   !  evec   : matrix array for transform
   !  vorten : vorticity tendency(correction) array of spectrum coefficients
   !  divten : divergence tendency(correction) array of spectrum coefficients
   !  phiten : geopotential tendency(correction)array of spectrum coefficients
   !  lev    : total vertical levels
   !  **** output ****
   !  vorten : vorticity tendency(correction) array of spectrum coefficients
   !  divten : divergence tendency(correction) array of spectrum coefficients
   !  phiten : geopotential tendency(correction)array of spectrum coefficients
   !  **** present on device ****
   !  evec, vorten, divten, phiten
   !  mlist, Llist
   !  jtrun, jtmax, lev, levp
   !---------------------------------------------------------------------------
   use index
   use const, only: RTYPE
   use openacc
   use cudafor
   use cublas

   implicit none

   integer, intent(in):: lev, jtrun, jtmax

   real(kind=RTYPE), intent(in):: evec(lev, lev)
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax), intent(inout):: &
      vorten, phiten, divten
   real(kind=RTYPE) vor(lev, 2, jtrun, jtmax), &
      div(lev, 2, jtrun, jtmax), &
      phe(lev, 2, jtrun, jtmax)

   !2dMPI
   REAL(kind=RTYPE), dimension(lev, 2, 3, jtrun, jtmax):: wrk
   REAL(kind=RTYPE), dimension(:), allocatable, device :: vars
   integer m, mf, n, j, l, k, KK, KL, llistnum_fj, idx
   REAL(kind=RTYPE) ONE, ZERO

   integer async_id, ierr
   logical cg_created
   integer(kind=cuda_stream_kind) stream, matmul_stream(jtmax)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle
   type(acc_graph_t) cg_graph
   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   ONE = 1.
   ZERO = 0.
   CALL mpe2d_unify_spec_lev_zx_gpu(wrk, vorten, divten, phiten, &
                                    lev, levp, jtrun, jtmax, mlistnum, &
                                    nsizex, nccl_row_comm)

   if (.not. cg_created) then
      CUDACHECK(cudaEventCreate(spread_event))
      CUDACHECK(cudaEventCreate(pack_event))
      handle = cublasGetHandle()
      do m = 1, mlistnum
         matmul_stream(m) = acc_get_cuda_stream(m + 1)
      end do

#ifdef USE_HIP
      ! Allocate `vars` before capture starts (not stream-ordered inside the
      ! graph): the `!$omp target ... has_device_addr(vars)` kernel below
      ! bakes in whatever address `vars` holds at capture time, and amdflang's
      ! libomptarget does not re-virtualize that address against a graph
      ! memory node the way CUDA Fortran's `device`-attributed cudaMallocAsync
      ! does on NVIDIA. Allocating first gives `vars` a real, stable address
      ! before anything is recorded, so the baked-in address stays valid for
      ! every future graph replay (root-caused 2026-09-08: without this,
      ! replay hit a GPU memory access fault reading a stale/placeholder
      ! address). Never freed -- it must outlive every replay of this graph.
      CUDACHECK(cudaMallocAsync(vars, lev*6*jtrun*jtmax, stream))
      CUDACHECK(cudaStreamSynchronize(stream))
#endif
      call accx_async_begin_capture(async_id)

#ifndef USE_HIP
      CUDACHECK(cudaMallocAsync(vars, lev*6*jtrun*jtmax, stream))
#endif

      CUDACHECK(cudaEventRecord(spread_event, stream))

      do m = 1, mlistnum
         CUDACHECK(cudaStreamWaitEvent(matmul_stream(m), spread_event, 0))
         mf = mlist(m)
         llistnum_fj = jtrun - mf + 1
         n = llistnum_fj*6
         idx = 1 + (mf - 1)*lev*6 + (m - 1)*lev*6*jtrun
         ierr = cublasSetStream(handle, matmul_stream(m))
         !$acc host_data use_device(evec, wrk)
#ifdef SP
         call sgemm('N', 'N', lev, n, lev, &
                    1._4, evec, lev, wrk(1, 1, 1, mf, m), lev, &
                    0._4, vars(idx), lev)
#else
         call dgemm('N', 'N', lev, n, lev, &
                    1., evec, lev, wrk(1, 1, 1, mf, m), lev, &
                    0., vars(idx), lev)
#endif
         !$acc end host_data
         CUDACHECK(cudaEventRecord(pack_event, matmul_stream(m)))
         CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
      end do

      !
#ifdef USE_HIP
      ! `vars` is a raw device pointer from cudaMallocAsync (CUDA Fortran
      ! `device` attribute stripped to a plain `pointer` by acc2omp -- see
      ! cmake/acc2omp.py). NVIDIA's OpenACC compiler natively recognizes a
      ! `device`-attributed array as already device-resident and emits no
      ! mapping for it. amdflang/libomptarget has no such information for
      ! the translated `pointer`, so an ordinary implicit-map
      ! `!$omp target teams distribute` here tries to treat `vars`'s (device)
      ! address as a host address needing a fresh host->device copy, which
      ! fails with "hsa_amd_memory_lock: HSA_STATUS_ERROR" (root-caused
      ! 2026-09-08). `has_device_addr(vars)` tells the compiler `vars`'s
      ! address is already a valid device address, so no copy is attempted.
      ! (`is_device_ptr` was tried first but amdflang requires that clause's
      ! argument to be literally type(c_ptr); `vars` is a Fortran pointer.)
      !$omp target teams distribute parallel do collapse(4) private(mf, kk, idx) has_device_addr(vars)
#else
      !$acc parallel loop collapse(4) private(mf, kk, idx) async(async_id)
#endif
      do m = 1, mlistnum
         do l = 1, jtrun
            do n = 1, 2
               do k = 1, levp
                  mf = mlist(m)
                  if (l .ge. mf) then
                     kk = Llist(k)
                     idx = kk + (n - 1)*lev + (l - 1)*lev*6 + (m - 1)*lev*6*jtrun
                     vorten(k, n, l, m) = vars(idx)
                     divten(k, n, l, m) = vars(idx + lev*2)
                     phiten(k, n, l, m) = vars(idx + lev*4)
                  end if
               end do
            end do
         end do
      end do
#ifndef USE_HIP
      CUDACHECK(cudaFreeAsync(vars, stream))
#endif

      call accx_async_end_capture(async_id, cg_graph)
      cg_created = .true.
      CUDACHECK(cudaEventDestroy(spread_event))
      CUDACHECK(cudaEventDestroy(pack_event))
   end if
   call accx_graph_launch(cg_graph, async_id)
   !
   return
end subroutine zx_gpu
