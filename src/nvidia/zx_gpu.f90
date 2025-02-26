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

      call accx_async_begin_capture(async_id)

      CUDACHECK(cudaMallocAsync(vars, lev*6*jtrun*jtmax, stream))

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
      !$acc parallel loop collapse(4) private(mf, kk, idx) async(async_id)
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
      CUDACHECK(cudaFreeAsync(vars, stream))

      call accx_async_end_capture(async_id, cg_graph)
      cg_created = .true.
      CUDACHECK(cudaEventDestroy(spread_event))
      CUDACHECK(cudaEventDestroy(pack_event))
   end if
   call accx_graph_launch(cg_graph, async_id)
   !
   return
end subroutine zx_gpu
