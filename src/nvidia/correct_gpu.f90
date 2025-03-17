#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine correct_gpu(vorten, divten, phiten, tmcor, pmcor, vornow, &
                       divnow, temnow, plnow, jtrun, jtmax, lev)
   !
   !  purpose :  add correction to variables for NNMI
   !-----------------------------------------------------------------------
   !  **** input ****
   !  vorten  :  vorticity correction array in spectrum domaon
   !  divten  :  divergence correction array in spectrum domain
   !  phiten  :  geopotantial correction array in spectrum domain
   !  pmcor   :  array for compute surface pressure correction
   !  tmcor   :  array for compute temperature correction
   !  lev     :  total vertical levels
   !  **** output ****
   !  vornow  : vorticity spectrum coefficient
   !  divnow  : divergence spectrum coefficient
   !  temnow  : temperature spectrum coefficient
   !  plnow   : surface pressure spectrum coefficient
   !
   ! modify to f90 in 2015 by C-H Lee
   ! modify 2015 by River Chen
   !-------------------------------------------------------------------------
   !
   use const, only: RTYPE
   use index, only: levp, mlistnum, mlist, Llist, &
                    nsizex, nccl_row_comm
   use openacc
   use cudafor
   use cublas
   implicit none
   !
   integer, intent(in) :: jtrun, jtmax, lev
   real(kind=RTYPE), intent(in):: vorten(levp, 2, jtrun, jtmax), &
                                  divten(levp, 2, jtrun, jtmax), &
                                  phiten(levp, 2, jtrun, jtmax)

   real(kind=RTYPE), intent(in):: tmcor(lev, lev), pmcor(lev)
   real(kind=RTYPE), intent(inout):: vornow(levp, 2, jtrun, jtmax), &
                                     divnow(levp, 2, jtrun, jtmax), &
                                     temnow(levp, 2, jtrun, jtmax), &
                                     plnow(jtrun, jtmax, 2)

   real(kind=RTYPE):: temnow1(lev, 2, jtrun, jtmax), &
                      phiten1(lev, 2, jtrun, jtmax), loc_sum
   integer m, mf, n, k, i, kk, llistnum_fj

   integer async_id, ierr
   logical cg_created
   integer(kind=cuda_stream_kind) stream, matmul_stream(jtmax)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle
   type(acc_graph_t) cg_graph
   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(phiten1, temnow1) async(async_id)
   !2dMPI >
   call mpe2d_unify_lev_gpu(phiten, phiten1, lev, levp, jtrun, jtmax, &
                            mlistnum, nsizex, nccl_row_comm)
   call mpe2d_unify_lev_gpu(temnow, temnow1, lev, levp, jtrun, jtmax, &
                            mlistnum, nsizex, nccl_row_comm)
   !2dMPI <

   !
   !  compute surface pressure correction from geopotential correction and
   !  add to surface spectrum array
   !
   !$acc parallel loop collapse(3) private(mf, loc_sum) async(async_id)
   do m = 1, mlistnum
      do n = 1, jtrun
         do i = 1, 2
            mf = mlist(m)
            if (n .ge. mf) then
               loc_sum = 0.
               !$acc loop vector reduction(+:loc_sum)
               do k = 1, lev
                  loc_sum = loc_sum + pmcor(k)*phiten1(k, i, n, m)
               end do
               plnow(n, m, i) = plnow(n, m, i) + loc_sum
            end if
         end do
      end do
   end do
   !
   !  compute temperature correction from geopotential correction and
   !  add to temperature spectrum array
   !
   CUDACHECK(cudaEventCreate(spread_event))
   CUDACHECK(cudaEventCreate(pack_event))
   handle = cublasGetHandle()
   do m = 1, mlistnum
      matmul_stream(m) = acc_get_cuda_stream(m + 1)
   end do

   CUDACHECK(cudaEventRecord(spread_event, stream))
   do m = 1, mlistnum
      CUDACHECK(cudaStreamWaitEvent(matmul_stream(m), spread_event, 0))
      mf = mlist(m)
      llistnum_fj = jtrun - mf + 1
      n = llistnum_fj*2
      ierr = cublasSetStream(handle, matmul_stream(m))
      !$acc host_data use_device(tmcor, phiten1, temnow1)
#ifdef SP
      call sgemm('N', 'N', lev, n, lev, &
                 1._4, tmcor, lev, phiten1(1, 1, mf, m), lev, &
                 1._4, temnow1(1, 1, mf, m), lev)
#else
      call dgemm('N', 'N', lev, n, lev, &
                 1., tmcor, lev, phiten1(1, 1, mf, m), lev, &
                 1., temnow1(1, 1, mf, m), lev)
#endif
      !$acc end host_data
      CUDACHECK(cudaEventRecord(pack_event, matmul_stream(m)))
      CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
   end do

   CUDACHECK(cudaEventDestroy(spread_event))
   CUDACHECK(cudaEventDestroy(pack_event))
   !

   !$acc parallel loop collapse(4) private(mf, kk) async(async_id)
   do m = 1, mlistnum
      do n = 1, jtrun
         do i = 1, 2
            do k = 1, levp
               mf = mlist(m)
               if (n .ge. mf) then
                  KK = Llist(k)
                  !
                  !  add vorticity and divergence correction to vorticity and divergence
                  !  spectrum array
                  !
                  vornow(k, i, n, m) = vornow(k, i, n, m) + vorten(k, i, n, m)
                  divnow(k, i, n, m) = divnow(k, i, n, m) + divten(k, i, n, m)

                  !2dMPI >
                  temnow(k, i, n, m) = temnow1(KK, i, n, m)
                  !2dMPI <
               end if
            end do
         end do
      end do
   end do
   !$acc exit data delete(phiten1, temnow1) async(async_id)

   return
end subroutine correct_gpu
