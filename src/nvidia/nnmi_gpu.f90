#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)
subroutine nnmi_gpu(x, eval, evec, no, bal, cutfreq, nnlist, &
                    cg_created, cg_graph, wrk)
   !
   !
   !  purpose : to get delta x vector at each iteration
   !-------------------------------------------------------------------
   !  **** input ****
   !  x    : vector matrix of variable array( tendency )
   !  eval : eigenval array of coeffienent matrix
   !  evec : eigenvector matrix array of coefficient matrix
   !  nn   : dimension of array
   !  c    : work array
   !  d    : work array
   !  cutfreq : cut of frequency
   !  **** output ****
   !  x    : vector matrix of variable array (correction)
   !  bal  : convergence indicator
   !
   !  modify to f90 by C-H Lee and sort by River Chen in 2005
   !
   !--------------------------------------------------------------------
   use rank
   use const, only: nnmivm
   use param, only: jtrun, jtmax
   use index, only: mlistnum, mlist
   use openacc
   use cudafor
   use cublas

   implicit none

   real, dimension(no*2, 2*jtmax*nnmivm), intent(inout):: x
   real, dimension(no*no, 2*jtmax*nnmivm), intent(in):: evec
   real, dimension(no, 2*jtmax*nnmivm), intent(in):: eval
   real, intent(inout):: bal(jtrun, nnmivm), cutfreq
   integer, intent(in):: no
   logical, intent(inout):: cg_created
   type(acc_graph_t), intent(inout):: cg_graph

   real, dimension(no, 2, 2*jtmax*nnmivm) :: wrk
   integer nnlist(2*jtmax*nnmivm)
   real bal_t, w1_t, w2_t, e_t, swap

   integer i, j, L, m, mf, nbig, nn, k
   integer ind
   integer async_id, ierr
   integer(kind=cuda_stream_kind) stream, subwrk_stream(2*jtmax*nnmivm)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle
   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   !
   ! set cut of frequency ( now , freq=1.0 means 24 hours )
   !
   !      data freq/1.0/
   if (.not. cg_created) then
      CUDACHECK(cudaEventCreate(spread_event))
      CUDACHECK(cudaEventCreate(pack_event))
      handle = cublasGetHandle()
      do L = 1, nnmivm
         do m = 1, mlistnum
            do k = 1, 2
               ind = k + (m - 1)*2 + (L - 1)*jtmax*2
               subwrk_stream(ind) = acc_get_cuda_stream(ind + 1)
            end do
         end do
      end do

      call accx_async_begin_capture(async_id)

      CUDACHECK(cudaEventRecord(spread_event, stream))
      do L = 1, nnmivm
         do m = 1, mlistnum
            do k = 1, 2
               ind = k + (m - 1)*2 + (L - 1)*jtmax*2
               CUDACHECK(cudaStreamWaitEvent(subwrk_stream(ind), spread_event, 0))
               ierr = cublasSetStream(handle, subwrk_stream(ind))
               !
               !   matrix multiplication to get alpha, (5.54)
               !
               nn = nnlist(ind)
               !  Inverse matrix array of eigenvector matrix is equal to
               !  the transpose of the eigenvector matrix.
               !$acc host_data use_device(evec, x, wrk)
               call dgemm('T', 'N', nn, 2, nn, &
                          1.0, evec(1, ind), nn, x(1, ind), no, &
                          0.0, wrk(1, 1, ind), no)
               !$acc end host_data
               CUDACHECK(cudaEventRecord(pack_event, subwrk_stream(ind)))
               CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
            end do
         end do
      end do
      !
      !   compute BAL, by (5.64)
      !
      !$acc parallel loop collapse(2) private(bal_t, mf) &
      !$acc& async(async_id)
      do L = 1, nnmivm
         do m = 1, mlistnum
            mf = mlist(m)
            bal_t = 0.
            !$acc loop vector collapse(2) reduction(+:bal_t) &
            !$acc& private(ind, nn, e_t, w1_t, w2_t, swap)
            do k = 1, 2
               do i = 1, no
                  ind = k + (m - 1)*2 + (L - 1)*jtmax*2
                  nn = nnlist(ind)
                  if (i .le. nn) then
                     e_t = eval(i, ind)
                     ! if (abs(e_t) .ge. cutfreq) then
                     if ((abs(e_t) .ge. cutfreq) &
                         .and. (nn .gt. 2)) then
                        !        if(abs(eval(i)).le.cutfreq)then
                        w1_t = wrk(i, 1, ind)
                        w2_t = wrk(i, 2, ind)
                        bal_t = bal_t + w1_t**2 + w2_t**2
                        swap = -w2_t/e_t
                        w2_t = w1_t/e_t
                        w1_t = swap
                     else
                        w1_t = 0.
                        w2_t = 0.
                     end if
                     wrk(i, 1, ind) = w1_t
                     wrk(i, 2, ind) = w2_t
                  end if
               end do
            end do
            bal(mf, L) = bal_t
         end do
      end do

      CUDACHECK(cudaEventRecord(spread_event, stream))
      do L = 1, nnmivm
         do m = 1, mlistnum
            do k = 1, 2
               ind = k + (m - 1)*2 + (L - 1)*jtmax*2
               CUDACHECK(cudaStreamWaitEvent(subwrk_stream(ind), spread_event, 0))
               nn = nnlist(ind)
               !
               !    transform alpha to get x, (5.55)
               !
               ierr = cublasSetStream(handle, subwrk_stream(ind))
               !$acc host_data use_device(evec, x, wrk)
               call dgemm('N', 'N', nn, 2, nn, &
                          1.0, evec(1, ind), nn, wrk(1, 1, ind), no, &
                          0., x(1, ind), no)
               !$acc end host_data
               CUDACHECK(cudaEventRecord(pack_event, subwrk_stream(ind)))
               CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
            end do
         end do
      end do
      call accx_async_end_capture(async_id, cg_graph)
      cg_created = .true.
      CUDACHECK(cudaEventDestroy(spread_event))
      CUDACHECK(cudaEventDestroy(pack_event))
   end if
   call accx_graph_launch(cg_graph, async_id)

   return
end subroutine nnmi_gpu
