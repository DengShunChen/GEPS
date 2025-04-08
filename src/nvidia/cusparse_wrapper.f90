subroutine cusparseDgtsvInterleavedBatch_async &
   (handle, algo, m, dl, d, du, x, batchCount, async_id)
   use cudafor
   use cusparse
   use openacc
   use machine, only: kind_phys
   implicit none
   real(kind=kind_phys), dimension(*) :: dl, d, du, x
   integer :: m, batchCount, algo
   integer(kind=cuda_stream_kind) :: stream
   integer :: istat, async_id
   integer(8), value :: buffer_size
   character(1), allocatable :: pbuffer(:)
   type(cusparseHandle) :: handle

   stream = acc_get_cuda_stream(async_id)
   istat = cusparseSetStream(handle, stream)

   !$acc host_data use_device(dl,d,du,x)
   !istat = cusparseDgtsv2StridedBatch_bufferSizeExt &
   !        (handle,m,dl,d,du,x,batchCount,batchStride,buffer_size)
   istat = cusparseDgtsvInterleavedBatch_bufferSizeExt &
           (handle, algo, m, dl, d, du, x, batchCount, buffer_size)
   !$acc wait(async_id)
   !$acc end host_data

   allocate (pbuffer(buffer_size))
   !$acc enter data create(pbuffer) async(async_id)
   !$acc host_data use_device(dl,d,du,x,pbuffer)
   !istat = cusparseDgtsv2StridedBatch &
   !        (handle,m,dl,d,du,x,batchCount,batchStride,pbuffer)
   istat = cusparseDgtsvInterleavedBatch &
           (handle, algo, m, dl, d, du, x, batchCount, pbuffer)
   !$acc end host_data
   !$acc exit data delete(pbuffer) async(async_id)
   !$acc wait(async_id)
   deallocate (pbuffer)

end subroutine cusparseDgtsvInterleavedBatch_async
