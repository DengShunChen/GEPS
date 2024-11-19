!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)
#define CUFFTCHECK(ierr) call cufft_check_helper(ierr, __FILE__, __LINE__)

subroutine fft_create_plan(plan, auto)
   use cufft
   implicit none
   integer(4) :: plan, auto

   CUFFTCHECK(cufftCreate(plan))
   if (auto .eq. 0) then
      CUFFTCHECK(cufftSetAutoAllocation(plan, 0))
   end if
end subroutine fft_create_plan

subroutine fft_make_plan(inc, jump, n, m, isign, plan, work_size)
   use cufft
   implicit none

   integer :: inc, jump, n, m, isign, nn

   integer :: isign
   integer(4) :: plan
   integer(4) :: rank
   integer :: inembed, onembed
   integer(4) :: istride, idist, ostride, odist
   integer(4) :: ffttype, batch
   integer(kind=int_ptr_kind()) :: work_size

   nn = jump/2

   if (isign .eq. 1) then
      rank = 1
      inembed = n
      istride = inc
      idist = nn
      onembed = n
      ostride = inc
      odist = jump
#ifdef SP
      ffttype = CUFFT_C2R
#else
      ffttype = CUFFT_Z2D
#endif
      batch = m
   else
      rank = 1
      inembed = n
      istride = inc
      idist = jump
      onembed = n
      ostride = inc
      odist = nn
#ifdef SP
      ffttype = CUFFT_R2C
#else
      ffttype = CUFFT_D2Z
#endif
      batch = m
   end if
   CUFFTCHECK(cufftMakePlanMany(plan, rank, n, inembed, istride, idist, onembed, ostride, odist, ffttype, batch, work_size))
end subroutine fft_make_plan

subroutine fft_destroy_plan()
   use cufft
   implicit none

   interface c_interface
      subroutine fft_plan_size(size) bind(C, name="fft_plan_size")
         use iso_c_binding
         implicit none
         integer(c_int), intent(out) :: size
      end subroutine fft_plan_size
      subroutine get_plan_list(plan_list) bind(C, name="get_plan_list")
         use iso_c_binding
         implicit none
         integer(c_int), intent(out), dimension(*) :: plan_list
      end subroutine get_plan_list
   end interface c_interface

   integer :: num_plan
   integer(4), allocatable :: plan_list(:)
   integer :: i

   call fft_plan_size(num_plan)
   allocate (plan_list(num_plan))
   call get_plan_list(plan_list)

   do i = 1, num_plan
      CUFFTCHECK(cufftDestroy(plan_list(i)))
   end do

end subroutine fft_destroy_plan

subroutine fft_set_workspace(plan, workspace)
   use cufft
   implicit none
   integer(4) :: plan
   integer :: workspace(*)

   !$acc host_data use_device(workspace)
   CUFFTCHECK(cufftSetWorkArea(plan, workspace))
   !$acc end host_data
end subroutine fft_set_workspace

subroutine fft_exec_async(a, inc, jump, n, m, isign, plan, pa, async_id)
  ! ------------------------------------------------------------
  ! The input and output varialbes depend on value of isign
  ! isign = 1:
  !   - input: pa
  !   - output: a
  ! ising = -1:
  !   - input: a
  !   - output: a
  ! ------------------------------------------------------------
   use const, only: RTYPE
   use cudafor
   use cufft
   use openacc
   implicit none
   real(8) :: scale
   integer :: inc, jump, n, m, isign, i, j
   real(kind=RTYPE), dimension(jump, m) :: pa
   real(kind=RTYPE), dimension(jump, m) :: a
   integer(4) :: plan
   integer :: async_id
   integer(kind=cuda_stream_kind) :: stream

   !$acc data copy(a) create(pa) async(async_id)
   if (mod(jump, 2) .eq. 0) then

      stream = acc_get_cuda_stream(async_id)
      CUFFTCHECK(cufftSetStream(plan, stream))

      if (isign .eq. 1) then

         !$acc update device(pa) async(async_id)
         ! !$acc parallel loop collapse(2) async(async_id)
         ! do j = 1, m
         ! do i = 1, n + 2
         !    pa(i, j) = a(i, j)
         ! end do
         ! end do

         !$acc host_data use_device(pa, a)
#ifdef SP
         CUFFTCHECK(cufftExecC2R(plan, pa, a))
#else
         CUFFTCHECK(cufftExecZ2D(plan, pa, a))
#endif
         !$acc end host_data
      else
         scale = 1.0/float(n)

         !$acc host_data use_device(a, pa)
#ifdef SP
         CUFFTCHECK(cufftExecR2C(plan, a, pa))
#else
         CUFFTCHECK(cufftExecD2Z(plan, a, pa))
#endif
         !$acc end host_data

         !$acc parallel loop collapse(2) async(async_id)
         do j = 1, m
         do i = 1, n + 2
            a(i, j) = pa(i, j)*scale
         end do
         end do
      end if
   else
      print *, 'fft jump is odd, CWB obsoleted, jump=', jump
   end if
   !$acc end data

#if DEBUG
   CUDACHECK(cudaStreamSynchronize(stream))
#endif

   return
end subroutine fft_exec_async

subroutine rfftmlt_gpu(a, work, trigs, ifax, inc, jump, n, m, isign)
   use const, only: RTYPE
   use iso_c_binding
   implicit none

   interface c_interface
      subroutine find_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="find_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign
         integer(c_int), intent(out) :: plan
      end subroutine find_fft_plan
      subroutine cache_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="cache_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign, plan
      end subroutine cache_fft_plan
   end interface c_interface

   integer(4) :: plan
   real(8) :: trigs(*)
   integer :: ifax(*), inc, jump, n, m, isign
   real(kind=RTYPE), dimension(jump, m) :: a
   real(kind=RTYPE), dimension(jump, m) :: work
   integer(kind=int_ptr_kind()) :: work_size
   integer(4) :: async_id

   call find_fft_plan(inc, jump, n, m, isign, plan)

   if (plan .eq. -1) then
      call fft_create_plan(plan, 1)
      call fft_make_plan(inc, jump, n, m, isign, plan, work_size)
      call cache_fft_plan(inc, jump, n, m, isign, plan)
   end if

   async_id = plan

   call fft_exec_async(a, inc, jump, n, m, isign, plan, work, async_id)
   !$acc wait(async_id)
end subroutine rfftmlt_gpu

subroutine fftfax_gpu(n, ifax, trigs)
   implicit none
   real*8 :: trigs(*)
   integer :: ifax(*), n
end subroutine fftfax_gpu

subroutine rfftmlt_loop(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, jump, m, isign)
   use const, only: RTYPE
   use cudafor
   use cufft
   use openacc
   use iso_c_binding
   implicit none

   interface c_interface
      subroutine find_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="find_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign
         integer(c_int), intent(out) :: plan
      end subroutine find_fft_plan
      subroutine cache_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="cache_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign, plan
      end subroutine cache_fft_plan
   end interface c_interface

   integer :: jlistnum, jump, m, isign, jlistnum
   integer :: jj, j, nxj
   real(kind=RTYPE), dimension(jump, m, *) :: cc, gwk1 ! Present on device
   real, dimension(4096, *) :: trigsj
   integer, dimension(19, *) :: ifaxj
   integer, dimension(*) :: jlist1, nxdef
   integer(4), dimension(jlistnum) :: plan_list
   integer(kind=int_ptr_kind()) :: work_size
   integer :: async_id
   integer(kind=cuda_stream_kind) :: stream
   integer :: i, k
   real(8) :: scale

   async_id = 1

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      call find_fft_plan(1, jump, nxj, m, isign, plan_list(jj))
      if (plan_list(jj) .eq. -1) then
         call fft_create_plan(plan_list(jj), 1)
         call fft_make_plan(1, jump, nxj, m, isign, plan_list(jj), work_size)
         call cache_fft_plan(1, jump, nxj, m, isign, plan_list(jj))
      end if

      if (mod(jump, 2) .eq. 0) then
         stream = acc_get_cuda_stream(async_id)
         CUFFTCHECK(cufftSetStream(plan_list(jj), stream))

         if (isign .eq. 1) then

            !$acc parallel loop collapse(2) async(async_id)
            do k = 1, m
            do i = 1, nxj + 2
               gwk1(i, k, jj) = cc(i, k, jj)
            end do
            end do

            !$acc host_data use_device(cc, gwk1)
#ifdef SP
            CUFFTCHECK(cufftExecC2R(plan_list(jj), gwk1(1, 1, jj), cc(1, 1, jj)))
#else
            CUFFTCHECK(cufftExecZ2D(plan_list(jj), gwk1(1, 1, jj), cc(1, 1, jj)))
#endif
            !$acc end host_data
         else
            scale = 1.0/dfloat(nxj)

            !$acc host_data use_device(cc, gwk1)
#ifdef SP
            CUFFTCHECK(cufftExecR2C(plan_list(jj), cc(1, 1, jj), gwk1(1, 1, jj)))
#else
            CUFFTCHECK(cufftExecD2Z(plan_list(jj), cc(1, 1, jj), gwk1(1, 1, jj)))
#endif
            !$acc end host_data

            !$acc parallel loop collapse(2) async(async_id)
            do k = 1, m
            do i = 1, nxj + 2
               cc(i, k, jj) = gwk1(i, k, jj)*scale
            end do
            end do
         end if

      else
         print *, 'fft jump is odd, CWB obsoleted, jump=', jump
      end if
   end do

end subroutine rfftmlt_loop

subroutine rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, jump, m, isign)
   ! Present on device: cc, gwk1, jlist1, nxdef
   use cudafor
   use cufft
   use openacc
   use iso_c_binding
   implicit none

   interface c_interface
      subroutine find_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="find_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign
         integer(c_int), intent(out) :: plan
      end subroutine find_fft_plan
      subroutine cache_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="cache_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign, plan
      end subroutine cache_fft_plan
   end interface c_interface

   integer :: jlistnum, jump, m, isign, jlistnum
   integer :: jj, j, nxj
   real(8), dimension(jump, m, *) :: cc, gwk1 ! Present on device
   real, dimension(4096, *) :: trigsj
   integer, dimension(19, *) :: ifaxj
   integer, dimension(*) :: jlist1, nxdef
   integer(kind=int_ptr_kind()) :: work_size
   integer :: async_id
   integer(kind=cuda_stream_kind) :: stream, plan_stream(jlistnum)
   integer :: i, k
   integer(4) :: base_plan, plan_id, plan_cur
   type(cudaEvent) :: spread_event, pack_event

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   j = jlist1(1)
   nxj = nxdef(j)
   call find_fft_plan(1, jump, nxj, m, isign, base_plan)
   if (base_plan .eq. -1) then
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         call fft_create_plan(plan_id, 1)
         call fft_make_plan(1, jump, nxj, m, isign, plan_id, work_size)
         call find_fft_plan(1, jump, nxj, m, isign, plan_cur)
         if (plan_cur .eq. -1) then
            call cache_fft_plan(1, jump, nxj, m, isign, plan_id)
         end if
      end do
      base_plan = plan_id - jlistnum + 1
   end if

   do jj = 1, jlistnum
      plan_id = base_plan + jj - 1
      plan_stream(jj) = acc_get_cuda_stream(plan_id)
   end do

   CUDACHECK(cudaEventCreate(spread_event))
   CUDACHECK(cudaEventCreate(pack_event))

   if (mod(jump, 2) .eq. 0) then
      if (isign .eq. 1) then
         !$acc parallel loop gang collapse(2) private(j, nxj) async(async_id)
         do jj = 1, jlistnum
            do k = 1, m
               j = jlist1(jj)
               nxj = nxdef(j)
               !$acc loop vector
               do i = 1, nxj + 2
                  gwk1(i, k, jj) = cc(i, k, jj)
               end do
            end do
         end do
         CUDACHECK(cudaEventRecord(spread_event, stream))
         do jj = 1, jlistnum
            plan_id = base_plan + jj - 1
            CUDACHECK(cudaStreamWaitEvent(plan_stream(jj), spread_event, 0))
            CUFFTCHECK(cufftSetStream(plan_id, plan_stream(jj)))
            !$acc host_data use_device(cc, gwk1)
            CUFFTCHECK(cufftExecZ2D(plan_id, gwk1(1, 1, jj), cc(1, 1, jj)))
            !$acc end host_data
            CUDACHECK(cudaEventRecord(pack_event, plan_stream(jj)))
            CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
         end do
      else
         do jj = 1, jlistnum
            plan_id = base_plan + jj - 1
            CUFFTCHECK(cufftSetStream(plan_id, stream))
            !$acc host_data use_device(cc, gwk1)
            CUFFTCHECK(cufftExecD2Z(plan_id, cc(1, 1, jj), gwk1(1, 1, jj)))
            !$acc end host_data
         end do
         !$acc parallel loop gang collapse(2) private(j, nxj) async(async_id)
         do jj = 1, jlistnum
            do k = 1, m
               j = jlist1(jj)
               nxj = nxdef(j)
               !$acc loop vector
               do i = 1, nxj + 2
                  cc(i, k, jj) = gwk1(i, k, jj)/dfloat(nxj)
               end do
            end do
         end do
      end if
   else
      print *, 'fft jump is odd, CWB obsoleted, jump=', jump
   end if

end subroutine rfftmlt_loop_identical

subroutine rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, jump, m, isign, graph)
   ! Present on device: cc, gwk1, jlist1, nxdef
  ! ------------------------------------------------------------
  ! The input and output varialbes depend on value of isign
  ! isign = 1:
  !   - input: gwk1
  !   - output: cc
  ! ising = -1:
  !   - input: cc
  !   - output: cc
  ! ------------------------------------------------------------
   use const, only: RTYPE
   use cudafor
   use cufft
   use openacc
   use iso_c_binding
   implicit none

   interface c_interface
      subroutine find_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="find_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign
         integer(c_int), intent(out) :: plan
      end subroutine find_fft_plan
      subroutine cache_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="cache_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign, plan
      end subroutine cache_fft_plan
   end interface c_interface

   integer :: jlistnum, jump, m, isign
   integer :: jj, j, nxj
   real(kind=RTYPE), dimension(jump, m, *) :: cc, gwk1 ! Present on device
   real, dimension(4096, *) :: trigsj
   integer, dimension(19, *) :: ifaxj
   integer, dimension(*) :: jlist1, nxdef
   type(cudaGraph) :: graph
   integer(kind=int_ptr_kind()) :: work_size
   integer :: async_id
   integer(kind=cuda_stream_kind) :: stream, plan_stream(jlistnum)
   integer :: i, k
   integer(4) :: base_plan, plan_id, plan_cur
   type(cudaEvent) :: spread_event, pack_event

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   j = jlist1(1)
   nxj = nxdef(j)
   call find_fft_plan(1, jump, nxj, m, isign, base_plan)
   if (base_plan .eq. -1) then
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         call fft_create_plan(plan_id, 1)
         call fft_make_plan(1, jump, nxj, m, isign, plan_id, work_size)
         call find_fft_plan(1, jump, nxj, m, isign, plan_cur)
         if (plan_cur .eq. -1) then
            call cache_fft_plan(1, jump, nxj, m, isign, plan_id)
         end if
      end do
      base_plan = plan_id - jlistnum + 1
   end if

   do jj = 1, jlistnum
      plan_id = base_plan + jj - 1
      plan_stream(jj) = acc_get_cuda_stream(plan_id)
   end do

   CUDACHECK(cudaEventCreate(spread_event))
   CUDACHECK(cudaEventCreate(pack_event))

   CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))

   if (mod(jump, 2) .eq. 0) then
      if (isign .eq. 1) then
         ! !$acc parallel loop gang collapse(2) private(j, nxj) async(async_id)
         ! do jj = 1, jlistnum
         !    do k = 1, m
         !       j = jlist1(jj)
         !       nxj = nxdef(j)
         !       !$acc loop vector
         !       do i = 1, nxj + 2
         !          gwk1(i, k, jj) = cc(i, k, jj)
         !       end do
         !    end do
         ! end do
         CUDACHECK(cudaEventRecord(spread_event, stream))
         do jj = 1, jlistnum
            plan_id = base_plan + jj - 1
            CUDACHECK(cudaStreamWaitEvent(plan_stream(jj), spread_event, 0))
            CUFFTCHECK(cufftSetStream(plan_id, plan_stream(jj)))
            !$acc host_data use_device(cc, gwk1)
#ifdef SP
            CUFFTCHECK(cufftExecC2R(plan_id, gwk1(1, 1, jj), cc(1, 1, jj)))
#else
            CUFFTCHECK(cufftExecZ2D(plan_id, gwk1(1, 1, jj), cc(1, 1, jj)))
#endif
            !$acc end host_data
            CUDACHECK(cudaEventRecord(pack_event, plan_stream(jj)))
            CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
         end do
      else
         CUDACHECK(cudaEventRecord(spread_event, stream))
         do jj = 1, jlistnum
            plan_id = base_plan + jj - 1
            CUDACHECK(cudaStreamWaitEvent(plan_stream(jj), spread_event, 0))
            CUFFTCHECK(cufftSetStream(plan_id, plan_stream(jj)))
            !$acc host_data use_device(cc, gwk1)
#ifdef SP
            CUFFTCHECK(cufftExecR2C(plan_id, cc(1, 1, jj), gwk1(1, 1, jj)))
#else
            CUFFTCHECK(cufftExecD2Z(plan_id, cc(1, 1, jj), gwk1(1, 1, jj)))
#endif
            !$acc end host_data
            CUDACHECK(cudaEventRecord(pack_event, plan_stream(jj)))
            CUDACHECK(cudaStreamWaitEvent(stream, pack_event, 0))
         end do
         !$acc parallel loop gang collapse(2) private(j, nxj) async(async_id)
         do jj = 1, jlistnum
            do k = 1, m
               j = jlist1(jj)
               nxj = nxdef(j)
               !$acc loop vector
               do i = 1, nxj + 2
                  cc(i, k, jj) = gwk1(i, k, jj)/float(nxj)
               end do
            end do
         end do
      end if
   else
      print *, 'fft jump is odd, CWB obsoleted, jump=', jump
   end if

   CUDACHECK(cudaStreamEndCapture(stream, graph))
   CUDACHECK(cudaEventDestroy(spread_event))
   CUDACHECK(cudaEventDestroy(pack_event))

end subroutine rfftmlt_loop_identical_cuda_graph
