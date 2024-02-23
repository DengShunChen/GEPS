!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fft_create_plan(plan, auto)
    use cufft
    implicit none
    integer(4) :: plan, istat, auto

    istat = cufftCreate(plan)
    if (istat .ne. CUFFT_SUCCESS) then
        print *, 'cufftCreate', istat
    end if
    if (auto .eq. 0) then
        istat = cufftSetAutoAllocation(plan, 0)
        if (istat .ne. CUFFT_SUCCESS) then
            print *, 'cufftSetAutoAllocation', istat
        end if
    end if
end subroutine fft_create_plan

subroutine fft_make_plan(inc, jump, n, m, isign, plan, work_size)
    use cufft
    implicit none

    integer :: inc, jump, n, m, isign, nn

    integer :: isign
    integer :: istat
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
        ffttype = CUFFT_Z2D
        batch = m
    else
        rank = 1
        inembed = n
        istride = inc
        idist = jump
        onembed = n
        ostride = inc
        odist = nn
        ffttype = CUFFT_D2Z
        batch = m
    end if
    istat = cufftMakePlanMany(plan, rank, n, inembed, istride, idist, onembed, ostride, odist, ffttype, batch, work_size)
    if (istat .ne. CUFFT_SUCCESS) then
        print *, 'cufftMakePlanMany', istat
    end if
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
    integer(4) :: istat
    integer :: i

    call fft_plan_size(num_plan)
    allocate (plan_list(num_plan))
    call get_plan_list(plan_list)

    do i = 1, num_plan
        istat = cufftDestroy(plan_list(i))
        if (istat .ne. CUFFT_SUCCESS) then
            print *, 'cufftDestroy', istat
        end if
    end do

end subroutine fft_destroy_plan

subroutine fft_set_workspace(plan, workspace)
    use cufft
    implicit none
    integer(4) :: plan, istat
    integer :: workspace(*)

    !$acc host_data use_device(workspace)
    istat = cufftSetWorkArea(plan, workspace)
    !$acc end host_data
    if (istat .ne. CUFFT_SUCCESS) then
        print *, 'cufftSetWorkArea', istat
    end if
end subroutine fft_set_workspace

subroutine fft_exec_async(a, inc, jump, n, m, isign, plan, pa, async_id)
    use cudafor
    use cufft
    use openacc
    implicit none
    real(8) :: scale
    integer :: inc, jump, n, m, isign, i, j
    real(8), dimension(jump, m) :: pa
    real(8), dimension(jump, m) :: a
    integer(4) :: plan, istat
    integer :: async_id
    integer(kind=cuda_stream_kind) :: stream

    !$acc data copy(a) create(pa) async(async_id)
    if (mod(jump, 2) .eq. 0) then

        stream = acc_get_cuda_stream(async_id)
        istat = cufftSetStream(plan, stream)

        if (isign .eq. 1) then

            !$acc parallel loop collapse(2) async(async_id)
            do j = 1, m
            do i = 1, n + 2
                pa(i, j) = a(i, j)
            end do
            end do

            !$acc host_data use_device(pa, a)
            istat = cufftExecZ2D(plan, pa, a)
            !$acc end host_data
            if (istat .ne. CUFFT_SUCCESS) then
                print *, 'cufftExecZ2D', istat
            end if
        else
            scale = 1.0/dfloat(n)

            !$acc host_data use_device(a, pa)
            istat = cufftExecD2Z(plan, a, pa)
            !$acc end host_data
            if (istat .ne. CUFFT_SUCCESS) then
                print *, 'cufftExecD2Z', istat
            end if

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
    istat = cudaStreamSynchronize(stream)
    if (istat .ne. cudaSuccess) then
        print *, 'cudaStreamSynchronize', istat
    end if
    istat = cudaGetLastError()
    if (istat .ne. cudaSuccess) then
        print *, 'cudaGetLastError', istat
    end if
#endif

    return
end subroutine fft_exec_async

subroutine rfftmlt_gpu(a, work, trigs, ifax, inc, jump, n, m, isign)
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
    real(8), dimension(jump, m) :: a
    real(8), dimension(jump, m) :: work
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
    real(8), dimension(jump, m, *) :: cc, gwk1
    real, dimension(4096, *) :: trigsj
    integer, dimension(19, *) :: ifaxj
    integer, dimension(*) :: jlist1, nxdef
    integer(4), dimension(jlistnum) :: plan_list
    integer(kind=int_ptr_kind()) :: work_size

    do jj = 1, jlistnum
        j = jlist1(jj)
        nxj = nxdef(j)
        call find_fft_plan(1, jump, nxj, m, isign, plan_list(jj))
        if (plan_list(jj) .eq. -1) then
            call fft_create_plan(plan_list(jj), 1)
            call fft_make_plan(1, jump, nxj, m, isign, plan_list(jj), work_size)
            call cache_fft_plan(1, jump, nxj, m, isign, plan_list(jj))
        end if
    end do

    ! To avoid one plan being accessed by two cufft exec, we set the stream id equal to plan id
    do jj = 1, jlistnum
        j = jlist1(jj)
        nxj = nxdef(j)
        call fft_exec_async(cc(1, 1, jj), 1, jump, nxj, m, isign, plan_list(jj), gwk1(1, 1, jj), plan_list(jj))
    end do

    do jj = 1, jlistnum
        !$acc wait(plan_list(jj))
    end do

end subroutine rfftmlt_loop


