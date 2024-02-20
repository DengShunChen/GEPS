!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine device_init(rank, size)
    use openacc
    implicit none
    integer :: rank, size
    integer :: num_device, device_id

    num_device = acc_get_num_devices(acc_device_nvidia)
    device_id = mod(rank, num_device)
    call acc_set_device_num(device_id, acc_device_nvidia)
end subroutine device_init

