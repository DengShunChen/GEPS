program lb7
   use omp_lib
   implicit none
   integer, parameter :: n1 = 769, n2 = 72, n3 = 64, nrep = 2000
   real(8), allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:), a4(:,:,:), a5(:,:,:), a6(:,:,:), a7(:,:,:), a8(:,:,:)
   real(8), allocatable :: b1(:,:), b2(:,:), b3(:,:), b4(:,:), b5(:,:), b6(:,:)
   integer, allocatable :: i1(:,:), i2(:,:), i3(:,:)
   integer :: ot, r
   real(8) :: t0, t1, sc
   allocate(a1(n1,n2,n3),a2(n1,n2,n3),a3(n1,n2,n3),a4(n1,n2,n3),a5(n1,n2,n3),a6(n1,n2,n3),a7(n1,n2,n3),a8(n1,n2,n3))
   allocate(b1(n2,n3),b2(n2,n3),b3(n2,n3),b4(n2,n3),b5(n2,n3),b6(n2,n3),i1(n2,n3),i2(n2,n3),i3(n2,n3))
   a1=1; a2=2; a3=3; a4=4; a5=5; a6=6; a7=7; a8=8; b1=1; b2=1; b3=1; b4=1; b5=1; b6=1; i1=1; i2=2; i3=3
   sc = 0.5d0
   !$omp target enter data map(to:a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)
   !$omp target data use_device_addr(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)
   do r = 1, 10
      !$omp target teams distribute parallel do firstprivate(sc) has_device_addr(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + b3(1,ot) + b4(1,ot) + b5(1,ot) + b6(1,ot) + i1(1,ot) + i2(1,ot) + i3(1,ot) &
            + a1(1,1,ot) + a2(1,1,ot) + a3(1,1,ot) + a4(1,1,ot) + a5(1,1,ot) + a6(1,1,ot) + a7(1,1,ot) + a8(1,1,ot)*sc
      end do
   end do
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do firstprivate(sc) has_device_addr(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + b3(1,ot) + b4(1,ot) + b5(1,ot) + b6(1,ot) + i1(1,ot) + i2(1,ot) + i3(1,ot) &
            + a1(1,1,ot) + a2(1,1,ot) + a3(1,1,ot) + a4(1,1,ot) + a5(1,1,ot) + a6(1,1,ot) + a7(1,1,ot) + a8(1,1,ot)*sc
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'has_device_addr, 17 arrays:    ', (t1-t0)/nrep*1e6, ' us/launch'
   !$omp end target data
   !$omp target update from(b1)
   print *, 'b1(1,1) =', b1(1,1), ' expected', 1 + (nrep+10)*(5+6+36.0d0/2+ 0.5d0*8)
end program lb7
