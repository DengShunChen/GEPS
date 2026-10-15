! launch-overhead microbenchmark: mirrors a translated GEPS kernel
! (collapse(2) over 64x72, ~20 mapped arrays, firstprivate scalars).
program lb
   use omp_lib
   implicit none
   integer, parameter :: n1 = 769, n2 = 72, n3 = 64, nrep = 2000
   real(8), allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:), a4(:,:,:), a5(:,:,:), a6(:,:,:), a7(:,:,:), a8(:,:,:)
   real(8), allocatable :: b1(:,:), b2(:,:), b3(:,:), b4(:,:), b5(:,:), b6(:,:)
   integer, allocatable :: i1(:,:), i2(:,:), i3(:,:)
   integer :: ot, k, i, r, nloc, ob
   real(8) :: t0, t1, sc
   allocate(a1(n1,n2,n3),a2(n1,n2,n3),a3(n1,n2,n3),a4(n1,n2,n3),a5(n1,n2,n3),a6(n1,n2,n3),a7(n1,n2,n3),a8(n1,n2,n3))
   allocate(b1(n2,n3),b2(n2,n3),b3(n2,n3),b4(n2,n3),b5(n2,n3),b6(n2,n3),i1(n2,n3),i2(n2,n3),i3(n2,n3))
   a1=1; a2=2; a3=3; a4=4; a5=5; a6=6; a7=7; a8=8; b1=1; b2=1; b3=1; b4=1; b5=1; b6=1; i1=1; i2=2; i3=3
   nloc = n3; ob = 1; sc = 0.5d0
   !$omp target enter data map(to:a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)
   ! warm up
   do r = 1, 10
      !$omp target teams distribute parallel do collapse(2) private(i) firstprivate(ob, nloc, sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + b1(k,ot)*i1(k,ot)
            end do
         end do
      end do
   end do
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do collapse(2) private(i) firstprivate(ob, nloc, sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + b1(k,ot)*i1(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'small kernel (4 arrays):        ', (t1-t0)/nrep*1e6, ' us/launch'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do collapse(2) private(i) firstprivate(ob, nloc, sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + a4(i,k,ot) + a5(i,k,ot) + a6(i,k,ot) + a7(i,k,ot) + a8(i,k,ot) &
                  + b1(k,ot)*i1(k,ot) + b2(k,ot)*i2(k,ot) + b3(k,ot)*i3(k,ot) + b4(k,ot) + b5(k,ot) + b6(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'kernel with 17 mapped arrays:   ', (t1-t0)/nrep*1e6, ' us/launch'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do firstprivate(nloc)
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + 1.0d0
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'trivial kernel (64 iters):      ', (t1-t0)/nrep*1e6, ' us/launch'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target update to(b1(1:1,1))
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'target update to (8 bytes):     ', (t1-t0)/nrep*1e6, ' us'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target enter data map(alloc:b6)
      !$omp target exit data map(delete:b6)
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'enter/exit data alloc 37 KB:    ', (t1-t0)/nrep*1e6, ' us'
   t0 = omp_get_wtime()
   do r = 1, 200
      !$omp target enter data map(alloc:a8)
      !$omp target exit data map(delete:a8)
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'enter/exit data alloc 28 MB:    ', (t1-t0)/200*1e6, ' us'
   deallocate(a8); allocate(a8(n1*30,n2,n3))
   t0 = omp_get_wtime()
   do r = 1, 50
      !$omp target enter data map(alloc:a8)
      !$omp target exit data map(delete:a8)
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'enter/exit data alloc 850 MB:   ', (t1-t0)/50*1e6, ' us'
end program lb
