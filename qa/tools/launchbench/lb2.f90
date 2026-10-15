! coalescing microbenchmark: collapse(2)+serial i (translated-code pattern)
! vs collapse(3) vs teams distribute + parallel do over i
program lb2
   use omp_lib
   implicit none
   integer, parameter :: n1 = 769, n2 = 72, n3 = 64, nrep = 500
   real(8), allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:), a4(:,:,:), a5(:,:,:), a6(:,:,:), a7(:,:,:), a8(:,:,:)
   real(8), allocatable :: b1(:,:), b2(:,:), b3(:,:)
   integer, allocatable :: i1(:,:), i2(:,:), i3(:,:)
   integer :: ot, k, i, r
   real(8) :: t0, t1, sc
   allocate(a1(n1,n2,n3),a2(n1,n2,n3),a3(n1,n2,n3),a4(n1,n2,n3),a5(n1,n2,n3),a6(n1,n2,n3),a7(n1,n2,n3),a8(n1,n2,n3))
   allocate(b1(n2,n3),b2(n2,n3),b3(n2,n3),i1(n2,n3),i2(n2,n3),i3(n2,n3))
   a1=1; a2=2; a3=3; a4=4; a5=5; a6=6; a7=7; a8=8; b1=1; b2=1; b3=1; i1=1; i2=2; i3=3; sc=0.5d0
   !$omp target enter data map(to:a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,i1,i2,i3)
   do r = 1, 5
      !$omp target teams distribute parallel do collapse(2) private(i) firstprivate(sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + a4(i,k,ot) + a5(i,k,ot) + a6(i,k,ot) + a7(i,k,ot) + a8(i,k,ot) &
                  + b1(k,ot)*i1(k,ot) + b2(k,ot)*i2(k,ot) + b3(k,ot)*i3(k,ot)
            end do
         end do
      end do
   end do
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do collapse(2) private(i) firstprivate(sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + a4(i,k,ot) + a5(i,k,ot) + a6(i,k,ot) + a7(i,k,ot) + a8(i,k,ot) &
                  + b1(k,ot)*i1(k,ot) + b2(k,ot)*i2(k,ot) + b3(k,ot)*i3(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'collapse(2) + serial i:            ', (t1-t0)/nrep*1e6, ' us'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do collapse(3) firstprivate(sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + a4(i,k,ot) + a5(i,k,ot) + a6(i,k,ot) + a7(i,k,ot) + a8(i,k,ot) &
                  + b1(k,ot)*i1(k,ot) + b2(k,ot)*i2(k,ot) + b3(k,ot)*i3(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'collapse(3):                       ', (t1-t0)/nrep*1e6, ' us'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute collapse(2) firstprivate(sc)
      do ot = 1, n3
         do k = 1, n2
            !$omp parallel do
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + a3(i,k,ot) + a4(i,k,ot) + a5(i,k,ot) + a6(i,k,ot) + a7(i,k,ot) + a8(i,k,ot) &
                  + b1(k,ot)*i1(k,ot) + b2(k,ot)*i2(k,ot) + b3(k,ot)*i3(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', 'teams distribute + parallel do i:  ', (t1-t0)/nrep*1e6, ' us'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do collapse(2) private(i) firstprivate(sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + b1(k,ot)*i1(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', '2-array collapse(2)+serial i:      ', (t1-t0)/nrep*1e6, ' us'
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do collapse(3) firstprivate(sc)
      do ot = 1, n3
         do k = 1, n2
            do i = 1, n1
               a1(i,k,ot) = a2(i,k,ot)*sc + b1(k,ot)*i1(k,ot)
            end do
         end do
      end do
   end do
   t1 = omp_get_wtime()
   print '(a,f8.1,a)', '2-array collapse(3):               ', (t1-t0)/nrep*1e6, ' us'
end program lb2
