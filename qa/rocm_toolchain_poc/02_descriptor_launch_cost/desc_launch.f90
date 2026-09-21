! PoC 02 - launch overhead grows ~10 us per allocatable / assumed-shape array
! referenced in a target region, because flang re-maps every array
! descriptor ("Copying data ... Size=88") on every launch even though the
! data is already present.  Explicit-shape dummies carry no descriptor and
! are almost free.
!
! The kernels are trivial (64 iterations) so the number is pure launch cost.
module kernels
   implicit none
contains
   ! 17 explicit-shape dummies: no descriptor, only a base pointer per array
   subroutine k_explicit(n1, n2, n3, sc, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6, i1,i2,i3)
      integer, intent(in) :: n1, n2, n3
      real(8), intent(in) :: sc
      real(8) :: a1(n1,n2,n3),a2(n1,n2,n3),a3(n1,n2,n3),a4(n1,n2,n3),a5(n1,n2,n3),a6(n1,n2,n3),a7(n1,n2,n3),a8(n1,n2,n3)
      real(8) :: b1(n2,n3),b2(n2,n3),b3(n2,n3),b4(n2,n3),b5(n2,n3),b6(n2,n3)
      integer :: i1(n2,n3),i2(n2,n3),i3(n2,n3)
      integer :: ot
      !$omp target teams distribute parallel do
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + b3(1,ot) + b4(1,ot) + b5(1,ot) + b6(1,ot) + i1(1,ot) + i2(1,ot) + i3(1,ot) &
            + a1(1,1,ot) + a2(1,1,ot) + a3(1,1,ot) + a4(1,1,ot) + a5(1,1,ot) + a6(1,1,ot) + a7(1,1,ot) + a8(1,1,ot)*sc
      end do
   end subroutine
   ! same kernel, 17 assumed-shape dummies: each carries a descriptor
   subroutine k_assumed(n3, sc, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6, i1,i2,i3)
      integer, intent(in) :: n3
      real(8), intent(in) :: sc
      real(8) :: a1(:,:,:),a2(:,:,:),a3(:,:,:),a4(:,:,:),a5(:,:,:),a6(:,:,:),a7(:,:,:),a8(:,:,:)
      real(8) :: b1(:,:),b2(:,:),b3(:,:),b4(:,:),b5(:,:),b6(:,:)
      integer :: i1(:,:),i2(:,:),i3(:,:)
      integer :: ot
      !$omp target teams distribute parallel do
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + b3(1,ot) + b4(1,ot) + b5(1,ot) + b6(1,ot) + i1(1,ot) + i2(1,ot) + i3(1,ot) &
            + a1(1,1,ot) + a2(1,1,ot) + a3(1,1,ot) + a4(1,1,ot) + a5(1,1,ot) + a6(1,1,ot) + a7(1,1,ot) + a8(1,1,ot)*sc
      end do
   end subroutine
   ! 2 explicit-shape dummies: the floor for one launch on this system
   subroutine k_explicit2(n2, n3, b1, b2)
      integer, intent(in) :: n2, n3
      real(8) :: b1(n2,n3), b2(n2,n3)
      integer :: ot
      !$omp target teams distribute parallel do
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + 1.0d0
      end do
   end subroutine
end module

program desc_launch
   use omp_lib
   use kernels
   implicit none
   integer, parameter :: n1 = 769, n2 = 72, n3 = 64
   integer :: nrep
   character(len=16) :: arg
   real(8), allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:), a4(:,:,:), a5(:,:,:), a6(:,:,:), a7(:,:,:), a8(:,:,:)
   real(8), allocatable :: b1(:,:), b2(:,:), b3(:,:), b4(:,:), b5(:,:), b6(:,:)
   integer, allocatable :: i1(:,:), i2(:,:), i3(:,:)
   integer :: ot, r
   real(8) :: t0, sc, us_expl2, us_expl17, us_assumed17, us_alloc17, us_alloc17_present
   allocate(a1(n1,n2,n3),a2(n1,n2,n3),a3(n1,n2,n3),a4(n1,n2,n3),a5(n1,n2,n3),a6(n1,n2,n3),a7(n1,n2,n3),a8(n1,n2,n3))
   allocate(b1(n2,n3),b2(n2,n3),b3(n2,n3),b4(n2,n3),b5(n2,n3),b6(n2,n3),i1(n2,n3),i2(n2,n3),i3(n2,n3))
   a1=1; a2=2; a3=3; a4=4; a5=5; a6=6; a7=7; a8=8; b1=1; b2=1; b3=1; b4=1; b5=1; b6=1; i1=1; i2=2; i3=3
   sc = 0.5d0
   nrep = 2000
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg); read(arg, *) nrep
   end if
   !$omp target enter data map(to:a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)

   do r = 1, 20   ! warm-up all four paths
      call k_explicit2(n2, n3, b1, b2)
      call k_explicit(n1, n2, n3, sc, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6, i1,i2,i3)
      call k_assumed(n3, sc, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6, i1,i2,i3)
   end do

   t0 = omp_get_wtime()
   do r = 1, nrep
      call k_explicit2(n2, n3, b1, b2)
   end do
   us_expl2 = (omp_get_wtime() - t0)/nrep*1e6

   t0 = omp_get_wtime()
   do r = 1, nrep
      call k_explicit(n1, n2, n3, sc, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6, i1,i2,i3)
   end do
   us_expl17 = (omp_get_wtime() - t0)/nrep*1e6

   t0 = omp_get_wtime()
   do r = 1, nrep
      call k_assumed(n3, sc, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6, i1,i2,i3)
   end do
   us_assumed17 = (omp_get_wtime() - t0)/nrep*1e6

   ! allocatables referenced directly from the program scope (same as k_assumed)
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do firstprivate(sc)
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + b3(1,ot) + b4(1,ot) + b5(1,ot) + b6(1,ot) + i1(1,ot) + i2(1,ot) + i3(1,ot) &
            + a1(1,1,ot) + a2(1,1,ot) + a3(1,1,ot) + a4(1,1,ot) + a5(1,1,ot) + a6(1,1,ot) + a7(1,1,ot) + a8(1,1,ot)*sc
      end do
   end do
   us_alloc17 = (omp_get_wtime() - t0)/nrep*1e6

   ! ... and with an explicit map(present,alloc:) - does not help
   t0 = omp_get_wtime()
   do r = 1, nrep
      !$omp target teams distribute parallel do firstprivate(sc) &
      !$omp&   map(present,alloc: a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,i1,i2,i3)
      do ot = 1, n3
         b1(1,ot) = b2(1,ot) + b3(1,ot) + b4(1,ot) + b5(1,ot) + b6(1,ot) + i1(1,ot) + i2(1,ot) + i3(1,ot) &
            + a1(1,1,ot) + a2(1,1,ot) + a3(1,1,ot) + a4(1,1,ot) + a5(1,1,ot) + a6(1,1,ot) + a7(1,1,ot) + a8(1,1,ot)*sc
      end do
   end do
   us_alloc17_present = (omp_get_wtime() - t0)/nrep*1e6

   print '(a)', 'us per launch (trivial kernel, data already resident):'
   print '(a,f8.1)', '   2 explicit-shape dummies            ', us_expl2
   print '(a,f8.1)', '  17 explicit-shape dummies            ', us_expl17
   print '(a,f8.1)', '  17 assumed-shape dummies (descriptor)', us_assumed17
   print '(a,f8.1)', '  17 allocatables in scope             ', us_alloc17
   print '(a,f8.1)', '  17 allocatables + map(present,alloc:)', us_alloc17_present
   print '(a,f6.1,a)', 'descriptor cost per array: ~', (us_assumed17 - us_expl17)/17, ' us'
end program
