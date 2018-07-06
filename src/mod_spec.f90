      module spec
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use index

      implicit none

      public

      integer,dimension(:),allocatable,save :: jtwv,jtwvp

      real,dimension(:,:),allocatable,save :: uzm

      real,dimension(:,:,:,:),allocatable,save :: vornow,divnow,temnow,qnow,             &
                                                  vorold,divold,temold,qold,trefs,       &
                                                  vorten,divten,temten,hldten
!!                                                vorten,divten,temten,qten,hldten

      real,dimension(:,:,:),  allocatable,save :: plnow,plold,dsqgeo,spgeo,plten

      real,dimension(:,:),  allocatable,save :: plnowL,ploldL,pltenL   !  for 2dMPI, allocated in cons.f90


      contains 

         subroutine allocate_spec_array

           integer  ierr

           allocate (vornow(levp,2,jtrun,jtmax),divnow(levp,     2,jtrun,jtmax), &
                     temnow(levp,2,jtrun,jtmax),  qnow(levp*ncld,2,jtrun,jtmax), &
                     vorold(levp,2,jtrun,jtmax),divold(levp,     2,jtrun,jtmax), &
                     temold(levp,2,jtrun,jtmax),  qold(levp*ncld,2,jtrun,jtmax), &
                      trefs(levp,2,jtrun,jtmax),                                 &
                     vorten(levp,2,jtrun,jtmax),divten(levp,     2,jtrun,jtmax), &
!!                   temten(levp,2,jtrun,jtmax),  qten(levp*ncld,2,jtrun,jtmax), &
                     temten(levp,2,jtrun,jtmax),                                 &
                     hldten(levp,2,jtrun,jtmax),                                 &
                      plnow(jtrun,jtmax,2),     plold(jtrun,jtmax,2),            &
                     dsqgeo(jtrun,jtmax,2),     spgeo(jtrun,jtmax,2),            &
                     plten(jtrun,jtmax,2), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_spec : allocate fail 1 '
               stop
           end if

!CWB2017 2dMPI
           allocate (jtwv(jtrun*jtmax), jtwvp((jtrun*jtmax/nsizex)+1), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_spec : allocate fail 2 '
               stop
           end if

!CWB2015
           vornow=0.
           divnow=0.
           temnow=0.
           qnow=0.
           vorold=0.
           divold=0.
           temold=0.
           qold=0.
           trefs=0.
           vorten=0.
           divten=0.
           temten=0.
!!         qten=0.
           hldten=0.
           plnow=0.
           plold=0.
           dsqgeo=0.
           spgeo=0.

           allocate (uzm(my,lev), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_spec : allocate fail 2 '
               stop
           end if

           return

         end subroutine

         subroutine deallocate_spec_array

           deallocate (vornow,divnow,temnow,qnow,             &
                       vorold,divold,temold,qold,trefs,       &
                       vorten,divten,temten,hldten,           &
!!                     vorten,divten,temten,qten,hldten,      &
                       plnow,plold,dsqgeo,spgeo,plten)

           deallocate (uzm)
           deallocate (jtwv,jtwvp)

           return

         end subroutine

      end module spec
