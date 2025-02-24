      module spec
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use const, only : RTYPE
      use param
      use index

      implicit none

      public

      integer,dimension(:),allocatable,save :: jtwv,jtwvp


!byl      real,dimension(:,:,:,:),allocatable,save :: vornow,divnow,temnow,qnow,             &
!byl                                                  vorold,divold,temold,qold,trefs,       &
      real(kind=RTYPE),dimension(:,:,:,:),allocatable,save :: vornow,divnow,temnow,         &
                                                              vormid,divmid,temmid,         &
                                                              vorold,divold,temold,trefs,   &
                                                              vorten,divten,temten,hldten
!!                                                            vorten,divten,temten,qten,hldten
!byl      real,dimension(:,:,:),  allocatable,save :: plnow,plold,dsqgeo,spgeo,plten
      real(kind=RTYPE),dimension(:,:,:),  allocatable,save :: plnow,plmid,plold,plten,spgeo

!ch   real,dimension(:,:),  allocatable,save :: plnowL,ploldL,pltenL   !  for 2dMPI, allocated in cons.f90
      real(kind=RTYPE),dimension(:,:),  allocatable,save :: plnowL,ploldL,pltenL   !  for 2dMPI, allocated in cons.f90


      contains 

         subroutine allocate_spec_array

           integer  ierr

           allocate (vornow(levp,2,jtrun,jtmax),divnow(levp,     2,jtrun,jtmax), &
                     temnow(levp,2,jtrun,jtmax),temold(levp     ,2,jtrun,jtmax), &
!byl                     qnow(levp*ncld,2,jtrun,jtmax),qold(levp*ncld,2,jtrun,jtmax),& 
                     vorold(levp,2,jtrun,jtmax),divold(levp,     2,jtrun,jtmax), &
                     vormid(levp,2,jtrun,jtmax),divmid(levp,     2,jtrun,jtmax), &
                     temmid(levp,2,jtrun,jtmax),                                 &
                      trefs(levp,2,jtrun,jtmax),                                 &
                     vorten(levp,2,jtrun,jtmax),divten(levp,     2,jtrun,jtmax), &
!!                   temten(levp,2,jtrun,jtmax),  qten(levp*ncld,2,jtrun,jtmax), &
                     temten(levp,2,jtrun,jtmax),                                 &
                     hldten(levp,2,jtrun,jtmax),                                 &
                      plnow(jtrun,jtmax,2),     plold(jtrun,jtmax,2),            &
                     plmid(jtrun,jtmax,2),                                       &
!byl                     dsqgeo(jtrun,jtmax,2),     spgeo(jtrun,jtmax,2),            &
                     spgeo(jtrun,jtmax,2),                                       &
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
           vormid=0.
           divmid=0.
           temmid=0.
!byl           qnow=0.
           vorold=0.
           divold=0.
           temold=0.
!byl           qold=0.
           trefs=0.
           vorten=0.
           divten=0.
           temten=0.
!!         qten=0.
           hldten=0.
           plnow=0.
           plold=0.
!byl           dsqgeo=0.
           spgeo=0.
           jtwv=0.
           plmid=0.   ! avoid undefine valuse

           return

         end subroutine

         subroutine deallocate_spec_array

!byl           deallocate (vornow,divnow,temnow,qnow,             &
!byl                       vorold,divold,temold,qold,trefs,       &
           deallocate (vornow,divnow,temnow,                  &
                       vormid,divmid,temmid,                  &
                       vorold,divold,temold,trefs,            &
                       vorten,divten,temten,hldten,           &
!!                     vorten,divten,temten,qten,hldten,      &
!byl                       plnow,plold,dsqgeo,spgeo,plten)
                       plnow,plmid,plold,plten,spgeo)

           deallocate (jtwv,jtwvp)
           deallocate (plnowL,ploldL,pltenL)

           return

         end subroutine

      end module spec
