      module grid
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use index
      use const, only : RTYPE

      implicit none

      public

      real, dimension(:,:),allocatable,save ::                      &
                           pdiff,t1000,tsave,std
      real(kind=RTYPE),dimension(:,:),allocatable,save :: pt,sgeo,  &
                           ptend,ptp

      real, dimension(:,:,:),allocatable,save :: plt
!! for Semi-Lagrangian
      real(kind=RTYPE), dimension(:,:,:),allocatable,save :: rdivm
      real(kind=RTYPE), dimension(:,:),allocatable,save :: dlpl,dtpl
      real(kind=RTYPE), dimension(:,:,:),allocatable,save ::        &
                           ut,vt,tt,qt,up,vp,ttp,qp,qm,sd,vvel,     &
                           ut_sl,vt_sl,phi,dlphi,dtphi,rvor,rdiv,   &
                           pk,pk2

!!                             ut_sl,vt_sl,uum_sl,vvm_sl,ttm_sl
!!      real, dimension(:,:,:,:),allocatable,save :: qm_sl
!!    real, dimension(:,:),allocatable,save :: pt_sl,ptp_sl

      integer   lonfull,lonhalf,lonpart,lonlenmax,mylonlen
      integer   latfull,lathalf,latpart,latlenmax,mylatlen
      integer   ndslhvar,ndslvvar,xy

      integer, allocatable :: lonstr(:),lonlen(:)
      integer, allocatable :: latstr(:),latlen(:)
      real(kind=RTYPE), allocatable :: gslati(:),gglati(:)
      real(kind=RTYPE), allocatable :: fa1(:),fa2(:),fa3(:),fa4(:)
      contains

         subroutine allocate_grid_array

           integer  ierr

           allocate (ut(nxp,lev,my_max),  &
                     vt(nxp,lev,my_max),  &
                     sd(nxp,lev,my_max),  &
                   vvel(nxp,lev,my_max),  &
                   rvor(nxp,lev,my_max),  &
                   rdiv(nxp,lev,my_max),  &
                     tt(nxp,lev,my_max),  &
                  qt(nxp,lev*ncld,my_max),&
                    phi(nxp,lev,my_max),  &
                    plt(nxp,lev,my_max),  &
                     pk(nxp,lev,my_max),  &
                    pk2(nxp,lev,my_max),  &
                     up(nxp,lev,my_max),  &
                     vp(nxp,lev,my_max),  &
                    ttp(nxp,lev,my_max),  &
                  rdivm(nxp,lev,my_max),  &
                  qp(nxp,lev*ncld,my_max),&
                  qm(nxp,lev*ncld,my_max),&
!! for Semi-Lagrangian
                  ut_sl(nx,levp,my_max),      &
                  vt_sl(nx,levp,my_max),      &
!!                 uum_sl(nx,levp,my_max),      &
!!                 vvm_sl(nx,levp,my_max),      &
!!                 ttm_sl(nx,levp,my_max),      &
!!                  qm_sl(nx,levp,ncld,my_max), &
!!                pt_sl(nx,my_max),           &
!!                  ptp_sl(nx,my_max),           &
                 stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_grid : allocate fail 1 '
               stop
           end if
           ut=0. ;  vt=0.;  sd=0.; rvor=0.; rdiv=0.
           tt=0. ;  qt=0.; phi=0.;  plt=0.;   pk=0.
           pk2=0.;  up=0.;  vp=0.;  ttp=0.;   qp=0.

           allocate (pt(nxp,my_max),dlpl(nxp,my_max),dtpl(nxp,my_max), &
                    sgeo(nxp,my_max), pdiff(nxp,my_max),&
                    ptend(nxp,my_max), t1000(nxp,my_max), tsave(nxp,my_max), &
                    std(nxp,my_max), ptp(nxp,my_max) ,stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_grid : allocate fail 2 '
               stop
           end if
           sgeo=0.; std=0.; pdiff=0.; tsave=0.; pt=0.

           allocate (gslati(my*2+1),gglati(my*2+1),                 &
                     lonstr(npe),lonlen(npe),                       &
                     latstr(npe),latlen(npe), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_grid for ndsl : allocate fail 3'
               stop
           end if

           allocate (dlphi(nxp,lev,my_max),  &
                     dtphi(nxp,lev,my_max),stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_grid for ndsl : allocate fail 4'
               stop
           end if

           allocate (fa1(my*2),fa2(my*2),fa3(my*2),fa4(my*2), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_grid for ndsl : allocate fail 5'
               stop
           end if
           dlphi=0.
           dtphi=0.

!CWB2015
           ptend=0.

!CWB2018
           ut=0.
           vt=0.
           up=0.
           vp=0.
           tt=0.
           ttp=0.
           qt=0.
           qp=0.
!! for Semi-Lagrangiain
           rdivm=0.
           ut_sl=0.
           vt_sl=0.
           qm=0.
           sd=0.
           vvel=0.
           fa1=0.
           fa2=0.
           fa3=0.
           fa4=0.

           return

         end subroutine

         subroutine deallocate_grid_array

           deallocate (ut,vt,sd,vvel,rvor,rdiv,tt,qt,phi,plt,pk,pk2,up,vp,ttp,qp,qm)
           deallocate (rdivm)
           deallocate ( pt,dlpl,dtpl,sgeo,pdiff, &
                ptend,t1000,tsave,std,ptp)
           ! for Semi-Lagrangian
           deallocate (gslati,gglati,lonstr,lonlen,latstr,latlen)
           deallocate (fa1,fa2,fa3,fa4)
           deallocate (dlphi,dtphi)
!!         deallocate (ut_sl,vt_sl,uum_sl,vvm_sl,ttm_sl,qm_sl,pt_sl,ptp_sl)
           deallocate (ut_sl,vt_sl)

           return

         end subroutine

      end module grid
