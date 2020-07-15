      module grid
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use index

      implicit none

      public

      real, dimension(:,:),allocatable,save :: pt,dlpl,dtpl,sgeo, &
                           pdiff,ptend,t1000,tsave,std,ptp

      real, dimension(:,:,:),allocatable,save :: ut,vt,sd,rvor,rdiv,&
                           tt,qt,phi,plt,pk,pk2,up,vp,ttp,qp
!! for Semi-Lagrangian
      real, dimension(:,:,:),allocatable,save :: dlphi,dtphi,     &
                             ut_sl,vt_sl
!!                             ut_sl,vt_sl,uum_sl,vvm_sl,ttm_sl
!!      real, dimension(:,:,:,:),allocatable,save :: qm_sl
!!    real, dimension(:,:),allocatable,save :: pt_sl,ptp_sl

      integer   lonfull,lonhalf,lonpart,lonlenmax,mylonlen
      integer   latfull,lathalf,latpart,latlenmax,mylatlen
      integer   ndslhvar,ndslvvar,xy

      integer, allocatable :: lonstr(:),lonlen(:)
      integer, allocatable :: latstr(:),latlen(:)
      real, allocatable :: gslati(:),gglati(:)

      contains 

         subroutine allocate_grid_array

           integer  ierr

           allocate (ut(nxp,lev,my_max),  &
                     vt(nxp,lev,my_max),  &
                     sd(nxp,lev,my_max),  &
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
                  qp(nxp,lev*ncld,my_max),&
!! for Semi-Lagrangian
                  ut_sl(nx,levp,my_max),      &
                  vt_sl(nx,levp,my_max),      &
!!                 uum_sl(nx,levp,my_max),      &
!!                 vvm_sl(nx,levp,my_max),      &
!!                 ttm_sl(nx,levp,my_max),      &
!!                  qm_sl(nx,levp,ncld,my_max), &
!!                pt_sl(nx,my_max),           &
!!               ptp_sl(nx,my_max),           &
                 stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_grid : allocate fail 1 '
               stop
           end if

           allocate (pt(nxp,my_max),dlpl(nxp,my_max),dtpl(nxp,my_max), &
                    sgeo(nxp,my_max), pdiff(nxp,my_max),&
                    ptend(nxp,my_max), t1000(nxp,my_max), tsave(nxp,my_max), &
                    std(nxp,my_max), ptp(nxp,my_max) ,stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_grid : allocate fail 2 '
               stop
           end if

           allocate (gslati(my*2+1),gglati(my*2+1),     &
                     lonstr(npe),lonlen(npe),           &
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

!CWB2015
           ptend=0.

!CWB2018
           sd=0.

           return

         end subroutine

         subroutine deallocate_grid_array

           deallocate (ut,vt,sd,rvor,rdiv,tt,qt,phi,plt,pk,pk2,up,vp,ttp,qp)
           deallocate ( pt,dlpl,dtpl,sgeo,pdiff, &
               ptend,t1000,tsave,std,ptp)
! for Semi-Lagrangian
           deallocate (gslati,gglati,lonstr,lonlen,latstr,latlen)
           deallocate (dlphi,dtphi)
!!         deallocate (ut_sl,vt_sl,uum_sl,vvm_sl,ttm_sl,qm_sl,pt_sl,ptp_sl)
           deallocate (ut_sl,vt_sl)

           return

         end subroutine

      end module grid
