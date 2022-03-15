      module phygrid
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use index

      implicit none

      public

      real, dimension(:,:),allocatable,save ::        &
             snr       ,   gwr       ,     tg       , &
                            ss       ,     rs       , &
           ustar       , tstar       ,  qstar       , &
           hflux       , qflux       ,raintot       , &
          raincu       ,rainlp       , totalp       , &
          curate       ,  plcl       , cumtop       , &
          tgclim       ,  gwet       ,     z0       , &
             alb       ,gwclim       ,   acld       , &
            ctot       ,  chig       ,   cmid       , &
            clow       ,  hpbl       ,   cosz       , & 
         rainlp6       ,raincu6      ,                &
         rainlp3       ,raincu3      ,                &
         rainlp1       ,raincu1      ,tsflw

      logical, allocatable,save :: land(:,:),ice(:,:),ocean(:,:)

      integer, allocatable,save :: il(:,:),ib(:,:)

      real, dimension(:,:),allocatable,save :: cof
      real, dimension(:,:),allocatable,save :: xlon
      real, dimension(:)  ,allocatable,save :: xlat

      real, dimension(:,:),allocatable,save :: u10,v10,t2,rh2,rh10,q2  &
                                              ,fm,fm10,fh,fh2,srflag
 
      real, dimension(:,:),allocatable,save :: fpsp,fpsp1

      real, dimension(:,:,:),allocatable,save :: e,eps,o3l,dtrad,asl,atl
      real, dimension(:,:,:),allocatable,save :: ftp,fqp,ftp1,fqp1
      real, dimension(:,:,:),allocatable,save :: deltaq,cnvwr,cnvcr

      contains 

         subroutine allocate_phygrid_array

           integer  ierr

           allocate (  e(nxp,lev,my_max),  eps(nxp,lev,my_max),  &
                     o3l(nxp,lev,my_max),dtrad(nxp,lev,my_max),  &
                     asl(nxp,lev,my_max),  atl(nxp,lev,my_max),  &
                     ftp(nxp,lev,my_max),  fqp(nxp,lev,my_max),  &
                    ftp1(nxp,lev,my_max), fqp1(nxp,lev,my_max),  &
                    deltaq(nxp,lev,my_max),cnvwr(nxp,lev,my_max),&
                    cnvcr(nxp,lev,my_max),  stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 1 '
               stop
           end if
           deltaq = 0.
           cnvcr  = 0.
           cnvwr  = 0.

           allocate (                                 &
             snr(nxp,my_max),   gwr(nxp,my_max),     tg(nxp,my_max), &
                            ss(nxp,my_max),     rs(nxp,my_max), &
           ustar(nxp,my_max), tstar(nxp,my_max),  qstar(nxp,my_max), &
           hflux(nxp,my_max), qflux(nxp,my_max),raintot(nxp,my_max), &
          raincu(nxp,my_max),rainlp(nxp,my_max), totalp(nxp,my_max), &
          curate(nxp,my_max),  plcl(nxp,my_max), cumtop(nxp,my_max), &
          tgclim(nxp,my_max),  gwet(nxp,my_max),     z0(nxp,my_max), &
             alb(nxp,my_max),gwclim(nxp,my_max),  acld(lev,my), &
            ctot(nxp,my_max),  chig(nxp,my_max),   cmid(nxp,my_max), &
            clow(nxp,my_max),  hpbl(nxp,my_max),   cosz(nxp,my_max), &
         rainlp6(nxp,my_max),raincu6(nxp,my_max),                    &
         rainlp3(nxp,my_max),raincu3(nxp,my_max),                    &
         rainlp1(nxp,my_max),raincu1(nxp,my_max), tsflw(nxp,my_max), &
                                          stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 2 '
               stop
           end if

!CWB2015
           ss=0.
           rs=0.

!CWB2015
           ctot=0.
           chig=0.
           cmid=0.
           clow=0.
           hpbl=0.

!CWB2016
           curate=0.
           tsflw=0.

           allocate (land(nxp,my_max),ice(nxp,my_max), &
                     ocean(nxp,my_max), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 3 '
               stop
           end if

           allocate (il(nxp,4),ib(nxp,4), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 4 '
               stop
           end if


           allocate (cof(nxp*3,4),xlon(nx,my_max),xlat(my), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 5 '
               stop
           end if

!CWB2015
           il=0
           ib=0
           cof=0.

           allocate (u10(nxp,my_max),v10(nxp,my_max),srflag(nxp,my_max) &
                     ,t2(nxp,my_max),rh2(nxp,my_max),rh10(nxp,my_max)   &
                     ,q2(nxp,my_max),fm(nxp,my_max),fm10(nxp,my_max)    &
                     ,fh(nxp,my_max),fh2(nxp,my_max), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 6 '
               stop
           end if
!
           rh2=0.
           rh10=0.
           gwr=0.
!
           allocate (fpsp(nxp,my_max),fpsp1(nxp,my_max),stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 7 '
               stop
           end if
!
           return

         end subroutine

         subroutine deallocate_phygrid_array

           deallocate (e,eps,o3l,dtrad,asl,atl,ftp,fqp,ftp1,fqp1)
           deallocate (deltaq,cnvwr,cnvcr)
           deallocate (                                           &
                       snr,gwr,tg,   ss,rs,                       &
             ustar,tstar,qstar,hflux,qflux,raintot,raincu,rainlp, &
             totalp,curate,plcl,cumtop,tgclim,gwet,z0,alb,        &
             gwclim,acld,ctot,chig,cmid,clow,hpbl,cosz)

           deallocate (land,ice,ocean)
           deallocate (il,ib)
           deallocate (cof,xlon,xlat)
           deallocate (u10,v10,t2,rh2,rh10,srflag,q2,fm,fm10,fh,fh2)
           deallocate (fpsp,fpsp1)
           deallocate (rainlp6,raincu6,rainlp3,raincu3,rainlp1,raincu1)
           deallocate (tsflw)

           return

         end subroutine

      end module phygrid
