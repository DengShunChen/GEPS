      module phygrid
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use index
      use const, only: RTYPE, naero

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
         rainlp6       ,raincu6      ,tg_diff       , &
         rainlp3       ,raincu3      ,tg_ocn        , &
         rainlp1       ,raincu1      ,tsflw

      logical, allocatable,save :: land(:,:),ice(:,:),ocean(:,:)

      integer, allocatable,save :: il(:,:),ib(:,:)

!helio>
      integer, dimension(:,:),allocatable,save :: ls_full(:,:)
      integer, dimension(:,:),allocatable,save :: ls_redu(:,:)
      real, dimension(:,:,:),allocatable,save :: outp(:,:,:)
!helio<

      real, dimension(:,:),allocatable,save :: cof
      real, dimension(:,:),allocatable,save :: xlon
      real, dimension(:)  ,allocatable,save :: xlat

      real, dimension(:,:),allocatable,save :: u10,v10,t2,rh2,rh10,q2  &
#ifdef TIMCOMCPL
                                              ,fm,fm10,fh,fh2,srflag   &
                                              ,ustress,vstress,ssu,ssv &
                                              ,ifrac,icedp,snodp
#else
                                              ,fm,fm10,fh,fh2,srflag
#endif
 
      real, dimension(:,:),allocatable,save :: fpsp,fpsp1

      real, dimension(:,:,:),allocatable,save :: e,eps,dtrad,asl,atl
      real, dimension(:,:,:),allocatable,save :: ftp,fqp,ftp1,fqp1
      real, dimension(:,:,:),allocatable,save :: deltaq,cnvwr,cnvcr
      real, dimension(:,:,:),allocatable,save :: dtcup,ducup,dvcup,    &
                                                 dtshl,dushl,dvshl,    &
                                                 dtlsp,dulsp,dvlsp
      real(kind=RTYPE), dimension(:,:,:),allocatable,save :: o3l
      real(kind=RTYPE), dimension(:,:,:),allocatable,save :: aeroclxm
      real(kind=RTYPE), dimension(:,:,:,:),allocatable,save ::         &
                                                 aerosave1,aerosave2

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
           e    = 0.
           eps  = 0.
           o3l  = 0.
           dtrad= 0.
           asl  = 0.
           atl  = 0.
           ftp  = 0.
           fqp  = 0.
           ftp1 = 0.
           fqp1 = 0

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
         rainlp6(nxp,my_max),raincu6(nxp,my_max),tg_diff(nxp,my_max), &
         rainlp3(nxp,my_max),raincu3(nxp,my_max),tg_ocn(nxp,my_max), &
         rainlp1(nxp,my_max),raincu1(nxp,my_max), tsflw(nxp,my_max), &
                                          stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 2 '
               stop
           end if

!CWB2015
           ss=0.
           rs=0.

!CWB2016
           curate=0.
           tsflw=0.
           cumtop=0.
           plcl  =0.
             snr=0.;        gwr=0.;          tg=0.
           ustar=0.;      tstar=0.;       qstar=0.
           hflux=0.;      qflux=0.;     raintot=0.
          raincu=0.;     rainlp=0.;      totalp=0.
          tgclim=0.;       gwet=0.;          z0=0.
             alb=0.;     gwclim=0.;        acld=0.
            ctot=0.;       chig=0.;        cmid=0.
            clow=0.;       hpbl=0.;        cosz=0.
         rainlp6=0.;    raincu6=0.
         rainlp3=0.;    raincu3=0.
         rainlp1=0.;    raincu1=0.

           allocate (land(nxp,my_max),ice(nxp,my_max), &
                     ocean(nxp,my_max), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 3 '
               stop
           end if
           land=.false.; ice=.false.; ocean=.false.

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
           xlon=0.
           xlat=0.

           allocate (u10(nxp,my_max),v10(nxp,my_max),srflag(nxp,my_max) &
                     ,t2(nxp,my_max),rh2(nxp,my_max),rh10(nxp,my_max)   &
                     ,q2(nxp,my_max),fm(nxp,my_max),fm10(nxp,my_max)    &
#ifdef TIMCOMCPL
                     ,fh(nxp,my_max),fh2(nxp,my_max)                    &
                     ,ustress(nxp,my_max),vstress(nxp,my_max)           &
                     ,ssu(nxp,my_max),ssv(nxp,my_max)                   &
                     ,ifrac(nxp,my_max),icedp(nxp,my_max),snodp(nxp,my_max) &
                     ,stat=ierr)
#else
                     ,fh(nxp,my_max),fh2(nxp,my_max), stat=ierr)
#endif
           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 6 '
               stop
           end if
!
           u10=0.;     v10=0.;     t2=0.;     rh2=0.
          rh10=0.;      q2=0.;     fm=0.;   fm10=0.
            fh=0.;     fh2=0.; srflag=0.
           gwr=0.;     
#ifdef TIMCOMCPL
           ssv=0.;    ssu=0.
           ustress=0.; vstress=0.
           ifrac=0.;   icedp=0.;  snodp=0.
#endif

!
           allocate (fpsp(nxp,my_max),fpsp1(nxp,my_max),stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 7 '
               stop
           end if
           fpsp=0.;    fpsp1=0.
!
           allocate (dtcup(nxp,lev,my_max),ducup(nxp,lev,my_max),    &
                     dvcup(nxp,lev,my_max),dtshl(nxp,lev,my_max),    &
                     dushl(nxp,lev,my_max),dvshl(nxp,lev,my_max),    &
                     dtlsp(nxp,lev,my_max),dulsp(nxp,lev,my_max),    &
                     dvlsp(nxp,lev,my_max),  stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 8 '
               stop
           end if
!
           dtcup = 0.
           ducup = 0.
           dvcup = 0.
           dtshl = 0.
           dushl = 0.
           dvshl = 0.
           dtlsp = 0.
           dulsp = 0.
           dvlsp = 0.
!

!helio>
           allocate (ls_full(nx,my_max),ls_redu(nx,my_max), stat=ierr)
           allocate (outp(nx,my_max,8), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_phygrid : allocate fail 9 '
               stop
           end if
           ls_full = 0.
           ls_redu = 0.
           outp = 0.
!helio<

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
           deallocate (dtcup,ducup,dvcup,dtshl,dushl,dvshl,dtlsp,dulsp,dvlsp)
!helio>
           deallocate (ls_full,ls_redu)
           deallocate (outp)
!helio<

           return

         end subroutine

         subroutine allocate_aerogrid_array
           integer  ierr
           allocate ( aeroclxm(nxp,naero*lev,my_max), &
                      aerosave1(nxp,lev,my_max,naero), &
                      aerosave2(nxp,lev,my_max,naero),stat=ierr )
           if (ierr/= 0) stop 'mod_phygrid : allocate aeroclxm'
           aeroclxm = 0.
           aerosave1 = 0.
           aerosave2 = 0.
           return
         end subroutine

         subroutine deallocate_aerogrid_array
           deallocate ( aeroclxm,aerosave1,aerosave2 )
           return
         end subroutine

      end module phygrid
