module mod_outflds_gpu
 implicit none

 real,allocatable::      tens(:)
 integer,parameter,private:: async_id = 1
 
contains
  subroutine divgout(nx,my,my_max,lpout,lev,itau,idtg, &
        plev,num,whtlev,pkout,pk,pklp,rdiv,rdivb,div,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param !for write grb2
      use const ,only:outdms,outgrb2 ,RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,ncnt

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max),       &
                rdivb(nxp,my_max),               &
                plev(lpout),whtlev(num)
      real       rdiv(nxp,lev,my_max)
      real(kind=RTYPE) div(nxp,my_max,lpout)
      real(kind=RTYPE) wk1(nx,my),pout(nx,my)

!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!
      integer    k,lpl,lenc,n,num,istat
      integer:: ptp0(9),ptp1(9)
  
      logical :: lwrite
!$acc wait(async_id)
!$acc enter data create(wk1,pout) async(async_id)

      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,rdiv,rdivb,pkout,div,tens)

      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'230'
      end do
      lrec(lpout) = 'h00230'
!
!
      lenc= nx*my
      ncnt= 0
!
      do 20 n=1,num
      do 10 k=1,lpout
      if(plev(k).eq.whtlev(n)) then
      call unify_reduceintp_gpu(nx,my,my_max,div(1,1,k),wk1)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
      call qmaxn3_w(wk1,1,1,1,nx,my,1)
      ptp0=(/0,2,11,6,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
      go to 20
      endif
   10 continue
   20 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!$acc wait(async_id)
!$acc exit data delete(wk1,pout) async(async_id)
!$acc wait(async_id)
!
      return
  end subroutine divgout

  subroutine dragout(nx,my,my_max,lpout,lev,itau,idtg, &
        plev,num,whtlev,pkout,pk,pklp,rdrag,rdragb,drag,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mpe
      use mod_grb2_param !for write grb2
      use const ,only:outdms,outgrb2 ,RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ncnt

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max),               &
                rdrag(nxp,lev,my_max),rdragb(nxp,my_max),                       &
                plev(lpout),whtlev(num)
      real(kind=RTYPE) wk1(nx,my),pout(nx,my),drag(nxp,my_max,lpout)
!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!
!      real      pk_lev_m1(nx,my)
!

!
      integer   k,lpl,lenc,j,nxj,ij,jj,i,n,istat
      integer:: ptp0(9),ptp1(9)

      logical lwrite

!$acc wait(async_id)
!$acc enter data create(wk1,pout) async(async_id)

      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,rdrag,rdragb,pkout,drag,tens)
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'290'
      end do
      lrec(lpout) = 'h00290'
!
!
      lenc= nx*my
      ncnt= 0
!
      do 20 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do 45 jj=1, jlistnum
      do 45 i=1,nxp
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      if(i<=nxj)then
        if(pkout(k).gt.pk(i,lev-1,jj)) drag(i,jj,k)= 0.
      endif
   45 continue
!
      call unify_reduceintp_gpu(nx,my,my_max,drag(1,1,k),wk1)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
      call qmaxn3_w(wk1,1,1,1,nx,my,1)
      ptp0=(/0,2,196,6,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
      go to 20
      endif
   10 continue
   20 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!$acc wait(async_id)
!$acc exit data delete(wk1,pout) async(async_id)
!$acc wait(async_id)
!
      return
  end subroutine dragout

  subroutine geopout(nx,my,my_max,lpout,lev,itau,idtg,plev, &
        num,whtlev,pkout,pk,pklp,phi,phib,phips,glob,ggdef,phistd,      &
        h850,h500,lwrite)

      use index
      use rank
      use mod_grb2_param !for write grb2
      use const ,only:outdms,outgrb2 ,RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ncnt
!
      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max),            &
                phi(nxp,lev,my_max),phib(nxp,my_max),plev(lpout),whtlev(num)

      real      phistd(lpout)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),slp(nx,my),phips(nxp,my_max,lpout)
      real(kind=RTYPE) h850(nxp,my_max),h500(nxp,my_max),tmp(nxp,my_max)
!

      integer*8    idtg

      character*6  lrec(lpout)
      character*4  ggdef
!
      integer      i,k,n,lenc,istat,lpl,jj,j,nxj
      integer:: ptp0(9),ptp1(9)

      logical lwrite

!$acc wait(async_id)
!$acc enter data create(tmp,slp,pout) async(async_id)
!
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,phi,phib,pkout,phips,tens)

      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'000'
      end do
      lrec(lpout) = 'h00000'
!
!
      lenc= nx*my
      ncnt= 0
!
      do k=1,lpout
        if(plev(k).eq.850.)then
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj=1, jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj)  h850(i,jj)= phips(i,jj,k)
          enddo
        enddo
        else if(plev(k).eq.500.)then
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj=1, jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj) h500(i,jj)= phips(i,jj,k)
          enddo
        enddo
        endif
      enddo
!
      do 20 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do  jj=1, jlistnum
       do  i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) tmp(i,jj)= phips(i,jj,k)+phistd(k)
       enddo
      enddo
   11 continue
      call unify_reduceintp_gpu(nx,my,my_max,tmp,slp)
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
      call smth9(nx,my,slp,glob,1)
      call smth9(nx,my,glob,slp,2)
!
      call qmaxn3_w(slp,1,1,1,nx,my,1)
      ptp0=(/0,3,5,1,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,slp,pout,ptp0,ptp1)
!
      go to 20
      endif
   10 continue
   20 continue

      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!$acc wait(async_id)
!$acc exit data delete(tmp,slp,pout) async(async_id)
!$acc wait(async_id)
!
      return
  end subroutine geopout

  subroutine shumout(nx,my,my_max,lpout,lev,itau,idtg   &
      , plev,num,whtlev,pkout,pk,pklp,dpd,dpdb,dew,glob,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:outdms,outgrb2 ,RTYPE,kflag

      implicit  none
      integer   nx,my,my_max,lpout,lev,itau,num,ncnt

      real      pkout(lpout),pklp(nxp,my_max)                            &
      , pk(nxp,lev,my_max),dpd(nxp,lev,my_max),dpdb(nxp,my_max)          &
      , plev(lpout),whtlev(num)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),tmp(nxp,my_max)           &
      ,                dew(nxp,my_max,lpout),ffx(nx,my_max)

      integer   i,k,lpl,n,lenc,istat,jj,j,nxj
      integer:: ptp0(9),ptp1(9)
!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!
      logical :: lwrite
!$acc wait(async_id)
!$acc enter data create(pout) async(async_id)

      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,dpd,dpdb,pkout,dew,tens)

      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'510'
      end do
      lrec(lpout) = 'h00510'
!
!
      lenc= nx*my
      ncnt= 0
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
! relative humidity must be smaller or equal 1.0
!
      call unify_reduceintp_gpu(nx,my,my_max,dew(1,1,k),glob)

!$acc parallel loop collapse(2) async(async_id)
        do j=1,my
          do i=1,nx
             glob(i,j)=min(100.,max(glob(i,j)*100.,0.0))
          enddo
        enddo
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,1,1,2,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
!$acc wait(async_id)
!$acc exit data delete(pout) async(async_id)
!$acc wait(async_id)
      return
  end subroutine shumout

  subroutine shumout2(nx,my,my_max,lpout,lev,itau,idtg   &
      , plev,num,whtlev,pkout,pk,pklp,dpd,dpdb,dew,glob,ggdef,ntrac,lwrite)
!
      use index
      use rank, only : myrank
      use radn, only : ntoz
      use param, only : ncld
      use mod_grb2_param
      use const ,only:outdms,outgrb2 , RTYPE,kflag,qmin,nmmiph

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ntrac,ncnt,ntrchk

      real      pkout(lpout),pklp(nxp,my_max)                        &
      , pk(nxp,lev,my_max),dpd(nxp,lev,my_max),dpdb(nxp,my_max)      &
      , plev(lpout),whtlev(num)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),tmp(nxp,my_max)       &
      , dew(nxp,my_max,lpout),ffx(nx,my_max)
!

      integer   i,k,lpl,n,lenc,istat,jj,j,nxj
      integer:: gtp0(9),gtp1(9)

      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!      character*3 cspec(8)
      character*3,dimension(:),allocatable :: cspec
      logical :: lwrite
      integer::Ptp0,Ptp1,Ptp2,Ptp3
!      integer,dimension(6)::cspe0,cspe1,cspe2,cspe3
      integer,dimension(:),allocatable ::cspe0,cspe1,cspe2,cspe3

!$acc wait(async_id)
!$acc enter data create(pout) async(async_id)

!key=556 for mixing ratio of hail
!key=571~575 for number concentration of cloud droplet, ice, rain, snow, and graupel
!key=572 : inc (ntinc=7)
!key=573 : rnc (ntrnc=8)
      if ( nmmiph .eq. 18 ) then
        IF (.NOT. ALLOCATED( cspec )) &
        allocate ( cspec(8),cspe0(8),cspe1(8),cspe2(8),cspe3(8) )
        cspec=(/'500','551','553','552','554','555','572','573'/)
        cspe0=(/  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  /)
        cspe1=(/  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  /)
        cspe2=(/  0  , 22  , 24  , 82  , 25  , 32  , 207 , 104 /)  !not sure of inc
        cspe3=(/  6  ,  8  ,  8  ,  8  ,  8  ,  8  ,  8  ,  8  /)
      elseif ( nmmiph .eq. 16 ) then
        IF (.NOT. ALLOCATED( cspec )) &
        allocate ( cspec(7),cspe0(7),cspe1(7),cspe2(7),cspe3(7) )
        cspec=(/'500','551','553','552','554','555','556'/)
        cspe0=(/  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  /)
        cspe1=(/  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  /)
        cspe2=(/  0  , 22  , 24  , 82  , 25  , 32  , 71  /)
        cspe3=(/  6  ,  8  ,  8  ,  8  ,  8  ,  8  ,  8  /)
      else
        IF (.NOT. ALLOCATED( cspec )) &
        allocate ( cspec(6),cspe0(6),cspe1(6),cspe2(6),cspe3(6) )
        cspec=(/'500','551','553','552','554','555'/)
        cspe0=(/  0  ,  0  ,  0  ,  0  ,  0  ,  0  /)
        cspe1=(/  1  ,  1  ,  1  ,  1  ,  1  ,  1  /)
        cspe2=(/  0  , 22  , 24  , 82  , 25  , 32  /)
        cspe3=(/  6  ,  8  ,  8  ,  8  ,  8  ,  8  /)
      endif
!
      if ( ntoz .gt. 0 ) then
        ntrchk = ntoz - 1
      else
        ntrchk = ncld
      endif
!
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,dpd,dpdb,pkout,dew,tens)

      if(ntrac.le.ntrchk)then ! all hydrometeors

      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,cspec(ntrac)
      end do
      write( lrec(lpout), '(a3,a3)' ) 'h00',cspec(ntrac)
      Ptp0=cspe0(ntrac) ;Ptp1=cspe1(ntrac)
      Ptp2=cspe2(ntrac) ;Ptp3=cspe3(ntrac)
!
      else if(ntrac.eq.ntoz)then
!
        do k = 1, lpout-1
          lpl = int(plev(k)+0.001)
          write( lrec(k), '(i3.3,a3)' ) lpl,'560'   ! ozone
        end do
        lrec(lpout) = 'h00560'
        Ptp0=0 ;Ptp1=14 ;Ptp2=1 ;Ptp3=8 !grib code
!
      else if(ntrac.eq.ncld+1)then
!
        do k = 1, lpout-1
          lpl = int(plev(k)+0.001)
          write( lrec(k), '(i3.3,a3)' ) lpl,'550'   ! combine all condensates together
        end do
        lrec(lpout) = 'h00550'
        Ptp0=0 ;Ptp1=1 ;Ptp2=235 ;Ptp3=9 !grib code 
      else
        goto 40
      endif
!
!
      lenc= nx*my
      ncnt= 0
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!
! relative humidity must be smaller or equal 1.0
!
      call unify_reduceintp_gpu(nx,my,my_max,dew(1,1,k),glob)
!$acc parallel loop collapse(2) async(async_id)
      do j=1,my
        do  i=1,nx
          glob(i,j)= max(glob(i,j),0.0)
        enddo
      enddo
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      gtp0=(/ptp0,ptp1,ptp2,ptp3,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,glob,pout,gtp0,gtp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,gtp1(1),gtp1(2),gtp1(3),gtp1(4),gtp1(5) &
            ,gtp1(6),gtp1(7),pout)
      endif
!
   40 continue
      deallocate ( cspec,cspe0,cspe1,cspe2,cspe3 )
!$acc wait(async_id)
!$acc exit data delete(pout) async(async_id)
!$acc wait(async_id)
      return
  end subroutine shumout2

! cloud fraction output
  subroutine cloudout(nx,my,my_max,lpout,lev,itau,idtg   &
      , plev,num,whtlev,pkout,pk,pklp,clds,cldb,cldfc,glob,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use param, only : ncld
      use mod_grb2_param
      use const ,only:outdms,outgrb2 , RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ntrac,ncnt,ntrchk

      real      pkout(lpout),pklp(nxp,my_max)                        &
      , pk(nxp,lev,my_max),clds(nxp,lev,my_max),cldb(nxp,my_max)      &
      , plev(lpout),whtlev(num)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),tmp(nxp,my_max)       &
      , cldfc(nxp,my_max,lpout),ffx(nx,my_max)
!
      integer   i,k,lpl,n,lenc,istat,jj,j,nxj

      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
      logical :: lwrite
      integer:: ptp0(9),ptp1(9)
!$acc wait(async_id)
!$acc enter data create(pout) async(async_id)
!
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,clds,cldb,pkout,cldfc,tens)

      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'770'
      end do
      lrec(lpout) = 'h00770'
!
      lenc= nx*my
      ncnt= 0
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      call unify_reduceintp_gpu(nx,my,my_max,cldfc,glob)
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,6,32,2,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
   40 continue

!$acc wait(async_id)
!$acc exit data delete(pout) async(async_id)
!$acc wait(async_id)
      return
  end subroutine cloudout


  subroutine surfout(nx,my,my_max,itau,idtg,taudir,ntau,pdiff  &
       ,pt,ptop,slp,ptend,glob,ggdef,lwrite)
!
      use index
      use mpe
      use rank
      use mod_grb2_param
      use const ,only:outdms,outgrb2 ,RTYPE,kflag ,domfc,ihdgo,ihdgo2
!
      implicit  none
      integer   nx,my,my_max,i,j,jj,kk,n,lev,nxj,itau,ntau,num,lenc,istat

      real      pdiff(nxp,my_max)
      real(kind=RTYPE) ptend(nxp,my_max),pt(nxp,my_max),glob(nx,my),   &
                       slp(nxp,my_max),tmp(nxp,my_max)
      character*18 taudir(ntau)
      character*4 ggdef
!
      real      ptop,tnshun
!
      integer*8 idtg
      character*6 label(ntau),labx
   
      logical :: lwrite
!
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
      do jj = 1, jlistnum
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj) slp(i,jj)= ( pt(i,jj)+pdiff(i,jj) )
        enddo
      enddo

!
      num= 0
      do 20 n=1,ntau
        read(taudir(n),'(a6,1x,i4)') labx,lev
        if(lev.eq.0) then
          num= num+1
          label(num)= labx
        endif
   20 continue

      if(num.eq.0) return
!
      tnshun= 1.0
      lenc= nx*my
!
!$acc wait(async_id)
!$acc enter data create(tmp) async(async_id)
      do 100 kk=1,num
!
!  sea surface level pressure
!

      if(label(kk).eq.'SSL010' .or. label(kk).eq.'ssl010') then
!$acc wait(async_id)
!$acc parallel loop collapse(2)  private( j,nxj ) async(async_id)
        do jj = 1, jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj) tmp(i,jj)= slp(i,jj) * 100.0 !hPa -> Pa
          enddo
        enddo
        call unify_reduceintp_gpu(nx,my,my_max,tmp,glob)
        call syslbl_w('ssl010',idtg,itau,ggdef)
        if( itau==0 .or. itau .gt. nint(domfc) )then
         call qmaxn3_w(glob,1,1,1,nx,my,1)
         if(outgrb2==1.and.myrank==0)then
           ihdgo2 = ihdgo
           call wrtgrb2_v2_gpu(itau,0,3,1,1,101,0,0,glob)
         endif
         if(outdms.gt.0)then
          !$acc update self(glob) async(async_id)
          !$acc wait(async_id)
           glob=glob/100.0
           if(lwrite) call dmswrit(nx,my,lenc,kflag,glob,istat)
         endif
        endif !itau .gt. domfc
!
!  terrain pressure
!
      else if(label(kk).eq.'B00010' .or. label(kk).eq.'b00010') then
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
        do jj = 1, jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj)  tmp(i,jj) = ( pt(i,jj) + ptop ) * 100.0
          enddo
        enddo
        call unify_reduceintp_gpu(nx,my,my_max,tmp,glob)
        call syslbl_w('b00010',idtg,itau,ggdef)
        call qmaxn3_w(glob,1,1,1,nx,my,1)
        if(outgrb2==1.and.myrank==0)then
          ihdgo2 = ihdgo
          call wrtgrb2_v2_gpu(itau,0,3,0,1,103,0,0,glob)
        endif
        if(outdms.gt.0)then
          !$acc update self(glob) async(async_id)
          !$acc wait(async_id)
          glob=glob/100.0
          if(lwrite) call dmswrit(nx,my,lenc,kflag,glob,istat)
        endif

!
!  terrain pressure tendency
!
!      else if(label(kk).eq.'B00011' .or. label(kk).eq.'b00011') then
!        $acc parallel loop collapse(2) private( j,nxj ) async(async_id)
!        do jj = 1, jlistnum
!          do i=1,nxp
!          j=jlist1(jj)
!          nxj=nxdef_2d(j)
!           if(i<=nxj)then
!            tmp(i,jj)= ptend(i,jj)*3600.0
!           endif
!          enddo
!        enddo
!        call unify_reduceintp_gpu(nx,my,my_max,tmp,glob)
!        call syslbl_w('b00011',idtg,itau,ggdef)
!        call qmaxn3_w(glob,1,1,1,nx,my,1)
!        if(lwrite) call dmswrit(nx,my,lenc,kflag,glob,istat)
!
      endif
!
  100 continue
!$acc wait(async_id)
!$acc exit data delete(tmp) async(async_id)
!$acc wait(async_id)
      return
  end subroutine surfout

  subroutine tempout(nx,my,my_max,lpout,lev,itau,idtg  &
      , plev,num,whtlev,pkout,pk,pklp,tt,ttbot,temp,ggdef,lwrite)
!
      use index
      use rank,   only : myrank
      use mod_grb2_param
      use const ,only:outdms,outgrb2 ,RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,k,i,j,n,istat
      integer   num,lpl,lenc,ncnt
      integer:: ptp0(9),ptp1(9)
      real      tnshun

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max)        &
      , tt(nxp,lev,my_max),ttbot(nxp,my_max),plev(lpout),whtlev(num)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),slp(nx,my),temp(nxp,my_max,lpout)

!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!
     logical :: lwrite

!$acc wait(async_id)
!$acc enter data create(glob,slp,pout) async(async_id)
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,tt,ttbot,pkout,temp,tens)
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'100'
      end do
      lrec(lpout) = 'h00100'
!
!
      lenc= nx*my
!
      tnshun= 1.0
!
      ncnt=0
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then  
!
!
      call unify_reduceintp_gpu(nx,my,my_max,temp(1,1,k),slp)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
      call smth9(nx,my,slp,glob,1)
      call smth9(nx,my,glob,slp,2)
!
      call qmaxn3_w(slp,1,1,1,nx,my,1)
      ptp0=(/0,0,0,2,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,slp,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!$acc wait(async_id)
!$acc exit data delete(glob,slp,pout) async(async_id)
!$acc wait(async_id)
!
      return
  end subroutine tempout

  subroutine vortout(nx,my,my_max,lpout,lev,itau,idtg     &
      , plev,num,whtlev,pkout,pk,pklp,rvor,rvorb,vor,ggdef,v850,v700,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:outdms,outgrb2 , RTYPE,kflag

      implicit  none

      integer   nx,my,my_max,lpout,lev,itau,num,ncnt

      real      pkout(lpout),pklp(nxp,my_max)                        &
      , pk(nxp,lev,my_max),rvorb(nxp,my_max)  &
      , plev(lpout),whtlev(num)
      real rvor(nxp,lev,my_max)
      real(kind=RTYPE) vor(nxp,my_max,lpout)
      real(kind=RTYPE) v850(nxp,my_max),v700(nxp,my_max)
      real(kind=RTYPE) wk1(nx,my),pout(nx,my)

!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef

      integer   k,lpl,lenc,i,n,istat,jj,j,nxj
      integer:: ptp0(9),ptp1(9)
!
      logical :: lwrite
!$acc wait(async_id)
!$acc enter data create(wk1,pout) async(async_id)
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,rvor,rvorb,pkout,vor,tens)
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'240'
      end do
      lrec(lpout) = 'h00240'
!
!
      lenc= nx*my
      ncnt= 0
!
      do k=1,lpout
        if(plev(k).eq.850.)then
!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
        do jj=1, jlistnum
          do i=1,nxp
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          if(i<=nxj) v850(i,jj)= vor(i,jj,k)
          enddo
        enddo
!
!  reduceintp has been done in voterp (2011/5)
!
        else if(plev(k).eq.700.)then
!$acc parallel loop collapse(2)  private(j,nxj) async(async_id)
        do jj=1, jlistnum
          do i=1,nxp
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          if(i<=nxj) v700(i,jj)= vor(i,jj,k)
          enddo
        enddo
!
!  reduceintp has been done in voterp (2011/5)
!
        endif
      enddo
!
      do 20 n=1,num
      do 10 k=1,lpout
      if(plev(k).eq.whtlev(n)) then
      call unify_reduceintp_gpu(nx,my,my_max,vor(1,1,k),wk1)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
      call qmaxn3_w(wk1,1,1,1,nx,my,1)
      ptp0=(/0,2,12,6,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
      go to 20
      endif
   10 continue
   20 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!$acc wait(async_id)
!$acc exit data delete(wk1,pout) async(async_id)
!$acc wait(async_id)
!
      return
  end subroutine vortout

  subroutine windout(nx,my,my_max,lpout,lev,itau,idtg      &
      , plev,num,whtlev,cosl,pkout,pk,pklp,ut,vt,sdhat,work3d,utb,vtb     &
      , wind,glob,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:outdms,outgrb2 , RTYPE,kflag
      use openacc

      implicit none

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max)          &
      , rdiv(nxp,lev,my_max),work3d(nxp,lev,my_max)                       &
      , plev(lpout),whtlev(num)
      real(kind=RTYPE) ut(nxp,lev,my_max),vt(nxp,lev,my_max)              &
      , sdhat(nxp,lev,my_max),cosl(my),wind(nxp,my_max,lpout)
      real      utb(nxp,my_max),vtb(nxp,my_max),wtb(nxp,my_max)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),tmp(nxp,my_max)
!
      integer   nx,my,my_max,lpout,lev,itau,jj,nxj,ncnt
      integer   num,k,lenc,lpl,n,i,j,istat
      integer:: ptp0(9),ptp1(9)

      real      rad,xxx
      integer, parameter:: async_id = 1

      integer*8 idtg
      character*6 lrec(lpout),krec(lpout),mrec(lpout)
      character*4 ggdef
!
      data rad/6.371e6/
!
     logical :: lwrite

!$acc wait(async_id)
!$acc enter data create(wtb,tmp,pout) async(async_id)
      lenc= nx*my
      ncnt= 0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'200'
       write( krec(k), '(i3.3,a3)' ) lpl,'210'
       write( mrec(k), '(i3.3,a3)' ) lpl,'220'
      end do
      lrec(lpout) = 'h00200'
      krec(lpout) = 'h00210'
      mrec(lpout) = 'h00220'
!
! first: do the u components
!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
     do jj = 1, jlistnum
       do k = 1, lev
        do i = 1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj) work3d(i,k,jj)=ut(i,k,jj)
        enddo
       enddo
     enddo
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,work3d,utb,pkout,wind,tens)
!

      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
!$acc parallel loop collapse(2) private(j,nxj,xxx) async(async_id)
        do jj = 1, jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj)then
            xxx= rad/cosl(j)
            if(pkout(k).gt.pklp(i,jj)) wind(i,jj,k)= utb(i,jj)
            tmp(i,jj)=wind(i,jj,k)*xxx
           endif
          enddo
        enddo
      call unify_reduceintp_gpu(nx,my,my_max,tmp,glob)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,2,2,2,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!----------------------------------------
!
!  now the v components
!
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
     do jj = 1, jlistnum
       do k = 1, lev
        do i = 1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj) work3d(i,k,jj)=vt(i,k,jj)
        enddo
       enddo
     enddo
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,work3d,vtb,pkout,wind,tens)
!
      ncnt= 0
      do 40 n=1,num
      do 20 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
!$acc parallel loop collapse(2) private(j,nxj,xxx) async(async_id)
        do jj = 1, jlistnum
          do i=1,nxp
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if(i<=nxj)then
            xxx= rad/cosl(j)
            if(pkout(k).gt.pklp(i,jj)) wind(i,jj,k)= vtb(i,jj)
            tmp(i,jj)=wind(i,jj,k)*xxx
           endif
          enddo
        enddo
      call unify_reduceintp_gpu(nx,my,my_max,tmp,glob)
!
      call syslbl_w(krec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,2,3,2,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 40
      endif
   20 continue
   40 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!----------------------------------------
!
!  now the w components
!
!$acc parallel loop collapse(2) private( j,nxj ) async(async_id)
     do jj = 1, jlistnum
       do i = 1,nxp
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       if(i<=nxj) wtb(i,jj)=0.0
       enddo
     enddo
!$acc parallel loop collapse(3) private(j,nxj) async(async_id)
     do jj = 1, jlistnum
       do k = 1, lev
        do i = 1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj) work3d(i,k,jj)=sdhat(i,k,jj) * 100.0
        enddo
       enddo
     enddo
      call voterp_gpu(nx,my,my_max,lev,lpout,pk,pklp,work3d,wtb,pkout,wind,tens)
!
      ncnt= 0
      do 42 n=1,num
      do 22 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      call unify_reduceintp_gpu(nx,my,my_max,wind(1,1,k),glob)
      call syslbl_w(mrec(k),idtg,itau,ggdef)
!
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,2,8,6,100,-2,nint(plev(k)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 42
      endif
   22 continue
   42 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) then
         !$acc update self(pout) async(async_id)
         !$acc wait(async_id)
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!$acc wait(async_id)
!$acc exit data delete(wtb,tmp,pout) async(async_id)
!$acc wait(async_id)
!
      return
  end subroutine windout
!
end module mod_outflds_gpu
