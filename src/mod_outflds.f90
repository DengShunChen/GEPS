module mod_outflds
 implicit none
 
 
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
                work3d(nxp,lev,my_max),rdivb(nxp,my_max),               &
                plev(lpout),whtlev(num)
      real(kind=RTYPE) rdiv(nxp,lev,my_max),div(nxp,my_max,lpout)
      real      tens(lev+1)
      real(kind=RTYPE) wk1(nx,my),pout(nx,my)

!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!
      integer    k,lpl,lenc,n,num,istat
      integer:: ptp0(9),ptp1(9)
  
      logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'230'
      end do
      lrec(lpout) = 'h00230'
      work3d=rdiv
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,work3d,rdivb,pkout,div,tens)
!
      lenc= nx*my
      ncnt= 0
!
      do 20 n=1,num
      do 10 k=1,lpout
      if(plev(k).eq.whtlev(n)) then
      call unify_reduceintp(nx,my,my_max,div(1,1,k),wk1)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,wk1,istat)
      call qmaxn3_w(wk1,1,1,1,nx,my,1)
      ptp0=(/0,2,11,6,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
      go to 20
      endif
   10 continue
   20 continue

      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
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
      real      tens(lev+1)
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

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'290'
      end do
      lrec(lpout) = 'h00290'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,rdrag,rdragb,pkout,drag,tens)
!
      lenc= nx*my
      ncnt= 0
!
!!      do jj =1, jlistnum
!!      j=jlist1(jj)
!!      nxj=nxdef_2d(j)
!!      do i=1,nxj
!!        pk_lev_m1(i,j)=pk(i,lev-1,jj)
!!      enddo
!!      enddo
!
!!      call mpe_unify(pk_lev_m1,nx,my,2,mpe_double)
!
      do 20 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      do 45 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 45 i=1,nxj
!      ij=(j-1)*nx+i
      if(pkout(k).gt.pk(i,lev-1,jj)) drag(i,jj,k)= 0.
   45 continue
!
      call unify_reduceintp(nx,my,my_max,drag(1,1,k),wk1)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,wk1,istat)
      call qmaxn3_w(wk1,1,1,1,nx,my,1)
      ptp0=(/0,2,196,6,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
      go to 20
      endif
   10 continue
   20 continue

      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
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

      real      tens(lev+1),phistd(lpout)
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

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'000'
      end do
      lrec(lpout) = 'h00000'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,phi,phib,pkout,phips,tens)
!
      lenc= nx*my
      ncnt= 0
!
      do k=1,lpout
        if(plev(k).eq.850.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            h850(i,jj)= phips(i,jj,k)
          enddo
        enddo
!byl          call unify_reduceintp(nx,my,my_max,phips(1,1,k),h850)
!!          do i=1,lenc
!!            h850(i,1)=phips(i,k)
!!          enddo
!
!          if(itau.le.72 .and. lreduce.eq.1)then
!byl            call smth9(nx,my,h850,glob,1)
!byl            h850=glob
!byl            call smth9(nx,my,h850,glob,2)
!byl            h850=glob
!          endif
        else if(plev(k).eq.500.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            h500(i,jj)= phips(i,jj,k)
          enddo
        enddo
!byl          call unify_reduceintp(nx,my,my_max,phips(1,1,k),h500)
!!          do i=1,lenc
!!            h500(i,1)=phips(i,k)
!!          enddo
!          if(itau.le.72 .and. lreduce.eq.1)then
!byl            call smth9(nx,my,h500,glob,1)
!byl            h500=glob
!byl            call smth9(nx,my,h500,glob,2)
!byl            h500=glob
!          endif
        endif
      enddo
!
      do 20 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      do 11 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 11 i=1,nxj
       tmp(i,jj)= phips(i,jj,k)+phistd(k)
   11 continue
      call unify_reduceintp(nx,my,my_max,tmp,glob)
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!      if(lreduce.eq.1 .and. itau.le.72)then
        call smth9(nx,my,glob,slp,1)
        glob=slp
        call smth9(nx,my,glob,slp,2)
!        glob=slp
!      endif
!
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,slp,istat)
      call qmaxn3_w(slp,1,1,1,nx,my,1)
      ptp0=(/0,3,5,1,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,slp,pout,ptp0,ptp1)
!
      go to 20
      endif
   10 continue
   20 continue

      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
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
      , plev(lpout),whtlev(num),tens(lev+1)
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

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'510'
      end do
      lrec(lpout) = 'h00510'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,dpd,dpdb,pkout,dew,tens)
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
!!      do 20 jj=1, jlistnum
!!      j=jlist1(jj)
!!      nxj=nxdef_2d(j)
!!      do 20 i=1,nxj
!!       tmp(i,jj)= min(100.,max(dew(i,jj,k)*100.,0.0))
!!   20 continue
      call mpe2d_unify_nx(ffx,dew(1,1,k))
      if( lreduce.eq.1 ) then
        do jj =1, jlistnum
          j=jlist1(jj)
          call reduceintp(ffx(1,jj),nxdef(j),nx,1)
          do i=1,nx
             ffx(i,jj)=min(100.,max(ffx(i,jj)*100.,0.0))
          enddo
        enddo
      endif
      call mpe2d_unify_my(glob,ffx)
!      call unify_reduceintp(nx,my,my_max,dew(1,1,k),glob)
!      do 20 i=1,lenc
!      glob(i,1)= min(100.,max(glob(i,1)*100.,0.0))
!   20 continue
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,1,1,2,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
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
      , plev(lpout),whtlev(num),tens(lev+1)
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
      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
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
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,dpd,dpdb,pkout,dew,tens)
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
!      do 20 jj=1, jlistnum
!      j=jlist1(jj)
!      nxj=nxdef_2d(j)
!      do 20 i=1,nxj
!       tmp(i,jj)= max(dew(i,jj,k),0.0)
!   20 continue
      call mpe2d_unify_nx(ffx,dew(1,1,k))
      if( lreduce.eq.1 ) then
        do jj =1, jlistnum
          j=jlist1(jj)
          call reduceintp(ffx(1,jj),nxdef(j),nx,1)
          do i=1,nx
             ffx(i,jj)=max(ffx(i,jj),qmin)
          enddo
        enddo
      endif
      call mpe2d_unify_my(glob,ffx)
!      do 20 i=1,lenc
!        glob(i,1)= max(glob(i,1),0.0)
!   20 continue
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      gtp0=(/ptp0,ptp1,ptp2,ptp3,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,glob,pout,gtp0,gtp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,gtp1(1),gtp1(2),gtp1(3),gtp1(4),gtp1(5) &
            ,gtp1(6),gtp1(7),pout)
      endif
!
   40 continue
      deallocate ( cspec,cspe0,cspe1,cspe2,cspe3 )
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
      , plev(lpout),whtlev(num),tens(lev+1) 
      real(kind=RTYPE) pout(nx,my),glob(nx,my),tmp(nxp,my_max)       &
      , cldfc(nxp,my_max,lpout),ffx(nx,my_max)
!
      integer   i,k,lpl,n,lenc,istat,jj,j,nxj

      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
      logical :: lwrite
      integer:: ptp0(9),ptp1(9)
!
      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'770'
      end do
      lrec(lpout) = 'h00770'
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,clds,cldb,pkout,cldfc,tens)
!
      lenc= nx*my
      ncnt= 0
!
      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      do 11 jj=1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 11 i=1,nxj
       tmp(i,jj)= cldfc(i,jj,k) * 100.0
   11 continue
      call unify_reduceintp(nx,my,my_max,tmp,glob)
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,6,32,2,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
   40 continue

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
      character*16 taudir(ntau)
      character*4 ggdef
!
      real      ptop,tnshun
!
      integer*8 idtg
      character*6 label(ntau),labx
   
      logical :: lwrite
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
!byl          slp(i,j)= pt(i,jj)+pdiff(i,jj)
          slp(i,jj)= pt(i,jj)+pdiff(i,jj)
          ptend(i,jj)= ptend(i,jj)*3600.0
        enddo
      enddo
!
!byl      call unify_reduceintp(nx,my,my_max,tmp,slp)
!byl      call mpe_unify(slp,nx,my,2,mpe_double)
!byl      if( lreduce.eq.1 ) call reduceintp (slp,nxdef,nx,my)
!byl      call smth9(nx,my,slp,glob,1)
!byl      slp=glob
!byl      call smth9(nx,my,slp,glob,2)
!byl      slp=glob
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
      do 100 kk=1,num
!
!  sea surface level pressure
!

      if(label(kk).eq.'SSL010' .or. label(kk).eq.'ssl010') then
        call unify_reduceintp(nx,my,my_max,slp,glob)
        call syslbl_w('ssl010',idtg,itau,ggdef)
        if( itau==0 .or. itau .gt. nint(domfc) )then
        if(outdms.gt.0)then
          if(lwrite) call dmswrit(nx,my,lenc,kflag,glob,istat)
        endif
        if(outgrb2==1.and.myrank==0)then
          ihdgo2 = ihdgo
          glob=glob*100.0
          call wrt_grb2_v2(itau,0,3,1,1,101,0,0,glob)
        endif
        call qmaxn3_w(glob,1,1,1,nx,my,1)
        endif !itau .gt. domfc
!
!  terrain pressure
!
      else if(label(kk).eq.'B00010' .or. label(kk).eq.'b00010') then
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
!byl            glob(i,j) = pt(i,jj) + ptop
            tmp(i,jj) = pt(i,jj) + ptop
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,tmp,glob)
!     if(myrank.eq.0) then
!       open(30,file='sfc_pres.bin',status='unknown', &
!           form='unformatted',access='direct',recl=262656)
!       write(30,rec=1) ((glob(i,j),i=432,647),j=410,561)
!       close(30)
!     endif
!byl        call mpe_unify(glob,nx,my,2,mpe_double)
        call syslbl_w('b00010',idtg,itau,ggdef)
!byl        if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
        if(outdms.gt.0)then
          if(lwrite) call dmswrit(nx,my,lenc,kflag,glob,istat)
        endif
        if(outgrb2==1.and.myrank==0)then
          ihdgo2 = ihdgo
          glob=glob*100.0
          call wrt_grb2_v2(itau,0,3,0,1,103,0,0,glob)
        endif
        call qmaxn3_w(glob,1,1,1,nx,my,1)

!
!  terrain pressure tendency
!
!      else if(label(kk).eq.'B00011' .or. label(kk).eq.'b00011') then
!        do i = 1, lenc
!          glob(i,1) = ptend(i,1)
!        end do
!        call syslbl('b00011',idtg,itau,ggdef,lrec)
!        if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!        if(lwrite) call dmswrit(nx,my,lrec,lenc,kflag,ifilout,glob,istat)
!        call qmaxn3(glob,lrec(1:14),lrec(15:26),1,1,1,nx,my,1)
!
      endif
!
  100 continue
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

      integer   nx,my,my_max,lpout,lev,itau,k,i,n,istat
      integer   num,lpl,lenc,ncnt
      integer:: ptp0(9),ptp1(9)
      real      tnshun

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max)        &
      , tt(nxp,lev,my_max),ttbot(nxp,my_max),plev(lpout),whtlev(num)
      real      tens(lev+1)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),slp(nx,my),temp(nxp,my_max,lpout)

!
      integer*8 idtg
      character*6 lrec(lpout)
      character*4 ggdef
!
     logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev+1)= 0.0
      tens(lev)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'100'
      end do
      lrec(lpout) = 'h00100'
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,tt,ttbot,pkout,temp,tens)
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
      call unify_reduceintp(nx,my,my_max,temp(1,1,k),glob)
!!      do 11 i=1,lenc
!!      glob(i,1)= temp(i,k)
!!   11 continue
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
        call smth9(nx,my,glob,slp,1)
        glob=slp
        call smth9(nx,my,glob,slp,2)
!!        glob=slp
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(slp,1,1,1,nx,my,1)
      ptp0=(/0,0,0,2,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,slp,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
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
      , pk(nxp,lev,my_max),work3d(nxp,lev,my_max),rvorb(nxp,my_max)  &
      , plev(lpout),whtlev(num)
      real(kind=RTYPE) rvor(nxp,lev,my_max),vor(nxp,my_max,lpout)
      real      tens(lev+1)
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

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
!
      do k = 1, lpout-1
       lpl = int(plev(k)+0.001)
       write( lrec(k), '(i3.3,a3)' ) lpl,'240'
      end do
      lrec(lpout) = 'h00240'
      work3d=rvor
!
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,work3d,rvorb,pkout,vor,tens)
!
      lenc= nx*my
      ncnt= 0
!
      do k=1,lpout
        if(plev(k).eq.850.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            v850(i,jj)= vor(i,jj,k)
          enddo
        enddo
!byl        call unify_reduceintp(nx,my,my_max,vor(1,1,k),v850)
!!          do i=1,lenc
!!            v850(i,1)=vor(i,k)
!!          enddo
!
!  reduceintp has been done in voterp (2011/5)
!
        else if(plev(k).eq.700.)then
        do jj=1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            v700(i,jj)= vor(i,jj,k)
          enddo
        enddo
!byl        call unify_reduceintp(nx,my,my_max,vor(1,1,k),v700)
!!          do i=1,lenc
!!            v700(i,1)=vor(i,k)
!!          enddo
!
!  reduceintp has been done in voterp (2011/5)
!
        endif
      enddo
!
      do 20 n=1,num
      do 10 k=1,lpout
      if(plev(k).eq.whtlev(n)) then
      call unify_reduceintp(nx,my,my_max,vor(1,1,k),wk1)
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,wk1,istat)
      call qmaxn3_w(wk1,1,1,1,nx,my,1)
      ptp0=(/0,2,12,6,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,wk1,pout,ptp0,ptp1)
      go to 20
      endif
   10 continue
   20 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) & 
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
      return
  end subroutine vortout

  subroutine windout(nx,my,my_max,lpout,lev,itau,idtg      &
      , plev,num,whtlev,cosl,pkout,pk,pklp,ut,vt,sdhat,utb,vtb     &
      , wind,glob,ggdef,lwrite)
!
      use index
      use rank, only : myrank
      use mod_grb2_param
      use const ,only:outdms,outgrb2 , RTYPE,kflag

      implicit none

      real      pkout(lpout),pklp(nxp,my_max),pk(nxp,lev,my_max)          &
      , rdiv(nxp,lev,my_max),work3d(nxp,lev,my_max)                       &
      , utb(nxp,my_max),vtb(nxp,my_max),plev(lpout),whtlev(num)
      real(kind=RTYPE) ut(nxp,lev,my_max),vt(nxp,lev,my_max)              &
      , sdhat(nxp,lev,my_max),cosl(my),wind(nxp,my_max,lpout)

      real      tens(lev+1),wtb(nxp,my_max)
      real(kind=RTYPE) pout(nx,my),glob(nx,my),tmp(nxp,my_max)
!
      integer   nx,my,my_max,lpout,lev,itau,jj,nxj,ncnt
      integer   num,k,lenc,lpl,n,i,j,istat
      integer:: ptp0(9),ptp1(9)

      real      rad,xxx
!lzl +add
!!      real      pklzl(nx,my)
!lzl -end

      integer*8 idtg
      character*6 lrec(lpout),krec(lpout),mrec(lpout)
      character*4 ggdef
!
      data rad/6.371e6/
!
     logical :: lwrite

      do k = 1, lev+1
       tens(k) = 1.0
      end do
      tens(lev)= 0.0
      tens(lev+1)= 0.0
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
      work3d=ut
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,work3d,utb,pkout,wind,tens)
!
!!      if( lreduce.eq.1 ) call reduceintp (utb,nxdef,nx,my)
!
!lzl +add======================================================
!!      do i =1,nx
!!      do j =1,my
!!         pklzl(i,j)=pklp(i,j)
!!      end do
!!      end do
!
!!      if( lreduce.eq.1 ) then
!!          call reduceintp(pklzl,nxdef,nx,my)
!!      endif

!lzl -end=====================================================          

      do 30 n=1,num
      do 10 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xxx= rad/cosl(j)
          do i=1,nxj
           if(pkout(k).gt.pklp(i,jj)) wind(i,jj,k)= utb(i,jj)
           tmp(i,jj)=wind(i,jj,k)*xxx
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,tmp,glob)
!!      do 45 i=1,nx*my
!lzl c      if(pkout(k).gt.pklp(i,1)) wind(i,1,k)= utb(i,1)
!!      if(pkout(k).gt.pklzl(i,1)) wind(i,1,k)= utb(i,1)  !lzl use full grid(pklzl)
!!   45 continue
!
!
!  deweight wind with cos latitude, earth radius
!
!!      do 50 j=1,my
!!      xxx= rad/cosl(j)
!!      do 50 i=1,nx
!!      glob(i,j)= wind(i,j,k)*xxx
!!   50 continue
!
      call syslbl_w(lrec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,2,2,2,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 30
      endif
   10 continue
   30 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
!  now the v components
!
      work3d=vt
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,work3d,vtb,pkout,wind,tens)
!
!!      if( lreduce.eq.1 ) call reduceintp (vtb,nxdef,nx,my)
!
!lzl +add===============================================================
!!      do i =1,nx
!!      do j =1,my
!!         pklzl(i,j)=pklp(i,j)
!!      end do
!!      end do
!
!!      if( lreduce.eq.1 ) then
!!          call reduceintp(pklzl,nxdef,nx,my)
!!      endif

!lzl -end================================================================
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
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xxx= rad/cosl(j)
          do i=1,nxj
           if(pkout(k).gt.pklp(i,jj)) wind(i,jj,k)= vtb(i,jj)
           tmp(i,jj)=wind(i,jj,k)*xxx
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,tmp,glob)
!
!!      do 55 i=1,nx*my
!lzl c      if(pkout(k).gt.pklp(i,1)) wind(i,1,k)= vtb(i,1)
!!      if(pkout(k).gt.pklzl(i,1)) wind(i,1,k)= vtb(i,1) !lzl use full grid (pklzl)
!!   55 continue
!
!  deweight wind with cos latitude, earth radius
!
!!      do 60 j=1,my
!!      xxx= rad/cosl(j)
!!      do 60 i=1,nx
!!      glob(i,j)= wind(i,j,k)*xxx
!!   60 continue
!
      call syslbl_w(krec(k),idtg,itau,ggdef)
!
!  reduceintp has been done in voterp (2011/5)
!
!!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      ptp0=(/0,2,3,2,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 40
      endif
   20 continue
   40 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) & 
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
!  now the w components
!
      wtb=0.
      work3d=sdhat
      call voterp(nx,my,my_max,lev,lpout,pk,pklp,work3d,wtb,pkout,wind,tens)
!
      ncnt= 0
      do 42 n=1,num
      do 22 k=1,lpout
!
      if(plev(k).eq.whtlev(n)) then
!
      call unify_reduceintp(nx,my,my_max,wind(1,1,k),glob)
!!      do j=1,my
!!      do i=1,nx
!!        glob(i,j)= wind(i,j,k)
!!      enddo
!!      enddo
!
      call syslbl_w(mrec(k),idtg,itau,ggdef)
!
!      if(lwrite) call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
      call qmaxn3_w(glob,1,1,1,nx,my,1)
      glob=glob*100.
      ptp0=(/0,2,8,6,100,-2,nint(plev(k)),-999,-999/)
      call split2(nx,my,lenc,ncnt,glob,pout,ptp0,ptp1)
      go to 42
      endif
   22 continue
   42 continue
      if(outdms.gt.0)then
      if(lwrite .and. myrank .lt. ncnt ) &
        call dmswrit_split(nx,my,lenc,kflag,pout,istat)
      endif
      if(outgrb2 == 1 )then
      if(lwrite .and. myrank .lt. ncnt ) &
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),pout)
      endif
!
      return
  end subroutine windout
!
end module mod_outflds
