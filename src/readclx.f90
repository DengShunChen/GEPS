      subroutine readclx(nx,my,my_max,julian,land,sea,ice,solt,wet,z0, &
                         alb,sst,bckfile,sigmaf,istyp,ivegtyp,ls       &
                        ,shdmax,shdmin,slopetyp,snoalb,ggdef,isot,ivegsrc)
!
!  read climt data from data base
!
!----------------------------------------------------------------
!  input :
!         nx     : dimension of e-w direction
!         my     : dimension of s-n direction
!         julian : julian day
!  output :
!         land   : index of bare soil or not (logical)
!         sea    : index of open water or not (logical)
!         ice    : index of ice cover or not (logical)
!         solt   : soil temperature ( i.e. ground temp.) (k)
!         wet    : ground wetness
!         z0     : surface roughness  (m)
!         alb    : surface albedo
!         sst    : sea surface temperature (k)
!         sigmaf : vegetation fraction
!         istyp  : soil type(1-9)
!         ivegtyp: vegetation type(1-13)
!-----------------------------------------------------------------
!
      use index
!
      implicit  none

      integer   nx,my,my_max,julian,isot,ivegsrc

      real      solt(nxp,my_max),wet(nxp,my_max),z0(nxp,my_max),alb(nxp,my_max),  &
                sst(nxp,my_max)                           &
!soil
               ,sigmaf(nxp,my_max)                        &
!noah
               ,shdmax(nxp,my_max),shdmin(nxp,my_max) &
               ,snoalb(nxp,my_max)

      integer slopetyp(nxp,my_max)

      integer istyp(nxp,my_max),ivegtyp(nxp,my_max)

!
      logical land(nxp,my_max),sea(nxp,my_max),ice(nxp,my_max)
!
!  working array as climt data base
!
      real      sstcl(nx,my,2),soltcl(nx,my,2),wetcl(nx,my,2),      &
                albcl(nx,my,2),z0cl(nx,my,2),                       &
!soil
                vfrcl(nx,my,2)
!soil
      integer   ls(nxp,my_max),icex(nxp,my_max),iglob(nx,my)

      character bckfile*80,lrec*26,blnk*1,ggdef*4
      integer   mon(12),mondy(13)
      data mon/15,46,74,105,135,166,196,227,258,288,319,349/
      data mondy/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data blnk/' '/

      integer   i,j,k,jul,nxj,mm,istat,monidex,lncrec,jj,ii
      real      coef1,coef2
!
      lncrec=nx*my
!
      jul=julian
      if(jul .ge. 366)jul=365
!
!  what month is it ?
!
      do k=1,12
      if(jul .gt. mondy(k) .and. jul .le. mondy(k+1))monidex=k
      end do
!
  11  format('W00100','gbck',a4,4x,i2.2,6x)  ! W10
  12  format('S9M100','gbck',a4,4x,i2.2,6x)  ! X10
  14  format('S00030','gbck',a4,4x,i2.2,6x)  ! S35
  15  format('S00040','gbck',a4,4x,i2.2,6x)  ! S44
  16  format('S00070','gbck',a4,11x,a1) ! S07
  17  format('S00090','gbck',a4,4x,i2.2,6x) ! S09
!soil
  13  format('S9M5B0','gbck',a4,4x,i2.2,6x)  ! for new soil
  18  format('S000ST','gbck',a4,11x,a1)      ! for new soil
  19  format('S000VF','gbck',a4,4x,i2.2,6x)  ! for new soil
  21  format('S000VT','gbck',a4,11x,a1)      ! for new soil
  22  format('S9Y100','gbck',a4,11x,a1)      ! for new soil
!soil
  23  format('S000VX','gbck',a4,11x,a1)      ! for new soil
  24  format('S000VN','gbck',a4,11x,a1)      ! for new soil
  25  format('S00063','gbck',a4,11x,a1)      ! for new soil
  26  format('S0003X','gbck',a4,11x,a1)      ! for new soil
!source 2(19-soil, 20-veg)
  28  format('S00XST','gbck',a4,11x,a1)      ! for new soil
  31  format('S00XVT','gbck',a4,11x,a1)      ! for new soil
!
!  to interpolat linearly based on julian day
!
!-- climat dataset has 12 months

      if(jul .le. mon(1))jul=jul+365
      if(jul .gt. mon(12))then
!                                -- read dmsfile --
      mm=12
      write(lrec,11)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,sstcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (sstcl(1,1,1),nxdef,nx,my)

!ch?  write(lrec,12)ggdef,mm
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,1),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,1),nxdef,nx,my)

      write(lrec,14)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,albcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (albcl(1,1,1),nxdef,nx,my)

      write(lrec,15)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,z0cl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (z0cl(1,1,1),nxdef,nx,my)

      write(lrec,13)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,wetcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (wetcl(1,1,1),nxdef,nx,my)

!soil
      write(lrec,19)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,vfrcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (vfrcl(1,1,1),nxdef,nx,my)
!soil
      mm=1
      write(lrec,11)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,sstcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (sstcl(1,1,2),nxdef,nx,my)

!ch?  write(lrec,12)ggdef,mm
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,2),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,2),nxdef,nx,my)

      write(lrec,14)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,albcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (albcl(1,1,2),nxdef,nx,my)

      write(lrec,15)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,z0cl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (z0cl(1,1,2),nxdef,nx,my)

      write(lrec,13)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,wetcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (wetcl(1,1,2),nxdef,nx,my)
!soil
      write(lrec,19)ggdef,mm
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,vfrcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (vfrcl(1,1,2),nxdef,nx,my)
!soil
!
      coef1=float(jul-mon(12))/float(380-mon(12))
      coef2=1.-coef1
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            sst(i,jj) =coef1*sstcl(ii,j,2) +coef2*sstcl(ii,j,1)
!ch?        solt(i,jj)=coef1*soltcl(ii,j,2)+coef2*soltcl(ii,j,1)
            alb(i,jj) =coef1*albcl(ii,j,2) +coef2*albcl(ii,j,1)
            z0(i,jj)=coef1*z0cl(ii,j,2)+coef2*z0cl(ii,j,1)
            wet(i,jj)=coef1*wetcl(ii,j,2)+coef2*wetcl(ii,j,1)
!soil
          sigmaf(i,jj)=coef1*vfrcl(ii,j,2)+coef2*vfrcl(ii,j,1)
!soil
            ii=ii+1
        enddo
      enddo
      end if
!
      do 30 k=2,12
      if(jul .gt. mon(k-1) .and. jul .le. mon(k))then
!                                -- read dmsfile --
      write(lrec,11)ggdef,k-1
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,sstcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (sstcl(1,1,1),nxdef,nx,my)

!ch?  write(lrec,12)ggdef,k-1
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,1),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,1),nxdef,nx,my)

      write(lrec,14)ggdef,k-1
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,albcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (albcl(1,1,1),nxdef,nx,my)

      write(lrec,15)ggdef,k-1
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,z0cl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (z0cl(1,1,1),nxdef,nx,my)

      write(lrec,13)ggdef,k-1
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,wetcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (wetcl(1,1,1),nxdef,nx,my)
!soil
      write(lrec,19)ggdef,k-1
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,vfrcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (vfrcl(1,1,1),nxdef,nx,my)
!soil
      write(lrec,11)ggdef,k
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,sstcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (sstcl(1,1,2),nxdef,nx,my)

!ch?  write(lrec,12)ggdef,k
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,2),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,2),nxdef,nx,my)
!
      write(lrec,14)ggdef,k
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,albcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (albcl(1,1,2),nxdef,nx,my)

      write(lrec,15)ggdef,k
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,z0cl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (z0cl(1,1,2),nxdef,nx,my)

      write(lrec,13)ggdef,k
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,wetcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (wetcl(1,1,2),nxdef,nx,my)
!soil
      write(lrec,19)ggdef,k
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,vfrcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (vfrcl(1,1,2),nxdef,nx,my)
!soil

!
      coef1=float(jul-mon(k-1))/float(mon(k)-mon(k-1))
      coef2=1.-coef1
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            sst(i,jj) =coef1*sstcl(ii,j,2) +coef2*sstcl(ii,j,1)
!ch?        solt(i,jj)=coef1*soltcl(ii,j,2)+coef2*soltcl(ii,j,1)
            alb(i,jj) =coef1*albcl(ii,j,2) +coef2*albcl(ii,j,1)
            z0(i,jj)=coef1*z0cl(ii,j,2)+coef2*z0cl(ii,j,1)
            wet(i,jj)=coef1*wetcl(ii,j,2)+coef2*wetcl(ii,j,1)
!soil
            sigmaf(i,jj)=coef1*vfrcl(ii,j,2)+coef2*vfrcl(ii,j,1)
!soil
            ii=ii+1
         enddo
      enddo
      end if
  30  continue
!
!-- climat dataset has 2 half years
!                                  -- read dmsfile --
!
!-- land, sea and ice table
!                                   -- read dmsflie --
      write(lrec,17)ggdef,monidex
      call dmsreadi(nx,my,lrec,lncrec,'I',bckfile,iglob,istat)
      if( lreduce.eq.1 ) call reducepicki(iglob,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            icex(i,jj)=iglob(ii,j)
            ii=ii+1
         enddo
      enddo

      write(lrec,16)ggdef,blnk
      call dmsreadi(nx,my,lrec,lncrec,'I',bckfile,iglob,istat)
      if( lreduce.eq.1 ) call reducepicki(iglob,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            ls(i,jj)=iglob(ii,j)
            ii=ii+1
         enddo
      enddo

!soil
!-- soiltyp
!      write(lrec,18)ggdef,blnk
      if(isot .eq. 0)write(lrec,18)ggdef,blnk
      if(isot .eq. 1)write(lrec,28)ggdef,blnk
!
      call dmsreadi(nx,my,lrec,lncrec,'I',bckfile,iglob,istat)
      if( lreduce.eq.1 ) call reducepicki(iglob,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            istyp(i,jj)=iglob(ii,j)
            ii=ii+1
         enddo
      enddo
!-- vegtyp
!      write(lrec,21)ggdef,blnk
      if(ivegsrc .eq. 0)write(lrec,21)ggdef,blnk
      if(ivegsrc .eq. 1)write(lrec,31)ggdef,blnk
!
      call dmsreadi(nx,my,lrec,lncrec,'I',bckfile,iglob,istat)
      if( lreduce.eq.1 ) call reducepicki(iglob,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            ivegtyp(i,jj)=iglob(ii,j)
            ii=ii+1
         enddo
      enddo
!-- annual mean Tg
      write(lrec,22)ggdef,blnk
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (soltcl(1,1,2),nxdef,nx,my)
!soil
!noah
!-- shdmax
      write(lrec,23)ggdef,blnk
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,sstcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (sstcl(1,1,1),nxdef,nx,my)
!      call dmsread(nx,my,lrec,lncrec,'H',bckfile,shdmax,istat)
!      if( lreduce.eq.1 ) call reducepick (shdmax,nxdef,nx,my)
!-- shdmin
      write(lrec,24)ggdef,blnk
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,sstcl(1,1,2),istat)
      if( lreduce.eq.1 ) call reducepick (sstcl(1,1,2),nxdef,nx,my)
!      call dmsread(nx,my,lrec,lncrec,'H',bckfile,shdmin,istat)
!      if( lreduce.eq.1 ) call reducepick (shdmin,nxdef,nx,my)
!-- slopetyp
      write(lrec,25)ggdef,blnk
      call dmsreadi(nx,my,lrec,lncrec,'I',bckfile,iglob,istat)
      if( lreduce.eq.1 ) call reducepicki(iglob,nxdef,nx,my)
      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            slopetyp(i,jj)=iglob(ii,j)
            ii=ii+1
         enddo
      enddo
!-- snoalb
      write(lrec,26)ggdef,blnk
      call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,1),istat)
      if( lreduce.eq.1 ) call reducepick (soltcl(1,1,1),nxdef,nx,my)
!      call dmsread(nx,my,lrec,lncrec,'H',bckfile,snoalb,istat)
!      if( lreduce.eq.1 ) call reducepick (snoalb,nxdef,nx,my)

      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            shdmax(i,jj)=sstcl(ii,j,1)
            shdmin(i,jj)=sstcl(ii,j,2)
            snoalb(i,jj)=soltcl(ii,j,1)
            solt(i,jj)=soltcl(ii,j,2)
            if ( (ls(i,jj).eq.0) .and. (icex(i,jj).eq.1) ) solt(i,jj)=271.2
            ii=ii+1
        enddo
      enddo

      do jj=1,jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            ice(i,jj) =.false.
            land(i,jj)=.false.
            sea(i,jj) =.false.
            if(ls(i,jj) .eq. 0)sea(i,jj) =.true.
            if(ls(i,jj) .eq. 1)land(i,jj)=.true.
            if(ls(i,jj) .eq. 1)slopetyp(i,jj)=max(1,slopetyp(i,jj))
!soil if(icex(i,j) .eq. 1)then
!soil sea(i,j)=.false.
!soil land(i,j)=.false.
!soil ice(i,j)=.true.
!soil end if
        enddo
      enddo
!
!  modify albedo base on ground wetness
!
!      do 51 j=1,my
!       nxj=nxdef(j)
!      do 51 i=1,nxj
!     xt=0.31-0.17*wet(i,1)
!     if(alb(i,j).lt.xt)alb(i,j)=xt
!51   continue
!
!  modify deep soil temp for sea-ice points (=271.2)
!
      do 52 jj=1,jlistnum
        j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 52 i=1,nxj
      if ( (ls(i,jj).eq.0) .and. (icex(i,jj).eq.1) ) then
!soil
! in new soil model, ice present sea-ice
      ice(i,jj)=.true.        !sea ice
      sea(i,jj)=.false.
!soil
      endif
  52  continue
!
      return
      end
