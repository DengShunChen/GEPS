      subroutine readclx(nx,my,my_max,julian,land,sea,ice,solt,wet,z0  &
                        ,alb,sst,sigmaf,istyp,ivegtyp,ls               &
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
      use namelist_soilveg
!helio>
      use mpe
      use phygrid, only : ls_full,ls_redu
!helio<
      use const, only: ihdgi,bckfile,idtg2
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
      real    glob(nx,my)
      real      sstcl(nxp,my_max,2),soltcl(nxp,my_max,2),    &
                wetcl(nxp,my_max,2),albcl(nxp,my_max,2),     &
                z0cl(nxp,my_max,2),                          &
!soil
                vfrcl(nxp,my_max,2)
!soil
      integer   ls(nxp,my_max),icex(nxp,my_max),iglob(nx,my)
      integer   wla(nxp,my_max)
      real      wmx(nxp,my_max),wmxcl(nxp,my_max,2)

      character blnk*1,ggdef*4
      character cdate*12,cmmdd*4,cyymmdd*8  ! jwhwu
      integer   mon(12),mondy(13)
      data mon/15,46,74,105,135,166,196,227,258,288,319,349/
      data mondy/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data blnk/' '/

      integer   i,j,k,jul,nxj,mm,istat,monidex,lncrec,jj,ii
      real      coef1,coef2
!
      lncrec=nx*my
!jwhwu> 
      wla=0
      wmx=0.
      wmxcl=0.
      cdate=' '
!
      write(cdate,'(i12.12)') idtg2
      cyymmdd=cdate(1:8)
      cmmdd=cdate(5:8)
!jwhwu<
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
#ifdef I38K
  11  format('W00100','  gbck',a4,4x,i2.2,6x)  ! W10
  12  format('S9M100','  gbck',a4,4x,i2.2,6x)  ! X10
  14  format('S00030','  gbck',a4,4x,i2.2,6x)  ! S35
  15  format('S00040','  gbck',a4,4x,i2.2,6x)  ! S44
  16  format('S00070','  gbck',a4,11x,a1) ! S07
  17  format('S00090','  gbck',a4,4x,i2.2,6x) ! S09
!soil
  13  format('S9M5B0','  gbck',a4,4x,i2.2,6x)  ! for new soil
  18  format('S000ST','  gbck',a4,11x,a1)      ! for new soil
  19  format('S000VF','  gbck',a4,4x,i2.2,6x)  ! for new soil
  21  format('S000VT','  gbck',a4,11x,a1)      ! for new soil
  22  format('S9Y100','  gbck',a4,11x,a1)      ! for new soil
!soil
  23  format('S000VX','  gbck',a4,11x,a1)      ! for new soil
  24  format('S000VN','  gbck',a4,11x,a1)      ! for new soil
  25  format('S00063','  gbck',a4,11x,a1)      ! for new soil
  26  format('S0003X','  gbck',a4,11x,a1)      ! for new soil
!source 2(19-soil, 20-veg)
  28  format('S00XST','  gbck',a4,11x,a1)      ! for new soil
  31  format('S00XVT','  gbck',a4,11x,a1)      ! for new soil
!source 3(19-soil, 20-veg)
  29  format('S00WST','  gbck',a4,11x,a1)      ! for new soil from MODIS
  32  format('S00WVT','  gbck',a4,11x,a1)      ! for new soil from MODIS
  33  format('S00WVF','  gbck',a4,4x,i2.2,6x)      ! for new soil from MODIS
#else
  11  format('W00100','gbck',a4,4x,i2.2,6x)  ! W10
  12  format('S9M100','gbck',a4,4x,i2.2,6x)  ! X10
  14  format('S00030','gbck',a4,4x,i2.2,6x)  ! S35
  15  format('S00040','gbck',a4,4x,i2.2,6x)  ! S44
  16  format('S00070','gbck',a4,11x,a1) ! S07
  17  format('S00090','gbck',a4,4x,i2.2,6x) ! S09
  161 format('S00WLA','gbck',a4,11x,a1) ! S00WLA
  111 format('WMX100','gbck',a4,4x,i2.2,6x)  ! WMX100
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
!source 3(19-soil, 20-veg)
  29  format('S00WST','gbck',a4,11x,a1)      ! for new soil from MODIS
  32  format('S00WVT','gbck',a4,11x,a1)      ! for new soil from MODIS
  33  format('S00WVF','gbck',a4,4x,i2.2,6x)      ! for new soil from MODIS
#endif

!
!  to interpolat linearly based on julian day
!
!-- climat dataset has 12 months

      if(jul .le. mon(1))jul=jul+365
      if(jul .gt. mon(12))then
!                                -- read dmsfile --
      mm=12
#if (!defined DYCLM) && (!defined DYANL)
      write(ihdgi,11)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,sstcl(1,1,1))
!jwhwu>
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,111)ggdef,mm
        call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
        call unify_reducepick(nx,my,my_max,glob,wmxcl(1,1,1))
      endif
#endif
!jwhwu<
#endif

!ch?  write(lrec,12)ggdef,mm
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,1),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,1),nxdef,nx,my)

      write(ihdgi,14)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,albcl(1,1,1))

      write(ihdgi,15)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,z0cl(1,1,1))

      write(ihdgi,13)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,wetcl(1,1,1))

!soil
      if(ivegsrc .le. 1)write(ihdgi,19)ggdef,mm
      if(ivegsrc .eq. 2)write(ihdgi,33)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,vfrcl(1,1,1))
!soil
      mm=1
#if (!defined DYCLM) && (!defined DYANL)
      write(ihdgi,11)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,sstcl(1,1,2))
!jwhwu>
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,111)ggdef,mm
        call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
        call unify_reducepick(nx,my,my_max,glob,wmxcl(1,1,2))
      endif
#endif
#endif
!jwhwu<

!ch?  write(lrec,12)ggdef,mm
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,2),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,2),nxdef,nx,my)

      write(ihdgi,14)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,albcl(1,1,2))

      write(ihdgi,15)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,z0cl(1,1,2))

      write(ihdgi,13)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,wetcl(1,1,2))
!soil
      if(ivegsrc .le. 1)write(ihdgi,19)ggdef,mm
      if(ivegsrc .eq. 2)write(ihdgi,33)ggdef,mm
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,vfrcl(1,1,2))
!soil
!
      coef1=float(jul-mon(12))/float(380-mon(12))
      coef2=1.-coef1
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
#if (!defined DYCLM) && (!defined DYANL)
            sst(i,jj) =coef1*sstcl(i,jj,2) +coef2*sstcl(i,jj,1)
#endif
!ch?        solt(i,jj)=coef1*soltcl(ii,j,2)+coef2*soltcl(ii,j,1)
            alb(i,jj) =coef1*albcl(i,jj,2) +coef2*albcl(i,jj,1)
            z0(i,jj)=coef1*z0cl(i,jj,2)+coef2*z0cl(i,jj,1)
            wet(i,jj)=coef1*wetcl(i,jj,2)+coef2*wetcl(i,jj,1)
!soil
          sigmaf(i,jj)=coef1*vfrcl(i,jj,2)+coef2*vfrcl(i,jj,1)
!soil
        enddo
#if (!defined DYCLM) && (!defined DYANL)
!jwhwu>
#ifdef WMX
        if(isot .eq. 2) then
          do i=1,nxj
            wmx(i,jj) =coef1*wmxcl(i,jj,2) +coef2*wmxcl(i,jj,1)
          enddo
        endif
#endif
!jwhwu<
#endif
      enddo
      end if
!
      do 30 k=2,12
      if(jul .gt. mon(k-1) .and. jul .le. mon(k))then
!                                -- read dmsfile --
#if (!defined DYCLM) && (!defined DYANL)
      write(ihdgi,11)ggdef,k-1
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,sstcl(1,1,1))
!jwhwu>
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,111)ggdef,k-1
        call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
        call unify_reducepick(nx,my,my_max,glob,wmxcl(1,1,1))
      endif
#endif
#endif
!jwhwu<

!ch?  write(lrec,12)ggdef,k-1
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,1),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,1),nxdef,nx,my)

      write(ihdgi,14)ggdef,k-1
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,albcl(1,1,1))

      write(ihdgi,15)ggdef,k-1
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,z0cl(1,1,1))

      write(ihdgi,13)ggdef,k-1
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,wetcl(1,1,1))
!soil
      if(ivegsrc .le. 1)write(ihdgi,19)ggdef,k-1
      if(ivegsrc .eq. 2)write(ihdgi,33)ggdef,k-1
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,vfrcl(1,1,1))
!soil
#if (!defined DYCLM) && (!defined DYANL)
      write(ihdgi,11)ggdef,k
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,sstcl(1,1,2))
!jwhwu>
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,111)ggdef,k
        call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
        call unify_reducepick(nx,my,my_max,glob,wmxcl(1,1,2))
      endif
#endif
!jwhwu<
#endif

!ch?  write(lrec,12)ggdef,k
!ch?  call dmsread(nx,my,lrec,lncrec,'H',bckfile,soltcl(1,1,2),istat)
!ch?  if( lreduce.eq.1 ) call reducepick (soltcl(1,1,2),nxdef,nx,my)
!
      write(ihdgi,14)ggdef,k
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,albcl(1,1,2))

      write(ihdgi,15)ggdef,k
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,z0cl(1,1,2))

      write(ihdgi,13)ggdef,k
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,wetcl(1,1,2))
!soil
      if(ivegsrc .le. 1)write(ihdgi,19)ggdef,k
      if(ivegsrc .eq. 2)write(ihdgi,33)ggdef,k
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,vfrcl(1,1,2))
!soil

!
      coef1=float(jul-mon(k-1))/float(mon(k)-mon(k-1))
      coef2=1.-coef1
      do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
#if (!defined DYCLM) && (!defined DYANL)
            sst(i,jj) =coef1*sstcl(i,jj,2) +coef2*sstcl(i,jj,1)
#endif
!ch?        solt(i,jj)=coef1*soltcl(i,jj,2)+coef2*soltcl(i,jj,1)
            alb(i,jj) =coef1*albcl(i,jj,2) +coef2*albcl(i,jj,1)
            z0(i,jj)=coef1*z0cl(i,jj,2)+coef2*z0cl(i,jj,1)
            wet(i,jj)=coef1*wetcl(i,jj,2)+coef2*wetcl(i,jj,1)
!soil
            sigmaf(i,jj)=coef1*vfrcl(i,jj,2)+coef2*vfrcl(i,jj,1)
!soil
         enddo
#if (!defined DYCLM) && (!defined DYANL)
#ifdef WMX
         if(isot .eq. 2) then
           do i=1,nxj
             wmx(i,jj) =coef1*wmxcl(i,jj,2) +coef2*wmxcl(i,jj,1)
           enddo
         endif
#endif
#endif
      enddo
      end if
  30  continue
#if (defined DYCLM) || (defined DYANL)
#ifdef DYCLM
      write(ihdgi,'(a6,a4,a4,4x,a4,4x)')"W00100","gbck",ggdef,cmmdd
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,sst)
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,'(a6,a4,a4,4x,a4,4x)')"WMX100","gbck",ggdef,cmmdd
        call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
        call unify_reducepick(nx,my,my_max,glob,wmx)
      endif
#endif
#else
      write(ihdgi,'(a6,a4,a4,a8,4x)')"W00100","gbck",ggdef,cyymmdd
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,sst)
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,'(a6,a4,a4,a8,4x)')"WMx100","gbck",ggdef,cyymmdd
        call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
        call unify_reducepick(nx,my,my_max,glob,wmx)
      endif
#endif
#endif
#endif

!
!-- climat dataset has 2 half years
!                                  -- read dmsfile --
!
!-- land, sea and ice table
!                                   -- read dmsflie --
#ifdef DYCLM
      write(ihdgi,'(a6,a4,a4,4x,a4,4x)')"S00090","gbck",ggdef,cmmdd
#elif defined DYANL
      write(ihdgi,'(a6,a4,a4,a8,4x)')"S00090","gbck",ggdef,cyymmdd
#else
      write(ihdgi,17)ggdef,monidex
#endif
      call dmsreadi(nx,my,lncrec,'I',bckfile,iglob,istat)
      call unify_reducepicki(nx,my,my_max,iglob,icex)

      write(ihdgi,16)ggdef,blnk
      call dmsreadi(nx,my,lncrec,'I',bckfile,iglob,istat)
!helio>

      do jj = 1, jlistnum
         j=jlist1(jj)
         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         ls_full(:,jj) = iglob(:,j)
         if(lreduce.eq.1) call reducepicki (iglob(1,j),nxdef(j),nx,1)
         do i = 1, nxj
            ls(i,jj) = iglob(ii,j)
            ii=ii+1
         enddo
         ls_redu(:,jj) = iglob(:,j)
      enddo
!helio<
!jwhwu>
#ifdef WMX
      if(isot .eq. 2) then
        write(ihdgi,161)ggdef,blnk
        call dmsreadi(nx,my,lncrec,'I',bckfile,iglob,istat)
        do jj = 1, jlistnum
           j=jlist1(jj)
           ii=nxjstart(j)
           nxj=nxdef_2d(j)
           if(lreduce.eq.1) call reducepicki (iglob(1,j),nxdef(j),nx,1)
           do i = 1, nxj
              wla(i,jj) = iglob(ii,j)
              ii=ii+1
           enddo
        enddo
      endif
#endif
!jwhwu<


!soil
!-- soiltyp
!      write(lrec,18)ggdef,blnk
      if(isot .eq. 0)write(ihdgi,18)ggdef,blnk
      if(isot .eq. 1)write(ihdgi,28)ggdef,blnk
      if(isot .eq. 2)write(ihdgi,29)ggdef,blnk
!
      call dmsreadi(nx,my,lncrec,'I',bckfile,iglob,istat)
      call unify_reducepicki(nx,my,my_max,iglob,istyp)
!-- vegtyp
!      write(lrec,21)ggdef,blnk
      if(ivegsrc .eq. 0)write(ihdgi,21)ggdef,blnk
      if(ivegsrc .eq. 1)write(ihdgi,31)ggdef,blnk
      if(ivegsrc .eq. 2)write(ihdgi,32)ggdef,blnk
!
      call dmsreadi(nx,my,lncrec,'I',bckfile,iglob,istat)
      call unify_reducepicki(nx,my,my_max,iglob,ivegtyp)
! the value of z0 is followed from new veg.type, read in lookup table.
       if (ivegsrc .ge. 1)then
       do jj=1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
             z0(i,jj)=z0_data(ivegtyp(i,jj))
          enddo
       enddo
       endif
!-- annual mean Tg
      write(ihdgi,22)ggdef,blnk
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,solt)

!soil
!noah
!-- shdmax
      write(ihdgi,23)ggdef,blnk
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,shdmax)
!      call dmsread(nx,my,lrec,lncrec,'H',bckfile,shdmax,istat)
!      if( lreduce.eq.1 ) call reducepick (shdmax,nxdef,nx,my)
!-- shdmin
      write(ihdgi,24)ggdef,blnk
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,shdmin)
!      call dmsread(nx,my,lrec,lncrec,'H',bckfile,shdmin,istat)
!      if( lreduce.eq.1 ) call reducepick (shdmin,nxdef,nx,my)
!-- slopetyp
      write(ihdgi,25)ggdef,blnk
      call dmsreadi(nx,my,lncrec,'I',bckfile,iglob,istat)
      call unify_reducepicki(nx,my,my_max,iglob,slopetyp)
!-- snoalb
      write(ihdgi,26)ggdef,blnk
      call dmsread(nx,my,lncrec,'H',bckfile,glob,istat)
      call unify_reducepick(nx,my,my_max,glob,snoalb)
!      call dmsread(nx,my,lrec,lncrec,'H',bckfile,snoalb,istat)
!      if( lreduce.eq.1 ) call reducepick (snoalb,nxdef,nx,my)


      do jj=1,jlistnum
         j=jlist1(jj)
!byl         ii=nxjstart(j)
         nxj=nxdef_2d(j)
         do i=1,nxj
            ice(i,jj) =.false.
            land(i,jj)=.false.
            sea(i,jj) =.false.
            if(ls(i,jj) .eq. 0)then !ocean
              sea(i,jj) =.true.
              if( (icex(i,jj).eq.1) )then !sea ice
                 solt(i,jj)=271.2 ! deep soil temp for sea-ice points (=271.2)
                 ice(i,jj)=.true.        
                 sea(i,jj)=.false.
              endif
            elseif(ls(i,jj) .eq. 1)then !land
              land(i,jj)=.true.
              slopetyp(i,jj)=max(1,slopetyp(i,jj))
            endif

!soil if(icex(i,j) .eq. 1)then
!soil sea(i,j)=.false.
!soil land(i,j)=.false.
!soil ice(i,j)=.true.
!soil end if
        enddo
!jwhwu>
! for lake points: wmx >  273.15 --> water
!                  wmx <= 273.15 --> ice
#ifdef WMX
        if(isot .eq. 2) then
          do i=1,nxj
            if(wla(i,jj) .gt. 2) then
              if(wmx(i,jj) .gt. 273.15) then
                sea(i,jj)=.true.
                sst(i,jj)=wmx(i,jj)
                ice(i,jj)=.false.
                alb(i,jj)=0.09
              else
                ice(i,jj)=.true.
                solt(i,jj)=271.2
                sea(i,jj)=.false.
                alb(i,jj)=0.55
              endif
            endif
          enddo
        endif
#endif
!jwhwu<
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
!      do 52 jj=1,jlistnum
!        j=jlist1(jj)
!       nxj=nxdef_2d(j)
!      do 52 i=1,nxj
!      if ( (ls(i,jj).eq.0) .and. (icex(i,jj).eq.1) ) then
!!soil
!! in new soil model, ice present sea-ice
!      ice(i,jj)=.true.        !sea ice
!      sea(i,jj)=.false.
!!soil
!      endif
!  52  continue
!
      return
      end
