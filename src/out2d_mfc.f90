      subroutine out2d_mfc (nx,lev,my,my_max,itau,idtg    &
                            ,raincu1,rainlp1,raintot,glob,t2,q2,rh2,rh10 &
                            ,u10,v10,tmax,tmin,rld,sld,ctot,pt,ggdef )
!
      use rank
      use mpe
      use index
      use const ,only : grav,ptop,rgas,cp ,outdms ,outgrb2 ,ifilout_grb, &
                        RTYPE,kflag,out_hp,ihdgo,ihdgo2
      use grid  ,only : tt,qt,plt,pk,pk2,sgeo
      use mod_grb2_param , only :ofdir,wrt_grb2_v2,wrt_grb2_accu_v2
      use phygrid ,only: raincu3, rainlp3
      implicit  none

      integer   nx,lev,my,my_max,itau,ntau,num,nc
      parameter (num=14)

      real      raintot(nxp,my_max),t2(nxp,my_max),u10(nxp,my_max),   &
                v10(nxp,my_max),ctot(nxp,my_max)
      real(kind=RTYPE) pt(nxp,my_max)

      real                   q2(nxp,my_max),rh2(nxp,my_max),          &
           rh10(nxp,my_max),tmax(nxp,my_max),tmin(nxp,my_max),        &
                          rld(nxp,my_max),sld(nxp,my_max),            &
           raincu1(nxp,my_max),rainlp1(nxp,my_max)
!
      character*4 ggdef
      integer*8 idtg
      character*6 dmskey(num)
!     --- local variable
      real(kind=RTYPE) mfcout(nxp,my_max,num)

      integer,dimension(num):: ptp0 ,ptp1 ,ptp2 ,ptp3 ,ptp4 ,ptp5
!
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
!
      integer   n,levz,lenc,lenc2,i,ia,kk,j,nxj,istat,jj,llts,k
      real      tnshun
      real(kind=RTYPE) alaps,rdg,ttb,ttp,ttt,ttt1,ttt2,anlslp,apha
      real(kind=RTYPE) phi(nxp,lev,my_max),hld1(nxp,my_max),hld2(nxp,my_max)
! for due point temperature
      real(kind=RTYPE),parameter:: TdAlpha = 17.27 ,TdBeta = 237.7
      real(kind=RTYPE)      TdGamma
!
      data dmskey/'b00621','b0062t','b02100','b02500','b02510', &
                  'b10200','b10210','b02171','b02181','b02150', &
                  's003x0','s003u0','x00770','ssl010'/

      integer:: gtp1(9),gtp0(9)
!     grib code 0,1,2:variable   3:order  4:layer  5:above_land_height
!                 1h  Tot                   T2M T2M  2M DW DW       
!                 p   p  T2 sh2 rh2 u10 v10 max min DPT LW SW CC SLP
      data ptp0/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0,  0 /
      data ptp1/  1,  1,  0,  1,  1,  2,  2,  0,  0,  0, 5, 4, 6,  3 /
      data ptp2/  8, 49,  0,  0,  1,  2,  3,  4,  5,  6, 3, 7, 1,  1 /
      data ptp3/  2,  1,  2,  6,  2,  2,  2,  2,  2,  2, 2, 2, 3,  1 /
      data ptp4/103,103,103,103,103,103,103,103,103,103, 1, 1,10,101 /
      data ptp5/  0,  0,  2,  2,  2, 10, 10,  2,  2,  2, 0, 0, 0,  0 /

!
      ntau=itau
      tnshun= 1.0
      lenc= nx*my
      lenc2= lev*my
      alaps = 0.0065
      rdg = rgas/grav
!     
        if( outgrb2 == 1)then
 134                      format( A  ,A ,I10.10 , i4.4       )
             write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau
             if(myrank==0) call system("mkdir -p "//trim(ofdir) )
        endif
 
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          mfcout(i,jj,1)=raincu1(i,jj) + rainlp1(i,jj) !rain1
          mfcout(i,jj,2)=raintot(i,jj)
          mfcout(i,jj,3)=t2  (i,jj)
          mfcout(i,jj,4)=q2  (i,jj)
          mfcout(i,jj,5)=rh2 (i,jj) * 100.0  ! ( % )
          mfcout(i,jj,6)=u10 (i,jj)
          mfcout(i,jj,7)=v10 (i,jj)
          mfcout(i,jj,8)=tmax(i,jj)
          mfcout(i,jj,9)=tmin(i,jj)

          TdGamma = TdAlpha * (t2(i,jj)-273.15)  / ( TdBeta + ( t2(i,jj) - 273.15 ) ) + &
                     log(rh2(i,jj))
          mfcout(i,jj,10)= ( TdBeta * TdGamma / ( TdAlpha - TdGamma ) ) + 273.15
          mfcout(i,jj,11)=rld (i,jj)
          mfcout(i,jj,12)=sld (i,jj)
          mfcout(i,jj,13)=ctot(i,jj) * 100.0 ! total cloud cover
        enddo
      enddo

!
!  hydrostatic equation
!
      do jj =1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          phi(i,lev,jj)= cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj)) &
                       + sgeo(i,jj)
        enddo
        do k=lev-1,1,-1
          do i=1,nxj
            phi(i,k,jj)=phi(i,k+1,jj)+cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj)) &
                       + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
          enddo
        enddo
      enddo
!
! Sea level pressure(hPa)
!
!  compute sea level pressure
!  The method is based on one used by ecmwf, reseach manual 2 (1988)
!
!  llts layer's temperature is used to derive an alternative
!  surface skin temperature
!
      llts = lev-5
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          ttb  = tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
          ttp  = tt(i,llts,jj)*pk(i,llts,jj)/(1.0+0.608*qt(i,llts,jj))
          ttt1 = ttb + alaps*rdg*ttb*   &
              ((pt(i,jj)+ptop)/plt(i,lev,jj)-1.0)
          ttt2 = ttp + alaps*(phi(i,llts,jj)-sgeo(i,jj))/grav
          hld1(i,jj) = 0.25*ttt1 + 0.75*ttt2
          hld2(i,jj) = hld1(i,jj) + alaps*sgeo(i,jj)/grav
          if( sgeo(i,jj) .lt. 0.1 ) then
            anlslp = pt(i,jj) + ptop
          else if( hld1(i,jj) .le. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
            apha = rgas*(290.5-hld1(i,jj))/sgeo(i,jj)
            ttt = sgeo(i,jj)/(rgas*hld1(i,jj))
            anlslp = (pt(i,jj)+ptop)*exp(ttt*(1.0-0.5*apha*ttt+0.333333*  &
                     apha*ttt*apha*ttt) )
          else if( hld1(i,jj) .gt. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
            hld1(i,jj) = (hld1(i,jj)+290.5)*0.5
            anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
          else if( hld1(i,jj) .lt. 255.0 .and. hld2(i,jj) .lt. 255.0 ) then
            hld1(i,jj) = (hld1(i,jj)+255.0)*0.5
            anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
          else
            apha = alaps * rdg
            ttt = sgeo(i,jj)/(rgas*hld1(i,jj))
            anlslp = (pt(i,jj)+ptop)*exp(ttt*(1.0-0.5*apha*ttt+0.333333*  &
                     apha*ttt*apha*ttt) )
          endif
          mfcout(i,jj,14) = anlslp * 100.0  !hPa to Pa
        enddo
      enddo
!
! Total Precp.
!byl      call mpe2d_unify(glob,raintot)
      nc=0
      do n=1,num
        call syslbl_w (dmskey(n),idtg,ntau,ggdef)
        call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),glob)
        call qmaxn3_w (glob,1,1,1,nx,my,1)
!        if ( myrank .eq. n-1 ) then
!          mout=glob
!          ihdgo2=ihdgo
          gtp0=(/ptp0(n),ptp1(n),ptp2(n),ptp3(n),ptp4(n),0,ptp5(n),-999,-999/)
          if(n==1)then
             gtp0(8:9)=(/1,1/) !1hr precip
          else if(n==8)then
             gtp0(8:9)=(/2,1/) !MaxT2m
          else if(n==9)then
             gtp0(8:9)=(/3,1/) !MinT2m
          endif
!        endif
          call split2(nx,my,lenc,nc,glob,mout,gtp0,gtp1)
      enddo
!
      if (myrank .lt. nc ) then

        if(outgrb2 == 1 )then
          if(gtp1(8)==-999)then
          call wrt_grb2_v2(itau,gtp1(1),gtp1(2),gtp1(3),gtp1(4),gtp1(5) &
              ,gtp1(6),gtp1(7),mout)
          else
          call wrt_grb2_accu_v2(itau,gtp1(1),gtp1(2),gtp1(3),gtp1(4),gtp1(5) &
              ,gtp1(6),gtp1(7),gtp1(8),gtp1(9),mout)
          endif
        endif !outgrb2

        if(outdms.gt.0)then
          if ( myrank .eq. (14-1) ) mout = mout / 100.0
          call dmswrit_split(nx,my,lenc,kflag,mout,istat)
        endif ! outdms .gt. 0

      endif
!
!! rh10
!!byl      call mpe2d_unify(glob,rh10)
!      call syslbl ('b10510',idtg,ntau,ggdef,ihdg)
!      call unify_reduceintp(nx,my,my_max,rh10,glob)
!!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!      glob=glob*100.0
!!     call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
!      call dmswrit_mfc(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)

! 3hr accum. precipitation for HomePageWifi
      if (out_hp)then
       do jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i=1,nxj
           raincu3(i,jj)= raincu3(i,jj) + raincu1(i,jj)
           rainlp3(i,jj)= rainlp3(i,jj) + rainlp1(i,jj)
         enddo
       enddo

       if( mod( itau , 3 ) == 0 )then

       if(outgrb2 == 1 )then
        if(myrank==0)then
        !convective precipitation
        call syslbl_w ('B00632',idtg,ntau,ggdef)
        glob=raincu3
        call unify_reduceintp(nx,my,my_max,glob,mout)
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,1,10,2,103,0,0,1,3,mout)
        !
        call syslbl_w ('B00642',idtg,ntau,ggdef)
        glob=rainlp3
        call unify_reduceintp(nx,my,my_max,glob,mout)
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,1,9,2,103,0,0,1,3,mout)

        call syslbl_w ('B00622',idtg,ntau,ggdef)
        glob=raincu3 + rainlp3
        call unify_reduceintp(nx,my,my_max,glob,mout)
        ihdgo2 = ihdgo
        call wrt_grb2_accu_v2(itau,0,1,9,2,103,0,0,1,3,mout)

        endif
       endif

        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            raincu3(i,jj)= 0.0
            rainlp3(i,jj)= 0.0
          enddo
        enddo

       endif ! mod(itau,3)==0

      endif !out_hp

      return
      end

