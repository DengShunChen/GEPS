      subroutine out2d_mfc (nx,lev,my,my_max,ifilout,itau,idtg    &
                            ,raincu1,rainlp1,raintot,glob,t2,q2,rh2,rh10 &
                            ,u10,v10,tmax,tmin,rld,sld,ctot,pt,ggdef )
!
      use rank
      use mpe
      use index
      use const ,only : grav,ptop,rgas,cp ,outdms ,outgrb2 ,ifilout_grb, &
                        RTYPE,kflag
      use grid  ,only : tt,qt,plt,pk,pk2,sgeo
      use mod_grb2_param  !for write grib2 data
!
      implicit  none

      integer   nx,lev,my,my_max,itau,ntau,num
      parameter (num=14)

      real      raintot(nxp,my_max),t2(nxp,my_max),u10(nxp,my_max),   &
                v10(nxp,my_max),ctot(nxp,my_max)
      real(kind=RTYPE) pt(nxp,my_max)

      real                   q2(nxp,my_max),rh2(nxp,my_max),          &
           rh10(nxp,my_max),tmax(nxp,my_max),tmin(nxp,my_max),        &
                          rld(nxp,my_max),sld(nxp,my_max),            &
           raincu1(nxp,my_max),rainlp1(nxp,my_max)
!
      real(kind=RTYPE) mfcout(nxp,my_max,num)
!
      character*4 ggdef
      integer*8 idtg
      character*6 dmskey(num)

      integer,dimension(num):: ptp0 ,ptp1 ,ptp2 ,ptp3 ,ptp4 ,ptp5
!
      real(kind=RTYPE) glob(nx,my),mout(nx,my)
!
      character*80 ifilout
      character*26 ihdg,ihdg2
!
      integer   n,levz,lenc,lenc2,i,ia,kk,j,nxj,istat,jj,llts,k
      real      tnshun,alaps,rdg,ttb,ttp,ttt,ttt1,ttt2,anlslp,apha
      real      phi(nxp,lev,my_max),hld1(nxp,my_max),hld2(nxp,my_max)
! for due point temperature
      real,parameter:: TdAlpha = 17.27 ,TdBeta = 237.7
      real      TdGamma
!
      data dmskey/'b00621','b0062t','b02100','b02500','b02510', &
                  'b10200','b10210','b02171','b02181','b02150', &
                  's003x0','s003u0','x00770','ssl010'/

!     grib code 0,1,2:variable   3:order  4:layer  5:above_land_height
!                 1h  Tot                   T2M T2M  2M DW DW       
!                 p   p  T2 sh2 rh2 u10 v10 max min DPT LW SW CC SLP
      data ptp0/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0,  0 /
      data ptp1/  1,  1,  0,  1,  1,  2,  2,  0,  0,  0, 5, 4, 6,  3 /
      data ptp2/  8, 49,  0,  0,  1,  2,  3,  4,  5,  6, 3, 7, 1,  1 /
      data ptp3/  2,  1,  2,  6,  2,  2,  2,  2,  2,  2, 2, 2, 3,  2 /
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
          mfcout(i,jj,13)=ctot(i,jj) ! total cloud cover
        enddo
      enddo

     !mfcout(:,:,1)=rain1(:,:)
     !mfcout(:,:,2)=raintot(:,:)
     !mfcout(:,:,3)=t2(:,:)
     !mfcout(:,:,4)=q2(:,:)
     !mfcout(:,:,5)=rh2(:,:) * 100.0
     !mfcout(:,:,6)=u10(:,:)
     !mfcout(:,:,7)=v10(:,:)
     !mfcout(:,:,8)=tmax(:,:)
     !mfcout(:,:,9)=tmin(:,:)
     !mfcout(:,:,10)=td(:,:)
     !mfcout(:,:,11)=rld(:,:)
     !mfcout(:,:,12)=sld(:,:)
     !mfcout(:,:,13)=ctot(:,:)
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
          mfcout(i,jj,14) = anlslp
        enddo
      enddo
!
! Total Precp.
!byl      call mpe2d_unify(glob,raintot)

!====== grib2 output
      if(outgrb2==1 )then
          n=1
          call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),mout)
          if(myrank==0)then
            call wrt_grb2_accu(ntau,ptp0(n),ptp1(n),ptp2(n),ptp3(n),ptp4(n),0,float(ptp5(n)),1,1,mout ) 
          endif
          do n=2,7
            call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),mout)
            if(myrank==0)then
              call wrt_grb2(ntau,ptp0(n),ptp1(n),ptp2(n),ptp3(n),ptp4(n),0,float(ptp5(n)) ,mout ) 
            endif
          enddo

          n=8          !Tmax2m 
          call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),mout)
          if(myrank==0)call wrt_grb2_accu(ntau,ptp0(n),ptp1(n),ptp2(n),ptp3(n),ptp4(n),0,float(ptp5(n)),2,1,mout ) 
          n=9          !Tmin2m
          call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),mout)
          if(myrank==0)call wrt_grb2_accu(ntau,ptp0(n),ptp1(n),ptp2(n),ptp3(n),ptp4(n),0,float(ptp5(n)),3,1,mout ) 

          do n=10,num
            call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),mout)
            if(myrank==0)then
              call wrt_grb2(ntau,ptp0(n),ptp1(n),ptp2(n),ptp3(n),ptp4(n),0,float(ptp5(n)) ,mout ) 
            endif
          enddo
      endif


    if(outdms.gt.0)then
      do n=1,num
        call syslbl (dmskey(n),idtg,ntau,ggdef,ihdg)
        call unify_reduceintp(nx,my,my_max,mfcout(1,1,n),glob)
        if ( myrank .eq. n-1 ) then
          mout=glob
          ihdg2=ihdg
        endif
      enddo
!
      if (myrank .lt. num ) call dmswrit_split(nx,my,ihdg2,lenc,kflag,ifilout,mout,istat)
    endif ! outdms .gt. 0
!
!! rh10
!!byl      call mpe2d_unify(glob,rh10)
!      call syslbl ('b10510',idtg,ntau,ggdef,ihdg)
!      call unify_reduceintp(nx,my,my_max,rh10,glob)
!!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
!      glob=glob*100.0
!!     call dmswrit(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)
!      call dmswrit_mfc(nx,my,ihdg,lenc,kflag,ifilout,glob,istat)


      return
      end

