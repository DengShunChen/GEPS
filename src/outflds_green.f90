      subroutine outflds_green( itau,nx,my,my_max,lev,ncld       &
             , idtg,cp,rgas,grav,t2,u10,v10,ss,pk                &
             , sgeo,pt,plt,ptop,ut,vt,tt,qt,cosl,raincu6,rainlp6)
!
!  output driver subroutine to process sigma level data to 40m & 100m
!
      use mpe
      use rank
      use index
      use const ,only:aki,bki ,outdms ,outgrb2 ,ifilout_grb , &
                      RTYPE,kflag,ggdef
      use mod_grb2_param , only :ofdir,wrt_grb2_v2,wrt_grb2_accu_v2

      use noah,only:runoff, cice ,zice 
      use phygrid,only:ice,ocean,land

      implicit  none

      integer   itau,nx,my,my_max,lev,ncld,istat,lenc,nc
      real      ptop,cp,rgas,grav
      real      rcp

      real      plt(nxp,lev,my_max)                                    &
              , u10(nxp,my_max),v10(nxp,my_max),t2(nxp,my_max)         &
              , ss(nxp,my_max)                                         &
              , tht(nxp,my_max),raincu6(nxp,my_max),rainlp6(nxp,my_max)
      real(kind=RTYPE) ut(nxp,lev,my_max),vt(nxp,lev,my_max),          &
                       tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max),     &
                       sgeo(nxp,my_max),pt(nxp,my_max),cosl(my),       &
                       pk(nxp,lev,my_max)

      integer*8 idtg
      integer:: ptp0(9),ptp1(9)
!
! local work arrays
!
      real(kind=RTYPE) glob(nx,my),mout(nx,my),wrk(nxp,my_max),      &
                       rh0(nxp,my_max)
!
      real      whtlev(100),whtlevq(100),whtlevz(100)
      character*6 labx

      integer   jj,j,nxj,k,i,n,nk,ntrac,ll,mm,la,kk
!
      real oqt(nxp,my_max),oqc(nxp,my_max),ou(nxp,my_max),  &
           ov(nxp,my_max),ot(nxp,my_max),pla(nxp,my_max),   &
           pp(nxp,my_max),p2(nxp,my_max),p10(nxp,my_max),   &
           rhtmp(nxp)
      real, parameter ::rad=6.371e6
      integer,parameter :: l= 4, m= 2
      real   avett,p(l),hm(m),xxx,temp,tepl(nxp,l,my_max)
      real(kind=RTYPE) akir(l),bkir(l)
      data hm/100.0,40.0/
!      data aki/   .00000,   .02193,   .26557,   .97701/ !M60~M57
!      data bki/.99058760,.98124505,.96996497,.95697164/
      character layer(2)*3,var(6)*3,wtemp*6
      data layer/'H10','B40'/
      data var/'010','500','200','210','100','550'/
!
      pla=0.
      oqt=0.
      oqc=0.
      ou=0.
      ov=0.
      ot=0.
      runoff=0.
      rcp=rgas/cp
      lenc = nx*my
      nc=0
      do k=1,l
        kk=lev-k+1
        akir(k)=aki(kk)
        bkir(k)=bki(kk)
      enddo

        if( outgrb2 == 1)then
 134                      format( A  ,A ,I10.10 , i4.4       )
             write(ofdir,134 )trim(ifilout_grb),'/',idtg/100 ,itau
             if(myrank==0) call system("mkdir -p "//trim(ofdir) )
        endif
!=======================================================================
      if(myrank .eq. 0) print*,'   in outflds_green for tau= ',itau
!----------------------------------------------------------------------
      do mm=1,m    ! 1: 100m, 2: 40m.
!-----------------------------------------------------------------------
! calculate PQUVT at new layer.
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
          do i = 1, nxj
! ave. temp.
          avett=(tt(i,lev,jj)+tt(i,lev-1,jj))/2.*                  &
                ((plt(i,lev,jj)+plt(i,lev-1,jj))/2./1000.)**rcp/     &
                (1.+0.608*(qt(i,lev-1,jj)+qt(i,lev-1,jj))/2.)
! calculate 40m&100m P, 2m P, 10m P.
        pp(i,jj)=(pt(i,jj)+ptop)*(exp(-(hm(mm)*grav)/(rgas*avett)))
        p2(i,jj)=(pt(i,jj)+ptop)*(exp(-(2.0*grav)/(rgas*avett)))
        p10(i,jj)=(pt(i,jj)+ptop)*(exp(-(10.0*grav)/(rgas*avett)))
        pla(i,jj)=pp(i,jj)
! transfer temperature from tt(virtual theta) to tepl(temperature).
        la=0
        do ll=lev,lev-3,-1
          la=la+1
          tepl(i,la,jj)=(tt(i,ll,jj)*(pla(i,jj)/1000.)**rcp)/(1.+0.608*qt(i,ll,jj))
        enddo

        call sigmap(l,akir(:),bkir(:),pt(i,jj),p)
! find pp pressure position
         if(pp(i,jj).gt.p(1))then
          temp=(pp(i,jj)-p(2))/(p(1)-p(2))-1.0
          oqt(i,jj)=qt(i,lev,jj)+temp*(qt(i,lev,jj)-qt(i,lev-1,jj))
          oqt(i,jj)=max(oqt(i,jj),1.e-8)
          oqc(i,jj)=qt(i,lev*2,jj)+temp*(qt(i,lev*2,jj)-qt(i,lev*2-1,jj))
          oqc(i,jj)=max(oqc(i,jj),1.e-8)

         xxx= rad/cosl(j)
         ou(i,jj)=u10(i,jj)/xxx
         ov(i,jj)=v10(i,jj)/xxx
          temp=(p(1)-pp(i,jj))/(p(1)-p10(i,jj))
          ou(i,jj)=ou(i,jj)*temp+ut(i,lev,jj)*(1.0-temp)
          ov(i,jj)=ov(i,jj)*temp+vt(i,lev,jj)*(1.0-temp)

          temp=(p(1)-pp(i,jj))/(p(1)-p2(i,jj))
          ot(i,jj)=t2(i,jj)*temp+tepl(i,1,jj)*(1.0-temp)

          goto 99
          endif

        do ll=1,l-1
         if((pp(i,jj).le.p(ll)).and.(pp(i,jj)).gt.p(ll+1))then
          temp=(p(ll)-pp(i,jj))/(p(ll)-p(ll+1))
          oqt(i,jj)=qt(i,lev-ll+1,jj)*temp+qt(i,lev-ll,jj)*(1.0-temp)
          oqt(i,jj)=max(oqt(i,jj),1.e-8)
          oqc(i,jj)=qt(i,lev*2-ll+1,jj)*temp+qt(i,lev*2-ll,jj)*(1.0-temp)
          oqc(i,jj)=max(oqc(i,jj),1.e-8)
          ou(i,jj)=ut(i,lev-ll+1,jj)*temp+ut(i,lev-ll,jj)*(1.0-temp)
          ov(i,jj)=vt(i,lev-ll+1,jj)*temp+vt(i,lev-ll,jj)*(1.0-temp)
          ot(i,jj)=tepl(i,ll,jj)*temp+tepl(i,ll+1,jj)*(1.0-temp)
          goto 99
         endif
        enddo

         if(pp(i,jj).le.p(l))then
          temp=(p(l-1)-pp(i,jj))/(p(l-1)-p(l))-1.0
          oqt(i,jj)=qt(i,lev-l+1,jj)+temp*(qt(i,lev-l+1,jj)-qt(i,lev-l+2,jj))
          oqt(i,jj)=max(oqt(i,jj),1.e-8)
          oqc(i,jj)=qt(i,lev*2-l+1,jj)+temp*(qt(i,lev*2-l+1,jj)-qt(i,lev*2-l+2,jj))
          oqc(i,jj)=max(oqc(i,jj),1.e-8)
          ou(i,jj)=ut(i,lev-l+1,jj)+temp*(ut(i,lev-l+1,jj)-ut(i,lev-l+2,jj))
          ov(i,jj)=vt(i,lev-l+1,jj)+temp*(vt(i,lev-l+1,jj)-vt(i,lev-l+2,jj))
          ot(i,jj)=tepl(i,l,jj)+temp*(tepl(i,l,jj)-tepl(i,l-1,jj))
          goto 99
         endif
!-----------------------------------------------------------------------
 99   continue
! transfer tt(theta tv) to T
!        tepl(i,j)=(ot(i,j)*(pla(i,j)/1000.)**rcp)/(1.+0.608*oqt(i,j))

        enddo  ! end (i)
      enddo  ! end (jj)

!      call mpe_unify(pla,nx,my,2,mpe_double)
!      call mpe_unify(oqt,nx,my,2,mpe_double)
!      call mpe_unify(oqc,nx,my,2,mpe_double)
!      call mpe_unify(tepl,nx,my,2,mpe_double)
!      call mpe_unify(ot,nx,my,2,mpe_double)
      do 50 jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      xxx= rad/cosl(j)
      do 50 i=1,nxj
       ou(i,jj)= ou(i,jj)*xxx
       ov(i,jj)= ov(i,jj)*xxx
   50 continue

!      call mpe_unify(globu,nx,my,2,mpe_double)
!      call mpe_unify(globv,nx,my,2,mpe_double)
!-----------------------------------------------------------------------

!output P
      write(wtemp,'(a3,a3)')layer(mm),var(1)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=pla * 100.0
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,3,0,1,103,0,nint(hm(mm)),-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output Q
      write(wtemp,'(a3,a3)')layer(mm),var(2)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=oqt
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,1,0,6,103,0,nint(hm(mm)),-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      write(wtemp,'(a3,a3)')layer(mm),var(6)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=oqc
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,1,235,8,103,0,nint(hm(mm)),-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output U,V
      write(wtemp,'(a3,a3)')layer(mm),var(3)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=ou
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,2,2,2,103,0,nint(hm(mm)),-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      write(wtemp,'(a3,a3)')layer(mm),var(4)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=ov
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,2,3,2,103,0,nint(hm(mm)),-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output T
      write(wtemp,'(a3,a3)')layer(mm),var(5)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=ot
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,0,0,2,103,0,nint(hm(mm)),-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!-----------------------------------------------------------------------
      enddo  ! end (mm)
!=======================================================================
!output S00310(net SW flux at the surface)
      write(wtemp,'(a6)')'S00310'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      wrk=ss
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,4,9,2,1,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output RH at bottom level
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
          do i = 1, nxj
          tht(i,jj) =tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
          enddo
       call qsatq(nxj,tht(1,jj),plt(1,lev,jj),rhtmp)
          do i = 1, nxj
           rh0(i,jj)=rhtmp(i)
           rh0(i,jj)=100.*(qt(i,lev,jj)/rh0(i,jj))
           rh0(i,jj)= min( 100., max( 1., rh0(i,jj) ) )
          enddo
      enddo
      write(wtemp,'(a6)')'B00510'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,rh0,glob)
      ptp0=(/0,1,1,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output b00010
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
          do i = 1, nxj
           wrk(i,jj)=pt(i,jj)+ptop
           wrk(i,jj)=wrk(i,jj) * 100.0
          enddo
      enddo
      write(wtemp,'(a6)')'B00010'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/0,3,0,1,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!=======================================================================
!output 6hr prec.
      if (mod(float(itau)+0.00001, 6. ) .lt. 0.01) then
      call syslbl_w ('b00633',idtg,itau,ggdef)
      wrk=raincu6
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,10,2,103,0,0,1,6/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!
      call syslbl_w ('b00643',idtg,itau,ggdef)
      wrk=rainlp6
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,47,2,103,0,0,1,6/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!
      call syslbl_w ('b00623',idtg,itau,ggdef)
      do 98 jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 98 i=1,nxj
       wrk(i,jj)=raincu6(i,jj)+rainlp6(i,jj)
 98   continue
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,8,2,103,0,0,1,6/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      !runoff
      wrk=runoff
      call syslbl_w ('runoff',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,5,2,103,0,0,1,6/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      runoff(:,:)=0.0 !6 hr zero out

      endif !mod(float(itau)+0.00001, 6. ) .lt. 0.01

!
      if (mod(float(itau)+0.00001, 24. ) .lt. 0.01) then
!      !!!sea ice fraction
!      wrk=cice
!      call syslbl_w ('w00091',idtg,itau,ggdef)
!      call unify_reduceintp(nx,my,my_max,wrk,glob)
!      ptp0=(/10,2,0,2,103,0,0,-999,-999/)
!      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      !!!sea ice thickness
      wrk=zice
      call syslbl_w ('w00092',idtg,itau,ggdef)
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      ptp0=(/10,2,1,2,103,0,0,-999,-999/)
      call split2(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
      endif

!
      if(outdms.gt.0)then
      if ( myrank .lt. nc )             &
         call dmswrit_split(nx,my,lenc,kflag,mout,istat)
      endif

      if(outgrb2 == 1 )then
       if ( myrank .lt. nc ) then
        if(ptp1(8)==-999)then
        call wrt_grb2_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),mout)
        else
        call wrt_grb2_accu_v2(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
        endif
       endif
      endif
!
!=======================================================================
      return
      end

!***********************************************************************
      subroutine sigmap(layer,aki,bki,psfc,p)
!
      use const, only: RTYPE
!
      implicit none
      integer i, layer
      real p(layer)
      real(kind=RTYPE) psfc,aki(layer), bki(layer)

      do i=1,layer
      p(i)=aki(i)+(bki(i)*psfc)
      enddo
      return
      end
