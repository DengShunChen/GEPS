      subroutine outflds_green( itau,nx,my,my_max,lev,ncld       &
             , idtg,ifilout,cp,rgas,grav,t2,u10,v10              &
             , sgeo,pt,plt,ptop,ut,vt,tt,qt,cosl,raincu6,rainlp6 &
             , ggdef,lwrite)
!
!  output driver subroutine to process sigma level data to 40m & 100m
!
      use mpe
      use rank
      use index

      implicit  none

      integer   itau,nx,my,my_max,lev,ncld,istat,lenc
      real      ptop,cp,rgas,grav
      real      rcp

      real      sgeo(nxp,my_max),pt(nxp,my_max),plt(nxp,lev,my_max)      &
              , ut(nxp,lev,my_max),vt(nxp,lev,my_max),tt(nxp,lev,my_max) &
              , qt(nxp,lev*ncld,my_max),cosl(my)                       &
              , u10(nxp,my_max),v10(nxp,my_max),t2(nxp,my_max)           &
              , raincu6(nxp,my_max),rainlp6(nxp,my_max)

      character ifilout*80, ggdef*4, ihdg*26
      integer*8 idtg
!
! local work arrays
!
      real      slp(nx,my),glob(nx,my),glob1(nx,my)                   &
                ,globu(nx,my),globv(nx,my)
!
      real      whtlev(100),whtlevq(100),whtlevz(100)
      character*6 labx

      integer   jj,j,nxj,k,i,n,nk,ntrac,ll,mm,la
!
      real oqt(nx,my),oqc(nx,my),ou(nx,my),ov(nx,my),ot(nx,my),pla(nx,my) &
           ,wrk(nx,1),pp(nx,my),p2(nx,my),p10(nx,my)
      real, parameter ::rad=6.371e6
      integer,parameter :: l= 4, m= 2
      real   avett,p(l),hm(m),aki(l),bki(l),xxx,temp,tepl(nx,l,my)
      data hm/100.0,40.0/
      data aki/   .00000,   .02193,   .26557,   .97701/ !M60~M57
      data bki/.99058760,.98124505,.96996497,.95697164/
      character layer(2)*3,var(6)*3,wtemp*6
      data layer/'H10','B40'/
      data var/'010','500','200','210','100','550'/
      logical :: lwrite
!
      rcp=rgas/cp
      lenc = nx*my
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
        pp(i,j)=(pt(i,jj)+ptop)*(exp(-(hm(mm)*grav)/(rgas*avett)))
        p2(i,j)=(pt(i,jj)+ptop)*(exp(-(2.0*grav)/(rgas*avett)))
        p10(i,j)=(pt(i,jj)+ptop)*(exp(-(10.0*grav)/(rgas*avett)))
         pla(i,j)=pp(i,j)
! transfer temperature from tt(virtual theta) to tepl(temperature).
          la=0
           do ll=lev,lev-3,-1
           la=la+1
        tepl(i,la,j)=(tt(i,ll,jj)*(pla(i,j)/1000.)**rcp)/(1.+0.608*qt(i,ll,jj))
           enddo

        call sigmap(l,aki(:),bki(:),pt(i,jj),p)
! find pp pressure position
         if(pp(i,j).gt.p(1))then
          temp=(pp(i,j)-p(2))/(p(1)-p(2))-1.0
          oqt(i,j)=qt(i,lev,jj)+temp*(qt(i,lev,jj)-qt(i,lev-1,jj))
          oqt(i,j)=max(oqt(i,j),1.e-8)
          oqc(i,j)=qt(i,lev*2,jj)+temp*(qt(i,lev*2,jj)-qt(i,lev*2-1,jj))
          oqc(i,j)=max(oqc(i,j),1.e-8)

         xxx= rad/cosl(j)
         ou(i,j)=u10(i,jj)/xxx
         ov(i,j)=v10(i,jj)/xxx
          temp=(p(1)-pp(i,j))/(p(1)-p10(i,j))
          ou(i,j)=ou(i,j)*temp+ut(i,lev,jj)*(1.0-temp)
          ov(i,j)=ov(i,j)*temp+vt(i,lev,jj)*(1.0-temp)

          temp=(p(1)-pp(i,j))/(p(1)-p2(i,j))
          ot(i,j)=t2(i,jj)*temp+tepl(i,1,j)*(1.0-temp)

          goto 99
          endif

        do ll=1,l-1
         if((pp(i,j).le.p(ll)).and.(pp(i,j)).gt.p(ll+1))then
          temp=(p(ll)-pp(i,j))/(p(ll)-p(ll+1))
          oqt(i,j)=qt(i,lev-ll+1,jj)*temp+qt(i,lev-ll,jj)*(1.0-temp)
          oqt(i,j)=max(oqt(i,j),1.e-8)
          oqc(i,j)=qt(i,lev*2-ll+1,jj)*temp+qt(i,lev*2-ll,jj)*(1.0-temp)
          oqc(i,j)=max(oqc(i,j),1.e-8)
          ou(i,j)=ut(i,lev-ll+1,jj)*temp+ut(i,lev-ll,jj)*(1.0-temp)
          ov(i,j)=vt(i,lev-ll+1,jj)*temp+vt(i,lev-ll,jj)*(1.0-temp)
          ot(i,j)=tepl(i,ll,j)*temp+tepl(i,ll+1,j)*(1.0-temp)
          goto 99
         endif
        enddo

         if(pp(i,j).le.p(l))then
          temp=(p(l-1)-pp(i,j))/(p(l-1)-p(l))-1.0
          oqt(i,j)=qt(i,lev-l+1,jj)+temp*(qt(i,lev-l+1,jj)-qt(i,lev-l+2,jj))
          oqt(i,j)=max(oqt(i,j),1.e-8)
          oqc(i,j)=qt(i,lev*2-l+1,jj)+temp*(qt(i,lev*2-l+1,jj)-qt(i,lev*2-l+2,jj))
          oqc(i,j)=max(oqc(i,j),1.e-8)
          ou(i,j)=ut(i,lev-l+1,jj)+temp*(ut(i,lev-l+1,jj)-ut(i,lev-l+2,jj))
          ov(i,j)=vt(i,lev-l+1,jj)+temp*(vt(i,lev-l+1,jj)-vt(i,lev-l+2,jj))
          ot(i,j)=tepl(i,l,j)+temp*(tepl(i,l,j)-tepl(i,l-1,j))
          goto 99
         endif
!-----------------------------------------------------------------------
 99   continue
! transfer tt(theta tv) to T
!        tepl(i,j)=(ot(i,j)*(pla(i,j)/1000.)**rcp)/(1.+0.608*oqt(i,j))

        enddo  ! end (i)
      enddo  ! end (jj)

      call mpe_unify(pla,nx,my,2,mpe_double)
      call mpe_unify(oqt,nx,my,2,mpe_double)
      call mpe_unify(oqc,nx,my,2,mpe_double)
!      call mpe_unify(tepl,nx,my,2,mpe_double)
      call mpe_unify(ot,nx,my,2,mpe_double)
      do 50 j=1,my
      xxx= rad/cosl(j)
      do 50 i=1,nx
      globu(i,j)= ou(i,j)*xxx
      globv(i,j)= ov(i,j)*xxx
   50 continue

      call mpe_unify(globu,nx,my,2,mpe_double)
      call mpe_unify(globv,nx,my,2,mpe_double)
!-----------------------------------------------------------------------

!output P
      write(wtemp,'(a3,a3)')layer(mm),var(1)
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (pla,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,pla,istat)

!output Q
      write(wtemp,'(a3,a3)')layer(mm),var(2)
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (oqt,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,oqt,istat)

      write(wtemp,'(a3,a3)')layer(mm),var(6)
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (oqc,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,oqc,istat)

!output U,V
      write(wtemp,'(a3,a3)')layer(mm),var(3)
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (globu,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,globu,istat)

      write(wtemp,'(a3,a3)')layer(mm),var(4)
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (globv,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,globv,istat)

!output T
      write(wtemp,'(a3,a3)')layer(mm),var(5)
      call syslbl(wtemp,idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (ot,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,ot,istat)
!-----------------------------------------------------------------------
      enddo  ! end (mm)
!=======================================================================
      call mpe2d_unify(glob,raincu6)
      call mpe2d_unify(glob1,rainlp6)
      call syslbl ('b00633',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!
      call syslbl ('b00643',idtg,itau,ggdef,ihdg)
      if( lreduce.eq.1 ) call reduceintp (glob1,nxdef,nx,my)
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob1,istat)
      call qmaxn3 (glob1,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)
!
      call syslbl ('b00623',idtg,itau,ggdef,ihdg)
      do 98 j=1,my
      do 98 i=1,nx
       glob(i,j)=glob(i,j)+glob1(i,j)
 98   continue
      if(lwrite) call dmswrit(nx,my,ihdg,lenc,'H',ifilout,glob,istat)
      call qmaxn3 (glob,ihdg(1:14),ihdg(15:26),1,1,1,nx,my,1)

!=======================================================================
      return
      end

!***********************************************************************
      subroutine sigmap(layer,aki,bki,psfc,p)
      implicit none
      integer i, layer
      real aki(layer), bki(layer), psfc, p(layer)

      do i=1,layer
      p(i)=aki(i)+(bki(i)*psfc)
      enddo
      return
      end
