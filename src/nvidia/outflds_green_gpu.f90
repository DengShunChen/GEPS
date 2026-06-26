      subroutine outflds_green_gpu( itau,nx,my,my_max,lev,ncld   &
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
      use mod_grb2_param , only :ofdir,wrtgrb2_v2_gpu,wrtgrb2_accu_v2_gpu

      use noah,only:runoff
      use mod_qsatq

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
      real      whtlev(100),whtlevq(100),whtlevz(100)
      character*6 labx

      integer   jj,j,nxj,k,i,n,nk,ntrac,ll,mm,la,kk
!
      real(kind=RTYPE) oqt(nxp,my_max),oqc(nxp,my_max), ou(nxp,my_max),  &
                        ov(nxp,my_max), ot(nxp,my_max),pla(nxp,my_max)
      real pp(nxp,my_max),p2(nxp,my_max),p10(nxp,my_max),   &
           rhtmp(nxp,my_max), plt_bt(nxp,my_max)
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
      real plyr(nxp,l,my_max)
      integer,parameter:: async_id = 1
!
!$acc wait(async_id)
!$acc enter data create(pla,oqt,oqc,ou,ov,ot,tepl,plyr) async(async_id)
!$acc enter data create(pp,p2,p10,tht) async(async_id)
!$acc enter data create(wrk,glob,mout) async(async_id)
!$acc enter data create(akir,bkir) async(async_id)
!$acc enter data create(plt_bt,rhtmp) async(async_id)
!$acc enter data copyin(hm) async(async_id)

!$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
           pla(i,jj)= 0.0
           oqt(i,jj)= 0.0
           oqc(i,jj)= 0.0
           ou (i,jj)= 0.0
           ov (i,jj)= 0.0
           ot (i,jj)= 0.0
         endif
        enddo
      enddo

      rcp=rgas/cp
      lenc = nx*my
      nc=0

!$acc parallel loop private(kk) async(async_id)
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
      call sigmap_gpu(nxp,my_max,l,akir,bkir,pt,plyr)
      do mm=1,m    ! 1: 100m, 2: 40m.
!-----------------------------------------------------------------------
! calculate PQUVT at new layer.
!$acc parallel loop collapse(2) private(j,nxj,avett,la,temp,xxx) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
! ave. temp.
          avett=(tt(i,lev,jj)+tt(i,lev-1,jj))/2.*                  &
                ((plt(i,lev,jj)+plt(i,lev-1,jj))/2./1000.)**rcp/     &
                (1.+0.608*(qt(i,lev-1,jj)+qt(i,lev-1,jj))/2.)
! calculate 40m&100m P, 2m P, 10m P.
        pp(i,jj)=(pt(i,jj)+ptop)*(exp(-(hm(mm)*grav)/(rgas*avett)))
        p2(i,jj)=(pt(i,jj)+ptop)*(exp(-(2.0*grav)/(rgas*avett)))
        p10(i,jj)=(pt(i,jj)+ptop)*(exp(-(10.0*grav)/(rgas*avett)))
        pla(i,jj)=pp(i,jj) * 100.0
! transfer temperature from tt(virtual theta) to tepl(temperature).
        la=0
        do ll=lev,lev-3,-1
          la=la+1
          tepl(i,la,jj)=(tt(i,ll,jj)*(pp(i,jj)/1000.)**rcp)/(1.+0.608*qt(i,ll,jj))
        enddo

        !call sigmap(l,akir(:),bkir(:),pt(i,jj),p)
! find pp pressure position
         if(pp(i,jj).gt.plyr(i,1,jj))then
          temp=(pp(i,jj)-plyr(i,2,jj))/(plyr(i,1,jj)-plyr(i,2,jj))-1.0
          oqt(i,jj)=qt(i,lev,jj)+temp*(qt(i,lev,jj)-qt(i,lev-1,jj))
          oqt(i,jj)=max(oqt(i,jj),1.e-8)
          oqc(i,jj)=qt(i,lev*2,jj)+temp*(qt(i,lev*2,jj)-qt(i,lev*2-1,jj))
          oqc(i,jj)=max(oqc(i,jj),1.e-8)

          xxx= rad/cosl(j)
          ou(i,jj)=u10(i,jj)/xxx
          ov(i,jj)=v10(i,jj)/xxx
          temp=(plyr(i,1,jj)-pp(i,jj))/(plyr(i,1,jj)-p10(i,jj))
          ou(i,jj)=ou(i,jj)*temp+ut(i,lev,jj)*(1.0-temp)
          ov(i,jj)=ov(i,jj)*temp+vt(i,lev,jj)*(1.0-temp)

          temp=(plyr(i,1,jj)-pp(i,jj))/(plyr(i,1,jj)-p2(i,jj))
          ot(i,jj)=t2(i,jj)*temp+tepl(i,1,jj)*(1.0-temp)

          goto 99
          endif

        do ll=1,l-1
         if((pp(i,jj).le.plyr(i,ll,jj)).and.(pp(i,jj)).gt.plyr(i,ll+1,jj))then
          temp=(plyr(i,ll,jj)-pp(i,jj))/(plyr(i,ll,jj)-plyr(i,ll+1,jj))
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

         if(pp(i,jj).le.plyr(i,l,jj))then
          temp=(plyr(i,l-1,jj)-pp(i,jj))/(plyr(i,l-1,jj)-plyr(i,l,jj))-1.0
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
        endif
        enddo  ! end (i)
      enddo  ! end (jj)

!$acc parallel loop collapse(2) private(j,nxj,xxx) async(async_id)
      do 50 jj = 1, jlistnum
      do 50 i=1,nxp
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      if(i<=nxj)then
      xxx= rad/cosl(j)
       ou(i,jj)= ou(i,jj)*xxx
       ov(i,jj)= ov(i,jj)*xxx
      endif
   50 continue

!-----------------------------------------------------------------------
!output P
      write(wtemp,'(a3,a3)')layer(mm),var(1)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,pla,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,3,0,1,103,0,nint(hm(mm)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output Q
      write(wtemp,'(a3,a3)')layer(mm),var(2)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,oqt,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,0,6,103,0,nint(hm(mm)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      write(wtemp,'(a3,a3)')layer(mm),var(6)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,oqc,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,235,8,103,0,nint(hm(mm)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output U,V
      write(wtemp,'(a3,a3)')layer(mm),var(3)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,ou,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,2,2,2,103,0,nint(hm(mm)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      write(wtemp,'(a3,a3)')layer(mm),var(4)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,ov,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,2,3,2,103,0,nint(hm(mm)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output T
      write(wtemp,'(a3,a3)')layer(mm),var(5)
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,ot,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,0,0,2,103,0,nint(hm(mm)),-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!-----------------------------------------------------------------------
      enddo  ! end (mm)
!=======================================================================
!output S00310(net SW flux at the surface)
      write(wtemp,'(a6)')'S00310'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) wrk(i,jj)=ss(i,jj)
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,4,9,2,1,0,0,-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output RH at bottom level
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          if(i<=nxj)then
           tht(i,jj) =tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
           plt_bt(i,jj) = plt(i,lev,jj)
          endif
        enddo
      enddo
      call qsatq_gpu (nx,my_max,tht,plt_bt,rhtmp)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          if(i<=nxj)then
           wrk(i,jj)=100.*(qt(i,lev,jj)/rhtmp(i,jj))
           wrk(i,jj)= min( 100., max( 1., wrk(i,jj) ) )
          endif
        enddo
      enddo
      write(wtemp,'(a6)')'B00510'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,1,2,103,0,0,-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!output b00010
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i = 1, nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) wrk(i,jj)=( pt(i,jj)+ptop ) * 100.0
        enddo
      enddo
      write(wtemp,'(a6)')'B00010'
      call syslbl_w(wtemp,idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,3,0,1,103,0,0,-999,-999/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)
!=======================================================================
!output 6hr prec.
      if (mod(float(itau)+0.00001, 6. ) .lt. 0.01) then
      call syslbl_w ('b00633',idtg,itau,ggdef)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) wrk(i,jj)=raincu6(i,jj)
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,10,2,103,0,0,1,6/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!
      call syslbl_w ('b00643',idtg,itau,ggdef)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) wrk(i,jj)=rainlp6(i,jj)
        enddo
      enddo
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,47,2,103,0,0,1,6/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

!
      call syslbl_w ('b00623',idtg,itau,ggdef)
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do 98 jj = 1, jlistnum
      do 98 i=1,nxp
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      if (i<=nxj) wrk(i,jj)=raincu6(i,jj)+rainlp6(i,jj)
 98   continue
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/0,1,8,2,103,0,0,1,6/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      !runoff
      !$acc parallel loop collapse(2) private(j,nxj) async(async_id)
      do jj = 1, jlistnum
        do i=1,nxp
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        if(i<=nxj) wrk(i,jj)=runoff(i,jj)
        enddo
      enddo
      call syslbl_w ('b00663',idtg,itau,ggdef)
      call unify_reduceintp_gpu(nx,my,my_max,wrk,glob)
      call qmaxn3_w (glob,1,1,1,nx,my,1)
      ptp0=(/2,0,5,2,103,0,0,1,6/)
      call split_v2_gpu(nx,my,lenc,nc,glob,mout,ptp0,ptp1)

      endif !mod(float(itau)+0.00001, 6. ) .lt. 0.01

!

!
      if(outdms.gt.0)then
      if ( myrank .lt. nc ) then
         !$acc update self(mout) async(async_id)
         !$acc wait(async_id)
         call dmswrit_split(nx,my,lenc,kflag,mout,istat)
      endif
      endif

      if(outgrb2 == 1 )then
       if ( myrank .lt. nc ) then
        if(ptp1(8)==-999)then
        call wrtgrb2_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),mout)
        else
        call wrtgrb2_accu_v2_gpu(itau,ptp1(1),ptp1(2),ptp1(3),ptp1(4),ptp1(5) &
            ,ptp1(6),ptp1(7),ptp1(8),ptp1(9),mout)
        endif
       endif
      endif
!$acc wait(async_id)
!$acc exit data delete(wrk,glob,mout) async(async_id)
!$acc exit data delete(pla,oqt,oqc,ou,ov,ot,tepl,plyr) async(async_id)
!$acc exit data delete(pp,p2,p10,tht) async(async_id)
!$acc exit data delete(akir,bkir,hm) async(async_id)
!$acc exit data delete(plt_bt,rhtmp) async(async_id)
!$acc wait(async_id)
!
!=======================================================================
      return
      end

!***********************************************************************
      subroutine sigmap_gpu(nxp,my_max,layer,aki,bki,psfc,p)
      use const, only: RTYPE
      use index,only:jlistnum,nxdef_2d,jlist1
      implicit none
      integer i,j,jj,k,nxj,nxp,my_max, layer
      real p(nxp,layer,my_max)
      real(kind=RTYPE) psfc(nxp,my_max),aki(layer), bki(layer)
      integer,parameter:: async_id = 1
      !$acc parallel loop collapse(3) private(j,nxj) async(async_id)
      do  jj = 1, jlistnum
       do k=1,layer
        do i=1,nxp
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         if(i<=nxj)then
          p(i,k,jj)=aki(k)+(bki(k)*psfc(i,jj))
         endif
        enddo
       enddo
      enddo

      return
      end

