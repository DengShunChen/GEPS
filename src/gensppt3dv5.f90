subroutine gensppt3dv5(nx,my,my_max,lev,sppt3d                &
      ,sppt2d500,sppt2d1000,sppt2d2000                        &
      ,facsppt500,facsppt1000,facsppt2000                     &
      ,de_corretime_500,de_corretime_1000,de_corretime_2000   &
      ,dt,itimestep)
!c
!c purpose: To generate 3D SPPT strucutre
!c output: sppt3d
!c

!      use fftcom
      use mpe
      use rank, only : myrank
      use index
      implicit none
!c
!c generate SPPT 3d structure following ECMW technical report
!c
      integer nx,my,dimtotal,ndim,idum,ii,i,jj,j,k,ncx,nxj,lev,my_max
!c
      real dt
      real sppt2d(nx,my)  
      real sppt3d(nxp,lev,my_max)   ! for plot figure
      integer itimestep
!c
      real,intent(out) :: sppt2d500(nxp,my_max)              &
                         ,sppt2d1000(nxp,my_max)             &
                         ,sppt2d2000(nxp,my_max)
!c
      real,dimension(:,:),allocatable,save :: rold500,rold1000,rold2000
!      real rold500(mlmax,2),rold1000(mlmax,2),rold2000(mlmax,2)
      integer mlmax500,mlmax1000,mlmax2000
      integer jtrun500,jtrun1000,jtrun2000
      real de_corretime_500,de_corretime_1000,de_corretime_2000
      real facsppt500,facsppt1000,facsppt2000

!c if you need to plot sppt3d structure, then turn it on
!c      nxmy4=nx*my*4
!c      open(11,file='strucsppt3dv5.dat',access='direct'
!c     &    ,form='unformatted',recl=nxmy4,status='unknown')
!c
!c  initialize ifax and trigs for fft991 routine
!c
!c      call fftfax (nx,ifax,trigs)
!c
!c  Please make sure your definition in init_block.f90 
!c
!      ncx=128 ! for 500km perturbation
      ncx=2.*3.14159*6371./500. 
      jtrun500= 2*((1+(ncx-1)/3)/2) ! for 500km perturbation
      mlmax500= jtrun500*(jtrun500+1)/2
!c      de_corretime_500=6. ! unit: hrs
      jtrun1000= 2*((1+(ncx/2-1)/3)/2) ! for 1000km perturbation
      mlmax1000= jtrun1000*(jtrun1000+1)/2
!c      de_corretime_1000=3.*24. ! unit: hrs
      jtrun2000= 2*((1+(ncx/4-1)/3)/2) ! for 2000km perturbation
      mlmax2000= jtrun2000*(jtrun2000+1)/2
!c      de_corretime_2000=30.*24. ! unit: hrs
!c
     if ( .not. allocated(rold500) ) allocate(rold500(mlmax500,2))
     if ( .not. allocated(rold1000)) allocate(rold1000(mlmax1000,2))
     if ( .not. allocated(rold2000)) allocate(rold2000(mlmax2000,2))

      if (itimestep.eq.1) then
          rold500=0.
          rold1000=0.
          rold2000=0.
      endif
!c
!c 500km perturbation
!c
      if (myrank .eq. 0 ) &
      call gensppt2d(nx,my,sppt2d,rold500,mlmax500,jtrun500     &
                          ,de_corretime_500,dt,itimestep)
      call mpe_bcast(sppt2d,nx*my,0,mpe_double)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) then
          call reducepick (sppt2d(1,j),nxdef(j),nx,1)
        endif
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        do i=1,nxj
          sppt2d500(i,jj)=sppt2d(ii,j)
          ii=ii+1
        enddo
      enddo
!c
!c 1000km perturbation
!c
      if (myrank .eq. 0 ) &
      call gensppt2d(nx,my,sppt2d,rold1000,mlmax1000,jtrun1000 &
                          ,de_corretime_1000,dt,itimestep)
      call mpe_bcast(sppt2d,nx*my,0,mpe_double)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) then
          call reducepick (sppt2d(1,j),nxdef(j),nx,1)
        endif
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        do i=1,nxj
          sppt2d1000(i,jj)=sppt2d(ii,j)
          ii=ii+1
        enddo
      enddo
!c
!c 2000km perturbation
!c
      if (myrank .eq. 0 ) &
      call gensppt2d(nx,my,sppt2d,rold2000,mlmax2000,jtrun2000 &
                          ,de_corretime_2000,dt,itimestep)
      call mpe_bcast(sppt2d,nx*my,0,mpe_double)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) then
          call reducepick (sppt2d(1,j),nxdef(j),nx,1)
        endif
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        do i=1,nxj
          sppt2d2000(i,jj)=sppt2d(ii,j)
          ii=ii+1
        enddo
      enddo
!
      do jj=1,jlistnum
          j=jlist1(jj)
        do k=1,lev
          nxj=nxdef_2d(j)
          do i=1,nxj
          sppt3d(i,k,jj)=( sppt2d500(i,jj)*facsppt500     &
                         +sppt2d1000(i,jj)*facsppt1000    &
                         +sppt2d2000(i,jj)*facsppt2000 )  &
!     &                      *exp(-(k-40.)*(k-40.)/30.)
                            *exp(-(k-50.)*(k-50.)/8000.)
!c     &                     *exp(-(k-10.)*(k-10.)/600.)
          enddo
        enddo
      enddo
!
      do jj=1,jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
        do i=1,nxj
        sppt3d(i,1,jj)=sppt3d(i,1,jj)*0.2
        sppt3d(i,2,jj)=sppt3d(i,2,jj)*0.4
        sppt3d(i,3,jj)=sppt3d(i,3,jj)*0.6
        enddo
      enddo
!c
      if (myrank.eq.0) then
      print*,'itimestep, sppt3d(1,10,1)= ' &
              ,itimestep,sppt3d(1,10,1)
      endif
!c
!c could try
!c      facsppt500=0.15
!c      facsppt1000=0.08
!c      facsppt2000=0.045
!c suggested setup
!c      facsppt500=0.15
!c      facsppt1000=0.05
!c      facsppt2000=0.015
!c
      !TOP=1 Boundry=60
!
!c
!c for plot, turn it on
!c
!c      do k=1,lev
!c      do j=1,my
!c      do i=1,nx
!c      sppt2d(i,j)=sppt3d(i,k,j)
!c      enddo
!c      enddo
!c      write(11,rec=nrec) sppt2d 
!c      print*,'it, sppt2d(100,50)= ',it,nrec,sppt2d(100,50)
!c      nrec=nrec+1
!c      enddo
!c
      return
      end
!c
!c
      subroutine gensppt2d(nx,my,sppt2d,rold,mlmax,jtrun &
                          ,de_corretime,dt,itimestep)
!c
!c purpose: To generate 2D SPPT strucutre
!c output: sppt3d
!c
!c SPPT 2D 500km example
!c      parameter (ncx=128,mcy=ncx/2)
!c      parameter (jtrun= 2*((1+(ncx-1)/3)/2), mlmax= jtrun*(jtrun+1)/2)
!c SPPT 2D 1000km example
!c      parameter (ncx=64,mcy=ncx/2)
!c      parameter (jtrun= 2*((1+(ncx-1)/3)/2), mlmax= jtrun*(jtrun+1)/2)
!c SPPT 2D 500km example
!c      parameter (ncx=32,mcy=ncx/2)
!c      parameter (jtrun= 2*((1+(ncx-1)/3)/2), mlmax= jtrun*(jtrun+1)/2)
!c
      use fftcom
!c      use mpe
      use rank, only : myrank
!c      use index
!c      include '../include60/fftcom_sppt.h'

!c
!c generate SPPT 2d structure following ECMW technical report
!c
      integer nx,my,lev,dimtotal,ndim,idum,j,my2
      real sppt2d(nx,my),ran1,gasdev
      integer itimestep
      integer timearray(3),iseed
!c
!c
      integer nperr,mlmax,jtrun,rl,rm,rlm,ml
      integer ran2
      real scaleval(nx,my)  &
          ,ran_Renumber(mlmax),ran_Imnumber(mlmax)
      real aves_Re,vars_Re,stds_Re
      real aves_Im,vars_Im,stds_Im
      real aves,vars,stds
      real f0,rkT,correL,cp,capa,rgas,pi,radsq,rad,one,onem,irad,dt
!c    &    ,ran_number(nx*my)
!c
!c  Spheric Harmonic Constants
!c
      dimension msort(mlmax),lsort(mlmax),mlsort(jtrun,jtrun)
      dimension poly(mlmax,my/2),dpoly(mlmax,my/2),eps4(mlmax) &
      ,         cim(mlmax)
      real cosl(my),onocos(my)
      real weight(my),sinl(my)
!c
      real de_corretime
      real totalwavenumber
      real rold(mlmax,2),rnow(mlmax,2),rten(mlmax,2),sigma_n(mlmax)
!c=====================================================================
!c
!c  prepare spheric harmonic coefficients
!c
      cp=1004.24
      capa= 1.0/3.5
      rgas= capa*cp
      pi  = 4.0*atan(1.0)
      rad = 6.371e6
      radsq= rad*rad
!c
!c  build pointer arrays for locating zonal and total wavenumber
!c  values in the one-dimensional spherical harmonic arrays.
!c
      call sortml (jtrun,mlmax,msort,lsort,mlsort)
!c
      do 150 ml = 1, mlmax
      rl = lsort(ml)
      rm = msort(ml)-1
      rlm= rl-1.0
      if (msort(ml).eq.1)  rm= 0.0
      if (lsort(ml).eq.1)  rlm= 0.0
      eps4(ml)= rl*rlm/radsq
      cim(ml)= rm
  150 continue
!c
!c  gaussian quadrature weights and latitudes
!c
      one = 1.0
      onem= -one
      call gausl3 (my,onem,one,weight,sinl)
!c
      my2= my/2
!cdir$ ivdep
      do 180 j = 1, my2
      sinl(my+1-j)  = -sinl(j)
      weight(my+1-j)= weight(j)
      onocos(j)     = 1.0/(1.0-sinl(j)*sinl(j))
      onocos(my+1-j)= onocos(j)
      cosl(j)       = 1.0/sqrt(onocos(j))
      cosl(my+1-j)  = cosl(j)
!c      print *,'j= ',j,'sinl(j)= ',sinl(j)
  180 continue
!c
!c  define associated legendre polynomials and their derivatives
!c
      call lgndr_sppt (my2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
!c
      correL=10.0*1.e3    ! correlation lenght 10 km
      rkT=0.5*correL*correL/radsq
      if (myrank.eq.0) print *,'rkT= ',rkT
      totalwavenumber=0.
      do ml=1,mlmax
      totalwavenumber=totalwavenumber                   &
                     +(2*ml+1)*exp(-1.0*rkT*ml*(ml+1))
      enddo
      if (myrank.eq.0) then
      print *,'mlmax= ',mlmax
      print *,'totalwavenumber= ',totalwavenumber
      endif
!c
!c 60 hrs e-folding decay
!c      de_phi=exp(-1.*(itimestep-1)*dt/(120.*3600.))
!c      de_phi=exp(-1.*itimestep*dt/(120.*3600.))
!c      de_phi=exp(-1.*dt/(12.*3600.))  ! 2015-20180102 setup
      de_phi=exp(-1.*dt/(de_corretime*3600.))  ! 20180122 setup
!c      de_phi=exp(-1.*dt/(6.*3600.))
!c      de_phi=0.96
!c      print *,'de_phi= ',de_phi
!c
      dimtotal=mlmax
      nperr=21+3*(itimestep-1)*1331
!c      idum=-1+nperr*13241*(-1.)
!c      iRe_dum=mod(-1+nperr*13241*(-1.),1956516643)
!c      iIm_dum=mod(1+(nperr-131)*13241*1.,2147483647)
      call itime(timearray)
      iseed=irand(0)
      iRe_dum=irand(timearray(1)+(i*11467-iseed)*timearray(2) &
              +(iseed-i*23)*timearray(3))
      call itime(timearray)
      iseed=irand(0)
      iIm_dum=irand(timearray(1)+(i*11467-iseed)*timearray(2) &
              +(iseed-i*23)*timearray(3))
!c
!c      print *,'nperr, iRe_dum, iIm_dum= ',nperr, iRe_dum, iIm_dum
!c
      do ndim=1,dimtotal
        ran_Renumber(ndim)=gasdev(iRe_dum)*1.0
        ran_Imnumber(ndim)=gasdev(iIm_dum)*1.0
      enddo
!c
!c normalized the ran_number
!c
      call avevar_sppt(ran_Renumber,dimtotal,aves_Re,vars_Re,stds_Re)
      call avevar_sppt(ran_Imnumber,dimtotal,aves_Im,vars_Im,stds_Im)
!c
      do i=1,dimtotal
      ran_Renumber(i)=(ran_Renumber(i)-aves_Re)/stds_Re
      ran_Imnumber(i)=(ran_Imnumber(i)-aves_Im)/stds_Im
      enddo
      do ml=1,mlmax
      rnow(ml,1)=ran_Renumber(ml)
      rnow(ml,2)=ran_Imnumber(ml)
      enddo
!c
      call transr_sppt(jtrun,mlmax,nx,my,1,poly,rnow,sppt2d)
!c calculate vars on real grid
      call avevar_sppt2d(sppt2d,nx,my,aves,vars,stds)
!c
!c set initial random value
!c
      if (itimestep.eq.1) then
!c
!c set the amplitude of noise
!c
      f0=sqrt( vars*(1-de_phi*de_phi)/(2.*totalwavenumber) )
      do ml=1,mlmax
      sigma_n(ml)=f0*exp(-rkT*ml*(ml+1)/2.)
      rold(ml,1)=sigma_n(ml)*rnow(ml,1)/(sqrt(1.-de_phi*de_phi))
      rold(ml,2)=sigma_n(ml)*rnow(ml,2)/(sqrt(1.-de_phi*de_phi))
      enddo
      call transr_sppt(jtrun,mlmax,nx,my,1,poly,rold,sppt2d)
!c
      call avevar_sppt2d(sppt2d,nx,my,aves,vars,stds)
!c      call removegt3std(sppt2d,nx,my,aves,stds)    !=========20180306

      else
!c
!c set the amplitude of noise
!c
      f0=sqrt( vars*(1-de_phi*de_phi)/(2.*totalwavenumber) )
!c      print *,'itimestep, f0= ',itimestep,f0
!c
      do ml=1,mlmax
      sigma_n(ml)=f0*exp(-rkT*ml*(ml+1)/2.)
      rten(ml,1)=de_phi*rold(ml,1)+sigma_n(ml)*rnow(ml,1)
      rten(ml,2)=de_phi*rold(ml,2)+sigma_n(ml)*rnow(ml,2)
      enddo
!c      print *,'sigma_n(10)*rnow(10,1),de_phi*rold(10,1)=',
!c     &         sigma_n(10)*rnow(10,1),de_phi*rold(10,1)
!c
      call transr_sppt(jtrun,mlmax,nx,my,1,poly,rten,sppt2d)
!c
      call avevar_sppt2d(sppt2d,nx,my,aves,vars,stds)
!c      call removegt3std(sppt2d,nx,my,aves,stds)   !=======20180306
!c
      do ml=1,mlmax
      rold(ml,1)=rten(ml,1)
      rold(ml,2)=rten(ml,2)
      enddo
      endif  ! end of different time step set up
!c
      return
      end
!c
!c
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
!CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
   11   v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 11
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
!c
!c
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
      NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
!c
!c
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,    &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!c
!c
      SUBROUTINE avevar_sppt(data,n,ave,var,std)
      INTEGER n
      REAL ave,var,data(n)
      INTEGER j
      REAL s,ep,std
      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
!c      print *,'s= ',s
!c      print *,'ep= ',ep
12    continue
!c      var=(var-ep**2/n)/(n-1)    ! sample numners < 30
      var=(var-ep**2/n)/n
      std=sqrt(var)
      return
      end
!c
!c
      SUBROUTINE avevar_sppt2d(data2d,n,m,ave,var,std)
      INTEGER n,m,nmdim
      REAL ave,var,data2d(n,m),data(n*m)
      INTEGER j
      REAL s,ep,std
!c
      nmdim=0
      do j=1,m
      do i=1,n
      nmdim=nmdim+1
      data(nmdim)=data2d(i,j)
      enddo
      enddo
!c
      
      ave=0.0
      do 11 j=1,nmdim
        ave=ave+data(j)
11    continue
      ave=ave/nmdim
      var=0.0
      ep=0.0
      do 12 j=1,nmdim
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
!c      print *,'s= ',s
!c      print *,'ep= ',ep
12    continue
!c      var=(var-ep**2/n)/(n-1)    ! sample numners < 30
      var=(var-ep**2/nmdim)/nmdim
      std=sqrt(var)
      return
      end
!c
!c
      subroutine removegt2std(scaleval,nx,my,aves,stds)
      integer nx,my
      real scaleval(nx,my),aves,stds
      real vcheck2p,vcheck2n
!c
      vcheck2p=2.0*stds
      vcheck2n=-2.0*stds
!c
!c      print *,'aves= ',aves
!c
      do j=1,my
      do i=1,nx
        if ((scaleval(i,j)-aves).gt.vcheck2p) scaleval(i,j)=vcheck2p
        if ((scaleval(i,j)-aves).lt.vcheck2n) scaleval(i,j)=vcheck2n 
!c        if ((scaleval(i,j)-aves).gt.vcheck2p) then
!c          scaleval(i,j)=vcheck2p
!c          print *,'remove overflow value'
!c        endif
!c        if ((scaleval(i,j)-aves).lt.vcheck2n) then
!c          scaleval(i,j)=vcheck2n
!c          print *,'remove underflow value'
!c        endif
      enddo
      enddo
!c
      return
      end
!c
!c
      subroutine removegt3std(scaleval,nx,my,aves,stds)
      integer nx,my
      real scaleval(nx,my),aves,stds
      real vcheck3p,vcheck3n
!c
      vcheck3p=3.0*stds
      vcheck3n=-3.0*stds
!c
!c      print *,'aves= ',aves
!c
      do j=1,my
      do i=1,nx
        if ((scaleval(i,j)-aves).gt.vcheck3p) scaleval(i,j)=vcheck3p
        if ((scaleval(i,j)-aves).lt.vcheck3n) scaleval(i,j)=vcheck3n 
!c        if ((scaleval(i,j)-aves).gt.vcheck3p) then
!c          scaleval(i,j)=vcheck3p
!c          print *,'remove overflow value'
!c        endif
!c        if ((scaleval(i,j)-aves).lt.vcheck3n) then
!c          scaleval(i,j)=vcheck3n
!c          print *,'remove underflow value'
!c        endif
      enddo
      enddo
!c
      return
      end
!c
!c
      subroutine lgndr_sppt (my2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
!c
!c  generate legendre polynomials and their derivatives on the
!c  gaussian latitudes
!c
!c ***input***
!c
!c  my2:  number of gaussian latitudes from south pole and equator
!c  jtrun:  zonal wavenumber truncation limit
!c  mlmax: total number of triangular truncation spherical harmonics
!c  mlsort: pointer array of 1-d indexs at functions of zonal and
!c          total wavenumbers
!c  sinl: sin of gaussian latitudes
!c
!c  ***output***
!c
!c  poly: associated legendre coefficients
!c  dpoly: d(poly)/d(sinl)
!c
!c ******************************************************************
!c
!c ref= belousov, s. l., 1962= tables of normalized associated
!c        legendre polynomials. pergamon press, new york
!c
      dimension poly(mlmax,my2),dpoly(mlmax,my2),sinl(my2) &
      , mlsort(jtrun,jtrun)
!c
!c      parameter (jtrunx= 100)
      dimension pnm(jtrun+1,jtrun+1),dpnm(jtrun+1,jtrun+1)
!c
!c sinl is sin(latitude) = cos(colatitude)
!c pnm(np,mp) is legendre polynomial p(n,m) with np=n+1, mp=m+1
!c pnm(mp,np+1) is x derivative of p(n,m) with np=n+1, mp=m+1
!c
      jtrunp= jtrun+1
      do 1001 j=1,my2
      xx= sinl(j)
      sn= sqrt(1.0-xx*xx)
	sn2i = 1.0/(1.0 - xx*xx)
      rt2= sqrt(2.0)
	c1 = rt2
!c
	pnm(1,1) = 1.0/rt2
      theta=-atan(xx/sqrt(1.0-xx*xx))+2.0*atan(1.0)
!c
      do 20 n=1,jtrun
	np = n + 1
      fn=n
	fn2 = fn + fn
	fn2s = fn2*fn2
!c eq 22
      c1= c1*sqrt(1.0-1.0/fn2s)
      c3= c1/sqrt(fn*(fn+1.0))
	ang = fn*theta
	s1 = 0.0
	s2 = 0.0
	c4 = 1.0
	c5 = fn
	a = -1.0
	b = 0.0
!c
      do 27 kp=1,np,2
	k = kp - 1
      s2= s2+c5*sin(ang)*c4
      if (k.eq.n) c4 = 0.5*c4
      s1= s1+c4*cos(ang)
	a = a + 2.0
	b = b + 1.0
      fk=k
	ang = theta*(fn - fk - 2.0)
	c4 = (a*(fn - b + 1.0)/(b*(fn2 - a)))*c4
	c5 = c5 - 2.0
   27 continue
!c eq 19
	pnm(np,1) = s1*c1
!c eq 21
	pnm(np,2) = s2*c3
   20 continue
!c
      do 4 mp=3,jtrunp
	m = mp - 1
      fm= m
	fm1 = fm - 1.0
	fm2 = fm - 2.0
	fm3 = fm - 3.0
      c6= sqrt(1.0+1.0/(fm+fm))
!c eq 23
	pnm(mp,mp) = c6*sn*pnm(m,m)
      if (mp - jtrunp) 3,4,4
    3 continue
	nps = mp + 1
!c
      do 41 np=nps,jtrunp
	n = np - 1
      fn= n
	fn2 = fn + fn
	c7 = (fn2 + 1.0)/(fn2 - 1.0)
	c8 = (fm1 + fn)/((fm + fn)*(fm2 + fn))
      c= sqrt((fn2+1.0)*c8*(fm3+fn)/(fn2-3.0))
      d= -sqrt(c7*c8*(fn-fm1))
      e= sqrt(c7*(fn-fm)/(fn+fm))
!c eq 17
	pnm(np,mp) = c*pnm(np-2,mp-2) &
                 + xx*(d*pnm(np-1,mp-2) + e*pnm(np - 1,mp))
   41 continue
    4 continue
!c
      do 50 mp=1,jtrun
      fm= mp-1.0
	fms = fm*fm
      do 50 np=mp,jtrun
      fnp= np
	fnp2 = fnp + fnp
	cf = (fnp*fnp - fms)*(fnp2 - 1.0)/(fnp2 + 1.0)
      cf= sqrt(cf)
!c der
      dpnm(np,mp)   = -sn2i*(cf*pnm(np+1,mp) - fnp*xx*pnm(np,mp))
   50 continue
!c
      do 71 m=1,jtrun
      do 71 l=m,jtrun
      ml= mlsort(m,l)
      poly(ml,j)= pnm(l,m)
      dpoly(ml,j)=dpnm(l,m)
   71 continue
      dpoly(1,j)= 0.0
 1001 continue
      return
      end
!c
!c
      subroutine transr_sppt (jtrun,mlmax,nx,my,ll,poly,s,r)
!c
!c  subroutine to transform a spectral coefficient field to
!c  grid point form
!c
!c *** const ***
!c
!c  jtrun: zonal wavenumber resolution limit
!c  mlmax: number of spectral coefficients (horizontal field)
!c  nx: e-w dimension no.
!c  my: n-s dimension no.
!c  ll: number of levels to transform
!c  poly: legendre polynomials
!c
!c *** input variable ***
!c
!c  s: spectral coefficient array to transform
!c
!c *** output variable ***
!c
!c  r: 3-d output grid point fields
!c
!c  **************************************
!c
      use fftcom 
!c
      dimension poly(mlmax,my/2),s(mlmax,2,ll),r(nx,ll,my)
!c      include '../include/fftcom.h'
!c      common/fft/ trigs(4096),ifax(19)
!csun  include '../include/paramt.h' .. change im,jm to nx,my
      dimension cc(nx+2,my),work(nx*my,2)
      integer mlmax,my,nx,ll,jtrun
      integer mlx,i,k,m,j,jj,ml,mm,mp
!c
      mlx= (jtrun/2)*((jtrun+1)/2)
      do 20 k=1,ll
      do 55 m=1,(nx+2)*my/2
      cc(m,1)= 0.0
      cc(m,my/2+1)= 0.0
   55 continue
!c
      do 5 j=1,my/2
      jj= my+1-j
      ml= 2*mlx
!cdir@ ivdep
      do 3 m=2,jtrun,2
      ml= ml+1
      mm= 2*m-1
      mp= mm+1
      cc(mm,j)= poly(ml,j)*s(ml,1,k)
      cc(mp,j)= poly(ml,j)*s(ml,2,k)
      cc(mm,jj)= cc(mm,j)
      cc(mp,jj)= cc(mp,j)
    3 continue
!c
      m1= 0
      do 5 l=jtrun-1,1,-2
!cdir@ ivdep
      do 6 m=1,l
      mm= 2*m-1
      mp= mm+1
      ml= m+m1
      mk= ml+mlx
      cc(mm,j)= cc(mm,j)+poly(ml,j)*s(ml,1,k)+poly(mk,j)*s(mk,1,k)
      cc(mm,jj)=cc(mm,jj)+poly(ml,j)*s(ml,1,k)-poly(mk,j)*s(mk,1,k)
      cc(mp,j)= cc(mp,j)+poly(ml,j)*s(ml,2,k)+poly(mk,j)*s(mk,2,k)
      cc(mp,jj)=cc(mp,jj)+poly(ml,j)*s(ml,2,k)-poly(mk,j)*s(mk,2,k)
    6 continue
      m1= m1+l
    5 continue
!c
!c      call fft991(cc,work,trigs,ifax,1,nx+3,nx,my,1)
      call rfftmlt(cc,work,trigs,ifax,1,nx+2,nx,my,1)
!c
      do 22 j=1,my
      do 22 i=1,nx
      r(i,k,j)= cc(i,j)
   22 continue
!c
   20 continue
!c
      return
      end
