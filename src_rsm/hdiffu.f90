      subroutine hdiffu ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld,amp   &
                        , rad,cosl,ut,vt,vornow,divnow,temnow,eps4     &
                        , trefs)
      use index
      use mpe
      use rank
      use const, only : hdk1,hdk2,radsq,hord,vd,RTYPE
      use param, only : octahedral,af

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      rad
      real      dta

      real(kind=RTYPE) vornow(levp,2,jtrun,jtmax),divnow(levp,2,jtrun,jtmax),  &
                       temnow(levp,2,jtrun,jtmax),trefs(levp,2,jtrun,jtmax),   &
                       ut(nxp,lev,my_max),vt(nxp,lev,my_max),eps4(jtrun,jtmax),&
                       cosl(my)
!
!     parameter ( ktop=4, ktop2=ktop/2 ) ! top "ktop" levels are inhenced
!
      real      wmax(lev)
      real      windmax1,windmax2,windmax3

      integer   jj,j,nxj,k,i,m,n,mf,nc,kk,KL
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu,finc
      real      hfilt,hfilt2,nf,dec,coefu,factop,powd,kfac
      real      c1,c2,c3
      logical   windchk

      data      windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!

!
      nf=jtrun-1
      wmax(1:lev)= 0.0
!
      do jj =1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        xx=rad/cosl(j)
        do k=1,lev
          do i=1,nxj
            wmax(k)= max(wmax(k),xx*sqrt(ut(i,k,jj)**2+vt(i,k,jj)**2))
          enddo
        enddo
      enddo
!
      call mpe_global_max(wmax,lev,mpe_double)
!
!!      do k=1,lev
!!        if( wmax(k) .gt. windmax3 ) then
!!          if(myrank.eq.0)print *,'wmax gt windmax at','k= ',k,         &
!!                                 ' windmax=',wmax(k)
!!        endif
!!      enddo
!
!
      finc = 10.*max(af,0.001)
      powd = float(hord) / 2.
      hfilt  = (radsq/(nf*(nf+1)))**powd
      hfilt2 = radsq/(nf*(nf+1))
      factop = 1.5
      coefu = factop/float(hdk2(3)-hdk2(2))
      if ( octahedral ) then
        hfilt  = hfilt/(6.*dta)
        hfilt2 = hfilt2/(6.*dta)
      else
        hfilt  = hfilt/dta
        hfilt2 = hfilt2/dta
      endif

      do 100 k=1,levp  ! levp -> lev
!
!  compute diffusion coefficients
!
        KL=Llist(k)
!
        kfac  = (1.+finc)*max(float(hdk2(1)-KL),0.)
        facd = max(amp,kfac)
        facv = max(min(amp,1.),kfac)
        fact = max(min(amp,1.),kfac)
!!        facd = 1. * amp * (kfac + 2.*max(float(hdk1-KL),0.))
!!        facv = 1. * (kfac + 1.*max(float(hdk1-KL),0.))
!!        fact = 1. * (kfac + 1.*max(float(hdk1-KL),0.))
!
!  difuse vorticity and divergence fields
!  diffuse moisture and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun

            c1=1.+dta*facv*hfilt*eps4(n,m)**powd
            c3=1.+dta*fact*hfilt*eps4(n,m)**powd

            if ( KL .le. hdk1 ) then
              c2=1.+dta*facd*hfilt2*eps4(n,m)+vd*exp(-0.5*(k-1))
            else
              c2=1.+dta*facd*hfilt*eps4(n,m)**powd
            endif

            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
!!            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
!!            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
            temnow(k,1,n,m)=temnow(k,1,n,m)/c3
            temnow(k,2,n,m)=temnow(k,2,n,m)/c3
          enddo
        enddo
 100  continue
!
!-------------------------------------------------------------------
!
!  filter layers near the upper bound if too strong wind speed 
!  happens at the top layer
!
!  2003/10/7 :
!  sometimes wind speed greater than 100m/s happens at k=2
!
      windchk=.false.
      do k=1,8
        if ( wmax(k) .gt. windmax3 ) windchk=.true.
      enddo
      if ( windchk ) then
       call filter_top(jtrun,jtmax,levp,hdk1,ncld,temnow,vornow,divnow)
      end if
!--------------------------------------------------------------------
      return
      end
!
!--------------------------------------------------------------------
      subroutine ohdiffu ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld,amp   &
                        , rad,cosl,ut,vt,vornow,divnow,temnow,eps4     &
                        , trefs)
!     horizontal diffusion for Cubic Truncation Octahedral Reduced
!     Gaussain grid (under testing for now) -- 03/03/2020 Pang-Yang Liu 
      use index
      use mpe
      use rank
      use const, only : hdk1,hdk2,radsq,RTYPE
      use param, only : octahedral,af

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      dta,rad

      real(kind=RTYPE) vornow(levp,2,jtrun,jtmax),divnow(levp,2,jtrun,jtmax),  &
                       temnow(levp,2,jtrun,jtmax),trefs(levp,2,jtrun,jtmax),   &
                       ut(nxp,lev,my_max),vt(nxp,lev,my_max),                  &
                       eps4(jtrun,jtmax),cosl(my)
!
!     parameter ( ktop=4, ktop2=ktop/2 ) ! top "ktop" levels are inhenced
!
      real      wmax(lev)
      real      windmax1,windmax2,windmax3

      integer   jj,j,nxj,k,i,m,n,mf,nc,kk,KL
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu
      real      hfilt,hfilt2,nf,ncut,ncor,kfac,finc
      real      c1,c2,c3
      logical   windchk

      data      windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!

!
      finc = 10.*max(af,0.001)
      nf=jtrun-1
      ncut=0.5*jtrun
      wmax(1:lev)= 0.0
!
      do jj =1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        xx=rad/cosl(j)
        do k=1,lev
          do i=1,nxj
            wmax(k)= max(wmax(k),xx*sqrt(ut(i,k,jj)**2+vt(i,k,jj)**2))
          enddo
        enddo
      enddo
!
      call mpe_global_max(wmax,lev,mpe_double)
!
      do k=1,lev
        if( wmax(k) .gt. windmax3 ) then
          if(myrank.eq.0)print *,'wmax gt windmax at','k= ',k,         &
                                 ' windmax=',wmax(k)
        endif
      enddo
!
!
      hfilt = (radsq/(nf*(nf+1)))**2.
      hfilt = hfilt/(6.*dta)
      hfilt2 = radsq/(nf*(nf+1))
      hfilt2 = hfilt2/(6.*dta)


      do 100 k=1,levp  ! levp -> lev
!
!  compute diffusion coefficients
!

         KL=Llist(k)
!
         kfac  = min((5.+finc)*max(float(hdk2(1)-KL),0.),60.+5.*finc)
         facd = amp + kfac
         facv = amp + kfac
         fact = amp + kfac
!         facd= 1.0 + max(float(hdk2(1)-KL),0.)
!         facv= 1.0 + max(float(hdk2(1)-KL),0.)
!         fact= 1.0 + max(float(hdk2(1)-KL),0.)

!         if ( KL .le. hdk1 ) facd=facd*(1.+float(hdk1-KL))
!
!          facd = amp * facd
!          facv = amp * facv
!          fact = amp * fact
!
          
!
!  difuse vorticity and divergence fields
!  diffuse moisture and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            ncor=1.
            if ( n .gt. ncut ) then
              ncor=exp(-0.5*( (n-jtrun)**2. / (n-ncut)**2. ))
            else
              ncor=0.
            endif
            c3=1.+dta*fact*hfilt*ncor*eps4(n,m)**2

            ncor=min(ncor+0.2*max(float(hdk2(1)-KL),0.),1.)
            c1=1.+dta*facv*hfilt*ncor*eps4(n,m)**2
            if ( KL .le. hdk1 ) then
!              c2=1.+dta*facd*hfilt2*ncor*eps4(n,m)
              c2=1.+dta*facd*hfilt2*eps4(n,m)
            else
              c2=1.+dta*facd*hfilt*ncor*eps4(n,m)**2
            endif

            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
!            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
!            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
            temnow(k,1,n,m)=temnow(k,1,n,m)/c3
            temnow(k,2,n,m)=temnow(k,2,n,m)/c3
          enddo
        enddo
 100  continue
!
!-------------------------------------------------------------------
!
!  filter layers near the upper bound if too strong wind speed 
!  happens at the top layer
!
!  2003/10/7 :
!  sometimes wind speed greater than 100m/s happens at k=2
!
      windchk=.false.
      do k=1,8
        if ( wmax(k) .gt. windmax3 ) windchk=.true.
      enddo
!
      if ( windchk ) then
       call filter_top(jtrun,jtmax,levp,hdk1,ncld,temnow,vornow,divnow)
      end if
!--------------------------------------------------------------------
      return
      end
!
!--------------------------------------------------------------------
      subroutine whdiffu ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld,amp   &
                        , rad,cosl,ut,vt,vornow,divnow,eps4,trefs) 
      use index
      use mpe
      use rank
      use const, only : hdk1,hdk2,radsq,vd,RTYPE
      use param, only : octahedral,af,mwhd

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      dta,rad

      real(kind=RTYPE) vornow(levp,2,jtrun,jtmax),divnow(levp,2,jtrun,jtmax),  &
                       temnow(levp,2,jtrun,jtmax),trefs(levp,2,jtrun,jtmax),   &
                       ut(nxp,lev,my_max),vt(nxp,lev,my_max),                  &
                       eps4(jtrun,jtmax),cosl(my)
!
!     parameter ( ktop=4, ktop2=ktop/2 ) ! top "ktop" levels are inhenced
!
      real      wmax(lev)
      real      windmax1,windmax2,windmax3

      integer   jj,j,nxj,k,i,m,n,mf,nc,kk,KL
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu
      real      hfilt2,hfilt4,hfilt6,nf,kfac,finc,fl
      real      c1,c2,c3,c4
      logical   windchk

      data      windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!
      temnow = 0.0
      wmax(1:lev)= 0.0
!
      do jj =1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        xx=rad/cosl(j)
        do k=1,lev
          do i=1,nxj
            wmax(k)= max(wmax(k),xx*sqrt(ut(i,k,jj)**2+vt(i,k,jj)**2))
          enddo
        enddo
      enddo
!
      call mpe_global_max(wmax,lev,mpe_double)
!
      nf=jtrun-1
!
      finc = 10.*max(af,0.001)
      fl   = 150./hdk2(2)-1.
      hfilt6 = (radsq/(nf*(nf+1)))**3.
      hfilt4 = (radsq/(nf*(nf+1)))**2.
      hfilt2 = radsq/(nf*(nf+1))
      if ( octahedral ) then
        hfilt6 = hfilt6/(6.*dta)
        hfilt4 = hfilt4/(6.*dta)
        hfilt2 = hfilt2/(6.*dta)
      else
        hfilt6 = hfilt6/dta
        hfilt4 = hfilt4/dta
        hfilt2 = hfilt2/dta
      endif

      do 100 k=1,levp  ! levp -> lev
!
!  compute diffusion coefficients
!

        KL=Llist(k)
!
        kfac = (fl+finc)*max(float(hdk2(2)-KL),0.)
        facd = mwhd * max(amp,kfac)
        facv = max(min(amp,1.),kfac)
!!        fact = amp * kfacv 
!          endif

!
!  difuse vorticity and divergence fields
!  diffuse moisture and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun

!!            c1=1.+dta*facv*hfilt6*eps4(n,m)**3.
            c1=1.+dta*facv*hfilt4*eps4(n,m)**2.
!!            c2=1.+dta*facd*hfilt6*eps4(n,m)**3.
!!            c2=1.+dta*facd*hfilt2*eps4(n,m)

            if ( KL .le. hdk1 ) then
              c2=1.+dta*facd*hfilt2*eps4(n,m)+vd*exp(-0.5*(k-1))
            else
              c2=1.+dta*facd*hfilt4*eps4(n,m)**2.
            endif

!!            c3=1.+dta*fact*hfilt6*eps4(n,m)**3.
!!            c3=1.+dta*fact*hfilt4*eps4(n,m)**2.


            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
!!            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
!!            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
          enddo
        enddo
 100  continue
!!
!
      windchk=.false.
      do k=1,8
        if ( wmax(k) .gt. windmax3 ) windchk=.true.
      enddo
      if ( windchk ) &
       call filter_top(jtrun,jtmax,levp,hdk1,ncld,temnow,vornow,divnow)
!
!--------------------------------------------------------------------
      return
      end
!
!--------------------------------------------------------------------
      subroutine filter_top(jtrun,jtmax,lev,ktop,ncld,temnow     &
                       ,vornow,divnow)
!
!  apply Lanczos filter to top "ktop" layers
!
      use index
      use mpe
      use const, only: RTYPE
!
      implicit  none

      integer   ktop,ktopm1
!      parameter ( ktop=10, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!     parameter ( ktop=4, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!     parameter ( ktop=6, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!
      integer   jtrun,jtmax,lev,ncld

      real(kind=RTYPE) temnow(lev,2,jtrun,jtmax),  &
                       vornow(lev,2,jtrun,jtmax),  &
                       divnow(lev,2,jtrun,jtmax)
!
      real      wvn_top(ktop+1),djt

      integer   k,mode,m,mf,n,nflt,IERR,KL
      real      pi,flt,fac
!

!2dMPI >
!     if(ktop.gt.levp)then
!        print *,'filter_top fatal: ktop greater than lev partial !'
!        call MPI_FINALIZE(IERR)
!        stop
!     endif
!2dMPI <

!      wvn_top(1) = jtrun*2./3.
      wvn_top(1) = max(min(jtrun/3.,155.),55.)
      wvn_top(ktop+1) = jtrun
!!      djt = ( wvn_top(ktop) - wvn_top(1) ) / ktopm1

      do k = 2, ktop
        djt = ( wvn_top(ktop+1) - wvn_top(1) ) * exp(-0.7*(k-1))
!        wvn_top(k) = (jtrun + wvn_top(k-1))*0.5
        wvn_top(k) = min( wvn_top(ktop+1) - djt , float(jtrun) )
      enddo
!
      pi = 3.141596
!
!  mode = 0 : just truncate into assigned wavenumbers without 
!             extra filtering
!  mode = 1 or other : add fitering along with truncating
!
      mode = 1
!
      if( mode .eq. 0 ) then
        do k = 1, lev
!2dMPI >
        KL=Llist(k)
        if( KL .le. ktop ) then
!2dMPI <
          do m = 1, mlistnum
            mf=max(2,mlist(m))
            do n = mf, jtrun
              nflt = int ( wvn_top(k) / float(n) )
              flt = min( 1.0, float(nflt) )
              vornow(k,1,n,m)= vornow(k,1,n,m)*flt
              divnow(k,1,n,m)= divnow(k,1,n,m)*flt
              temnow(k,1,n,m)= temnow(k,1,n,m)*flt
              vornow(k,2,n,m)= vornow(k,2,n,m)*flt
              divnow(k,2,n,m)= divnow(k,2,n,m)*flt
              temnow(k,2,n,m)= temnow(k,2,n,m)*flt
            enddo
          enddo
!          do m = 1, mlistnum
!            mf=max(2,mlist(m))
!            do n = mf, jtrun
!              nflt = int ( wvn_top(k) / float(n) )
!              flt = min( 1.0, float(nflt) )
!              do nc=1,ncld
!                kk=(nc-1)*lev+k
!                qnow(kk,1,n,m)  = qnow(kk,1,n,m)*flt
!                qnow(kk,2,n,m)  = qnow(kk,2,n,m)*flt
!              enddo
!            enddo
!          enddo
!2dMPI >
        endif
!2dMPI <
        enddo
      else
        do k = 1, lev
!2dMPI >
        KL=Llist(k)
        if( KL .le. ktop ) then
!2dMPI <
          do m = 1, mlistnum
            mf=max(2,mlist(m))
            do n = mf, jtrun
              fac = min(wvn_top(k),float(n-1)) * pi / wvn_top(k)
              flt = sin(fac)/fac
              vornow(k,1,n,m)= vornow(k,1,n,m)*flt
              divnow(k,1,n,m)= divnow(k,1,n,m)*flt
              temnow(k,1,n,m)= temnow(k,1,n,m)*flt
              vornow(k,2,n,m)= vornow(k,2,n,m)*flt
              divnow(k,2,n,m)= divnow(k,2,n,m)*flt
              temnow(k,2,n,m)= temnow(k,2,n,m)*flt
            enddo
          enddo
!          do m = 1, mlistnum
!            mf=max(2,mlist(m))
!            do n = mf, jtrun
!              fac = min(wvn_top(k),float(n-1)) * pi / wvn_top(k)
!              flt = sin(fac)/fac
!              do nc=1,ncld
!                kk=(nc-1)*lev+k
!                qnow(kk,1,n,m)  = qnow(kk,1,n,m)*flt
!                qnow(kk,2,n,m)  = qnow(kk,2,n,m)*flt
!              enddo
!            enddo
!          enddo
!2dMPI >
        endif
!2dMPI <
        enddo
      endif
!  
      return
      end
