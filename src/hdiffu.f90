      subroutine hdiffu_3tl ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld,amp&
                        , rad,cosl,ut,vt,vornow,divnow,temnow,eps4     &
                        , trefs)
      use index
      use mpe
      use rank
      use const, only : hdk1,hdk2,radsq,vd,RTYPE
      use param, only : octahedral

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      rad

      real(kind=RTYPE) dta
      real(kind=RTYPE) cosl(my),ut(nxp,lev,my_max),vt(nxp,lev,my_max),  &
                vornow(levp,2,jtrun,jtmax),divnow(levp,2,jtrun,jtmax),   &
                temnow(levp,2,jtrun,jtmax),eps4(jtrun,jtmax),            &
                trefs(levp,2,jtrun,jtmax)
!
!     parameter ( ktop=4, ktop2=ktop/2 ) ! top "ktop" levels are inhenced
!
      real      wmax(lev)
      real      windmax1,windmax2,windmax3

      integer   jj,j,nxj,k,i,m,n,mf,nc,kk,KL
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu
      real      hfilt,hfilt2,nf,kfac
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
      do k=1,lev
        if( wmax(k) .gt. windmax3 ) then
          if(myrank.eq.0)print *,'wmax gt windmax at','k= ',k,         &
                                 ' windmax=',wmax(k)
        endif
      enddo
!
!
        hfilt = (radsq/(nf*(nf+1)))**2.
        hfilt2 = radsq/(nf*(nf+1))
      if ( octahedral ) then
        hfilt = hfilt/(6.*dta)
        hfilt2 = hfilt2/(6.*dta)
      else
        hfilt = 16.*hfilt/dta
        hfilt2 =16.*hfilt2/dta
      endif

      do 100 k=1,levp  ! levp -> lev
!
!       if( wmax(k) .gt. windmax2 ) then
!         if(myrank.eq.0)print *,'wmax gt windmax at','k= ',k,         &
!                                ' windmax=',wmax(k)
!
!  compute diffusion coefficients
!

        KL=Llist(k)
!
!            facd= 1.
!            facv= 1.
!            fact= 1.

!          if ( KL .le. hdk2 ) then
        kfac = 1.0 + 1. * min(max(float(hdk2(1)-KL),0.),60.)
        kfac = kfac*(1.+vd*exp(-0.5*max(float(KL-hdk1),0.)))
!!        facd = 2.* amp * (kfac + 1.*max(float(hdk1-KL),0.))
!!        facv = amp * (kfac + 0.5*max(float(hdk1-KL),0.))
!!        fact = amp * (kfac + 0.5*max(float(hdk1-KL),0.))
        facd = 2.* amp * kfac
        facv = amp * kfac
        fact = amp * kfac
!          endif


!        if ( KL .le. hdk1 ) then
!          ddiffu =facd*hfilt2
!        else
!          ddiffu =facd*hfilt
!        endif

!        tdiffu =fact*hfilt
!        vdiffu =facv*hfilt


!
!  difuse vorticity and divergence fields
!  diffuse moisture and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun

            c1=1.+dta*facv*hfilt*eps4(n,m)**2
            c3=1.+dta*fact*hfilt*eps4(n,m)**2

            if ( KL .le. hdk1 ) then
              c2=1.+dta*facd*hfilt2*eps4(n,m)
            else
              c2=1.+dta*facd*hfilt*eps4(n,m)**2
            endif

            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
!            temnow(k,1,n,m)=temnow(k,1,n,m)/c3
!            temnow(k,2,n,m)=temnow(k,2,n,m)/c3
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
!      if (wmax(1).gt.windmax3 .or. wmax(2).gt.windmax3 .or. &
!          wmax(3).gt.windmax3)then
      if ( windchk ) then
       call filter_top_3tl(jtrun,jtmax,levp,ncld,temnow,vornow,divnow)
      end if
!--------------------------------------------------------------------
      return
      end            
!      
      subroutine hdiffu ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld,amp   &
                        , rad,cosl,ut,vt,vornow,divnow,temnow,eps4     &
                        , trefs)
      use index
      use mpe
      use rank
      use const, only : hdk1,hdk2,radsq,doskeb,onocos,wcfac,wdfac      &
                      , poly,dpoly,hord,vd,factop,RTYPE
      use param, only : octahedral

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      rad

      real(kind=RTYPE) dta
      real(kind=RTYPE) vordiss(levp,2,jtrun,jtmax),divdiss(levp,2,jtrun,jtmax)

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
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu
      real      hfilt,nf,dec,coefu,powd,kfac,dect
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
      powd = float(hord) / 2.
      hfilt  = (radsq/(nf*(nf+1)))**powd
      coefu = factop/float(hdk2(1)-hdk1)
      if ( octahedral ) then
        hfilt  = hfilt/(6.*dta)
      else
        hfilt  = hfilt/dta
      endif

      do 100 k=1,levp  ! levp -> lev
!
!  compute diffusion coefficients
!
        KL=Llist(k)
!
        kfac = min(coefu*max(float(hdk2(1)-KL),0.),factop)!  & 
!             + min(1.*max(float(hdk2(3)-KL),0.),4.)
!        if ( KL .le. hdk1 ) kfac = kfac*(1.+vd*exp(-0.5*KL))
!        kfac = kfac*(1.+vd*exp(-0.2*max(float(KL-hdk1),0.)))
        dect = float(min(max(hdk1-KL,1-hdk1),hdk1-1))/float((hdk1-1))
        if ( dect .ge. 0. ) then
          dec  = 0.5*(1.+dect**(1./3.))
        else
          dect = -1.*dect
          dec  = 0.5*(1.-dect**(1./3.))
        endif
        kfac = kfac*(1.+vd*dec)
        facd = max(1.,kfac)*amp
        facv = max(1.,kfac)*amp
        fact = max(1.,kfac)*amp
!!        facd = 1. * amp * (kfac + 2.*max(float(hdk1-KL),0.))
!!        facv = 1. * (kfac + 1.*max(float(hdk1-KL),0.))
!!        fact = 1. * (kfac + 1.*max(float(hdk1-KL),0.))
!


!
!  difuse vorticity, divergence and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun

!            if ( KL .le. hdk1 ) then
!              c1=1.+dta*facv*hfilt*eps4(n,m)**powd+vd*exp(-0.7*k)
!              c2=1.+dta*facd*hfilt*eps4(n,m)**powd+vd*exp(-0.7*k)
!              c3=1.+dta*fact*hfilt*eps4(n,m)**powd+vd*exp(-0.7*k)
!            else
              c1=1.+dta*facv*hfilt*eps4(n,m)**powd
              c2=1.+dta*facd*hfilt*eps4(n,m)**powd
              c3=1.+dta*fact*hfilt*eps4(n,m)**powd
!            endif

!  if doskeb = .true. estimate the dissipation of kinectic energy for SKEB
            vordiss(k,1,n,m)=(1.-1./c1)*vornow(k,1,n,m)
            vordiss(k,2,n,m)=(1.-1./c1)*vornow(k,2,n,m)
            divdiss(k,1,n,m)=(1.-1./c2)*divnow(k,1,n,m)
            divdiss(k,2,n,m)=(1.-1./c2)*divnow(k,2,n,m)
!
            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
!            temnow(k,1,n,m)=temnow(k,1,n,m)/c3
!            temnow(k,2,n,m)=temnow(k,2,n,m)/c3
          enddo
        enddo
 100  continue
!
!       estimate the dissipation of kinetic energy
      if ( doskeb ) then
         ! transfer the dissipation to grid point from spectrum
         call tranuv (jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac &
                    , poly,dpoly,vordiss,divdiss,ut,vt,nsizey)
      endif

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
      do k=1,hdk1
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
      use param, only : octahedral

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      rad
 
      real(kind=RTYPE) dta
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
      real      hfilt,hfilt2,nf,ncut,ncor,kfac
      real      c1,c2,c3
      logical   windchk

      data      windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!

!
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
         kfac  = min(5.*max(float(hdk2(1)-KL),0.),60.)
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
            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
!            temnow(k,1,n,m)=temnow(k,1,n,m)/c3
!            temnow(k,2,n,m)=temnow(k,2,n,m)/c3
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
      do k=1,hdk1
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
                        , rad,cosl,ut,vt,vornow,divnow,temnow          &
                        , eps4,trefs) 
      use index
      use mpe
      use rank
      use const, only : hdk1,hdk2,radsq,vd,factop,RTYPE,hord
      use param, only : octahedral

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      rad

      real(kind=RTYPE) dta
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
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu,dec,dect
      real      hfilt,nf,kfac,fl,powd
      real      c1,c2,c3,c4
      logical   windchk

      data      windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!
!      temnow = 0.0
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
      powd = float(hord) / 2.
      fl   = factop/float(hdk2(1)-hdk1)
      hfilt = (radsq/(nf*(nf+1)))**powd
      if ( octahedral ) then
        hfilt = hfilt/(6.*dta)
      else
        hfilt = hfilt/dta
      endif

      do 100 k=1,levp  ! levp -> lev
!
!  compute diffusion coefficients
!

        KL=Llist(k)
!
        kfac = min(fl*max(float(hdk2(1)-KL),0.),factop)!    &
!              + min(1.*max(float(hdk2(3)-KL),0.),4.)
!        if ( KL .le. hdk1 ) kfac = kfac*(1.+vd*exp(-0.5*KL))
!        kfac = kfac*(1.+vd*exp(-0.5*max(float(KL-hdk1),0.)))
        dect = float(min(max(hdk1-KL,1-hdk1),hdk1-1))/float((hdk1-1))
        if ( dect .ge. 0. ) then
          dec  = 0.5*(1.+dect**(1./3.))
        else
          dect = -1.*dect
          dec  = 0.5*(1.-dect**(1./3.))
        endif
        kfac = kfac*(1.+vd*dec)

!        facd = mwhd * max(amp,kfac)
!        facv = max(min(amp,1.),kfac)
        facd = max(1.,kfac)*amp
        facv = max(1.,kfac)*amp
        fact = max(1.,kfac)*amp
!!        fact = amp * kfacv 
!          endif

!
!  difuse vorticity and divergence fields
!  diffuse moisture and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun

!            if ( KL .le. hdk1 ) then
!              c1=1.+dta*facv*hfilt4*eps4(n,m)**2.+vd*exp(-0.7*k)
!              c2=1.+dta*facd*hfilt4*eps4(n,m)**2.+vd*exp(-0.7*k)
!              c3=1.+dta*fact*hfilt4*eps4(n,m)**2.+vd*exp(-0.7*k)
!            else
              c1=1.+dta*facv*hfilt*eps4(n,m)**powd
              c2=1.+dta*facd*hfilt*eps4(n,m)**powd
              c3=1.+dta*fact*hfilt*eps4(n,m)**powd
!            endif





            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
!            temnow(k,1,n,m)=temnow(k,1,n,m)/c3
!            temnow(k,2,n,m)=temnow(k,2,n,m)/c3
          enddo
        enddo
 100  continue
!!
!
      windchk=.false.
      do k=1,hdk1
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
        djt = ( wvn_top(ktop+1) - wvn_top(1) ) * exp(-0.5*(k-1))
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

!
!--------------------------------------------------------------------
      subroutine filter_skeb(jtrun,jtmax,lev,dissest,wvn_top)
!
!  apply Lanczos filter to estimation of kinetic energy dissipation.
!
      use index
      use mpe
      use const, only : RTYPE
!
      implicit  none

!
      integer   jtrun,jtmax,lev,ncld

      real(kind=RTYPE) dissest(lev,2,jtrun,jtmax)
!
      real      wvn_top

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
          do m = 1, mlistnum
            mf=max(2,mlist(m))
            do n = mf, jtrun
              nflt = int ( wvn_top / float(n) )
              flt = min( 1.0, float(nflt) )
              dissest(k,1,n,m)= dissest(k,1,n,m)*flt
              dissest(k,2,n,m)= dissest(k,2,n,m)*flt
            enddo
          enddo
!2dMPI <
        enddo
      else
        do k = 1, lev
!2dMPI >
        KL=Llist(k)
          do m = 1, mlistnum
            mf=max(2,mlist(m))
            do n = mf, jtrun
              fac = min(wvn_top,float(n-1)) * pi / wvn_top
              flt = sin(fac)/fac
              dissest(k,1,n,m)= dissest(k,1,n,m)*flt
              dissest(k,2,n,m)= dissest(k,2,n,m)*flt
            enddo
          enddo
!2dMPI <
        enddo
      endif
!  
      return
      end
!
!--------------------------------------------------------------------
      subroutine filter_top_3tl(jtrun,jtmax,lev,ncld,temnow     &
                       ,vornow,divnow)
!
!  apply Lanczos filter to top "ktop" layers
!
      use index
      use mpe
      use const, only : RTYPE
!
      implicit  none

      integer   ktop,ktopm1
      parameter ( ktop=10, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!     parameter ( ktop=4, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!     parameter ( ktop=6, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!
      integer   jtrun,jtmax,lev,ncld

      real(kind=RTYPE) temnow(lev,2,jtrun,jtmax), &
                       vornow(lev,2,jtrun,jtmax), &
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

!      wvn_top(1) = jtrun*1./3.
      wvn_top(1) = 155
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
