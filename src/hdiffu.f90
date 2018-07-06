      subroutine hdiffu ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld,hfilt  &
                        , rad,cosl,ut,vt,vornow,divnow,temnow,eps4     &
                        , trefs)
      use index
      use mpe
      use rank
      use const, only : hdktop,hdk1,hdk2,hdk3,factop,coefu

      implicit  none

      integer   my,my_max,nx,jtrun,jtmax,lev,ncld
      real      dta,rad,hfilt

      real      cosl(my),ut(nxp,lev,my_max),vt(nxp,lev,my_max),  &
                vornow(levp,2,jtrun,jtmax),divnow(levp,2,jtrun,jtmax),   &
                temnow(levp,2,jtrun,jtmax),eps4(jtrun,jtmax),            &
                trefs(levp,2,jtrun,jtmax)
!
!     parameter ( ktop=4, ktop2=ktop/2 ) ! top "ktop" levels are inhenced
!
      real      wmax(lev)
      real      windmax1,windmax2,windmax3

      integer   jj,j,nxj,k,i,m,n,mf,nc,kk,KL
      real      xx,facd,facv,fact,amp,ddiffu,vdiffu,tdiffu,dec,dect
      real      c1,c2,c3

      data      windmax1/80./, windmax2/100./, windmax3/130./
!!      data      windmax1/70./, windmax2/100./, windmax3/130./
!
!cc    diffu=1.0/(24.*3600.*(jtrun*(jtrun+1)/radsq)**2)
!
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
        if( wmax(k) .gt. windmax2 ) then
          if(myrank.eq.0)print *,'wmax gt windmax at','k= ',k,         &
                                 ' windmax=',wmax(k)
        endif
      enddo

      do 100 k=1,levp  ! levp -> lev
!
!       if( wmax(k) .gt. windmax2 ) then
!         if(myrank.eq.0)print *,'wmax gt windmax at','k= ',k,         &
!                                ' windmax=',wmax(k)
!
!  compute diffusion coefficients
!
!ch      dec=1.4*min(max(k-k1,0),k2-k1)                     &
!ch         +0.8*min(max(k-k3,0),ktop-k3)+0.3*min(19-k,0)

         KL=Llist(k)
         dec=max(min(coefu*(hdktop-KL),factop),0.)               &
!!            -max(min(0.4*(hdk3-KL),1.2),0.) 
            -max(min(1.5*coefu*(hdk3-KL),factop),0.) 

!!!         dect=max(min(coefu*(hdktop-5-KL),factop),0.)            &
!!!             -max(min(1.5*coefu*(hdk3-5-KL),factop),0.)


!
!          facd= factop - dec
!!!          facv= 1.0 + dec
          fact= 1.0 + dec


!ch       amp = max(min(0.91*(22-k),15),1)
          amp = min(1.+1.0*max(hdk1-KL,0.),8.)

          facd = max(fact*amp,1.0)
          facv = max(fact*amp,1.0)
          fact = 0.25*max(fact*amp,1.0)
!!!          fact = 0.2*amp
!
        ddiffu =hfilt*facd
        tdiffu =hfilt*fact
        vdiffu =hfilt*facv
!
!  difuse vorticity and divergence fields
!  diffuse moisture and temperature fields
!
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            c1=1.+dta*vdiffu*eps4(n,m)**2
            c2=1.+dta*ddiffu*eps4(n,m)**2
            c3=1.+dta*tdiffu*eps4(n,m)**2
            vornow(k,1,n,m)=vornow(k,1,n,m)/c1
            vornow(k,2,n,m)=vornow(k,2,n,m)/c1
            divnow(k,1,n,m)=divnow(k,1,n,m)/c2
            divnow(k,2,n,m)=divnow(k,2,n,m)/c2
            temnow(k,1,n,m)=(temnow(k,1,n,m)+(c3-1.)*trefs(k,1,n,m))/c3
            temnow(k,2,n,m)=(temnow(k,2,n,m)+(c3-1.)*trefs(k,2,n,m))/c3
          enddo
        enddo
!!        do m=1,mlistnum
!!          mf=mlist(m)
!!          do n=mf,jtrun
!!            c3=1.+dta*qdiffu*eps4(n,m)**2
!!            do nc=1,ncld
!!              kk=(nc-1)*lev+k
!!              qnow(kk,1,n,m)=(qnow(kk,1,n,m)+(c3-1.)*qrefs(kk,1,n,m))/c3
!!              qnow(kk,2,n,m)=(qnow(kk,2,n,m)+(c3-1.)*qrefs(kk,2,n,m))/c3
!!            enddo
!!          enddo
!!        enddo
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
      if (wmax(1).gt.windmax3 .or. wmax(2).gt.windmax3 .or. &
          wmax(3).gt.windmax3)then
       call filter_top(jtrun,jtmax,levp,ncld,temnow,vornow,divnow)
      end if
!--------------------------------------------------------------------
      return
      end

      subroutine filter_top(jtrun,jtmax,lev,ncld,temnow     &
                       ,vornow,divnow)
!
!  apply Lanczos filter to top "ktop" layers
!
      use index
      use mpe
!
      implicit  none

      integer   ktop,ktopm1
      parameter ( ktop=6, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!     parameter ( ktop=4, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!     parameter ( ktop=6, ktopm1=ktop-1 ) ! top "ktop" levels are filtered
!
      integer   jtrun,jtmax,lev,ncld

      real      temnow(lev,2,jtrun,jtmax),                           &
                vornow(lev,2,jtrun,jtmax),divnow(lev,2,jtrun,jtmax)
!
      real      wvn_top(ktop+1),djt

      integer   k,mode,m,mf,n,nflt,IERR
      real      pi,flt,fac
!

!2dMPI >
!     if(ktop.gt.levp)then
!        print *,'filter_top fatal: ktop greater than lev partial !'
!        call MPI_FINALIZE(IERR)
!        stop
!     endif
!2dMPI <

!      wvn_top(1) = jtrun*0.5
      wvn_top(1) = 155
      wvn_top(ktop) = jtrun
      djt = ( wvn_top(ktop) - wvn_top(1) ) / ktopm1
      do k = 2, ktopm1
!        wvn_top(k) = (jtrun + wvn_top(k-1))*0.5
        wvn_top(k) = min( wvn_top(k-1) + djt , jtrun )
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
        do k = 1, ktop
!2dMPI >
        if((K.ge.Lstart) .and. (K.le.Lend))then
!2dMPI <
          do m = 1, mlistnum
            mf=max(2,mlist(m))
            do n = mf, jtrun
              nflt = int ( wvn_top(k) / float(n) )
              flt = min( 1.0, float(nflt) )
              vornow(k,1,n,m)= vornow(k,1,n,m)*flt
              divnow(k,1,n,m)= divnow(k,1,n,m)*flt
!!              temnow(k,1,n,m)= temnow(k,1,n,m)*flt
              vornow(k,2,n,m)= vornow(k,2,n,m)*flt
              divnow(k,2,n,m)= divnow(k,2,n,m)*flt
!!             temnow(k,2,n,m)= temnow(k,2,n,m)*flt
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
        do k = 1, ktop
!2dMPI >
        if((K.ge.Lstart) .and. (K.le.Lend))then
!2dMPI <
          do m = 1, mlistnum
            mf=max(2,mlist(m))
            do n = mf, jtrun
              fac = min(wvn_top(k),float(n-1)) * pi / wvn_top(k)
              flt = sin(fac)/fac
              vornow(k,1,n,m)= vornow(k,1,n,m)*flt
              divnow(k,1,n,m)= divnow(k,1,n,m)*flt
!!              temnow(k,1,n,m)= temnow(k,1,n,m)*flt
              vornow(k,2,n,m)= vornow(k,2,n,m)*flt
              divnow(k,2,n,m)= divnow(k,2,n,m)*flt
!!              temnow(k,2,n,m)= temnow(k,2,n,m)*flt
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
