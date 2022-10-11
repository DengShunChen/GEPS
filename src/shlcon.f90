      subroutine shlcon ( nxj,nx,lev,ktshl,dt,g,r,cp,xkapa,hltm,ptop,dsigma &
                        , tg,pok,pst,topo,phi,po,tt,qt,nshl,rcup,ncld )
!
!     shallow convection parameterization, vertical mixing of t and q
!     (following hogen et al., 1990, tiedtke, et al., 1988)
!
!     c.-s. liou  april 1992
!     modify to f90 by C-H Lee and sort by River Chen in 2015
!
!     conditions for shallow convection:
!     (1) tg       > t surface air (ts)
!     (2) rh(lev)   > 0.8
!     (3) surface lcl is within lowest 200 mb of the atmosphere
!     (4) the atmosphere is moist unstable somewhere in lowest 200 mb
!
!     mixing coef.:
!     xsc = 10.0 m2/s,  from surface to active cloud layer top
!     xsc = 2.0  m2/s,  at layer just above the cloud top (inversion)
!     xsc = 0.0      ,  elsewhere
!
!     parameters:
!
!     nx : horizontal dimension
!    lev : vertical dimension
!  ktshl : highest sigma level allowed for shallow convection
!     dt : time step for shallow convection adjustment
!      g : gravity
!      r : gas constant
!     cp : heat specific of air
!  xkapa : r/cp
!   hltm : constant of latent heat release for water vaper
!   ptop : model top pressure (mb),
!   dsigma : sigma level thickness (lev),
!     tg : ground surface temperature (nx),                 (k)
!    pst : (terrain pressure - ptop) (nx),                  (mb)
!    pok : odd level pressure**kapa (nx,lev),               (mb)
!   topo : terrain geopotential (nx),                       (m2/s2)
!    phi : geopotential on odd levels                       (m2/s2)
!     po : odd level pressure (nx,lev),                     (mb)
!     tt : real temperature on odd levels (nx,lev),         (k)
!     qt : specific humidity on odd levels (nx,lev),        (kg/kg)
!   nshl : total shallow convection active points
!
!  local work arrays
!
!    xsc : mixing coef. due to shallow convection (nx,lev)  (m2/s2)
!

      use paramt
      use const, only: RTYPE

      implicit  none

      integer   nxj,nx,lev,ktshl,ncld

      real      tg(nx),po(nx,lev)                                   &
              , ql(nx,lev),ttmp(nx)
      real(kind=RTYPE) tt(nx,lev),qt(nx,lev*ncld),phi(nx,lev)       &
              ,        topo(nx),pst(nx),dsigma(lev,2),pok(nx,lev)
!
!     local work arrays
!
      logical cld(im)
!
      real      xsc(im,lm),ps(im),aa(im,lm,4)                       &
              , ts(im),wk(im),pp(im),zz(im),hh(im)

      integer   indx0(im),indx(im),ktop0(im),ktop(im)
!
      real      pox(im,lm),ttx(im,lm),phix(im,lm)                   &
               ,tlclx(im),plclx(im),hhx(im,lm),pokx(im,lm)
! for new kcu
      real      tlcl(im),plcl(im),zlcl(im),telcl(im),rcup(im)
      integer   klcl(im),klfc(im),klcl0(im)
!
      real      tmp1(im)
!
      integer nshl,kmax,nxlev,i,k,itt,ii,kclb,ktpmin,ktp,ka,kc,ntrac
      real    ptop,hltm,xkapa,cp,r,g,dt,hlcl,xkapa1,rh,xlog10,fac
      real    tdp,tsx,aln10,p_lcl,ddd,cc,zhig,dhh,www,xsc1,xsc2,tp1,hhu
      real    hhl,wkl,zzl,ddp

      data hlcl / 200.0 /
!
      nshl = 0
      xkapa1 = 1.0/xkapa
      kmax = lev - ktshl + 1
      nxlev = nx * lev
      do 50 i = 1, nxlev
      xsc(i,1) = 0.0
  50  continue
      do 60 i = 1, nxj
      ps(i) = pst(i) + ptop
  60  continue
!
!     check condition (1), assume surface to lev-level has same theda
!
!...modification note: conditional processing has no meaning
!                      if (ps(i)/po(i,lev).lt.0.) do 120 still produces NaNQ
!
!     tmp_min = 1.0
!     do i=1, nxj
!     tmp_min=min(tmp_min,ps(i)/po(i,lev))
!     enddo
!
!     if(tmp_min .gt. 0.) then

!     do 100 i = 1, nxj
!     ts(i) = tt(i,lev) * exp(xkapa*log(ps(i)/po(i,lev)))
!     if ( tg(i) .gt. ts(i) )  then
!        cld(i) = .true.
!     else
!        cld(i) = .false.
!     endif
! 100 continue
!
!     else
!     do 120 i = 1, nxj
!     ts(i) = tt(i,lev) * (ps(i)/po(i,lev))**xkapa
!     if ( tg(i) .gt. ts(i) )  then
!        cld(i) = .true.
!     else
!        cld(i) = .false.
!     endif
! 120 continue
!
!     endif
!---
      do i=1,nxj
        ts(i)=ps(i)/po(i,lev)
      enddo
      call vlog(ts,ts,nxj)
      do i=1,nxj
        ts(i)=xkapa*ts(i)
      enddo
      call vexp(ts,ts,nxj)

      do 100 i = 1, nxj
      ts(i) = tt(i,lev) * ts(i)
      if ( tg(i) .gt. ts(i) )  then
         cld(i) = .true.
      else
         cld(i) = .false.
      endif
  100 continue
!
!     check condition (2)  (use wk to store qsat)
!     besure qt>0.
      do k=1,lev
        do i=1,nxj
          qt(i,k)=max(qt(i,k),1.0e-10)
        enddo
      enddo
!
      ttmp(:)=tt(:,lev)
      call qsatq (nxj, ttmp, po(1,lev), wk )
!
      do 200 i = 1, nxj
      rh = qt(i,lev)/wk(i)
      cld(i) = cld(i) .and. (rh.gt.0.7)
      cld(i) = cld(i) .and. (rcup(i).lt. 1.0e-7)
  200 continue
!
!     check condition (3), use formula given by r.l. inman for t(lcl)
!                          see jam, feb 1969, p155-158, equation 21.
!
!     do 300 i = 1, nxj
!     fac  = log(ps(i)*qt(i,lev)/(6.11*(qt(i,lev)+0.62197)))/log(10.0)
!     tdp  = 237.3*fac/(7.5-fac)
!     tsx  = ts(i) - 273.16
!     tlcl(i) =tdp-(0.212+0.001571*tdp-0.000436*tsx)*(tsx-tdp) + 273.16
!     plcl(i) = po(i,lev)*((tlcl(i)/tt(i,lev))**xkapa1)
!     cld(i) = cld(i) .and. ((ps(i)-plcl(i)).lt.hlcl)
!     cld(i) = cld(i) .and. (po(i,lev).gt.plcl(i))
! 300 continue
!------
      do i=1,nxj
      tmp1(i) = (ps(i)*qt(i,lev)/(6.11*(qt(i,lev)+0.62197)))
      enddo
      call vlog(tmp1,tmp1,nxj)
      xlog10=log(10.0)

      do i = 1, nxj
      fac  = tmp1(i)/xlog10
      tdp  = 237.3*fac/(7.5-fac)
      tsx  = ts(i) - 273.16
      tlcl(i) =tdp-(0.212+0.001571*tdp-0.000436*tsx)*(tsx-tdp) + 273.16
      plcl(i) =tlcl(i)/tt(i,lev)
      enddo
      call vlog(plcl,plcl,nxj)
      do i = 1,nxj
        plcl(i)=xkapa1*plcl(i)
      enddo
      call vexp(plcl,plcl,nxj)

      do 300 i = 1, nxj
      plcl(i) = po(i,lev)*plcl(i)
      cld(i) = cld(i) .and. ((ps(i)-plcl(i)).lt.hlcl)
      cld(i) = cld(i) .and. (po(i,lev).gt.plcl(i))
  300 continue
!
!     check condition (4), use saturated moist static energy to check
!                          moist instability ( hh(k+1) > hh(k) )
!     pack possible shallow convection points for further calculation
!
      itt = 0
      do 400 i = 1, nxj
      if ( cld(i) )  then
         itt = itt + 1
         indx0(itt) = i
      endif
  400 continue
      if ( itt .eq. 0 )  return
!
!  pack array
!
      do k = ktshl-1, lev
      do i = 1, itt
       ii = indx0(i)
       pox(i,k) = po(ii,k)
       ttx(i,k) = tt(ii,k)
       phix(i,k)= phi(ii,k)/g
      end do
      end do
      do i = 1, itt
       ii = indx0(i)
       plclx(i) = plcl(ii)
       tlclx(i) = tlcl(ii)
      end do
!
! search for lcl height and k-level
!
      aln10 = log(10.)
!
      do k = lev, ktshl, -1
       do i = 1, itt
        p_lcl=plclx(i)
        ddd = (p_lcl-pox(i,k))*(p_lcl-pox(i,k-1))
        if( ddd .lt. 0. )then
!ibm        cc=(log(p_lcl)-log(pox(i,k)))
!ibm    +      /(log(pox(i,k-1))-log(pox(i,k)))
         cc=log(p_lcl/pox(i,k))/log(pox(i,k-1)/pox(i,k))
!                              !environment temp at lcl lev
         telcl(i)=(ttx(i,k-1)-ttx(i,k))*cc+ttx(i,k)
         if(cc.ge.0.5)then
          klcl(i)=k-1
          else
          klcl(i)=k
          endif
         endif
       enddo
      enddo

!
      do i = 1, itt
       cld(i) = .false.
      end do
!
!-----------------------------------------------
! search for lfc level
!-----------------------------------------------
      call qsatq ( itt, tlclx, plclx, wk )
      call theta_e(itt,1,plclx,tlclx,wk,hh)
!
      do i = 1, itt
       klfc(i) = klcl(i)
       ktop0(i) = klcl(i)
       cld(i) = tlclx(i).ge.telcl(i)
      end do
!
      do k = ktshl, lev
       call qsatq ( itt, ttx(1,k), pox(1,k), pokx(1,k) )
       call theta_e(itt,1,pox(1,k),ttx(1,k),pokx(1,k),hhx(1,k))
      end do
!
      do k = lev, ktshl, -1
      do i = 1, itt
       if( .not. cld(i) .and. k .le. (klcl(i)-1) )then
         if( hhx(i,k) .lt. hh(i) ) then
           klfc(i) = k+1
           cld(i) = .true.
         end if
       end if
      end do
      end do

!--------------------------------------------------
! find cloud top(ktop)
!--------------------------------------------------
!
      do i = 1, itt
       cld(i) = .false.
      end do
!
      do i = 1, itt
       ii = indx0(i)
       wk(i) = topo(ii)/g
      end do
!
      do k = lev, ktshl, -1
      do i = 1, itt
       if( .not. cld(i) .and. k .le. (klfc(i)-1) )then
       zhig = phix(i,k+1) - wk(i)
       dhh = hh(i) - hhx(i,k)
!                             set shl cu an upper bound on 4000m
        if( zhig .le. 4000. .and. dhh .lt. 0. ) then
         ktop0(i) = k+1
         cld(i) = .true.
        end if
       end if
      end do
      end do

!
!     re-pack indx to select unstable points only
!
      do i=1,itt
      cld(i)=.false.
      enddo
!
      do 480 i = 1, itt
      cld(i)=(ktop0(i).lt.klcl(i)).and.(ktop0(i).lt.klfc(i))
      if ( cld(i) )  then
         nshl = nshl + 1
         indx(nshl) = indx0(i)
         ktop(nshl) = ktop0(i)
         klcl0(nshl)=klcl(i)
      endif
  480 continue
      if ( nshl .eq. 0 )  return

!
!  re-pack array
!
      do k = ktshl-1, lev
      do i = 1, nshl
       ii = indx(i)
       pox(i,k) = po(ii,k)
       pokx(i,k) = pok(ii,k)
       ttx(i,k) = tt(ii,k)
       phix(i,k)= phi(ii,k)/g
      end do
      end do
      do i = 1, nshl
       pp(i) = ps(indx(i))*100.0
      end do
!
!.........................................
! define kc
!.........................................
!
      kclb=lev
      ktpmin=lev
!
      do i = 1, nshl
       ktp=ktop(i)
       kclb=(klcl0(i)+lev)*0.5
       ktpmin= min( ktpmin, ktp )
       www = (phix(i,kclb)+phix(i,ktp)) * 0.5
       hh(i) = 1.0/((www-phix(i,kclb))*(www-phix(i,kclb)))
      end do
!
      do k = ktpmin, lev
      do i = 1, nshl
      ktp=ktop(i)
      kclb=(klcl0(i)+lev)*0.5
      if( k .ge. ktp .and. k .le. kclb ) then
!                               local bomex version
        xsc1 = -10.*hh(i)*(phix(i,k)-phix(i,kclb)) * (phix(i,k)-phix(i,ktp))
!           nonlocal bomex version k=10.*( 1-1/(a+b*F(z)) )   a=1. b=1.5
        xsc2 = -xsc1/(1.+1.5*(-hh(i))*(phix(i,k)-phix(i,kclb)) &
                  *(phix(i,k)-phix(i,ktp)) )
        xsc(i,k) = xsc1+xsc2
      end if
      end do
      end do

!
!     update t and q by vertical mixing due to shallow convection
!
!      (construct coef for tri-diagonal equation system)
!      (use hh ,wk to store air density ro)
!
      tp1 = g*dt
      do 600 i = 1, nshl
      hhu   = pox(i,ktshl-1)*100.0/(r*ttx(i,ktshl-1))
      hh(i) = pox(i,ktshl)*100.0/(r*ttx(i,ktshl))
!     wk(i) = tp1*0.5*(hhu+hh(i))*xsc(i,ktshl)
      wk(i) = tp1*0.5*(hhu+hh(i))*xsc(i,ktshl-1)
      zz(i) = phix(i,ktshl-1) - phix(i,ktshl)
  600 continue
      do 630 k = ktshl, lev
      ka = k - ktshl + 1
      do 610 i = 1, nshl
      if ( k .lt. lev )  then
         hhl = pox(i,k+1)*100.0/(r*ttx(i,k+1))
!        wkl = tp1*0.5*(hh(i)+hhl)*xsc(i,k+1)
         wkl = tp1*0.5*(hh(i)+hhl)*xsc(i,k)
         zzl = phix(i,k) - phix(i,k+1)
      else
         hhl = 0.0
         wkl = 0.0
         zzl = 1.0
      endif
      ddp = pp(i)*dsigma(k,1) + dsigma(k,2)
      aa(i,ka,1) = -wk(i)/(ddp*zz(i))
      aa(i,ka,3) = -wkl  /(ddp*zzl)
      aa(i,ka,2) = 1.0 - aa(i,ka,1) - aa(i,ka,3)
      aa(i,ka,4) = ttx(i,k)/pokx(i,k)
      hh(i) = hhl
      wk(i) = wkl
      zz(i) = zzl
  610 continue
  630 continue
!     fix upper b.c
!
      do 640 i = 1, nshl
      aa(i,1,4)= aa(i,1,4)-aa(i,1,1)*ttx(i,ktshl-1) / pokx(i,ktshl-1)
  640 continue
!
      call trigau ( aa, nx, nshl, lev, kmax, 0 )
!
!     update tt,   and prepare for qt update
!
      do 660 k = ktshl, lev
      ka = k - ktshl + 1
      do 670 i = 1, nshl
      tt(indx(i),k) = aa(i,ka,4)*pokx(i,k)
      aa(i,ka,4) = qt(indx(i),k)
  670 continue
  660 continue
!
      do 700 i = 1, nshl
      aa(i,1,4) = aa(i,1,4) - aa(i,1,1)*qt(indx(i),ktshl-1)
  700 continue
!
      call trigau ( aa, nx, nshl, lev, kmax, 1 )
!
      do 720 k = ktshl, lev
      ka = k - ktshl + 1
      do 710 i = 1, nshl
      qt(indx(i),k) = aa(i,ka,4)
  710 continue
  720 continue
!
! prepare and update ql--------------------------
!
      if(ncld .ge. 2)then
      ntrac=2
      do k=1,lev
      kc=k+(ntrac-1)*lev
      do i=1,nshl
      ql(indx(i),k)=qt(indx(i),kc)
      enddo
      enddo
!
      do 800 k = ktshl, lev
      ka = k - ktshl + 1
      do 820 i = 1, nshl
      aa(i,ka,4) = ql(indx(i),k)
  820 continue
  800 continue
!
      do 900 i = 1, nshl
      aa(i,1,4) = aa(i,1,4) - aa(i,1,1)*ql(indx(i),ktshl-1)
  900 continue
!
      call trigau ( aa, nx, nshl, lev, kmax, 1 )
!
      do 920 k = ktshl, lev
      ka = k - ktshl + 1
      do 910 i = 1, nshl
      ql(indx(i),k) = aa(i,ka,4)
  910 continue
  920 continue

      do k=1,lev
      kc=k+(ntrac-1)*lev
      do i=1,nshl
      qt(indx(i),kc)=ql(indx(i),k)
      enddo
      enddo
!
      endif
!
!xx   print 900, itt, nshl
!x900 format( 1x,'# of shallow convection points (itt,nshl) = ', 2i7 )
      return
      end


        subroutine theta_e(nx,lev,p,t,q,thetae)
!***********************************************************************
! file name        : theta_e.f                                         *
!                                                                      *
! description      : To derive equivalent potential temperature.       *
!                                                                      *
! data file        : none                                              *
!                                                                      *
! I/O parameters   :                                                   *
!                                                                      *
!  input:
!          p(nx,lev)  real        pressure(hPa).                       *
!          t(nx,lev)  real        temperature(K).                      *
!          q(nx,lev) real         specific humility (g/g)              *
!                                                                      *
!  output:
!          thetae(nx,lev) real    equivalent potential temperature(K). *
!
! inter. parameter :                                                   *
!          p0         real        1000. hPa pressure.                  *
!          theta      real        potential temperature(K).            *
!          rgas       real        gas constant for dry air.            *
!          cp         real        specific heat of dry air at constant *
!                                 pressure.                            *
!          aslon      real        constant, = 0.622                    *
!                                                                      *
! called fun./sub. : none                                              *
!                                                                      *
!                                                                      *
! reference        : Bolton,D., 1980: The computation of equivalent    *
!                    potential temperature. MOn. Wea. Rev., 108,       *
!                    1046-1053.                                        *
!                                                                      *
! auther           : Deng Shiung Ming                                  *
!                    - Institude for Information Industry              *
!
! create date      : Jul. 31, 1994                                     *
!----------------------------------------------------------------------*
! modified by      :  Chin-Tzu Fong                                    *
! modified date    :  Apr. 25, 1995                                    *
! modified reason  :  Fit into the structure of GFS                    *
! modify to f90 by C-H Lee and sort by River Chen
!***********************************************************************

        parameter( p0=1000.,rgas=287.,cp=1004.,aslon=0.622 )
        dimension p(nx,lev),t(nx,lev),q(nx,lev),thetae(nx,lev)
        dimension td(nx,lev)
!
        dimension tmp1(nx,lev),tmp2(nx,lev),tmp3(nx,lev)
!
!* To compute dew-point temperature

!       do 50 i = 1, nx*lev

!* To compute vapor pressure from specific humidity.

!       p_td=q(i,1)*p(i,1)/(aslon+q(i,1)*(1.-aslon))

!* To derive dew-point temperature from vapor pressure.

!       rlog=log(p_td/6.112)
!       td(i,1)=243.5*rlog/(17.67-rlog) + 273.15
! 50    continue
!
!       do 100 k = 1, lev
!       do 100 i = 1, nx

!* To compute potential temperature.

!       theta=t(i,k)*(p0/p(i,k))**(rgas/cp)

!* To compute isentropic condensation temperature.

!       tc=1./(1./(td(i,k)-56.)+log(t(i,k)/td(i,k))/800.)+56.

!* To compute mixing ratio.

!       tdz=td(i,k)-273.15
!       p_w=6.112*exp(17.67*tdz/(tdz+243.5))
!       w=aslon*p_w/(p(i,k)-p_w)

!* To compute equivalent potential temperature.

!       if( p(i,k).gt.100. )then
!         thetae(i,k)=theta*exp(2675.*w/tc)
!       else
!         thetae(i,k)=theta
!       endif
! 100   continue
!---

!* To compute dew-point temperature
      do i=1,nx*lev
        p_td=q(i,1)*p(i,1)/(aslon+q(i,1)*(1.-aslon))
        td(i,1)=p_td/6.112
      enddo
      call vlog(td,td,nx*lev)

      do 50 i = 1, nx*lev
!       p_td=q(i,1)*p(i,1)/(aslon+q(i,1)*(1.-aslon))
!       rlog=log(p_td/6.112)
        rlog=td(i,1)
        td(i,1)=243.5*rlog/(17.67-rlog) + 273.15
  50  continue

!* To compute potential temperature.
      do i=1,nx*lev
        tmp1(i,1)=p0/p(i,1)
       enddo
      call vlog(tmp1,tmp1,nx*lev)
      do i=1,nx*lev
        tmp1(i,1)=(rgas/cp)*tmp1(i,1)
        tmp2(i,1)=t(i,1)/td(i,1)
        tdz=td(i,1)-273.15
        tmp3(i,1)=(17.67*tdz/(tdz+243.5))
      enddo
      call vexp(tmp1,tmp1,nx*lev)
      call vlog(tmp2,tmp2,nx*lev)
      call vexp(tmp3,tmp3,nx*lev)

      do 100 i = 1, nx*lev
!       theta=t(i,1)*(p0/p(i,1))**(rgas/cp)
!       tc=1./(1./(td(i,1)-56.)+log(t(i,1)/td(i,1))/800.)+56.
!       tdz=td(i,1)-273.15
!       p_w=6.112*exp(17.67*tdz/(tdz+243.5))
        theta=t(i,1)*tmp1(i,1)
        tc=1./(1./(td(i,1)-56.)+tmp2(i,1)/800.)+56.
        p_w=6.112*tmp3(i,1)
        w=aslon*p_w/(p(i,1)-p_w)
        if( p(i,1).gt.100. )then
          thetae(i,1)=theta*exp(2675.*w/tc)
        else
          thetae(i,1)=theta
        endif
  100 continue
!
        return
        end
