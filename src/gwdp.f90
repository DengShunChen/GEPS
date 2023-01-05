      subroutine gwdp(j,nxj,nx,lev,u,v,t,q,p,pk,pk2,hi,var,dt,g,r,cp,drag, &
                      ugws,vgws,cgw)
!=======================================================================
!
!  1. function description
!
!     perform the gravity wave drag parameterization
!     following the method of palmer et al. 1986
!
!
!  2. common block specification
!     none
!
!
!  3. parameter specification
!  
!  3-1. input and output variable
!
!     j     : latitude index
!     nx    : x-dimension of model grid
!     lev   : total vertical level of model
!     u,v   : wind in e-w & n-s direction (m/s)
!     t     : temperature  (k)
!     q     : specific humidity (kg/kg)
!     p     : pressure at odd-level (mb)
!     pk    : p**kcapa at odd-level ((=p/1000.)**0.286)
!     pk2   : p**kcapa at even-level ((=p/1000.)**0.286)
!     hi    : geopotential height (m2/s2)
!     var   : variance of the terrain field (m2)
!     dt    : time interval (=2.*dt)
!     g     : gravitational constant (=9.806)
!     r     : dry gas constant (=287.)
!     cp    : specific heat constant of the dry air (=1004.)
!     drag  : gravity wave drag force
!     ugws  : zonal component of surface gravity wave stress 
!     vgws  : meridional component of surface gravity wave stress 
!
!  3-2. local variable
!
!     p2    : pressure at even-level (mb)
!     pt    : potential temperature (k)
!     den   : density (kg/m3)
!     sn2   : brunt-vaisala frequency (/s2)
!     ri    : bulk richardson number
!     stress: magnitude of the gravity wave stress (kg/s2*m)
!     wss   : reference level wind speed (m/s)                   
!     wsu   : reference level wind vector in e-w direction 
!     wsv   : reference level wind vector in n-s direction 
!     wsd   : reference level wind direction 
!     snm   : mean brunt-vaisala frequency of reference level
!     rim   : mean richardson number of reference level
!     denm  : mean density of reference level
!
!
!  4. vertical structure description
!
!           ------------------------------------------
!              .......p,u,v,hi,t,q......drag............ pk.1 (odd-level) 
!   pk2.1   --p2,sn,ri,stress-------------------------
!              .......p,u,v,hi,t,q......drag............ pk.2                
!   pk2.2   --p2,sn,ri,stress-------------------------
!              .......p,u,v,hi,t,q......drag............ pk.3                
!                               .
!                               .
!                               .
!              .......p,u,v,hi,t,q......drag............ pk.k-1              
!   pk2.k-1 --p2,sn,ri,stress-------------------------
!              .......p,u,v,hi,t,q......drag............ pk.k                
!   pk2.k   --p2,sn,ri,stress-------------------------
!              .......p,u,v,hi,t,q......drag............ pk.k+1              
!                               .
!                               .
!                               .
!              .......p,u,v,hi,t,q......drag............ pk.LEV-1              
! pk2.LEV-1 --p2,sn,ri,stress-------------------------
!              .......p,u,v,hi,t,q......drag............ pk.LEV              
! pk2.LEV   --p2,sn,ri,stress-------------------------
!  (even-level)
!
!
!  5. calling modules
!     diabat
!
!
!  6. usage
!     call gwdp(j,nx,lev,u,v,t,q,p,pk,pk2,hi,var,dt,g,r,cp)
!
!
!  7. modules called
!     none
!
!
!  8. date
!     created  sep. 1994
!
!
!  9. author
!     c.-h. shiao 
!
!
! 10. reference
!     palmer, t. n., g. j. shutts and r. swinbank, 1986:  alleviation
!          of a systematic westerly bias in circulation and numberical
!          weather prediction models through an orographic gravity-
!          wave drag parameterization.  quart. j. roy. meteor. soc.,
!          112, 1001-1039.
!
!=======================================================================
!
      use paramt
      use const, only : RTYPE

      implicit  none

! ... input and output variable
!
      integer   j,nxj,nx,lev

      real      dt,g,r,cp

      real      p(nx,lev)
      real      drag(nx,lev)
      real      ugws(nx),vgws(nx)
      real(kind=RTYPE) u(nx,lev),v(nx,lev),t(nx,lev),q(nx,lev),   &
                       hi(nx,lev)
      real(kind=RTYPE) pk(nx,lev),pk2(nx,lev),var(nx)
!
! ... local variable
!
      real      p2(im,lm)
      real      pt(im,lm),den(im,lm),sn2(im,lm),ri(im,lm)
      real      stress(im,lm),wri(im,lm)
      real      wss(im),wsu(im),wsv(im),wsd(im)
      real      snm(im),rim(im),denm(im)

      integer   k,i,levm,levref,levt

      real     cgw,varmax,dz1,s1,s2,u1,v1,su2,pp,wsp,sn,rr,aa,bb, & 
               delth,delth2,um,vm,umnew,denx,alp,stmin,dpk,wsold, &
               wsnew


!
! ***** preprocess for gwdp
!
! for L18
!cc   cgw=2.5e-5
! for L30 (2003/01/20)
!cc   cgw=1.25e-4
!
! tuned by phon-lu for better performance (2011/5)
!
!      cgw=7.0e-5
! for L60
!      cgw=2.5e-5
!      cgw=1.0e-4
!      cgw=1.25e-4
!
!for T511      varmax=400.*400.
      varmax=1200.*1200.    !TCo639L72
!byl add
      levt=10 !TCo639L72
!
! ... variance of the terrain field
!
      do 70 k=1,lev
      do 70 i=1,nxj
      stress(i,k)=0.
      drag(i,k)=0.
      wri(i,k)=9999.
   70 continue
!
      do 80 i=1,nxj
      if(var(i).ge.varmax) var(i)=varmax
   80 continue
!
! ... pressure at even-level (p2)
! ... potential temperature (pt)
! ... density (den)
!
      do 100 k=1,lev
      do 100 i=1,nxj
      p2(i,k)=1000.0*pk2(i,k)*pk2(i,k)*pk2(i,k)*sqrt(pk2(i,k))
      pt(i,k)=t(i,k)*(1.0+0.608*q(i,k))/pk(i,k)
      den(i,k)=p(i,k)*100./(r*t(i,k))
  100 continue
!
! ... ri and sn                       
!
      do 130 k=1,lev-1
      do 130 i=1,nxj
!
      dz1=(hi(i,k)-hi(i,k+1))/g
!
! ... brunt-vaisala frequency
!
         s1=pt(i,k)-pt(i,k+1)
         s2=(pt(i,k)+pt(i,k+1))/2.
      sn2(i,k)=g*s1/(s2*dz1)
!
! ... richardson number
!
         u1=u(i,k)-u(i,k+1)
         v1=v(i,k)-v(i,k+1)
      su2=u1*u1+v1*v1
         if(su2.le.1.e-6) su2=1.e-6
!
      ri(i,k)=sn2(i,k)*dz1*dz1/su2
!
      if(sn2(i,k).le.0.) then
      sn2(i,k)=0.
      ri(i,k)=0.
      endif
!
  130 continue
!
      do 140 i=1,nxj
      sn2(i,lev)=sn2(i,lev-1)
      ri(i,lev)=ri(i,lev-1)
  140 continue
!
!======================================================================
!
! L18 ... reference level = lowest 3-level mean
! L30 ... reference level = lowest 5-level mean
! reference level is about at sig=0.9
!
!L18  levref=lev-3    ! L18
!L30  levref=lev-5    ! L30
!L40      levref=lev-7    ! L40
      levref=lev-10    ! L72
!
! fix an error (2011/5)
!err  levm=levref*3/4   ! supposed to be at sig=0.45
!cc   levm=int((levref*0.45)+0.51)
! tuned by phon-lu (2011/5)
!      levm=int((levref*0.4)+0.51)
! for L60(op7)
      levm=int((levref*0.4)+0.51)
!
      do 200 i=1,nxj
      wsu(i)=0.              
      wsv(i)=0.              
  200 continue
      do 210 k=lev-1,levref,-1
      do 210 i=1,nxj
      pp=p2(i,k)-p2(i,k-1)
      wsu(i)=wsu(i)+u(i,k)*pp
      wsv(i)=wsv(i)+v(i,k)*pp
  210 continue
      do 220 i=1,nxj
      wsp=p2(i,lev-1)-p2(i,levref-1)
      wsu(i)=wsu(i)/wsp
      wsv(i)=wsv(i)/wsp
      wss(i)=sqrt(wsu(i)*wsu(i)+wsv(i)*wsv(i))
      if(wss(i).le.1.e-3) then
      wss(i)=0.
      wsu(i)=0.
      wsv(i)=0.
      else      
      wsu(i)=wsu(i)/wss(i)
      wsv(i)=wsv(i)/wss(i)
      endif
         
  220 continue
!
      do 230 i=1,nxj
      snm(i)=0.
      rim(i)=0.
  230 continue
      do 238 k=lev-1,levref-1,-1
      do 235 i=1,nxj
      pp=p(i,k+1)-p(i,k)
      snm(i)=snm(i)+sn2(i,k)*pp 
      rim(i)=rim(i)+ri(i,k)*pp 
  235 continue
  238 continue
      do 240 i=1,nxj
      wsp=p(i,lev)-p(i,levref-1)
      snm(i)=snm(i)/wsp
      rim(i)=rim(i)/wsp
  240 continue
!
      do 245 i=1,nxj
      denm(i)=0.
  245 continue
      do 250 k=lev-1,levref,-1
      do 248 i=1,nxj
      pp=p2(i,k)-p2(i,k-1)
      denm(i)=denm(i)+den(i,k)*pp
  248 continue
  250 continue
      do 255 i=1,nxj
      wsp=p2(i,lev-1)-p2(i,levref-1)
      denm(i)=denm(i)/wsp
  255 continue
!
! ... low-level stress (at levref layer)
!
      do 260 i=1,nxj
!shiao      if(snm(i).le.0. .or. rim(i).le.0.0 ) then
      if(snm(i).le.0. .or. rim(i).le.0.25 ) then
      stress(i,levref)=0.
      else  
        sn=sqrt(snm(i))
        rr=sqrt(rim(i))
        aa=1./rr+2.  
        bb=2.*sqrt(aa)-aa
        delth=(wss(i)*bb)**2./snm(i)
        delth2=min(var(i),delth)
      stress(i,levref)=cgw*denm(i)*sn*wss(i)*delth2
      endif
  260 continue            
      do 265 k=lev,levref+1,-1
      do 265 i=1,nxj
      stress(i,k)=stress(i,levref)
  265 continue
!
!======================================================================
!
! ... surface gravity wave stress
!
      do 300 i=1,nxj
      ugws(i)=stress(i,levref)*wsu(i)
      vgws(i)=stress(i,levref)*wsv(i)
  300 continue
!
!======================================================================
!
! ... stress at any level (from reference level to level 1)
!
      do 610 k=levref-1,1,-1
      do 600 i=1,nxj
!
! ... eliassen-palm theorem
!
      stress(i,k)=stress(i,k+1)
!
! ... any level wind vector projected in the surface wind direction
!
      um=(u(i,k)+u(i,k+1))/2.
      vm=(v(i,k)+v(i,k+1))/2.
      umnew=um*wsu(i)+vm*wsv(i)
!
!shiao      if( (umnew.le.0.0) .or. (ri(i,k).le.0.0) 
      if( (umnew.le.0.0) .or. (ri(i,k).le.0.25)  &
                         .or. (sn2(i,k).le.0.0) ) then
        stress(i,k)=0.0
        wri(i,k)=0.0
      else
!
! .. wave amplitude ( delth )
!
        denx=(den(i,k)+den(i,k+1))/2.
        sn=sqrt(sn2(i,k))
        bb=cgw*denx*sn*umnew
        delth=sqrt(stress(i,k)/bb)
!
! .. wave richardson number ( wri )
!
        rr=sqrt(ri(i,k))
        alp=sn*delth/umnew
        bb=1.+rr*alp
!c        wri=ri(i,k)*(1.-alp)/(bb*bb)
        wri(i,k)=ri(i,k)*(1.-alp)/(bb*bb)
!
! .. saturation hypothesis
!
!c          if(wri.le.0.25) then     
          if(wri(i,k).le.0.25) then     
            aa=1./rr+2.  
            bb=2.*sqrt(aa)-aa
            delth=(umnew*bb)**2./sn2(i,k)
            delth2=min(var(i),delth)
          stress(i,k)=cgw*denx*sn*umnew*delth2
          stress(i,k)=min(stress(i,k+1),stress(i,k)) 
          endif        
      endif
!
! ... modification to the lower boundary       
!
      if(k.ge.levm) then
        stmin=0.5*stress(i,k+1)
        if(stress(i,k).le.stmin) then               
          stress(i,k)=stress(i,k+1)*0.5
        endif
      endif
!
  600 continue
  610 continue
!
! ... modification to the upper boundary
!
!cc      do 700 i=1,nxj
!cc      stress(i,1)=stress(i,2)
!cc  700 continue
!
! ...............................................
! ... new winds which are adjusted by gravity wave drag
!
!byl      do 830 k=3,levref
      do 830 k=levt+1,levref
      do 840 i=1,nxj
!
        dpk=(p2(i,k-1)-p2(i,k))*100.
      drag(i,k)=-g*(stress(i,k-1)-stress(i,k))/dpk
!
         wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
         u(i,k)=u(i,k)+dt*drag(i,k)*wsu(i)
         v(i,k)=v(i,k)+dt*drag(i,k)*wsv(i)
         wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
         t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
!c      drag(i,k)=drag(i,k)*wsu(i) 
!c      drag(i,k)=drag(i,k)*wsv(i) 
!
  840 continue
  830 continue
!
      k=1
      do 850 i=1,nxj
!
      drag(i,k)=0.                
!
         wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
         u(i,k)=u(i,k)+dt*drag(i,k)*wsu(i)
         v(i,k)=v(i,k)+dt*drag(i,k)*wsv(i)
         wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
         t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
!c      drag(i,k)=drag(i,k)*wsu(i) 
!c      drag(i,k)=drag(i,k)*wsv(i) 
!
  850 continue
!
      k=2
!byl add
      do 860 k=levt,2,-1
      do 860 i=1,nxj
!
!byl      drag(i,k)=drag(i,k+1)/2.    
      drag(i,k)=0.5*drag(i,k+1)
!
         wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
         u(i,k)=u(i,k)+dt*drag(i,k)*wsu(i)
         v(i,k)=v(i,k)+dt*drag(i,k)*wsv(i)
         wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
         t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
!c      drag(i,k)=drag(i,k)*wsu(i) 
!c      drag(i,k)=drag(i,k)*wsv(i) 
!
  860 continue
!
      return
      end
