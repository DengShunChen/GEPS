      subroutine xyintpo(flag0,im,jm,flag,nx,my,ain,aout,iwnd  &
                       ,xr,yr,first)
!c
!c   process the interpolation from gg(or ga, gx) to gg(or ga)
!c   gg : the gaussion grid system
!c   ga : the equal distant grid system
!c   gx : the equal distant grid system, but the first point at half
!c        a grid's distance
!c
      use const, only: RTYPE
      real(kind=RTYPE) ain(im,jm)
      real(kind=RTYPE) weight0(jm),sinl0(jm),weight(my),sinl(my)
      dimension aout(nx,my),xr(nx),yr(my)
      character*2 flag0, flag
      integer iwnd
      logical first
!c
!c  working arrays
!c
      dimension alat0(jm+2),alon0(im+1),alat(my),alon(nx)   
      dimension aaa(im+1,jm+2),bbb(nx*my),wkx(nx*my),wky(nx*my)
      integer jend
!c
      pi= 4.0*atan(1.0)

!c-------------------------
      if( first ) then
!c
!c  gaussian quadrature weights and latitudes
!c----
      if( flag0.eq.'gg' .or. flag0.eq.'GG' )then       !--- gg
!c
      call gausl3 (jm,-1.,1.,weight0,sinl0)
      do j = 1, jm
       alat0(j+1) = asin(sinl0(j))
      end do
      alat0(1)    = -pi*0.5
      alat0(jm+2) = pi*0.5
      jend = jm+1
!c
      elseif( flag0.eq.'ga' .or. flag0.eq.'GA' )then   !--- ga
!c
      yy= pi/float(jm-1)
      do j = 1, jm
       alat0(j) = -0.5*pi+(j-1)*yy
      end do
      jend = jm-1
!c
!c  gx : equal distance, but the first point beginning at half delta
!c
      elseif( flag0.eq.'gx' .or. flag0.eq.'GX' )then   !--- gx
      yy= pi/float(jm)
      do j = 1, jm 
       alat0(j+1) = -0.5*pi+0.5*yy+(j-1)*yy
      end do
      alat0(1)    = -pi*0.5
      alat0(jm+2) = pi*0.5
      jend = jm+1
!c
      else
!c
      print*, flag0,' No such grid system --- error '
      stop
!c
      end if
!c
      xx= (2.0*pi)/float(im)
      do i = 1, im
       alon0(i) = (i-1)*xx
      end do
      alon0(im+1) = pi*2.0
!c----
      if( flag.eq.'gg' .or. flag.eq.'GG' )then
!c
      call gausl3 (my,-1.,1.,weight,sinl)
!c
      do j = 1, my
       alat(j) = asin(sinl(j))
      end do
!c
      elseif( flag.eq.'ga' .or. flag.eq.'GA' )then   !--- else
!c
      yy= pi/float(my-1)
      do j = 2, my-1
       alat(j) = -0.5*pi+(j-1)*yy
      end do
       alat(1)  = -0.5*pi+0.0001
       alat(my) =  0.5*pi-0.0001
!c
      else
!c
      print*, flag,' No such grid system for outfields'
      stop
!c
      end if
!c
      if( flag0.eq.'gx' .or. flag0.eq.'GX' )then
!c
!c  shift alon with half grid's distance
!c
      xx= (2.0*pi)/float(nx)
      alon(1) = 2.0*pi-0.5*xx
      do i = 2, nx
       alon(i) = 0.5*xx+(i-1)*xx
      end do
!c
      else
!c
      xx= (2.0*pi)/float(nx)
      do i = 1, nx
       alon(i) = (i-1)*xx
      end do
!c
      end if
!c----
      do ix = 1, nx
      do i  = 1, im
      if( alon(ix).ge.alon0(i) .and. alon(ix).lt.alon0(i+1) )then
         xr(ix) = float(i) + (alon(ix)-alon0(i))/(alon0(i+1)-alon0(i))
      end if
      end do
      end do
!c
      do jy = 1, my
      do j  = 1, jend
      if( alat(jy).ge.alat0(j) .and. alat(jy).lt.alat0(j+1) )then
         yr(jy) = float(j) + (alat(jy)-alat0(j))/(alat0(j+1)-alat0(j))
      end if
      end do
      end do
!c
      end if
!c------------------
!c
!c  check if the fields is of wind
!c
!c----
      if( flag0.eq.'gg' .or. flag0.eq.'GG' .or. flag0.eq.'gx' & 
            .or. flag0.eq.'GX')then
!c
      do j = 1, jm
      do i = 1, im
       aaa(i,j+1) = ain(i,j)
      end do
      end do
!c
      do j = 2, jm+1
       aaa(im+1,j) = aaa(1,j)
      end do
!c
      if( iwnd .eq. 1 ) then
       do i = 1, im+1
        aaa(i,1) = 0.
        aaa(i,jm+2) = 0.
       end do
      else   ! ---- else
       ps= ain(1,1)
       pn= ain(1,jm)
       do i = 2, im
        ps= ps + ain(i,1)
        pn= pn + ain(i,jm)
       end do
        ps = ps/float(im)
        pn = pn/float(im)
       do i = 1, im+1
        aaa(i,1) = ps
        aaa(i,jm+2) = pn
       end do
      end if
!c
      else   !--- else
!c
      ps = ain(1,1)
      pn = ain(1,jm)
      do j = 1, jm
      do i = 1, im
       aaa(i,j) = ain(i,j)
      end do
      end do
!c
      do j = 1, jm
       aaa(im+1,j) = aaa(1,j)
      end do
!c
      end if
!c----
!c
      mn = 0
      do j = 1, my
      do i = 1, nx
       mn = mn+1
       wkx(mn) = xr(i)
       wky(mn) = yr(j)
      end do
      end do
!c
      nxmy = nx*my
!c
      if( flag0.eq.'gg' .or. flag0.eq.'GG' .or. flag0.eq.'gx'  & 
              .or. flag0.eq.'GX' )then
       call bcubij(aaa,im+1,jm+2,bbb,nxmy,wkx,wky)
      else
       call bcubij(aaa,im+1,jm  ,bbb,nxmy,wkx,wky)
      end if
!c
      mn = 0
      do j = 1, my
      do i = 1, nx
       mn = mn+1
       aout(i,j) = bbb(mn)
      end do
      end do
!c
      if( flag.eq.'ga' .or. flag.eq.'GA' ) then
       if( iwnd .eq. 1 ) then
        do i = 1, nx
         aout(i,1) = 0.
         aout(i,my) = 0.
        end do
       else
        do i = 1, nx
         aout(i,1) = ps
         aout(i,my) = pn
        end do
       end if
      end if
!c
      return
      end

      subroutine bcubij(f,ix,jy,dout,mn,xin,yin)
!c         
!c  parameter list   
!c         
!c  f,ix,jy,dout,mn,xin,yin - see below  
!c                
      dimension f(ix,jy),dout(mn),xin(mn),yin(mn) 
!c         
!c
!c  working array
      dimension fxx(ix,jy),fyy(ix,jy),pix(mn,4),pjy(mn,4)  &
     , tp1(mn,4),tp2(mn,4),tp3(ix,jy)     
      dimension ipt(mn) 
!c         
!c          a bicubic spline interpolator to interpolate from a grid   
!c          with constant grid spacing to a grid with constant or      
!c          variable grid spacing. all grids are assumed to have point 
!c          (1,1) in the lower left corner with i increasing to the right        
!c          and j increasing upward.     
!c         
!c          arguments:         
!c         
!c          f(ix,jy): fwa of data array to be interpolated from (given 
!c                    by user) 
!c         
!c          ix: first (i) dimension of f (given by user)     
!c         
!c          jy: second (j) dimension of f (given by user)    
!c         
!c          dout(mn): fwa of array of interpolated values (given on    
!c                     output) 
!c         
!c         
!c         mn: number of points (dimension) of output grid   
!c         
!c          xin(mn): x-coordinates of points in dout relative to the   
!c                   x-coordinates of f. a "1" refers to the leftmost  
!c                   boundary of f (given by user) 
!c         
!c          yin(mn): y-coordinates of points in dout relative to the   
!c                   y-corrdinates of f. a "1" refers to the bottom    
!c                   row of f (given by user)      
!c         
!c   scratch work arrays
!c
!c          pix(mn,4): array to hold coefficients for interpolation in 
!c                     x-direction     
!c         
!c          pjy(mn,4): array to hold coefficients for interpolation in 
!c                     y-direction 
!c         
!c          tp1(mn,4): work space        
!c         
!c          tp2(mn,4): work space        
!c         
!c          tp3(ix,jy): work space       
!c         
!c          fxx(ix,jy): array to hold cubic spline values (computed    
!c                      internally)      
!c         
!c          fyy(ix,jy): array to hold cubic spline values (computed    
!c                      internally)      
!c         
!c          ipt(mn): array that holds 2-d coordinate, relative to f grid,
!c                   of each point in dout 
!c         
!c          compute ipt,jpt,pix and pjy  
!c         
      call stupij(xin,yin,mn,ix,pix,pjy,ipt)
!c         
      ijm2=ix*jy-2  
      ixjym2=ix*(jy-2)
!c         
!c          interpolate in x-direction   
!c         
!c          compute fyy        
!c         
      do 100 i=1,ixjym2       
      fyy(i,2)=(f(i,1)-2.0*f(i,2)+f(i,3))         
  100 continue      
      call trdih(ix,jy-2,fyy(1,2))     
!c         
      do 105 i=1,ix
      fyy(i,1)= 0.0
      fyy(i,jy)= 0.0
  105 continue
!c         
!c          compute fxxyy      
!c         
      do 130 i=1,ijm2         
      fxx(i+1,1)= fyy(i,1)+(fyy(i+2,1)-2.0*fyy(i+1,1))      
  130 continue      
!c         
      do 205 j=1,jy 
      fxx(1,j)= fyy(ix,j)-2.0*fyy(1,j)+fyy(2,j)   
      fxx(ix,j)= fyy(ix-1,j)-2.0*fyy(ix,j)+fyy(1,j)         
  205 continue      
!c
      call tpose(fxx,ix,jy,tp3)        
      call trdiph(jy,ix,tp3)      
      call tpose(tp3,jy,ix,fxx)        
!c         
!c          fxx holds fxxyy and fyy holds fyy      
!c         
      call gathij(1,mn,ix,jy,ipt,fxx,fyy,tp2)    
!c
      do 550 i=1,mn 
      tp2(i,1)= pix(i,1)*tp2(i,1)    
      tp2(i,2)= pix(i,2)*tp2(i,2)    
      tp2(i,3)= pix(i,3)*tp2(i,3)    
      tp2(i,4)= pix(i,4)*tp2(i,4)    
      tp1(i,1)=tp2(i,1)+tp2(i,2)+tp2(i,3)+tp2(i,4)
  550 continue      
!c
      call gathij(0,mn,ix,jy,ipt,fxx,fyy,tp2)    
!c
      do 555 i=1,mn 
      tp2(i,1)= pix(i,1)*tp2(i,1)    
      tp2(i,2)= pix(i,2)*tp2(i,2)    
      tp2(i,3)= pix(i,3)*tp2(i,3)    
      tp2(i,4)= pix(i,4)*tp2(i,4)    
      tp1(i,2)=tp2(i,1)+tp2(i,2)+tp2(i,3)+tp2(i,4)
  555 continue      
!c         
!c          compute fxx        
!c         
      do 170 i=1,ijm2         
      fxx(i+1,1)= f(i,1)+f(i+2,1)-2.0*f(i+1,1)  
  170 continue      
!c
      do 155 j=1,jy 
      fxx(1,j)= f(ix,j)-2.0*f(1,j)+f(2,j)         
      fxx(ix,j)= f(ix-1,j)-2.0*f(ix,j)+f(1,j)     
  155 continue      
!c
      call tpose(fxx,ix,jy,tp3)        
      call trdiph(jy,ix,tp3)
      call tpose(tp3,jy,ix,fxx)        
!c         
      call gathij(1,mn,ix,jy,ipt,fxx,f,tp2)      
!c
      do 560 i=1,mn 
      tp2(i,1)= pix(i,1)*tp2(i,1)    
      tp2(i,2)= pix(i,2)*tp2(i,2)    
      tp2(i,3)= pix(i,3)*tp2(i,3)    
      tp2(i,4)= pix(i,4)*tp2(i,4)    
      tp1(i,3)=tp2(i,1)+tp2(i,2)+tp2(i,3)+tp2(i,4)
  560 continue      
!c
      call gathij(0,mn,ix,jy,ipt,fxx,f,tp2)      
!c
      do 570 i=1,mn 
      tp2(i,1)= pix(i,1)*tp2(i,1)    
      tp2(i,2)= pix(i,2)*tp2(i,2)    
      tp2(i,3)= pix(i,3)*tp2(i,3)    
      tp2(i,4)= pix(i,4)*tp2(i,4)    
      tp1(i,4)=tp2(i,1)+tp2(i,2)+tp2(i,3)+tp2(i,4)
!c         
!c          interpolate in y-direction   
!c
      tp1(i,1)= tp1(i,1)*pjy(i,1)    
      tp1(i,2)= tp1(i,2)*pjy(i,2)    
      tp1(i,3)= tp1(i,3)*pjy(i,3)    
      tp1(i,4)= tp1(i,4)*pjy(i,4)    
      dout(i)=tp1(i,1)+tp1(i,2)+tp1(i,3)+tp1(i,4) 
  570 continue      
      return        
      end
      subroutine tpose(x,im,jm,y)       
      dimension x(im,jm),y(jm,im)       
      do 1 i=1,im   
      do 1 j=1,jm   
    1 y(j,i)= x(i,j)
      return        
      end 
      subroutine trdiph (m,n,y)         
!c         
!c  vectorized periodic gaussian elimination solver
!c         
!      dimension y(m,n),work(10000)
      dimension y(m,n),work(m+3*n)  ! wei 20231019 
!c         
!c gaussian elimination        
!c         
      nt2= 2*n      
      work(n+1)= 0.25         
      v = 1.0       
      work(1) = work(n+1)     
      work(nt2+1) = work(n+1) 
      bn = -v*work(nt2+1)+4.0
      do 101 j=2,n-2
        ne = j+n   
        work(ne) = 1.0/(4.0-work(j-1))
        work(j) = work(ne)   
        nu = j+nt2 
        work(nu) = -work(nu-1)*work(ne)
        v = -v*work(j-1)     
        bn = bn-v*work(nu)   
  101 continue      
!c
      v = 1.0-v*work(n-2)     
      ne = nt2      
      work(ne-1) = 1.0/(4.0-work(n-2)) 
      nu = nt2+n    
      work(n-1) = (1.0-work(nu-2))*work(ne-1)      
      work(ne) = 1.0/(bn-v*work(n-1))    
!c
      v= 1.0        
!DIR$ IVDEP
!ocl  novrec
      do 201 i=1,m  
      y(i,1)= y(i,1)*work(n+1)
      work(i+3*n)=y(i,n)-v*y(i,1)       
  201 continue      
!c
      do 103 j=2,n-2
      v= -v*work(j-1)         
      ne = j+n   
!DIR$ IVDEP
!ocl  novrec
      do 113 i=1,m  
      y(i,j)= (y(i,j)-y(i,j-1))*work(ne)
      work(i+3*n)= work(i+3*n)-v*y(i,j) 
  113 continue      
  103 continue      
      v = 1.0-v*work(n-2)     
      ne = nt2      
!DIR$ IVDEP
!ocl  novrec
      do 203 i=1,m  
      y(i,n-1)=(y(i,n-1)-y(i,n-2))*work(ne-1)       
      work(i+3*n)= work(i+3*n)-v*y(i,n-1)
!c         
!c backwards substitution      
!c         
      y(i,n)= work(i+3*n)*work(ne)      
      y(i,n-1)= y(i,n-1)-work(n-1)*y(i,n)  
  203 continue      
!c
      do 104 k=n-2,1,-1
      nu = k+nt2 
!DIR$ IVDEP
!ocl  novrec
      do 104 i=1,m  
      y(i,k)= y(i,k)-work(k)*y(i,k+1)-work(nu)*y(i,n)       
  104 continue      
      return        
      end 
      subroutine gathij(kp,mn,ix,jy,ijpt,fxx,f,tp2)
!c
      dimension fxx(ix*jy),f(ix*jy),tp2(mn,4)
      dimension ijpt(mn)
!c
      if(kp.eq.0) then
!c
!c  lower left corner of interpolation square
!c
      do 10 i=1,mn
      inx= ijpt(i)-ix
      tp2(i,2)= fxx(inx)
      tp2(i,4)= f  (inx)
   10 continue
!c
      else
!c
!c  upper left corner of interpolation square
!c
      do 20 i=1,mn
      tp2(i,2)= fxx(ijpt(i))
      tp2(i,4)= f  (ijpt(i))
   20 continue
      endif
!c
      if(kp.eq.0) then
!c  lower right corner of interpolation square
      i1= 1-ix
      else
!c  upper right corner of interpolation square
      i1= 1
      endif
!c
      do 30 i=1,mn
      inx= ijpt(i)+i1
      if(mod(inx,ix).eq.1) inx= inx-ix
      tp2(i,1)= fxx(inx)
      tp2(i,3)= f  (inx)
   30 continue
!c
      return
      end


      subroutine prexp ( nx,lev,ptop,sig,pt,pk,pk2,plt )
!c
!c  subroutine to compute p to the kapa on odd and even levels
!c
!c *** input ****
!c
!c  nx: e-w dimension no.
!c  lev: number of vertical levels
!c  sig: sigma levels
!c  pt: terrain pressure
!c
!c *** output ***
!c
!c  pk: odd (full) level p**capa
!c  pk2: even( half) level p**capa
!c
!c **************************************************
!c
!c
      dimension pt(nx),pk2(nx,lev),pk(nx,lev),sig(lev+1) &
      ,         plt(nx,lev)
!c
!csun  include '../include/paramt.h'  .. change im to nx
      dimension pl2(nx,2)
!c
!c  compute  pressure variables
!c
      capa  = 1.0/3.5
      capap1= 1.0+capa
      opok  = 1.0/1000.0**capa
!c
      kbot= 1
      k= 1
      ptopk= ptop*opok*ptop**capa
      do 80 i=1, nx
      pl2(i,kbot)= sig(2)*pt(i)+ptop
      pk2(i,1)   = opok*pl2(i,kbot)**capa
      pk(i,k)    = (pl2(i,kbot)*pk2(i,k)-ptopk)  &
                / (capap1*(pl2(i,kbot)-ptop))
      plt(i,k)   = 1000.0*pk(i,k)*pk(i,k)*pk(i,k)*sqrt(pk(i,k))
   80 continue
!c
      do 90 k = 2, lev
      ktop= kbot
      kbot= 3 - ktop
      do 90 i = 1, nx
      pl2(i,kbot)= sig(k+1)*pt(i)+ptop
      pk2(i,k)   = opok*pl2(i,kbot)**capa
      pk(i,k)    = (pl2(i,kbot)*pk2(i,k)-pl2(i,ktop)*pk2(i,k-1))  &
                / (capap1*(pl2(i,kbot)-pl2(i,ktop)))
      plt(i,k)   = 1000.0*pk(i,k)*pk(i,k)*pk(i,k)*sqrt(pk(i,k))
   90 continue
!c
      return
      end

      subroutine phi2pt(nx,my,lmax,zz,anlslp,t1000,puvphi,sgeo,pt)
!c
      use const, only: RTYPE
!c
      dimension sgeo(nx,my),pt(nx,my),zz(nx,my,lmax)
      dimension fld1(nx,lmax+2,my),fld2(nx,lmax+2,my),t1000(nx,my) &
             , anlslp(nx,my)                                       & 
             , hld1(nx,my),hld2(nx,my)                             &
             , tens(lmax+2)
      dimension phistd(lmax),puvphi(lmax),pdiff(nx,my)
      real(kind=RTYPE) presp(nx,lmax+2,my)
!c
      data cp/1004.24/, grav/9.80616/
!c
      capa= 1.0/3.5
      rgas= capa*cp
!c
      lmaxp1= lmax + 1
      lmaxp2= lmax + 2
      do 20 k = 1, lmaxp2
      tens(k) = 1.0
  20  continue
      tens(lmaxp2) = 0.0
      tens(1) = 0.0
!c
      call geostd ( lmax,puvphi,phistd )
!c
!c  read geopotential height
!c
      do 60 k=1,lmax
      pstdt= log(puvphi(k))
      do 50 j=1,my
      do 50 i=1,nx
      presp(i,k+1,j)= pstdt
      fld2(i,k+1,j)= -grav*zz(i,j,k)
      fld1(i,k+1,j)= -grav*phistd(k)-fld2(i,k+1,j)
  50  continue
  60  continue
!c
!c  initialization of terrain pressure
!c
!c  define bottom & top boundary condition for pt interpolation
!c
      l1= 1
      alaps= 6.5e-4
      eps = 0.001
      call ttstd ( 1, 1.0, tp00 )
      call ttstd ( 1, puvphi(1), tp01 )
!c
      do 130 j=1,my
      do 132 i=1,nx
      hld1(i,j)= anlslp(i,j)
      pt(i,j)  = 0.1
!c
!c dealing with inconsistency between anlslp and 1000hPa height
!c 2009/12/10
      if( fld2(i,lmaxp1,j).gt.0 .and. anlslp(i,j).ge.1000.)then
        pxx = log(hld1(i,j))
        pt(i,j) = fld2(i,lmaxp1,j)  &
               +rgas*t1000(i,j)*(pxx-presp(i,lmaxp1,j))
      endif
!c
      if(sgeo(i,j).le.-1.0) then
        px= sgeo(i,j)/(rgas*(t1000(i,j)-alaps*0.5*sgeo(i,j)))
        hld1(i,j)= hld1(i,j)*(1.0-px)
        pt(i,j)  =-sgeo(i,j) + 0.1
      endif
  132 continue
      call geostd (nx,hld1(1,j),hld2(1,j))
!c
      do 133 i=1,nx
      hld1(i,j) = log(hld1(i,j))
!ccc   if ( hld1(i,j) .le. (presp(i,lmaxp1,j)+eps) )  then
      if ( hld1(i,j) .le. presp(i,lmaxp1,j) )  then
        presp(i,lmaxp2,j)= 2.0*presp(i,lmaxp1,j) - presp(i,lmax,j)
        fld2(i,lmaxp2,j) = 2.0*fld2(i,lmaxp1,j) - fld2(i,lmax,j)
        fld1(i,lmaxp2,j) = fld1(i,lmaxp1,j)
      else
        presp(i,lmaxp2,j)= hld1(i,j)
        fld2(i,lmaxp2,j) = pt(i,j)
        fld1(i,lmaxp2,j) =-fld2(i,lmaxp2,j)-grav*hld2(i,j)
      endif
      presp(i,1,j)= log(1.0)
      fld2(i,1,j) = fld2(i,2,j) - rgas*tp00*(presp(i,2,j)-presp(i,1,j))
      fld1(i,1,j) = fld1(i,2,j)
!c
      pdiff(i,j)  = -sgeo(i,j)
      if ((sgeo(i,j).le.0.0).and.(sgeo(i,j).ge.(-1.0+eps)))  &
       pdiff(i,j)= 0.0
!c
  133 continue
!c
!c  interpolate to get terrain pressure
!c  interpolation is of log p - cubic as function of geopotential
!c  output is log p in array pt
!c
      call vterpj( nx,lmaxp2,l1,fld2(1,1,j),presp(1,1,j),pdiff(1,j) &
                , pt(1,j),tens)
!c
      do 135 i=1,nx
      pt(i,j)= exp(pt(i,j))
  135 continue
  130 continue
!c
      return
      end
