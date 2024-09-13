!----------------------------------------------------------------------------

      subroutine ndslfv_init(lonf,latg,ntrac, &
                             coslat,wgt)
!
! initialize some constant variable/arrays for ndsl advection -hmhj
!
! program log
! 2011 02 20 : henry jaung, initiated for ndsl advection
! 2013 06 20 : henry jaung, enhance cyclic condition for extra strong wind in wam
! 2013 09 30 : henry jaung, add option to advect between regular and
!                           reduced grids.
!
      use physcons
      use rank
      use index
      use grid
      use const, only : ndsladvh2, RTYPE
!
      implicit none

      real(kind=RTYPE),intent(in):: coslat(latg)
!      real(kind=RTYPE)   ,intent(in):: colrad(latg/2)
      real(kind=RTYPE),intent(in):: wgt   (latg/2)
!      integer,intent(in):: lats_nodes_a(nsize), lonf, latg, ntrac
      integer,intent(in)::  lonf, latg, ntrac
!
      integer	jm2,jm,jmh,i,j
      real(kind=RTYPE) 	pi,hfpi,twopi
      real(kind=RTYPE), dimension(:), allocatable ::  gglat,ggfact

!      logical   lprint
!
      integer	n,nm,nr,lonp
!
!      lprint = .false.
      pi = 4.0*atan(1.0)
      hfpi  = pi * 0.5
      twopi = pi * 2.0
      jm  = latg
      jm2 = jm*2
      jmh = jm/2
! ------------------------------------------------------------
! ------------------- gaussian latitude
! cosglat   cos(lat)
! gglat     gaussian latitude count from 0 to 2 pi
!           from north pole to south pole back to north pole
! gglati    gaussian latitude at edges for gglat
! gslati    cos(lat) at gglati
! ggfact    interpolation factor from gglat to gglati for wind
!
! --------------------- for parallel --------------------
!      allocate ( gslati(jm2+1) )
      allocate ( gglat(jm2),ggfact(jm2) )
!
      do j=1,jmh
        gglat(      j) = 0.5*pi - acos(coslat(j))
        gglat(jm +1-j) = 0.5*pi + acos(coslat(j))
        gglat(jm +  j) = 1.5*pi - acos(coslat(j))
        gglat(jm2+1-j) = 1.5*pi + acos(coslat(j))
      enddo
!
      gslati(    1)=0.0
      gslati(jm +1)=2.0
      gslati(jm2+1)=4.0
      do j=1,jmh
        gslati(    j+1) = gslati(    j  ) + wgt(j)
        gslati(jm -j+1) = gslati(jm -j+2) - wgt(j)
        gslati(jm +j+1) = gslati(jm +j  ) + wgt(j)
        gslati(jm2-j+1) = gslati(jm2-j+2) - wgt(j)
      enddo

!
! real(kind=RTYPE) latitude values first at edge, use temporary ggfact
      gglati(1) = hfpi
      do j=1,jmh
        gglati(j+1) = asin( sin(gglati(j)) - wgt(j) )
        ggfact(j) = gglati(j) - gglati(j+1)
      enddo
!
! then relative latitude from 0 to 2pi
      gglati(    1)=0.0
      gglati(jm +1)=pi
      gglati(jm2+1)=twopi
      do j=1,jmh
        gglati(    j+1) = gglati(    j  ) + ggfact(j)
        gglati(jm -j+1) = gglati(jm -j+2) - ggfact(j)
        gglati(jm +j+1) = gglati(jm +j  ) + ggfact(j)
        gglati(jm2-j+1) = gglati(jm2-j+2) - ggfact(j)
      enddo

!
! coefficient of Gaussian latitude interpolation
      do j=3,jm2-1
        fa1(j)=((gglati(j)-gglat(j-1))*(gglati(j)-gglat(j))       & 
              *(gglati(j)-gglat(j+1)))/((gglat(j-2)-gglat(j-1))   &
              *(gglat(j-2)-gglat(j))*(gglat(j-2)-gglat(j+1)))
        fa2(j)=((gglati(j)-gglat(j-2))*(gglati(j)-gglat(j))       &
              *(gglati(j)-gglat(j+1)))/((gglat(j-1)-gglat(j-2))   &
              *(gglat(j-1)-gglat(j))*(gglat(j-1)-gglat(j+1)))
        fa3(j)=((gglati(j)-gglat(j-2))*(gglati(j)-gglat(j-1))     &
              *(gglati(j)-gglat(j+1)))/((gglat(j)-gglat(j-2))     &
              *(gglat(j)-gglat(j-1))*(gglat(j)-gglat(j+1)))
        fa4(j)=((gglati(j)-gglat(j-2))*(gglati(j)-gglat(j-1))     &
              *(gglati(j)-gglat(j)))/((gglat(j+1)-gglat(j-2))     &
              *(gglat(j+1)-gglat(j-1))*(gglat(j+1)-gglat(j)))
      enddo
      ! over pole
      fa1(2)=((gglati(2)-gglat(1))*(gglati(2)-gglat(2))           &
            *(gglati(2)-gglat(3)))/((-1.*gglat(1)-gglat(1))       &
            *(-1.*gglat(1)-gglat(2))*(-1.*gglat(1)-gglat(3)))
      fa4(jm2)=fa1(2)

      fa2(2)=((gglati(2)+gglat(1))*(gglati(2)-gglat(2))           &
            *(gglati(2)-gglat(3)))/((gglat(1)+gglat(1))           &
            *(gglat(1)-gglat(2))*(gglat(1)-gglat(3)))
      fa3(jm2)=fa2(2)

      fa3(2)=((gglati(2)+gglat(1))*(gglati(2)-gglat(1))           &
            *(gglati(2)-gglat(3)))/((gglat(2)+gglat(1))           &
            *(gglat(2)-gglat(1))*(gglat(2)-gglat(3)))
      fa2(jm2)=fa3(2)

      fa4(2)=((gglati(2)+gglat(1))*(gglati(2)-gglat(1))           &
            *(gglati(2)-gglat(2)))/((gglat(3)+gglat(1))           &
            *(gglat(3)-gglat(1))*(gglat(3)-gglat(2)))
      fa1(jm2)=fa4(2)

      fa1(1)=((gglati(1)+gglat(1))*(gglati(1)-gglat(1))           &
            *(gglati(1)-gglat(2)))/((gglat(jm2-1)-gglat(jm2))     &
            *(-1.*gglat(2)-gglat(1))*(-1.*gglat(2)-gglat(2)))
      fa2(1)=((gglati(1)+gglat(2))*(gglati(1)-gglat(1))           &
            *(gglati(1)-gglat(2)))/((gglat(jm2)-gglat(jm2-1))     &
            *(-1.*gglat(1)-gglat(1))*(-1.*gglat(1)-gglat(2)))
      fa3(1)=((gglati(1)+gglat(2))*(gglati(1)+gglat(1))           &
            *(gglati(1)-gglat(2)))/((gglat(1)+gglat(2))           &
            *(gglat(1)+gglat(1))*(gglat(1)-gglat(2)))
      fa4(1)=((gglati(1)+gglat(2))*(gglati(1)+gglat(1))           &
            *(gglati(1)-gglat(1)))/((gglat(2)+gglat(2))           &
            *(gglat(2)+gglat(1))*(gglat(2)-gglat(1)))

      deallocate(gglat,ggfact)
!
!      allocate( lonstr(nsize), lonlen(nsize) )
!      allocate( latstr(nsize), latlen(nsize) )
!
      lonfull = lonf
      lonhalf = lonf / 2	! lonf has to be even
      lonpart = (lonhalf-1)/nsizey+1
!
      lonp = lonhalf / nsizey
      do n=1,nsizey
        lonlen(n)=lonp
      enddo
      nr=mod(lonhalf,nsizey)
      if( nr.ne.0 ) then
        do n=1,nr
          lonlen(n)=lonlen(n)+1
        enddo
      endif
      nm=1
      do n=1,nsizey
        lonstr(n) = nm
        nm = nm + lonlen(n)
!       print *,' node lonstr lonlen ',n,lonstr(n),lonlen(n)
      enddo
!!
      latfull = my * 2
      lathalf = my
      latpart = 0
      do n=1,nsizey
        latpart=max(latpart,jlistnum_sl(n))
        latlen(n) = jlistnum_sl(n)
      enddo
!      if ( myrank .eq. 0 ) then
!        print *,'latpart & my_max : ',latpart,' & ',my_max
!        do n=1,nsize
!         print *,'node=',n,' jlistnum_sl=',jlistnum_sl(n)
!        enddo
!      endif

      nm=1
      do n=1,nsizey
        latstr(n) = nm
        nm = nm + latlen(n)
!       print *,' node latstr latlen ',n,latstr(n),latlen(n)
      enddo
!
      lonlenmax=0
      latlenmax=0
      do n=1,nsizey
        lonlenmax = max(lonlenmax,lonlen(n))
        latlenmax = max(latlenmax,latlen(n))
      enddo
!
      mylonlen = lonlen(col_rank+1)
      mylatlen = latlen(col_rank+1)
!
!!      ndslhvar = 4 + ntrac 	! u,v,t,ps,tracers
      ndslhvar = 3 + ntrac 	! u,v,t,tracers
      ndslvvar = 3 + ntrac	! u,v,t,tracers
!      print *,' ndslfv_init: total variables for ndsl hadv =',ndslhvar
!      print *,' ndslfv_init: total variables for ndsl vadv =',ndslvvar
!  define sequence for first integration in NDSL
      if ( ndsladvh2 ) then
        xy = 0
      else
        xy = 1
      endif
!!
      return
      end subroutine ndslfv_init

! -------------------------------------------------------------------------------
      subroutine cyclic_cell_massadvx(im,imf,levs,nvars,delt,uc,qq,mass)
!
! compute local positive advection with mass conservation
! qq is advected by uc from past to next position
!
! author: hann-ming henry juang 2008
!
!      use grid
!      use gfs_dyn_layout1
!
      use const, only : RTYPE

      implicit none
!
      integer	im,imf,levs,nvars,mass
      real(kind=RTYPE)	delt,pi
      real(kind=RTYPE)	uc(imf,levs)
      real(kind=RTYPE)	qq(imf,levs,nvars)
!
      real(kind=RTYPE)	past(im,nvars),next(im,nvars),da(im,nvars)
      real(kind=RTYPE)	dxfact(im)
      real(kind=RTYPE)	xreg(im+1),xpast(im+1),xnext(im+1)
      real(kind=RTYPE)	uint(im+1)
      real(kind=RTYPE) 	dist(im+1),sc,ds(im+1),step(10),dist_step
      real(kind=RTYPE), parameter :: fa1 = 9./16.
      real(kind=RTYPE), parameter :: fa2 = 1./16.

      integer  	i,k,n,nn,nf,nv,nst,nstep

!     sc = ggloni(im+1)-ggloni(1)
      pi = 4.0 * atan(1.0)
      sc = 2.0 * pi
      do i=1,im+1
        ds(i) = sc / float(im)
      enddo
      do i=1,im+1
        xreg(i) = (i-1.5)*ds(i)
      enddo
      nv = nvars
!
      do k=1,levs
!
! 4th order interpolation from mid point to cell interfaces
!
        do i=3,im-1
          uint(i)=fa1*(uc(i,k)+uc(i-1,k))-fa2*(uc(i+1,k)+uc(i-2,k))
        enddo
        uint(2)=fa1*(uc(2,k)+uc(1 ,k))-fa2*(uc(3,k)+uc(im  ,k))
        uint(1)=fa1*(uc(1,k)+uc(im,k))-fa2*(uc(2,k)+uc(im-1,k))
        uint(im+1)=uint(1)
        uint(im  )=fa1*(uc(im,k)+uc(im-1,k))                            & 
                       -fa2*(uc(1,k)+uc(im-2,k))

!
! compute past and next positions of cell interfaces
!
        do i=1,im+1
          dist(i)  = uint(i) * delt
        enddo
!cflx   call def_cfl_max (im+1,dist,ds,nstep)
        call def_cfl_step(im+1,dist,ds,step,nstep,levs+1-k,'advx')
!
!  mass positive advection
!
       do nst = 1, nstep
!
        do i=1,im+1
          dist_step = dist(i)*step(nst)
          xpast(i) = xreg(i) - dist_step
          xnext(i) = xreg(i) + dist_step
        enddo
        if( mass.eq.1 ) then
         do i=1,im
          dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
         enddo
        endif
!
        do n=1,nv
          past(1:im,n) = qq(1:im,k,n)
        enddo
        call cyclic_cell_ppm_intp(xreg,past,xpast,da,im,nv,im,im,sc)
!        call cyclic_cell_plm_intp(xreg,past,xpast,da,im,nv,im,im,sc)
        if( mass.eq.1) then
          do n=1,nv
            da(1:im,n) = da(1:im,n) * dxfact(1:im)
          enddo
        endif
        call cyclic_cell_ppm_intp(xnext,da,xreg,next,im,nv,im,im,sc)
!        call cyclic_cell_plm_intp(xnext,da,xreg,next,im,nv,im,im,sc)
        do n=1,nv
          qq(1:im,k,n) = next(1:im,n)
        enddo
!
       enddo
!
      enddo

      return
      end subroutine cyclic_cell_massadvx
!
!
! -------------------------------------------------------------------------------
      subroutine cyclic_cell_massadvxl(im,imf,levs,nvars,delt,uc,qq,mass,forward)
!
! compute local positive advection with mass conservation
! qq is advected by uc from past to next position
!
! author: hann-ming henry juang 2008
!
!      use grid
!      use gfs_dyn_layout1
!
      use const, only : RTYPE

      implicit none
!
      integer	im,imf,levs,nvars,mass
      real(kind=RTYPE)	delt,pi
      real(kind=RTYPE)	uc(imf,levs)
      real(kind=RTYPE)	qq(imf,levs,nvars)
!
      real(kind=RTYPE)	past(im,nvars),next(im,nvars),da(im,nvars)
      real(kind=RTYPE)	dxfact(im)
      real(kind=RTYPE)	xreg(im+1),xpast(im+1),xnext(im+1)
      real(kind=RTYPE)	uint(im+1)
      real(kind=RTYPE) 	dist(im+1),sc,ds(im+1),step(10),dist_step
      real(kind=RTYPE), parameter :: fa1 = 9./16.
      real(kind=RTYPE), parameter :: fa2 = 1./16.

      integer  	i,k,n,nn,nf,nv,nst,nstep
      logical   forward

!     sc = ggloni(im+1)-ggloni(1)
      pi = 4.0 * atan(1.0)
      sc = 2.0 * pi
      do i=1,im+1
        ds(i) = sc / float(im)
      enddo
      do i=1,im+1
        xreg(i) = (i-1.5)*ds(i)
      enddo
      nv = nvars
!
      do k=1,levs
!
! 4th order interpolation from mid point to cell interfaces
!
        do i=3,im-1
          uint(i)=fa1*(uc(i,k)+uc(i-1,k))-fa2*(uc(i+1,k)+uc(i-2,k))
        enddo
        uint(2)=fa1*(uc(2,k)+uc(1 ,k))-fa2*(uc(3,k)+uc(im  ,k))
        uint(1)=fa1*(uc(1,k)+uc(im,k))-fa2*(uc(2,k)+uc(im-1,k))
        uint(im+1)=uint(1)
        uint(im  )=fa1*(uc(im,k)+uc(im-1,k))                            & 
                       -fa2*(uc(1,k)+uc(im-2,k))

!
! compute past and next positions of cell interfaces
!
        do i=1,im+1
          dist(i)  = uint(i) * delt
        enddo
!cflx   call def_cfl_max (im+1,dist,ds,nstep)
        call def_cfl_step(im+1,dist,ds,step,nstep,levs+1-k,'advx')
!
!  mass positive advection
!
       do nst = 1, nstep
!
        if ( forward ) then
          do i=1,im+1
            dist_step = dist(i)*step(nst)
            xpast(i) = xreg(i)
            xnext(i) = xreg(i) + dist_step
          enddo
        else
          do i=1,im+1
            dist_step = dist(i)*step(nst)
            xpast(i) = xreg(i) - dist_step
            xnext(i) = xreg(i) + dist_step
          enddo
        endif
        if( mass.eq.1 ) then
         do i=1,im
          dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
         enddo
        endif
!
        if ( forward ) then
          do n=1,nv
            da(1:im,n) = qq(1:im,k,n)
          enddo
        else
          do n=1,nv
            past(1:im,n) = qq(1:im,k,n)
          enddo
          call cyclic_cell_ppm_intp(xreg,past,xpast,da,im,nv,im,im,sc)
!          call cyclic_cell_plm_intp(xreg,past,xpast,da,im,nv,im,im,sc)
        endif
        if( mass.eq.1) then
          do n=1,nv
            da(1:im,n) = da(1:im,n) * dxfact(1:im)
          enddo
        endif
        call cyclic_cell_ppm_intp(xnext,da,xreg,next,im,nv,im,im,sc)
!        call cyclic_cell_plm_intp(xnext,da,xreg,next,im,nv,im,im,sc)
        do n=1,nv
          qq(1:im,k,n) = next(1:im,n)
        enddo
!
       enddo
!
      enddo

      return
      end subroutine cyclic_cell_massadvxl
!
!
!-------------------------------------------------------------------
      subroutine cyclic_cell_massadvy(jm,lev,nvars,delt,vc,qq,mass)
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
!
      use grid     , only : gglati,fa1,fa2,fa3,fa4
!
      use const, only : RTYPE

      implicit none
!
      integer   jm,lev,nvars,mass
      real(kind=RTYPE)      delt
      real(kind=RTYPE)      vc(jm,lev)
      real(kind=RTYPE)      qq(jm,lev,nvars)
!
      real(kind=RTYPE)      var(jm)
      real(kind=RTYPE)      past(jm,nvars),da(jm,nvars),next(jm,nvars)
      real(kind=RTYPE)      dyfact(jm)
      real(kind=RTYPE)      ypast(jm+1),ynext(jm+1)
      real(kind=RTYPE)      dist (jm+1), ds(jm), step(10), dist_step
      real(kind=RTYPE)      sc

      integer   n,k,j,jmh,nv,nst,nstep
!
! preparations ---------------------------
!
      jmh  = jm / 2
      sc = gglati(jm+1)-gglati(1)
      do j=1,jm
        ds(j) = gglati(j+1) - gglati(j)
      enddo
      nv   = nvars
!
      do k=1,lev
!
        do j=1,jmh
          var(j)     =  vc(j    ,k) * delt
          var(j+jmh) = -vc(j+jmh,k) * delt
        enddo
! for Gaussian latitude
        do j=3,jm-1
          dist(j)=fa1(j)*var(j-2)+fa2(j)*var(j-1)+fa3(j)*var(j)+fa4(j)*var(j+1)
        enddo
        ! over pole
        dist(2)=fa1(2)*var(jm)+fa2(2)*var(1)+fa3(2)*var(2)+fa4(2)*var(3)
        dist(jm)=fa1(jm)*var(jm-2)+fa2(jm)*var(jm-1)+fa3(jm)*var(jm)+fa4(jm)*var(1)

        dist(1)=fa1(1)*var(jm-1)+fa2(1)*var(jm)+fa3(1)*var(1)+fa4(1)*var(2)
        dist(jm+1)=dist(1)
! for reguler grid
!        do j=3,jm-1
!          dist(j)=fa1*(var(j)+var(j-1))-fa2*(var(j+1)+var(j-2))
!        enddo
!        ! over pole
!        dist(2)=fa1*(var(2)+var(1 ))-fa2*(var(3)+var(jm  ))
!        dist(1)=fa1*(var(1)+var(jm))-fa2*(var(2)+var(jm-1))
!        dist(jm+1)=dist(1)
!        dist(jm  )=fa1*(var(jm)+var(jm-1))-fa2*(var(1)+var(jm-2))
!cflx   call def_cfl_max (jm+1,dist,ds,nstep)
        call def_cfl_step(jm+1,dist,ds,step,nstep,lev+1-k,'advy')

!
! advection all in y
!
       do nst = 1, nstep
!
        do j=1,jm+1
          dist_step = dist(j)*step(nst)
          ypast(j) = gglati(j) - dist_step
          ynext(j) = gglati(j) + dist_step
        enddo
        if( mass.eq.1 ) then
         do j=1,jm
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
         enddo
        endif

        do n=1,nv
          past(1:jm,n) = qq(1:jm,k,n)
        enddo
        call cyclic_cell_ppm_intp(gglati,past,ypast,da,jm,nv,jm,jm,sc)
!        call cyclic_cell_plm_intp(gglati,past,ypast,da,jm,nv,jm,jm,sc)

        if( mass.eq.1 ) then
          do n=1,nv
            da(1:jm,n) = da(1:jm,n) * dyfact(1:jm)
          enddo
        endif
        call cyclic_cell_ppm_intp(ynext,da,gglati,next,jm,nv,jm,jm,sc)
!        call cyclic_cell_plm_intp(ynext,da,gglati,next,jm,nv,jm,jm,sc)

        do n=1,nv
          qq(1:jm,k,n) = next(1:jm,n)
        enddo
!
       enddo
!
      enddo

      return
      end subroutine cyclic_cell_massadvy
!
!-------------------------------------------------------------------
      subroutine cyclic_cell_massadvyl(jm,lev,nvars,delt,vc,qq,mass,forward)
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
!
      use grid     , only : gglati,fa1,fa2,fa3,fa4
!
      use const, only : RTYPE

      implicit none
!
      integer   jm,lev,nvars,mass
      real(kind=RTYPE)      delt
      real(kind=RTYPE)      vc(jm,lev)
      real(kind=RTYPE)      qq(jm,lev,nvars)
!
      real(kind=RTYPE)      var(jm)
      real(kind=RTYPE)      past(jm,nvars),da(jm,nvars),next(jm,nvars)
      real(kind=RTYPE)      dyfact(jm)
      real(kind=RTYPE)      ypast(jm+1),ynext(jm+1)
      real(kind=RTYPE)      dist (jm+1), ds(jm), step(10), dist_step
      real(kind=RTYPE)      sc

      integer   n,k,j,jmh,nv,nst,nstep
      logical   forward
!
! preparations ---------------------------
!
      jmh  = jm / 2
      sc = gglati(jm+1)-gglati(1)
      do j=1,jm
        ds(j) = gglati(j+1) - gglati(j)
      enddo
      nv   = nvars
!
      do k=1,lev
!
        do j=1,jmh
          var(j)     =  vc(j    ,k) * delt
          var(j+jmh) = -vc(j+jmh,k) * delt
        enddo
! for Gaussian latitude
        do j=3,jm-1
          dist(j)=fa1(j)*var(j-2)+fa2(j)*var(j-1)+fa3(j)*var(j)+fa4(j)*var(j+1)
        enddo
        ! over pole
        dist(2)=fa1(2)*var(jm)+fa2(2)*var(1)+fa3(2)*var(2)+fa4(2)*var(3)
        dist(jm)=fa1(jm)*var(jm-2)+fa2(jm)*var(jm-1)+fa3(jm)*var(jm)+fa4(jm)*var(1)

        dist(1)=fa1(1)*var(jm-1)+fa2(1)*var(jm)+fa3(1)*var(1)+fa4(1)*var(2)
        dist(jm+1)=dist(1)
! for equal latitude
!        do j=3,jm-1
!          dist(j)=fa1*(var(j)+var(j-1))-fa2*(var(j+1)+var(j-2))
!        enddo
! over pole
!        dist(2)=fa1*(var(2)+var(1 ))-fa2*(var(3)+var(jm  ))
!        dist(1)=fa1*(var(1)+var(jm))-fa2*(var(2)+var(jm-1))
!        dist(jm+1)=dist(1)
!        dist(jm  )=fa1*(var(jm)+var(jm-1))-fa2*(var(1)+var(jm-2))
!cflx   call def_cfl_max (jm+1,dist,ds,nstep)
        call def_cfl_step(jm+1,dist,ds,step,nstep,lev+1-k,'advy')

!
! advection all in y
!
       do nst = 1, nstep
!
        if ( forward ) then
          do j=1,jm+1
            dist_step = dist(j)*step(nst)
            ypast(j) = gglati(j)
            ynext(j) = gglati(j) + dist_step
          enddo
        else
          do j=1,jm+1
            dist_step = dist(j)*step(nst)
            ypast(j) = gglati(j) - dist_step
            ynext(j) = gglati(j) + dist_step
          enddo
        endif
        if( mass.eq.1 ) then
         do j=1,jm
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
         enddo
        endif

        if ( forward ) then
          do n=1,nv
            da(1:jm,n) = qq(1:jm,k,n)
          enddo
        else
          do n=1,nv
            past(1:jm,n) = qq(1:jm,k,n)
          enddo
          call cyclic_cell_ppm_intp(gglati,past,ypast,da,jm,nv,jm,jm,sc)
!          call cyclic_cell_plm_intp(gglati,past,ypast,da,jm,nv,jm,jm,sc)
        endif

        if( mass.eq.1 ) then
          do n=1,nv
            da(1:jm,n) = da(1:jm,n) * dyfact(1:jm)
          enddo
        endif
        call cyclic_cell_ppm_intp(ynext,da,gglati,next,jm,nv,jm,jm,sc)
!        call cyclic_cell_plm_intp(ynext,da,gglati,next,jm,nv,jm,jm,sc)

        do n=1,nv
          qq(1:jm,k,n) = next(1:jm,n)
        enddo
!
       enddo
!
      enddo

      return
      end subroutine cyclic_cell_massadvyl
!
!-------------------------------------------------------------------
      subroutine fixend_cell_massadvy(jm,jmh,levs,nvars,delt,vc,qq,mass)
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
      use grid     , only : gslati
!      use gfs_dyn_layout1
!
      use const, only : RTYPE

      implicit none
!
      integer   jm,jmh,levs,nvars,mass
      real(kind=RTYPE)	delt
      real(kind=RTYPE)	vc(jm,levs)
      real(kind=RTYPE)	qq(jm,levs,nvars)
!
      real(kind=RTYPE)	var(jmh)
      real(kind=RTYPE)	past(jmh,nvars),da(jmh,nvars),next(jmh,nvars)
      real(kind=RTYPE)	dyfact(jmh)
      real(kind=RTYPE)	ypast(jmh+1),ynext(jmh+1)
      real(kind=RTYPE)	dist (jmh+1), ds(jmh)
      real(kind=RTYPE)	step(10), dist_step
      real(kind=RTYPE) 	hfpi,pi
      real(kind=RTYPE), parameter :: fa1 = 9./16.
      real(kind=RTYPE), parameter :: fa2 = 1./16.

      integer  	n,k,j,nv,nst,nstep
!
! preparations ---------------------------
!
      do j=1,jmh
        ds(j) = gslati(j+1) - gslati(j)
      enddo
      pi = 4.0 * atan(1.0)
      hfpi = pi * 0.5
      nv   = nvars
!
      do k=1,levs
!
! first hemisphere
        do j=1,jmh
          var(j)     = vc(j    ,k) * delt
        enddo
        do j=2,jmh
          dist(j)=(var(j)*ds(j-1)+var(j-1)*ds(j)) / (ds(j-1)+ds(j))
        enddo
        dist(1)=0.0
        dist(jmh+1)=0.0
!cflx   call def_cfl_max (jmh+1,dist,ds,nstep)
        call def_cfl_step(jmh+1,dist,ds,step,nstep,levs+1-k,'advy')
!
        do nst = 1, nstep
!
        do j=1,jmh+1
          dist_step = dist(j)*step(nst)
          ypast(j) = gslati(j) - dist_step
          ynext(j) = gslati(j) + dist_step
        enddo
        if( mass.eq.1 ) then
         do j=1,jmh
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
         enddo
        endif
        do n=1,nv
          past(1:jmh,n) = qq(1:jmh,k,n)
        enddo
        call fixend_cell_plm_intp(gslati,past,ypast,da,jmh,nv)
        if( mass.eq.1 ) then
          do n=1,nv
            da(1:jmh,n) = da(1:jmh,n) * dyfact(1:jmh)
          enddo
        endif
        call fixend_cell_plm_intp(ynext,da,gslati,next,jmh,nv)	
        do n=1,nv
          qq(1:jmh,k,n) = next(1:jmh,n)
        enddo

        enddo   ! for nst
!
! second hemisphere
        do j=1,jmh
          var(j) =  -vc(j+jmh,k) * delt
        enddo
        do j=2,jmh
          dist(j)=(var(j)*ds(j-1)+var(j-1)*ds(j)) / (ds(j-1)+ds(j))
        enddo
        dist(1)=0.0
        dist(jmh+1)=0.0
!cflx   call def_cfl_max (jmh+1,dist,ds,nstep)
        call def_cfl_step(jmh+1,dist,ds,step,nstep,levs+1-k,'advy')
!
        do nst = 1, nstep
!
        do j=1,jmh+1
          dist_step = dist(j)*step(nst)
          ypast(j) = gslati(j) - dist_step
          ynext(j) = gslati(j) + dist_step
        enddo
        if( mass.eq.1 ) then
         do j=1,jmh
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
         enddo
        endif

        do n=1,nv
          past(1:jmh,n) = qq(jmh+1:jm,k,n)
        enddo
        call fixend_cell_plm_intp(gslati,past,ypast,da,jmh,nv)
        if( mass.eq.1 ) then
          do n=1,nv
            da(1:jmh,n) = da(1:jmh,n) * dyfact(1:jmh)
          enddo
        endif
        call fixend_cell_plm_intp(ynext,da,gslati,next,jmh,nv)	
        do n=1,nv
          qq(jmh+1:jm,k,n) = next(1:jmh,n)
        enddo
!
        enddo    ! for nst
!
      enddo ! for k

      return
      end subroutine fixend_cell_massadvy
!
!-------------------------------------------------------------------------------
      subroutine cyclic_cell_intpx(levs,imp,imf,qq)
!
! do  mass conserving interpolation from different grid at given latitude
!
! author: hann-ming henry juang 2008
!
      use grid      , only : lonfull
!     use gfs_dyn_layout1
      use const, only : RTYPE

      implicit none
!
      integer	 levs, imp, imf
      real(kind=RTYPE)	 qq(lonfull,levs)
!
      real(kind=RTYPE)	old(lonfull,levs),new(lonfull,levs)
      real(kind=RTYPE)	xpast(lonfull+1),xnext(lonfull+1)
      real(kind=RTYPE)	two_pi,dxp,dxf,hfdxp,hfdxf,sc,pi
!
      integer  	i,k,im
!
      im = lonfull

! ..................................
      if( imp.ne.imf ) then
! ..................................
        pi  = 4.0 * atan(1.0)
        two_pi = 2.0 * pi
        dxp = two_pi / imp
        dxf = two_pi / imf
        hfdxp = 0.5 * dxp
        hfdxf = 0.5 * dxf

        do i=1,imp+1
          xpast(i) = (i-1) * dxp - hfdxp
        enddo

        do i=1,imf+1
          xnext(i) = (i-1) * dxf - hfdxf
        enddo

        sc=two_pi

        old(1:imp,1:levs)=qq(1:imp,1:levs)
        call cyclic_cell_ppm_intp(xpast,old,xnext,new,im,levs,imp,imf,sc)
!        call cyclic_cell_plm_intp(xpast,old,xnext,new,im,levs,imp,imf,sc)
      	
        qq(1:imf,1:levs)=new(1:imf,1:levs)

! .................
      endif
! .................

      return
      end subroutine cyclic_cell_intpx
!
!-------------------------------------------------------------------------------
      subroutine cyclic_cell_intpxl(levs,imp,imf,qq)
!
! do  mass conserving interpolation from different grid at given latitude
!
! author: hann-ming henry juang 2008
!
      use grid      , only : lonfull
!      use gfs_dyn_layout1
      use const      , only : RTYPE

      implicit none
!
      integer	 levs, imp, imf
      real(kind=RTYPE)	 qq(lonfull,levs)
!
      real(kind=RTYPE)	old(lonfull,levs),new(lonfull,levs)
      real(kind=RTYPE)	xpast(lonfull+1),xnext(lonfull+1)
      real(kind=RTYPE)	two_pi,dxp,dxf,hfdxp,hfdxf,sc,pi
!
      integer  	i,k,im
!
      im = lonfull

! ..................................
      if( imp.ne.imf ) then
! ..................................
        pi  = 4.0 * atan(1.0)
        two_pi = 2.0 * pi
        dxp = two_pi / imp
        dxf = two_pi / imf
        hfdxp = 0.5 * dxp
        hfdxf = 0.5 * dxf

        do i=1,imp+1
          xpast(i) = (i-1) * dxp - hfdxp
        enddo

        do i=1,imf+1
          xnext(i) = (i-1) * dxf - hfdxf
        enddo

        sc=two_pi

        old(1:imp,1:levs)=qq(1:imp,1:levs)
!        call cyclic_cell_ppm_intp(xpast,old,xnext,new,im,levs,imp,imf,sc)
        call cyclic_cell_plm_intp(xpast,old,xnext,new,im,levs,imp,imf,sc)
      
        qq(1:imf,1:levs)=new(1:imf,1:levs)

! .................
      endif
! .................

      return
      end subroutine cyclic_cell_intpxl
!
! -------------------------------------------------------------------------
      subroutine cyclic_cell_plm_intp(pp,qq,pn,qn,lons,nv,lonp,lonn,sc)
!
! linear quick version of cyclic cell interpolation
!
! pp    location at interfac point as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
! lonp  numer of cells for input
! lonn  numer of cells for output
! sc    length of the entire circle
!
! author : henry.juang@noaa.gov
!
!      use gfs_dyn_layout1
!
      use const      , only : RTYPE

      implicit none
!
      real(kind=RTYPE)      pp(lons+1)
      real(kind=RTYPE)      qq(lons  ,nv)
      real(kind=RTYPE)      pn(lons+1)
      real(kind=RTYPE)      qn(lons  ,nv)
      integer   lons,lonp,lonn,nv
      real(kind=RTYPE)      sc
!
      real(kind=RTYPE)      px(lonp+2),ps(lonp+2)
      real(kind=RTYPE)      hfds(lonp+1),rdsi(lonp+1)
      real(kind=RTYPE)      locs(lonp+lonn+3),dp(lonp+lonn+2),dt(lonp+lonn+2)
      real(kind=RTYPE)      qc(0:lonp+2),dq(lonp+1),qmi(lonp+1),qpi(lonp+1)
      real(kind=RTYPE)      vals(lonp+lonn+2)
      real(kind=RTYPE)      rdd(lonn)
      real(kind=RTYPE)      shift,hfsc,dd,check,ss
      integer   js(lonn+1),ix4i(lonp+2),i4j(lonp+lonn+3)
      integer   i,j,n,ip,ix,is,in
!
      hfsc = 0.5 * sc
!
! arrange input array cover output location with cyclic boundary
! condition
!
! bring input location close to output location
! and cover output starting point.
      if( pp(1).ne.pn(1) ) then
        if( pp(1).gt.pn(1) ) then
          shift = (int((pn(1)-pp(1))/sc)-1) * sc
        else
          shift =  int((pn(1)-pp(1))/sc) * sc
        endif
        do i=1,lonp+1
          px(i)=pp(i) + shift                   !bring pp to left of pn
        enddo
        px(lonp+2) = px(2) + sc
        if( pn(1)-px(1) .lt. hfsc ) then        !search start point
          forward: do i=1,lonp+1
            if( px(i).gt.pn(1) ) then
              ix = i - 1
              exit forward
            endif
          enddo forward
        else
          backward: do i=lonp,1,-1
            if( px(i).le.pn(1) ) then
              ix = i
              exit backward
            endif
          enddo backward
        endif
        i = 0
        do j=ix,lonp
          i=i+1
          ps(i) = px(j)
          ix4i(i) = j
        enddo
        do j=1,ix-1
          i=i+1
          ps(i) = px(j) + sc
          ix4i(i) = j
        enddo
        ps  (lonp+1) = ps(1) + sc
        ix4i(lonp+1) = ix
        ps  (lonp+2) = ps(2) + sc
        ix4i(lonp+2) = ix+1
      else
        do i=1,lonp+1
          ps(i) = pp(i)
        enddo
        ps(lonp+2) = ps(2) + sc
        do i=1,lonp
          ix4i(i) = i
        enddo
        ix4i(lonp+1) = 1
        ix4i(lonp+2) = 2
      endif

      locs(1) = ps(1)
      is=2
      i4j(1)=1
      in=1
      do j=2,lonp+lonn+1
        if( ps(is).lt.pn(in) ) then
          locs(j) = ps(is)
          is = is + 1
        else
          locs(j) = pn(in)
          js(in) = j
          in = in + 1
        endif
        i4j(j) = is - 1
      enddo
      locs(lonp+lonn+2) = pn(lonn+1)
      js(lonn+1)        = lonp+lonn+2
      i4j(lonp+lonn+2)  = lonp+1
      locs(lonp+lonn+3) = ps(lonp+2)
!     if( is.gt.lonp+2) then
!        print *,' ================================================= '
!        print *,' --- last j and lonp+lonn+3 sc ',j,lonp+lonn+3,sc
!        print *,' locs ',(locs(j),j=lonp+lonn,lonp+lonn+3)
!        print *,' --- last is i4j ',is,i4j(lonp+lonn+2)
!        print *,' ps ',(ps(i),i=lonp,lonp+2)
!        print *,' --- last in js ',in,js(in-1)
!        print *,' pn ',(pn(i),i=lonn,lonn+1)
!      endif
!
! interpolation coefficient
!
      do i=1,lonp+1
        hfds(i) = 0.5*(ps(i+1)-ps(i))
      enddo
      rdsi(1) = 1.0 / (hfds(1)+hfds(lonp))
      do i=2,lonp+1
        rdsi(i) = 1.0 / (hfds(i)+hfds(i-1))
      enddo
      do i=1,lonn
        dd = 0.0
        do j=js(i),js(i+1)-1
          dp(j) = locs(j+1) - locs(j)
          dd    = dd    + dp(j)
        enddo
        rdd(i) = 1.0 / dd
      enddo
      do j=1,lonp+lonn+2
        i = i4j(j)
        dt(j) = (locs(j) + 0.5*dp(j)) - (ps(i) + hfds(i))
      enddo
!
! start interpolation by integral of ppm
!
      do n=1,nv


        qc(0) = qq(ix4i(lonp),n)
        do i=1,lonp+1
          ix = ix4i(i)
          qc(i) = qq(ix,n)
        enddo

        do i=1,lonp+1
          qmi(i) = (qc(i)-qc(i-1))*rdsi(i)
        enddo
        do i=1,lonp
          qpi(i) = qmi(i+1)
        enddo
        qpi(lonp+1) = qmi(2)
        do i=1,lonp+1
          check = qmi(i)*qpi(i)
          if( check.lt.0.0 ) then
            dq(i) = 0.0
          else
            dq(i)  = 0.5*(qmi(i)+qpi(i))
          endif
        enddo

        do j=1,lonp+lonn+2
          i = i4j(j)
          vals(j) = qc(i) + dq(i)*dt(j)
        enddo

        do i=1,lonn
          ss = 0.0
          do j=js(i),js(i+1)-1
            ss = ss + vals(j)*dp(j)
          enddo
          qn(i,n) = ss * rdd(i)
        enddo

      enddo
!
      return
      end subroutine cyclic_cell_plm_intp
!
! -------------------------------------------------------------------------
      subroutine fixend_cell_plm_intp(pp,qq,pn,qn,lons,nv)
!
! linear quick version of the same two end points cell interpolation
! so pp(1) = pn(1) and pp(lonp+1) = pn(lonn+1)
!
! pp    location at interfac point as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
!
! author : henry.juang@noaa.gov
!
!      use gfs_dyn_layout1
!
      use index, only : col_rank
      use const, only : RTYPE

      implicit none
!
      real(kind=RTYPE)      pp(lons+1)
      real(kind=RTYPE)      qq(lons  ,nv)
      real(kind=RTYPE)      pn(lons+1)
      real(kind=RTYPE)      qn(lons  ,nv)
      integer   lons,nv
!
      real(kind=RTYPE)      hfdp(lons),rdsi(lons)
      real(kind=RTYPE)      locf(lons+lons),df(lons+lons-1),dt(lons+lons-1)
      real(kind=RTYPE)      qc(lons),dq(lons),qmi(lons),qpi(lons)
      real(kind=RTYPE)      valf(lons+lons-1)
      real(kind=RTYPE)      rdd(lons)
      real(kind=RTYPE)      dd,check,ss
      integer   js(lons+1),i4j(lons+lons)
      integer   i,j,n,ip,in
!
! arrange input array cover output location with cyclic boundary
! condition
!
      locf(1) = pp(1)
      i4j(1)=1
      js(1)=1
      ip=2
      in=2
      do j=2,lons+lons-1
        if( pp(ip).le.pn(in) ) then
          locf(j) = pp(ip)
          ip = ip + 1
        else
          locf(j) = pn(in)
          js(in) = j
          in = in + 1
        endif
        i4j(j) = ip - 1
      enddo
      locf(lons+lons) = pp(lons+1)
      js(lons+1) = lons+lons
!
! interpolation coefficient
!
      do i=1,lons
        hfdp(i) = 0.5*(pp(i+1)-pp(i))
      enddo
      do i=2,lons
        rdsi(i) = 1.0 / (hfdp(i)+hfdp(i-1))
      enddo
      do i=1,lons
        dd = 0.0
        do j=js(i),js(i+1)-1
          df(j) = locf(j+1) - locf(j)
          dd    = dd    + df(j)
        enddo
        rdd(i) = 1.0 / dd
      enddo
      do j=1,lons+lons-1
        i = i4j(j)
        dt(j) = (locf(j) + 0.5*df(j)) - (pp(i) + hfdp(i))
      enddo
!
! start interpolation by integral of ppm
!
      dq(1) = 0.0
      dq(lons) = 0.0

      do n=1,nv

        do i=1,lons
          qc(i) = qq(i,n)
        enddo
        do i=2,lons
          qmi(i) = (qc(i)-qc(i-1))*rdsi(i)
        enddo
        do i=2,lons-1
          qpi(i) = qmi(i+1)
        enddo
        do i=2,lons-1
          check = qmi(i)*qpi(i)
          if( check.lt.0.0 ) then
            dq(i) = 0.0
          else
            dq(i)  = 0.5*(qmi(i)+qpi(i))
          endif
        enddo

        do j=1,lons+lons-1
          i = i4j(j)
          valf(j) = qc(i) + dq(i)*dt(j)
        enddo

        do i=1,lons
          ss = 0.0
          do j=js(i),js(i+1)-1
            ss = ss + valf(j)*df(j)
          enddo
          qn(i,n) = ss * rdd(i)
        enddo

      enddo
!
      return
      end subroutine fixend_cell_plm_intp
!
! ------------------------------------------------------------------------
      subroutine vertical_cell_advect(lons,londim,levs,nvars,           & 
                                      deltim,ssi,wwi,qql,mass,forward)
!
      use const, only : RTYPE

      implicit none

      integer 	londim,levs,nvars,lons,mass
      real(kind=RTYPE) 	deltim
      real(kind=RTYPE)	ssi(londim,levs+1)
      real(kind=RTYPE)	wwi(londim,levs+1)
      real(kind=RTYPE)	qql(londim,levs,nvars)

      real(kind=RTYPE) ssii(levs+1)
      real(kind=RTYPE) ssid(levs+1),ssia(levs+1)
      real(kind=RTYPE) dd(levs+1),ds(levs),step(10),dd_step
      real(kind=RTYPE) dsfact(levs), sstmp, dpdt, check
      real(kind=RTYPE) rqmm(levs,nvars),rqnn(levs,nvars),rqda(levs,nvars)
      integer km,i,k,n,nst,nstep
      logical forward

      do i=1,lons

        do k=1,levs
          ds(k)=ssi(i,k)-ssi(i,k+1)
        enddo
        dd(1) = 0.0
        do k=2,levs
          dd(k) = wwi(i,k) * deltim
        enddo
        dd(levs+1) = 0.0
!cflx   call def_cfl_max (levs+1,dd,ds,nstep)
        call def_cfl_step(levs+1,dd,ds,step,nstep,levs+1,'advv')

!
       do nst = 1, nstep
!
!       do k=1,levs+1
!         ssii(k)=-ssi(i,k)
!         ssid(k)=-ssi(i,k)+dd(k)
!         ssia(k)=-ssi(i,k)-dd(k)
!       enddo
!hmhj give direction for value larger with k larger
        if ( forward ) then
          do k=1,levs+1
            dd_step= dd(k)*step(nst)
! for ppm interpolation
            ssii(k)=ssi(i,k)
            ssid(k)=ssi(i,k)
            ssia(k)=ssi(i,k)+dd_step
! for plm interpolation
!            ssii(k)=-ssi(i,k)
!            ssid(k)=-ssi(i,k)+dd_step
!            ssia(k)=-ssi(i,k)-dd_step
          enddo
        else
          do k=1,levs+1
            dd_step= dd(k)*step(nst)
! for ppm interpolation
            ssii(k)=ssi(i,k)
            ssid(k)=ssi(i,k)-dd_step
            ssia(k)=ssi(i,k)+dd_step
! for plm interpolation
!            ssii(k)=-ssi(i,k)
!            ssid(k)=-ssi(i,k)+dd_step
!            ssia(k)=-ssi(i,k)-dd_step
          enddo
        endif
        if( mass.eq.1 ) then
          do k=1,levs
            dsfact(k)=(ssid(k)-ssid(k+1))/(ssia(k)-ssia(k+1))
          enddo
        endif
!
        if ( forward ) then
          do n=1,nvars
            do k=1,levs
              rqda(k,n) = qql(i,k,n)
            enddo
          enddo
        else
          do n=1,nvars
            do k=1,levs
              rqmm(k,n) = qql(i,k,n)
            enddo
          enddo
          call vertical_cell_ppm_intp(ssii,rqmm,ssid,rqda,levs,nvars,i)
!          call fixend_cell_plm_intp(ssii,rqmm,ssid,rqda,levs,nvars)
        endif
        if( mass.eq.1 ) then
          do n=1,nvars
            do k=1,levs
              rqda(k,n) = rqda(k,n) * dsfact(k)
            enddo
          enddo
        endif
        call vertical_cell_ppm_intp(ssia,rqda,ssii,rqnn,levs,nvars,i)
!        call fixend_cell_plm_intp(ssia,rqda,ssii,rqnn,levs,nvars)
        do n=1,nvars
          do k=1,levs
            qql(i,k,n)=rqnn(k,n)
          enddo
        enddo
!
       enddo
!
      enddo

      return
      end subroutine vertical_cell_advect
!
!
! -------------------------------------------------------------------------
      subroutine cyclic_cell_ppm_intp(pp,qq,pn,qn,lons,nv,lonp,lonn,sc)
!
! mass conservation in cyclic bc interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in horizontal with mass conservation is under the condition that
! variable value at pp(1)= pp(lons+1)=pn(lons+1)
!
! pp    location at interfac point as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
! lonp  numer of cells for input
! lonn  numer of cells for output
! levs  number of vertical layers
! mono  monotonicity o:no, 1:yes
!
! author : henry.juang@noaa.gov
!
!
      use const, only : RTYPE
      implicit none
!
      real(kind=RTYPE)      pp(lons+1)
      real(kind=RTYPE)      qq(lons  ,nv)
      real(kind=RTYPE)      pn(lons+1)
      real(kind=RTYPE)      qn(lons  ,nv)
      integer   lons,lonp,lonn,nv
      real(kind=RTYPE)      sc
!
      integer   ik,le,kstr,kend
      integer   i,k, kl, kh, kk, kkl, kkh, n
      integer, parameter :: mono=1

      real(kind=RTYPE) locs  (3*lonp)
      real(kind=RTYPE) mass  (3*lonp,nv)
      real(kind=RTYPE) hh    (3*lonp)
      real(kind=RTYPE) fm    (3*lonp)
      real(kind=RTYPE) fn    (3*lonp)
      real(kind=RTYPE) dqmono(3*lonp,nv)
      real(kind=RTYPE) qmi   (3*lonp,nv)
      real(kind=RTYPE) qpi   (3*lonp,nv)
      real(kind=RTYPE) cyclic_length
      real(kind=RTYPE) pnmin,pnmax,locbndmin,locbndmax
      real(kind=RTYPE) dqi,dqimax,dqimin
      real(kind=RTYPE) tl,tl2,tl3,qql,tlp,tlm,tlc
      real(kind=RTYPE) th,th2,th3,qqh,thp,thm,thc
      real(kind=RTYPE) dql(nv),dqh(nv)
      real(kind=RTYPE) dpp,dqq,c1,c2,cc,r3,r6
      real(kind=RTYPE) rdthtl
!
!     cyclic_length = pp(lonp+1) - pp(1)
      cyclic_length = sc
!
! arrange input array cover output location with cyclic boundary
! condition
!
      locs(lonp+1:2*lonp) = pp(1:lonp)
      do i=1,lonp
        locs(i) = locs(i+lonp) - cyclic_length
        locs(i+2*lonp) = locs(i+lonp) + cyclic_length
      enddo
      mass(1       :  lonp,1:nv) = qq(1:lonp,1:nv)
      mass(1+  lonp:2*lonp,1:nv) = qq(1:lonp,1:nv)
      mass(1+2*lonp:3*lonp,1:nv) = qq(1:lonp,1:nv)

      pnmin = pn(1)
      pnmax = pn(lonn+1)
!!    do i=2,lonn
!!      pnmin = min( pnmin, pn(i) )
!!      pnmax = max( pnmax, pn(i) )
!!    enddo

      locbndmin=locs(  lonp+4)
      locbndmax=locs(2*lonp-4)
      if( pnmin.lt.locbndmin-sc ) then
        do i=1,lonn+1
          pn(i)=pn(i)+int((locbndmin-pnmin)/sc)*sc
        enddo
      else if( pnmin.gt.locbndmax ) then
        do i=1,lonn+1
          pn(i)=pn(i)+(int((locbndmax-pnmin)/sc)-1)*sc
        enddo
      endif

      pnmin = pn(1)
      pnmax = pn(lonn+1)
!!    do i=2,lonn
!!      pnmin = min( pnmin, pn(i) )
!!      pnmax = max( pnmax, pn(i) )
!!    enddo

      if( pnmin.lt.locs(lonp+1) ) then
        do i=lonp,1,-1
          if( pnmin.ge.locs(i) .and. pnmin.lt.locs(i+1) ) then
            kstr = i
            go to 10
          endif
        enddo
      else
        do i=lonp+1,2*lonp
          if( pnmin.ge.locs(i) .and. pnmin.lt.locs(i+1) ) then
            kstr = i
            go to 10
          endif
        enddo
      endif
      print *,' Error: can not find kstr: pnmin locs(1) locs(2*lonp) ',&
                                          pnmin,locs(1),locs(2*lonp)
      print *,' Error: pn(1) pn(2) pn(3) ',pn(1),pn(2),pn(3)

 10   kstr=max(3,kstr)

      if( pnmax.lt.locs(2*lonp+1) ) then
        do i=2*lonp,lonp,-1
          if( pnmax.ge.locs(i) .and. pnmax.lt.locs(i+1) ) then
            kend = i+1
            go to 20
          endif
        enddo
      else
        do i=2*lonp+1,3*lonp-1
          if( pnmax.ge.locs(i) .and. pnmax.lt.locs(i+1) ) then
            kend = i+1
            go to 20
          endif
        enddo
      endif
      print *,' Error: cannot get kend: pnmax locs(lonp) locs(3*lonp)',&
                                  kend, pnmax,locs(lonp),locs(3*lonp)
      print *,' Error: pn(lonn-1) pn(lonn) pn(lonn+1) ',               &
                       pn(lonn-1),pn(lonn),pn(lonn+1)

 20   kend=min(3*lonp-2,kend)
!
! prepare grid spacing
!
      do i=kstr-2,kend+2
        hh(i) = locs(i+1)-locs(i)
      enddo
      do i=kstr-1,kend+2
       cc = 1./(hh(i)+hh(i-1))
       fm(i) = hh(i  ) * cc
       fn(i) = hh(i-1) * cc
      enddo
!
! prepare location with monotonic concerns
!
      do n=1,nv
      do i=kstr-2,kend+2
        dqi = 0.25*(mass(i+1,n)-mass(i-1,n))
        dqimax = max(mass(i-1,n),mass(i,n),mass(i+1,n)) - mass(i,n)
        dqimin = mass(i,n) - min(mass(i-1,n),mass(i,n),mass(i+1,n))
        dqmono(i,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
      enddo
      enddo
!
! compute value at interface with monotone
!
      r3 = 1./3.
      do n=1,nv
      do i=kstr-1,kend+2
        qmi(i,n)=mass(i-1,n)*fm(i)+mass(i,n)*fn(i)                 &
             +(dqmono(i-1,n)-dqmono(i,n))*r3
      enddo
      enddo
      qpi(kstr-1:kend+1,1:nv) = qmi(kstr:kend+2,1:nv)
!
! do less diffusive
!
!!      do n=1,nv
!!      do i=kstr-1,kend+2
!!        qmi(i,n)=mass(i,n)-sign(min(abs(2.*dqmono(i,n)),           &
!!                            abs(qmi(i,n)-mass(i,n))),              &
!!                            2.*dqmono(i,n))
!!        qpi(i,n)=mass(i,n)+sign(min(abs(2.*dqmono(i,n)),           &
!!                            abs(qpi(i,n)-mass(i,n))),              &
!!                            2.*dqmono(i,n))
!!      enddo
!!      enddo
!
! do monotonicity within cell
!
      r6 = 1./6.
      if( mono.eq.1 ) then
        do n=1,nv
        do i=kstr-1,kend+1
          c1=qpi(i,n)-mass(i,n)
          c2=mass(i,n)-qmi(i,n)
          if( c1*c2.le.0.0 ) then
            qmi(i,n)=mass(i,n)
            qpi(i,n)=mass(i,n)
          else
            cc=qpi(i,n)-qmi(i,n)
            c1=cc*(mass(i,n)-0.5*(qpi(i,n)+qmi(i,n)))
            c2=cc*cc*r6
            if( c1.gt.c2 ) then
              qmi(i,n)=3.*mass(i,n)-2.*qpi(i,n)
            else if( c1.lt.-c2 ) then
              qpi(i,n)=3.*mass(i,n)-2.*qmi(i,n)
            endif
          endif
        enddo
        enddo
      endif
!
! start interpolation by integral of ppm
!
      kkl = kstr
      tl=(pn(1)-locs(kkl))/hh(kkl)
      tl2=tl*tl
      tl3=tl2*tl
      tlp = tl3-tl2
      tlm = tl3-2.*tl2+tl
      tlc = -2.*tl3+3.*tl2
      do n=1,nv
        dql(n)=tlp*qpi(kkl,n)+tlm*qmi(kkl,n)+tlc*mass(kkl,n)
      enddo

      do i=1,lonn

        kl = i
        kh = i + 1
! find kkh
        do kk=kkl+1,kend+2
          if( pn(kh).lt.locs(kk) ) then
            kkh = kk-1
            go to 100
          endif
        enddo

        print *,' Error in cyclic_cell_ppm_intp location not found '
        print *,' lons=',lons,' lonp=',lonp,' lonn=',lonn
        print *,' pnmin=',pnmin,' pnmax=',pnmax
        print *,' pn(1)=',pn(1),' pn(lonn+1)=',pn(lonn+1)
        print *,' kstr =',kstr ,' kend =',kend
        print *,' kh=',kh,' pn(kh)=',pn(kh)
        print *,' kkl +1=',kkl +1,' locs(kkl +1)=',locs(kkl +1)
        print *,' kend+1=',kend+1,' locs(kend+1)=',locs(kend+1)
        call abort

 100    continue
! mass interpolate
        th=(pn(kh)-locs(kkh))/hh(kkh)
        th2=th*th
        th3=th2*th
        thp = th3-th2
        thm = th3-2.*th2+th
        thc = -2.*th3+3.*th2
        do n=1,nv
          dqh(n)=thp*qpi(kkh,n)+thm*qmi(kkh,n)+thc*mass(kkh,n)
        enddo
        if( kkh.eq.kkl ) then
          rdthtl = th - tl
          if ( rdthtl.ne.0.0 ) rdthtl = 1. / rdthtl
          do n=1,nv
!hmhj            qn(i,n) = (dqh(n)-dql(n))/(th-tl)
            qn(i,n) = (dqh(n)-dql(n))*rdthtl
          enddo
        else if( kkh.gt.kkl ) then
          dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
          do kk=kkl+1,kkh-1
            dpp = dpp + hh(kk)
          enddo
          do n=1,nv
            dql(n) = mass(kkl,n)-dql(n)
            dqq  = dql(n)*hh(kkl) + dqh(n)*hh(kkh)
            do kk=kkl+1,kkh-1
              dqq = dqq + mass(kk,n)*hh(kk)
            enddo
            qn(i,n) = dqq / dpp
          enddo
        else
          print *,' Error in cyclic_cell_ppm_intp location messed up '
          print *,' kkl=',kkl,' kkh=',kkh
          print *,' kh=',kh,' pn(kh)=',pn(kh)
          print *,' kkl-1=',kkl-1,' locs(kkl-1)=',locs(kkl-1)
          call abort
        endif

! next one
        kkl = kkh
        tl = th
        do n=1,nv
          dql(n) = dqh(n)
        enddo

      enddo
!
      return
      end subroutine cyclic_cell_ppm_intp
!
! ------------------------------------------------------------------------
      subroutine vertical_cell_ppm_intp(pp,qq,pn,qn,levs,nvars,i)
!
! mass conservation in vertical interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in vertical with mass conservation is under the condition that
! pp(1)=pn(1), pp(levs+1)=pn(levs+1)
!
! pp    pressure at interfac level as input
! qq    quantity at layer as input
! pn    pressure at interface of new grid structure as input
! qn    quantity at layer as output
! levs  numer of verical layers
!
! author : henry.juang@noaa.gov
!
      use const, only : RTYPE

      implicit none
!
      real(kind=RTYPE)      pp(levs+1)
      real(kind=RTYPE)      qq(levs,nvars)
      real(kind=RTYPE)      pn(levs+1)
      real(kind=RTYPE)      qn(levs,nvars)
      integer   levs,nvars
!
      real(kind=RTYPE)      massm,massc,massp,massbot,masstop
      real(kind=RTYPE)      qmi(levs,nvars),qpi(levs,nvars)
      real(kind=RTYPE)      dql(nvars),dqh(nvars)
      real(kind=RTYPE)      hh(levs)
      real(kind=RTYPE)      dqi,dqimax,dqimin,dqmono(levs,nvars)
      real(kind=RTYPE)      tl,tl2,tl3,tlp,tlm,tlc
      real(kind=RTYPE)      th,th2,th3,thp,thm,thc
      real(kind=RTYPE)      dpp,dqq,c1,c2,rdthtl
      integer   i,k, kl, kh, kk, kkl, kkh,n
      integer, parameter :: mono=1
!
      if( pp(1).ne.pn(1) .or. pp(levs+1).ne.pn(levs+1) ) then
        print *,' Error in vertical_cell_ppm_intp for domain values '
        print *,' i pp1 pn1 ppt pnt ',i,                            &
                pp(1),pn(1),pp(levs+1),pn(levs+1)
        call abort
      endif
!
! prepare thickness for grid
!
      do k=1,levs
        hh(k) = pp(k+1)-pp(k)
      enddo
!
! prepare location with monotonic concerns
!
      do n=1,nvars
        massbot=(3.*hh(1)+hh(2))*qq(1,n)-2.*hh(1)*qq(2,n)
        massm = massbot/(hh(1)+hh(2))
        massc = qq(1  ,n)
        massp = qq(1+1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(1,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        do k=2,levs-1
          massp = qq(k+1,n)
          massc = qq(k  ,n)
          massm = qq(k-1,n)
          dqi = 0.25*(massp-massm)
          dqimax = max(massm,massc,massp) - massc
          dqimin = massc - min(massm,massc,massp)
          dqmono(k,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        enddo
        masstop=(3.*hh(levs)+hh(levs-1))*qq(levs,n)                  &
                   -2.*hh(levs)*qq(levs-1,n)
        massp = masstop/(hh(levs)+hh(levs-1))
        massc = qq(levs  ,n)
        massm = qq(levs-1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(levs,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
!
! compute value at interface with momotone
!
        do k=2,levs
          qmi(k,n)=(qq(k-1,n)*hh(k)+qq(k,n)*hh(k-1))/(hh(k)+hh(k-1)) &
             +(dqmono(k-1,n)-dqmono(k,n))/3.0
        enddo
        do k=1,levs-1
          qpi(k,n)=qmi(k+1,n)
        enddo
        qmi(1,n)=qq(1,n)
        qpi(1,n)=qq(1,n)
        qmi(levs,n)=qq(levs,n)
        qpi(levs,n)=qq(levs,n)
      enddo
!
! do monotonicity
!
      if( mono.eq.1 ) then
        do n=1,nvars
        do k=1,levs
          c1=qpi(k,n)-qq(k,n)
          c2=qq(k,n)-qmi(k,n)
          if( c1*c2.le.0.0 ) then
            qmi(k,n)=qq(k,n)
            qpi(k,n)=qq(k,n)
          endif
        enddo
        do k=1,levs
          c1=(qpi(k,n)-qmi(k,n))*(qq(k,n)-0.5*(qpi(k,n)+qmi(k,n)))
          c2=(qpi(k,n)-qmi(k,n))*(qpi(k,n)-qmi(k,n))/6.
          if( c1.gt.c2 ) then
            qmi(k,n)=3.*qq(k,n)-2.*qpi(k,n)
          else if( c1.lt.-c2 ) then
            qpi(k,n)=3.*qq(k,n)-2.*qmi(k,n)
          endif
        enddo
        enddo
      endif
!
! start interpolation by integral of ppm spline
!
      kkl = 1
      tl=0
      do n=1,nvars
        dql(n)=0.0
      enddo

      do k=1,levs

        kl = k
        kh = k + 1
! find kkh
        do kk=kkl+1,levs+1
          if( pn(kh).ge.pp(kk) ) then
            kkh = kk-1
            go to 100
          endif
        enddo
        print *,' Error in vertical_cell_ppm_intp for no lev found '
        print *,' i kh kl ',i,kh,kl
        print *,' pn ',(pn(kk),kk=1,levs+1)
        print *,' pp ',(pp(kk),kk=1,levs+1)
        call abort
 100    continue
        th=(pn(kh)-pp(kkh))/hh(kkh)
        th2=th*th
        th3=th2*th
        thp = th3-th2
        thm = th3-2.*th2+th
        thc = -2.*th3+3.*th2
        do n=1,nvars
          dqh(n)=thp*qpi(kkh,n)+thm*qmi(kkh,n)+thc*qq(kkh,n)
        enddo
! mass interpolate
        if( kkh.eq.kkl ) then
          rdthtl = th - tl
          if ( rdthtl.ne.0.0 ) rdthtl = 1. / rdthtl
          do n=1,nvars
!hmhj            qn(k,n) = (dqh(n)-dql(n))/(th-tl)
            qn(k,n) = (dqh(n)-dql(n))*rdthtl
          enddo
        else if( kkh.gt.kkl ) then
          dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
          do kk=kkl+1,kkh-1
            dpp = dpp + hh(kk)
          enddo
          do n=1,nvars
            dql(n) = qq(kkl,n)-dql(n)
            dqq  = dql(n)*hh(kkl) + dqh(n)*hh(kkh)
            do kk=kkl+1,kkh-1
              dqq = dqq + qq(kk,n)*hh(kk)
            enddo
            qn(k,n) = dqq / dpp
          enddo
        else
          print *,' Error in vertical_cell_ppm_intp for lev messed up '
          print *,' i kh kl ',i,kh,kl
          print *,' pn ',(pn(kk),kk=1,levs+1)
          print *,' pp ',(pp(kk),kk=1,levs+1)
          call abort
        endif
! next one
        kkl = kkh
        tl  = th
        do n=1,nvars
          dql(n) = dqh(n)
        enddo

      enddo     ! end of k loop
!
      return
      end subroutine vertical_cell_ppm_intp
!!
!!
      subroutine def_cfl_step (im,dist,del,step,nstep,k,job)
!
! compute the deformation cfl condition
! select the maxima value of the deformation CFL and provide step to
! avoid it.
!
!
      use const, only : RTYPE

      implicit none
!
      integer im,nstep
      real(kind=RTYPE) dist(im),del(im-1),step(10)
! local
      integer   n,k,nchk
      real(kind=RTYPE)      rstep,check,check_max,check_point
      real(kind=RTYPE)      safe_step,last_step
      character*4 job
!
      check_point=1.00
      safe_step=0.99
      nstep = 1
      step(1) = 1.0

      check_max = 0.0
      check_loop: do n=1,im-1
        check = abs ( (dist(n+1)-dist(n))/del(n) )
        if( check.lt.check_point ) then
          cycle check_loop
        else
          if ( check .gt. check_max ) then
             nchk=n
             check_max=check
          endif
!!          check_max = max( check_max, check )
        endif
      enddo check_loop
      if(check_max.ge.check_point) then
        nstep = int(check_max/safe_step) + 1
        if ( job .eq. 'advv' ) then
          print *,' max def_cfl ',check_max,' needs ',nstep,    &
          'steps at level',im+1-nchk,'of',im,'in ',job,' processing'
        else if ( job .eq. 'advx' .and. im .gt. 25 ) then
          print *,' max def_cfl ',check_max,' needs ',nstep,    &
          'steps in',nchk,'of',im,'at level',k,'in ',job,' processing'
        else if ( job .eq. 'advy' ) then
          print *,' max def_cfl ',check_max,' needs ',nstep,    &
          'steps in',nchk,'of',im,'at level',k,'in ',job,' processing'
        endif
        rstep =  safe_step / check_max
        do n=1,nstep-1
          step(n) = rstep
        enddo
        last_step = 1. - ( nstep - 1 ) * rstep
        step(nstep) = last_step
      endif
!
      return
      end subroutine def_cfl_step
!!
      subroutine mymaxmin(a,im,ix,kx,ch)
      use const, only : RTYPE
        implicit none
        real(kind=RTYPE) :: fmax
        real(kind=RTYPE) :: fmin
        integer :: i
        integer :: im
        integer :: ix
        integer :: k
        integer :: kx

      real(kind=RTYPE) a(ix,kx)
      character*(*) ch
      do k=1,kx
        fmin=a(1,k)
        fmax=a(1,k)
        do i=1,im
          fmin=min(fmin,a(i,k))
          fmax=max(fmax,a(i,k))
        enddo
        print *,' max=',fmax,' min=',fmin,' at k=',k,' for ',ch
      enddo
      return
      end subroutine mymaxmin
!
! ------------------------------------------------------------------------
!
      subroutine ndslfv_update (lonsperlat,vdzonl,vdmerd,vdzonlr,vdmerdr,deltim,forward)

!  update all horizontal components into momentum eqs
!  for Semi-Lagrangian vertical advection

      use index
      use param, only : nx,my,lev,my_max
      use const, only : onocos,radsq,RTYPE
      use grid , only : dlphi,dtphi

      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)    vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max)
      real(kind=RTYPE)    vdmerdr(nxp,lev,my_max),vdzonlr(nxp,lev,my_max)
      integer i,ii,k,lan,lat,lons_lat
      integer dt2
      logical forward

      if ( forward ) then
        dt2 = deltim
      else
        dt2 = 2. * deltim
      endif
!
!$omp parallel do                                                   &
!$omp private(lan,lat,lons_lat,i,k)                                 &
!$omp schedule(dynamic)
      do lan=1,jlistnum
!
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
        do k=1,lev
         do i=1,lons_lat
!ttl           vdzonl(i,k,lan) = (vdzonlr(i,k,lan)-dlphi(i,k,lan)/radsq) &
!ttl                             * dt2 + vdzonl(i,k,lan)
!ttl           vdmerd(i,k,lan) = (vdmerdr(i,k,lan)-dtphi(i,k,lan)/radsq  &
!ttl                             / onocos(lat))*dt2 + vdmerd(i,k,lan)
           vdzonl(i,k,lan) = vdzonlr(i,k,lan) * dt2 + vdzonl(i,k,lan)
           vdmerd(i,k,lan) = vdmerdr(i,k,lan) * dt2 + vdmerd(i,k,lan)
         enddo
        enddo
      enddo
!$omp end parallel do
!
! ===============================
!
      return
      end subroutine ndslfv_update
!
      subroutine ndslfv_update_3tl (lonsperlat,vdzonl,vdmerd,vdzonlr,vdmerdr,deltim)

!  update all horizontal components into momentum eqs
!  for Semi-Lagrangian vertical advection for 3tl

      use index
      use param, only : nx,my,lev,my_max
      use const, only : onocos,radsq,RTYPE
      use grid , only : dlphi,dtphi

      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

!ch   real    vdmerd(nx+3,lev,my_max),vdzonl(nx+3,lev,my_max)
      real(kind=RTYPE) vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max)
      real(kind=RTYPE) vdmerdr(nxp,lev,my_max),vdzonlr(nxp,lev,my_max)
      integer i,ii,k,lan,lat,lons_lat
      integer dt2

      dt2 = 2. * deltim
!
!$omp parallel do                                                   &
!$omp private(lan,lat,lons_lat,i,k)                                 &
!$omp schedule(dynamic)
      do lan=1,jlistnum
!
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
        do k=1,lev
         do i=1,lons_lat
           vdzonl(i,k,lan) = (vdzonlr(i,k,lan)-dlphi(i,k,lan)/radsq) &
                             * dt2 + vdzonl(i,k,lan)
           vdmerd(i,k,lan) = (vdmerdr(i,k,lan)-dtphi(i,k,lan)/radsq  &
                             / onocos(lat))*dt2 + vdmerd(i,k,lan)
         enddo
        enddo
      enddo
!$omp end parallel do
!
! ===============================
!
      return
      end subroutine ndslfv_update_3tl

!-------------------------------------------------------------------------
!CWB2021 note, the double precision is called by reducepick/reduceintp
!-------------------------------------------------------------------------
      subroutine cyclic_cell_ppm_intp_dp(pp,qq,pn,qn,lons,nv,lonp,lonn,sc)
!
! mass conservation in cyclic bc interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in horizontal with mass conservation is under the condition that
! variable value at pp(1)= pp(lons+1)=pn(lons+1)
!
! pp    location at interfac point as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
! lonp  numer of cells for input
! lonn  numer of cells for output
! levs  number of vertical layers
! mono  monotonicity o:no, 1:yes
!
! author : henry.juang@noaa.gov
!
!
      implicit none
!
      real      pp(lons+1)
      real      qq(lons  ,nv)
      real      pn(lons+1)
      real      qn(lons  ,nv)
      integer   lons,lonp,lonn,nv
      real      sc
!
      real      locs  (3*lonp)
      real      mass  (3*lonp,nv)
      real      hh    (3*lonp)
      real      fm    (3*lonp)
      real      fn    (3*lonp)
      real      dqmono(3*lonp,nv)
      real      qmi   (3*lonp,nv)
      real      qpi   (3*lonp,nv)
      real      cyclic_length
      integer   ik,le,kstr,kend
      real      pnmin,pnmax,locbndmin,locbndmax
      real      dqi,dqimax,dqimin
      real      tl,tl2,tl3,qql,tlp,tlm,tlc
      real      th,th2,th3,qqh,thp,thm,thc
      real      dql(nv),dqh(nv)
      real      dpp,dqq,c1,c2,cc,r3,r6
      real      rdthtl
      integer   i,k, kl, kh, kk, kkl, kkh, n
      integer, parameter :: mono=1
!
!     cyclic_length = pp(lonp+1) - pp(1)
      cyclic_length = sc
!
! arrange input array cover output location with cyclic boundary
! condition
!
      locs(lonp+1:2*lonp) = pp(1:lonp)
      do i=1,lonp
        locs(i) = locs(i+lonp) - cyclic_length
        locs(i+2*lonp) = locs(i+lonp) + cyclic_length
      enddo
      mass(1       :  lonp,1:nv) = qq(1:lonp,1:nv)
      mass(1+  lonp:2*lonp,1:nv) = qq(1:lonp,1:nv)
      mass(1+2*lonp:3*lonp,1:nv) = qq(1:lonp,1:nv)

      pnmin = pn(1)
      pnmax = pn(lonn+1)
!!    do i=2,lonn
!!      pnmin = min( pnmin, pn(i) )
!!      pnmax = max( pnmax, pn(i) )
!!    enddo

      locbndmin=locs(  lonp+4)
      locbndmax=locs(2*lonp-4)
      if( pnmin.lt.locbndmin-sc ) then
        do i=1,lonn+1
          pn(i)=pn(i)+int((locbndmin-pnmin)/sc)*sc
        enddo
      else if( pnmin.gt.locbndmax ) then
        do i=1,lonn+1
          pn(i)=pn(i)+(int((locbndmax-pnmin)/sc)-1)*sc
        enddo
      endif

      pnmin = pn(1)
      pnmax = pn(lonn+1)
!!    do i=2,lonn
!!      pnmin = min( pnmin, pn(i) )
!!      pnmax = max( pnmax, pn(i) )
!!    enddo

      if( pnmin.lt.locs(lonp+1) ) then
        do i=lonp,1,-1
          if( pnmin.ge.locs(i) .and. pnmin.lt.locs(i+1) ) then
            kstr = i
            go to 10
          endif
        enddo
      else
        do i=lonp+1,2*lonp
          if( pnmin.ge.locs(i) .and. pnmin.lt.locs(i+1) ) then
            kstr = i
            go to 10
          endif
        enddo
      endif
      print *,' Error: can not find kstr: pnmin locs(1) locs(2*lonp) ',&
                                          pnmin,locs(1),locs(2*lonp)
      print *,' Error: pn(1) pn(2) pn(3) ',pn(1),pn(2),pn(3)

 10   kstr=max(3,kstr)

      if( pnmax.lt.locs(2*lonp+1) ) then
        do i=2*lonp,lonp,-1
          if( pnmax.ge.locs(i) .and. pnmax.lt.locs(i+1) ) then
            kend = i+1
            go to 20
          endif
        enddo
      else
        do i=2*lonp+1,3*lonp-1
          if( pnmax.ge.locs(i) .and. pnmax.lt.locs(i+1) ) then
            kend = i+1
            go to 20
          endif
        enddo
      endif
      print *,' Error: cannot get kend: pnmax locs(lonp) locs(3*lonp)',&
                                  kend, pnmax,locs(lonp),locs(3*lonp)
      print *,' Error: pn(lonn-1) pn(lonn) pn(lonn+1) ',               &
                       pn(lonn-1),pn(lonn),pn(lonn+1)

 20   kend=min(3*lonp-2,kend)
!
! prepare grid spacing
!
      do i=kstr-2,kend+2
        hh(i) = locs(i+1)-locs(i)
      enddo
      do i=kstr-1,kend+2
       cc = 1./(hh(i)+hh(i-1))
       fm(i) = hh(i  ) * cc
       fn(i) = hh(i-1) * cc
      enddo
!
! prepare location with monotonic concerns
!
      do n=1,nv
      do i=kstr-2,kend+2
        dqi = 0.25*(mass(i+1,n)-mass(i-1,n))
        dqimax = max(mass(i-1,n),mass(i,n),mass(i+1,n)) - mass(i,n)
        dqimin = mass(i,n) - min(mass(i-1,n),mass(i,n),mass(i+1,n))
        dqmono(i,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
      enddo
      enddo
!
! compute value at interface with monotone
!
      r3 = 1./3.
      do n=1,nv
      do i=kstr-1,kend+2
        qmi(i,n)=mass(i-1,n)*fm(i)+mass(i,n)*fn(i)                 &
             +(dqmono(i-1,n)-dqmono(i,n))*r3
      enddo
      enddo
      qpi(kstr-1:kend+1,1:nv) = qmi(kstr:kend+2,1:nv)
!
! do less diffusive
!
!!      do n=1,nv
!!      do i=kstr-1,kend+2
!!        qmi(i,n)=mass(i,n)-sign(min(abs(2.*dqmono(i,n)),           &
!!                            abs(qmi(i,n)-mass(i,n))),              &
!!                            2.*dqmono(i,n))
!!        qpi(i,n)=mass(i,n)+sign(min(abs(2.*dqmono(i,n)),           &
!!                            abs(qpi(i,n)-mass(i,n))),              &
!!                            2.*dqmono(i,n))
!!      enddo
!!      enddo
!
! do monotonicity within cell
!
      r6 = 1./6.
      if( mono.eq.1 ) then
        do n=1,nv
        do i=kstr-1,kend+1
          c1=qpi(i,n)-mass(i,n)
          c2=mass(i,n)-qmi(i,n)
          if( c1*c2.le.0.0 ) then
            qmi(i,n)=mass(i,n)
            qpi(i,n)=mass(i,n)
          else
            cc=qpi(i,n)-qmi(i,n)
            c1=cc*(mass(i,n)-0.5*(qpi(i,n)+qmi(i,n)))
            c2=cc*cc*r6
            if( c1.gt.c2 ) then
              qmi(i,n)=3.*mass(i,n)-2.*qpi(i,n)
            else if( c1.lt.-c2 ) then
              qpi(i,n)=3.*mass(i,n)-2.*qmi(i,n)
            endif
          endif
        enddo
        enddo
      endif
!
! start interpolation by integral of ppm
!
      kkl = kstr
      tl=(pn(1)-locs(kkl))/hh(kkl)
      tl2=tl*tl
      tl3=tl2*tl
      tlp = tl3-tl2
      tlm = tl3-2.*tl2+tl
      tlc = -2.*tl3+3.*tl2
      do n=1,nv
        dql(n)=tlp*qpi(kkl,n)+tlm*qmi(kkl,n)+tlc*mass(kkl,n)
      enddo

      do i=1,lonn

        kl = i
        kh = i + 1
! find kkh
        do kk=kkl+1,kend+2
          if( pn(kh).lt.locs(kk) ) then
            kkh = kk-1
            go to 100
          endif
        enddo

        print *,' Error in cyclic_cell_ppm_intp location not found '
        print *,' lons=',lons,' lonp=',lonp,' lonn=',lonn
        print *,' pnmin=',pnmin,' pnmax=',pnmax
        print *,' pn(1)=',pn(1),' pn(lonn+1)=',pn(lonn+1)
        print *,' kstr =',kstr ,' kend =',kend
        print *,' kh=',kh,' pn(kh)=',pn(kh)
        print *,' kkl +1=',kkl +1,' locs(kkl +1)=',locs(kkl +1)
        print *,' kend+1=',kend+1,' locs(kend+1)=',locs(kend+1)
        call abort

 100    continue
! mass interpolate
        th=(pn(kh)-locs(kkh))/hh(kkh)
        th2=th*th
        th3=th2*th
        thp = th3-th2
        thm = th3-2.*th2+th
        thc = -2.*th3+3.*th2
        do n=1,nv
          dqh(n)=thp*qpi(kkh,n)+thm*qmi(kkh,n)+thc*mass(kkh,n)
        enddo
        if( kkh.eq.kkl ) then
          rdthtl = th - tl
          if ( rdthtl.ne.0.0 ) rdthtl = 1. / rdthtl
          do n=1,nv
!            qn(i,n) = (dqh(n)-dql(n))/(th-tl)
            qn(i,n) = (dqh(n)-dql(n))*rdthtl
          enddo
        else if( kkh.gt.kkl ) then
          dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
          do kk=kkl+1,kkh-1
            dpp = dpp + hh(kk)
          enddo
          do n=1,nv
            dql(n) = mass(kkl,n)-dql(n)
            dqq  = dql(n)*hh(kkl) + dqh(n)*hh(kkh)
            do kk=kkl+1,kkh-1
              dqq = dqq + mass(kk,n)*hh(kk)
            enddo
            qn(i,n) = dqq / dpp
          enddo
        else
          print *,' Error in cyclic_cell_ppm_intp location messed up '
          print *,' kkl=',kkl,' kkh=',kkh
          print *,' kh=',kh,' pn(kh)=',pn(kh)
          print *,' kkl-1=',kkl-1,' locs(kkl-1)=',locs(kkl-1)
          call abort
        endif

! next one
        kkl = kkh
        tl = th
        do n=1,nv
          dql(n) = dqh(n)
        enddo

      enddo
!
      return
      end subroutine cyclic_cell_ppm_intp_dp
