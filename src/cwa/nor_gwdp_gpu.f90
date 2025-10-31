      subroutine nor_gwdp_gpu(dump, nxjp, nx, lev, u, v, t, q, plt, pk, pk2, hi, dt, grav, rgas, sinl &
                          , cosl, avgdrag_u, avgdrag_v, cp, ptop)
!
! reference : Scinocca(2002,2003) JAS and Andrew Orr(2010) J. of Climate
!
!---- input:
!---- nxj    : number of point in lontitude of reduced gauslian grid.
!---- nx     : dimension of lontitude.
!---- lev    : dimension of vertical.
!---- dt     : 2.*dt dt: time step of model integration
!---- grav   : gravity of earth
!     rgas   : dry gas constant (=287.)
!---- sinl   : sin of latitude
!     t     : temperature  (k)
!     q     : specific humidity (kg/kg)
!     plt   : pressure at odd-level (mb)
!     pk    : p**kcapa at odd-level ((=p/1000.)**0.286)
!     pk2   : p**kcapa at even-level ((=p/1000.)**0.286)
!     hi    : geopotential height (m2/s2)
!---- u      : velocity of x-component
!---- v      : velocity of y-component
!---- output:
!---  u,v
!-----------------------------------------
!---- den    : density of atmosphere
!---- sn2    : Brunt-Vasisla frequency
!---- p2     : pressure at even-level (mb)
!---- nc        : number of spectral elements of the gravity wave phase speed.
!---- nphi      : number of azimuths of the gravity waves.
!---- cmin      : minimun phase speed of the gravity waves(m/s).
!---- cmax      : the maximun phase speed of the gravity waves(m/s).
!---- taul      : constant=0.6 s/m
!---- cstar     : constant=1.
!---- pw        : constant=1.5 exponent power of gravity wave frequency
!---- s         : constant=1. low-m exponent
!---- lunch_lev     : lunch level of gravity wave
!---- wavelenth : characteric vertical wave-length
!---- fluxtoatl : total lunch E-P flux in each nc direction
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!

!ch>
         use index, only: nxdef, jlist1, jlistnum
         use const, only: RTYPE
         use param, only: my, my_max, ncld

!ch<

         implicit none

         integer j, nx, lev, jj, nxjp(my), myim(my_max), dump
         real dt, grav, rgas, sinl(my), cosl(my), cp, sinljj(my_max), cosljj(my_max)

         real plt(nx, lev, my_max), &
            avgdrag_u(lev, my_max), avgdrag_v(lev, my_max)
         real(kind=RTYPE) u(nx, lev, my_max), v(nx, lev, my_max), t(nx, lev, my_max), &
            q(nx, lev*ncld, my_max), hi(nx, lev, my_max), pk(nx, lev, my_max), &
            pk2(nx, lev, my_max)

! local array
         integer, parameter :: nc = 20, nphi = 4
         integer, parameter :: ltop = 1, lunch_lev = 44
         real p2, den, sn2(nx, ltop:lunch_lev, my_max), &
            pt
         real fluxsat(nc, nx, ltop:lunch_lev, nphi, my_max), &
            flux(nc, nx, ltop:lunch_lev, nphi, my_max), &
            dflux(nx, ltop:lunch_lev, nphi, my_max), &
            fluxsat_nc(nx, ltop:lunch_lev, nphi, my_max), &
            flux_nc(nx, ltop:lunch_lev, nphi, my_max), &
            dflux_u(nx, ltop:lunch_lev - 1, my_max), &
            dflux_v(nx, ltop:lunch_lev - 1, my_max), &
            a(nx, my_max), cn, uv_lunch, uv, uv1, &
            fxin(nc), fxout(nc), drag_u(nx, ltop:lunch_lev, my_max), &
            drag_v(nx, ltop:lunch_lev, my_max)
         real mstar
!
         integer ilon, jlat, kk, i, k, nc1, n, nn
         real cmin, cmax, taul, cstar, pw, s, wavelenth, fluxtotal, ptop
         real pi, omega, acor(my_max), cor(my_max), wtcosl(my_max), dz1, s1, s2, &
              p, corp(my_max), xmin, xmax, ndir
         real b1, b2, dc, xbar, x, s3, uhat, fxtotal, chat, flux1, flux2, dflux1
         real uhat_z0, uhat_z1, dpk, wsold, wsnew, chata(nc, nphi)

!      data cmin/0.25/, cmax/2000./ taul/0.6/
!      data cstar/1./, pw/1./, s/1./, wavelenth/2000./
!      data fluxtotal/3.75e-4/, lunch_lev/44/, ltop/1/
         data cmin/0.25/, cmax/1000./ taul/0.25/
         data cstar/1./, wavelenth/2000./
         data fluxtotal/3.75e-4/, pw/1.5/, s/0./    !CG1
!
         !subroutine normalize inlined
         real     newsum, ave, var, ss, ep, chattmp
         integer :: async_id = -1
         !GPU register
         real :: flux_ncr, fluxr, fluxsat_ncr, sn2r, denr, ar, dfluxr, fluxsatr, &
            fluxtmp, p21, pk2r, pk2r1, pt1, fluxs(nc)
         
         write(*,*) 3340
         ilon = 1142
         jlat = nx/4
         kk = 31
         pi = 4.*atan(1.)
         omega = 7.292e-5
         mstar = 2.*pi/wavelenth
         p = 2.-pw
         !$acc data create(myim, sinljj, cosljj, cor, acor, wtcosl, corp, &
         !$acc&     sn2, a, dflux, fluxsat_nc, flux, fluxsat, &
         !$acc&     flux_nc, dflux_u, dflux_v, drag_u, &
         !$acc&     drag_v, chata) async(async_id)
         !$acc parallel loop async(async_id) private(j)
         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
            sinljj(jj) = sinl(j)
            cosljj(jj) = cosl(j)
            
   !
            cor(jj) = 2.*omega*sinljj(jj)
            acor(jj) = abs(cor(jj))
   !
   !ch   wtcosl=nxj*cosl
            wtcosl(jj) = nxdef(j)*cosljj(jj)
            corp(jj) = acor(jj)**p
         end do
   !
   ! ... pressure at even-level (p2)
   ! ... potential temperature (pt)
   ! ... density (den)
   !
         
   !
   ! ... brunt-vaisala frequency (sn2)
   !
         !$acc parallel loop async(async_id) collapse(2) gang private(dz1, s1, s2)
         do jj = 1, jlistnum
            do k = 1, lunch_lev
               !$acc loop vector
               do i = 1, myim(jj)
                  pt = t(i, k, jj)*(1.0 + 0.608*q(i, k, jj))/pk(i, k, jj)
                  pt1 = t(i, k + 1, jj)*(1.0 + 0.608*q(i, k + 1, jj))/pk(i, k + 1, jj)
                  dz1 = (hi(i, k, jj) - hi(i, k + 1, jj))/grav
                  s1 = pt - pt1
                  s2 = (pt + pt1)/2.
                  sn2(i, k, jj) = grav*s1/(s2*dz1)
                  if (sn2(i, k, jj) .le. acor(jj)) sn2(i, k, jj) = acor(jj)
               end do
            end do
         end do
   !      do i=1,nxj
   !        sn2(i,lev)=sn2(i,lev-1)
   !      enddo
   !      if(j.eq.jlat)then
   !          call findmaxf(nxj,nx,lev,sn2,vmax,iv,kv,vmin,ivm,kvm)
   !          print *,' sn2=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
   !          call findmaxf(nxj,nx,lev,den,vmax,iv,kv,vmin,ivm,kvm)
   !          print *,' den=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
   !      endif
   !
   !  coefficient of A(indepedent of height) in eq-(4) (Orr,2010)
   !
   !      if(j.eq.jlat)print *,'mstar,p,cor,corp=',mstar,p,cor,corp
   !
    
         !$acc parallel loop async(async_id) gang vector collapse(2)
         do jj = 1, jlistnum
            do i = 1, nx
               a(i, jj) = 0.
               if (i .le. myim(jj)) then
                  ar = cstar*mstar**3*(sn2(i, lunch_lev, jj)**p - corp(jj))/p
                  if (ar .gt. 0.) a(i, jj) = ar
               end if
            end do
         end do
   !      if(j.eq.jlat)then
   !          call findmaxf(nxj,nx,1,a,vmax,iv,kv,vmin,ivm,kvm)
   !          print *,' a=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
   !      endif
   !
         k = lunch_lev
         
   !      if(j.eq.jlat)then
   !        do nn=1,nphi
   !          print *,' lunch nn,u,v=',nn,uv_lunch(ilon,nn)
   !        enddo
   !      endif
   !
   ! do coordinate stretch, eq-(28),(29),(30) (Scinocca,2003)
   ! cn: transformed wave phase speed
   !
         xmin = 1./cmin
         xmax = 1./cmax
         b1 = (xmax - xmin)/(exp((xmax - xmin)/taul) - 1)
         b2 = xmin - b1
!
         dc = xmax - xmin
         nc1 = nc - 1
         !$acc parallel loop async(async_id) private(xbar, x)
         do n = 1, nc
            xbar = xmin + dc*(n - 1)/nc1
            x = b1*exp((xbar - xmin)/taul) + b2
            cn = 1./x
            chata(n, 1) = cn
            chata(n, 2) = -cn
            chata(n, 3) = cn
            chata(n, 4) = -cn
         end do
   !
   !      if(j.eq.jlat)then
   !          print *,' b1,b2=',b1,b2
   !          print *,' cn=',(cn(n),n=1,nc)
   !      endif
   !
         s3 = s + 3
   !
   !  for lunch level, assign fluxtotal value to spectral flux density and
   !  saturate flux density and assign u,v at lunch  level to variable 'uv'
   !  for pw=1,p=2-pw=1 and uhat=0, the formula of flux density is simply rewriten.
   !
   !
   !CWB2015 turn off for reproducible
   !!$omp  parallel do default(none) &
   !!$omp  private(nn,ndir,k,i,uv,uhat_z1,uhat_z0,n,chat,flux1,flux2,   &
   !!$omp  fxout,dflux1,fxtotal) &
   !!$omp  shared(ltop,lunch_lev,nxj,p,fxin, &
   !!$omp  u,v,uv_lunch,cn,cstar,den,a,sn2,s3,mstar,fluxsat,flux,dflux, &
   !!$omp  fluxsat_nc,flux_nc)
         uhat = 0.
         !$acc parallel loop async(async_id) collapse(3) gang vector private(ndir, &
         !$acc&         chat, flux1, flux2, newsum, ave, var, ep, ss, uv, uv1, &
         !$acc&         fxin, k, uv_lunch, flux_ncr, den, &
         !$acc&         fluxr, fluxsat_ncr, sn2r, denr, ar)
         do jj = 1, jlistnum
            do nn = 1, nphi
               do i = 1, nx
                  if (i .le. myim(jj)) then
                     ndir = 1.
                     k = lunch_lev
                     uv_lunch = u(i, k, jj)
                     if (nn .eq. 3 .or. nn .eq. 4) uv_lunch = v(i, k, jj)
                     if (nn .eq. 2 .or. nn .eq. 4) ndir = -1.
                     newsum = fluxtotal*ndir
                     fluxsat_ncr = 0.
                     sn2r = sn2(i, k, jj)
                     den = plt(i, k, jj)*100./(rgas*t(i, k, jj))
                     !$acc loop seq
                     do n = 1, nc
                        chat = chata(n, nn)
                        flux1 = cstar*den*a(i, jj)*(chat - uhat)/sn2r &
                                !    &           *abs((chat-uhat)/chat)**p
                                *abs((chat - uhat)/chat)
                        flux2 = 1./(1.+(mstar*(chat - uhat)/sn2r)**s3)
                        fxin(n) = flux1*flux2
                        fluxsat_ncr = fluxsat_ncr + flux1
                     end do
                     flux_ncr = 0.
      !
      !          if(j.eq.jlat .and. i.eq.ilon)then
      !          if(nn.eq.1)then
      !            print *,'bef norm fxin at lunch=',(fxin(n,nn),n=1,nc)
      !          endif
      !          endif
      !
                     !call normalize(fxtotal, nc, fxin(1, nn, jj), flux(1, i, k, nn, jj))
                     !
                     ! new sum should be added
                     !
               !
                     !call reassignsum(x,nx,newsum,y)
                     !call avevar(x,n,ave,var)
                     ave = 0.0
                     !$acc loop seq
                     do n = 1, nc
                       ave = ave + fxin(n)
                     end do
                     ave = ave/nc
                     var = 0.0
                     ep = 0.0
                     newsum = newsum/nc
                     !$acc loop seq
                     do n = 1, nc
                       ss = fxin(n) - ave
                       ep = ep + ss
                       var = var + ss*ss
                     end do
                     var = (var - ep**2/nc)/(nc - 1)
                     ! end call avevar
                     !$acc loop seq
                     do n = 1, nc
                       fluxr = fxin(n) - ave
                       fluxr = fluxr + newsum
                       flux(n, i, k, nn, jj) = fluxr
                     !end call reassignsum
                     !end call normalize
      !
      !          if(j.eq.jlat .and. i.eq.ilon)then
      !          if(nn.eq.1)then
      !            print *,'aft norm fxin at lunch=',(flux(n,i,k,nn),n=1,nc)
      !          endif
      !          endif
      !
                        flux_ncr = flux_ncr + fluxr
                     end do
                     flux_nc(i, lunch_lev, nn, jj) = flux_ncr
                     !!!!!!!!!!!!!!!!!!!!!
                  end if
               end do
            end do
         end do
         
         !$acc parallel loop async(async_id) collapse(3) gang vector private(ndir, &
         !$acc&         chat, flux1, flux2, newsum, ave, uv, uv1, den, &
         !$acc&         uhat_z1, uhat_z0, dflux1, fxin, k, uv_lunch, flux_ncr, &
         !$acc&         fluxr, fluxsat_ncr, sn2r, dfluxr, fluxtmp, chattmp, fluxs)
         do jj = 1, jlistnum
            do nn = 1, nphi ! 4
               do i = 1, nx
                  !$acc cache(chata)
                  if (i .le. myim(jj)) then
                     ndir = 1.
                     k = lunch_lev
                     uv_lunch = u(i, lunch_lev, jj)
                     if (nn .eq. 3 .or. nn .eq. 4) uv_lunch = v(i, lunch_lev, jj)
                     if (nn .eq. 2 .or. nn .eq. 4) ndir = -1.
                     uv1 = u(i, lunch_lev, jj)
                     if (nn .eq. 3 .or. nn .eq. 4) uv1 = v(i, lunch_lev, jj)
                     flux_ncr = flux_nc(i, lunch_lev, nn, jj)
                     ar = a(i, jj)
                     !$acc loop seq
                     do k = lunch_lev - 1, ltop, -1 ! 44 -> 1
      !
      ! compute spectral flux(eq-4) and saturated spectral flux(eq-6) Orr(2010)
      ! uhat_z1: upper level, uhat_z0: below level.
      ! for pw=1,p=2-pw=1, the formula of flux density is simply rewriten.
      !
                        if (nn .eq. 1 .or. nn .eq. 2) then
                           uv = u(i, k, jj)
                        else
                           uv = v(i, k, jj)
                        end if
                        uhat_z1 = uv - uv_lunch
                        uhat_z0 = uv1 - uv_lunch
                        uv1 = uv
                        newsum = flux_ncr/nc
                        flux_ncr = 0.
                        ave = 0.0
                        !var = 0.0
                        !ep = 0.0
                        sn2r = sn2(i, k, jj)
                        den = plt(i, k, jj)*100./(rgas*t(i, k, jj))
                        fluxtmp = cstar*den*ar/sn2r
                        dfluxr = 0.
                        !$acc loop seq
                        do n = 1, nc ! 20
                           chat = chata(n, nn)
                           chattmp = chat - uhat_z1
                           flux1 = fluxtmp*(chattmp) &
                                   !    &             *abs((chat-uhat_z1)/chat)**p
                                   *abs((chattmp)/chat)
                           flux2 = 1./(1.+(mstar*(chattmp)/sn2r)**s3)
                           fluxs(n) = flux1
                           fxin(n) = flux1*flux2
      !
                          ave = ave + fxin(n)
                        end do
                        ave = ave/nc
                        !$acc loop seq
                        do n = 1, nc ! 20
                          fluxr = fxin(n) - ave
                          fluxr = fluxr + newsum
                           
      !
      ! encounter critical level and do critical level filtering
      !
                           chat = chata(n, nn)
                           if (chat .le. uhat_z1 .and. chat .ge. uhat_z0) then
                              dflux1 = fluxr - flux(n, i, k + 1, nn, jj)
                              dfluxr = dfluxr + dflux1
                              fluxs(n) = fluxs(n) - dflux1
                              fluxr = fluxr - dflux1
                           else
                              fluxr = flux(n, i, k + 1, nn, jj)
                           end if
                           flux(n, i, k, nn, jj) = fluxr
                           fluxsat(n, i, k, nn, jj) = fluxs(n)
                        end do
                        dflux(i, k, nn, jj) = dfluxr
                     end do
                  end if
               end do
            end do
         end do
         !$acc parallel loop async(async_id) collapse(3) gang private(flux_ncr, &
         !$acc&         fluxsat_ncr, fluxr, dfluxr, ndir, dflux1, fluxsatr)
         do jj = 1, jlistnum
            do nn = 1, nphi
               do k = ltop, lunch_lev - 1
                  ndir = 1.
                  if (nn .eq. 2 .or. nn .eq. 4) ndir = -1.
                  !$acc loop vector
                  do i = 1, myim(jj)
                     flux_ncr = flux_nc(i, k, nn, jj)
                     fluxsat_ncr = fluxsat_nc(i, k, nn, jj)
                     dfluxr = dflux(i, k, nn, jj)
                     !$acc loop seq 
                     do n = 1, nc
      !
      ! nonlinear dispersion mechanism (saturation theorey)
      !
                        fluxr = flux(n, i, k, nn, jj)
                        fluxsatr = fluxsat(n, i, k, nn, jj)
                        dflux1 = fluxr - fluxsatr
                        if (ndir*dflux1 > 0) then
                           dfluxr = dfluxr + dflux1
                           fluxr = fluxsatr
                        end if
                        fluxsat_ncr = fluxsat_ncr + fluxsatr
                        flux_ncr = flux_ncr + fluxr
                     end do
                     flux_nc(i, k, nn, jj) = flux_ncr
                     fluxsat_nc(i, k, nn, jj) = fluxsat_ncr
                     dflux(i, k, nn, jj) = dfluxr
                  end do
               end do
            end do
         end do
   !
         
         
   !!$omp end parallel do
   !
         !$acc parallel loop async(async_id) collapse(2) gang 
         do jj = 1, jlistnum
            do k = ltop, lunch_lev - 1
               !$acc loop vector
               do i = 1, myim(jj)
                  dflux_u(i, k, jj) = dflux(i, k, 1, jj) + dflux(i, k, 2, jj)
                  dflux_v(i, k, jj) = dflux(i, k, 3, jj) + dflux(i, k, 4, jj)
               end do
            end do
         end do
   !      if(j.eq.jlat)then
   !          call findmaxf(nxj,nx,lev,dflux_u,vmax,iv,kv,vmin,ivm,kvm)
   !          print *,' dflux_u=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
   !          call findmaxf(nxj,nx,lev,dflux_v,vmax,iv,kv,vmin,ivm,kvm)
   !          print *,' dflux_v=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
   !      endif
   !
   !      if(j.eq.jlat)then
   !        print *,' orig u,v'
   !        call findmaxf(nxj,nx,lev,u,umax,iu,ku,umin,ium,kum)
   !        call findmaxf(nxj,nx,lev,v,vmax,iv,kv,vmin,ivm,kvm)
   !        print *,' max_u=',umax,iu,ku,' umin=',umin,ium,kum
   !        print *,' max_v=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
   !      endif
   !
         !$acc parallel loop async(async_id) collapse(2) gang vector
         do jj = 1, jlistnum
            do k = 1, lev
               avgdrag_u(k, jj) = 0.
               avgdrag_v(k, jj) = 0.
            end do
         end do
         !$acc parallel loop async(async_id) collapse(2) gang private(dpk, &
         !$acc&         wsold, wsnew, pk2r, pk2r1, p2, p21)
         do jj = 1, jlistnum
            do k = ltop + 1, lunch_lev - 1
               !$acc loop vector
               do i = 1, myim(jj)
                  pk2r = pk2(i, k, jj)
                  pk2r1 = pk2(i, k - 1, jj)
                  p2 = 1000.0*pk2r*pk2r*pk2r*sqrt(pk2r)
                  p21 = 1000.0*pk2r1*pk2r1*pk2r1*sqrt(pk2r1)
                  dpk = (p21 - p2)*100.
                  drag_u(i, k, jj) = -grav*(dflux_u(i, k - 1, jj) - dflux_u(i, k, jj))/dpk
                  drag_v(i, k, jj) = -grav*(dflux_v(i, k - 1, jj) - dflux_v(i, k, jj))/dpk
                  wsold = 0.5*(u(i, k, jj)*u(i, k, jj) + v(i, k, jj)*v(i, k, jj))
      !        wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
                  u(i, k, jj) = u(i, k, jj) + dt*drag_u(i, k, jj)
                  v(i, k, jj) = v(i, k, jj) + dt*drag_v(i, k, jj)
      !        wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
      !        t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
                  wsnew = 0.5*(u(i, k, jj)*u(i, k, jj) + v(i, k, jj)*v(i, k, jj))
                  t(i, k, jj) = t(i, k, jj) + (wsold - wsnew)/cp
                  !$acc atomic
                  avgdrag_u(k, jj) = avgdrag_u(k, jj) + drag_u(i, k, jj)*cosljj(jj)
                  !$acc atomic
                  avgdrag_v(k, jj) = avgdrag_v(k, jj) + drag_v(i, k, jj)*cosljj(jj)
               end do
            avgdrag_u(k, jj) = avgdrag_u(k, jj)/wtcosl(jj)
            avgdrag_v(k, jj) = avgdrag_v(k, jj)/wtcosl(jj)
            end do
         end do
         k = 1
         !$acc parallel loop async(async_id) gang private(dpk, wsold, wsnew, pk2r, p2)
         do jj = 1, jlistnum
            !$acc loop vector
            do i = 1, myim(jj)
               pk2r = pk2(i, k, jj)
               p2 = 1000.0*pk2r*pk2r*pk2r*sqrt(pk2r)
               wsold = 0.5*(u(i, k, jj)*u(i, k, jj) + v(i, k, jj)*v(i, k, jj))
   !        wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
               dpk = (ptop - p2)*100.
               drag_u(i, k, jj) = grav*dflux_u(i, k, jj)/dpk
               drag_v(i, k, jj) = grav*dflux_v(i, k, jj)/dpk
               u(i, k, jj) = u(i, k, jj) + dt*drag_u(i, k, jj)
               v(i, k, jj) = v(i, k, jj) + dt*drag_v(i, k, jj)
   !        u(i,k)=u(i,k)+dt*drag_u(i,k+1)*2.
   !        v(i,k)=v(i,k)+dt*drag_v(i,k+1)*2.
   !        wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
   !        t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
               wsnew = 0.5*(u(i, k, jj)*u(i, k, jj) + v(i, k, jj)*v(i, k, jj))
               t(i, k, jj) = t(i, k, jj) + (wsold - wsnew)/cp
               !$acc atomic
               avgdrag_u(k, jj) = avgdrag_u(k, jj) + drag_u(i, k + 1, jj)*2.*cosljj(jj)
               !$acc atomic
               avgdrag_v(k, jj) = avgdrag_v(k, jj) + drag_v(i, k + 1, jj)*2.*cosljj(jj)
            end do
            avgdrag_u(k, jj) = avgdrag_u(k, jj)/wtcosl(jj)
            avgdrag_v(k, jj) = avgdrag_v(k, jj)/wtcosl(jj)
         end do
         !$acc end data
!
!      if(j.eq.jlat)then
!          print *,' drag u,v'
!          call findmaxf(nxj,nx,lev,u,umax,iu,ku,umin,ium,kum)
!          call findmaxf(nxj,nx,lev,v,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' max_u=',umax,iu,ku,' umin=',umin,ium,kum
!          print *,' max_v=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!      endif
!
         return
      end
!
