      subroutine pbl_noah_gpu(nxjp, nx, lev, ktpbl, dt, g, r, cp, xkapa &
                              , hltm, ptop, tice, hice, tg, z0, land &
                              , topo, phi, phii, pss, u, v, t, q, ut, vt &
                              , tt, qt, pk, pk2, ustar, tstar, qstar, e &
                              , eps, hflux, qflux, itstp, gwclim, tgclim &
                              , ocean, ice, sheleg, totalp, ss, rs, alb &
                              , imx, xkmx, idg, xkmd, itype, t2, q2, rh2 &
                              , rh10, u10, v10, fm, fh, fm10, fh2, srflag &
                              , rld, stbo, km, smc, stc, canopy, runoff &
                              , sigmaf, istyp, ivegtyp, ncld, dsigma &
                              , islopetyp, slc, sncover, snwdph, shdmax &
                              , shdmin, snoalb, albedo2, sld, zice, cice &
                              , xtice, hpbl, asl, atl, xmu, gfx, kpbl &
                              , nmpbl, nmmiph, isot, ivegsrc, sfemis_g &
                              , dudt, dvdt, dtdt, dqdt)
!
!#######################################################################
!                     subroutine description
!
! 1. function description
!
!    this version is based on the Troen and Mahrt (1986) to deal with
!    the vertical diffusion with account of non-local flux.
!    the code is origially from ncepgfs
!
!
! 2. parameter specification
!
!     nx : dimension of horizontal-direction
!    lev : dimension of z-direction
!  ktpbl : starting level to consider vertical mixing by pbl
!     dt : time step for ut, vt, tt, qt updates                (s)
!      g : gravity                                             (m/s2)
!      r : air gas constant                                    (j/kg)
!     cp : air specific heat constant                          (j/kg)
!  xkapa : r/cp
!   hltm : water vapor latent heat constant                    (j/kg)
!   ptop : model top pressure                                  (mb)
!   tice : temperature of melting point                        (k)
!   hice : latent heat of snow melting                         (j/kg)
!     tg : ground real temp. (sea or land)  (nx)               (k)
!     z0 : surface roughness length    (nx)                    (m)
!   gwet : ground wetness              (nx)
!   land : logic for bare soil land    (nx)  (=true for bare land)
!   topo : terrain geopotential   (nx)                         (m2/s2)
!    phi : odd level geopotential (nx)                         (m2/s2)
!    pss : (terrain pres -ptop)   (nx)     at current time lvl (mb)
!      u : x-direction velocity   (nx,lev) at current time lvl (m/s)
!      v : y-direction velocity   (nx,lev) at current time lvl (m/s)
!      t : potent temp on odd lvl (nx,lev) at current time lvl (k)
!      q : moisture    on odd lvl (nx,lev) at current time lvl (kg/kg)
!     ut : x-direction velocity   (nx,lev) at past/future  lvl (m/s)
!     vt : y-direction velocity   (nx,lev) at past/future  lvl (m/s)
!     tt : real temp.  on odd lvl (nx,lev) at past/future  lvl (k)
!     qt : moisture    on odd lvl (nx,lev) at past/future  lvl (kg/kg)
!     pk : p**capa on odd levels  (exner func)
!    pk2 : p**capa on even levels (exner func)
!   ustar: friction velocity         (nx)                      (m/s)
!   tstar: potential temp scale      (nx)                      (k)
!   qstar: moisture scale            (nx)                      (kg/kg)
!     e  : turbulence kinetic energy (nx,lev) at even lvls     (m2/s2)
!    eps : tke dissipation           (nx,lev)                  (m2/s2)
!   hflux: upward surface heat flux         (nx)               (w/m2)
!   qflux: upward surface moisture flux     (nx)               (w/m2)
!   itstp: indicated which time step it is
!    gwr : soil water amount           (nx)                    (mm)
!  gwclim: soil water climate          (nx)                    (mm)
!  tgclim: climate ground temperature  (nx)                    (k)
!   ocean: logic for open water      (nx)  (=true for open water)
!     ice: logic for ice cover       (nx)  (=true for ice covered)
!     snr: snow depth in water state (nx)                      (mm)
!  totalp: total precipitation       (nx)                      (mm)
!      ss: net short wave radiat. flux down to ground (nx)     (wat/m2)
!      rs: net long  wave radiat. flux up from ground (nx)     (wat/m2)
!     alb: diffusion albedo           (nx)
!     imx: address of the two points with max xkm at the two lowest lvl
!    xkmx: max mixing coef xkm at two lowest levels            (m2/s)
!     idg: i-index address of selected point to print profile
!    xkmd: xkm profile at the selected idg point               (m2/s)
!   itype: surface layer update type, =0 use surface stress as b.c.
!                                     =1 directly update lowest levels
!      t2: 2m air temp derived by similarity theory (nx)     (K)
!     u10: 10m air x-direction wind derived by simil. theory(nx)  (m/s)
!     v10: 10m air y-direction wind derived by simil. theory(nx)  (m/s)
!
!soil io
!
!     smc: volumetric soil moisture content (nx,km_soil)
!     stc: soil temperature                 (nx,km_soil)
!  canopy: canopy moisture content(<0.5mm)          (nx)
!  sigmaf: green vegetation fraction        (nx)
!   istyp: soil type(1-9)                   (nx)
! ivegtyp: vegetation type(1-13)            (nx)
!  runoff: accumulate run off water(runof+drain) (nx)             (mm)
!     rld: long wave radiat. flux down to ground (nx)         (wat/m2)
!   wlt  : wiltling point of 9 soil types  (9)
!   ref  : field capacity of 9 soil types  (9)
! tsat   : satuation point of 9 soil types (9)
! dfkt   : soil thermal diffusivity     (22,9)     m2/s
! xktk   : soil hydraulic conductivity  (22,9)     m/s
!  dfk   : soil hydraulic diffusivity   (22,9)     m2/s
!
! 4. local variable
!
!    phi : geopential on odd lvl  (nx,lev) at current time lvl (m2/s2)
!    xkm : eddy mixing coef. for u,v (nx,lev)                  (m2/s)
!    xkh : eddy mixing coef. for t,q (nx,lev)                  (m2/s)
!      ts: surface air temperature    (nx)                     (k)
!    qsfc: surface moisture (specific humidity)      (nx)      (kg/kg)
!     zl : height(lev) over monin-obukhov scale height.
!          a stability indicator.                    (nx)
!     czh: ground thermal capacity                   (nx)      (j/kg.m2)
!    snow: logic for snow cover                      (nx)
!
!soil
!
!  rhscnpy: righ hand side terms of canopy prediction(nx)
!    rhsmc: righ hand side terms of smc prediction(nx,km)
!      aim: working array smc prediction           (nx,km)
!      bim: working array smc prediction           (nx,km)
!      cim: working array smc prediction           (nx,km)
!   drain: drainage water from bottom soil layer   (nx)        (mm/s)
!     snomt: snow melt rate                        (nx)        (mm/s)
!     zsoil: soil depth of soil layer(-0.1,-2.)    (nx,km)        (m)
!     runof: run off water from surface and soil   (nx)        (mm/s)
!      gfx: ground heat flx                        (nx)        (wat/m2)
!     t850: 850mb air temp                         (nx)         (K)
!     gesh: control number,control precipitation   (nx)
!           conver to snow(0) or not(1)
!
!
! 5. calling modules
!
!     main program
!
!     created       november 1991
!     modified      july     1992,  april   1993
!     modified      jan      2005  by f. wang
!     modified      may      2009  by f. wang
!
!            ed., american meteorology socity, 67-98.
!  mixpbl:
!       troen, i., and l. mahrt, 1986: a simple model of the atmospheric
!          boundary layer: sensitivity to surface evaporation. bound.
!          layer meteor., 37, 129-148.
!
!  land process:
!       Mahrt, L., and H. -L. Pan,1984: Atwo-layer model of soil hydrology.
!             Boundary Layer Meteorol., 29, 1-20.
!       Pan, H.-L., and L.Mahrt, 1987: Interaction between soil hydrology
!             and boundary-layer development. Boundary Layer Meteorol.,
!             38,185-202.
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#####################################################################
!
         use mpe
         use rank
         use index
         use radn, only: ntcw, ntiw, ntinc, ntrnc, ntoz, ntrw, ntsw, ntgl
         use const, only: RTYPE
         use param, only: my, my_max
         use cusparse
!ch   use paramt

!
         implicit none
!
! input & output variable
!
         integer nxj, nxjp(my), nx, lev, ktpbl, itype, idg(my), ncld, km, &
            nmpbl, nmmiph, nc, myim(my_max)
         real dt, g, r, cp, xkapa, hltm, ptop, stbo, hice, tice

         integer imx(2, my), itstp

         real tg(nx, my_max), z0(nx, my_max), &
            u(nx, lev, my_max), v(nx, lev, my_max), t(nx, lev, my_max), &
            ustar(nx, my_max), tstar(nx, my_max), qstar(nx, my_max), &
            e(nx, lev, my_max), eps(nx, lev, my_max), hflux(nx, my_max), &
            qflux(nx, my_max), pkd(nx, my_max), pk2d(nx, my_max), gwclim(nx, my_max), &
            tgclim(nx, my_max), snr(nx), totalp(nx, my_max), &
            ss(nx, my_max), rs(nx, my_max), alb(nx, my_max), xkmx(2, my), &
            xkmd(lev), t2(nx, my_max), u10(nx, my_max), v10(nx, my_max)
         real(kind=RTYPE) qt(nx, lev*ncld, my_max), q(nx, lev*ncld, my_max), &
            phi(nx, lev, my_max), topo(nx, my_max), pss(nx, my_max), &
            ut(nx, lev, my_max), vt(nx, lev, my_max), tt(nx, lev, my_max), &
            pk(nx, lev, my_max), pk2(nx, lev, my_max)
!soil
         real smc(nx, km, my_max), stc(nx, km, my_max), canopy(nx, my_max), &
            sigmaf(nx, my_max), rld(nx, my_max), runoff(nx, my_max)

         integer istyp(nx, my_max), ivegtyp(nx, my_max)
!
         logical land(nx, my_max), ocean(nx, my_max), ice(nx, my_max)
!
!  local work arrays
!
!ch   real     hgt(im,lm),xkm(im,lm),xkh(im,lm),ts(im),                    &
!ch            qsfc(im),zl(im),czh(im),ps(im),sfcw(im),                    &
!ch            dhgt(im,lm),ro2(im,lm),dhgtz(im,lm)
         real hgt(nx, my_max), xkm(nx, lev), xkh(nx, lev), &
            qsfc(nx), zl(nx), czh(nx), ps(nx, my_max), sfcw(nx, my_max), &
            dhgt, ro2(nx, my_max), dhgtz(nx, lev)
!soil
         real rhscnpy(nx), rhsmc(nx, km), aim(nx, km), bim(nx, km), &
            cim(nx, km), drain(nx, my_max), snomt(nx, my_max), zsoil(nx, km), &
            runof(nx, my_max), gfx(nx, my_max), t850(nx, my_max), gesh(nx)
!
!ch   logical snow(im)
         logical snow(nx)
!
! ncep_pbl
!ch   integer   kpbl(im)
!ch   real      rb(im),fm(im),fh(im),hpbl(im),                          &
!ch             heat(im),evap(im),stress(im),                           &
!ch             prsl(im,lm),prslk(im,lm),phil(im,lm),del(im,lm),        &
!ch             prsi(im,lm+1),phi2(im,lm+1),phii(im,lm+1),              &
!ch             dsigma(lm,2),rcl(im),                                   &
         integer kpbl(nx, my_max)
         real rb(nx, my_max), fm(nx, my_max), fh(nx, my_max), hpbl(nx, my_max), &
            heat(nx, my_max), evap(nx, my_max), stress(nx, my_max), &
            prsl(nx, lev, my_max), prslk(nx, lev, my_max), &
            phil(nx, lev, my_max), del(nx, lev, my_max), &
            !byl                prsi(nx,lev+1),phi2(nx,lev+1),phii(nx,lev+1),              &
            prsi(nx, lev + 1, my_max), phii(nx, lev + 1, my_max), &
            rcl(nx, my_max), &
            u1(nx, lev, my_max), v1(nx, lev, my_max), t1(nx, lev, my_max)
         real(kind=RTYPE) dsigma(lev, 2)
         real, dimension(:, :, :, :), allocatable :: q1
!
         real pk2x(nx, lev, my_max), pkx(nx, lev, my_max)
!
! noah  ------------
! input and output
         integer islopetyp(nx, my_max), islmsk(nx, my_max), io(my_max), jo(my_max)
         real sld(nx, my_max), &
            slc(nx, km, my_max), &
            sheleg(nx, my_max), snwdph(nx, my_max), sncover(nx, my_max), &
            zice(nx, my_max), cice(nx, my_max), xtice(nx, my_max), &
            rh2(nx, my_max), rh10(nx, my_max), &
            shdmin(nx, my_max), shdmax(nx, my_max), snoalb(nx, my_max), &
            albedo2(nx, my_max), sfemis_g(nx, my_max)
! local
         real z0rl(nx, my_max), prslki(nx, my_max), cd(nx, my_max), cdq(nx, my_max), &
            tsurf(nx, my_max), psi(nx, my_max), prsl1(nx, my_max), &
            fm10(nx, my_max), fh2(nx, my_max), fh10(nx, my_max), &
            qsurf(nx, my_max), evapc(nx), cmm(nx), chh(nx), ep1d(nx, my_max), &
            !                radsl(nx) ,tprcp(nx),                                   &
            tprcp(nx, my_max), &
            !                phy_f2d(nx),q2(nx)
            ddvel(nx, my_max), q2(nx, my_max)
!
         real srflag(nx, my_max), sfemis(nx, my_max)
!orig logical   flag_guess(nx),flag_iter(nx),mom4ice(nx)
         logical flag_guess(nx, my_max), flag_iter(nx, my_max), mom4ice(my_max)
!
!2018 new
         logical redrag(my_max)
         integer ivegsrc, isot
!----------------------------------------------------
! for new pbl: asl,atl,xmu
         real asl(nx, lev, my_max), atl(nx, lev, my_max), swh(nx, lev, my_max), &
            hlw(nx, lev, my_max), xmu(nx, my_max)

!----------------------------------------------------
! for fractional step:
         real dudt(nx, lev, my_max), dvdt(nx, lev, my_max), &
            dtdt(nx, lev, my_max), dqdt(nx, lev, my_max)

         integer lsm(my_max), i, k, iter, kc, ntrac, jj, j
         real ppd, ppp, ttt, ppu, dth(my_max), p850, ddd, cc, qqq
         real time1, time2
         integer :: istat, async_id = 1
         type(cusparseHandle) :: sparsehandle
!
! noah mode
         istat = cusparseCreate(sparsehandle)
         do jj = 1, jlistnum
            lsm(jj) = 1
            mom4ice(jj) = .false.
            redrag(jj) = .false.
            io(jj) = 480
            jo(jj) = 240
            dth(jj) = dt
         end do
!------------------------------------------------------------
!     io=67
!     jo=251
!     io=589
!     jo=148
!     io=705
!     jo=236
!
         ntrac = ncld
         if (nmmiph .eq. 6) ntrac = ncld - 3
         if (nmmiph .eq. 8) ntrac = ncld - 4
         if (nmmiph .eq. 18) ntrac = 5
         if (nmmiph .eq. 11 .or. nmmiph .eq. 12 .or. nmmiph .eq. 13) ntrac = 7
         if (nmmiph .eq. 15) ntrac = 5
         if (nmmiph .eq. 16) ntrac = 5
!      if ( nmmiph .eq.16 ) ntrac=8

!
! --- ensure ktpbl selection is greater than 2
!
         if (ktpbl .lt. 2) ktpbl = 2
!
! --- compute surface pres and surface air temp at current time level
!
         allocate (q1(nx, lev, ntrac, my_max))

         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
         end do

         do jj = 1, jlistnum
            do k = 1, lev
               pkd(:, jj) = pk(:, k, jj)
               pk2d(:, jj) = pk2(:, k, jj)
               call vlog(pk2x(1, k, jj), pk2d(1, jj), myim(jj))
               call vlog(pkx(1, k, jj), pkd(1, jj), myim(jj))
               do i = 1, myim(jj)
                  pk2x(i, k, jj) = pk2x(i, k, jj)*(cp/r)
                  pkx(i, k, jj) = pkx(i, k, jj)*(cp/r)
               end do
               call vexp(pk2x(1, k, jj), pk2x(1, k, jj), myim(jj))
               call vexp(pkx(1, k, jj), pkx(1, k, jj), myim(jj))
            end do
         end do
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               ps(i, jj) = pss(i, jj) + ptop
               hgt(i, jj) = (phi(i, lev, jj) - topo(i, jj))/g
            end do
         end do
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               qqq = max(qt(i, lev, jj), 1.0e-8)
               ttt = tt(i, lev, jj)*(1.+0.608*qqq)
               ppp = pkx(i, lev, jj)*1000.
               ro2(i, jj) = 100.0*ppp/(r*ttt)
            end do
         end do
         !
         !      do 130 k = 1, lev-1
         !      do 130 i = 1, nxj
         !      ppu = pkx(i,k)*1000.
         !      ppd = pkx(i,k+1)*1000.
         !      ro2(i,k) = ( ppd - ppu ) * 100. / ( phi(i,k) - phi(i,k+1) )
         !  130 continue
         !

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !                                                                      c
         !     determine ground condition and predict ghround temp and wetness  c
         !                                                                      c
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !
         !
         !
         !                                                                      c
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !
         ! NOAH land soil model ....................................
         !
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !
         !..........................................
         ! define some variable needed in noah sub
         !..........................................
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               ! build slimsk table for input
               if (ocean(i, jj)) islmsk(i, jj) = 0
               if (land(i, jj)) islmsk(i, jj) = 1
               if (ice(i, jj)) islmsk(i, jj) = 2
               ! build roughness from ustar over ocean for input
               if (itstp .eq. 1 .and. ocean(i, jj)) &
                  z0(i, jj) = ustar(i, jj)*ustar(i, jj)*0.014/g ! get z0 from ustar
               !
               prslki(i, jj) = pk2(i, lev, jj)/pk(i, lev, jj)   !(ps/p1)**r/cp
               !     sfemis   - real, sfc lw emissivity (fractional)
               sfemis(i, jj) = sfemis_g(i, jj)
               !        radsl(i)           =-ss(i)-rld(i)  ! snet + rld  upward
               !
               tprcp(i, jj) = totalp(i, jj)/1000.  ! dth precip (m)
               !
               !        radsl(i)            =-ss(i)-rld(i)  ! snet + rld  upward
               !
               !       shdmin(i)           =0.01
               !       shdmax(i)           =0.99
               !       snoalb(i)           =0.6
               !        phy_f2d(i)          =0.
               ddvel(i, jj) = 0.
               rcl(i, jj) = 1.
            end do
         end do
         !..............................................................
         !  some variable pass from input and may be changed in this sub
         !  so these variables should be put back to input
         !.............................................................
         !      do k=1,lev
         !       do i=1,nx
         !       t(i,k)              =t(i,k)*pk(i,k) ! transfer t to temp
         !       enddo
         !      enddo

         !     io=67
         !     if(jj.eq.jo)then
         !     print*,'before noah t(in  tmpt)-----------------'
         !     do k=1,lev
         !     print*,k,tt(io,k),pk(io,k)
         !     enddo
         !     endif
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               z0rl(i, jj) = z0(i, jj)*100.         ! transport from m to cm
               !       sheleg(i)           =snr(i)    !snow equiv water depth (mm)
               !       snwdph(i)           =snr(i)*10. !snow real depth (mm)
               !        psi(i)              =ps(i)*0.1  ! surface pressure (mb to cb)
               !        prsl(i,1)           =pkx(i,lev)*100.     !cb
               !noah new
               psi(i, jj) = ps(i, jj)*100.  ! surface pressure (mb to pa)
               prsl1(i, jj) = pkx(i, lev, jj)*1000.*100.     !pa
            end do
         end do

         !...........................
         ! initialize  variable
         !..........................
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               tsurf(i, jj) = tg(i, jj)
               flag_guess(i, jj) = .false.
               flag_iter(i, jj) = .true.
               drain(i, jj) = 0.
               ep1d(i, jj) = 0.
               runof(i, jj) = 0.
               hflux(i, jj) = 0.
               qflux(i, jj) = 0.
            end do
         end do
         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !        to get environment temp at 850mb
         !
         p850 = 850.
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               t850(i, jj) = tt(i, lev, jj)
            end do
         end do
         !
         do jj = 1, jlistnum
            do k = lev, 2, -1
               do i = 1, myim(jj)
                  ppd = pkx(i, k, jj)*1000.
                  ppu = pkx(i, k - 1, jj)*1000.
                  ddd = (p850 - ppd)*(p850 - ppu)
                  if (ddd .lt. 0.) then
                     cc = (log(p850/ppd))/(log(ppu/ppd))
                     t850(i, jj) = (tt(i, k - 1, jj) - tt(i, k, jj))*cc &
                                   + tt(i, k, jj)
                  end if
               end do
            end do
         end do
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               srflag(i, jj) = 0.
               if (t850(i, jj) .le. 273.16) then
                  srflag(i, jj) = 1.
               end if
            end do
         end do

         !.............................................
         ! loop needed to remove unstable in calm situation (sfcw < 2m/s)
         !............................................
         do iter = 1, 2
            do jj = 1, jlistnum
               !
               !surface exchange coefficient
               !
               call sfc_diff(myim(jj), nx, psi(1, jj), ut(1, lev, jj), &
                             vt(1, lev, jj), tt(1, lev, jj), &
                             qt(1, lev, jj), &
                             !                    tg,z0rl,cd,cdq,rb,                                  &
                             !                    rcl,prsl(1,1),work1,slmsk,                          &
                             !                    ustar,sfcw,phy_f2d,fm10,fh2,                        &
                             !                    sigmaf,ivegtyp,shdmax,                              &
                             !                    tsurf,flag_iter)
                             hgt(1, jj), snwdph(1, jj), tg(1, jj), z0rl(1, jj), cd(1, jj), &
                             cdq(1, jj), rb(1, jj), prsl1(1, jj), prslki(1, jj), islmsk(1, jj), &
                             stress(1, jj), fm(1, jj), fh(1, jj), &
                             ustar(1, jj), sfcw(1, jj), ddvel(1, jj), fm10(1, jj), &
                             fh2(1, jj), fh10(1, jj), sigmaf(1, jj), &
                             ivegtyp(1, jj), shdmax(1, jj), ivegsrc, &
                             tsurf(1, jj), flag_iter(1, jj), redrag(jj))
               !
               !     print*,'pblnoah,diff'
            end do
            do jj = 1, jlistnum
               do i = 1, myim(jj)
                  if ((iter .eq. 1) .and. (sfcw(i, jj) .lt. 2.)) &
                     flag_guess(i, jj) = .true.
               end do
            end do
            !
            ! surface flux over ocean
            !
            !      call sfc_ocean(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev), &
            !                     tg,qsurf,evapc,gfx,cd,cdq,rcl,prsl(1,1),work1,     &
            !                     slmsk,qflux,hflux,ep1d,phy_f2d,                    &
            !                     flag_iter)
            do jj = 1, jlistnum
               call sfc_ocean(myim(jj), nx, psi(1, jj), ut(1, lev, jj), &
                              vt(1, lev, jj), tt(1, lev, jj), &
                              qt(1, lev, jj), tg(1, jj), cd(1, jj), cdq(1, jj), &
                              prsl1(1, jj), prslki(1, jj), islmsk(1, jj), &
                              ddvel(1, jj), flag_iter(1, jj), &
                              qsurf(1, jj), gfx(1, jj), qflux(1, jj), &
                              hflux(1, jj), ep1d(1, jj))
            end do
            !
            !     print*,'pblnoah,ocean'
            !
            !     if(jj.eq.jo)then
            !     print*,'before sfcdrv-----------------'
            !     print*,'  k     u     v      t        q'
            !     do k=1,lev
            !     print*,k,ut(io,k),vt(io,k),tt(io,k),qt(io,k)
            !     enddo
            !     print*,'sheleg sncover snwdph=',sheleg(io),sncover(io),snwdph(io)
            !     print*,'istyp,ivegtyp,sigmaf=',istyp(io),ivegtyp(io),sigmaf(io)
            !     print*,'-----------------------------'
            !     print*,'ss,rld,sld,radsl=',ss(io),rld(io),sld(io),radsl(io)
            !     print*,'prsl,work1,slmsk=',prsl(io,1),work1(io),slmsk(io)
            !     print*,'tg,tsurf=',tg(io),tsurf(io)
            !     print*,'smc     slc      stc='
            !     do k=km,1,-1
            !     print*,k,smc(io,k),slc(io,k),stc(io,k)
            !     enddo
            !     endif
            !
            !      call sfc_drv(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),&
            !                     sheleg,sncover,snwdph,tg,qsurf,tprcp,SRFLAG,       &
            !                     smc,stc,slc,evapc,istyp,sigmaf,                    &
            !                     ivegtyp,canopy,rld,sld,                            &
            !                     radsl,dth,tgclim,gfx,cd,cdq,                       &
            !                     rcl,prsl(1,1),work1,slmsk,                         &
            !                     drain,qflux,hflux,ep1d,phy_f2d,                    &
            !                     runof,SLOPETYP,SHDMIN,SHDMAX,SNOALB,alb,           &
            !                     tsurf, flag_iter, flag_guess,albedo2,jj,io,jo)
            do jj = 1, jlistnum
               j = jlist1(jj)
               call sfc_drv(myim(jj), nx, km, psi(1, jj), ut(1, lev, jj), &
                            vt(1, lev, jj), tt(1, lev, jj), &
                            qt(1, lev, jj), istyp(1, jj), &
                            ivegtyp(1, jj), sigmaf(1, jj), sfemis(1, jj), &
                            rld(1, jj), sld(1, jj), ss(1, jj), dth(jj), &
                            tgclim(1, jj), cd(1, jj), cdq(1, jj), prsl1(1, jj), prslki(1, jj), &
                            hgt(1, jj), islmsk(1, jj), ddvel(1, jj), islopetyp(1, jj), &
                            shdmin(1, jj), shdmax(1, jj), &
                            snoalb(1, jj), alb(1, jj), flag_iter(1, jj), &
                            flag_guess(1, jj), isot, ivegsrc, sheleg(1, jj), &
                            snwdph(1, jj), tg(1, jj), tprcp(1, jj), &
                            srflag(1, jj), smc(1, 1, jj), &
                            stc(1, 1, jj), slc(1, 1, jj), &
                            canopy(1, jj), tsurf(1, jj), z0rl(1, jj), sncover(1, jj), &
                            qsurf(1, jj), gfx(1, jj), drain(1, jj), qflux(1, jj), &
                            hflux(1, jj), ep1d(1, jj), runof(1, jj), &
                            albedo2(1, jj), j, io(jj), jo(jj))
            end do
            !
            !       call sfc_sice(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),   &
            !                      zice,cice,xtice,sld,        &    ! FOR SEA-ICE - XW Nov04
            !                      sheleg,snwdph,tg,qsurf,tprcp,SRFLAG,stc,evapc,    &
            !                      rld,radsl,SNOMT,dth,gfx,cd,cdq,                   &
            !                      rcl,prsl(1,1),work1,slmsk,                        &
            !                      qflux,hflux,ep1d,phy_f2d,flag_iter,               &
            !                      mom4ice,lsm)
            do jj = 1, jlistnum
               call sfc_sice(myim(jj), nx, km, psi(1, jj), ut(1, lev, jj), &
                             vt(1, lev, jj), tt(1, lev, jj), &
                             qt(1, lev, jj), dth(jj), sfemis(1, jj), rld(1, jj), &
                             ss(1, jj), sld(1, jj), srflag(1, jj), &
                             cd(1, jj), cdq(1, jj), prsl1(1, jj), prslki(1, jj), &
                             islmsk(1, jj), ddvel(1, jj), &
                             flag_iter(1, jj), mom4ice(jj), lsm(jj), zice(1, jj), &
                             cice(1, jj), xtice(1, jj), sheleg(1, jj), &
                             tg(1, jj), tprcp(1, jj), stc(1, 1, jj), ep1d(1, jj), &
                             snwdph(1, jj), qsurf(1, jj), snomt(1, jj), gfx(1, jj), &
                             qflux(1, jj), hflux(1, jj))
            end do
            !
            do jj = 1, jlistnum
               do i = 1, myim(jj)
                  flag_iter(i, jj) = .False.
                  flag_guess(i, jj) = .False.
                  if ((islmsk(i, jj) .eq. 1) .and. (iter .eq. 1)) then
                     if (sfcw(i, jj) .lt. 2.) flag_iter(i, jj) = .true.
                  end if
               end do
            end do                                 !!!!! <---- Clu_q2m_iter
         end do
         !

         !** update near surface fields
         !
         !      call sfc_diag(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev), &
         !                    tg,qsurf,u10,v10,t2,q2,rcl,work1,slmsk,             &
         !                    qflux,fm,fh,fm10,fh2,fh10,rh2,rh10)
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxjp(j)
            call sfc_diag(myim(jj), nx, psi(1, jj), ut(1, lev, jj), vt(1, lev, jj), &
                          tt(1, lev, jj), qt(1, lev, jj), tg(1, jj), &
                          qsurf(1, jj), u10(1, jj), v10(1, jj), t2(1, jj), &
                          q2(1, jj), prslki(1, jj), qflux(1, jj), fm(1, jj), &
                          fh(1, jj), fm10(1, jj), fh2(1, jj), fh10(1, jj), &
                          rh2(1, jj), rh10(1, jj))
         end do

         !

         !     if(jj.eq.jo)then
         !     print*,'after noah -----------------'
         !     print*,'k       tt       t2         q2        rh2-----------------'
         !     print*,k,tt(io,lev),t2(io),q2(io),rh2(io)
         !     print*,'qflux,hflux=',qflux(io),hflux(io)
         !     endif
         !      do k=1,lev
         !      do i=1,nx
         !       t(i,k)          =t(i,k)/pk(i,k)   ! transfer t back to poten temp
         !      enddo
         !      enddo
         !

         do jj = 1, jlistnum
            do i = 1, myim(jj)
               z0(i, jj) = z0rl(i, jj)/100.         ! transport from cm to m
            end do
         end do
         !soil
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !                                                                      c
         !     mixing layer parameterization:                                   c
         !     compute vertical eddy mixing coeff. and update all level         c
         !     u, v, t and q according to the mixing and surface drag           c
         !                                                                      c
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !
         !
         if (nmmiph .eq. 6) then  !WSM6
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     q1(i, kc, 1, jj) = qt(i, k, jj)
                     q1(i, kc, 2, jj) = qt(i, lev*(ntcw - 1) + k, jj)
                     q1(i, kc, 3, jj) = qt(i, lev*(ntiw - 1) + k, jj)
                     q1(i, kc, 4, jj) = qt(i, lev*(ntoz - 1) + k, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 8) then ! Thompson
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     q1(i, kc, 1, jj) = qt(i, k, jj)
                     q1(i, kc, 2, jj) = qt(i, lev*(ntcw - 1) + k, jj)
                     q1(i, kc, 3, jj) = qt(i, lev*(ntiw - 1) + k, jj)
                     q1(i, kc, 4, jj) = qt(i, lev*(ntinc - 1) + k, jj)
                     q1(i, kc, 5, jj) = qt(i, lev*(ntoz - 1) + k, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 18) then ! 2M Thompson
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     q1(i, kc, 1, jj) = qt(i, k, jj)
                     q1(i, kc, 2, jj) = qt(i, lev*(ntcw - 1) + k, jj)
                     q1(i, kc, 3, jj) = qt(i, lev*(ntiw - 1) + k, jj)
                     q1(i, kc, 4, jj) = qt(i, lev*(ntinc - 1) + k, jj)
                     q1(i, kc, 5, jj) = qt(i, lev*(ntoz - 1) + k, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 11 .or. nmmiph .eq. 12 .or. nmmiph .eq. 13) then ! GFDL MP
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     q1(i, kc, 1, jj) = qt(i, k, jj)
                     q1(i, kc, 2, jj) = qt(i, lev*(ntcw - 1) + k, jj)
                     q1(i, kc, 3, jj) = qt(i, lev*(ntiw - 1) + k, jj)
                     q1(i, kc, 4, jj) = qt(i, lev*(ntrw - 1) + k, jj)
                     q1(i, kc, 5, jj) = qt(i, lev*(ntsw - 1) + k, jj)
                     q1(i, kc, 6, jj) = qt(i, lev*(ntgl - 1) + k, jj)
                     q1(i, kc, 7, jj) = qt(i, lev*(ntoz - 1) + k, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 15 .or. nmmiph .eq. 16) then ! Goddard MP
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     q1(i, kc, 1, jj) = qt(i, k, jj)
                     q1(i, kc, 2, jj) = qt(i, lev*(ntcw - 1) + k, jj)
                     q1(i, kc, 3, jj) = qt(i, lev*(ntiw - 1) + k, jj)
                     q1(i, kc, 4, jj) = qt(i, lev*(ntrw - 1) + k, jj)
                     q1(i, kc, 5, jj) = qt(i, lev*(ntoz - 1) + k, jj)
                  end do
               end do
            end do
         else
            do jj = 1, jlistnum
               do nc = 1, ntrac
                  do k = 1, lev
                     kc = lev - k + 1
                     do i = 1, myim(jj)
                        q1(i, kc, nc, jj) = qt(i, lev*(nc - 1) + k, jj)
                     end do
                  end do
               end do
            end do
         end if
         !
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               do i = 1, myim(jj)
                  t1(i, kc, jj) = tt(i, k, jj)
                  u1(i, kc, jj) = ut(i, k, jj)
                  v1(i, kc, jj) = vt(i, k, jj)
                  !
                  prslk(i, kc, jj) = pk(i, k, jj)
                  prsl(i, kc, jj) = pkx(i, k, jj)*100000.     !Pa
                  prsi(i, kc, jj) = pk2x(i, k, jj)*100000.    !Pa
                  del(i, kc, jj) = (pss(i, jj)*dsigma(k, 1) + dsigma(k, 2))*100. !Pa
                  phil(i, kc, jj) = phi(i, k, jj) - topo(i, jj)
                  !byl           phii(i,kc+1)=phi2(i,k)
               end do
            end do
         end do
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               prsi(i, lev + 1, jj) = ptop*100.     !Pa
               !byl       phii(i,1)=phi2(i,lev+1)
            end do
         end do
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               heat(i, jj) = hflux(i, jj)
               evap(i, jj) = qflux(i, jj)
               !      stress(i)=stress(i)
               rcl(i, jj) = 1.
            end do
         end do
         !
         !     if(jj.eq.jo)then
         !     print*,'before mixpbl uvtq-----------------'
         !     do k=lev,1,-1
         !     print*,k,u1(io,k),v1(io,k),t1(io,k),q1(io,k,1),q1(io,k,2)
         !     enddo
         !     print*,'before  phii,p2,phil,p-----------------'
         !     do k=lev,1,-1
         !     print*,k,phii(io,k),prsi(io,k),phil(io,k),prsl(io,k)
         !     enddo
         !     print*,'before  del prslk -----------------'
         !     do k=lev,1,-1
         !     print*,k,del(io,k),prslk(io,k)
         !     enddo
         !     print*,'sfc tg,sfcw,heat,evap,stress--='
         !     print*,tg(io),sfcw(io),heat(io),evap(io),stress(io)
         !     print*,'sfc rb,fm,fh--='
         !     print*,rb(io),fm(io),fh(io)
         !     endif
         !
         if (nmpbl .eq. 1) then
            if (myrank .eq. 0) &
               print *, ' warning !!!, nmpbl can not be 1, reassign nmpbl=2'
            nmpbl = 2
         end if
         if (nmpbl .eq. 2) then
            do jj = 1, jlistnum
               j = jlist1(jj)
               call mixpbl_n(nx, myim(jj), lev, ntrac, u1(1, 1, jj), v1(1, 1, jj), &
                             t1(1, 1, jj), q1(1, 1, 1, jj), &
                             pk2(1, lev, jj), rb(1, jj), fm(1, jj), fh(1, jj), &
                             tg(1, jj), heat(1, jj), evap(1, jj), stress(1, jj), sfcw(1, jj), &
                             kpbl(1, jj), prsi(1, 1, jj), del(1, 1, jj), &
                             prsl(1, 1, jj), prslk(1, 1, jj), &
                             phii(1, 1, jj), phil(1, 1, jj), rcl(1, jj), dt, &
                             hpbl(1, jj), ice(1, jj), g, cp, hltm, &
                             r, j)
            end do
         end if
         !
         !---
         ! new_sas pbl
         if (nmpbl .eq. 3) then
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     swh(i, kc, jj) = asl(i, k, jj)/86400.
                     hlw(i, kc, jj) = atl(i, k, jj)/86400.
                  end do
               end do
            end do
            do jj = 1, jlistnum
               j = jlist1(jj)
               call mixpbl_new(nx, myim(jj), lev, ntrac, ntcw, u1(1, 1, jj), &
                               v1(1, 1, jj), t1(1, 1, jj), &
                               q1(1, 1, 1, jj), swh(1, jj), hlw(1, jj), &
                               xmu(1, jj), pk2(1, lev, jj), &
                               rb(1, jj), fm(1, jj), fh(1, jj), tg(1, jj), heat(1, jj), &
                               evap(1, jj), stress(1, jj), sfcw(1, jj), kpbl(1, jj), prsi(1, 1, jj), &
                               del(1, 1, jj), prsl(1, 1, jj), prslk(1, 1, jj), &
                               phii(1, 1, jj), phil(1, 1, jj), &
                               rcl(1, jj), dt, hpbl(1, jj), ice(1, jj), g, &
                               cp, hltm, r, j)
            end do
            !
         end if

         !---
         ! YSU pbl scheme
         if (nmpbl .eq. 5) then
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     swh(i, kc, jj) = asl(i, k, jj)/86400.
                     hlw(i, kc, jj) = atl(i, k, jj)/86400.
                  end do
               end do
            end do
            do jj = 1, jlistnum
               j = jlist1(jj)
               call ysu2d(u1(1, 1, jj), v1(1, 1, jj), t1(1, 1, jj), q1(1, 1, 1, jj), &
                          prsl(1, 1, jj), prsi(1, 1, jj), prslk(1, 1, jj), &
                          pk2(1, lev, jj), nx, myim(jj), lev, ntrac, del(1, 1, jj), cp, &
                          g, r, hltm, phii(1, 1, jj), phil(1, 1, jj), psi(1, jj), &
                          z0(1, jj), stress(1, jj), hpbl(1, jj), kpbl(1, jj), &
                          fm(1, jj), fh(1, jj), islmsk(1, jj), heat(1, jj), evap(1, jj), &
                          sfcw(1, jj), rb(1, jj), dt, rcl(1, jj), u10(1, jj), v10(1, jj), &
                          swh(1, jj), hlw(1, jj), xmu(1, jj), j)
            end do
            !
         end if

         !
         ! NCEP GFS moninedmf
         if (nmpbl .eq. 4) then
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     swh(i, kc, jj) = asl(i, k, jj)/86400.
                     hlw(i, kc, jj) = atl(i, k, jj)/86400.
                  end do
               end do
            end do

            !$acc enter data copyin(swh, hlw, xmu, rb, z0rl, u10, v10, &
            !$acc&      fm, fh, tg, heat, evap, stress, sfcw, prsi, &
            !$acc&      del, prsl, prslk, phil) async(async_id)
            !$acc enter data copyin(u1, v1, t1, q1, kpbl, hpbl) async(async_id)
            !call cpu_time(time1)
            call moninedmf_gpu(nx, nxjp, lev, ntrac, ntcw, u1, &
                               v1, t1, q1, &
                               swh, hlw, xmu, pk2, &
                               rb, z0rl, u10, v10, fm, &
                               fh, tg, heat, evap, stress, &
                               sfcw, kpbl, prsi, &
                               del, prsl, &
                               prslk, phii, phil, dt, &
                               hpbl, sparsehandle)
            !call cpu_time(time2)
            !$acc exit data copyout(u1, v1, t1, q1, kpbl, hpbl) async(async_id)
            !$acc exit data delete(swh, hlw, xmu, rb, z0rl, u10, v10, &
            !$acc&     fm, fh, tg, heat, evap, stress, sfcw, prsi, &
            !$acc&     del, prsl, prslk, phil) async(async_id)
            !$acc wait(async_id)
            !write(*,*) 'aaaaa', time2-time1
            !
         end if
         !
         if (nmmiph .eq. 6) then  !WSM6
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     qt(i, lev*(ntcw - 1) + k, jj) = q1(i, kc, 2, jj)
                     qt(i, lev*(ntiw - 1) + k, jj) = q1(i, kc, 3, jj)
                     qt(i, lev*(ntoz - 1) + k, jj) = q1(i, kc, 4, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 8) then ! Thompson
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     qt(i, lev*(ntcw - 1) + k, jj) = q1(i, kc, 2, jj)
                     qt(i, lev*(ntiw - 1) + k, jj) = q1(i, kc, 3, jj)
                     qt(i, lev*(ntinc - 1) + k, jj) = q1(i, kc, 4, jj)
                     qt(i, lev*(ntoz - 1) + k, jj) = q1(i, kc, 5, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 18) then ! 2M Thompson
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     qt(i, lev*(ntcw - 1) + k, jj) = q1(i, kc, 2, jj)
                     qt(i, lev*(ntiw - 1) + k, jj) = q1(i, kc, 3, jj)
                     qt(i, lev*(ntinc - 1) + k, jj) = q1(i, kc, 4, jj)
                     qt(i, lev*(ntoz - 1) + k, jj) = q1(i, kc, 5, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 11 .or. nmmiph .eq. 12 .or. nmmiph .eq. 13) then ! GFDL MP
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     qt(i, lev*(ntcw - 1) + k, jj) = q1(i, kc, 2, jj)
                     qt(i, lev*(ntiw - 1) + k, jj) = q1(i, kc, 3, jj)
                     qt(i, lev*(ntrw - 1) + k, jj) = q1(i, kc, 4, jj)
                     qt(i, lev*(ntsw - 1) + k, jj) = q1(i, kc, 5, jj)
                     qt(i, lev*(ntgl - 1) + k, jj) = q1(i, kc, 6, jj)
                     qt(i, lev*(ntoz - 1) + k, jj) = q1(i, kc, 7, jj)
                  end do
               end do
            end do
         else if (nmmiph .eq. 15 .or. nmmiph .eq. 16) then ! Goddard MP
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     qt(i, lev*(ntcw - 1) + k, jj) = q1(i, kc, 2, jj)
                     qt(i, lev*(ntiw - 1) + k, jj) = q1(i, kc, 3, jj)
                     qt(i, lev*(ntrw - 1) + k, jj) = q1(i, kc, 4, jj)
                     qt(i, lev*(ntoz - 1) + k, jj) = q1(i, kc, 5, jj)
                  end do
               end do
            end do
         else
            do jj = 1, jlistnum
               do nc = 2, ntrac
                  do k = 1, lev
                     kc = lev - k + 1
                     do i = 1, myim(jj)
                        qt(i, lev*(nc - 1) + k, jj) = q1(i, kc, nc, jj)
                     end do
                  end do
               end do
            end do
         end if
         !
         do jj = 1, jlistnum
            do k = 1, lev
               !jh       do k=ktpbl,lev
               kc = lev - k + 1
               do i = 1, myim(jj)
                  ! time split
                  !            tt(i,k) = t1(i,kc)
                  !            ut(i,k) = u1(i,kc)
                  !            vt(i,k) = v1(i,kc)
                  ! no time split
                  dtdt(i, kc, jj) = (t1(i, kc, jj) - tt(i, k, jj))/dt
                  dudt(i, kc, jj) = (u1(i, kc, jj) - ut(i, k, jj))/dt
                  dvdt(i, kc, jj) = (v1(i, kc, jj) - vt(i, k, jj))/dt
                  dqdt(i, kc, jj) = (q1(i, kc, 1, jj) - qt(i, k, jj))/dt
               end do
            end do
         end do
         !

         !     call maxp ( nx,xkm(1,lev)  ,xkmx(1),imx(1) )
         !     call maxp ( nx,xkm(1,lev-1),xkmx(2),imx(2) )
         !     if ( idg .ne. 0 )  then
         !       do 260 k = 1, lev
         !       xkmd(k) = xkm(idg,k)
         ! 260   continue
         !     endif
         !
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               hflux(i, jj) = hflux(i, jj)*ro2(i, jj)*cp          ! transfer to W/m2
               qflux(i, jj) = qflux(i, jj)*ro2(i, jj)*hltm        ! transfer to W/m2
               tstar(i, jj) = -heat(i, jj)/ustar(i, jj)
               qstar(i, jj) = -evap(i, jj)/ustar(i, jj)
               runoff(i, jj) = (drain(i, jj) + runof(i, jj))*dth(jj) + runoff(i, jj) !wei add at 20231012
            end do
         end do
         deallocate (q1)
         istat = cusparseDestroy(sparsehandle)
!     if(jj.eq.jo)then
!     print*,'after mixpbl -----------------'
!     do k=lev,1,-1
!     print*,k,u1(io,k),v1(io,k),t1(io,k),q1(io,k,1),q1(io,k,2)
!     enddo
!     print*,'sfc kpbl,hpbl=',kpbl(io),hpbl(io)
!     print*,'hflux,qflux(w/m2)=',hflux(io),qflux(io)
!     endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     change tt and tg back to temp                                    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     do 300 k = 1, lev
!     do 300 i = 1, nx
!     tt(i,k) = tt(i,k)*pk(i,k)
! 300 continue
!soil do 320 i = 1, nx
!     tg(i) = tg(i)*pk2(i,lev)
! 320 continue
!soil
         return
      end
