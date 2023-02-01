      subroutine grdtpw ( mn,dt,hltm,tice,hice,czh,land,ocean,snow,ts &
                        , tgclim,tg,gwclim,gwr,gwet,ss,rs,hflux,qflux &
                        , snr,totalp )
!
!#######################################################################
!           subroutine description
!
! 1. function description
!
!    this subroutine predicts ground temperature for land points with
!    energy budget concept including newton cooling term, and ground
!    water for land and snow points with water mass budget concept.
!
! 2. parameter specification
!
!    mn : dimension of horizontal-direction
!    dt : time step for tg updates                         (s)     i
!  hltm : water vapor latent heat constant                 (j/kg)  i
!  tice : temperature of melting point               (mn)  (k)     i
!  hice : latent heat of snow melting                      (j/kg)  i
!   czh : ground thermal capacity                          (j/k m2)i
! land  : logic for bare soil land                   (mn)          i
! ocean : logic for open water (=ture for open water)(mn)          i
! snow  : logic for snow covered area                (mn)          i
!   ts  : surface air temperature                    (mn)   (k)    i
! tgclim: ground real temp. climatolagy              (mn)  (k)     i
!    tg : ground real temp. (sea or land)            (mn)  (k)     i/o
! gwclim: soil water climatology                     (mn)  (mm)    i
!   gwr : soil water for land and ice points         (mn)  (mm)    i/o
!  gwet : ground wetness ( = gwr/20.0 )              (mn)          i/o
!    ss : net short wave flux down to ground surface (mn)  (wat/m2)i
!    rs : net long  wave flux up from ground surface (mn)  (wat/m2)i
! hflux : upward surface heat flux                   (mn)  (wat/m2)i
! qflux : upward surface moisture flux               (mn)  (wat/m2)i
!  snr  : snow depth in water state                  (mn)   (mm)   i/o
! totalp: total precipitation during dt              (mn)   (mm)   i
!
! 4. local veriable
!
!  betatg : newton cooling restorine rate (=1/(100 hr))     (1/s)
!     fhs : net budget of total flux                        (wat/m2)
!
! 5. calling modules
!
!     pbltke
!
! 6. usage
!
!     call grdtpw ( mn,dt,hltm,tice,hice,czh,land,ocean,snow,ts
!    1            , tgclim,tg,gwclim,gwr,gwet,ss,rs,hflux,qflux
!    2            , snr,totalp )
!
! 7. modules called
!
! 8. limitation
!
! 9. date
!
!    created       jan. 1992
!    modified      may  1993
!
! 10. author
!
!     f. wang,  c.s. liou
!     modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#######################################################################
      implicit  none
!
! input & output variable
      integer   mn
      real      dt,hltm,tice,hice
      real      czh(mn),ts(mn),tgclim(mn),                            &
                tg(mn),gwclim(mn),gwr(mn),gwet(mn),ss(mn),rs(mn),     &
                hflux(mn),qflux(mn),snr(mn),totalp(mn)

      logical   land(mn),ocean(mn),snow(mn)
!-----------------------------------------------------------------------
!
      integer   ij
      real      belnd,beice,alp,deltg,betg,fhs,tmp,fgtg
      real      dtywtr,gwrcc,rainr,water,cbhice,snomlt


      belnd= 6.2832/(3600.*100.)
      beice= 6.2832/(3600.*360.)
      alp  = 0.5
      deltg= 5.0
!
      do 100 ij = 1, mn
        if ( .not. ocean(ij))  then
          if ( land(ij) ) then
             betg = belnd
          else
             betg = beice
          endif
          fhs = ss(ij) - rs(ij) - qflux(ij) - hflux(ij)
          tmp = (tg(ij)+ dt*(fhs/czh(ij)+betg*tgclim(ij)))/(1.+betg*dt)
!
!         apply time filter when fhs is large (kalnay & kanamitsu,1988)
!         criterion: alpha*(p-1) >= 1.,  p=3.5
!
          fgtg = dt*fhs/(czh(ij)*deltg)
          if ( abs(fgtg) .ge. 0.4 )  then
            tg(ij) = (1.0-alp)*tmp + alp*tg(ij)
          else
            tg(ij) = tmp
          endif
!
! snr: mm
! totalp: mm
! dtywtr: kg/m3
! snomlt: mm
! cbhice: (j/m2k)/(j/kg)=kg/m2k
!
           dtywtr = 1000.0
           gwrcc  = 20.0
!
           if (land(ij).or.snow(ij))  then
             if (ts(ij).lt.tice)  then
               snr(ij) = snr(ij) + totalp(ij)
               rainr   = 0.0
             else
               rainr   = totalp(ij)
             endif
             if (tg(ij).lt.tice)  then
               water = rainr
             else
               cbhice  = czh(ij)/hice
               snomlt  = min(cbhice*(tg(ij)-tice)/dtywtr*1000.,snr(ij))
               tg(ij)  = tg(ij) - snomlt/1000.*dtywtr/cbhice
               snr(ij) = snr(ij) - snomlt
               water   = rainr + snomlt
             endif
             gwr(ij) = (gwr(ij) + dt*(water/dt- qflux(ij)/hltm + betg  &
                    * gwclim(ij))) / (1.0 + dt*(water/dt/gwrcc+betg))
             gwr(ij) = max( 0.01, min( gwrcc, gwr(ij) ) )
             gwet(ij)= gwr(ij)/gwrcc
           else
             tg(ij) = min( tg(ij),tice )
           endif
        endif
 100  continue
!
      return
      end
