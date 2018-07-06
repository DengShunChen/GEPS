      subroutine grdcon ( mn,tg,gwet,czh,land,ocean,snow,ice,tice &
                        , as,snr )
!
!***********************************************************************
!            subroutine description
!
! 1. function description
!
!    this subroutine deterime ground condition, then depending on
!    the condition to calculate thermal capacity (czh) and diffusion
!    albedo (as).
!
! 2. parameter specification
!
!    mn : dimension of horizontal-direction                          i
!    tg : ground temperature                            (mn)   (k)   i/o
!  gwet : ground wetness(water volume/soil space volume (mn)         i/o
!   czh : soil thermal capacity                         (mn)  (j/m2k)o
!  land : logic for bare soil land                      (mn)         i
!  ocean: logic for open water   (=true for open water) (mn)         i
!  snow : logic for snow covered area                   (mn)         o
!   ice : logic for ice covered area                    (mn)         i
!  tice : temperature of melting point                  (mn)   (k)   i
!    as : diffusion albedo                              (mn)         o
!   snr : snow depth in water state                     (mn)   (mm)  i/o
!
! 3. local variable
!
! 4. calling modules
!
!    pbltke
!
! 5. usage
!
!    call grdcon ( mn,tg,gwet,czh,land,ocean,snow,ice,tice,as,snr )
!
! 6. modules called
!
! 7. limitation
!
! 8. date
!
!    created    feb. 1992
!    modified   may  1993
!
! 9. author
!
!    f. wang,  c.s. liou
!    modify to f90 by C-H Lee and sort by River Chen in 2015
!
!10. reference
!#######################################################################
      implicit  none
!
!  input & ouput variable
      integer   mn
      real      tice
      real      tg(mn),gwet(mn),czh(mn),as(mn),snr(mn)
      logical   land(mn),ocean(mn),snow(mn),ice(mn)
!
!----------------------------------------------------------------------
      integer   i
      real      pi,frday,wet
!
      pi    = 3.141592654
      frday = 2.0*pi/86400.0
!
!----------------------------------------------------------------------
! determine ground condition by
!  from data base:
!    ocean            for open water
!    land             for bare land
!    ice              for ice (seaice & land ice)
! from snow depth:
!    snow             check if (snr .gt. 0.) on land and ice
!                     when snr >0.   land(or ice) become snow
!
!----------------------------------------------------------------------
!
      do 400 i = 1, mn
        snow(i) = land(i) .and. (snr(i).gt.0.0)
!xx     snow(i) = (land(i).or.ice(i)).and.(snr(i).gt.0.)
 400  continue
!
! ---------------------------------------------------------------------
!   determine czh and as ,depending on ground condiction.
!        czh is soil thermal capacity (j/m2/k)
!        as is albedo
!
!              czh : rho(kg/m3) *cp(j/kgk) *h(m)]
!              h   : effective thickness (=sqrt(2ks/w))
!              ks  : thermal diffusivity of soil
!              w   : frequency (=2pi/86400)
!----------------------------------------------------------------------
!
      do 500 i = 1, mn
        if (snow(i)) then
          as(i)  = 0.7
          czh(i) = 4.2e4*2.3
        else if (ice(i)) then
          czh(i) = 4.2e4*5.1
        else if (land(i)) then
          wet = gwet(i)
          if (wet.le.1.e-10)  wet = 0.0
!xxx      if (tg(i).le.tice)  wet = min( wet,0.01 )
          czh(i)= 4.2e4*sqrt( (0.387+0.15*wet)*(1.+wet)*2.e-3/frday )
        endif
 500  continue
!
      return
      end
