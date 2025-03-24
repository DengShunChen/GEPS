!-----------------------------------
      subroutine sfc_ocean                                              &
!...................................
!  ---  inputs: &
           ( imj,im, ps, u1, v1, t1, q1, tskin, cm, ch,                 &  
             prsl1, prslki, islimsk, ddvel, flag_iter,                  &
!  ---  outputs:
!    &       qsurf, cmm, chh, gflux, evap, hflx, ep                     &  
#ifdef TIMCOMCPL
             qsurf, gflux, evap, hflx, ep, ssu, ssv)
#else
             qsurf, gflux, evap, hflx, ep)
#endif

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_ocean                                                     !
!       inputs:                                                         !
!          ( im, ps, u1, v1, t1, q1, tskin, cm, ch,                     !
!            prsl1, prslki, islimsk, ddvel, flag_iter,                  !
!       outputs:                                                        !
!            qsurf, cmm, chh, gflux, evap, hflx, ep )                   !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!         2005  -- created from the original progtm to account for      !
!                  ocean only                                           !
!    oct  2006  -- h. wei      added cmm and chh to the output          !
!    apr  2009  -- y.-t. hou   modified to match the modified gbphys.f  !
!                  reformatted the code and added program documentation !
!    sep  2009  -- s. moorthi removed rcl and made pa as pressure unit  !
!                  and furthur reformatted the code                     !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     ps       - real, surface pressure                            im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     islimsk  - integer, sea/land/ice mask (=0/1/2)               im   !
!     ddvel    - real, wind enhancement due to convection (m/s)    im   !
!     flag_iter- logical,                                          im   !
!                                                                       !
!  outputs:                                                             !
!     qsurf    - real, specific humidity at sfc                    im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     gflux    - real, ground heat flux (zero for ocean)           im   !
!     evap     - real, evaporation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
!     use funcphys, only : fpvs
      use physcons, only : cp => con_cp, rd => con_rd, eps => con_eps,    &
                           epsm1 => con_epsm1, hvap => con_hvap,          &
                           rvrdm1 => con_fvirt
      use const,    only : RTYPE
!
      implicit none
!
!  ---  constant parameters:
      real (kind=kind_phys), parameter :: cpinv  = 1.0/cp                 &
      ,                                   hvapi  = 1.0/hvap               &
      ,                                   elocp  = hvap/cp

!  ---  inputs:
      integer, intent(in) :: im, imj

      real (kind=RTYPE), dimension(im), intent(in) :: u1, v1, t1, q1
      real (kind=kind_phys), dimension(im), intent(in) :: ps,             &
#ifdef TIMCOMCPL
              tskin, cm, ch, prsl1, prslki, ddvel, ssu, ssv
#else
              tskin, cm, ch, prsl1, prslki, ddvel
#endif
      integer, dimension(im), intent(in):: islimsk

      logical, intent(in) :: flag_iter(im)

!  ---  outputs:
!     real (kind=kind_phys), dimension(im), intent(out) :: qsurf,         &
      real (kind=kind_phys), dimension(im) :: qsurf,         &
             cmm, chh, gflux, evap, hflx, ep, cpl_u1, cpl_v1

!  ---  locals:

      real (kind=kind_phys) :: q0, qss, rch, rho, wind, tem, fpvs

      integer :: i

      logical :: flag(im)
#ifdef TIMCOMCPL
      cpl_u1=0.
      cpl_v1=0.
      cpl_u1=u1
      cpl_v1=v1
#endif
!
!===> ...  begin here
!
!  --- ...  flag for open water
      do i = 1, imj
        flag(i) = ( islimsk(i) == 0 .and. flag_iter(i) )

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specified
!           ps is in pascals, wind is wind speed,
!           rho is density, qss is sat. hum. at surface

        if ( flag(i) ) then

#ifdef TIMCOMCPL
          wind     = max(sqrt((cpl_u1(i)-ssu(i))**2 + (cpl_v1(i)-ssv(i))**2)          &
                       + max( 0.0, min( ddvel(i), 30.0 ) ), 1.0)
#else
          wind     = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))                  &
                       + max( 0.0, min( ddvel(i), 30.0 ) ), 1.0)
#endif          
          q0       = max( q1(i), 1.0e-8 )
          rho      = prsl1(i) / (rd*t1(i)*(1.0 + rvrdm1*q0))

          qss      = fpvs( tskin(i) )
          qss      = eps*qss / (ps(i) + epsm1*qss)

          evap(i)  = 0.0
          hflx(i)  = 0.0
          ep(i)    = 0.0
          gflux(i) = 0.0

!  --- ...    rcp  = rho cp ch v

          rch      = rho * cp * ch(i) * wind
          cmm(i)   = cm(i) * wind
          chh(i)   = rho * ch(i) * wind

!  --- ...  sensible and latent heat flux over open water

          hflx(i)  = rch * (tskin(i) - t1(i) * prslki(i))

          evap(i)  = elocp*rch * (qss - q0)
          qsurf(i) = qss

          tem      = 1.0 / rho
          hflx(i)  = hflx(i) * tem * cpinv
          evap(i)  = evap(i) * tem * hvapi
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean
!-----------------------------------
