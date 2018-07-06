 
      SUBROUTINE SFC_OCEAN(imj,IM,KM,PS,U1,V1,T1,Q1,  &
                        TSKIN,QSURF,DM,GFLUX,CM,CH,   &
                        RCL,PRSL1,PRSLKI,SLIMSK,      &
!                       CMM,CHH, 
                        EVAP,HFLX,EP,DDVEL,flag_iter) 
! 
      USE MACHINE , ONLY : kind_phys 
!     USE FUNCPHYS, ONLY : fpvs 
      USE PHYSCONS, HVAP => con_HVAP                                  ,&
                    CP => con_CP, HFUS => con_HFUS, JCAL => con_JCAL  ,&
                    EPS => con_eps, EPSM1 => con_epsm1                ,&
                    RVRDM1 => con_FVirt, RD => con_RD 
      implicit none 
! 
!     include 'constant.h' 
! 
      integer              IM, km,imj 
! 
      real(kind=kind_phys), parameter :: cpinv=1.0/cp, HVAPI=1.0/HVAP 
      real(kind=kind_phys) PS(IM),       U1(IM),      V1(IM),    &
                           T1(IM),       Q1(IM),                 &
                           TSKIN(IM),    QSURF(IM),   DM(IM),    &
                           GFLUX(IM),                            &
                           CM(IM),       CH(IM),      RCL(IM),   &
                           PRSL1(IM),    PRSLKI(IM),  SLIMSK(IM),& 
                           EVAP(IM),     HFLX(IM),               &
                           EP(IM),       DDVEL(IM),              &
                           CMM(IM),      CHH(IM) 
 
      logical              flag_iter(im), FLAG(IM) 
 
! 
!     Locals 
! 
      integer              k,i 
! 
      real(kind=kind_phys) PSURF(IM),   PS1(IM),     Q0(IM),   &
                           QSS(IM),     RCH(IM),     RHO(IM),  &
                           THETA1(IM),  TV1(IM),     XRCL(IM), &
                           WIND(IM) 
! 
      real(kind=kind_phys) elocp,  tem 
!c 
      PARAMETER (ELOCP=HVAP/CP) 
 
! 
! FLAG for open water 
! 
      DO I = 1, IMj 
         FLAG(I) = SLIMSK(I).EQ. 0. .AND. flag_iter(i) 
      ENDDO 
 
! 
!  INITIALIZE VARIABLES. ALL UNITS ARE SUPPOSEDLY M.K.S. UNLESS SPECIFIED 
!  PSURF IS IN PASCALS 
!  WIND IS WIND SPEED, THETA1 IS ADIABATIC SURFACE TEMP FROM LEVEL 1 
!  RHO IS DENSITY, QSS IS SAT. HUM. AT SURFACE 
! 
!! 
      DO I=1,IMj 
        IF(FLAG(I)) THEN 
        XRCL(I)  = SQRT(RCL(I)) 
        PSURF(I) = 1000. * PS(I) 
        PS1(I)   = 1000. * PRSL1(I) 
        WIND(I) = XRCL(I) * SQRT(U1(I) * U1(I) + V1(I) * V1(I))  &
                    + MAX(0.0, MIN(DDVEL(I), 30.0)) 
        WIND(I) = MAX(WIND(I),1.) 
        Q0(I) = MAX(Q1(I),1.E-8) 
        THETA1(I) = T1(I) * PRSLKI(I) 
        TV1(I) = T1(I) * (1. + RVRDM1 * Q0(I)) 
        RHO(I) = PS1(I) / (RD * TV1(I)) 
!       qss(i) = fpvs(tskin(i)) 
!       QSS(I) = EPS * QSS(I) / (PSURF(I) + EPSM1 * QSS(I)) 
!
        call qsatq(1, tskin(i), psurf(i)*0.01, qss(i))
!
        ENDIF 
      ENDDO 
 
      DO I=1,IMj 
       IF(FLAG(I)) THEN 
        EVAP(I) = 0. 
        HFLX(I) = 0. 
        GFLUX(I) = 0. 
        EP(I) = 0. 
       ENDIF 
      ENDDO 
 
! 
!  RCP = RHO CP CH V 
! 
      DO I = 1, IMj 
       IF(FLAG(I)) THEN 
        RCH(I) = RHO(I) * CP * CH(I) * WIND(I) 
!wei added 10/24/2006 
        CMM(I)=CM(I)* WIND(I) 
        CHH(I)=RHO(I)*CH(I)* WIND(I) 
       ENDIF 
      ENDDO 
! 
!  SENSIBLE AND LATENT HEAT FLUX OVER OPEN WATER 
! 
      DO I = 1, IMj 
        IF(FLAG(I)) THEN 
          EVAP(I) = ELOCP * RCH(I) * (QSS(I) - Q0(I)) 
          DM(I) = 1. 
          QSURF(I) = QSS(I) 
        ENDIF 
      ENDDO 
!! 
! 
!  CALCULATE SENSIBLE HEAT FLUX 
! 
      DO I = 1, IMj 
        IF(FLAG(I)) THEN 
        HFLX(I) = RCH(I) * (TSKIN(I) - THETA1(I)) 
        ENDIF 
      ENDDO 
! 
!  THE REST OF THE OUTPUT   <----  Note: Redundant calculation  
! 
!*      DO I = 1, IM 
!*        IF(SLIMSK(I).EQ.0.) THEN 
!*        QSURF(I) = Q1(I) + EVAP(I) / (ELOCP * RCH(I)) 
!*        DM(I) = 1. 
!*      ENDDO 
! 
      do i=1,imj 
        IF(FLAG(I)) THEN 
        tem     = 1.0 / rho(i) 
        hflx(i) = hflx(i) * tem * cpinv 
        evap(i) = evap(i) * tem * hvapi 
        ENDIF 
      enddo 
! 
      RETURN 
      END 
