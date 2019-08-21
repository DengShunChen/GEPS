      SUBROUTINE SFC_DIAG(imj,IM,KM,PS,U1,V1,T1,Q1, &
                        TSKIN,QSURF,                &
                        U10M,V10M,T2M,Q2M,          &
                        RCL,PRSLKI,SLIMSK,          &
                        EVAP,FM,FH,FM10,FH2,FH10,rh2,rh10) 
! 
      USE MACHINE , ONLY : kind_phys 
!     USE FUNCPHYS, ONLY : fpvs 
      USE PHYSCONS, grav => con_g, SBC => con_sbc, HVAP => con_HVAP  , &
                    CP => con_CP, HFUS => con_HFUS, JCAL => con_JCAL , &
                    EPS => con_eps, EPSM1 => con_epsm1               , &
                    RVRDM1 => con_FVirt, RD => con_RD 
      implicit none 
! 
!     include 'constant.h' 
! 
      integer              IM, km ,imj,j
! 
      real(kind=kind_phys) PS(IM),       U1(IM),      V1(IM),          &
                           T1(IM),       Q1(IM),                       &
                           TSKIN(IM),    QSURF(IM),                    &
                           F10M(IM),     U10M(IM),                     &
                           V10M(IM),     T2M(IM),     T10M(IM),        &
                           Q2M(IM),      Q10M(IM),                     &
                           RCL(IM),      PRSL1(IM),   PRSLKI(IM),      &
                           SLIMSK(IM),   EVAP(IM),                     &
                           FM(IM),       FH(IM),                       &
                           FM10(IM),     FH2(IM),   FH10(IM),          &
                           rh2(IM),rh10(IM),p2(IMj),p10(IMj),rosfc(im)
! 
!     Locals 
! 
      real (kind=kind_phys), parameter :: qmin=1.0e-8 
      integer              k,i 
! 
      real(kind=kind_phys)                        &
                           PSURF(IM),   QSS(IM),  &
                           THETA1(IM),  XRCL(IM)
! 
      real(kind=kind_phys) g,    sig2k 
! 
!c 
      PARAMETER (G=grav) 
! 
      LOGICAL FLAG(IM), FLAGSNW(IM) 
      real(kind=kind_phys) KT1(IM),       KT2(IM),      KTSOIL,    &
                           ET(IM,KM),                              &
                           STSOIL(IM,KM), AI(IM,KM),    BI(IM,KM), &
                           CI(IM,KM),     RHSTC(IM,KM) 
! 
! 
!     ESTIMATE SIGMA ** K AT 2 M 
! 
      SIG2K = 1. - 4. * G * 2. / (CP * 280.) 
! 
!  INITIALIZE VARIABLES. ALL UNITS ARE SUPPOSEDLY M.K.S. UNLESS SPECIFIE 
!  PSURF IS IN PASCALS 
!  THETA1 IS ADIABATIC SURFACE TEMP FROM LEVEL 1 
! 
!! 
      DO I=1,IMj 
        XRCL(I)  = SQRT(RCL(I)) 
        PSURF(I) = 1000. * PS(I) 
        THETA1(I) = T1(I) * PRSLKI(I) 
!
        rosfc(i)=psurf(i)/(rd*tskin(i))
        p2(i)=psurf(i)-rosfc(i)*g*2.
        p10(i)=psurf(i)-rosfc(i)*g*10.
!
      ENDDO 
!! 
! 
      DO I = 1, IMj 
        F10M(I) = FM10(I) / FM(I) 
        F10M(I) = min(F10M(I),1.) 
        U10M(I) = F10M(I) * XRCL(I) * U1(I) 
        V10M(I) = F10M(I) * XRCL(I) * V1(I) 
         T2M(I) = TSKIN(I) * (1. - FH2(I) / FH(I)) &
                + THETA1(I) * FH2(I) / FH(I) 
         T2M(I) = T2M(I) * SIG2K 
        T10M(I) = TSKIN(I) * (1. - FH10(I) / FH(I)) &
                + THETA1(I) * FH10(I) / FH(I) 
        T10M(I) = T10M(I) * SIG2K 
!        Q2M(I) = QSURF(I) * (1. - FH2(I) / FH(I))
!     &         + Q1(I) * FH2(I) / FH(I) 
!       T2M(I) = T1 
!       Q2M(I) = Q1 
        IF(EVAP(I).GE.0.) THEN 
! 
!  IN CASE OF EVAPORATION, USE THE INFERRED QSURF TO DEDUCE Q2M 
! 
          Q2M(I) = QSURF(I) * (1. - FH2(I) / FH(I)) &
               + max(qmin,Q1(I)) * FH2(I) / FH(I)      !  Moorthi 
         Q10M(I) = QSURF(I) * (1. - FH10(I) / FH(I)) &
               + max(qmin,Q1(I)) * FH10(I) / FH(I)      !  Moorthi 
!!   &         + Q1(I) * FH2(I) / FH(I) 
        ELSE 
! 
!  FOR DEW FORMATION SITUATION, USE SATURATED Q AT TSKIN 
! 
!jfe      QSS(I) = 1000. * FPVS(TSKIN(I)) 
!         qss(I) = fpvs(tskin(I)) 
!         QSS(I) = EPS * QSS(I) / (PSURF(I) + EPSM1 * QSS(I)) 
!
          call qsatq(1, tskin(i), psurf(i)*0.01, qss(i))
!
          Q2M(I) = QSS(I) * (1. - FH2(I) / FH(I))  &
               + max(qmin,Q1(I)) * FH2(I) / FH(I)      ! Moorthi 
!!   &         + Q1(I) * FH2(I) / FH(I) 
         Q10M(I) = QSS(I) * (1. - FH10(I) / FH(I))  &
               + max(qmin,Q1(I)) * FH10(I) / FH(I)      ! Moorthi 
        ENDIF 
!       QSS(I) = fpvs(t2m(I)) 
!       QSS(I) = EPS * QSS(I) / (PSURF(I) + EPSM1 * QSS(I)) 
!
          call qsatq(1, t2m(i), psurf(i)*0.01, qss(i))
!
        Q2M(I) = MIN(Q2M(I),QSS(I)) 
          call qsatq(1, t10m(i), psurf(i)*0.01, qss(i))
        Q10M(I) = MIN(Q10M(I),QSS(I)) 
      ENDDO 
!    
        call qsatq(imj,t2m,p2*0.01,qss)
!
        do i=1,imj
        rh2(i)=max(min(q2m(i)/qss(i),1.),0.)
        enddo
!
        call qsatq(imj,t10m,p10*0.01,qss)
!
        do i=1,imj
        rh10(i)=max(min(q10m(i)/qss(i),1.),0.)
        enddo
      RETURN 
      END 
