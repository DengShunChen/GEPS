      SUBROUTINE SFC_DIFF(imj,IM,PS,U1,V1,T1,Q1,      &
                          TSKIN,Z0RL,CM,CH,RB,        &
                          RCL,PRSL1,PRSLKI,SLIMSK,    &
                          STRESS,FM,FH,               &
!lu_q2m_iter [-1L/+2L]: add tsurf, flag_iter 
!*   &                    USTAR,WIND,DDVEL,FM10,FH2) 
                          USTAR,WIND,DDVEL,FM10,FH2,  &
                          FH10,SIGMAF,VEGTYPE,SHDMAX, &
                          tsurf,flag_iter) 
! 
      USE MACHINE , ONLY : kind_phys 
!     USE FUNCPHYS, ONLY : fpvs     
      USE PHYSCONS, grav => con_g, SBC => con_sbc     ,&
                    CP => con_CP, HFUS => con_HFUS    ,&
                    RVRDM1 => con_FVirt, RD => con_RD ,&
                    EPS => con_eps, EPSM1 => con_epsm1 
 
      implicit none 
! 
!     include 'constant.h' 
! 
      integer              IM, km, ipr,j,io,jo,imj
! 
      real(kind=kind_phys) PS(IM),       U1(IM),      V1(IM),     &
                           T1(IM),       Q1(IM),                  &
                           TSKIN(IM),    Z0RL(IM),                &
                           CM(IM),       CH(IM),      RB(IM),     &
                           RCL(IM),      PRSL1(IM),   PRSLKI(IM), &
                           SLIMSK(IM),   STRESS(IM),              &
                           FM(IM),       FH(IM),      USTAR(IM),  &
                           WIND(IM),     DDVEL(IM),               &
                           FM10(IM),FH2(IM),FH10(IM),SIGMAF(IM),  &
                           SHDMAX(IM) 
      integer VEGTYPE(IM)
 
!lu_q2m_iter [+1L]: add flag_iter 
      logical              flag_iter(im) 
 
! 
!     Locals 
! 
      integer              k,i 
! 
      real(kind=kind_phys) DTV(IM),     HL1(IM),     HL12(IM), &
                           HLINF(IM),   PH(IM),                &
                           PH2(IM),PH10(IM),  PM(IM),PM10(IM), &
                           PSURF(IM),   Q0(IM),      RAT(IM),  &
                           THETA1(IM),  THV1(IM),              &
                           TSURF(IM),   TV1(IM),               &
                           TVS(IM),     XRCL(IM),              &
                           Z0(IM),      Z0MAX(IM),   Z1(IM),   &
                           ZTMAX(IM),   PS1(IM),     QS1(IM) 
 
! 
      real(kind=kind_phys) a0,    a0p,    a1,    a1p,   aa,   aa0,    &
                           aa1,   adtv,   alpha, arnu,  b1,   b1p,    &
                           b2,    b2p,    bb,    bb0,   bb1,  bb2,    &
                           ca,    cc,     cc1,   cc2,   charnock,     &
                           cq,    fms,    fhs,   g,     hl0,  hl0inf, &
                           hl110, hlt,    hltinf,OLINF,               &
                           restar, rnu,   vis,   DELP
! 
!  
      PARAMETER (CHARNOCK=.014,CA=.4)!C CA IS THE VON KARMAN CONSTANT 
!     PARAMETER (CHARNOCK=.018,CA=.4)!C CA IS THE VON KARMAN CONSTANT 
      PARAMETER (G=grav) 
      PARAMETER (ALPHA=5.,A0=-3.975,A1=12.32,B1=-7.755,B2=6.041) 
      PARAMETER (A0P=-7.941,A1P=24.75,B1P=-8.705,B2P=7.899,VIS=1.4E-5) 
      PARAMETER (AA1=-1.076,BB1=.7045,CC1=-.05808) 
      PARAMETER (BB2=-.1954,CC2=.009999) 
!     PARAMETER (RNU=1.51E-5,ARNU=.135*RNU) 
! MBEK -- TOGA CORE flux algorithm
      PARAMETER (RNU=1.51E-5,ARNU=.11*RNU) 
! 
!  INITIALIZE VARIABLES. ALL UNITS ARE SUPPOSEDLY M.K.S. UNLESS SPECIFIED 
!  PSURF IS IN PASCALS 
!  WIND IS WIND SPEED, THETA1 IS ADIABATIC SURFACE TEMP FROM LEVEL 1 
!  SURFACE ROUGHNESS LENGTH IS CONVERTED TO M FROM CM 
! 
      DO I=1,IMj 
        if(flag_iter(i)) then  
        XRCL(I)  = SQRT(RCL(I)) 
        PSURF(I) = 1000. * PS(I) 
!**     TSURF(I) = TSKIN(I)                 !! <---- Clu_q2m_iter [-1L] 
        PS1(I)   = 1000. * PRSL1(I) 
        WIND(I) = XRCL(I) * SQRT(U1(I) * U1(I) + V1(I) * V1(I))  &
                    + MAX(0.0, MIN(DDVEL(I), 30.0)) 
        WIND(I) = MAX(WIND(I),1.) 
        Q0(I) = MAX(Q1(I),1.E-8) 
        THETA1(I) = T1(I) * PRSLKI(I) 
        TV1(I) = T1(I) * (1. + RVRDM1 * Q0(I)) 
        THV1(I) = THETA1(I) * (1. + RVRDM1 * Q0(I)) 
!lu_q2m_iter[-1L/+2L]: TVS is computed from avg(tsurf,tskin) 
!**     TVS(I) = TSURF(I) * (1. + RVRDM1 * Q0(I)) 
        TVS(I) = 0.5 * (TSURF(I)+TSKIN(I)) *   &
                 (1. + RVRDM1 * Q0(I)) 
!       qs1(i) = fpvs(t1(i)) 
!       QS1(I) = EPS * QS1(I) / (PS1(I) + EPSM1 * QS1(I)) 
        call qsatq(1,t1(i),ps1(i)*0.01,qs1(i))
        QS1(I) = MAX(QS1(I), 1.E-8) 
        Q0(I) = min(QS1(I),Q0(I)) 
 
        Z0(I) = .01 * Z0RL(i) 
!byl        Z1(I) = -RD * TV1(I) * LOG(PS1(I)/PSURF(I)) / G 
        DELP  = (PSURF(I) - PS1(I)) / G
        Z1(I) = DELP * RD * TV1(I) / PS1(I)
        endif 
      ENDDO 
!! 
! 
!  COMPUTE STABILITY DEPENDENT EXCHANGE COEFFICIENTS 
! 
!  THIS PORTION OF THE CODE IS PRESENTLY SUPPRESSED 
!
 
      DO I=1,IMj 
       if(flag_iter(i)) then  
        IF(SLIMSK(I).EQ.0.) THEN 
          USTAR(I) = SQRT(G * Z0(I) / CHARNOCK) 
        ENDIF 
! 
! phon 2012/2/23
!      ustar just follow last step; because it is saved during integral.
!                 ustar is updated in the end of this sub 
! 
!  COMPUTE STABILITY INDICES (RB AND HLINF) 
! 
        Z0MAX(I) = MIN(Z0(I),1. * Z1(I)) 
!
! **  test XUBN's new z0 for thermal roughness
!
        IF(SLIMSK(I).NE.0.) THEN 
!

!CWB2015, change all alog to dlog
        Z0MAX(I) = exp( ((1.-SHDMAX(I))**2)*dlog(0.01)+    & 
                 (1.-((1.-SHDMAX(I))**2))*dlog(Z0MAX(I)) )
!
          if(VEGTYPE(I).eq.7)then
          Z0MAX(I) = exp( ((1.-SHDMAX(I))**2)*dlog(0.01)+  &
                (1.-((1.-SHDMAX(I))**2))*dlog(0.07) )
          endif
!
          if(VEGTYPE(I).eq.8)then
          Z0MAX(I) = exp( ((1.-SHDMAX(I))**2)*dlog(0.01)+  &
                (1.-((1.-SHDMAX(I))**2))*dlog(0.05) )
          endif
!
          if(VEGTYPE(I).eq.9)then
          Z0MAX(I) = exp( ((1.-SHDMAX(I))**2)*dlog(0.01)+  &
                (1.-((1.-SHDMAX(I))**2))*dlog(0.01) )
          endif
!
          if(VEGTYPE(I).eq.11)then
          Z0MAX(I) = exp( ((1.-SHDMAX(I))**2)*dlog(0.01)+  &
                (1.-((1.-SHDMAX(I))**2))*dlog(0.01) )
          endif
       ENDIF
!
        ZTMAX(I)= Z0MAX(I)*exp( - ((1.-SIGMAF(I))**2)      &
                  *0.8*CA*sqrt(USTAR(I)*0.01/(1.5e-05)))
!
!***       ZTMAX(I) = Z0MAX(I) 
!
        IF(SLIMSK(I).EQ.0.) THEN 
          RESTAR = USTAR(I) * Z0MAX(I) / VIS 
          RESTAR = MAX(RESTAR,.000001) 
!  Rat taken from Zeng, Zhao and Dickinson 1997 
          RAT(I) = 2.67 * restar ** .25 - 2.57 
          RAT(I) = min(RAT(I),7.) 
          ZTMAX(I) = Z0MAX(I) * EXP(-RAT(I)) 
        ENDIF 
! *Z0 over sea ** refer to ECMWF 
!     
!       IF(SLIMSK(I).EQ. 0.)THEN
!       Z0MAX(I)= 0.11*RNU/USTAR(I) + 0.018*USTAR(I)*USTAR(I)/g
!       ZTMAX(I)= 0.4*RNU/USTAR(I) 
!       ZQMAX(I)= 0.62*RNU/USTAR(I)
!       ENDIF
!**
       endif 
      ENDDO 
!
!
!      if(j.eq.jo)then
!      print*,'flag_iter=',flag_iter(io)
!      print*,'ustar=',ustar(io)
!      print*,'z0,z1=',z0(io),z1(io)
!      print*, 'z0max=',z0max(io)
!      print*, 'ztmax=',ztmax(io)
!      print*, 'ztmax_ec=',0.4*rnu/ustar(io)
!      endif
!
!##DG  IF(LAT.EQ.LATD) THEN 
!##DG    PRINT *, ' z0max, ztmax, restar, RAT(I) =',  
!##DG &   z0max, ztmax, restar, RAT(I) 
!##DG  ENDIF 
      DO I = 1, IMj 
       if(flag_iter(i)) then  
        DTV(I) = THV1(I) - TVS(I) 
        ADTV = ABS(DTV(I)) 
        ADTV = MAX(ADTV,.001) 
        DTV(I) = SIGN(1.,DTV(I)) * ADTV 
        RB(I) = G * DTV(I) * Z1(I) / (.5 * (THV1(I) + TVS(I)) & 
                * WIND(I) * WIND(I)) 
        RB(I) = MAX(RB(I),-5000.) 
        FM(I) = LOG((Z0MAX(I)+Z1(I)) / Z0MAX(I)) 
        FH(I) = LOG((ZTMAX(I)+Z1(I)) / ZTMAX(I)) 
        HLINF(I) = RB(I) * FM(I) * FM(I) / FH(I) 
        FM10(I) = LOG((Z0MAX(I)+10.) / Z0MAX(I)) 
        FH2(I) = LOG((ZTMAX(I)+2.) / ZTMAX(I)) 
        FH10(I) = LOG((ZTMAX(I)+10.) / ZTMAX(I)) 
       endif 
      ENDDO 
!##DG  IF(LAT.EQ.LATD) THEN 
!##DG    PRINT *, ' DTV, RB(I), FM(I), FH(I), HLINF =', 
!##DG &   dtv, rb, FM(I), FH(I), hlinf 
!##DG  ENDIF 
! 
!  STABLE CASE 
! 
      DO I = 1, IMj 
       if(flag_iter(i)) then  
        IF(DTV(I).GE.0.) THEN 
          HL1(I) = HLINF(I) 
        ENDIF 
        IF(DTV(I).GE.0..AND.HLINF(I).GT..25) THEN 
          HL0INF = Z0MAX(I) * HLINF(I) / Z1(I) 
          HLTINF = ZTMAX(I) * HLINF(I) / Z1(I) 
          AA = SQRT(1. + 4. * ALPHA * HLINF(I)) 
          AA0 = SQRT(1. + 4. * ALPHA * HL0INF) 
          BB = AA 
          BB0 = SQRT(1. + 4. * ALPHA * HLTINF) 
          PM(I) = AA0 - AA + LOG((AA + 1.) / (AA0 + 1.)) 
          PH(I) = BB0 - BB + LOG((BB + 1.) / (BB0 + 1.)) 
          FMS = FM(I) - PM(I) 
          FHS = FH(I) - PH(I) 
          HL1(I) = FMS * FMS * RB(I) / FHS 
        ENDIF 
       endif 
      ENDDO 
! 
!  SECOND ITERATION 
! 
      DO I = 1, IMj 
       if(flag_iter(i)) then  
        IF(DTV(I).GE.0.) THEN 
          HL0 = Z0MAX(I) * HL1(I) / Z1(I) 
          HLT = ZTMAX(I) * HL1(I) / Z1(I) 
          AA = SQRT(1. + 4. * ALPHA * HL1(I)) 
          AA0 = SQRT(1. + 4. * ALPHA * HL0) 
          BB = AA 
          BB0 = SQRT(1. + 4. * ALPHA * HLT) 
          PM(I) = AA0 - AA + LOG((AA + 1.) / (AA0 + 1.)) 
          PH(I) = BB0 - BB + LOG((BB + 1.) / (BB0 + 1.)) 
          HL110 = HL1(I) * 10. / Z1(I) 
          AA = SQRT(1. + 4. * ALPHA * HL110) 
          PM10(I) = AA0 - AA + LOG((AA + 1.) / (AA0 + 1.)) 
          HL12(I) = HL1(I) * 2. / Z1(I) 
!         AA = SQRT(1. + 4. * ALPHA * HL12(I)) 
          BB = SQRT(1. + 4. * ALPHA * HL12(I)) 
          PH2(I) = BB0 - BB + LOG((BB + 1.) / (BB0 + 1.)) 
          BB = SQRT(1. + 4. * ALPHA * HL110) 
          PH10(I) = BB0 - BB + LOG((BB + 1.) / (BB0 + 1.)) 
        ENDIF 
       endif 
      ENDDO 
!! 
!##DG  IF(LAT.EQ.LATD) THEN 
!##DG    PRINT *, ' HL1(I), PM, PH =', 
!##DG &   HL1(I),  pm, ph 
!##DG  ENDIF 
! 
!  UNSTABLE CASE 
! 
! 
!  CHECK FOR UNPHYSICAL OBUKHOV LENGTH 
! 
      DO I=1,IMj 
       if(flag_iter(i)) then  
        IF(DTV(I).LT.0.) THEN 
          OLINF = Z1(I) / HLINF(I) 
          IF(ABS(OLINF).LE.50. * Z0MAX(I)) THEN 
            HLINF(I) = -Z1(I) / (50. * Z0MAX(I)) 
          ENDIF 
        ENDIF 
       endif 
      ENDDO 
! 
!  GET PM AND PH 
! 
      DO I = 1, IMj 
       if(flag_iter(i)) then  
        IF(DTV(I).LT.0..AND.HLINF(I).GE.-.5) THEN 
          HL1(I) = HLINF(I) 
          PM(I) = (A0 + A1 * HL1(I)) * HL1(I)   &
                  / (1. + B1 * HL1(I) + B2 * HL1(I) * HL1(I)) 
          PH(I) = (A0P + A1P * HL1(I)) * HL1(I) &
                  / (1. + B1P * HL1(I) + B2P * HL1(I) * HL1(I)) 
          HL110 = HL1(I) * 10. / Z1(I) 
          PM10(I) = (A0 + A1 * HL110) * HL110   &
                  / (1. + B1 * HL110 + B2 * HL110 * HL110) 
          HL12(I) = HL1(I) * 2. / Z1(I) 
          PH2(I) = (A0P + A1P * HL12(I)) * HL12(I)  &
                  / (1. + B1P * HL12(I) + B2P * HL12(I) * HL12(I)) 
          PH10(I) = (A0P + A1P * HL110) * HL110  &
                  / (1. + B1P * HL110 + B2P * HL110 * HL110) 
        ENDIF 
        IF(DTV(I).LT.0.AND.HLINF(I).LT.-.5) THEN 
          HL1(I) = -HLINF(I) 
          PM(I) = LOG(HL1(I)) + 2. * HL1(I) ** (-.25) - .8776 
          PH(I) = LOG(HL1(I)) + .5 * HL1(I) ** (-.5) + 1.386 
          HL110 = HL1(I) * 10. / Z1(I) 
          PM10(I) = LOG(HL110) + 2. * HL110 ** (-.25) - .8776 
          HL12(I) = HL1(I) * 2. / Z1(I) 
          PH2(I) = LOG(HL12(I)) + .5 * HL12(I) ** (-.5) + 1.386 
          PH10(I) = LOG(HL110) + .5 * HL110 ** (-.5) + 1.386 
        ENDIF 
       endif 
      ENDDO 
! 
!  FINISH THE EXCHANGE COEFFICIENT COMPUTATION TO PROVIDE FM AND FH 
! 
      DO I = 1, IMj 
       if(flag_iter(i)) then  
        FM(I) = FM(I) - PM(I) 
        FH(I) = FH(I) - PH(I) 
        FM10(I) = FM10(I) - PM10(I) 
        FH2(I) = FH2(I) - PH2(I) 
        FH10(I) = FH10(I) - PH10(I) 
        CM(I) = CA * CA / (FM(I) * FM(I)) 
        CH(I) = CA * CA / (FM(I) * FH(I)) 
        CQ = CH(I) 
        STRESS(I) = CM(I) * WIND(I) * WIND(I) 
        USTAR(I)  = SQRT(STRESS(I)) 
!       USTAR(I) = SQRT(CM(I) * WIND(I) * WIND(I)) 
       endif 
      ENDDO 
!##DG  IF(LAT.EQ.LATD) THEN 
!##DG    PRINT *, ' FM, FH, CM, CH(I), USTAR =', 
!##DG &   FM, FH, CM, ch, USTAR 
!##DG  ENDIF 
! 
!  UPDATE Z0 OVER OCEAN 
! 
      DO I = 1, IMj 
       if(flag_iter(i)) then  
        IF(SLIMSK(I).EQ.0.) THEN 
          Z0(I) = (CHARNOCK / G) * USTAR(I) ** 2 
! MBEK -Toga coare flux algorithem
!         Z0(I) = (CHARNOCK / G) * USTAR(I) ** 2 + arnu/ustar(i) 
!  NEW IMPLEMENTATION OF Z0 
!         CC = USTAR(I) * Z0 / RNU 
!         PP = CC / (1. + CC) 
!         FF = G * ARNU / (CHARNOCK * USTAR(I) ** 3) 
!         Z0 = ARNU / (USTAR(I) * FF ** PP) 
          Z0(I) = MIN(Z0(I),.1) 
          Z0(I) = MAX(Z0(I),1.E-7) 
          Z0RL(I) = 100. * Z0(I) 
        ENDIF 
       endif 
      ENDDO 
 
!      if(j.eq.jo)then
!      print*,'final-----'
!      print*,'ustar=',ustar(io)
!      print*,'z0(step1)=',charnock/g*ustar(io)*ustar(io)
!      print*, 'z0(step2)=',z0(io)
!      print*,'-----------'
!      endif
      RETURN 
      END
