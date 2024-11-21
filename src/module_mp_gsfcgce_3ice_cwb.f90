!WRF:MODEL_LAYER:PHYSICS
!
      MODULE module_mp_gsfcgce_3ice_cwb

!      USE module_wrf_error
!      USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm
!      USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep

      REAL, PRIVATE :: al,cp,rd1,rd2,grvt,cd1,cd2,t0,t00,alv,alf,als,  &
                       avc,afc,asc,rw,cw,ci,c76,c358,c172,c409,c218,   &
                       c580,c610,c149,c879,c141,as,bs,aw,bw,roqr
      INTEGER, PRIVATE, PARAMETER :: SPECTRUM = 1                       ! 0.MP;1.GAMMA+MAPPING;2.MP+MAPPING SPECTRUM
      INTEGER, PRIVATE, PARAMETER :: VTRX = 1                           ! 0.FIXED;1.NUMS;2.ANAL;3.GCE-4ICE RAINDROPS FALL SPEED
      REAL, PRIVATE, PARAMETER :: THRD = 1./3.
      REAL, PRIVATE, PARAMETER :: KMIN = 0.2222, KMAX = 0.8396          ! KMAX = 0.9441; KMAX =0.8396 for 15
      REAL, PRIVATE, PARAMETER :: cmin1 = 1.E-12, cmin2 = 1.E-8
      REAL, PRIVATE, PARAMETER :: MDR1 = 33.999045, MDR2 = -13.813047
      REAL, PRIVATE, PARAMETER :: MDR3 = 0.29834969, MDR4 = 1.6622025
      REAL, PRIVATE, PARAMETER :: MDR5 = 1.2480129E-3, MDR6 = -0.104641
      REAL, PRIVATE, PARAMETER :: MDS1 = 87.049912, MDS2 = -41.711453
      REAL, PRIVATE, PARAMETER :: MDS3 = 2.1772518, MDS4 = 4.9657494
      REAL, PRIVATE, PARAMETER :: MDS5 = 3.5452471E-3, MDS6 = -0.45465148
      REAL, PRIVATE, PARAMETER :: MDG1 = 231.32157, MDG2 = -93.387843
      REAL, PRIVATE, PARAMETER :: MDG3 = 1.810578, MDG4 = 9.5907234
      REAL, PRIVATE, PARAMETER :: MDG5 = 1.683093E-3, MDG6 = -0.37378561
      REAL, PRIVATE, PARAMETER :: MDH1 = 12.676969, MDH2 = -2.8835883
      REAL, PRIVATE, PARAMETER :: MDH3 = -3.8803716E-2, MDH4 = 0.43033064
      REAL, PRIVATE, PARAMETER :: MDH5 = 6.53219E-3, MDH6 = -6.4118002E-2
      REAL, PRIVATE, PARAMETER :: ROS1 = -87.375922, ROS2 = 40.320759
      REAL, PRIVATE, PARAMETER :: ROS3 = -1.432207, ROS4 = -4.2850883
      REAL, PRIVATE, PARAMETER :: ROS5 = -4.4996627E-3, ROS6 = 0.29658535
      REAL, PRIVATE, PARAMETER :: EFD1 = 1.6855156, EFD2 = 0.84147302
      REAL, PRIVATE, PARAMETER :: EFD3 = -0.46783369, EFD4 = 1.005305E-2
      REAL, PRIVATE, PARAMETER :: EFD5 = 3.4833332E-2, EFD6 = 1.4727223E-2
      REAL, PRIVATE, PARAMETER :: VTR1 = 8.07534, VTR2 = -1.7939129E-3
      REAL, PRIVATE, PARAMETER :: VTR3 = -4.1755472, VTR4 = 0.53826913
      REAL, PRIVATE, PARAMETER :: VTR5 = -5.1104358E-2, VTR6 = -0.37871377
      REAL, PRIVATE, PARAMETER :: VTR7 = 5.3815649E-2
      REAL, PARAMETER :: xnor = 8.E6
      REAL, PARAMETER :: xnos = 1.6E7
      REAL, PARAMETER :: xnoh = 2.E5
      REAL, PARAMETER :: xnog = 4.E6
      REAL, PARAMETER :: rhograul = 400.
      REAL, PARAMETER :: rhohail = 900.
!+---+-----------------------------------------------------------------+
CONTAINS
!-------------------------------------------------------------------
!  NASA/GSFC GCE
!  Tao et al, 2001, Meteo. & Atmos. Phy., 97-137
!-------------------------------------------------------------------
      SUBROUTINE gsfcgce_3ice_cwb                                      &
                        (th,qv,ql,qr,qi,qs,qg,rho,pii,p,dt_in,z,ht,    &
                         dz8w,itimestep,                               &
!                         ids,ide,jds,jde,kds,kde,                      &
                         ims,ime,jms,jme,kms,kme,                      &
                         its,ite,jts,jte,kts,kte,                      &
                         sr,rainnc,rainncv,snownc,snowncv,graupelnc,   &
                         graupelncv,xland)
!-------------------------------------------------------------------
      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
      INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
      INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte,itimestep
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th,qv,&
                                                ql,qr,qi,qs,qg
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: rho,pii, &
                                                              p,dz8w,z
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: rainnc,       &
            rainncv,snownc,snowncv,sr,graupelnc,graupelncv
      REAL, DIMENSION(its:ite,kts:kte,jts:jte) :: SMLF,GMLF
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: ht,xland
      REAL, INTENT(IN) :: dt_in
      INTEGER :: i,j,k,i24h,ihail,ice2

      ihail = 2
      ice2 = 0
      i24h = nint(86400./dt_in)
      if (mod(itimestep,i24h).eq.1) then
         write(6,*) 'ihail=',ihail,'  ice2=',ice2
         if (ice2.eq.0) then
            write(6,*) 'Running 3-ice scheme in GSFCGCE with'
            if (ihail.eq.0) then 
               write(6,*) '     ice, snow and graupel'
            elseif (ihail.eq.1) then
               write(6,*) '     ice, snow and hail'
            elseif (ihail.eq.2) then
               write(6,*) '     ice, snow (graupel or hail)'
            else
               write(6,*) 'ihail has to be either 0,1,or 2'
               stop
            endif !ihail
         elseif (ice2.eq.1) then
            write(6,*) 'Running 2-ice scheme in GSFCGCE with'
            write(6,*) '     ice and snow'
         elseif (ice2.eq.2) then
            write(6,*) 'Running 2-ice scheme in GSFCGCE with'
            write(6,*) '     ice and graupel'
         elseif (ice2.eq.3) then
            write(6,*) 'Running warm rain only scheme without any ice'
         else
            write(6,*) 'gsfcgce_2ice has to be 0, 1, 2, or 3'
            stop
         endif !ice2
      endif !itimestep
      CALL consat_s

      ! microphysics in GCE
      CALL SATICEL_S(dt_in,ihail,ice2,th,qv,ql,qr,qi,qs,qg,rho,pii,p,  &
                     itimestep,SMLF,GMLF,xland,                        &
!                     ids,ide,jds,jde,kds,kde,                          &
                     ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte) 

      CALL fall_flux(dt_in,ihail,ice2,th,pii,ql,qr,qi,qs,qg,p,rho,z,   &
                     dz8w,ht,itimestep,SMLF,GMLF,sr,rainnc,rainncv,    &
                     snownc,snowncv,graupelnc,graupelncv,ims,ime,jms,  &
                     jme,kms,kme,its,ite,jts,jte,kts,kte)

      END SUBROUTINE gsfcgce_3ice_cwb
!=======================================================================
!-----------------------------------------------------------------------
      SUBROUTINE fall_flux (dt,ihail,ice2,TH,PII,qc,qr,qi,qs,qg,p,rho, &
                            z,dz8w,topo,itimestep,SMLF3D,GMLF3D,sr,    &
                            rainnc,rainncv,snownc,snowncv,graupelnc,   &
                            graupelncv,ims,ime,jms,jme,kms,kme,its,ite,&
                            jts,jte,kts,kte)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ihail,ice2,ims,ime,jms,jme,kms,kme,its,   &
                             ite,jts,jte,kts,kte,itimestep
      REAL, INTENT(IN) :: dt
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: topo
      REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: sr,rainnc,    &
            rainncv,snownc,snowncv,graupelnc,graupelncv
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: qr,qi,&
                                                                 qs,qg
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: qc,TH,   &
                                                PII,rho,z,dz8w,p
      REAL, DIMENSION(its:ite,kts:kte,jts:jte), INTENT(INOUT) ::       &
                                                SMLF3D,GMLF3D
      REAL, DIMENSION(kts:kte) :: TK1D,TC1D,qcz,DSLL,RHOS,RHOG,SMLF,   &
            GMLF,N0R,N0S,N0G,N0H,MVDR,EFDR,MVDS,EFDS,MVDG,EFDG,MVDH,   &
            EFDH,KDXR,KDXS,KDXG,KDXH,AFAR,AFAS,AFAG,LAMR,LAMS,LAMG,    &
            MVDM,AVR,BVR,LTK,LTK2,LQR,LQR2,LQS,LQS2,LQG,LQG2,MVG0
! temperary vars
      INTEGER :: ihail2,k,i,j,min_q,max_q
      REAL :: tmp1,term0,pptrain,pptsnow,pptgraul,pptice,dtb,pi,consta,&
              constc,t_del_tv,del_tv,flux,fluxin,fluxout,tmpqrz
      REAL, DIMENSION(kts:kte) :: sqrhoz,qrz,qiz,qsz,qgz,zz,dzw,prez,  &
                                  rhoz,orhoz,vtr,vts,vtg,vti
      REAL, PARAMETER :: constb = 0.8, constd = 0.41, o6 = 1./6.,      &
                         cdrag = 0.6, abar = 19.3, bbar = 0.37,        &
                         p0 = 1.0e5, rhoe_s = 1.29
! for terminal velocity flux
      LOGICAL :: notlast
!-----------------------------------------------------------------------
!  This program calculates precipitation fluxes due to terminal velocities.
!-----------------------------------------------------------------------
      dtb = dt
      pi = ACOS(-1.)
      consta = 841.9967                        ! 2115.0*0.01**(1-constb)
      constc = 11.72
!***********************************************************************
! Calculate precipitation fluxes due to terminal velocities.
!***********************************************************************
!- Calculate termianl velocity (vt?)  of precipitation q?z
!- Find maximum vt? to determine the small delta t
      ihail2 = 0
      DO j = jts,jte   ! j_loop
      DO i = its,ite   ! i_loop:
      pptrain = 0.
      pptsnow = 0.
      pptgraul = 0.
      pptice  = 0.
      DO k = kts,kte
         qcz(k) = qc(i,k,j)
         qrz(k) = qr(i,k,j)
         qiz(k) = qi(i,k,j)
         qsz(k) = qs(i,k,j)
         rhoz(k) = rho(i,k,j)
         orhoz(k) = 1./rhoz(k)
         prez(k) = p(i,k,j)
         sqrhoz(k) = SQRT(rhoe_s/rhoz(k))
         zz(k) = z(i,k,j)
         dzw(k) = dz8w(i,k,j)
      ENDDO ! end of k-loop 

      IF (ice2.eq.0) THEN
         DO k = kts,kte
            qgz(k) = qg(i,k,j)
         ENDDO
      ELSE
         DO k = kts,kte
            qgz(k) = 0.
         ENDDO
      ENDIF
      !-------------------------------Tsai-----------------------------------------------
      DO k = kts, kte
         TK1D(k) = TH(i,k,j)*PII(i,k,j)
         TC1D(k) = TK1D(k)-t0
         RHOS(k) = 100.
         RHOG(k) = rhograul
         MVG0(k) = 0.
         AFAR(k) = 0.
         AFAS(k) = 0.
         AFAG(k) = 0.
         LAMR(k) = 0.
         LAMS(k) = 0.
         LAMG(k) = 0.
         vtr(k)  = 0.
         SMLF(k) = MAX(0.,SMLF3D(i,k,j))
         GMLF(k) = MAX(0.,GMLF3D(i,k,j))
      ENDDO
      !-------------------------------Tsai-----------------------------------------------
      !-- rain
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         IF (qrz(k).GE.cmin2) THEN
            min_q = MIN0(min_q,k)
            max_q = MAX0(max_q,k)
            IF (SPECTRUM.EQ.0) THEN
               AFAR(k) = 0.
               N0R(k) = LOG(xnor)
               LAMR(k) = EXP((N0R(k)+LOG(pi*1000./rhoz(k)/qrz(k)/6.)+  &
                         LGAMMA(AFAR(k)+4.))/(AFAR(k)+4.))
            ELSEIF (SPECTRUM.EQ.1) THEN
               LTK(k)  = LOG(TK1D(k))
               LQR(k)  = -1.*LOG(rhoz(k)*qrz(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQR2(k) = LQR(k)*LQR(k)
               MVDR(k) = EXP(MDR1+MDR2*LTK(k)+MDR3*LQR(k)+MDR4*        &
                         LTK2(k)+MDR5*LQR2(k)+MDR6*LTK(k)*LQR(k))
               EFDR(k) = EXP(EFD1+EFD2*LOG(MVDR(k))+EFD3*LOG(1.E3)+    &
                         EFD4*LOG(MVDR(k))**2.+EFD5*LOG(1.E3)**2.+     &
                         EFD6*LOG(MVDR(k))*LOG(1.E3))
               EFDR(k) = MAX(MIN(EFDR(k),MVDR(k)/KMIN**THRD),MVDR(k)/  &
                         KMAX**THRD)
               KDXR(k) = MAX(KMIN,MIN(KMAX,(MVDR(k)/EFDR(k))**3.))
               AFAR(k) = MAX(0.,(6.*KDXR(k)-3.+SQRT(8.*KDXR(k)+        &
                         1.))/(2.-2.*KDXR(k)))
               LAMR(k) = (AFAR(k)+3.)/EFDR(k)*1.E6
               N0R(k) = LOG(6.*rhoz(k)*qrz(k)/pi/1000.)+(4.+AFAR(k))*  &
                        LOG(LAMR(k))-LGAMMA(AFAR(k)+4.)
            ELSEIF (SPECTRUM.EQ.2) THEN
               MVDR(k) = EXP(8.425904+0.25*LOG(rhoz(k)*qrz(k)))
               AFAR(k) = 0.
               LAMR(k) = (EXP(LGAMMA(AFAR(k)+4.)-LGAMMA(AFAR(k)+1.)))**&
                         THRD/MVDR(k)*1.E6
               N0R(k) = LOG(6.*rhoz(k)*qrz(k)/pi/1000.)+(4.+AFAR(k))*  &
                        LOG(LAMR(k))-LGAMMA(AFAR(k)+4.)
            ENDIF

            IF (VTRX.EQ.0) THEN
               vtr(k) = consta*EXP(LGAMMA(constb+AFAR(k)+4.)-LGAMMA(   &
                        AFAR(k)+4.)-constb*LOG(LAMR(k)))*sqrhoz(k)
            ELSEIF (VTRX.EQ.1) THEN
               LTK(k)  = LOG(TK1D(k))
               LQR(k)  = LOG(-1.*LOG(qrz(k)*rhoz(k)))
               LQR2(k) = LQR(k)*LQR(k)
               vtr(k)  = MIN(9.1,EXP((VTR1+VTR2*LTK(k)+VTR3*LQR(k)+    &
                         VTR4*LQR2(k))/(1.+VTR5*LTK(k)+VTR6*LQR(k)+    &
                         VTR7*LQR2(k)))/1.E4)*sqrhoz(k)
            ELSEIF (VTRX.EQ.2) THEN
               vtr(k) = CHEN(rhoz(k),AFAR(k)+3.,LAMR(k))
            ELSEIF (VTRX.EQ.3) THEN
               vtr(k) = MAX(0.,-26.7+2.06E6/LAMR(k)-2.045E9/LAMR(k)**  &
                        2.+9.06E11/LAMR(k)**3.)*1.E-2*sqrhoz(k)
            ENDIF
            if (k.eq.1) then
!               del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtr(k))
               del_tv = DMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtr(k))
            else
!               del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtr(k))
               del_tv = DMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtr(k))
            endif
         else
            vtr(k) = 0.
         endif  ! qrz
      ENDDO ! end of k-loop
      if (max_q.ge.min_q) then
          !- Check if the summation of the small delta t >=  big delta t
          !             (t_del_tv)          (del_tv)             (dtb)
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
          ! use small delta t to calculate the qrz flux
          ! termi is the qrz flux pass in the grid box through the upper boundary
          ! termo is the qrz flux pass out the grid box through the lower boundary
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vtr(k)*qrz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qrz(k) = qrz(k)+del_tv*flux
!            qrz(k) = AMAX1(0.,qrz(k))
            qrz(k) = DMAX1(0.,qrz(k))
            qr(i,k,j) = qrz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptrain = pptrain+fluxin*del_tv
         else
            k = min_q-1
            qrz(k) = qrz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qr(i,k,j) = qrz(k)
         endif
      else
           notlast = .false.
      endif
      ENDDO  ! eno of DO while (notlast)
      
      !-- snow
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
      min_q = kte
      max_q = kts-1
      DO k = kts,kte
         IF (qsz(k).GE.cmin2) THEN
            min_q = MIN0(min_q,k)
            max_q = MAX0(max_q,k)
            IF (SPECTRUM.EQ.0) THEN
               RHOS(k) = 100.
               AFAS(k) = 0.
               N0S(k) = LOG(xnos)
               LAMS(k) = EXP((N0S(k)+LOG(pi*RHOS(k)/rhoz(k)/qsz(k)/6.)+&
                         LGAMMA(AFAS(k)+4.))/(AFAS(k)+4.))
               SMLF(k) = MIN(SMLF(k),1.-MIN(1.,RHOS(k)/1.E3))
            ELSEIF (SPECTRUM.EQ.1) THEN
               LTK(k) = LOG(TK1D(k))
               LQS(k) = -1.*LOG(rhoz(k)*qsz(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQS2(k) = LQS(k)*LQS(k)
               RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*  &
                         LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)* &
                         LQS(k))))
               SMLF(k) = MIN(SMLF(k),1.-MIN(1.,RHOS(k)/1.E3))
               RHOS(k) = 1.E3*SMLF(k)+(1.-SMLF(k))*RHOS(k)
               RHOS(k) = MIN(MAX(RHOS(k),100.),910.)
               MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+&
                         MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
               EFDS(k) = EXP(EFD1+EFD2*LOG(MVDS(k))+EFD3*LOG(RHOS(k))+ &
                         EFD4*LOG(MVDS(k))**2.+EFD5*LOG(RHOS(k))**2.+  &
                         EFD6*LOG(MVDS(k))*LOG(RHOS(k)))
               EFDS(k) = MAX(MIN(EFDS(k),MVDS(k)/KMIN**THRD),MVDS(k)/  &
                         KMAX**THRD)
               KDXS(k) = MAX(KMIN,MIN(KMAX,(MVDS(k)/EFDS(k))**3.))
               AFAS(k) = MAX(0.,(6.*KDXS(k)-3.+SQRT(8.*KDXS(k)+        &
                         1.))/(2.-2.*KDXS(k)))
               LAMS(k) = (AFAS(k)+3.)/EFDS(k)*1.E6
               N0S(k) = LOG(6.*rhoz(k)*qsz(k)/pi/RHOS(k))+(4.+AFAS(k))*&
                        LOG(LAMS(k))-LGAMMA(AFAS(k)+4.)
            ELSEIF (SPECTRUM.EQ.2) THEN
               LTK(k) = LOG(TK1D(k))
               LQS(k) = -1.*LOG(rhoz(k)*qsz(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQS2(k) = LQS(k)*LQS(k)
               RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*  &
                         LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)* &
                         LQS(k))))
               SMLF(k) = MIN(SMLF(k),1.-MIN(1.,RHOS(k)/1.E3))
               RHOS(k) = 1.E3*SMLF(k)+(1.-SMLF(k))*RHOS(k)
               RHOS(k) = MIN(MAX(RHOS(k),100.),910.)
               MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+&
                         MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
               AFAS(k) = 0.
               LAMS(k) = (EXP(LGAMMA(AFAS(k)+4.)-LGAMMA(AFAS(k)+1.)))**&
                         THRD/MVDS(k)*1.E6
               N0S(k) = LOG(6.*rhoz(k)*qsz(k)/pi/RHOS(k))+(4.+AFAS(k))*&
                        LOG(LAMS(k))-LGAMMA(AFAS(k)+4.)
            ENDIF
            vts(k) = constc*sqrhoz(k)*EXP(LGAMMA(constd+4.+AFAS(k))-   &
                     LGAMMA(AFAS(k)+4.)-constd*LOG(LAMS(k)))
            if (k.eq.1) then
!               del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vts(k))
               del_tv = DMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vts(k))
            else
!               del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vts(k))
               del_tv = DMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vts(k))
            endif
         else
            vts(k) = 0.
         endif
      enddo
      if (max_q.ge.min_q) then
          !- Check if the summation of the small delta t >=  big delta t
          !             (t_del_tv)          (del_tv)             (dtb)
         t_del_tv = t_del_tv+del_tv
         if (t_del_tv.ge.dtb) then
            notlast = .false.
            del_tv = dtb+del_tv-t_del_tv
         endif
          ! use small delta t to calculate the qsz flux
          ! termi is the qsz flux pass in the grid box through the upper boundary
          ! termo is the qsz flux pass out the grid box through the lower boundary
         fluxin = 0.
         do k = max_q,min_q,-1
            fluxout = rhoz(k)*vts(k)*qsz(k)
            flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
            qsz(k) = qsz(k)+del_tv*flux
!            qsz(k) = AMAX1(0.,qsz(k))
            qsz(k) = DMAX1(0.,qsz(k))
            qs(i,k,j) = qsz(k)
            fluxin = fluxout
         enddo
         if (min_q.eq.1) then
            pptsnow = pptsnow+fluxin*del_tv
         else
            k = min_q-1
            qsz(k) = qsz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
            qs(i,k,j) = qsz(k)
         endif
      else
         notlast=.false.
      endif
      ENDDO  ! end of DO while (notlast)
      if (ice2.eq.0) then 
         t_del_tv = 0.
         del_tv = dtb
         notlast = .true.
         DO while (notlast)
         min_q = kte
         max_q = kts-1
         do k = kts,kte
            IF (qgz(k).GE.cmin2) THEN
               MVG0(k) = EXP(8.8283+0.25*LOG(rhoz(k)*qgz(k)))/1.E6
               DSLL(k) = MAX(1.E-3,MIN(1.,2.*1.E-2*(EXP(MIN(20.,       &
                         -TC1D(k)/MAX(0.1,(1.1E4*rhoz(k)*(qcz(k)+      &
                         qrz(k))-1.3E3*rhoz(k)*qgz(k)+1.))))-1.)))
               if (ihail.eq.2) then
                  IF (MVG0(k).GT.DSLL(k)) THEN
                     ihail2 = 1
                  ELSE
                     ihail2 = 0
                  ENDIF
               ENDIF
               if (ihail.eq.1.or.ihail2.eq.1) then ! Hail
                  min_q = MIN0(min_q,k)
                  max_q = MAX0(max_q,k)
                  RHOG(k) = rhohail
                  GMLF(k) = MIN(GMLF(k),1.-MIN(1.,RHOG(k)/1.E3))
                  IF (SPECTRUM.EQ.0) THEN
                     AFAG(k) = 0.
                     N0H(k) = LOG(xnoh)
                     LAMG(k) = EXP((N0H(k)+LOG(pi*rhohail/rhoz(k)/     &
                               qgz(k)/6.)+LGAMMA(AFAG(k)+4.))/(AFAG(k)+&
                               4.))
                  ELSEIF (SPECTRUM.EQ.1) THEN
                     LTK(k) = LOG(TK1D(k))
                     LQG(k) = -1.*LOG(rhoz(k)*qgz(k))
                     LTK2(k) = LTK(k)*LTK(k)
                     LQG2(k) = LQG(k)*LQG(k)
                     MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*  &
                               LTK2(k)+MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
                     EFDH(k) = EXP(EFD1+EFD2*LOG(MVDH(k))+EFD3*        &
                               LOG(9.E2)+EFD4*LOG(MVDH(k))**2.+EFD5*   &
                               LOG(9.E2)**2.+EFD6*LOG(MVDH(k))*        &
                               LOG(9.E2))
                     EFDH(k) = MAX(MIN(EFDH(k),MVDH(k)/KMIN**THRD),    &
                               MVDH(k)/KMAX**THRD)
                     KDXH(k) = MAX(KMIN,MIN(KMAX,(MVDH(k)/EFDH(k))**3.))
                     AFAG(k) = MAX(0.,(6.*KDXH(k)-3.+SQRT(8.*KDXH(k)+  &
                               1.))/(2.-2.*KDXH(k)))
                     LAMG(k) = (AFAG(k)+3.)/EFDH(k)*1.E6
                     N0H(k) = LOG(6.*rhoz(k)*qgz(k)/pi/rhohail)+(4.+   &
                              AFAG(k))*LOG(LAMG(k))-LGAMMA(AFAG(k)+4.)
                  ELSEIF (SPECTRUM.EQ.2) THEN
                     LTK(k) = LOG(TK1D(k))
                     LQG(k) = -1.*LOG(rhoz(k)*qgz(k))
                     LTK2(k) = LTK(k)*LTK(k)
                     LQG2(k) = LQG(k)*LQG(k)
                     MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*  &
                               LTK2(k)+MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
                     AFAG(k) = 0.
                     LAMG(k) = (EXP(LGAMMA(AFAG(k)+4.)-LGAMMA(AFAG(k)+ &
                               1.)))**THRD/MVDH(k)*1.E6
                     N0H(k) = LOG(6.*rhoz(k)*qgz(k)/pi/rhohail)+(4.+   &
                              AFAG(k))*LOG(LAMG(k))-LGAMMA(AFAG(k)+4.)
                  ENDIF
                  vtg(k) = SQRT(4.*9.81*rhohail/3./rhoz(k)/cdrag)*     &
                           EXP(LGAMMA(4.5+AFAG(k))-LGAMMA(AFAG(k)+4.)- &
                           0.5*LOG(LAMG(k)))
                  if (k.eq.1) then
!                     del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
                     del_tv = DMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
                  else
!                      del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
                      del_tv = DMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
                  endif !k
               elseif (ihail.eq.0.or.ihail2.eq.0) then  ! Graupel
                  min_q = MIN0(min_q,k)
                  max_q = MAX0(max_q,k)
                  IF (SPECTRUM.EQ.0) THEN
                     AFAG(k) = 0.
                     RHOG(k) = rhograul
                     N0G(k) = LOG(xnog)
                     LAMG(k) = EXP((N0G(k)+LOG(pi*rhograul/rhoz(k)/    &
                               qgz(k)/6.)+LGAMMA(AFAG(k)+4.))/(AFAG(k)+&
                               4.))
                     GMLF(k) = MIN(GMLF(k),1.-MIN(1.,RHOG(k)/1.E3))
                  ELSEIF (SPECTRUM.EQ.1) THEN
                     LTK(k) = LOG(TK1D(k))
                     LQG(k) = -1.*LOG(rhoz(k)*qgz(k))
                     LTK2(k) = LTK(k)*LTK(k)
                     LQG2(k) = LQG(k)*LQG(k)
                     GMLF(k) = MIN(GMLF(k),1.-MIN(1.,RHOG(k)/1.E3))
                     RHOG(k) = MIN(900.,1.E3*GMLF(k)+(1.-GMLF(k))*MAX( &
                               250.,EXP(-0.056433*LQG(k)+6.608048)))
                     MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*  &
                               LTK2(k)+MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
                     EFDG(k) = EXP(EFD1+EFD2*LOG(MVDG(k))+EFD3*        &
                               LOG(RHOG(k))+EFD4*LOG(MVDG(k))**2.+     &
                               EFD5*LOG(RHOG(k))**2.+EFD6*LOG(MVDG(k))*&
                               LOG(RHOG(k)))
                     EFDG(k) = MAX(MIN(EFDG(k),MVDG(k)/KMIN**THRD),    &
                               MVDG(k)/KMAX**THRD)
                     KDXG(k) = MAX(KMIN,MIN(KMAX,(MVDG(k)/EFDG(k))**3.))
                     AFAG(k) = MAX(0.,(6.*KDXG(k)-3.+SQRT(8.*KDXG(k)+  &
                               1.))/(2.-2.*KDXG(k)))
                     LAMG(k) = (AFAG(k)+3.)/EFDG(k)*1.E6
                     N0G(k) = LOG(6.*rhoz(k)*qgz(k)/pi/RHOG(k))+(4.+   &
                               AFAG(k))*LOG(LAMG(k))-LGAMMA(AFAG(k)+4.)
                  ELSEIF (SPECTRUM.EQ.2) THEN
                     LTK(k) = LOG(TK1D(k))
                     LQG(k) = -1.*LOG(rhoz(k)*qgz(k))
                     LTK2(k) = LTK(k)*LTK(k)
                     LQG2(k) = LQG(k)*LQG(k)
                     GMLF(k) = MIN(GMLF(k),1.-MIN(1.,RHOG(k)/1.E3))
                     RHOG(k) = MIN(900.,1.E3*GMLF(k)+(1.-GMLF(k))*MAX( &
                               250.,EXP(-0.056433*LQG(k)+6.608048)))
                     MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*  &
                               LTK2(k)+MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
                     AFAG(k) = 0.
                     LAMG(k) = (EXP(LGAMMA(AFAG(k)+4.)-LGAMMA(AFAG(k)+ &
                               1.)))**THRD/MVDG(k)*1.E6
                     N0G(k) = LOG(6.*rhoz(k)*qgz(k)/pi/RHOG(k))+(4.+   &
                              AFAG(k))*LOG(LAMG(k))-LGAMMA(AFAG(k)+4.)
                  ENDIF
                  vtg(k) = abar*(p0/prez(k))**0.4*EXP(LGAMMA(bbar+4.+  &
                           AFAG(k))-LGAMMA(AFAG(k)+4.)-bbar*           &
                           LOG(LAMG(k)))
                  if (k.eq.1) then
!                     del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
                     del_tv = DMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
                  else
!                     del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
                     del_tv = DMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
                  endif !k
               endif !ihail
            else
               vtg(k) = 0.
            endif !qgz
         enddo ! end of k-loop
         if (max_q.ge.min_q) then
            !- Check if the summation of the small delta t >=  big delta t
            !             (t_del_tv)          (del_tv)             (dtb)
            t_del_tv = t_del_tv+del_tv
            if (t_del_tv.ge.dtb) then
               notlast = .false.
               del_tv = dtb+del_tv-t_del_tv
            endif
            ! use small delta t to calculate the qgz flux
            ! termi is the qgz flux pass in the grid box through the upper boundary
            ! termo is the qgz flux pass out the grid box through the lower boundary
            !
            fluxin = 0.
            do k = max_q,min_q,-1
               fluxout = rhoz(k)*vtg(k)*qgz(k)
               flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
               qgz(k) = qgz(k)+del_tv*flux
!               qgz(k) = AMAX1(0.,qgz(k))
               qgz(k) = DMAX1(0.,qgz(k))
               qg(i,k,j) = qgz(k)
               fluxin = fluxout
            enddo
            if (min_q.eq.1) then
               pptgraul = pptgraul+fluxin*del_tv
            else
               k = min_q-1
               qgz(k) = qgz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
               qg(i,k,j) = qgz(k)
            endif
         else
            notlast = .false.
         endif
         ENDDO  ! end of DO while (notlast)
      ENDIF !ice2
      !
      !-- cloud ice  (03/21/02) follow Vaughan T.J. Phillips at GFDL
      !
      t_del_tv = 0.
      del_tv = dtb
      notlast = .true.
      DO while (notlast)
         min_q = kte
         max_q = kts-1
         do k = kts,kte
            if (qiz(k).ge.1.E-12) then
               min_q = MIN0(min_q,k)
               max_q = MAX0(max_q,k)
               vti(k) = 3.29*(rhoz(k)*qiz(k))**0.16  ! Heymsfield and Donner
               if (k.eq.1) then
!                  del_tv = AMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vti(k))
                  del_tv = DMIN1(del_tv,0.9*(zz(k)-topo(i,j))/vti(k))
               else
!                  del_tv = AMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vti(k))
                  del_tv = DMIN1(del_tv,0.9*(zz(k)-zz(k-1))/vti(k))
               endif
            else
               vti(k) = 0.
            endif
         enddo
         if (max_q.ge.min_q) then
          !- Check if the summation of the small delta t >=  big delta t
          !             (t_del_tv)          (del_tv)             (dtb)
            t_del_tv = t_del_tv+del_tv
            if (t_del_tv.ge.dtb) then
               notlast = .false.
               del_tv = dtb+del_tv-t_del_tv
            endif
            ! use small delta t to calculate the qiz flux
            ! termi is the qiz flux pass in the grid box through the upper boundary
            ! termo is the qiz flux pass out the grid box through the lower boundary
            fluxin = 0.
            do k = max_q,min_q,-1
               fluxout = rhoz(k)*vti(k)*qiz(k)
               flux = (fluxin-fluxout)/rhoz(k)/dzw(k)
               qiz(k) = qiz(k)+del_tv*flux
!               qiz(k) = AMAX1(0.,qiz(k))
               qiz(k) = DMAX1(0.,qiz(k))
               qi(i,k,j) = qiz(k)
               fluxin = fluxout
            enddo
            if (min_q.eq.1) then
               pptice = pptice+fluxin*del_tv
            else
               k = min_q-1
               qiz(k) = qiz(k)+del_tv*fluxin/rhoz(k)/dzw(k)
               qi(i,k,j) = qiz(k)
            endif
         else
            notlast = .false.
         endif
         DO k = kts,kte
            SMLF3D(i,k,j) = SMLF(k)
            GMLF3D(i,k,j) = GMLF(k)
         ENDDO
      ENDDO !notlast

      snowncv(i,j) = pptsnow
      snownc(i,j) = snownc(i,j)+pptsnow
      graupelncv(i,j) = pptgraul
      graupelnc(i,j) = graupelnc(i,j)+pptgraul
      RAINNCV(i,j) = pptrain+pptsnow+pptgraul+pptice
      RAINNC(i,j)  = RAINNC(i,j)+pptrain+pptsnow+pptgraul+pptice
      sr(i,j) = 0.
      if (RAINNCV(i,j).gt.0.) sr(i,j) = (pptsnow+pptgraul+pptice)/     &
                                        RAINNCV(i,j) 

      ENDDO                                                             ! i_loop
      ENDDO                                                             ! j_loop

      END SUBROUTINE fall_flux
!----------------------------------------------------------------
!
!=======================================================================
      SUBROUTINE consat_s
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      al  = 2.5e10
      cp  = 1.00546e7
      rd1 = 1.E-3
      rd2 = 2.2
      grvt = 980.
      cd1 = 6.e-1
      cd2 = 4.*grvt/(3.*cd1)
      t0 = 273.16
      t00 = 238.16
      alv = 2.5e10
      alf = 3.336e9
      als = 2.8336e10
      avc = alv/cp
      afc = alf/cp
      asc = als/cp
      rw = 4.615E6
      cw = 4.187e7
      ci = 2.093e7
      c76 = 7.66
      c358 = 35.86
      c172 = 17.26939
      c409 = 4098.026
      c218 = 21.87456
      c580 = 5807.695
      c610 = 6.1078e3
      c149 = 1.496286e-5
      c879 = 8.794142
      c141 = 1.4144354e7
      as = 177.3894
      bs = 0.41
      aw = 2115.
      bw = 0.8
      roqr = 1.
!
      return
      END SUBROUTINE consat_s
!======================================================================
!
!======================================================================
      SUBROUTINE saticel_s(dt,ihail,ice2,ptwrf,qvwrf,qlwrf,qrwrf,qiwrf,&
                           qswrf,qgwrf,rho_mks,pi_mks,p0_mks,itimestep,&
                           SMLF,GMLF,xland,                            &
!                           ids,ide,jds,jde,kds,kde,                    &
                           ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
!-----------------------------------------------------------------------
      IMPLICIT NONE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!                                                                      c
!   Tao, W.-K., J. J. Shi,  S. Lang, C. Peters-Lidard, A. Hou, S.      c
!   Braun, and J. Simpson, 2007: New, improved bulk-microphysical      c
!   schemes for studying precipitation processes in WRF. Part I:       c
!   Comparisons with other schemes. to appear on Mon. Wea. Rev.        C
!                                                                      c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      integer :: ice2,ihail,itimestep
!      integer :: ids,ide,jds,jde,kds,kde
      integer :: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte
      integer :: i,j,k,ihail2
      real :: a1,a2,ami,ami50,del,ee1,r2ice,rt0,dt,d2t,r2is,r2ig,erw,  &
              eri,esr,egr,rn10,rn17,idt
      real :: cpi,cpi2,cpi4,ga3d,ag,bg,eiw,ui50,ri50
      real, dimension (its:ite,kts:kte,jts:jte), intent(inout) :: SMLF,&
                                                                  GMLF
      real, dimension (kts:kte) :: SMLF1D,GMLF1D,rp0,pi0,pii,pr0,r00,  &
            rr0,rrq,fv0,tns,tnw,roqs,roqg,p0,RHOS,RHOG,esi1,esi2,rhois,&
            tng,cpm,MVDS,MVDH,EFDS,EFDG,EFDH,KDXR,KDXS,KDXG,KDXH,MVDR, &
            EFDR,N0R,N0S,N0G,N0H,LAMR,AVR,BVR,rho,LTK,LTK2,LQR,  &
            LQR2,LQS,LQS2,LQG,LQG2,AFAR,AFAS,AFAG,LAMG,MVDG,DSLL,LZR,  &
            LZS,LZG,GR1,GR2,GR3,GR4,GR5,GR6,GR7,GBR3,GBR4,GBR6,GBR25,  &
            GS1,GS2,GS3,GS4,GS5,GS6,GBS3,GBS4,GBS25,GG1,GG2,GG3,GG4,   &
            GBG3,GBG4,GBG25,ff,ecs,ecg,mvrc,egi,MVG0,lan1D,ncloud,LLMR,&
            LAR4
      real, dimension(ims:ime,jms:jme) :: xland
      real, dimension(ims:ime,kms:kme,jms:jme) :: ptwrf,qvwrf,qlwrf,   &
            qrwrf,qiwrf,qswrf,qgwrf,rho_mks,pi_mks,p0_mks
      real, dimension(kts:kte) :: pt,qv,qc,qr,qi,qs,qg,tair,tairc,egs, &
            rtair,dep,dd,dd1,qvs,dm,rsub1,cnd,dlt2,dlt3,zr,vr,zs,vs,vg,&
            zg,pg,qsi,ssi,esi,qsw,ssw,y1,y2,y3,y4,y5,y6,y7,r_nci,dwv,  &
            qcsink,qisink,qrsoce,qrsink,qssoce,qssink,qgsoce,scv,tca
      real, dimension(kts:kte) :: praut,pracw,dgacs,dgacw,dgaci,dgacr, &
            pgacs,wgacs,qgacw,wgaci,qgacr,pgwet,pgdry,pgaut,pracs,     &
            psacr,qsacr,pgfr,psmlt,pgmlt,psdep,piacr,pevrv,pwacs,wgacr,&
            pidep,pint,pihom,pifr,piim,pimlt,psaut,qracs,psaci,psacw,  &
            piacw,qsacw,praci,pmlts,pmltg,pssub,pgsub,nraut
!+---+-----------------------------------------------------------------+
      ihail2 = 0
      ag = 0.
      bg = 0.
      r2ig = 1.
      r2is = 1.
      if (ice2.eq.1) then
         r2ig = 0.
         r2is = 1.
      endif
      if (ice2.eq.2) then
         r2ig = 1.
         r2is = 0.
      endif
      rt0  = 1./(t0-t00)
      cpi  = 4.*ATAN(1.)
      cpi2 = cpi*cpi
      cpi4 = 0.25*cpi
      idt  = 1./dt
      DO j = jts,jte   ! j_loop
      DO i = its,ite   ! i_loop:
      ! -----------------------------------------------------------------!
      ! put 3D array data to 1D array
      DO k = kts,kte
         lan1D(k) = xland(i,j)
         r00(k) = rho_mks(i,k,j)*0.001
         rho(k) = rho_mks(i,k,j)
         p0(k)  = p0_mks(i,k,j)*10.
         pi0(k) = pi_mks(i,k,j)
         rp0(k) = 3.799052e3/p0(k)
         pii(k) = 1./pi0(k)
         pr0(k) = 1./p0(k)
         rr0(k) = 1./r00(k)
         rrq(k) = SQRT(SQRT(rr0(k)))
         fv0(k) = SQRT(rho_mks(i,2,j)*0.001/r00(k))
         pt(k)  = ptwrf(i,k,j)
         qv(k)  = qvwrf(i,k,j)
         qc(k)  = qlwrf(i,k,j)
         qr(k)  = qrwrf(i,k,j)
         qi(k)  = qiwrf(i,k,j)
         qs(k)  = qswrf(i,k,j)
         qg(k)  = qgwrf(i,k,j)
      ENDDO

      ! main k-loop
      DO k = kts,kte
         ! initializing 
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qr(k).lt.cmin1) qr(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
         if (qs(k).lt.cmin1) qs(k) = 0.
         if (qg(k).lt.cmin1) qg(k) = 0.
         tair(k) = pt(k)*pi0(k)
         tairc(k) = tair(k)-t0
         cpm(k) = cp*(1.+0.887*qv(k))
         ncloud(k) = 0.
         tnw(k) = 0.
         zr(k) = 0.
         LZR(k) = 0.
         vr(k) = 0.
         AFAR(k) = 0.
         AVR(k) = aw
         BVR(k) = bw
         tns(k) = 0.
         roqs(k) = 0.
         zs(k) = 0.
         LZS(k) = 0.
         vs(k) = 0.
         AFAS(k) = 0.
         MVG0(k) = 0.
         DSLL(k) = 1.
         tng(k) = 0.
         roqg(k) = 0.
         zg(k) = 0.
         LZG(k) = 0.
         vg(k) = 0.
         AFAG(k) = 0.
         SMLF1D(k) = 0.
         GMLF1D(k) = 0.
        
         ! for Clould 
         IF (qc(k).GE.cmin1) THEN
            if (lan1D(k).eq.1.) then     ! land
               ncloud(k) = 250.*1.E6*rho(k)
            elseif (lan1D(k).eq.2.) then     ! ocean
               ncloud(k) = 50.*1.E6*rho(k)
            else
               ncloud(k) = 100.*1.E6*rho(k)
            endif
         ENDIF
         
         ! for Rain
         IF (qr(k).GE.cmin1) THEN
            ! Spectrum chosen
            IF (SPECTRUM.EQ.0) THEN
               AFAR(k) = 0.
               tnw(k) = LOG(xnor)-LOG(10.**(8.+2.*AFAR(k)))
               dd(k) = r00(k)*qr(k)
               zr(k) = EXP((tnw(k)+LOG(cpi*roqr/dd(k)/6.)+LGAMMA(      &
                       AFAR(k)+4.))/(AFAR(k)+4.))
            ELSEIF (SPECTRUM.EQ.1) THEN
               LTK(k)  = LOG(tair(k))
               LQR(k)  = -1.*LOG(qr(k)*rho(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQR2(k) = LQR(k)*LQR(k)
               MVDR(k) = EXP(MDR1+MDR2*LTK(k)+MDR3*LQR(k)+MDR4*LTK2(k)+&
                         MDR5*LQR2(k)+MDR6*LTK(k)*LQR(k))
               EFDR(k) = EXP(EFD1+EFD2*LOG(MVDR(k))+EFD3*LOG(1.E3)+    &
                         EFD4*LOG(MVDR(k))**2.+EFD5*LOG(1.E3)**2.+     &
                         EFD6*LOG(MVDR(k))*LOG(1.E3))
               EFDR(k) = MAX(MIN(EFDR(k),MVDR(k)/KMIN**THRD),MVDR(k)/  &
                         KMAX**THRD)
               KDXR(k) = MAX(KMIN,MIN(KMAX,(MVDR(k)/EFDR(k))**3.))
               AFAR(k) = MAX(0.,(6.*KDXR(k)-3.+SQRT(8.*KDXR(k)+        &
                         1.))/(2.-2.*KDXR(k)))
               zr(k) = (AFAR(k)+3.)/EFDR(k)*1.E6
               tnw(k) = LOG(6.*qr(k)*rho(k)/cpi/1.E3)+(4.+AFAR(k))*    &
                        LOG(zr(k))-LGAMMA(AFAR(k)+4.)-LOG(10.**(8.+2.* &
                        AFAR(k)))
               zr(k) = zr(k)/1.E2
            ELSEIF (SPECTRUM.EQ.2) THEN
               MVDR(k) = EXP(8.425904+0.25*LOG(rho(k)*qr(k)))
               AFAR(k) = 0.
               zr(k) = (EXP(LGAMMA(AFAR(k)+4.)-LGAMMA(AFAR(k)+1.)))**  &
                       THRD/MVDR(k)*1.E6
               tnw(k) = LOG(6.*qr(k)*rho(k)/cpi/1.E3)+(4.+AFAR(k))*    &
                        LOG(zr(k))-LGAMMA(AFAR(k)+4.)-LOG(10.**(8.+2.* &
                        AFAR(k)))
               zr(k) = zr(k)/1.E2
            ENDIF

            LZR(k) = LOG(zr(k))
            GR4(k) = LGAMMA(AFAR(k)+4.)
            GR5(k) = LGAMMA(AFAR(k)+5.)
            GR6(k) = LGAMMA(AFAR(k)+6.)
            IF (AFAR(k).LT.0.1) THEN
               BVR(k) = bw
               AVR(k) = aw
            ELSE
               LTK(k)  = LOG(tair(k))
               LQR(k)  = -1.*LOG(qr(k)*rho(k))
               LQR2(k) = LQR(k)*LQR(k)
               BVR(k)  = MAX(0.5,MIN(2.,1.2218728-0.1281004*LTK(k)+    &
                         2.6088596E-2*LQR(k)-7.4467639E-3*LQR2(k)+     &
                         7.7592532E-4*LQR(k)*LQR2(k)-1.7056075E-5*     &
                         LQR2(k)*LQR2(k)))
               AVR(k)  = EXP(7.6004532-0.7990953*LTK(k)+1.0281818*     &
                         LQR(k)-0.16595505*LQR2(k)+1.110037E-2*LQR(k)* &
                         LQR2(k)-2.0925743E-4*LQR2(k)*LQR2(k))  
               AVR(k)  = MAX(2.E2,MIN(3.E7,AVR(k)))*10.**(2.-2.*BVR(k))
            ENDIF

            IF (VTRX.EQ.0) THEN
               GBR4(k) = LGAMMA(BVR(k)+AFAR(k)+4.)
               vr(k) = AVR(k)*EXP(GBR4(k)-GR4(k)-BVR(k)*LZR(k))*fv0(k)
            ELSEIF (VTRX.EQ.1) THEN
               LTK(k)  = LOG(tair(k))
               LQR(k)  = LOG(-1.*LOG(qr(k)*rho(k)))
               LQR2(k) = LQR(k)*LQR(k)
               vr(k)   = MIN(9.1,EXP((VTR1+VTR2*LTK(k)+VTR3*LQR(k)+    &
                         VTR4*LQR2(k))/(1.+VTR5*LTK(k)+VTR6*LQR(k)+    &
                         VTR7*LQR2(k)))/1.E4)*1.E2*fv0(k)
            ELSEIF (VTRX.EQ.2) THEN
               LAMR(k) = zr(k)*1.E2
               vr(k) = CHEN(rho(k),AFAR(k)+3.,LAMR(k))*1.E2*fv0(k)
            ELSEIF (VTRX.EQ.3) THEN
               LAMR(k) = zr(k)*1.E2
               vr(k) = MAX(0.,-26.7+2.06E6/LAMR(k)-2.045E9/LAMR(k)**   &
                       2.+9.06E11/LAMR(k)**3.)*fv0(k)! GCE-4ICE
            ENDIF
         ENDIF
         
         ! for Snow
         IF (qs(k).GE.cmin1) THEN
            ! Spectrum chosen
            IF (SPECTRUM.EQ.0) THEN
               roqs(k) = 0.1
               AFAS(k) = 0.
               tns(k) = LOG(xnos)-LOG(10.**(8.+2.*AFAS(k)))
               dd(k) = r00(k)*qs(k)
               zs(k) = EXP((tns(k)+LOG(cpi*roqs(k)/dd(k)/6.)+LGAMMA(   &
                       AFAS(k)+4.))/(AFAS(k)+4.))
            ELSEIF (SPECTRUM.EQ.1) THEN
               LTK(k) = LOG(tair(k))
               LQS(k) = -1.*LOG(qs(k)*rho(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQS2(k) = LQS(k)*LQS(k) 
               RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*  &
                         LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)* &
                         LQS(k))))
               roqs(k) = RHOS(k)/1.E3
               MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+&
                         MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
               EFDS(k) = EXP(EFD1+EFD2*LOG(MVDS(k))+EFD3*LOG(RHOS(k))+ &
                         EFD4*LOG(MVDS(k))**2.+EFD5*LOG(RHOS(k))**2.+  &
                         EFD6*LOG(MVDS(k))*LOG(RHOS(k)))
               EFDS(k) = MAX(MIN(EFDS(k),MVDS(k)/KMIN**THRD),MVDS(k)/  &
                         KMAX**THRD)
               KDXS(k) = MAX(KMIN,MIN(KMAX,(MVDS(k)/EFDS(k))**3.))
               AFAS(k) = MAX(0.,(6.*KDXS(k)-3.+SQRT(8.*KDXS(k)+        &
                         1.))/(2.-2.*KDXS(k)))
               zs(k) = (AFAS(k)+3.)/EFDS(k)*1.E6
               tns(k) = LOG(6.*qs(k)*rho(k)/cpi/RHOS(k))+(4.+AFAS(k))* &
                        LOG(zs(k))-LGAMMA(AFAS(k)+4.)-LOG(10.**(8.+2.* &
                        AFAS(k)))
               zs(k) = zs(k)/1.E2
            ELSEIF (SPECTRUM.EQ.2) THEN
               LTK(k) = LOG(tair(k))
               LQS(k) = -1.*LOG(qs(k)*rho(k))
               LTK2(k) = LTK(k)*LTK(k)
               LQS2(k) = LQS(k)*LQS(k)
               RHOS(k) = MAX(100.,MIN(910.,EXP(ROS1+ROS2*LTK(k)+ROS3*  &
                         LQS(k)+ROS4*LTK2(k)+ROS5*LQS2(k)+ROS6*LTK(k)* &
                         LQS(k))))
               roqs(k) = RHOS(k)/1.E3
               MVDS(k) = EXP(MDS1+MDS2*LTK(k)+MDS3*LQS(k)+MDS4*LTK2(k)+&
                         MDS5*LQS2(k)+MDS6*LTK(k)*LQS(k))
               AFAS(k) = 0.
               zs(k) = (EXP(LGAMMA(AFAS(k)+4.)-LGAMMA(AFAS(k)+1.)))**  &
                       THRD/MVDS(k)*1.E6
               tns(k) = LOG(6.*qs(k)*rho(k)/cpi/RHOS(k))+(4.+AFAS(k))* &
                        LOG(zs(k))-LGAMMA(AFAS(k)+4.)-LOG(10.**(8.+2.* &
                        AFAS(k)))
               zs(k) = zs(k)/1.E2
            ENDIF
            LZS(k) = LOG(zs(k))
            GS2(k) = LGAMMA(AFAS(k)+2.)
            GS4(k) = LGAMMA(AFAS(k)+4.)
            GS5(k) = LGAMMA(AFAS(k)+5.)
            GS6(k) = LGAMMA(AFAS(k)+6.)
            GBS3(k) = LGAMMA(bs+AFAS(k)+3.)
            GBS4(k) = LGAMMA(bs+AFAS(k)+4.)
            vs(k) = as*fv0(k)*EXP(GBS4(k)-GS4(k)-bs*LZS(k))
         ENDIF

         ! for Graupel
         if (qg(k).ge.cmin1) then
            MVG0(k) = EXP(8.8283+0.25*LOG(rho(k)*qg(k)))*1.E-6
            DSLL(k) = MAX(1.E-3,MIN(1.,2.*1.E-2*(EXP(MIN(20.,-tairc(k)/&
                      MAX(0.1,(1.1E4*(qc(k)+qr(k))*rho(k)-1.3E3*qg(k)* &
                      rho(k)+1.))))-1.)))
         ENDIF

         if (ihail.eq.2) then
            IF (MVG0(k).GT.DSLL(k)) THEN
               ihail2 = 1
            ELSE
               ihail2 = 0
            ENDIF
         ENDIF
         ! for Graupel
         if (qg(k).ge.cmin1) then
            if (ihail.eq.1.or.ihail2.eq.1) then
               roqg(k) = .9
               ag = SQRT(cd2*roqg(k))
               bg = .5
               IF (SPECTRUM.EQ.0) THEN
                  AFAG(k) = 0.
                  tng(k) = LOG(xnoh)-LOG(10.**(8.+2.*AFAG(k)))
                  dd(k) = r00(k)*qg(k)
                  zg(k) = EXP((tng(k)+LOG(cpi*roqg(k)/dd(k)/6.)+       &
                          LGAMMA(AFAG(k)+4.))/(AFAG(k)+4.))
               ELSEIF (SPECTRUM.EQ.1) THEN
                  LTK(k) = LOG(tair(k))
                  LQG(k) = -1.*LOG(qg(k)*rho(k))
                  LTK2(k) = LTK(k)*LTK(k)
                  LQG2(k) = LQG(k)*LQG(k)
                  MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*     &
                            LTK2(k)+MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
                  EFDH(k) = EXP(EFD1+EFD2*LOG(MVDH(k))+EFD3*LOG(9.E2)+ &
                            EFD4*LOG(MVDH(k))**2.+EFD5*LOG(9.E2)**2.+  &
                            EFD6*LOG(MVDH(k))*LOG(9.E2))
                  EFDH(k) = MAX(MIN(EFDH(k),MVDH(k)/KMIN**THRD),       &
                            MVDH(k)/KMAX**THRD)
                  KDXH(k) = MAX(KMIN,MIN(KMAX,(MVDH(k)/EFDH(k))**3.))
                  AFAG(k) = MAX(0.,(6.*KDXH(k)-3.+SQRT(8.*KDXH(k)+     &
                            1.))/(2.-2.*KDXH(k)))
                  zg(k) = (AFAG(k)+3.)/EFDH(k)*1.E6
                  tng(k) = LOG(6.*qg(k)*rho(k)/cpi/900.)+(4.+AFAG(k))* &
                           LOG(zg(k))-LGAMMA(AFAG(k)+4.)-LOG(10.**(8.+ &
                           2.*AFAG(k)))
                  zg(k) = zg(k)/1.E2
               ELSEIF (SPECTRUM.EQ.2) THEN
                 LTK(k) = LOG(tair(k))
                 LQG(k) = -1.*LOG(qg(k)*rho(k))
                 LTK2(k) = LTK(k)*LTK(k)
                 LQG2(k) = LQG(k)*LQG(k)
                 MVDH(k) = EXP(MDH1+MDH2*LTK(k)+MDH3*LQG(k)+MDH4*      &
                           LTK2(k)+MDH5*LQG2(k)+MDH6*LTK(k)*LQG(k))
                 AFAG(k) = 0.
                 zg(k) = (EXP(LGAMMA(AFAG(k)+4.)-LGAMMA(AFAG(k)+       &
                         1.)))**THRD/MVDH(k)*1.E6
                 tng(k) = LOG(6.*qg(k)*rho(k)/cpi/900.)+(4.+AFAG(k))*  &
                          LOG(zg(k))-LGAMMA(AFAG(k)+4.)-LOG(10.**(8.+  &
                          2.*AFAG(k)))
                 zg(k) = zg(k)/1.E2
               ENDIF
               LZG(k) = LOG(zg(k))
               GG4(k) = LGAMMA(AFAG(k)+4.)
               GBG3(k) = LGAMMA(bg+AFAG(k)+3.)
               GBG4(k) = LGAMMA(bg+AFAG(k)+4.)
               GBG25(k) = LGAMMA(bg*0.5+AFAG(k)+2.5)
               vg(k) = ag*SQRT(rr0(k))*EXP(GBG4(k)-GG4(k)-bg*LZG(k))
            elseif (ihail.eq.0.or.ihail2.eq.0) then
               ag = 351.2
               bg = .37
               IF (SPECTRUM.EQ.0) THEN
                  AFAG(k) = 0.
                  roqg(k) = 0.4
                  tng(k) = LOG(xnog)-LOG(10.**(8.+2.*AFAG(k)))
                  dd(k) = r00(k)*qg(k)
                  zg(k) = EXP((tng(k)+LOG(cpi*roqg(k)/dd(k)/6.)+LGAMMA(&
                          AFAG(k)+4.))/(AFAG(k)+4.))
               ELSEIF (SPECTRUM.EQ.1) THEN
                  LTK(k) = LOG(tair(k))
                  LQG(k) = -1.*LOG(qg(k)*rho(k))
                  LTK2(k) = LTK(k)*LTK(k)
                  LQG2(k) = LQG(k)*LQG(k)
                  RHOG(k) = MIN(900.,MAX(250.,EXP(-0.056433*LQG(k)+    &
                            6.608048)))
                  roqg(k) = RHOG(k)/1.E3
                  MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*     &
                            LTK2(k)+MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
                  EFDG(k) = EXP(EFD1+EFD2*LOG(MVDG(k))+EFD3*           &
                            LOG(RHOG(k))+EFD4*LOG(MVDG(k))**2.+EFD5*   &
                            LOG(RHOG(k))**2.+EFD6*LOG(MVDG(k))*        &
                            LOG(RHOG(k)))
                  EFDG(k) = MAX(MIN(EFDG(k),MVDG(k)/KMIN**THRD),       &
                            MVDG(k)/KMAX**THRD)
                  KDXG(k) = MAX(KMIN,MIN(KMAX,(MVDG(k)/EFDG(k))**3.))
                  AFAG(k) = MAX(0.,(6.*KDXG(k)-3.+SQRT(8.*KDXG(k)+     &
                            1.))/(2.-2.*KDXG(k)))
                  zg(k) = (AFAG(k)+3.)/EFDG(k)*1.E6
                  tng(k) = LOG(6.*qg(k)*rho(k)/cpi/RHOG(k))+(4.+       &
                           AFAG(k))*LOG(zg(k))-LGAMMA(AFAG(k)+4.)-     &
                           LOG(10.**(8.+2.*AFAG(k)))
                  zg(k) = zg(k)/1.E2
               ELSEIF (SPECTRUM.EQ.2) THEN
                  LTK(k) = LOG(tair(k))
                  LQG(k) = -1.*LOG(qg(k)*rho(k))
                  LTK2(k) = LTK(k)*LTK(k)
                  LQG2(k) = LQG(k)*LQG(k)
                  RHOG(k) = MIN(900.,MAX(250.,EXP(-0.056433*LQG(k)+    &
                            6.608048)))
                  roqg(k) = RHOG(k)/1.E3
                  MVDG(k) = EXP(MDG1+MDG2*LTK(k)+MDG3*LQG(k)+MDG4*     &
                            LTK2(k)+MDG5*LQG2(k)+MDG6*LTK(k)*LQG(k))
                  AFAG(k) = 0.
                  zg(k) = (EXP(LGAMMA(AFAG(k)+4.)-LGAMMA(AFAG(k)+      &
                          1.)))**THRD/MVDG(k)*1.E6
                  tng(k) = LOG(6.*qg(k)*rho(k)/cpi/RHOG(k))+(4.+       &
                           AFAG(k))*LOG(zg(k))-LGAMMA(AFAG(k)+4.)-     &
                           LOG(10.**(8.+2.*AFAG(k)))
                  zg(k) = zg(k)/1.E2
               ENDIF
               LZG(k) = LOG(zg(k))
               GG4(k) = LGAMMA(AFAG(k)+4.)
               GBG3(k) = LGAMMA(bg+AFAG(k)+3.)
               GBG4(k) = LGAMMA(bg+AFAG(k)+4.)
               GBG25(k) = LGAMMA(bg*0.5+AFAG(k)+2.5)
               vg(k) = ag*fv0(k)*EXP(GBG4(k)-GG4(k)-bg*LZG(k))
            endif          ! ihail
            GG1(k) = LGAMMA(AFAG(k)+1.)
            GG2(k) = LGAMMA(AFAG(k)+2.)
            GG3(k) = LGAMMA(AFAG(k)+3.)
         endif          ! qg min

!     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
!     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
!     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
!     ***   Y2 : KINETIC VISCOSITY (V)
         y1(k) = c149*tair(k)**1.5/(tair(k)+120.)
         dwv(k) = c879*pr0(k)*tair(k)**1.81
         tca(k) = c141*y1(k)
         scv(k) = 1./((rr0(k)*y1(k))**.1667*dwv(k)**.333)
         praut(k) = 0.
         pracw(k) = 0.
         psaut(k) = 0.
         psaci(k) = 0.
         praci(k) = 0.
         piacr(k) = 0.
         psacw(k) = 0.
         piacw(k) = 0.
         qsacw(k) = 0.
         pracs(k) = 0.
         psacr(k) = 0.
         qsacr(k) = 0.
         dgaci(k) = 0.
         wgaci(k) = 0.
         pgwet(k) = 0.
         pgdry(k) = 0.
         pgaut(k) = 0.
         pgfr(k) = 0.
         dgacw(k) = 0.
         pgacs(k) = 0.
         dgacs(k) = 0.
         wgacs(k) = 0.
         qgacw(k) = 0.
         dgacr(k) = 0.
         qgacr(k) = 0.
         piim(k) = 0.
         SMLF1D(k)  = 0.
         GMLF1D(k)  = 0.
!-------------------------------------------------------------------------------
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
!* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
         if (qc(k).gt.1.E-5) then
            mvrc(k) = MIN(5.E-5,(6.E-3*qc(k)/ncloud(k)/cpi)**THRD/2.)
            nraut(k) = MIN(ncloud(k),ncloud(k)*ncloud(k)*EXP(-40.731+  &
                       5.372E5*mvrc(k)-2.0139E-5/mvrc(k)))
            praut(k) = nraut(k)*EXP(-21.37+1.9899E9*qc(k)/ncloud(k))
            praut(k) = MIN(qc(k)*idt,praut(k))
         endif
         if (qc(k).ge.cmin2.and.qr(k).ge.cmin2) then
            erw = 1.
            GBR3(k) = LGAMMA(BVR(k)+AFAR(k)+3.)
            y1(k) = EXP(tnw(k)+GBR3(k)-(BVR(k)+AFAR(k)+3.)*LZR(k))
            pracw(k) = MIN(qc(k)*idt,cpi4*erw*AVR(k)*fv0(k)*qc(k)*y1(k))
         endif
!*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
!*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
!*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
!*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
!*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
         if (tair(k).lt.t0) then
            esi1(k) = MIN(1.,0.05*EXP(0.1*tairc(k)))
            if (qi(k).gt.2.E-3) then
               psaut(k) = r2is*MAX(1.E-3*esi1(k)*(qi(k)-2.E-3),0.)
               psaut(k) = MIN(qi(k)*idt,psaut(k))
            endif
!*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
!* 18 * PGFR : FREEZING OF QR TO QG                               **18**
            if (qr(k).ge.cmin2) then
               GR7(k) = LGAMMA(AFAR(k)+7.)
               y1(k) = cpi2*1.E-4*roqr*rr0(k)/36.
               y2(k) = EXP(0.66*(t0-tair(k)))-1.
               y3(k) = EXP(tnw(k)+GR7(k)-(AFAR(k)+7.)*LZR(k))
               pgfr(k) = MIN(qr(k)*idt,r2ig*y1(k)*y2(k)*y3(k))
            endif
            if (qs(k).gt.2.E-3) then
               y1(k) = EXP(0.09*tairc(k))
               pgaut(k) = r2is*MAX(1.E-3*y1(k)*(qs(k)-2.E-3),0.)
            endif
            if (qi(k).ge.cmin2.and.qs(k).ge.cmin2) then
               IF (SPECTRUM.EQ.0) THEN
                  esi2(k) = 0.
               ELSEIF (SPECTRUM.EQ.1.OR.SPECTRUM.EQ.2) THEN
                  rhois(k) = (qi(k)*0.9+qs(k)*roqs(k))/(qi(k)+qs(k))
                  esi2(k) = 1.-SQRT(MIN(MAX(rhois(k),0.1),0.91)/0.91)
               ENDIF
               esi(k) = MIN(MAX(esi1(k),esi2(k),0.),1.)
               y1(k) = EXP(tns(k)+GBS3(k)-(bs+AFAS(k)+3.)*LZS(k))
               psaci(k) = r2is*cpi4*as*fv0(k)*esi(k)*qi(k)*y1(k)
               psaci(k) = MIN(qi(k)*idt,psaci(k))
            endif
!            if (qc(k).ge.cmin2.and.qi(k).ge.cmin2) then
!               eiw = 0.5
!               ui50 = 100.
!               ri50 = 5.e-3
!               ami50 = 3.76e-8
!               piacw(k) = MIN(qc(k)/dt,r2is*cpi*eiw*ui50*ri50*ri50*    &
!                          qc(k)*r00(k)*qi(k)*1.e-9/ami50)
!            endif
            if (qc(k).ge.cmin2.and.qs(k).ge.cmin2) then
               ecs(k) = EXP(0.1*tairc(k))
               y1(k) = EXP(tns(k)+GBS3(k)-(bs+AFAS(k)+3.)*LZS(k))
               psacw(k) = r2is*cpi4*ecs(k)*as*fv0(k)*qc(k)*y1(k)
               psacw(k) = MIN(qc(k)*idt,psacw(k))
            endif
            if (qi(k).ge.cmin2.and.qr(k).ge.cmin2) then
               eri = 1.
               ami = 1./(24.*6.e-9)
               GBR3(k) = LGAMMA(BVR(k)+AFAR(k)+3.)
               GBR6(k) = LGAMMA(BVR(k)+AFAR(k)+6.)
               y1(k) = EXP(tnw(k)+GBR3(k)-(BVR(k)+AFAR(k)+3.)*LZR(k))
               y2(k) = EXP(tnw(k)+GBR6(k)-(BVR(k)+AFAR(k)+6.)*LZR(k))
               praci(k) = r2is*cpi4*eri*AVR(k)*fv0(k)*qi(k)*y1(k)
               piacr(k) = r2is*cpi2*eri*AVR(k)*roqr*ami*fv0(k)*qi(k)*  &
                          y2(k)
               praci(k) = MIN(qi(k)*idt,praci(k))
               piacr(k) = MIN(qr(k)*idt,piacr(k))
            endif
!*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
!*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
            if (qr(k).ge.cmin2.and.qs(k).ge.cmin2) then
               esr = 1.
               GR1(k) = LGAMMA(AFAR(k)+1.)
               GR2(k) = LGAMMA(AFAR(k)+2.)
               GR3(k) = LGAMMA(AFAR(k)+3.)
               GS1(k) = LGAMMA(AFAS(k)+1.)
               GS3(k) = LGAMMA(AFAS(k)+3.)
               y1(k) = ABS(vr(k)-vs(k))
               y2(k) = EXP(tns(k)+tnw(k)+GS6(k)+GR1(k)-(AFAS(k)+6.)*   &
                       LZS(k)-(AFAR(k)+1.)*LZR(k))
               y3(k) = EXP(tns(k)+tnw(k)+GS5(k)+GR2(k)-(AFAS(k)+5.)*   &
                       LZS(k)-(AFAR(k)+2.)*LZR(k))
               y4(k) = EXP(tns(k)+tnw(k)+GS4(k)+GR3(k)-(AFAS(k)+4.)*   &
                       LZS(k)-(AFAR(k)+3.)*LZR(k))
               y5(k) = EXP(tns(k)+tnw(k)+GR6(k)+GS1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAS(k)+1.)*LZS(k))
               y6(k) = EXP(tns(k)+tnw(k)+GR5(k)+GS2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAS(k)+2.)*LZS(k))
               y7(k) = EXP(tns(k)+tnw(k)+GR4(k)+GS3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAS(k)+3.)*LZS(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               dd1(k) = y1(k)*(y5(k)+2.*y6(k)+y7(k))/24.
               pracs(k) = r2ig*r2is*cpi2*esr*roqs(k)*rr0(k)*dd(k)
               psacr(k) = r2is*r2ig*cpi2*esr*roqr*rr0(k)*dd1(k)
               pracs(k) = MIN(qs(k)*idt,pracs(k))
               psacr(k) = MIN(qr(k)*idt,psacr(k))
            endif
!*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
!* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
!* 17 * PGWET : WET GROWTH OF QG                                  **17**
            if (qi(k).ge.cmin2.and.qg(k).ge.cmin2) then
!               egi = 0.05
               egi(k) = MIN(1.,0.01*EXP(0.1*tairc(k)))
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y1(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  y2(k) = cpi4*egi(k)*ag*SQRT(rr0(k))
                  dgaci(k) = MIN(qi(k)*idt,r2ig*y2(k)*qi(k)*y1(k))
                  wgaci(k) = dgaci(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y3(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  dgaci(k) = 0.
                  wgaci(k) = r2ig*cpi4*egi(k)*ag*fv0(k)*qi(k)*y3(k)
                  wgaci(k) = MIN(qi(k)*idt,wgaci(k))
               endif
            endif
!TTT***** QG=QG+MIN(PGDRY,PGWET)
!*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
!* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
!* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
            if (qg(k).ge.cmin2.and.qs(k).ge.cmin2) then
               if (qc(k)+qr(k).lt.1.E-4) then
                  ee1 = 0.01
               else
                  ee1 = 1.
               endif
               egs(k) = ee1*EXP(0.09*tairc(k))
               y1(k) = ABS(vg(k)-vs(k))
               y2(k) = EXP(tns(k)+tng(k)+GS6(k)+GG1(k)-(AFAS(k)+6.)*   &
                       LZS(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tns(k)+tng(k)+GS5(k)+GG2(k)-(AFAS(k)+5.)*   &
                       LZS(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tns(k)+tng(k)+GS4(k)+GG3(k)-(AFAS(k)+4.)*   &
                       LZS(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               if (ihail.eq.1.or.ihail2.eq.1) then  ! hail
                  dgacs(k) = MIN(qs(k)*idt,r2ig*r2is*cpi2*roqs(k)*     &
                             rr0(k)*egs(k)*dd(k))
               elseif (ihail.eq.0.or.ihail2.eq.0) then ! graupel
                  dgacs(k) = 0.
               endif
               wgacs(k) = r2ig*r2is*cpi2*roqs(k)*rr0(k)*dd(k)
               wgacs(k) = MIN(qs(k)*idt,wgacs(k))
            endif
            if (qc(k).ge.cmin2.and.qg(k).ge.cmin2) then
               ecg(k) = EXP(0.05*tairc(k))
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y1(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  dgacw(k) = r2ig*cpi4*ecg(k)*ag*SQRT(rr0(k))*qc(k)*   &
                             y1(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y2(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  dgacw(k) = r2ig*cpi4*ecg(k)*ag*fv0(k)*qc(k)*y2(k)
               endif
               dgacw(k) = MIN(qc(k)*idt,dgacw(k))
            endif
            if (qr(k).ge.cmin2.and.qg(k).ge.cmin2) then
               egr = 1.
               y1(k) = ABS(vg(k)-vr(k))
               y2(k) = EXP(tng(k)+tnw(k)+GR6(k)+GG1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tng(k)+tnw(k)+GR5(k)+GG2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tng(k)+tnw(k)+GR4(k)+GG3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               dgacr(k) = MIN(qr(k)*idt,r2ig*cpi2*egr*roqr*rr0(k)*dd(k))
            endif
            if (qg(k).ge.cmin2) then
               y1(k) = 1./(alf+cw*tairc(k))
               y2(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
               y3(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*LZG(k))
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*rrq(k)*scv(k)*y3(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*SQRT(fv0(k))*scv(k)*&
                          y3(k)
               endif
               rn17 = cw-ci
               y7(k) = r00(k)*alv*dwv(k)*(rp0(k)-qv(k))-tca(k)*tairc(k)
               pgwet(k) = r2ig*MAX(y1(k)*(2.*cpi*rr0(k)*y7(k)*y6(k)+(  &
                          wgaci(k)+wgacs(k))*(alf+rn17*tairc(k))),0.)
            endif
         else                  ! t >= t0
            if (qc(k).ge.cmin2.and.qs(k).ge.cmin2) then
               y1(k) = EXP(tns(k)+GBS3(k)-(bs+AFAS(k)+3.)*LZS(k))
               qsacw(k) = r2is*cpi4*1.*as*fv0(k)*qc(k)*y1(k)
               qsacw(k) = MIN(qc(k)*idt,qsacw(k))
            endif
            if (qg(k).ge.cmin2.and.qs(k).ge.cmin2) then
               y1(k) = ABS(vg(k)-vs(k))
               y2(k) = EXP(tns(k)+tng(k)+GS6(k)+GG1(k)-(AFAS(k)+6.)*   &
                       LZS(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tns(k)+tng(k)+GS5(k)+GG2(k)-(AFAS(k)+5.)*   &
                       LZS(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tns(k)+tng(k)+GS4(k)+GG3(k)-(AFAS(k)+4.)*   &
                       LZS(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               pgacs(k) = r2ig*r2is*cpi2*roqs(k)*rr0(k)*1.*dd(k)
               pgacs(k) = MIN(qs(k)*idt,pgacs(k))
            endif
            if (qr(k).ge.cmin2.and.qs(k).ge.cmin2) then
               GS1(k) = LGAMMA(AFAS(k)+1.)
               GS3(k) = LGAMMA(AFAS(k)+3.)
               y1(k) = ABS(vr(k)-vs(k))
               y5(k) = EXP(tns(k)+tnw(k)+GR6(k)+GS1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAS(k)+1.)*LZS(k))
               y6(k) = EXP(tns(k)+tnw(k)+GR5(k)+GS2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAS(k)+2.)*LZS(k))
               y7(k) = EXP(tns(k)+tnw(k)+GR4(k)+GS3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAS(k)+3.)*LZS(k))
               dd(k) = y1(k)*(y5(k)+2.*y6(k)+y7(k))/24.
               qsacr(k) = r2is*r2ig*cpi2*1.*roqr*rr0(k)*dd(k)
               qsacr(k) = MIN(qr(k)*idt,qsacr(k))
            endif
            if (qc(k).ge.cmin2.and.qg(k).ge.cmin2) then
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y1(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  qgacw(k) = r2ig*cpi4*1.*ag*SQRT(rr0(k))*qc(k)*y1(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y2(k) = EXP(tng(k)+GBG3(k)-(bg+AFAG(k)+3.)*LZG(k))
                  qgacw(k) = r2ig*cpi4*1.*ag*fv0(k)*qc(k)*y2(k)
               endif
               qgacw(k) = MIN(qc(k)*idt,qgacw(k))
            endif
            if (qr(k).ge.cmin2.and.qg(k).ge.cmin2) then
               y1(k) = ABS(vg(k)-vr(k))
               y2(k) = EXP(tng(k)+tnw(k)+GR6(k)+GG1(k)-(AFAR(k)+6.)*   &
                       LZR(k)-(AFAG(k)+1.)*LZG(k))
               y3(k) = EXP(tng(k)+tnw(k)+GR5(k)+GG2(k)-(AFAR(k)+5.)*   &
                       LZR(k)-(AFAG(k)+2.)*LZG(k))
               y4(k) = EXP(tng(k)+tnw(k)+GR4(k)+GG3(k)-(AFAR(k)+4.)*   &
                       LZR(k)-(AFAG(k)+3.)*LZG(k))
               dd(k) = y1(k)*(y2(k)+2.*y3(k)+y4(k))/24.
               qgacr(k) = MIN(qr(k)*idt,r2ig*cpi2*1.*roqr*rr0(k)*dd(k))
            endif
         endif    ! tair
!********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
         qcsink(k) = (psacw(k)+praut(k)+pracw(k)+dgacw(k)+qsacw(k)+    &
                     qgacw(k)+piacw(k))*dt
         qc(k) = MAX(0.,qc(k)-qcsink(k))
         if (qc(k).lt.0.) then
            a1 = 1.
            if (qcsink(k).ne.0.) a1 = qc(k)/qcsink(k)+1.
            psacw(k) = psacw(k)*a1
            piacw(k) = piacw(k)*a1
            praut(k) = praut(k)*a1
            pracw(k) = pracw(k)*a1
            dgacw(k) = dgacw(k)*a1
            qsacw(k) = qsacw(k)*a1
            qgacw(k) = qgacw(k)*a1
            qc(k) = 0.
         endif
!******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
         wgacr(k) = pgwet(k)-dgacw(k)-wgaci(k)-wgacs(k)
         wgacr(k) = MIN(qr(k)*idt,wgacr(k))
         wgacr(k) = MAX(-qg(k)*idt,wgacr(k))
         pgdry(k) = dgacw(k)+dgaci(k)+dgacr(k)+dgacs(k)
         if (pgwet(k).ge.pgdry(k)) then
            wgacr(k) = 0.
            wgaci(k) = 0.
            wgacs(k) = 0.
         else
            dgacr(k) = 0.
            dgaci(k) = 0.
            dgacs(k) = 0.
         endif
!********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
         qisink(k) = (psaut(k)+psaci(k)+praci(k)+dgaci(k)+wgaci(k))*dt
         qi(k) = MAX(0.,qi(k)-qisink(k))
         if (qi(k).lt.0.) then
            a2 = 1.
            if (qisink(k).ne.0.) a2 = qi(k)/qisink(k)+1.
            psaut(k) = psaut(k)*a2
            psaci(k) = psaci(k)*a2
            praci(k) = praci(k)*a2
            dgaci(k) = dgaci(k)*a2
            wgaci(k) = wgaci(k)*a2
            qi(k) = 0.
         endif
         dlt3(k) = 0.
         dlt2(k) = 0.
         if (tair(k).lt.t0) then
            if (qr(k).lt.1.E-4) then
               dlt3(k) = 1.
               dlt2(k) = 1.
            endif
            if (qs(k).ge.1.E-4) then
               dlt2(k) = 0.
            endif
         endif
         if (ice2.eq.1) then
            dlt3(k) = 1.
            dlt2(k) = 1.
         endif
!********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
         del = 0.
         if (wgacr(k).lt.0.) del = 1.
         qrsink(k) = (piacr(k)+dgacr(k)+(1.-del)*wgacr(k)+psacr(k)+    &
                     pgfr(k))*dt
         qrsoce(k) = (praut(k)+pracw(k))*dt
         qr(k) = MAX(0.,qr(k)+qrsoce(k)-qrsink(k)-del*wgacr(k)*dt)
         if (qr(k).lt.0.) then
            a1 = 1.
            if (qrsink(k).ne.0.) a1 = qr(k)/qrsink(k)+1.
            piacr(k) = piacr(k)*a1
            dgacr(k) = dgacr(k)*a1
            if (wgacr(k).gt.0.) wgacr(k) = wgacr(k)*a1
            psacr(k) = psacr(k)*a1
            pgfr(k) = pgfr(k)*a1
            qr(k) = 0.
         endif
!********   HANDLING THE NEGATIVE SNOW (QS) ***************************
         qssoce(k) = (psaut(k)+psaci(k)+psacw(k)+dlt3(k)*praci(k)+     &
                     dlt3(k)*piacr(k)+dlt2(k)*psacr(k))*dt
         qssink(k) = (pgacs(k)+dgacs(k)+wgacs(k)+pgaut(k)+(1.-dlt2(k))*&
                     pracs(k))*dt
         qs(k) = MAX(0.,qs(k)+qssoce(k)-qssink(k))
         if (qs(k).lt.0.) then
            a2 = 1.
            if (qssink(k).ne.0.) a2 = qs(k)/qssink(k)+1.
            pgacs(k) = pgacs(k)*a2
            dgacs(k) = dgacs(k)*a2
            wgacs(k) = wgacs(k)*a2
            pgaut(k) = pgaut(k)*a2
            pracs(k) = pracs(k)*a2
            qssink(k) = qssink(k)*a2
            qs(k) = 0.
         endif
         if ((psacw(k)+dgacw(k)+psacr(k)+dgacr(k)).ge.cmin2) then
            if (tair(k).ge.265.16.and.tair(k).le.270.16) then
               ff(k) = 0.
               if (tair(k).ge.265.16.and.tair(k).lt.268.16)            &
               ff(k) = MAX(0.,3.5E8*(tairc(k)+8.)*THRD)
               if (tair(k).ge.268.16.and.tair(k).le.270.16)            &
               ff(k) = MAX(0.,3.5E8*(-3.-tairc(k))*0.5)
               piim(k) = ff(k)*(psacw(k)+dgacw(k)+psacr(k)+dgacr(k))*  &
                         1.E-12            ! mio = 1.e-12
            endif
         endif
!********   HANDLING THE NEGATIVE GRAUPEL (QG) *************************
         qgsoce(k) = ((1.-dlt3(k))*praci(k)+dgaci(k)+wgaci(k)+         &
                     dgacw(k)+(1.-dlt3(k))*piacr(k)+dgacr(k)+pgfr(k)+  &
                     wgacr(k)+(1.-dlt2(k))*psacr(k))*dt
         y2(k) = (psacw(k)+dgacw(k)+piacr(k)+dgacr(k)+wgacr(k)+        &
                 psacr(k)+pgfr(k))*dt
         pt(k) = pt(k)+y2(k)*alf*pii(k)/cpm(k)
         qi(k) = MAX(0.,qi(k)+piim(k)+piacw(k)*dt)
         qg(k) = MAX(0.,qg(k)+qgsoce(k)+qssink(k))
!-------------------------------------------------------------------------
!* 11 * PSMLT : MELTING OF QS                                     **11**
!* 19 * PGMLT : MELTING OF QG TO QR                               **19**
         psmlt(k) = 0.
         pgmlt(k) = 0.
         tair(k) = pt(k)*pi0(k)
         tairc(k) = tair(k)-t0
         y1(k) = tca(k)*tairc(k)-r00(k)*alv*dwv(k)*(rp0(k)-qv(k))
         if (tair(k).ge.t0.and.qs(k).ge.cmin1) then
            GBS25(k) = LGAMMA(bs*0.5+AFAS(k)+2.5)
            y2(k) = EXP(tns(k)+GS2(k)-(AFAS(k)+2.)*LZS(k))
            y3(k) = EXP(tns(k)+GBS25(k)-(0.5*bs+AFAS(k)+2.5)*LZS(k))
            y4(k) = 0.78*y2(k)+0.31*SQRT(as)*SQRT(fv0(k))*scv(k)*y3(k)
            dd(k) = 2.*cpi*dt/alf*rr0(k)*y1(k)*y4(k)+cw/alf*dt*        &
                    tairc(k)*(qsacw(k)+qsacr(k))
            psmlt(k) = r2is*MAX(0.,MIN(dd(k),qs(k)))
            SMLF1D(k) = MIN(1.,MAX(0.,psmlt(k)/qs(k)))
         endif
         if (tair(k).ge.t0.and.qg(k).ge.cmin1) then
            if (ihail.eq.1.or.ihail2.eq.1) then
               y2(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
               y3(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*LZG(k))
               y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*rrq(k)*scv(k)*y3(k)
            elseif (ihail.eq.0.or.ihail2.eq.0) then
               y4(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
               y5(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*LZG(k))
               y6(k) = 0.78*y4(k)+0.31*SQRT(ag)*SQRT(fv0(k))*scv(k)*   &
                       y5(k)
            endif
            dd1(k) = 2.*cpi/alf*rr0(k)*dt*y1(k)*y6(k)+cw/alf*dt*       &
                     tairc(k)*(qgacw(k)+qgacr(k))
            pgmlt(k) = r2ig*MAX(0.,MIN(dd1(k),qg(k)))
            GMLF1D(k) = MIN(1.,MAX(0.,pgmlt(k)/qg(k)))
         endif
         pt(k) = pt(k)-(psmlt(k)+pgmlt(k))*alf*pii(k)/cpm(k)
         qr(k) = MAX(0.,qr(k)+psmlt(k)+pgmlt(k))
         qs(k) = MAX(0.,qs(k)-psmlt(k))
         qg(k) = MAX(0.,qg(k)-pgmlt(k))
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
!* 25 * PIFR : IMMERSION-FREEZING OF QC TO QI ( T0-5 > T)
!* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
         pihom(k) = 0.
         pimlt(k) = 0.
         pifr(k) = 0.
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
         tair(k) = pt(k)*pi0(k)
         if (tair(k).le.t00.and.qc(k).ge.cmin1) then
            pihom(k) = qc(k)
         else
            pihom(k) = 0.
         endif
         if (tair(k).ge.t0.and.qi(k).ge.cmin1) then
            pimlt(k) = qi(k)
         else
            pimlt(k) = 0.
         endif
         if (tair(k).le.(t0-4.).and.qc(k).ge.cmin1) then
            pifr(k) = MIN(qc(k)*idt,1.E-4*EXP(0.66*(t0-tair(k))-1.)*   &
                      qc(k)*qc(k)*1.E6/ncloud(k))
         endif
         y1(k) = pihom(k)-pimlt(k)+pifr(k)
         pt(k) = pt(k)+y1(k)*alf*pii(k)/cpm(k)
         qc(k) = MAX(0.,qc(k)-y1(k))
         qi(k) = MAX(0.,qi(k)+y1(k))
!     CALCULATION OF PINT USES DIFFERENT VALUES OF THE INTERCEPT AND SLOPE FOR
!     THE FLETCHER EQUATION. ALSO, ONLY INITIATE MORE ICE IF THE NEW NUMBER
!     CONCENTRATION EXCEEDS THAT ALREADY PRESENT.
!* 31 * pint  : initiation of qi                                  **31**
!* 32 * pidep : deposition of qi                                  **32**
         pint(k) = 0.
         pidep(k) = 0.
         tair(k) = pt(k)*pi0(k)
         if (tair(k).lt.t0) then
            if (qi(k).lt.cmin1) qi(k) = 0.
            tairc(k) = tair(k)-t0
            rtair(k) = 1./(tair(k)-c76)
            y2(k) = EXP(c218-c580*rtair(k))
            qsi(k) = rp0(k)*y2(k)
            esi(k) = c610*y2(k)
            ssi(k) = qv(k)/qsi(k)-1.
            r_nci(k) = 1.E-4
            if (ssi(k).ge.1.E-5) then
              r_nci(k) = MIN(1.E-3*EXP(12.96*ssi(k)-0.639),1.)
            endif
            if (qi(k).ge.cmin1) then
               y1(k) = 1./tair(k)
               rn10 = 8.314E7/(.226*18.016)
               dd(k) = y1(k)*(als*als/rw*y1(k)-alv/2.43E3)+rn10*       &
                       tair(k)/esi(k)
               pidep(k) = MAX(4.*51.545E-4*dt*SQRT(rr0(k))*ssi(k)*     &
                          r_nci(k)*SQRT(qi(k))/dd(k),0.)
            endif
            ami50 = 3.76E-8
            dd(k) = MAX(1.E-9*r_nci(k)/r00(k)-qi(k)*1.E-9/ami50,0.)
            dm(k) = MAX((qv(k)-qsi(k)),0.)
            rsub1(k) = c580*asc*qsi(k)*rtair(k)*rtair(k)
            dep(k) = dm(k)/(1.+rsub1(k))
            pint(k) = MAX(MIN(dd(k),dm(k)),0.)
            pint(k) = MIN(pint(k)+pidep(k),dep(k))
            pt(k) = pt(k)+pint(k)*als*pii(k)/cpm(k)
            qv(k) = MAX(0.,qv(k)-pint(k))
            qi(k) = MAX(0.,qi(k)+pint(k))
         endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
!* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
         psdep(k) = 0.
         pssub(k) = 0.
         pgsub(k) = 0.
         tair(k) = pt(k)*pi0(k)
         if (tair(k).lt.t0) then
            rtair(k) = 1./(tair(k)-c76)
            qsi(k) = rp0(k)*EXP(c218-c580*rtair(k))
            ssi(k) = qv(k)/qsi(k)-1.
            y1(k) = als*als/rw*r00(k)/(tca(k)*tair(k)*tair(k))+        &
                    1./(dwv(k)*qsi(k))
            if (qs(k).ge.cmin1) then
               GBS25(k) = LGAMMA(bs*0.5+AFAS(k)+2.5)
               y2(k) = EXP(tns(k)+GS2(k)-(AFAS(k)+2.)*LZS(k))
               y3(k) = EXP(tns(k)+GBS25(k)-(0.5*bs+AFAS(k)+2.5)*LZS(k))
               y4(k) = 0.78*y2(k)+0.31*SQRT(as)*SQRT(fv0(k))*scv(k)*   &
                       y3(k)
               psdep(k) = 2.*cpi*dt*ssi(k)*y4(k)/y1(k)
               pssub(k) = psdep(k)
               psdep(k) = r2is*MAX(psdep(k),0.)
               pssub(k) = r2is*MAX(-qs(k),MIN(pssub(k),0.))
            endif
            if (qg(k).ge.cmin1) then
               if (ihail.eq.1.or.ihail2.eq.1) then
                  y2(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
                  y3(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*    &
                          LZG(k))
                  y6(k) = 0.78*y2(k)+0.31*SQRT(ag)*rrq(k)*scv(k)*y3(k)
               elseif (ihail.eq.0.or.ihail2.eq.0) then
                  y4(k) = EXP(tng(k)+GG2(k)-(AFAG(k)+2.)*LZG(k))
                  y5(k) = EXP(tng(k)+GBG25(k)-(0.5*bg+AFAG(k)+2.5)*    &
                          LZG(k))
                  y6(k) = 0.78*y4(k)+0.31*SQRT(ag)*SQRT(fv0(k))*       &
                          scv(k)*y5(k)
               endif
               pgsub(k) = -2.*r2ig*cpi*dt*ssi(k)*y6(k)/y1(k)
            endif
            dm(k) = qv(k)-qsi(k)
            rsub1(k) = c580*asc*qsi(k)*rtair(k)*rtair(k)
!     ********   DEPOSITION OR SUBLIMATION OF QS  **********************
            y1(k) = dm(k)/(1.+rsub1(k))
            psdep(k) = r2is*MIN(psdep(k),MAX(y1(k),0.))
            y2(k) = MIN(y1(k),0.)
            pssub(k) = r2is*MAX(pssub(k),y2(k))
!     ********   SUBLIMATION OF QG   ***********************************
            dd(k) = MAX((-y2(k)-qs(k)),0.)
            pgsub(k) = r2ig*MIN(dd(k),qg(k),MAX(pgsub(k),0.))
            pt(k) = pt(k)+(psdep(k)+pssub(k)-pgsub(k))*als*pii(k)/cpm(k)
            qv(k) = MAX(0.,qv(k)+pgsub(k)-pssub(k)-psdep(k))
            qs(k) = MAX(0.,qs(k)+psdep(k)+pssub(k))
            qg(k) = MAX(0.,qg(k)-pgsub(k))
         endif
!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
         pevrv(k) = 0.
         tair(k) = pt(k)*pi0(k)
         rtair(k) = 1./(tair(k)-c358)
         qsw(k) = rp0(k)*EXP(c172-c409*rtair(k))
         ssw(k) = qv(k)/qsw(k)-1.
         if (qr(k).ge.cmin1.and.qv(k).lt.qsw(k)) then
            dm(k) = qv(k)-qsw(k)
            rsub1(k) = c409*avc*qsw(k)*rtair(k)*rtair(k)
            dd1(k) = MAX(-dm(k)/(1.+rsub1(k)),0.)
            GR2(k) = LGAMMA(AFAR(k)+2.)
            GBR25(k) = LGAMMA(BVR(k)*0.5+AFAR(k)+2.5)
            y3(k) = EXP(tnw(k)+GR2(k)-(AFAR(k)+2.)*LZR(k))
            y4(k) = EXP(tnw(k)+GBR25(k)-(BVR(k)*0.5+AFAR(k)+2.5)*LZR(k))
            y5(k) = alv*alv/rw*r00(k)
            y1(k) = 0.78*y3(k)+0.31*SQRT(AVR(k))*SQRT(fv0(k))*scv(k)*  &
                    y4(k)
            y2(k) = y5(k)/(tca(k)*tair(k)*tair(k))+1./(dwv(k)*qsw(k))
            pevrv(k) = -2.*cpi*dt*ssw(k)*y1(k)/y2(k)
            pevrv(k) = MIN(dd1(k),qr(k),MAX(pevrv(k),0.))
            pt(k) = pt(k)-pevrv(k)*alv*pii(k)/cpm(k)
            qv(k) = MAX(0.,qv(k)+pevrv(k))
            qr(k) = MAX(0.,qr(k)-pevrv(k))
         endif
!*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
         tair(k) = pt(k)*pi0(k)
         cnd(k) = rt0*(tair(k)-t00)
         dep(k) = rt0*(t0-tair(k))
         y1(k) = 1./(tair(k)-c358)
         y2(k) = 1./(tair(k)-c76)
         qsw(k) = rp0(k)*EXP(c172-c409*y1(k))
         qsi(k) = rp0(k)*EXP(c218-c580*y2(k))
         dd(k) = c409*pi0(k)*y1(k)*y1(k)
         dd1(k) = c580*pi0(k)*y2(k)*y2(k)
         if (qc(k).lt.cmin1) qc(k) = cmin1
         if (qi(k).lt.cmin1) qi(k) = cmin1
         if (tair(k).ge.t0) then
            dep(k) = 0.
            cnd(k) = 1.
            qi(k) = 0.
         endif
         if (tair(k).lt.t00) then
            cnd(k) = 0.
            dep(k) = 1.
            qc(k) = 0.
         endif
         y5(k) = cnd(k)*alv*pii(k)/cpm(k)+dep(k)*als*pii(k)/cpm(k)
         y1(k) = qc(k)*qsw(k)/(qc(k)+qi(k))
         y2(k) = qi(k)*qsi(k)/(qc(k)+qi(k))
         y4(k) = dd(k)*y1(k)+dd1(k)*y2(k)
         qvs(k) = y1(k)+y2(k)
         rsub1(k) = (qv(k)-qvs(k))/(1.+y4(k)*y5(k))
         cnd(k) = cnd(k)*rsub1(k)
         dep(k) = dep(k)*rsub1(k)
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
!c    ******   condensation or evaporation of qc  ******
         cnd(k) = MAX(-qc(k),cnd(k))
!c    ******   deposition or sublimation of qi    ******
         dep(k) = MAX(-qi(k),dep(k))
         pt(k) = pt(k)+cnd(k)*alv*pii(k)/cpm(k)+dep(k)*als*pii(k)/cpm(k)
         qv(k) = MAX(0.,qv(k)-cnd(k)-dep(k))
         qc(k) = MAX(0.,qc(k)+cnd(k))
         qi(k) = MAX(0.,qi(k)+dep(k))
!=========================================================================
         if (qc(k).lt.cmin1) qc(k) = 0.
         if (qr(k).lt.cmin1) qr(k) = 0.
         if (qi(k).lt.cmin1) qi(k) = 0.
         if (qs(k).lt.cmin1) qs(k) = 0.
         if (qg(k).lt.cmin1) qg(k) = 0.
      ENDDO

      ! put back to 3D array
      DO k = kts,kte
         ptwrf(i,k,j) = pt(k)
         qvwrf(i,k,j) = qv(k)
         qlwrf(i,k,j) = qc(k)
         qrwrf(i,k,j) = qr(k)
         qiwrf(i,k,j) = qi(k)
         qswrf(i,k,j) = qs(k)
         qgwrf(i,k,j) = qg(k)
         SMLF(i,k,j)  = SMLF1D(k)
         GMLF(i,k,j)  = GMLF1D(k)
      ENDDO
      ! -----------------------------------------------------------------!
 
      ENDDO  ! end of i_loop
      ENDDO  ! end of j_loop

      return
      END SUBROUTINE saticel_s
!=======================================================================
!
!========================================================================
      REAL FUNCTION CHEN(RHO,AFAR,LAMR)                                 ! Raindrops terminal velocity
!========================================================================
      IMPLICIT NONE
      REAL, INTENT(IN) :: RHO,AFAR,LAMR
      REAL :: RHOQ,VRA1,VRA2,VRA3,VRB1,VRB2,VRB3,VRC1,VRC2,VRC3,LAMM

      RHOQ  = EXP(0.115231*RHO)
      VRA1  = 0.044612*RHOQ
      VRA2  = -0.263166*RHOQ
      VRA3  = 4.7178*RHOQ*RHO**(-0.47335)
      VRB1  = 2.2955-0.038465*RHO
      VRB2  = 2.2955-0.038465*RHO
      VRB3  = 1.1451-0.038465*RHO
      VRC1  = 0.
      VRC2  = 0.184325
      VRC3  = 0.184325
      LAMM  = LAMR/1.E3
      CHEN  = VRA1*EXP((AFAR+1.)*LOG(LAMM)+LGAMMA(VRB1+AFAR+1.)-(VRB1+ &
              AFAR+1.)*LOG(LAMM+VRC1)-LGAMMA(AFAR+1.))+VRA2*EXP((AFAR+ &
              1.)*LOG(LAMM)+LGAMMA(VRB2+AFAR+1.)-(VRB2+AFAR+1.)*       &
              LOG(LAMM+VRC2)-LGAMMA(AFAR+1.))+VRA3*EXP((AFAR+1.)*      &
              LOG(LAMM)+LGAMMA(VRB3+AFAR+1.)-(VRB3+AFAR+1.)*LOG(LAMM+  &
              VRC3)-LGAMMA(AFAR+1.))

      END FUNCTION CHEN
!======================================================================

!======================================================================
      END MODULE  module_mp_gsfcgce_3ice_cwb
!======================================================================
