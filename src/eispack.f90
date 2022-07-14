      SUBROUTINE RG(NM,N,A,WR,WI,MATZ,Z,IV1,FV1,IERR)
!!
!
      INTEGER N,NM,IS1,IS2,IERR,MATZ
      REAL A(NM,N),WR(N),WI(N),Z(NM,N),FV1(N)
      INTEGER IV1(N)
!
!
      IF (N.LE.NM) GO TO 10
      IERR=10*N
      GO TO 50
!
   10 CALL BALANC(NM,N,A,IS1,IS2,FV1)
      CALL ELMHES(NM,N,IS1,IS2,A,IV1)
      IF (MATZ.NE.0) GO TO 20
! ************* FIND EIGEN VALUES ONLY *********************
      CALL HQR(NM,N,IS1,IS2,A,WR,WI,IERR)
      GO TO 50
! ************* FIND BOTH EIGENVALUES AND EIGENVECTORS *****
   20 CALL ELTRAN(NM,N,IS1,IS2,A,IV1,Z)
      CALL HQR2(NM,N,IS1,IS2,A,WR,WI,Z,IERR)
      IF (IERR.NE.0) GO TO 50
      CALL BALBAK(NM,N,IS1,IS2,FV1,N,Z)
   50 RETURN
      END
!
!
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
!
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      REAL ABS
      LOGICAL NOCONV
!
! ***************** RADIX is a MACHINE DEPENDENT PARAMETER SPACIFYING
!                   THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION
! *************************************************************************
!
      RADIX=2.
!
      B2=RADIX*RADIX
      K=1
      L=N
      GO TO 100
! ***************IN-LINE PROCEDURE FOR ROW AND COLUMN EXCHANG *******
   20 SCALE(M)=J
      IF (J.EQ.M) GO TO 50
!
      DO 30 I=1,L
      F=A(I,J)
      A(I,J)=A(I,M)
      A(I,M)=F
   30 CONTINUE
!
      DO 40 I=K,N
      F=A(J,I)
      A(J,I)=A(M,I)
      A(M,I)=F
   40 CONTINUE
!
   50 GO TO (80,130), IEXC
! *************** SEARCH FOR ROWS ISOLATING AN EIGENVALUE AND PUSH THEM DOWN
   80 IF (L.EQ.1) GO TO 280
      L=L-1
! *************** FOR J+L STEP -1 UNTIL 1 do -- ********************
  100 DO 120 JJ=1, L
      J=L+1-JJ
!
      DO 110 I=1,L
        IF (I.EQ.J) GO TO 110
        IF (A(J,I).NE.0.) GO TO 120
  110 CONTINUE
!
      M=L
      IEXC=1
      GO TO 20
  120 CONTINUE
      GO TO 140
! *************** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE AND
!                 PUSH THEM LEFT **********************
  130 K=K+1
!
  140 DO 170 J=K,L
!
        DO 150 I=K,L
          IF(I.EQ.J) GO TO 150
          IF(A(I,J).NE.0.) GO TO 170
  150   CONTINUE
!
        M=K
        IEXC=2
        GO TO 20
  170 CONTINUE
! ***************** NOW BALANCE THE SUBMATRIX IN ROWS K TO L *************
      DO 180 I=K,L
  180 SCALE(I)=1.0
! ***************** ITERATIVE LOOP FOR NORM REDUCTION **************
  190 NOCONV=.FALSE.
!
      DO 270 I=K,L
        C=0.
        R=0.
!
      DO 200 J=K,L
        IF (J.EQ.I) GO TO 200
        C=C+ABS(A(J,I))
        R=R+ABS(A(I,J))
  200 CONTINUE
! ***************** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW *********
      IF (C.EQ.0.OR.R.EQ.0.) GO TO 270
      G=R/RADIX
      F=1.
      S=C+R
  210 IF(C.GE.G) GO TO 220
      F=F*RADIX
      C=C*B2
      GO TO 210
  220 G=R*RADIX
  230 IF (C.LT.G) GO TO 240
      F=F/RADIX
      C=C/B2
      GO TO 230
! *********************** NOW BALANCE *******************************
  240 IF ((C+R)/F.GE.0.95*S) GO TO 270
      G=1./F
      SCALE(I)=SCALE(I)*F
      NOCONV=.TRUE.
!
      do 250 J=K,N
  250 A(I,J)=A(I,J)*G
!
      DO 260 J=1,L
  260 A(J,I)=A(J,I)*F
!
  270 CONTINUE
!
      IF(NOCONV) GO TO 190
!
  280 LOW=K
      IGH=L
      RETURN
      END
!
!
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
!
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),Z(NM,M)
      REAL S
!
      IF(M.EQ.0) GO TO 200
      IF(IGH.EQ.LOW) GO TO 120
!
      DO 110 I=LOW,IGH
        S=SCALE(I)
! **************** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!                  IF THE FOREGOING STATEMENT IS REPLACED BY
!                  S=1.0/SCALE(i). *******************************
      DO 100 J=1,M
  100 Z(I,J)=Z(I,J)*S
!
  110 CONTINUE
! *****************_ FOR I=LOW-1 UNTIL 1,
!                    IGH+1 STEP 1 until N DO __ ********************
  120 DO 140 II=1,N
        I=II
        IF(I.GE.LOW.AND.I.LE.IGH) GO TO 140
        IF(I.LT.LOW) I=LOW-II
        K=SCALE(I)
        IF(K.EQ.I) GO TO 140
!
        DO 130 J=1,M
          S=Z(I,J)
          Z(I,J)=Z(K,J)
          Z(K,J)=S
  130   CONTINUE
!
  140 CONTINUE
!
  200 RETURN
      END
!
!
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
!
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
!      REAL ABS
      INTEGER INT(IGH)
!
      LA=IGH-1
      KP1=LOW+1
      IF(LA.LT.KP1) GO TO 200
!
      DO 180 M=KP1,LA
        MM1=M-1
        X=0.
        I=M
!
        DO 100 J=M,IGH
          IF (ABS(A(J,MM1)).LE.ABS(X)) GO TO 100
            X=A(J,MM1)
            I=J
  100   CONTINUE
!
        INT(M)=I
        IF (I.EQ.M) GO TO 130
! ***************** INTERCHANGE ROWS AND CLOUMNS OF A *****************
        DO 110 J=MM1,N
          Y=A(I,J)
          A(I,J)=A(M,J)
          A(M,J)=Y
  110   CONTINUE
!
        DO 120 J=1,IGH
          Y=A(J,I)
          A(J,I)=A(J,M)
          A(J,M)=Y
  120   CONTINUE
! ****************** END INTERCHANGE ***********************
  130   IF(X.EQ.0.) GO TO 180
        MP1=M+1
!
        DO 160 I=MP1,IGH
          Y=A(I,MM1)
          IF(Y.EQ.0.) GO TO 160
          Y=Y/X
          A(I,MM1)=Y
!
          DO 140 J=M,N
  140     A(I,J)=A(I,J)-Y*A(M,J)
!
          DO 150 J=1,IGH
  150     A(J,M)=A(J,M)+Y*A(J,I)
!
  160   CONTINUE
  180 CONTINUE
!
  200 RETURN
      END
!
!
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
!
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
!
! ***************** INITIALIZE Z TO IDENTITY MATRIX **************
      DO 80 I=1,N
!
        DO 60 J=1,N
   60   Z(I,J)=0.
!
        Z(I,I)=1.
   80 CONTINUE
!
      KL=IGH-LOW-1
      IF (KL.LT.1) GO TO 200
! ***************** FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO __ ********
      DO 140 MM=1,KL
        MP=IGH-MM
        MP1=MP+1
!
        DO 100 I=MP1,IGH
  100   Z(I,MP)=A(I,MP-1)
!
        I=INT(MP)
        IF(I.EQ.MP) GO TO 140
!
        DO 130 J=MP,IGH
          Z(MP,J)=Z(I,J)
          Z(I,J)=0.
  130   CONTINUE
!
        Z(I,MP)=1.
  140 CONTINUE
!
  200 RETURN
      END
!
!
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
!
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N)
      REAL P,Q,R,S,T,W,X,Y,ZZ,NORM,MACHEP
!      REAL SQRT,ABS,SIGN
!      INTEGER MIN0
      LOGICAL NOLATS
!
! **************** MACHEP IS A MACHINE DEPENDENT PARAMETER SPACIFYING
!                  THE RELATIVE PRECISION OF FLOATING POINT ARITHEMETIC
! ************************************************************************
      MACHEP=2.**(-47)
!
      IERR=0
      NORM=0.
      K=1
! **************** STORE ROOTS OSOLATED BY BALANC
!                  AND COMPUTE MATRIX NORM ***********************
      DO 50 I=1,N
!
        DO 40 J=K,N
   40   NORM=NORM+ABS(H(I,J))
!
        K=I
        IF (I.GE.LOW.AND.I.LE.IGH) GO TO 50
        WR(I)=H(I,J)
        WI(I)=0.
   50 CONTINUE
!
      EN=IGH
      T=0.
! ***************** SEARCH FOR NEXT EIGENVALUES *******************
   60 IF (EN.LT.LOW) GO TO 1001
      ITS=0
      NA=EN-1
      ENM2=NA-1
! ***************** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                   FOR L=en STEP -1 UNTIL LOW DO __ ***************
   70 DO 80 LL=LOW,EN
        L=EN+LOW-LL
        IF (L.EQ.LOW) GO TO 100
        S=ABS(H(L-1,L-1))+ABS(H(L,L))
        IF (S.EQ.0.) S=NORM
        IF (ABS(H(L,L-1)).LE.MACHEP*S) GO TO 100
   80 CONTINUE
! ***************** FORM SHIFT **************************************
  100 X=H(EN,EN)
      IF(L.EQ.EN) GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF (L.EQ.NA) GO TO 280
      IF (ITS.EQ.30) GO TO 1000
      IF (ITS.NE.10.AND.ITS.NE.20) GO TO 130
! ***************** FORM EXCEPTIONAL SHIFT *************************
      T=T+X
!
      DO 120 I=LOW,EN
  120 H(I,I)=H(I,I)-X
!
      S=ABS(H(EN,NA))+ABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
  130 ITS=ITS+1
! ***************** LOOK FOR TWO CONSECUTIVE SMALL
!                   SUB-DIAGONAL ELEMENTS.
!                   FOR m=EN-2 STEP-1 UNTIL L DO__ ******************
      DO 140 MM=L,ENM2
        M=ENM2+L-MM
        ZZ=H(M,M)
        R=X-ZZ
        S=Y-ZZ
        P=(R*S-W)/H(M+1,M)+H(M,M+1)
        Q=H(M+1,M+1)-ZZ-R-S
        R=H(M+2,M+1)
        S=ABS(P)+ABS(Q)+ABS(R)
        P=P/S
        Q=Q/S
        R=R/S
        IF (M.EQ.L) GO TO 150
        IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.MACHEP*ABS(P)       &
         * (ABS(H(M-1,M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
!
  150 MP2=M+2
!
      DO 160 I=MP2,EN
        H(I,I-2)=0.
        IF (I.EQ.MP2) GO TO 160
        H(I,I-3)=0.
  160 CONTINUE
! ******************** DOUBLE QR STEP INVOLVING ROWS L TO EN AND
!                      COLUMNS M TO EN **************************
      DO 260 K=M,NA
        NOLATS=K.NE.NA
        IF (K.EQ.M) GO TO 170
        P=H(K,K-1)
        Q=H(K+1,K-1)
        R=0.
        IF (NOLATS) R=H(K+2,K-1)
        X=ABS(P)+ABS(Q)+ABS(R)
        IF (X.EQ.0.) GO TO 260
        P=P/X
        Q=Q/X
        R=R/X
  170   S=SIGN(SQRT(P*P+Q*Q+R*R),P)
        IF (K.EQ.M) GO TO 180
        H(K,K-1)=-S*X
        GO TO 190
  180   IF(L.NE.M) H(K,K-1)=-H(K,K-1)
  190   P=P+S
        X=P/S
        Y=Q/S
        ZZ=R/S
        Q=Q/P
        R=R/P
! ******************** ROW MODIFICATION ***************************
        DO 210 J=K,EN
          P=H(K,J)+Q*H(K+1,J)
          IF (.NOT.NOLATS) GO TO 200
          P=P+R*H(K+2,J)
          H(K+2,J)=H(K+2,J)-P*ZZ
  200     H(K+1,J)=H(K+1,J)-P*Y
          H(K,J)=H(K,J)-P*X
  210   CONTINUE
!
        J=MIN0(EN,K+3)
! ******************** COLUMN MODIFICATION *************************
        DO 230 I=L,J
          P=X*H(I,K)+Y*H(I,K+1)
          IF (.NOT.NOLATS) GO TO 220
          P=P+ZZ*H(I,K+2)
          H(I,K+2)=H(I,K+2)-P*R
  220     H(I,K+1)=H(I,K+1)-P*Q
          H(I,K)=H(I,K)-P
  230   CONTINUE
!
  260 CONTINUE
!
      GO TO 70
! ******************** ONE ROOT FOUND ******************************
  270 WR(EN)=X+T
      WI(EN)=0.
      EN=NA
      GO TO 60
! ******************** TWO ROOTS FOUND *****************************
  280 P=(Y-X)/2.0
      Q=P*P+W
      ZZ=SQRT(ABS(Q))
      X=X+T
      IF (Q.LT.0.) GO TO 320
! ******************** REAL PAIR ***********************************
      ZZ=P+SIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF (ZZ.NE.0.) WR(EN)=X-W/ZZ
      WI(NA)=0.
      WI(EN)=0.
      GO TO 330
! ******************** COMPLEX PAIR *******************************
  320 WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
  330 EN=ENM2
      GO TO 60
! ******************** SET ERROR -- NO CONVERGENCE TO AN
!                      EIGENVALUE AFTER 30 ITERATIONS *************
 1000 IERR=EN
 1001 RETURN
      END
!
!
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
!
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN, &
              IGH,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
!      REAL SQRT,ABS,SIGN
!      INTEGER MIN0
      LOGICAL NOLATS
      COMPLEX Z3
!      COMPLEX CMPLX
!      REAL REAL,AIMAG
!
! ********************** MACHEP IS A MACHINE DEPENDENT PARAMETER SPACIFYING
!                        THE RELATIVE PRECISION OF FLOATING POINT ARITHEMETIC
      MACHEP=2.**(-47)
!
      IERR=0
      NORM=0.
      K=1
! ********************* STORE ROOTS ISOLATED BY BALANC
!                       AND COMPUTE MATRIX NORM **********************
      DO 50 I=1,N
!
        DO 40 J=K,N
  40    NORM=NORM+ABS(H(I,J))
!
        K=I
        IF (I.GE.LOW.AND.I.LE.IGH) GO TO 50
        WR(I)=H(I,I)
        WI(I)=0.
   50 CONTINUE
!
      EN=IGH
      T=0.
! ********************* SEARCH FOR NEXT EIGENVALUES ******************
   60 IF (EN.LT.LOW) GO TO 340
      ITS=0
      NA=EN-1
      ENM2=NA-1
! ********************* LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                       FOR L=EN STEP -1 UNTIL LOW DO __ **************
   70 DO 80 LL=LOW,EN
       L=EN+LOW-LL
       IF (L.EQ.LOW) GO TO 100
       S=ABS(H(L-1,L-1))+ABS(H(L,L))
       IF (S.EQ.0.) S=NORM
       IF (ABS(H(L,L-1)).LE.MACHEP*S) GO TO 100
   80 CONTINUE
! ******************** FORM SHIFT ***********************************
  100 X=H(EN,EN)
      IF (L.EQ.EN) GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF (L.EQ.NA) GO TO 280
      IF (ITS.EQ.30) GO TO 1000
      IF (ITS.NE.10.AND.ITS.NE.20) GO TO 130
! ******************** FORM EXCEPTIONAL SHIFT ***********************
      T=T+X
!
      DO 120 I=LOW,EN
  120 H(I,I)=H(I,I)-X
!
      S=ABS(H(EN,NA))+ABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
  130 ITS=ITS+1
! ******************** LOOK FOR TWO CONSECUTIVE SMALL
!                      SUB-DIAGONAL ELEMENTS.
!                      FOR m=EN-2 STEP -1 UNTIL L Do __ ************
      DO 140 MM=L,ENM2
      M=ENM2+L-MM
      ZZ=H(M,M)
      R=X-ZZ
      S=Y-ZZ
      P=(R*S-W)/H(M+1,M)+H(M,M+1)
      Q=H(M+1,M+1)-ZZ-R-S
      R=H(M+2,M+1)
      S=ABS(P)+ABS(Q)+ABS(R)
      P=P/S
      Q=Q/S
      R=R/S
      IF (M.EQ.L) GO TO 150
      IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)).LE.MACHEP*ABS(P)      &
       * (ABS(H(M-1,M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
!
  150 MP2=M+2
!
      DO 160 I=MP2,EN
        H(I,I-2)=0.
        IF (I.EQ.MP2) GO TO 160
        H(I,I-3)=0.
  160 CONTINUE
! ****************** DOUBLE QR STEP INVOLVING ROWS L TO EN AND
!                    COLUMNS M TO EN *****************************
      DO 260 K=M,NA
        NOLATS=K.NE.NA
        IF (K.EQ.M) GO TO 170
        P=H(K,K-1)
        Q=H(K+1,K-1)
        R=0.
        IF (NOLATS) R=H(K+2,K-1)
        X=ABS(P)+ABS(Q)+ABS(R)
        IF (X.EQ.0.) GO TO 260
        P=P/X
        Q=Q/X
        R=R/X
  170   S=SIGN(SQRT(P*P+Q*Q+R*R),P)
        IF (K.EQ.M) GO TO 180
        H(K,K-1)=-S*X
        GO TO 190
  180   IF (L.NE.M) H(K,K-1)=-H(K,K-1)
  190   P=P+S
        X=P/S
        Y=Q/S
        ZZ=R/S
        Q=Q/p
        R=R/P
! ****************** ROW MODIFICATION *************************
        DO 210 J=K,N
          P=H(K,J)+Q*H(K+1,J)
          IF (.NOT.NOLATS) GO TO 200
          P=P+R*H(K+2,J)
          H(K+2,J)=H(K+2,J)-P*ZZ
  200     H(K+1,J)=H(K+1,J)-P*Y
          H(K,J)=H(K,J)-P*X
  210   CONTINUE
!
        J=MIN0(EN,K+3)
! ****************** COLUMN MODIFICATION **********************
        DO 230 I=1,J
          P=X*H(I,K)+Y*H(I,K+1)
          IF (.NOT.NOLATS) GO TO 220
          P=P+ZZ*H(I,K+2)
          H(I,K+2)=H(I,K+2)-P*R
  220     H(I,K+1)=H(I,K+1)-P*Q
          H(I,K)=H(I,K)-P
  230   CONTINUE
! ****************** ACCUMULATE TRANSFORMATIONS ***************
        DO 250 I=LOW,IGH
          P=X*Z(I,K)+Y*Z(I,K+1)
          IF (.NOT.NOLATS) GO TO 240
          P=P+ZZ*Z(I,K+2)
          Z(I,K+2)=Z(I,K+2)-P*R
  240     Z(I,K+1)=Z(I,K+1)-P*Q
          Z(I,K)=Z(I,K)-P
  250   CONTINUE
!
  260 CONTINUE
!
      GO TO 70
! ***************** ONE ROOT FOUND ******************************
  270 H(EN,EN)=X+T
      WR(EN)=H(EN,EN)
      WI(EN)=0.
      EN=NA
      GO TO 60
! ***************** TWO ROOTS FOUND *****************************
  280 P=(Y-X)/2.
      Q=P*P+W
      ZZ=SQRT(ABS(Q))
      H(EN,EN)=X+T
      X=H(EN,EN)
      H(NA,NA)=Y+T
      IF (Q.LT.0.) GO TO 320
! ****************** REAL PAIR ***********************************
      ZZ=P+SIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF (ZZ.NE.0.) WR(EN)=X-W/ZZ
      WI(NA)=0.
      WI(EN)=0.
      X=H(EN,NA)
      S=ABS(X)+ABS(ZZ)
      P=X/S
      Q=ZZ/S
      R=SQRT(P*P+Q*Q)
      P=P/R
      Q=Q/R
! ****************** ROW MODIFICATION ***************************
      DO 290 J=NA,N
        ZZ=H(NA,J)
        H(NA,J)=Q*ZZ+P*H(EN,J)
        H(EN,J)=Q*H(EN,J)-P*ZZ
  290 CONTINUE
! ****************** COLUMN MODIFICATION ************************
      DO 300 I=1,EN
        ZZ=H(I,NA)
        H(I,NA)=Q*ZZ+P*H(I,EN)
        H(I,EN)=Q*H(I,EN)-P*ZZ
  300 CONTINUE
! ****************** ACCUMULATE TRANSFORMATION *******************
      DO 310 I=LOW,IGH
        ZZ=Z(I,NA)
        Z(I,NA)=Q*ZZ+P*Z(I,EN)
        Z(I,EN)=Q*Z(I,EN)-P*ZZ
  310 CONTINUE
!
      GO TO 330
! ***************** COMPLEX PAIR *********************************
  320 WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
  330 EN=ENM2
      GO TO 60
! ***************** ALL ROOTS FOUND. BACKSUBSTITUDE TO FIND
!                   VECTORS OF UPPER TRIANGULAR FORM *************
  340 IF (NORM.EQ.0.) GO TO 1001
! ***************** FOR EN=N STEP -1 UNTIL 1 DO __ ***************
      DO 800 NN=1,N
        EN=N+1-NN
        P=WR(EN)
        Q=WI(EN)
        NA=EN-1
        IF (Q) 710, 600, 800
! ***************** REAL VECTOR **********************************
  600 M=EN
      H(EN,EN)=1.
      IF (NA.EQ.0) GO TO 800
! ***************** FOR i=EN-1 STEP -1 UNTIL 1 DO __ *************
      DO 700 II=1,NA
        I=EN-II
        W=H(I,I)-P
        R=H(I,EN)
        IF (M.GT.NA) GO TO 620
!
        DO 610 J=M,NA
  610   R=R+H(I,J)*H(J,EN)
!
  620   IF (WI(I).GE.0.) GO TO 630
        ZZ=W
        S=R
        GO TO 700
  630   M=I
        IF (WI(I).NE.0.) GO TO 640
        T=W
        IF (W.EQ.0.) T=MACHEP*NORM
        H(I,EN)=-R/T
        GO TO 700
! *************** SOLVE REAL EQUATION *****************************
  640   X=H(I,I+1)
        Y=H(I+1,I)
        Q=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
        T=(X*S-ZZ*R)/Q
        H(I,EN)=T
        IF (ABS(X).LE.ABS(ZZ)) GO TO 650
        H(I+1,EN)=(-R-W*T)/X
        GO TO 700
  650   H(I+1,EN)=(-S-Y*T)/ZZ
  700 CONTINUE
! *************** END REAL VECTOR *********************************
      GO TO 800
! ************** COMPLEX VECOTR ***********************************
  710 M=NA
! ************** LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
!                EIGEN VECTOR MATRIX IS TRINAUGULAR ***************
      IF (ABS(H(EN,NA)).LE.ABS(H(NA,EN))) GO TO 720
      H(NA,NA)=Q/H(EN,NA)
      H(NA,EN)=-(H(EN,EN)-P)/H(EN,NA)
      GO TO 730
  720 Z3=CMPLX(0.0,-H(NA,EN))/CMPLX(H(NA,NA)-P,Q)
      H(NA,NA)=REAL(Z3)
      H(NA,EN)=AIMAG(Z3)
  730 H(EN,NA)=0.0
      H(EN,EN)=1.
      ENM2=NA-1
      IF (ENM2.EQ.0) GO TO 800
! ********************* FOR i=EN-2 STEP -1 UNTIL 1 DO __ ************
      DO 790 II=1,ENM2
        I=NA-II
        W=H(I,I)-P
        RA=0.
        SA=H(I,EN)
!
        DO 760 J=M,NA
          RA=RA+H(I,J)*H(J,NA)
          SA=SA+H(I,J)*H(J,EN)
  760   CONTINUE
!
        IF (WI(I).GE.0.) GO TO 770
        ZZ=W
        R=RA
        S=SA
        GO TO 790
  770   M=I
        IF (WI(I).NE.0.) GO TO 780
        Z3=CMPLX(-RA,-SA)/CMPLX(W,Q)
        H(I,NA)=REAL(Z3)
        H(I,EN)=AIMAG(Z3)
        GO TO 790
! ******************** SOLVE COMPLEX EQUATIONS *******************
  780   X=H(I,I+1)
        Y=H(I+1,I)
        VR=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
        VI=(WR(I)-P)*2.0*Q
        IF (VR.EQ.0.AND.VI.EQ.0.) VR=MACHEP*NORM    &
         * (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
        Z3=CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/CMPLX(VR,VI)
        H(I,NA)=REAL(Z3)
        H(I,EN)=AIMAG(Z3)
        IF (ABS(X).LE.ABS(ZZ)+ABS(Q)) GO TO 785
        H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X
        H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X
        GO TO 790
  785   Z3=CMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN))/CMPLX(ZZ,Q)
        H(I+1,NA)=REAL(Z3)
        H(I+1,EN)=AIMAG(Z3)
  790 CONTINUE
! ******************** END COMPLEX VECTOR ***************************
  800 CONTINUE
! ******************** END BACK SUBSTITUTION.
!                      VECTORS OF ISOLATED ROOTS ********************
      DO 840 I=1,N
        IF (I.GE.LOW.AND.I.LE.IGH) GO TO 840
!
        DO 820 J=1,N
  820   Z(I,J)=H(I,J)
!
  840 CONTINUE
! ******************** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!                      VECTORS OF ORIGINAL FULL MATRIX.
!                      FOR J=N STEP -1 UNTIL LOW DO __ **************
      DO 880 JJ=LOW,N
        J=N+LOW-JJ
        M=MIN0(J,IGH)
!
        DO 880 I=LOW,IGH
          ZZ=0.
!
          DO 860 K=LOW,M
  860     ZZ=ZZ+Z(I,K)*H(K,J)
!
          Z(I,J)=ZZ
  880 CONTINUE
!
      GO TO 1001
! ************** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS *********************
 1000 IERR=EN
 1001 RETURN
      END



