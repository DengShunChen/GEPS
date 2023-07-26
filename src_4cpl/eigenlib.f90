      SUBROUTINE RSB(NM,N,MB,A,W,MATZ,Z,FV1,FV2,IERR)
!
      INTEGER N,MB,NM,IERR,MATZ
      REAL A(NM,MB),W(N),Z(NM,N),FV1(N),FV2(N)
      LOGICAL TF
!
!     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
!     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
!     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
!     OF A REAL SYMMETRIC BAND MATRIX.
!
!     ON INPUT
!
!        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
!        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!        DIMENSION STATEMENT.
!
!        N  IS THE ORDER OF THE MATRIX  A.
!
!        MB  IS THE HALF BAND WIDTH OF THE MATRIX, DEFINED AS THE
!        NUMBER OF ADJACENT DIAGONALS, INCLUDING THE PRINCIPAL
!        DIAGONAL, REQUIRED TO SPECIFY THE NON-ZERO PORTION OF THE
!        LOWER TRIANGLE OF THE MATRIX.
!
!        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
!        BAND MATRIX.  ITS LOWEST SUBDIAGONAL IS STORED IN THE
!        LAST  N+1-MB  POSITIONS OF THE FIRST COLUMN, ITS NEXT
!        SUBDIAGONAL IN THE LAST  N+2-MB  POSITIONS OF THE
!        SECOND COLUMN, FURTHER SUBDIAGONALS SIMILARLY, AND
!        FINALLY ITS PRINCIPAL DIAGONAL IN THE  N  POSITIONS
!        OF THE LAST COLUMN.  CONTENTS OF STORAGES NOT PART
!        OF THE MATRIX ARE ARBITRARY.
!
!        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
!        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
!        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
!
!     ON OUTPUT
!
!        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
!
!        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
!
!        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
!           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
!           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
!
!        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (MB .GT. 0) GO TO 10
      IERR = 12 * N
      GO TO 50
   10 IF (MB .LE. N) GO TO 15
      IERR = 12 * N
      GO TO 50
!
   15 IF (MATZ .NE. 0) GO TO 20
!     .......... FIND EIGENVALUES ONLY ..........
      TF = .FALSE.
      CALL  BANDR(NM,N,MB,A,W,FV1,FV2,TF,Z)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 TF = .TRUE.
      CALL  BANDR(NM,N,MB,A,W,FV1,FV1,TF,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END

      SUBROUTINE BANDR(NM,N,MB,A,D,E,E2,MATZ,Z)
      INTEGER J,K,L,N,R,I1,I2,J1,J2,KR,MB,MR,M1,NM,N2,R1,UGL,MAXL,MAXR
      REAL A(NM,MB),D(N),E(N),E2(N),Z(NM,N)
      REAL G,U,B1,B2,C2,F1,F2,S2,DMIN,DMINRT
      LOGICAL MATZ
      DMIN = 2.0E0**(-64)
      DMINRT = 2.0E0**(-32)
      DO 30 J = 1, N
   30 D(J) = 1.0E0
      IF (.NOT. MATZ) GO TO 60
      DO 50 J = 1, N
         DO 40 K = 1, N
   40    Z(J,K) = 0.0E0
         Z(J,J) = 1.0E0
   50 CONTINUE
   60 M1 = MB - 1
      IF (M1 - 1) 900, 800, 70
   70 N2 = N - 2
      DO 700 K = 1, N2
         MAXR = MIN0(M1,N-K)
         DO 600 R1 = 2, MAXR
            R = MAXR + 2 - R1
            KR = K + R
            MR = MB - R
            G = A(KR,MR)
            A(KR-1,1) = A(KR-1,MR+1)
            UGL = K
            DO 500 J = KR, N, M1
               J1 = J - 1
               J2 = J1 - 1
               IF (G .EQ. 0.0E0) GO TO 600
               B1 = A(J1,1) / G
               B2 = B1 * D(J1) / D(J)
               S2 = 1.0E0 / (1.0E0 + B1 * B2)
               IF (S2 .GE. 0.5E0 ) GO TO 450
               B1 = G / A(J1,1)
               B2 = B1 * D(J) / D(J1)
               C2 = 1.0E0 - S2
               D(J1) = C2 * D(J1)
               D(J) = C2 * D(J)
               F1 = 2.0E0 * A(J,M1)
               F2 = B1 * A(J1,MB)
               A(J,M1) = -B2 * (B1 * A(J,M1) - A(J,MB)) - F2 + A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J,MB) + F1) + A(J1,MB)
               A(J,MB) = B1 * (F2 - F1) + A(J,MB)
               DO 200 L = UGL, J2
                  I2 = MB - J + L
                  U = A(J1,I2+1) + B2 * A(J,I2)
                  A(J,I2) = -B1 * A(J1,I2+1) + A(J,I2)
                  A(J1,I2+1) = U
  200          CONTINUE
               UGL = J
               A(J1,1) = A(J1,1) + B2 * G
               IF (J .EQ. N) GO TO 350
               MAXL = MIN0(M1,N-J1)
               DO 300 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = A(I1,I2) + B2 * A(I1,I2+1)
                  A(I1,I2+1) = -B1 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2) = U
  300          CONTINUE
               I1 = J + M1
               IF (I1 .GT. N) GO TO 350
               G = B2 * A(I1,1)
  350          IF (.NOT. MATZ) GO TO 500
               DO 400 L = 1, N
                  U = Z(L,J1) + B2 * Z(L,J)
                  Z(L,J) = -B1 * Z(L,J1) + Z(L,J)
                  Z(L,J1) = U
  400          CONTINUE
               GO TO 500
  450          U = D(J1)
               D(J1) = S2 * D(J)
               D(J) = S2 * U
               F1 = 2.0E0 * A(J,M1)
               F2 = B1 * A(J,MB)
               U = B1 * (F2 - F1) + A(J1,MB)
               A(J,M1) = B2 * (B1 * A(J,M1) - A(J1,MB)) + F2 - A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J1,MB) + F1) + A(J,MB)
               A(J,MB) = U
               DO 460 L = UGL, J2
                  I2 = MB - J + L
                  U = B2 * A(J1,I2+1) + A(J,I2)
                  A(J,I2) = -A(J1,I2+1) + B1 * A(J,I2)
                  A(J1,I2+1) = U
  460          CONTINUE
               UGL = J
               A(J1,1) = B2 * A(J1,1) + G
               IF (J .EQ. N) GO TO 480
               MAXL = MIN0(M1,N-J1)
               DO 470 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = B2 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2+1) = -A(I1,I2) + B1 * A(I1,I2+1)
                  A(I1,I2) = U
  470          CONTINUE
               I1 = J + M1
               IF (I1 .GT. N) GO TO 480
               G = A(I1,1)
               A(I1,1) = B1 * A(I1,1)
  480          IF (.NOT. MATZ) GO TO 500
               DO 490 L = 1, N
                  U = B2 * Z(L,J1) + Z(L,J)
                  Z(L,J) = -Z(L,J1) + B1 * Z(L,J)
                  Z(L,J1) = U
  490          CONTINUE
  500       CONTINUE
  600    CONTINUE
         IF (MOD(K,64) .NE. 0) GO TO 700
         DO 650 J = K, N
            IF (D(J) .GE. DMIN) GO TO 650
            MAXL = MAX0(1,MB+1-J)
            DO 610 L = MAXL, M1
  610       A(J,L) = DMINRT * A(J,L)
            IF (J .EQ. N) GO TO 630
            MAXL = MIN0(M1,N-J)
            DO 620 L = 1, MAXL
               I1 = J + L
               I2 = MB - L
               A(I1,I2) = DMINRT * A(I1,I2)
  620       CONTINUE
  630       IF (.NOT. MATZ) GO TO 645
            DO 640 L = 1, N
  640       Z(L,J) = DMINRT * Z(L,J)
  645       A(J,MB) = DMIN * A(J,MB)
            D(J) = D(J) / DMIN
  650    CONTINUE
  700 CONTINUE
  800 DO 810 J = 2, N
  810 E(J) = SQRT(D(J))
!fj--begin
!fj  When MATZ is true, E and E2 share the same address
!fj  and produces incorrect solution.  Following
!fj  modification avoids such aliasing issue.
!fj
      IF ( MATZ ) THEN
        DO 830 J = 1, N
          DO 820 K = 2, N
  820      Z(J,K) = E(K) * Z(J,K)
  830   CONTINUE
        U = 1.0E0
        DO 850 J = 2, N
          A(J,M1) = U * E(J) * A(J,M1)
          U = E(J)
          A(J,MB) = D(J) * A(J,MB)
          D(J) = A(J,MB)
          E(J) = A(J,M1)
  850   CONTINUE
        D(1) = A(1,MB)
        E(1) = 0.0E0
      ELSE
        U = 1.0E0
        DO 855 J = 2, N
          A(J,M1) = U * E(J) * A(J,M1)
          U = E(J)
          E2(J) = A(J,M1) ** 2
          A(J,MB) = D(J) * A(J,MB)
          D(J) = A(J,MB)
          E(J) = A(J,M1)
  855   CONTINUE
        D(1) = A(1,MB)
        E(1) = 0.0E0
        E2(1) = 0.0E0
      ENDIF
!fj--end
      GO TO 1001
  900 DO 950 J = 1, N
         D(J) = A(J,MB)
         E(J) = 0.0E0
         E2(J) = 0.0E0
  950 CONTINUE
 1001 RETURN
      END

      SUBROUTINE TQLRAT(N,D,E2,IERR)
!
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL D(N),E2(N)
      REAL B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
!     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
!
!     ON INPUT
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
!          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!          THE SMALLEST EIGENVALUES.
!
!        E2 HAS BEEN DESTROYED.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
!
      F = 0.0E0
      T = 0.0E0
      E2(N) = 0.0E0
!
      DO 290 L = 1, N
         J = 0
         H = ABS(D(L)) + SQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * S)
         R = PYTHAG(P,1.0E0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
!
         DO 140 I = L1, N
  140    D(I) = D(I) - H
!
         F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0E0) G = B
         H = G
         S = 0.0E0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0E0) G = B
            H = G * P / R
  200    CONTINUE
!
         E2(L) = S * G
         D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0E0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0E0) GO TO 130
  210    P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
!
  250    I = 1
  270    D(I) = P
  290 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END

      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
!
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1.
!
!        E HAS BEEN DESTROYED.
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      F = 0.0E0
      TST1 = 0.0E0
      E(N) = 0.0E0
!
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (TST1 .LT. H) TST1 = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = PYTHAG(P,1.0E0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
!
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + ABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END

      REAL FUNCTION EPSLON (X)
      REAL X
!
!     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
!
      REAL A,B,C,EPS
!
!     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
!     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
!        1.  THE BASE USED IN REPRESENTING FLOATING POINT
!            NUMBERS IS NOT A POWER OF THREE.
!        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
!            THE ACCURACY USED IN FLOATING POINT VARIABLES
!            THAT ARE STORED IN MEMORY.
!     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
!     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
!     ASSUMPTION 2.
!     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
!            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
!            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
!            C  IS NOT EXACTLY EQUAL TO ONE,
!            EPS  MEASURES THE SEPARATION OF 1.0 FROM
!                 THE NEXT LARGER FLOATING POINT NUMBER.
!     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
!     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
!
!     THIS VERSION DATED 4/6/83.
!
      A = 4.0E0/3.0E0
   10 B = A - 1.0E0
      C = B + B + B
      EPS = ABS(C-1.0E0)
      IF (EPS .EQ. 0.0E0) GO TO 10
      EPSLON = EPS*ABS(X)
      RETURN
      END

      REAL FUNCTION PYTHAG(A,B)
      REAL A,B
!
!     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
!
      REAL P,R,S,T,U
!     P = AMAX1(ABS(A),ABS(B))
      P = DMAX1(ABS(A),ABS(B))        ! CWB2015
      IF (P .EQ. 0.0E0) GO TO 20
!     R = (AMIN1(ABS(A),ABS(B))/P)**2
      R = (DMIN1(ABS(A),ABS(B))/P)**2 ! CWB2015
   10 CONTINUE
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         U = 1.0E0 + 2.0E0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
