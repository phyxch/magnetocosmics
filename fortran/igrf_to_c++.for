
c
      SUBROUTINE IGRF_TO_CC_(R,T,F,BR,BT,BF)
c
C     CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN SPHERICAL
C     GEOGRAPHICAL COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE
C     MODEL COEFFICIENTS (e.g., http://www.ngdc.noaa.gov/IAGA/wg8/igrf2000.html)
C
C     UPDATING THE COEFFICIENTS TO A GIVEN EPOCH IS MADE AUTOMATICALLY UPON THE FIRST
C     CALL AND AFTER EVERY CHANGE OF THE PARAMETER IY.
C
C-----INPUT PARAMETERS:
C
C     IY  -  YEAR NUMBER (FOUR-DIGIT; 1965 &LE IY &LE 2005)
C     NM  -  HIGHEST ORDER OF SPHERICAL HARMONICS IN THE SCALAR POTENTIAL (NM &LE 10)
C     R,T,F -  SPHERICAL COORDINATES (RADIUS R IN UNITS RE=6371.2 KM, GEOGRAPHIC
C                COLATITUDE  T  AND LONGITUDE  F  IN RADIANS)
C
C-----OUTPUT PARAMETERS:
C
C     BR,BT,BF - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C
C     LAST MODIFICATION:  JANUARY 5, 2001.
C     THE CODE WAS MODIFIED TO ACCEPT DATES THROUGH 2005.
C     IT HAS ALSO BEEN SLIGHTLY SIMPLIFIED BY TAKING OUT SOME REDUNDANT STATEMENTS,
C     AND A "SAVE" STATEMENT WAS ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
C     FORTRAN COMPILERS.
c
C     WRITTEN BY: N. A. TSYGANENKO
C    SAVE MA,IYR,IPR
C
      DIMENSION A(14),B(14)
c
      COMMON /IGRFCC/ G(105),H(105),REC(105),NM

130   PP=1./R
      P=PP
      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N
      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.
      U=T
      CF=COS(F)
      SF=SIN(F)
      C=COS(U)
      S=SIN(U)
      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.
         Y=1.
170      Q=P
         Z=D
         BI=0.
         P2=0.
         D2=0.
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.E-5) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BT=BBT
      IF(S.LT.1.E-5) GOTO 210
      BF=BBF/S
      RETURN
210   IF(C.LT.0.) BBF=-BBF
      BF=BBF

      RETURN
C
999   FORMAT(//1X,
     * 'IGRF WARNS:**** YEAR IS OUT OF INTERVAL 1965-2005: IY =',I5,/,
     *',        CALCULATIONS WILL BE DONE FOR IYR =',I5,' ****'//)
      END
