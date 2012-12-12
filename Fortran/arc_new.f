*     REM FINITE DIFFERENCE PROGRAM  CONVECTION MODEL   
      COMMON SFILE,TITLE,SUFILE
      COMMON IFCV
      COMMON /FIRST/ X,Y,NR,NC,T,TI,NOD
      COMMON /SECOND/ D,CHP,DTC,CHF,TIME
      COMMON /THIRD/ NODCV,RVSH,FVSH,TPMIN,TIC,VEL,TOTIM,KSAVE,TP0

      REAL X(120),Y(120),T(120,120),TMSRC
      REAL TI,TIME,TPMIN(9),VEL(8),TNEW(5),STYRS,TOTIM
      REAL CHP(9),D(9),DTC,CHF,TP0,TYCK,RVSH(8),FVSH(8)

      INTEGER NOD(120,120),NODCV(120,120),INCLM,YEAR1,ZFNUM,IFCV
      INTEGER MTI(2), MTJ(2), MTI1(2), MTJ1(2), MTI2(2), MTJ2(2)
      INTEGER LLL, LLK, KSAVE, MSRC, MSRR, ITK, BLST

      CHARACTER*8 YEARS, LOOPS,LPUT,LCY
      CHARACTER*80 AFILE, SFILE, SUFILE, NULLS
      CHARACTER*80 TITLE
      LOGICAL    exists
      DATA TM,TN /0,15/
      DATA YEARS,LOOPS,LPUT,LCY /' Years',' Loops', ' Update','every '/    

      TIME = 0.0
      TP0 = 100.0
      ZFNUM = 2
      NS = 0
      YEAR1 = 1 
      INCLM = 1
      DTC = 0.25
      WRITE (*,*) '                   FINITE DIFFERENCE HEAT FLOW'
      WRITE (*,*) '                     ENTER FILE NAME OF MODEL'
      READ (*,15) AFILE      
      WRITE (*,*)   'TO SAVE THE FINAL DATA FILE ENTER 1 AT THE PROMPT'
   3  READ(*,*, ERR = 4) KSAVE
      GO TO 5
   4  WRITE (*,*) ' INVALID VALUE FOR KSAVE - ENTER A NEW VALUE'
      GO TO 3
   5  IF (KSAVE) 7,7,6
   6  WRITE (*,*) '          ENTER FILE NAME TO SAVE AS FINAL MODEL'
      READ (*,15)  SFILE
   7  WRITE (*,*) '          ENTER FILE NAME TO SAVE AS A SURFER FILE'
      READ (*,15)  SUFILE

   8  CONTINUE 
      INQUIRE (FILE = AFILE, EXIST = exists)
      IF (.NOT. exists) THEN
         WRITE(*,*) 'INVALID MODEL FILE NAME - ENTER A DIFFERENT NAME'
         READ (*,15) AFILE
      END IF
   16 FORMAT(F6.2)
      OPEN (10,FILE=AFILE, STATUS = 'OLD')
   10 READ (10,*) NR, NC, IFCV
      READ (10,*) CHF,TIME
       TOTIM = TIME
  
  11  FORMAT(2I4)
  12  FORMAT(F5.1,F12.0)
      DTC = .25
      ACHF = CHF
      CHF = CHF*.001

      READ (10,15) TITLE
   15 FORMAT(A80)
C
      WRITE(*,15) TITLE
      WRITE(*,*)'NUMBER OF ROWS =    ',NR
      WRITE(*,*)'NUMBER OF COLUMNS = ',NC
      WRITE(*,*)'NUMBER OF LOOPS =   ',NT
      IF (IFCV.EQ.0) THEN  
          WRITE (*,*) 'NO CONVECTION'
          ELSE 
          WRITE(*,*) 'CONVECTION SUBROUTINE IS ACCESSED'
      ENDIF
      WRITE(*,*)'CONSTANT HEAT FLOW AT BASE OF MODEL=',ACHF,' mW M^2'
      WRITE(*,*) 'MODEL TIME ELAPSED = ',TIME
      TMAX = 0
      TMIN = 1000
      DO 30, I = 1, NR
       READ(10,*) (T(I,J),J=1,NC)
   30 CONTINUE
   21 FORMAT(3F9.3)
      WRITE(*,*) 'READ TEMPS IN A ',NC,'COL ',NR,'ROW ARRAY'

      DO 40, I = 1, NR
       READ(10,*) (NOD(I,J),J=1,NC)

   40 CONTINUE
   22 FORMAT(3I6)
      WRITE(*,*) 'READ CONDUCTION CODES '

      IF (IFCV.EQ.1) THEN
        DO I = 1,NR
          READ(10,*) (NODCV(I,J),J=1,NC)
        END DO
      WRITE(*,*) 'READ CONVECTION CODES'
      END IF
      READ(10,*) (X(I),I=1,NC)
   51 CONTINUE 
   23 FORMAT(11F7.3)
      READ(10,*) (Y(I),I=1,NR)
      WRITE(*,*) ' READ DIMENSIONS OF A ',NC,'COL ',NR,'ROW ARRAY'
      READ(10,*) (CHP(I),I=1,8)

      DO 53, I = 1,8
        CHP(I) = CHP(I)/1E6
   53 CONTINUE
   24 FORMAT(8F6.3)
      WRITE(*,*) ' READ 8 HEAT PRODUCTION VALUES'

      READ(10,*) (D(I),I=1,8)
      DO I = 1,8
        D(I) = D(I)*14.33
      END DO
      WRITE(*,*) 'CONVERTED 8 THERMAL CONDUCTIVITIES TO DIFF. IN m^2/y'
        DMAX = 0
        DO I = 1,8
           IF (D(I).GT.DMAX) THEN
              DMAX = D(I)
           END IF
        END DO
       WRITE(*,*) (D(I),I=1,8)
 760   CONTINUE

      IF (IFCV.EQ.1) THEN
        READ(10,*) (FVSH(I),I=1,8)
        READ(10,*) (RVSH(I),I=1,8)
        READ(10,*) (TPMIN(I),I=1,8)
        READ(10,*) (VEL(I),I=1,8)
        WRITE(*,*) ' READ 8 VELOCITIES IN m/yr'
        write(*,761) (VEL(I),I=1,8)
  761 FORMAT (8F10.2)
        V1=VEL(1)
        DO K1 = 2,8
          IF (V1.LT.VEL(K1)) V1 = VEL(K1)
        END DO
      END IF
      CLOSE (10)

      WRITE(*,*) 'PRESS ENTER TO CONTINUE'
      
      TIME = 0.0
       T1 = D(1)

       DO 1030, M = 1,8
          IF (T1.LT.D(M))  T1 = D(M)
 1030  CONTINUE

      T2 = X(1)
      T3 = Y(1)

      DO 1080, I = 1, NC
          IF (T2.GT.X(I))  T2 = X(I)
 1080 CONTINUE
       DO 1120, J = 1, NR
          IF (T3.GT.Y(J))  T3 = Y(J)
 1120    CONTINUE
C
C      CALCULATE TIME INCREMENT (TI) FOR MODEL
C      
        IF (IFCV) 1180, 1180, 1140
 1140   IF (T2.GT.T3) THEN
           TIC = T3/V1
        ELSE
           TIC = T2/V1
        ENDIF
 1180     IF (T2.LT.T3) THEN 
 1190       TI = T2 * T2 / (5 * T1)
       ELSE
 1210       TI = T3 * T3 / (5 * T1)
       ENDIF
      
      READ(*,*) ! PRESS ENTER TO CONTINUE
       
      WRITE (*,*) 'TO CHANGE TEMP. ON A BLOCK, ENTER 1, ELSE 0 '                                                  
      READ (*,*)  CGT
      IF (CGT.EQ.1) THEN
        WRITE (*,*) ' ENTER NUMBER OF RECTANGULAR BLOCKS TO CHANGE'
        READ(*,*) NBTC
       DO IB = 1, NBTC
        WRITE (*,*) 'BLOCK',IB
        WRITE (*,*)'ENTER COORDINATES OF UPPER LEFT CORNER row, column' 
        READ(*,*) I,J
        MTI1(IB)=I
        MTJ1(IB)=J
        WRITE (*,*)'ENTER COORDINATES OF LOWER RIGHT CORNER row, column'
        READ(*,*) II,JJ
        MTI(IB)= I+(II-I)/2
        MTJ(IB)= J+(JJ-I)/2
        MTI2(IB)=II
        MTJ2(IB)=JJ
        WRITE (*,*) 'CURRENT TEMPS FOR THE COORDINATES ARE'
        WRITE (*,*) I,J,T(I,J)
        WRITE (*,*) II,JJ,T(II,JJ)
        WRITE(*,*) 'ENTER NEW TEMPERATURE FOR THE BLOCK '
        READ(*,*) TNEW(IB)
          
          DO K = I, II
            DO KK = J, JJ
              T(K,KK) = TNEW(IB)
            END DO
          END DO
       
       END DO
        
        TM = TNEW(IB) 
      ENDIF
      
      
      DO I = 1,NR
       DO J = 1,NC
        IF (TM.LT.T(I,J)) TM=T(I,J)
        IF (TN.GT.T(I,J)) TN=T(I,J)
       END DO
      END DO
      RANGE = TM-TN
      IF (RANGE.EQ.0) RANGE = 1.0
      SPAN = 55/RANGE
      
      WRITE (*,*) 'TO START A MOVING SOURCE ENTER 1, ELSE ENTER 0'
      READ (*,*) MSRC
      IF (MSRC.NE.0) THEN
        WRITE (*,*) 'ENTER STARTING AND ENDING COLUMNS'
        READ (*,*) ITK,BLST
        WRITE (*,*) ' ENTER DEPTH TO SOURCE IN METERS'
        DEPMAX = INT(NR*Y(1))
        WRITE (*,*) ' MAXIMUM DEPTH IS ',DEPMAX, 'meters'
        READ (*,*) MSRR
        MSRR = INT(MSRR/Y(1))
        WRITE (*,*) 'ENTER TEMPERATURE OF MOVING SOURCE'
        READ (*,*) TMSRC
      ENDIF
 895   OPEN (12,FILE=SUFILE,STATUS='UNKNOWN')

      WRITE(12,1222) 'DSAA'
 1222 FORMAT(A4)
      WRITE(12,1223) NR,NC
      WRITE(12,1223) 0,NR
      WRITE(12,1223) 0,NC
 1223 FORMAT(2I4)
      DO I = 1,NR
       DO J = 1,NC
        IF (TM.LT.T(I,J)) TM=T(I,J)
        IF (TN.GT.T(I,J)) TN=T(I,J)
       END DO
      END DO
      WRITE(12,1224) TN,TM
 1224 FORMAT(2F7.2)

        TOUT = TI
       WRITE(*,*) 'EACH ITERATION IN TIME SPANS',TOUT,'YEARS'
       WRITE(*,*) 'ENTER A SHORTER ITERATION TIME IN YEARS IF DESIRED'
       READ(*,*) STYRS

        IF (STYRS.LT.TOUT) TI = STYRS
 1220     NT = 0
C      CALL clearscreen( $GCLEARSCREEN )

      NP = 0
      DLST =100.0
      DO I = 1,NR
        IF(DLST.GT.Y(I)) DLST = Y(I)
      END DO
      THTCST = DLST*DLST/DMAX
      WRITE (*,*) 'The thermal time constant for the vertical dimension 
     +is',THTCST, ' years' 
      WRITE (*,*) 'Enter time duration for calculation in YEARS'
      Read (*,*) TP0

      WRITE (*,*) 'Enter number of loops between screen updates - '
      READ (*,*) LLK

 2136 CONTINUE
      WRITE (*,*) 'Press ENTER to begin - '
      READ (*,15) NULLS

       NT = 0
       LLL = LLK
 1890  NT = NT + 1
       TOUT = TIME
       TYCK = TIME
       TOTIM = TOTIM + TOUT
       LLL = LLL - 1
        
       IF (LLL.LE.0) THEN 
           CALL SAVIT
           LLL = LLK

		   IF (NT.GT.LLK) THEN
			  WRITE(*,1892) NT,LOOPS,TOUT,YEARS
		   ELSE     
			  WRITE(*,1894) LPUT, LCY, LLK, LOOPS
		   ENDIF    
         END IF
 2215  FORMAT(I3)
 2216  FORMAT(F5.2)

 2103 FORMAT(10F8.2)
 1892  FORMAT(I6,A8,5X,F11.2,A8) 
 1893 FORMAT(I6,5X,I6)
 1894 FORMAT(1XA8,A8,I4,A8)
 2144  CONTINUE

C       THIS IS WHERE THE CALL TO THE MOVING SOURCE IS MADE        
c       The only direction at present is left to right and the move        
c       is made with each iteration.  Determine velocity and set 
c       iteration time according to the column spacing.  This part of
c       the routine is new and will be revised soon to allow for 
c       multidirectional movement and better velocity control
        
        IF (MSRC.NE.0) THEN
           IF (ITK.LT.BLST) THEN
             ITK = ITK+1
             T(MSRR,ITK) = TMSRC
           END IF
        END IF
        IF (IFCV.NE.0) THEN        
C       'CALL to CONVECTION ROUTINE'      
           CALL CONV
        END IF
        
C      'CALL to CONDUCTION ROUTINE'
       CALL TEMP

       TIME = TIME + TI
       IF (TIME-TP0) 1890,1890,1000
 1000  continue
       DO I = 1,NR
            WRITE(12,2103) (T(I,J), J=1,NC)
       END DO
       
       WRITE (*,*) 'Program finished - Press enter to clear screen'
       CLOSE (12)
       READ (*,*)
C***** SAVE MODEL RESULT TO FILE NAMED SFILE
       CALL SAVIT     
       END
*****************************
       SUBROUTINE TEMP
      
      COMMON /FIRST/ X,Y,NR,NC,T,TI,NOD
      COMMON /SECOND/ D,CHP,DTC,CHF,TIME
       
       REAL T(120,120), W(120,120),TI,TIME,XB,YB,AD,CT
      REAL X(120),Y(120)
      REAL CHP(9),D(9),QFAC,DTC,CHF
      INTEGER NOD(120,120)
      INTEGER KD, KD1, KKD, KKD1, L1, L2, L3, M1, M2, M3, KHP   
      
C      WRITE(*,*)'ENTERED SUBROUTINE TEMP'
c     Alas, a magic number appears here.  See comments below for
c     explanation.  I need to add a data column for density in the
c     conduction codes.  Then this calculation can be done on the fly.
 2720 QFAC = 14.33
C       DHF = CHF * QFAC * TI 
C       Temperature difference due to constant heat flow into the model
c       from the bottom is calculated from DHF.  QFAC is TI (seconds/yr)
c       divided by the product of density and heat capacity.  The value
c       for density is 2200 kg/m^3 and heat capacity is 1000 kJ/kg K.
       DHF = CHF * QFAC *TI
       DO 3760, I = 1, NR
C      WRITE(*,*) 'ROW',I
      DO 3750, J = 1, NC
      XB = 0.
      YB = 0.
*       KKD is the code for thermal conductivity (1-8)
*       NNN is the code for direction of heat flow (1-8)
*       KHP is the code for radioactive heat production (1-8)
*       the following statements determine KKD, KHP and NNN from NOD 
 2780     KKD = NOD(I, J) / 1000
      NNN = NOD(I, J) - KKD * 1000
      KHP = NNN / 100
      NNN = NNN - KHP * 100

c       The following computed GO TO has been kept to facilitate
c       converting the code to FORTRAN 77 if necessary.  When time
c       permits the nested IF's will be converted to SELECT CASE

C 2830  GO TO(2850,3320,2910,2880,2850,2850,2910,2880,2910,
C     + 2880,2850,2910,2880) NNN
      
      IF (NNN.EQ.1) THEN
      J1 = J + 1
      NN = 2
      ELSEIF (NNN.EQ.2) THEN
      GOTO 3320
      ELSEIF (NNN.EQ.3) THEN
      J1 = J + 1
      NN = 1
      ELSEIF (NNN.EQ.4) THEN
      J1 = J - 1
      NN = 1
      ELSEIF (NNN.EQ.5) THEN
      J1 = J + 1
      NN = 2
      ELSEIF (NNN.EQ.6) THEN
      J1 = J + 1
      NN = 2
      ELSEIF (NNN.EQ.7) THEN
      J1 = J + 1
      NN = 1
      ELSEIF (NNN.EQ.8) THEN
      J1 = J - 1
      NN = 1
      ELSEIF (NNN.EQ.9) THEN
      J1 = J + 1
      NN = 1
      ELSEIF (NNN.EQ.10) THEN
      J1 = J - 1
      NN = 1
      ELSEIF (NNN.EQ.11) THEN
      J1 = J + 1
      NN = 2
      ELSEIF (NNN.EQ.12) THEN
      J1 = J + 1
      NN = 1
      ELSEIF (NNN.EQ.13) THEN
      J1 = J - 1
      NN = 1
      END IF

      GO TO 2890

       
       
 2880  J1 = J - 1
       NN = 1

 2890  IF(J1.LE.0) THEN
            GOTO 3290
       ENDIF
       IF(J1.GT.NC) THEN
            GOTO 3290
       ENDIF
*       routines for diagonally partitioned cells - GONE

 2930  KKD1 = NOD(I, J1) / 1000
 2950  KD = KKD
 3120  KD1 = KKD1
 3280  CT = T(I, J1) - T(I, J)
       AD = X(J1) / D(KD1) + X(J) / D(KD)
       XB = XB + 2 * CT / (AD * X(J))
 3290  GO TO (3340,2880) NN      
 3320  W(I, J) = T(I, J)
 3340  GO TO(3360,3750,3360,3360,3420,3390,3420,3420,3390,3390,3390,
     + 3390,3390)NNN
 3360  I1 = I + 1
       NN = 2
       GO TO 3440
 3390  I1 = I - 1
       NN = 1
       GO TO 3440
 3420  I1 = I + 1
       NN = 1
 3440  KKD1 = NOD(I1, J) / 1000

C       WRITE(*,*) 'KKD1 AND KD AT THIS POINT ARE', KKD1,KD

 3510  KD1 = KKD1
 3670  CT = T(I1, J) - T(I, J)
C        write(*,*) 'number  D(KD1)  D(KD)  KD1  KD  NN' 
C     write(*,*) '3670', D(KD1), D(KD), KD1, KD,NN
      IF (KD1.LE.1) THEN KD1 = 1 
      AD = Y(I1) / D(KD1) + Y(I) / D(KD)
C     WRITE(*,*) AD,Y(I1),Y(I)
      YB = YB + 2 * CT / (AD * Y(I))
C     WRITE(*,*) YB,' VALUE FOR YB'
       GO TO(3710,3390) NN
 3710  W(I, J) = T(I, J) + TI * (XB + YB)
       IF (NNN-10) 3740,3740, 3730
 3730  W(I, J) = W(I, J) + DHF / Y(I)
 3740  W(I, J) = W(I, J) + CHP(KHP) * TI / DTC
 3750  CONTINUE
 3760  CONTINUE
     
       DO 3810, I = 1, NR
      DO 3800, J = 1, NC
      T(I, J) = W(I, J)
 3800  CONTINUE  
 3810  CONTINUE
      RETURN
      END
*********************************************************
      SUBROUTINE CONV
      COMMON /FIRST/ X,Y,NR,NC,T,TI,NOD
      COMMON /THIRD/ NODCV,RVSH,FVSH,TPMIN,TIC,VEL,TOTIM,KSAVE,TP0
      REAL W(120,120),T(120,120),FVSH(8),RVSH(8),X(120)
      REAL TPMIN(9),Y(120),VEL(8)
      REAL TI,TA,TIC,XB,YB,CT,AMT,DIST,RATIO
      INTEGER NODCV(120,120), NOD(120,120)
      INTEGER I,J,NR,NC,J1,I1,LCONV,HC,IHCVN,MT,NN,KF,KV,KR,I2
      I2 = 1
      LCONV = INT(TI / (10 * TIC))
      IF (LCONV.GT.5) THEN 
         LCONV = 5
      ENDIF
 5030  IF (LCONV.LE.0) THEN 
          LCONV = 1
       ENDIF
C       WRITE(*,*) LCONV
 5050 TA = TI / REAL(LCONV)
      DO 5850 IHCVN = 1, LCONV
C        HC = 0
C        NCC = 0
C        I2 = I2+1
 5080    DO  I = 1, NR
           DO  J = 1, NC

C               WRITE(*,*) 'ROW',I,'COLUMN',J
                        IF (NODCV(I, J)) 5110, 5110, 5130
 5110                   W(I, J) = T(I, J)
                        GO TO 5700
 5130                   KK = NOD(I, J) / 100
                        KK = NOD(I, J) - 100 * KK
                        IF (KK - 2) 5160, 5110, 5160
 5160                   XB = 0.0
                        YB = 0.0
                        MT = NODCV(I, J) / 10000
                        NN = NODCV(I, J) - 10000 * MT
                        KV = NN / 1000
                        NN = NN - KV * 1000
                        KF = NN / 100
                        NN = NN - KF * 100
                        KR = NN / 10
                        NN = NN - KR * 10
C        WRITE(*,*) NODCV(I,J),MT,KV,KF,KR,NN

c      GO TO(5270,5320,5340,5370,5410,5420,5450,5500,5520)NN
        SELECT CASE (NN)
                CASE (1)
 5270                    I1 = I - 1
                         J1 = J - 1
                CASE (2)
 5320                    I1 = I - 1
                         J1 = J
                CASE (3)
 5340                    I1 = I - 1
 5370                    J1 = J + 1
                CASE (4)
 5420                    I1 = I 
 5450                    J1 = J - 1
                CASE (5)
                        GO TO 5615
                CASE (6)
 5470                   I1 = I 
                        J1 = J + 1
                CASE (7)
 5500                    I1 = I + 1
 5520                    J1 = J - 1
                CASE (8)
 5560                    I1 = I + 1
                         J1 = J
 5570           CASE (9)                 
                         I1 = I + 1
                         J1 = J + 1 
                CASE DEFAULT 
                         GO TO 5615
                END SELECT
                
                IF (I1.GT.NR) GOTO 5615
                IF (J1.GT.NC) GOTO 5615
                IF (I1.LE.0) GOTO 5615
                IF (J1.LE.0) GOTO 5615

 
 5610          IF (T(I1, J1) - TPMIN(MT)) 5615, 5620, 5620
 5615          W(I, J) = T(I, J)
               GO TO 5700
 5620          XB = X(J1) + X(J)
               YB = Y(I1) + Y(I)
               CT = T(I1, J1) - T(I, J)
               AMT = (VEL(KV) * FVSH(KF) / RVSH(KR)) * TA
               DIST = SQRT((XB / 2.)**2 + (YB / 2.)**2)
               RATIO = AMT / DIST
               IF (RATIO-1) 5690,5690, 5670 
 5670          HC = 1
               RATIO = .999999
 5690          W(I, J) = T(I, J) + RATIO * CT
 5700      END DO
 5775   END DO
        DO  I = 1, NR
           DO  J = 1, NC
                 T(I, J) = W(I, J)
 5776      END DO
 5777   END DO   
C   The following was a column counting scheme left over from Chuck 
c   Brott's original FORTRAN 2 code.  I have kept it for no particular
c   reason other than I was never certain exactly what it did.  Thus,
c   when in doubt, don't change anything.
C         IF (HC) 5860, 5860, 5780 
C 5780       IF (NCC-1) 5790, 5820, 5850
C 5790         IF (I2.LE.NC) NCC = 1
C              GOTO 5080
C 5820   NCC = 2
C        GOTO 5860
 5850    CONTINUE
 5860 RETURN
      END
C*********************************       
       SUBROUTINE SAVIT
       COMMON SFILE,TITLE,SUFILE
       COMMON IFCV
       COMMON /FIRST/ X,Y,NR,NC,T,TI,NOD
       COMMON /SECOND/ D,CHP,DTC,CHF,TIME
       COMMON /THIRD/ NODCV,RVSH,FVSH,TPMIN,TIC,VEL,TOTIM,KSAVE,TP0
       
       REAL X(120),Y(120),T(120,120)
       REAL TI,TIME,TPMIN(9),VEL(8),TIC,TM,TN,RNCX,RNRY
       REAL CHP(9),D(9),DTC,CHF,RVSH(8),FVSH(8),TP0
       INTEGER NOD(120,120),NODCV(120,120),IFCV,KSAVE
       
       CHARACTER*80 SFILE,SUFILE
       CHARACTER*80 TITLE

       IF (KSAVE) 1220,1220,5
    5  OPEN (11,FILE=SFILE, STATUS = 'UNKNOWN')
   16  FORMAT(F6.2)
       WRITE (11,12) NR, NC, IFCV
       WRITE (11,11) CHF*1000,TIME
  11   FORMAT(F6.2,F15.3)
  12   FORMAT(3I4)
       WRITE (11,15) TITLE
   15  FORMAT(A80)

       DO I = 1, NR
         WRITE(11,21) (T(I,J),J=1,NC)
       END DO
   21  FORMAT(120F20.10)
     
       DO I = 1, NR
         WRITE(11,22) (NOD(I,J),J=1,NC)
       END DO
   22  FORMAT(120I6)
   
       IF (IFCV.EQ.1) THEN
        DO I = 1,NR
          WRITE(11,22) (NODCV(I,J),J=1,NC)
        END DO
       END IF

       WRITE(11,23) (X(I),I=1,NC)

   51  CONTINUE 
   23  FORMAT(50F10.2)
       WRITE(11,23) (Y(I),I=1,NR)
       WRITE(11,24) ((CHP(I)*1E6),I=1,8)
   24  FORMAT(8G10.2)
       WRITE(11,24) ((D(I)/14.33),I=1,8)
       IF (IFCV.EQ.1) THEN
        WRITE(11,24) (FVSH(I),I=1,8)
        WRITE(11,24) (RVSH(I),I=1,8)
        WRITE(11,24) (TPMIN(I),I=1,8)
        WRITE(11,24) (VEL(I),I=1,8)
       END IF
       CLOSE (11)

 1220 IF (TIME.LT.TP0) THEN              
         RETURN
      ENDIF
      OPEN (12,FILE=SUFILE,STATUS='UNKNOWN')
      WRITE(12,1222) 'DSAA'
 1222 FORMAT(A4)
      WRITE(12,1223) NC,NR
      
      RNCX=X(1)*real(NC)
      RNRY=Y(1)*real(NR)

      IF (X(1).LT.0.01) THEN
        RNCX = RNCX*1000.
      ELSE IF (X(1).LT.0.1) THEN
        RNCX = RNCX*100.
      ELSE IF (X(1).LT.1) THEN
        RNCX = RNCX*10.
      END IF
      IF (Y(1).LT.0.01) THEN
        RNRY = RNRY*1000.
      ELSE IF (Y(1).LT.0.1) THEN
        RNRY = RNRY*100.
      ELSE IF (Y(1).LT.1) THEN
        RNRY = RNRY*10.
      END IF
      WRITE(12,1224) 0.0,RNCX
      WRITE(12,1224) -RNRY,0.0
 1223 FORMAT(2I4)
      TN = 1000.0
      TM = -1000.0
      DO I = 1,NR
       DO J = 1,NC
        IF (T(I,J).GT.TM) TM=T(I,J)
        IF (TN.GT.T(I,J)) TN=T(I,J)
       END DO
      END DO
      WRITE(12,1224) TN,TM
 1224 FORMAT(2F10.2)
       DO I =  NR,1,-1
         WRITE(12,21) (T(I,J),J=1,NC)
       END DO

      RETURN
      END
