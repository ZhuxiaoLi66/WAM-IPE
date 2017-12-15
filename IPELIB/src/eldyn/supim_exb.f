!20111031: originally copied from supim2bc.f
!*
!      program test_supim_exb
      SUBROUTINE supim_EXBV(utime,IFL,LT_SEC,r_apex,GLON_deg,VEXBup_out)
      USE module_precision
      USE module_input_parameters,ONLY:NYEAR,NDAY,F107AV_ipe,F107D_ipe
     &,start_time
      IMPLICIT NONE       
!
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec),INTENT(IN) :: IFL
      REAL(KIND=real_prec),INTENT(IN) :: LT_SEC
      REAL(KIND=real_prec),INTENT(IN) :: r_apex
      REAL(KIND=real_prec),INTENT(IN) :: GLON_deg
      REAL(KIND=real_prec),INTENT(OUT) :: VEXBup_out
! local
      INTEGER,PARAMETER :: NPTS=601
      INTEGER,PARAMETER :: NFL=170
!      REAL,DIMENSION(NFL) :: GREQ
      INTEGER :: IYEAR,IDAY,F107A,F107
      REAL :: AP
      REAL :: PYE,R0,RAD
      INTEGER :: JDJ,JDA
      REAL,DIMENSION(NFL) :: VPE, VZONE
      COMMON
     +/BLK1/PYE,R0,RAD
     +/BLK9/VPE
!,VP(NPTS,NFL),DVP(NPTS,NFL)
     +     ,VZONE
!,VZON(NPTS,NFL)
     +/BLK15/IYEAR,IDAY,F107A,F107,AP

!parameters
      PYE=4.*ATAN(1.)
      RAD=PYE/180.
      R0=6370.E3 !meter

! input parameters
      IYEAR=NYEAR  !1998
      IDAY=NDAY    !74
      F107A=F107AV_ipe !187.
      F107=F107D_ipe   !184.
      AP=4.
! LOW ALTITUDE DRIFT PATTERN
      JDJ=10
! HIGH ALTITUDE DRIFT PATTERN
      JDA=1

!      DT=900. ![sec]
!*  MAIN CALCULATIONS
!      time_loop:  DO i=1,97
!      ZTSEC=REAL(i-1)*DT !ZTSEC+DT
      write(3000,*) LT_SEC, LT_SEC/3600.

!      dht=300.E3
!      ht0=200.E3
!      ifl_loop: DO IFL=1,14 !NFL  !410

!        INEB=IEBN(IFL)
!        ISEB=IEBS(IFL)
!      INEB=1
!      ISEB=NPTS

! APEX height: GREQ(IFL) meter
!      ht=ht0 + dht*REAL(IFL-1)
!      GREQ(IFL)=R0+ht
!      if(i==1) write(3001,*)'ht',ht
!GLOND(NM,IFL) degrees
!      GLOND_IFL=288.
      if(utime==start_time)
     &  write(3001,*)IFL,' GLOND',GLON_deg,' GREQ',r_apex

      CALL EXBV(JDJ,JDA,LT_SEC,r_apex,GLON_deg
!     &              ,INEB,ISEB
     &,IFL)
      write(3000,fmt=*) IFL,VPE(IFL)
      VEXBup_out=VPE(IFL)
!      END DO ifl_loop !410 IFL=1,2 !NFL
!      END DO time_loop !
!      end program test_supim_exb
      END SUBROUTINE supim_EXBV
!*-------------------------------------------------------------------------------
!*
      SUBROUTINE EXBV(IDJ,IDA,ZTSEC,GRE,GLONB
!     &,INEB,ISEB
     &,IFL)
!*
!********************************************************************************
!*  ROUTINE TO EVALUATE THE GEOMAGNETIC EQUATORIAL EXB DRIFT VELOCITY
!*  AND THE E-REGION ZONAL DRIFT VELOCITY
!*
!*  DRIFT ZERO                                       APEX ALTITUDE < ZE0
!*  LOW ALTITUDE DRIFT (INPUT, EG, JICAMARCA)  ZE1 < APEX ALTITUDE < ZE2
!*  MID ALTITUDE DRIFT (EG, ARECIBO)           ZE3 < APEX ALTITUDE < ZE4
!*  HIGH ALTITUDE DRIFT (ZERO PROBABLY)        ZE5 < APEX ALTITUDE
!*  INTERPOLATION AT INTERMEDIATE APEX ALTITUDES
!********************************************************************************
!*
      PARAMETER(NPTS=601,NFL=170)
      SAVE RZE0,RZE1,RZE2,RZE3,RZE4,RZE5,JVPE,AVPE
      REAL JVPE
      COMMON
     +/BLK1/PYE,R0,RAD
     +/BLK9/VPE(NFL)
!,VP(NPTS,NFL),DVP(NPTS,NFL)
     +     ,VZONE(NFL)
!,VZON(NPTS,NFL)
     +/BLK15/IYEAR,IDAY,F107A,F107,AP
     +/EXBVEL/EXBJT(97,10),EXBJD(97,10),EXBAT(97,1),EXBAD(97,1)
      DATA ZE0/150.E3/,ZE1/200.E3/,ZE2/600.E3/,ZE3/2000.E3/
     &    ,ZE4/3000.E3/,ZE5/4000.E3/
!*
      IF(IFL.GT.1) GOTO 40
!*
      RZE0=R0+ZE0+1.E3
      RZE1=R0+ZE1+1.E3
      RZE2=R0+ZE2+1.E3
      RZE3=R0+ZE3+1.E3
      RZE4=R0+ZE4+1.E3
      RZE5=R0+ZE5+1.E3
!*
!*     LOW-ALTITUDE DRIFT
!*
      IF(IDJ.EQ.0) THEN
!* NO DRIFT
         JVPE=0.
      ELSEIF(IDJ.GE.1.AND.IDJ.LE.9) THEN
!* DRIFT FROM BLOCK DATA EXBVD
         ZTSEC1=ZTSEC/3600.
         IF(ZTSEC1.LT.EXBJT(1,IDJ)) ZTSEC1=ZTSEC1+24.
         II=0
   20    II=II+1
         IJ=II+1
         IF(ABS(EXBJT(II,IDJ)-ZTSEC1).LT.0.001) THEN
            JVPE=EXBJD(II,IDJ)
         ELSE
            IF(EXBJT(II,IDJ).LT.ZTSEC1.AND.EXBJT(IJ,IDJ).GT.ZTSEC1) THEN
               W=(EXBJD(IJ,IDJ)-EXBJD(II,IDJ))
     $          /(EXBJT(IJ,IDJ)-EXBJT(II,IDJ))
               JVPE=EXBJD(II,IDJ)+W*(ZTSEC1-EXBJT(II,IDJ))
            ELSE
               GOTO 20
            ENDIF
         ENDIF
      ELSEIF(IDJ.EQ.10) THEN
!*  DRIFT FROM SCHERLIESS AND FEJER (JGR, 104, 6828-6842, 1999)
         SLT=ZTSEC/3600.
         CALL SFVDM(SLT,GLONB,IDAY,F107,JVPE)
      ELSE
         STOP 'INCORRECT LOW-ALTITUDE DRIFT PARAMETER'
      ENDIF
!*
!*     MID-ALTITUDE DRIFT
!*
      IF(IDA.EQ.0) THEN
!* NO DRIFT
         AVPE=0.
      ELSEIF(IDA.EQ.1) THEN
!* DRIFT FROM BLOCK DATA EXBVD
         ZTSEC1=ZTSEC/3600.
         IF(ZTSEC1.LT.EXBAT(1,IDA)) ZTSEC1=ZTSEC1+24.
         II=0
   30    II=II+1
         IJ=II+1
         IF(ABS(EXBAT(II,IDA)-ZTSEC1).LT.0.001) THEN
            AVPE=EXBAD(II,IDA)
         ELSE
            IF(EXBAT(II,IDA).LT.ZTSEC1.AND.EXBAT(IJ,IDA).GT.ZTSEC1) THEN
               W=(EXBAD(IJ,IDA)-EXBAD(II,IDA))
     $          /(EXBAT(IJ,IDA)-EXBAT(II,IDA))
               AVPE=EXBAD(II,IDA)+W*(ZTSEC1-EXBAT(II,IDA))
            ELSE
               GOTO 30
            ENDIF
         ENDIF
      ELSE
         STOP 'INCORRECT MID-ALTITUDE DRIFT PARAMETER'
      ENDIF
!*
!*     MULTIPLIED BY 2 TO GIVE APEX VALUE
         AVPE=AVPE*2.
!*
!*     CALCULATION OF VPE
!*
   40 IF(IFL.EQ.1.OR.IFL.EQ.NFL) THEN
!*        BOTH VPE(1) AND VPE(NFL) MUST BE ZERO
         VPE(IFL)=0.
      ELSEIF(GRE.LT.RZE0) THEN
         VPE(IFL)=0.
      ELSEIF(GRE.LT.RZE1) THEN
         VPE(IFL)=JVPE*(GRE-RZE0)/(RZE1-RZE0)
      ELSEIF(GRE.LT.RZE2) THEN
         VPE(IFL)=JVPE
      ELSEIF(GRE.LT.RZE3) THEN
         VPE(IFL)=JVPE+(AVPE-JVPE)*(GRE-RZE2)/(RZE3-RZE2)
      ELSEIF(GRE.LT.RZE4) THEN
         VPE(IFL)=AVPE
      ELSEIF(GRE.LT.RZE5) THEN
         VPE(IFL)=AVPE+(0.-AVPE)*(GRE-RZE4)/(RZE5-RZE4)
      ELSE
         VPE(IFL)=0.
      ENDIF

!*
!*  E-REGION EQUATORIAL ZONAL DRIFT (POSITIVE IN EASTWARD DIRECTION)
!*
      VZONE(IFL)=0.
!*
      RETURN
!*
      END
!*
!*-------------------------------------------------------------------------------
!*
      BLOCK DATA EXBVD
!*
!*  EXB DRIFT VELOCITY DATA FOR SUBROUTINE EXBV
!*
      COMMON
     +/EXBVEL/EXBJT(97,10),EXBJD(97,10),EXBAT(97,1),EXBAD(97,1)
!*
!*  LOW-ALTITUDE DRIFT PATTERNS
!*
!*DRIFT 1: EXB DRIFT AT JICAMARCA, EQUINOX, F107<100 (FEJER ET AL., 1991)
      DATA (EXBJT(I,1),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,1),I=51,97)/47*-1./
      DATA (EXBJD(I,1),I=1,50)
     $/  7.17,  9.10, 13.41, 15.47, 19.43, 21.15, 21.50, 22.55, 21.03
     $, 20.22, 17.86, 13.51, 10.76,  8.87,  6.71,  4.95,  2.42,  1.29
     $,  0.72,  0.73,  3.60, 10.11, 10.17,  2.76, -7.27,-13.50,-13.68
     $,-16.19,-15.53,-13.66,-12.73,-13.64,-14.17,-13.09,-15.09,-12.91
     $,-16.11,-14.12,-13.76,-14.40,-10.44, -7.70, -9.44, -3.61, -2.86
     $, -5.43, -3.77, -0.63,  5.23,  7.17/
     $,(EXBJD(I,1),I=51,97)/47*-1./
!*
!*DRIFT 2: EXB DRIFT AT JICAMARCA, EQUINOX, 100<F107<150 (FEJER ET AL., 1991)
      DATA (EXBJT(I,2),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,2),I=51,97)/47*-1./
      DATA (EXBJD(I,2),I=1,50)
     $/ 11.49, 13.33, 17.81, 22.99, 23.62, 23.75, 25.45, 23.04, 24.00
     $, 20.03, 17.03, 14.60, 11.81, 10.08,  8.69,  8.23,  6.98,  7.09
     $,  6.71,  8.38, 12.95, 23.25, 25.80,  5.55, -8.94,-20.70,-20.81
     $,-24.93,-26.36,-25.44,-25.99,-22.68,-19.75,-21.10,-19.74,-17.22
     $,-16.25,-15.77,-15.31,-13.87,-13.01,-11.37,-10.25,-10.91, -6.86
     $, -6.10, -1.23,  3.83,  9.65, 11.49/
     $,(EXBJD(I,2),I=51,97)/47*-1./
!*
!*DRIFT 3: EXB DRIFT AT JICAMARCA, EQUINOX, F107>150 (FEJER ET AL., 1991)
      DATA (EXBJT(I,3),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,3),I=51,97)/47*-1./
      DATA (EXBJD(I,3),I=1,50)
     $/  8.24, 10.18, 15.00, 18.62, 21.78, 24.95, 26.20, 26.72, 25.72
     $, 23.89, 22.39, 20.59, 19.96, 17.02, 16.33, 14.70, 13.22, 12.37
     $, 12.68, 15.60, 19.75, 33.44, 50.75, 38.96,  8.25,-22.53,-33.88
     $,-33.35,-32.39,-28.99,-28.14,-29.08,-28.32,-24.57,-25.76,-22.73
     $,-23.15,-21.74,-23.73,-21.53,-21.58,-25.09,-24.27,-21.01,-13.94
     $, -8.46, -2.79,  1.83,  6.32,  8.24/
     $,(EXBJD(I,3),I=51,97)/47*-1./
!*
!*DRIFT 4: EXB DRIFT AT JICAMARCA, WINTER, F107<100 (FEJER ET AL., 1991)
      DATA (EXBJT(I,4),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,4),I=51,97)/47*-1./
      DATA (EXBJD(I,4),I=1,50)
     $/  7.43, 11.71, 14.10, 15.56, 15.87, 19.57, 20.17, 19.55, 20.42
     $, 19.66, 19.29, 17.37, 15.72, 14.37, 12.61, 11.37,  8.83,  7.42
     $,  5.02,  3.24, -0.25, -2.34, -5.50,-11.05,-15.24,-16.50,-17.48
     $,-16.00,-16.23,-14.70,-15.31,-16.70,-15.58,-12.92,-11.56, -7.42
     $, -8.18, -9.87, -5.67, -4.33, -7.00, -7.53, -7.56, -7.40,-12.74
     $, -6.15, -0.42, -3.34,  3.14,  7.43/
     $,(EXBJD(I,4),I=51,97)/47*-1./
!*
!*DRIFT 5: EXB DRIFT AT JICAMARCA, WINTER, 100<F107<150 (FEJER ET AL., 1991)
      DATA (EXBJT(I,5),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,5),I=51,97)/47*-1./
      DATA (EXBJD(I,5),I=1,50)
     $/ 13.97, 16.25, 17.23, 18.61, 16.33, 18.21, 19.27, 18.12, 19.26
     $, 18.28, 17.75, 16.43, 15.20, 14.65, 11.92,  9.89,  9.82,  8.12
     $,  7.30,  7.38,  8.53, 11.61,  2.13, -6.27,-16.23,-18.60,-20.65
     $,-23.17,-25.46,-22.53,-21.49,-21.22,-20.16,-19.00,-18.56,-17.79
     $,-18.57,-17.79,-15.21,-14.30,-11.97, -8.29, -3.79, -2.34,  1.54
     $,  9.97,  8.60,  9.11, 11.68, 13.97/
     $,(EXBJD(I,5),I=51,97)/47*-1./
!*
!*DRIFT 6: EXB DRIFT AT JICAMARCA, WINTER, F107>150 (FEJER ET AL., 1991)
      DATA (EXBJT(I,6),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,6),I=51,97)/47*-1./
      DATA (EXBJD(I,6),I=1,50)
     $/ 11.62, 13.25, 15.92, 17.59, 19.45, 20.52, 21.70, 21.64, 21.67
     $, 21.29, 18.76, 18.83, 18.06, 17.27, 16.71, 15.13, 13.22, 12.76
     $, 12.19, 14.22, 17.03, 20.23, 11.34, -6.35,-15.32,-20.46,-23.60
     $,-26.24,-26.38,-26.34,-25.42,-25.90,-23.97,-23.18,-21.33,-19.64
     $,-19.35,-19.30,-16.56,-16.30,-16.44,-11.16,-12.06,-12.54, -4.08
     $,  7.78,  6.64,  8.57,  9.98, 11.62/
     $,(EXBJD(I,6),I=51,97)/47*-1./
!*
!*DRIFT 7: EXB DRIFT AT JICAMARCA, SUMMER, F107<100 (FEJER ET AL., 1991)
      DATA (EXBJT(I,7),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,7),I=51,97)/47*-1./
      DATA (EXBJD(I,7),I=1,50)
     $/  7.35,  9.55, 12.02, 14.59, 17.76, 18.28, 18.40, 15.60, 11.04
     $,  9.50,  6.17,  3.49,  2.30,  1.42, -1.01, -0.85, -1.05,  0.21
     $,  1.82,  1.99,  1.67,  1.55,  4.23,  3.56,  0.08, -1.90, -2.13
     $, -2.94, -5.79, -9.96, -9.49,-19.39,-17.31,-16.34,-11.95, -6.97
     $, -8.56,-10.60, -5.77, -7.82, -4.18, -8.54, -7.70, -4.58, -4.35
     $,  4.03,  1.19,  3.79,  5.15,  7.35/
     $,(EXBJD(I,7),I=51,97)/47*-1./
!*
!*DRIFT 8: EXB DRIFT AT JICAMARCA, SUMMER, 100<F107<150 (FEJER ET AL., 1991)
      DATA (EXBJT(I,8),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,8),I=51,97)/47*-1./
      DATA (EXBJD(I,8),I=1,50)
     $/  8.32,  9.71,  9.69, 11.15, 14.42, 14.31, 12.77, 13.82, 11.66
     $, 11.63,  9.25,  9.17,  9.29,  8.63,  5.82,  4.29,  4.95,  6.46
     $,  9.45,  9.74, 12.85, 15.79, 28.10, 25.50, 16.73,  7.78,  1.48
     $, -6.48,-12.02,-20.74,-18.56,-16.40,-17.73,-19.05,-14.32,-16.60
     $,-20.27,-21.33,-22.05,-22.28,-22.44,-23.38,-20.06,-12.42, -2.86
     $,  0.77,  2.19,  4.67,  6.93,  8.32/
     $,(EXBJD(I,8),I=51,97)/47*-1./
!*
!*DRIFT 9: EXB DRIFT AT JICAMARCA, SUMMER, F107>150 (FEJER ET AL., 1991)
      DATA (EXBJT(I,9),I=1,50)
     $/  8.00,  8.25,  8.75,  9.25,  9.75, 10.25, 10.75, 11.25, 11.75
     $, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25
     $, 16.75, 17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25, 20.75
     $, 21.25, 21.75, 22.25, 22.75, 23.25, 23.75, 24.25, 24.75, 25.25
     $, 25.75, 26.25, 26.75, 27.25, 27.75, 28.25, 28.75, 29.25, 29.75
     $, 30.25, 30.75, 31.25, 31.75, 32.00/
     $,(EXBJT(I,9),I=51,97)/47*-1./
      DATA (EXBJD(I,9),I=1,50)
     $/  8.11, 11.57, 15.98, 13.86, 16.79, 16.01, 15.92, 14.20, 12.95
     $, 11.16,  8.86,  8.93,  9.92,  7.11,  6.04,  4.01,  6.09,  5.99
     $,  5.96,  6.50,  9.49, 14.25, 23.02, 33.77, 25.17,  8.48, -0.51
     $, -2.33, -7.96,-14.89,-14.13,-11.73,-12.57,-13.68,-16.13,-19.55
     $,-17.94,-27.96,-30.04,-32.28,-35.30,-28.51,-19.88,-13.78, -7.90
     $,  2.75,  3.60,  4.26,  4.65,  8.11/
     $,(EXBJD(I,9),I=51,97)/47*-1./
!*
!*  MID-ALTITUDE DRIFT PATTERNS
!*
!*DRIFT 1: EXB DRIFT AT ARECIBO (FEJER, 1993)
      DATA (EXBAT(I,1),I=1,49)
     $/  0.00,  0.50,  1.00,  1.50,  2.00,  2.50,  3.00,  3.50,  4.00
     $,  4.50,  5.00,  5.50,  6.00,  6.50,  7.00,  7.50,  8.00,  8.50
     $,  9.00,  9.50, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00
     $, 13.50, 14.00, 14.50, 15.00, 15.50, 16.00, 16.50, 17.00, 17.50
     $, 18.00, 18.50, 19.00, 19.50, 20.00, 20.50, 21.00, 21.50, 22.00
     $, 22.50, 23.00, 23.50, 24.00/
     $,(EXBAT(I,1),I=50,97)/48*-1./
      DATA (EXBAD(I,1),I=1,49)
     $/ -8.00, -8.00, -8.00, -8.00, -7.00, -5.00, -4.00, -1.00,  3.00
     $,  5.00,  5.00,  2.00, -4.00,-10.00,-10.00, -6.00,  3.00,  9.00
     $, 12.00, 16.00, 18.00, 19.00, 20.00, 20.00, 18.00, 12.00, 10.00
     $,  8.00,  4.00,  1.00, -2.00, -3.00, -4.00, -5.00, -4.00, -3.00
     $, -2.00, -2.00, -2.00, -3.00, -4.00, -7.00,-10.00,-11.00,-12.00
     $,-12.00,-11.00, -9.00, -8.00/
     $,(EXBAD(I,1),I=50,97)/48*-1./
!*
      END
!*-------------------------------------------------------------------------------
!*
      SUBROUTINE SFVDM(SLT,GL,IDAY,F107,VPE)
!*
!********************************************************************************
!*
!*     ROUTINE TO DETERMINE VERTICAL EXB DRIFT
!*     (SCHERLIESS, L., AND B.G. FEJER,
!*      J. GEOPHYS. RES., 104, 6829-6842, 1999)
!*
!********************************************************************************
!*
!*    SLT:   SOLAR LOCAL TIME
!*    GL:    GEOGRAPHIC LONGITUDE (+ EAST)
!*    IDAY:  DAY OF YEAR
!*    F107:  F10.7
!*    VPE:   EQUATORIAL VERTICAL DRIFT
!*
!********************************************************************************
!*
      DIMENSION COEFF(624),FUNCT(6)

      DATA (COEFF(I),I=1,60)/
     $ -10.80592, -9.63722,-11.52666, -0.05716, -0.06288,  0.03564,
     $  -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138,
     $   2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171,
     $ -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656,
     $  -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419,
     $ -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984,
     $ -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062,
     $ -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297,
     $   1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023,
     $   5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166/
      DATA (COEFF(I),I=61,120)/
     $   3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456,
     $   7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008,
     $  -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507,
     $  -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477,
     $ -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528,
     $ -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075,
     $  14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723,
     $  12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484,
     $  18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249,
     $   4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394/
      DATA (COEFF(I),I=121,180)/
     $  14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739,
     $   7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626,
     $   7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039,
     $  10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177,
     $  21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715,
     $  19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434,
     $  26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631,
     $  21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903,
     $  28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622,
     $  22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206/
      DATA (COEFF(I),I=181,240)/
     $  31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322,
     $  46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982,
     $  13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970,
     $  21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589,
     $  16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570,
     $  18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991,
     $  10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775,
     $  12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851,
     $  -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433,
     $ -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385/
      DATA (COEFF(I),I=241,300)/
     $   2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932,
     $   3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084,
     $  17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516,
     $   0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331,
     $  15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934,
     $   4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308,
     $   9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124,
     $  13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465,
     $   5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222,
     $   9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418/
      DATA (COEFF(I),I=301,360)/
     $   9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919,
     $  13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801,
     $  10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208,
     $  10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751,
     $   6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705,
     $   5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151,
     $  29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325,
     $  17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398,
     $   8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374,
     $  -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187/
      DATA (COEFF(I),I=361,420)/
     $   8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156,
     $  14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061,
     $  14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129,
     $   6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993,
     $   7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631,
     $  -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188,
     $  -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214,
     $ -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655,
     $   2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650,
     $   5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376/
      DATA (COEFF(I),I=421,480)/
     $  13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223,
     $  -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585,
     $ -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463,
     $ -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679,
     $ -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355,
     $ -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593,
     $ -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839,
     $ -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226,
     $ -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004,
     $ -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699/
      DATA (COEFF(I),I=481,540)/
     $ -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796,
     $ -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959,
     $ -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792,
     $ -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535,
     $ -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143,
     $ -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688,
     $  -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422,
     $  -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978,
     $ -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144,
     $ -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133/
      DATA (COEFF(I),I=541,600)/
     $ -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949,
     $ -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583,
     $ -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213,
     $ -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053,
     $ -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085,
     $ -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726,
     $ -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529,
     $ -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544,
     $  -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807,
     $ -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214/
      DATA (COEFF(I),I=601,624)/
     $ -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652,
     $ -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235,
     $ -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214,
     $ -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/
!*
      CALL G(IDAY,F107,FUNCT,GL)
!*
      VPE=0.
      DO 11 I=1,13
        DO 12 J=1,8
          KK=8*(I-1)+J
            DO 13 K=1,6
              IND=6*(KK-1)+K
              BSPL4=BSPL4T(I,SLT)*BSPL4L(J,GL)
              VPE=VPE+BSPL4*FUNCT(K)*COEFF(IND)
   13       CONTINUE
   12   CONTINUE
   11 CONTINUE
!*
      RETURN
!*
      END
!*
!*-------------------------------------------------------------------------------
!*
      FUNCTION BSPL4T(I,T1)
!*
      DIMENSION TT(0:39),B(20,20)
!*
      DATA TT/ 0.00, 2.75, 4.75, 5.50, 6.25,
     *         7.25,10.00,14.00,17.25,18.00,
     *        18.75,19.75,21.00,24.00,26.75,
     *        28.75,29.50,30.25,31.25,34.00,
     *        38.00,41.25,42.00,42.75,43.75,
     *        45.00,48.00,50.75,52.75,53.50,
     *        54.25,55.25,58.00,62.00,65.25,
     *        66.00,66.75,67.75,69.00,72.00/
!*
      T=T1
      IF(I.GE.0.AND.T.LT.TT(I)) THEN
         T=T+24.
      ENDIF
      DO 11 J=I,I+4-1
        IF(T.GE.TT(J).AND.T.LT.TT(J+1)) THEN
           B(J,1)=1.
        ELSE
           B(J,1)=0.
        ENDIF
   11 CONTINUE
      DO 12 J=2,4
        DO 13 K=I,I+4-J
          B(K,J)=(T-TT(K))/(TT(K+J-1)-TT(K))*B(K,J-1)
          B(K,J)=B(K,J)+(TT(K+J)-T)/(TT(K+J)-TT(K+1))*B(K+1,J-1)
   13   CONTINUE
   12 CONTINUE
!*
      BSPL4T=B(I,4)
!*
      RETURN
!*
      END
!*
!*-------------------------------------------------------------------------------
!*
      FUNCTION BSPL4L(I,T1)
!*
      DIMENSION TL(0:24),B(20,20)
!*
      DATA TL/  0, 10,100,190,200,250,280,310,
     *        360,370,460,550,560,610,640,670,
     *        720,730,820,910,920,970,1000,1030,1080/
!*
      T=T1
      IF(I.GE.0.AND.T.LT.TL(I)) THEN
         T=T+360.
      ENDIF
      DO 11 J=I,I+4-1
        IF(T.GE.TL(J).AND.T.LT.TL(J+1)) THEN
           B(J,1)=1.
        ELSE
           B(J,1)=0.
        ENDIF
   11 CONTINUE
!*
      DO 12 J=2,4
        DO 13 K=I,I+4-J
          B(K,J)=(T-TL(K))/(TL(K+J-1)-TL(K))*B(K,J-1)
          B(K,J)=B(K,J)+(TL(K+J)-T)/(TL(K+J)-TL(K+1))*B(K+1,J-1)
   13   CONTINUE
   12 CONTINUE
!*
      BSPL4L=B(I,4)
!*
      RETURN
!*
      END
!*
!*-------------------------------------------------------------------------------
!*
      SUBROUTINE G(IDAY,F107,FUNCT,T)
!*
      DIMENSION FUNCT(6)
!*
      IF(F107.LT.75.) THEN
         FLUX=75.
      ELSEIF(F107.GT.230.) THEN
         FLUX=230.
      ELSE
         FLUX=F107
      ENDIF
      CFLUX=FLUX
!*
      A=0.
      IF(IDAY.GE.120.AND.IDAY.LE.240) THEN
          A=170.
          SIGMA=60.
      ENDIF
      IF(IDAY.GE.60.OR.IDAY.GE.300.) THEN
          A=170.
          SIGMA=40.
      ENDIF
      IF(FLUX.LT.95..AND.A.GT.0.) THEN
         GAUSS=EXP(-0.5*((T-A)**2)/SIGMA**2)
         CFLUX=GAUSS*95.+(1.-GAUSS)*FLUX
      ENDIF
!*
      DO 11 I=1,6
        FUNCT(I)=0.
   11 CONTINUE
!*
      IF(IDAY.GE.135.AND.IDAY.LE.230) THEN
          FUNCT(1)=1.
      ENDIF
      IF(IDAY.LE.45.OR.IDAY.GE.320) THEN
          FUNCT(2)=1.
      ENDIF
      IF(IDAY.GE.75.AND.IDAY.LE.105) THEN
          FUNCT(3)=1.
      ENDIF
      IF(IDAY.GE.260.AND.IDAY.LE.290) THEN
          FUNCT(3)=1.
      ENDIF
!*
      IF(IDAY.GE.45.AND.IDAY.LE.75) THEN    ! W-E
          FUNCT(2)=1.-(IDAY-45)/30.
          FUNCT(3)=1.-FUNCT(2)
      ENDIF
      IF(IDAY.GE.105.AND.IDAY.LE.135) THEN  ! E-S
          FUNCT(3)=1.-(IDAY-105)/30.
          FUNCT(1)=1.-FUNCT(3)
      ENDIF
      IF(IDAY.GE.230.AND.IDAY.LE.260) THEN  ! S-E
          FUNCT(1)=1.-(IDAY-230)/30.
          FUNCT(3)=1.-FUNCT(1)
      ENDIF
      IF(IDAY.GE.290.AND.IDAY.LE.320) THEN  ! E-W
          FUNCT(3)=1.-(IDAY-290)/30.
          FUNCT(2)=1.-FUNCT(3)
      ENDIF
!*
      FUNCT(4)=(CFLUX-140.)*FUNCT(1)
      FUNCT(5)=(CFLUX-140.)*FUNCT(2)
      FUNCT(6)=(CFLUX-140.)*FUNCT(3)
!*
      RETURN
!*
      END
