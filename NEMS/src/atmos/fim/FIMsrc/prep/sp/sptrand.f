C-----------------------------------------------------------------------
      SUBROUTINE SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &                   IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,
     &                   JBEG,JEND,JCPU,
     &                   WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  SPTRAND    PERFORM A GRADIENT SPHERICAL TRANSFORM
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-02-29
C
C ABSTRACT: THIS SUBPROGRAM PERFORMS A SPHERICAL TRANSFORM
C           BETWEEN SPECTRAL COEFFICIENTS OF SCALAR FIELDS
C           AND THEIR MEANS AND GRADIENTS ON A GLOBAL CYLINDRICAL GRID.
C           THE WAVE-SPACE CAN BE EITHER TRIANGULAR OR RHOMBOIDAL.
C           THE GRID-SPACE CAN BE EITHER AN EQUALLY-SPACED GRID
C           (WITH OR WITHOUT POLE POINTS) OR A GAUSSIAN GRID.
C           THE WAVE AND GRID FIELDS MAY HAVE GENERAL INDEXING,
C           BUT EACH WAVE FIELD IS IN SEQUENTIAL 'IBM ORDER',
C           I.E. WITH ZONAL WAVENUMBER AS THE SLOWER INDEX.
C           TRANSFORMS ARE DONE IN LATITUDE PAIRS FOR EFFICIENCY;
C           THUS GRID ARRAYS FOR EACH HEMISPHERE MUST BE PASSED.
C           IF SO REQUESTED, JUST A SUBSET OF THE LATITUDE PAIRS
C           MAY BE TRANSFORMED IN EACH INVOCATION OF THE SUBPROGRAM.
C           THE TRANSFORMS ARE ALL MULTIPROCESSED OVER LATITUDE EXCEPT
C           THE TRANSFORM FROM FOURIER TO SPECTRAL IS MULTIPROCESSED
C           OVER ZONAL WAVENUMBER TO ENSURE REPRODUCIBILITY.
C           TRANSFORM SEVERAL FIELDS AT A TIME TO IMPROVE VECTORIZATION.
C           SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C PROGRAM HISTORY LOG:
C   96-02-29  IREDELL
C 1998-12-15  IREDELL  OPENMP DIRECTIVES INSERTED
C
C USAGE:    CALL SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
C    &                   IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,
C    &                   JBEG,JEND,JCPU,
C    &                   WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)
C   INPUT ARGUMENTS:
C     IROMB    - INTEGER SPECTRAL DOMAIN SHAPE
C                (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C     MAXWV    - INTEGER SPECTRAL TRUNCATION
C     IDRT     - INTEGER GRID IDENTIFIER
C                (IDRT=4 FOR GAUSSIAN GRID,
C                 IDRT=0 FOR EQUALLY-SPACED GRID INCLUDING POLES,
C                 IDRT=256 FOR EQUALLY-SPACED GRID EXCLUDING POLES)
C     IMAX     - INTEGER EVEN NUMBER OF LONGITUDES.
C     JMAX     - INTEGER NUMBER OF LATITUDES.
C     KMAX     - INTEGER NUMBER OF FIELDS TO TRANSFORM.
C     IPRIME   - INTEGER LONGITUDE INDEX FOR THE PRIME MERIDIAN.
C                (DEFAULTS TO 1 IF IPRIME=0)
C     ISKIP    - INTEGER SKIP NUMBER BETWEEN LONGITUDES
C                (DEFAULTS TO 1 IF ISKIP=0)
C     JNSKIP   - INTEGER SKIP NUMBER BETWEEN N.H. LATITUDES FROM NORTH
C                (DEFAULTS TO IMAX IF JNSKIP=0)
C     JSSKIP   - INTEGER SKIP NUMBER BETWEEN S.H. LATITUDES FROM SOUTH
C                (DEFAULTS TO -IMAX IF JSSKIP=0)
C     KWSKIP   - INTEGER SKIP NUMBER BETWEEN WAVE FIELDS
C                (DEFAULTS TO (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C     KGSKIP   - INTEGER SKIP NUMBER BETWEEN GRID FIELDS
C                (DEFAULTS TO IMAX*JMAX IF KGSKIP=0)
C     JBEG     - INTEGER LATITUDE INDEX (FROM POLE) TO BEGIN TRANSFORM
C                (DEFAULTS TO 1 IF JBEG=0)
C                (IF JBEG=0 AND IDIR<0, WAVE IS ZEROED BEFORE TRANSFORM)
C     JEND     - INTEGER LATITUDE INDEX (FROM POLE) TO END TRANSFORM
C                (DEFAULTS TO (JMAX+1)/2 IF JEND=0)
C     JCPU     - INTEGER NUMBER OF CPUS OVER WHICH TO MULTIPROCESS
C     WAVE     - REAL (*) WAVE FIELDS IF IDIR>0
C     GRIDMN   - REAL (KMAX) GLOBAL MEANS IF IDIR<0
C     GRIDXN   - REAL (*) N.H. X-GRADIENTS (STARTING AT JBEG) IF IDIR<0
C     GRIDXS   - REAL (*) S.H. X-GRADIENTS (STARTING AT JBEG) IF IDIR<0
C     GRIDYN   - REAL (*) N.H. Y-GRADIENTS (STARTING AT JBEG) IF IDIR<0
C     GRIDYS   - REAL (*) S.H. Y-GRADIENTS (STARTING AT JBEG) IF IDIR<0
C     IDIR     - INTEGER TRANSFORM FLAG
C                (IDIR>0 FOR WAVE TO GRID, IDIR<0 FOR GRID TO WAVE)
C   OUTPUT ARGUMENTS:
C     WAVE     - REAL (*) WAVE FIELDS IF IDIR<0
C     GRIDMN   - REAL (KMAX) GLOBAL MEANS IF IDIR>0
C     GRIDXN   - REAL (*) N.H. X-GRADIENTS (STARTING AT JBEG) IF IDIR>0
C     GRIDXS   - REAL (*) S.H. X-GRADIENTS (STARTING AT JBEG) IF IDIR>0
C                [GRIDX=(D(WAVE)/DLAM)/(CLAT*RERTH)]
C     GRIDYN   - REAL (*) N.H. Y-GRADIENTS (STARTING AT JBEG) IF IDIR>0
C     GRIDYS   - REAL (*) S.H. Y-GRADIENTS (STARTING AT JBEG) IF IDIR>0
C                [GRIDY=(D(WAVE)/DPHI)/RERTH]
C
C SUBPROGRAMS CALLED:
C   SPWGET       GET WAVE-SPACE CONSTANTS
C   SPLAPLAC     COMPUTE LAPLACIAN IN SPECTRAL SPACE
C   SPTRANV      PERFORM A VECTOR SPHERICAL TRANSFORM
C
C REMARKS: MINIMUM GRID DIMENSIONS FOR UNALIASED TRANSFORMS TO SPECTRAL:
C   DIMENSION                    LINEAR              QUADRATIC
C   -----------------------      ---------           -------------
C   IMAX                         2*MAXWV+2           3*MAXWV/2*2+2
C   JMAX (IDRT=4,IROMB=0)        1*MAXWV+1           3*MAXWV/2+1
C   JMAX (IDRT=4,IROMB=1)        2*MAXWV+1           5*MAXWV/2+1
C   JMAX (IDRT=0,IROMB=0)        2*MAXWV+3           3*MAXWV/2*2+3
C   JMAX (IDRT=0,IROMB=1)        4*MAXWV+3           5*MAXWV/2*2+3
C   JMAX (IDRT=256,IROMB=0)      2*MAXWV+1           3*MAXWV/2*2+1
C   JMAX (IDRT=256,IROMB=1)      4*MAXWV+1           5*MAXWV/2*2+1
C   -----------------------      ---------           -------------
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      REAL WAVE(*),GRIDMN(KMAX),GRIDXN(*),GRIDXS(*),GRIDYN(*),GRIDYS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SET PARAMETERS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      KW=KWSKIP
      IF(KW.EQ.0) KW=2*MX
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO GRID
      IF(IDIR.GT.0) THEN
C$OMP PARALLEL DO PRIVATE(KWS)
        DO K=1,KMAX
          KWS=(K-1)*KW
          GRIDMN(K)=WAVE(KWS+1)/SQRT(2.)
          CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),1)
          WZ(1:2*MX,K)=0.
        ENDDO
        CALL SPTRANV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &               IPRIME,ISKIP,JNSKIP,JSSKIP,MDIM,KGSKIP,
     &               JBEG,JEND,JCPU,
     &               WD,WZ,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM GRID TO WAVE
      ELSE
C$OMP PARALLEL DO
        DO K=1,KMAX
          WD(1:2*MX,K)=0.
          WZ(1:2*MX,K)=0.
        ENDDO
        CALL SPTRANV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &               IPRIME,ISKIP,JNSKIP,JSSKIP,MDIM,KGSKIP,
     &               JBEG,JEND,JCPU,
     &               WD,WZ,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)
        IF(JBEG.EQ.0) THEN
C$OMP PARALLEL DO PRIVATE(KWS)
          DO K=1,KMAX
            KWS=(K-1)*KW
            CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),-1)
            WAVE(KWS+1)=GRIDMN(K)*SQRT(2.)
          ENDDO
        ELSE
C$OMP PARALLEL DO PRIVATE(KWS)
          DO K=1,KMAX
            KWS=(K-1)*KW
            CALL SPLAPLAC(IROMB,MAXWV,ENN1,WZ(1,K),WD(1,K),-1)
            WAVE(KWS+1:KWS+2*MX)=WAVE(KWS+1:KWS+2*MX)+WZ(1:2*MX,K)
            WAVE(KWS+1)=GRIDMN(K)*SQRT(2.)
          ENDDO
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
