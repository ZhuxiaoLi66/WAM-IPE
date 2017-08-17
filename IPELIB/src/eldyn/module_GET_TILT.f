!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_GET_TILT
!--------------------------------------------------------------------- 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
! - low/midlatitudes electric potential is from an empirical model from
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! - the transition zone is smoothed
! - output is the horizontal global electric field in magnetic coordinates direction
!  at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real:: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use physconst,     only: pi
c     use abortutils,    only: endrun
c     use cam_logfile,   only: iulog
   
      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: GET_TILT

      contains
!-----------------------------------------------------------------------
!*********************** Copyright 1996, Dan Weimer/MRC ***********************

!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following routines (translib.for) were added to return the dipole tilt. C
!  GET_TILT was initially a procedure (TRANS), here it has been changed into   C
!  a function which returns the dipole tilt.                                   C
! Barbara Emery (emery@ncar.ucar.edu) and William Golesorkhi, HAO/NCAR (3/96)  C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! COORDINATE TRANSFORMATION UTILITIES
!**********************************************************************        
	FUNCTION GET_TILT(YEAR,MONTH,DAY,HOUR)
!
!-----------------------------------------------------------------------
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following line initially was:                                           C
!       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)                            C
!  It has been changed to return the dipole tilt from this function call.      C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         
!      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
!      TRANSFORMATIONS, IDENTIFIED BY K.
!          K=1 TRANSFORMS GSE to GEO
!          K=2     "      GEO to MAG
!          K=3     "      GSE to MAG
!          K=4     "      GSE to GSM
!          K=5     "      GEO to GSM
!          K=6     "      GSM to MAG
!          K=7     "      GSE to GEI
!          K=8     "      GEI to GEO
!          K=9     "      GSM to SM 
!	   K=10    "      GEO to SM 
!	   K=11    "      MAG to SM 
!
!      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
!      FILE UNIT=IDBUG
!
!      The formal names of the coordinate systems are:
!	GSE - Geocentric Solar Ecliptic
!	GEO - Geographic
!	MAG - Geomagnetic
!	GSM - Geocentric Solar Magnetospheric
!	SM  - Solar Magnetic
!	
!      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
!      ST(I) AND CT(I) ARE SINES & COSINES.       
!
!      Program author:  D. R. Weimer
!
!      Some of this code has been copied from subroutines which had been
!      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
!      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
!
!      The formulas for the calculation of Greenwich mean sidereal time (GMST)
!      and the sun's location are from "Almanac for Computers 1990",
!      U.S. Naval Observatory.
!
!-----------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        USE module_ADJUST, ONLY: ADJUST
        USE module_JULDAY, ONLY: JULDAY
        implicit none 
!
!-----------------------------Return Value--------------------------
!
        real get_tilt
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)

        real epoch, th0, ph0, dipole
        COMMON/MFIELD/EPOCH,TH0,PH0,DIPOLE
c       DATA EPOCH,TH0,PH0,DIPOLE/1980.,11.19,-70.76,.30574/
!
!------------------------------Arguments--------------------------------
!
	INTEGER YEAR, MONTH, DAY
	REAL HOUR
!
!-----------------------------Parameters------------------------------
!
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
!
!---------------------------Local variables-----------------------------
!
        integer IDBUG
        integer j, k, jd, iyr, i, mjd

        REAL UT, T0, GMSTD, GMSTH, ECLIP, MA, LAMD, SUNLON, pi
        real b32, b33, b3
!
!-------------------------External Functions----------------------------
!
!nm20121012        integer julday
!nm20121012        external julday
!
!-----------------------------------------------------------------------
!
        EPOCH=1980.
        TH0=11.19
        PH0=-70.76
        DIPOLE=.30574
	pi=2.*ASIN(1.)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! IDBUG=0 to prevent printing data to the screen or writing data to a file.    C
        IDBUG = 0
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
c       T0=(real(MJD,r8)-51544.5)/36525.0
	T0=(float(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 +36000.770*T0 +3.87933E-4*T0*T0 + 
     &       15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
        MA=357.528 + 35999.050*T0 + 0.041066678*UT
        CALL ADJUST(MA)
        LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
        CALL ADJUST(LAMD)
        SUNLON=LAMD + (1.915-0.0048*T0)*SIN(MA*pi/180.) + 0.020* 
     &     SIN(2.*MA*pi/180.)
        CALL ADJUST(SUNLON)
c         IF(IDBUG.NE.0)THEN
c         WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
c         WRITE(IDBUG,*) 'MJD=',MJD
c         WRITE(IDBUG,*) 'T0=',T0
c         WRITE(IDBUG,*) 'GMSTH=',GMSTH
c         WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
c         WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
c         WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
c         WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
c         ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
! Derived later:
!       CX(6) = Dipole tilt angle  
!       CX(7) = Angle between sun and magnetic pole
!       CX(8) = Subsolar point latitude
!       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIN(CX(I)*pi/180.)
	  CT(I) = COS(CX(I)*pi/180.)
	ENDDO
!         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
!         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) + 
     &AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
!         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
!         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) + 
     &AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
!         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3  
!         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) + 
     &AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + 
     &AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) + 
     &AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + 
     &AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO 
!
	ST(6) = AM(3,1,GSEMAG)        
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIN(ST(6)*pi/180.)  

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0.
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0.
        AM(2,2,GSMSM) = 1.
        AM(2,3,GSMSM) = 0.
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0.
        AM(3,3,GSMSM) = CT(6)  
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +  
     &AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  
     &AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +  
     &  AM(I,2,GSMSM)*AM(J,2,GSMMAG) + 
     &AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
      
!
      CX(7)=ATAN2( AM(2,1,11) , AM(1,1,11) )
      
      CX(7)=CX(7)*180./pi
      CX(8)=ASIN( AM(3,1,1)*pi/180. )
      CX(9)=ATAN2( AM(2,1,1) , AM(1,1,1) )
      CX(9)=CX(9)*180./pi

      IF(IDBUG.NE.0)THEN
c     WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
c     WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',
c    &CX(7)
c     WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
c     WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
c        WRITE(IDBUG,1001) K
         DO I=1,3
c          WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The next line was added to return the dipole tilt from this function call.  C
      GET_TILT = CX(6)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END FUNCTION GET_TILT
!-----------------------------------------------------------------------
      end module module_GET_TILT
