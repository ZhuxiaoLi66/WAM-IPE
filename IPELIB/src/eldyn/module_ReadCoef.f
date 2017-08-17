!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_ReadCoef
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
      public :: ReadCoef

      contains
!*********************** Copyright 1996, Dan Weimer/MRC ***********************

	SUBROUTINE ReadCoef (wei96_file)
!
!-----------------------------------------------------------------------
!
! Read in the data file with the model coefficients
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!
! NCAR addition (Jan 97):  initialize constants used in GECMP
!-----------------------------------------------------------------------
!
c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use ioFileMod,     only : getfil
c     use units,         only : getunit, freeunit
c     use abortutils,    only : endrun
c     use cam_logfile,   only : iulog
      implicit none 
!
!-------------------------------Commons---------------------------------
!
      real alamn, alamx, alamr, stpd, stp2, cstp, sstp
      COMMON /CECMP/ ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP
!            ALAMN = Absolute min latitude (deg) of model
!            ALAMX = Absolute max latitude (deg) for normal gradient calc.
!            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
!            STP2  = Denominator in gradient calc

!
!------------------------------Arguments--------------------------------
!
      character(len=*), intent(in) :: wei96_file
!
!-----------------------------Parameters------------------------------
!
      real d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,  
     &            R2D = 57.2957795130823208767981548147)
!
!---------------------------Local variables-----------------------------
!
      INTEGER udat,unit,ios
      integer ll,mm,k,m,klimit,kk,nn,ii,i,n,ilimit,mlimit,l

      REAL C(0:3)
      real stpr, step

      CHARACTER*15 skip

      INTEGER MaxL,MaxM,MaxN,iulog
      REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
      COMMON /AllCoefs/Cn,MaxL,MaxM,MaxN

      character(len=256) :: locfn
!
!-----------------------------------------------------------------------
      iulog=14
      STEP = 10.
      STPR = STEP/6671.
      STPD = STPR*R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMN = 45.
      ALAMX = 90. - STPD
      ALAMR = ALAMN*D2R
!          End NCAR addition
! 
!  get coeff_file  
c     unit= getunit()
      unit= 200
c     write(iulog,*) 'Weimer: getting file ',trim(wei96_file),
c    &' unit ',unit
c     call getfil( wei96_file, locfn, 0 )
      locfn= wei96_file
!      
c     write(iulog,*) 'Weimer: opening file ',trim(locfn),
c    &' unit ',unit	
      OPEN(unit=unit,file=trim(locfn),  
     &     status = 'old',iostat = ios)
c     if(ios.gt.0) then
c      write(iulog,*) 'Weimer: error in opening wei96.cofcnts',
c    &' unit ',unit
c       call endrun
c     endif
      
  900 FORMAT(A15)
 1000 FORMAT(3I8)
 2000 FORMAT(3I2)
 3000 FORMAT(2I2,4E15.6)

!     READ(udat,900) skip
c     write(iulog,*) 'Weimer: reading file ',trim(locfn),
c    &' unit ',unit	
      READ(unit,1000,iostat = ios) MaxL,MaxM,MaxN
c     if(ios.gt.0) then
c     write(iulog,*) 
c    &'ReadCoef: error in reading wei96.cofcnts file',
c    &' unit ',unit	
c       call endrun
c     endif
      DO l=0,MaxL
        IF(l.LT.MaxM)THEN
          mlimit=l
        ELSE
          mlimit=MaxM
        ENDIF
        DO m=0,mlimit
          IF(m.LT.1)THEN
            klimit=0
          ELSE
            klimit=1
          ENDIF
          DO k=0,klimit
            READ(unit,2000,iostat = ios) ll,mm,kk
c           if(ios.gt.0) then
c     	      write(iulog,*) 
c    &'ReadCoef: error in reading wei96.cofcnts file',' unit ',
c    &unit	
c             call endrun
c           endif
c           IF(ll.NE.l .OR. mm.NE.m .OR. kk.NE.k)THEN
c             WRITE(IULOG,*)'Data File Format Error'
c             CALL ENDRUN
c           ENDIF
            DO n=0,MaxN
              IF(n.LT.1)THEN
        	ilimit=0
              ELSE
        	ilimit=1
              ENDIF
              DO i=0,ilimit
        	READ(unit,3000,iostat = ios) nn,ii,C
c               if(ios.gt.0) then
c     	          write(iulog,*) 'ReadCoef: error in reading',  
c    &                 ' wei96.cofcnts file',' unit ',unit	
c                 call endrun
c               endif
c       	IF(nn.NE.n .OR. ii.NE.i)THEN
c       	  WRITE(IULOG,*)'Data File Format Error'
c       	  CALL ENDRUN
c       	ENDIF
        	Cn(0,i,n,k,l,m)=C(0)
        	Cn(1,i,n,k,l,m)=C(1)
        	Cn(2,i,n,k,l,m)=C(2)
        	Cn(3,i,n,k,l,m)=C(3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!      
      close(unit)
c     call freeunit(unit)
!    
      RETURN
      END SUBROUTINE ReadCoef
!-----------------------------------------------------------------------
      end module module_ReadCoef
