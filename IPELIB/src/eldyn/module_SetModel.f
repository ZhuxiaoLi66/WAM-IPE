!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_SetModel
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
      public :: SetModel

      contains
!-----------------------------------------------------------------------
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!================================================================================================

	SUBROUTINE SetModel(angle,Bt,Tilt,SWVel)
!
!-----------------------------------------------------------------------
! Calculate the complete set of spherical harmonic coefficients,
! given an arbitrary IMF angle (degrees from northward toward +Y),
! magnitude Bt (nT), dipole tilt angle (degrees),
! and solar wind velocity (km/sec).
! Returns the Coef in the common block SetCoef.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        USE module_FSVal, ONLY: FSVal
        implicit none 
!
!-------------------------------Commons---------------------------------
!
	INTEGER MaxL,MaxM,MaxN
	REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	COMMON /AllCoefs/Cn,MaxL,MaxM,MaxN

	INTEGER ML,MM
	REAL Coef(0:1,0:8,0:3),pi
	COMMON/SetCoef/Coef,pi,ML,MM
!
!------------------------------Arguments--------------------------------
!
	REAL angle,Bt,Tilt,SWVel
!
!---------------------------Local variables-----------------------------
!
        integer n, k, ilimit, i, klimit, l, m, mlimit
	REAL FSC(0:1,0:4),  omega, sintilt !nm20121012
!
!-----------------------------------------------------------------------
!
	pi=2.*ASIN(1.)
	ML=MaxL
	MM=MaxM
	SinTilt=SIN(Tilt*pi/180.)
!	SinTilt=SIND(Tilt)

	omega=angle*pi/180.

        fsc(1,0) = 0.
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
! Retrieve the regression coefficients and evaluate the function
! as a function of Bt,Tilt,and SWVel to get each Fourier coefficient.
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  FSC(i,n)=Cn(0,i,n,k,l,m) + Bt*Cn(1,i,n,k,l,m) +  
     &	   SinTilt*Cn(2,i,n,k,l,m) + SWVel*Cn(3,i,n,k,l,m)
		ENDDO
	      ENDDO
! Next evaluate the Fourier series as a function of angle.
      	      Coef(k,l,m)=FSVal(omega,MaxN,FSC)
	    ENDDO
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE SetModel

!-----------------------------------------------------------------------
      end module module_SetModel
