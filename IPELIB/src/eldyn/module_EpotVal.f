!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_EpotVal
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
      public :: EpotVal

      contains
!-----------------------------------------------------------------------
	FUNCTION EpotVal(gLAT,gMLT)
!
!-----------------------------------------------------------------------
! Return the value of the electric potential in kV at
! corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
!
! Must first call ReadCoef and SetModel to set up the model coeficients for
! the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        USE module_LEGENDRE, ONLY: LEGENDRE
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real EpotVal
!
!-------------------------------Commons---------------------------------
!
	INTEGER ML,MM
	REAL Coef(0:1,0:8,0:3),pi
	COMMON/SetCoef/Coef,pi,ML,MM
!
!------------------------------Arguments--------------------------------
!
	REAL gLAT,gMLT
!
!---------------------------Local variables-----------------------------
!
        integer limit,l,m

	Real Theta,Phi,Z,ct,Phim
        real r
	REAL Plm(0:20,0:20)
!
!-----------------------------------------------------------------------
!
	r=90.-gLAT
	IF(r .LT. 45.)THEN
	  Theta=r*pi/45.
          Phi=gMLT*pi/12.
	  Z=Coef(0,0,0)
	  ct=COS(Theta)
	  CALL Legendre(ct,ML,MM,Plm)
	  DO l=1,ML
	    Z=Z + Coef(0,l,0)*Plm(l,0)
	    IF(l.LT.MM)THEN
	      limit=l
	    ELSE
	      limit=MM
	    ENDIF
	    DO m=1,limit
	      phim=phi*m
              Z=Z + Coef(0,l,m)*Plm(l,m)*COS(phim) +  
     &	   Coef(1,l,m)*Plm(l,m)*SIN(phim)
	    ENDDO
	  ENDDO
	ELSE
	  Z=0.
	ENDIF
	EpotVal=Z
	RETURN
	END FUNCTION EpotVal

!-----------------------------------------------------------------------
      end module module_EpotVal
