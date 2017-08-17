!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_constants
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
      public :: constants

      contains
!-----------------------------------------------------------------------
      subroutine constants
!---------------------------------------------------------------
! Purpose: set up constant values (e.g. magnetic grid, convertion
!      constants etc)
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!--------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
!-------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------
      integer  :: i,j
      real :: fac,lat

      rtd     = 180./pi 	        ! radians -> deg
      dtr     = pi/180.	        ! deg -> radians
      sqr2    = sqrt(2.e0)
      hr2rd   = pi/12.	        ! pi/12 hrs
      dy2rd   = 2.*pi/dy2yr          ! 2*pi/365.24  average year
      deg2mlt = 24./360.          ! convert degrees to MLT hours
      mlt2deg = 360./24.          ! for mlt to mlon       

!-------------------------------------------------------------------
! Set grid deltas:
!-------------------------------------------------------------------
      dlatm = 180./nmlat
      dlonm = 360./nmlon

!-------------------------------------------------------------------
! Set magnetic latitude array 
!-------------------------------------------------------------------
      do j = 0,nmlat
        ylatm(j) = j*dlatm
        lat = (ylatm(j) - 90.)*dtr
	fac = cos(lat)    ! sinIm = 2*sin(lam_m)/sqrt[4-3*cos^2(lam_m)]
	fac = 4. - 3.*fac*fac
	fac = 2./sqrt( fac )
	sinIm_mag(j) = fac*sin( lat )
      end do 

!------------------------------------------------------------------
! Set magnetic longitude array
!------------------------------------------------------------------
      do i = 0,nmlon
        ylonm(i) = i*dlonm
      end do ! i=1,nmlonp1

!-----------------------------------------------------------------
! find boundary index for weimer
!------------------------------------------------------------------
      do j = 0,nmlat
        nmlat_wei = j
        if( bnd_wei <= ylatm(j) ) then
           exit
        end if
      end do 

!-------------------------------------------------------------------
! find latitudinal shift
!-------------------------------------------------------------------
      do j = 0,nmlat
        ilat_sft = j
        if( lat_sft <= ylatm(j) ) then
           exit
        end if
      end do 

!------------------------------------------------------------------
! find index for linear interpolation of ed2 at mag.equator 
! use 12 deg - same as in TIEGCM      
!------------------------------------------------------------------
      do j = 0,nmlat
        lat = ylatm(j) - 90.
        if( lat <= -12. ) then
	  jmin = j
        else if( lat > 12. ) then
	  jmax = j
	  exit
       end if
      end do

      end subroutine constants
!-----------------------------------------------------------------------
      end module module_constants
