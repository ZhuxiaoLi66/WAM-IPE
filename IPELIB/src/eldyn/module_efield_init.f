!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_efield_init
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

!     use shr_kind_mod,  only: r8 => shr_kind_r8
!     use physconst,     only: pi
!     use abortutils,    only: endrun
!     use cam_logfile,   only: iulog
   
      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: efield_init   ! interface routine


      contains

      subroutine efield_init(efield_lflux_file, efield_hflux_file,      &
     &                       efield_wei96_file)
      USE efield !,ONLY:
      USE module_prep_pnm,ONLY:prep_pnm
      USE module_index_quiet ,ONLY:index_quiet
      USE module_read_acoef ,ONLY: read_acoef
      USE module_constants ,ONLY:constants
      USE module_prep_fk ,ONLY:prep_fk
      USE module_ReadCoef ,ONLY:ReadCoef
      implicit none
!--------------------------------------------------------------------
! Purpose: read in and set up coefficients needed for electric field
!          calculation (independent of time & geog. location)
!
! Method:
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-------------------------------------------------------------------
      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file
      character(len=*), intent(in) :: efield_wei96_file

      call constants	 ! calculate constants
!-----------------------------------------------------------------------
! low/midlatitude potential from Scherliess model
!-----------------------------------------------------------------------
      call read_acoef (efield_lflux_file, efield_hflux_file)	! read in A_klnm for given S_aM
      call index_quiet  ! set up index for f_m(mlt),f_l(UT),f_-k(d)
      call prep_fk	! set up the constant factors for f_k
      call prep_pnm	! set up the constant factors for P_n^m & dP/d phi
!-----------------------------------------------------------------------
!following part should be independent of time & location if IMF constant
!-----------------------------------------------------------------------
      call ReadCoef (efield_wei96_file)

      end subroutine efield_init


      end module module_efield_init
