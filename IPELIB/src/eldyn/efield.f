!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module efield
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

      public
!      public :: ed1,           ! zonal electric field Ed1  [V/m] 
!     &          ed2,           ! meridional electric field Ed2 [V/m] 
!     &          potent,        ! electric potential [V]
!     &	        nmlon, nmlat,  ! dimension of mag. grid 
!     &          dlatm, dlonm,  ! grid spacing of mag. grid 
!     &	        ylonm, ylatm   ! magnetic longitudes/latitudes (deg)
!     &,iday,iyear,iday_m,imo,f107d,by,bz,ut
!nm20121003:
!     &,day,dy2yr,date,nmlath,nmlat_wei,deg2mlt
!     &,rtd,sqr2,nm,mm,mmp,pmopmmo,nm,
!     private

      integer :: iday     ! day number of year
      integer :: iyear    ! year
      integer :: iday_m   ! day of month
      integer :: imo      ! month
      real    :: ut       ! universal time  

!---------------------------------------------------------------------- 
! solar parameters
!---------------------------------------------------------------------- 
      real ::   f107d           ! 10.7 cm solar flux
      real ::   by              ! By component of IMF [nT]
      real ::   bz              ! Bz component of IMF [nT]
!nm20121003      private

!---------------------------------------------------------------------- 
! mag. grid dimensions (assumed resolution of 2deg)
!---------------------------------------------------------------------- 
      integer, parameter :: nmlon   = 180     ! mlon 
      integer, parameter :: nmlat   = 90      ! mlat
      integer, parameter :: nmlath  = nmlat/2 ! mlat/2
      integer, parameter :: nmlonh  = nmlon/2 ! mlon/2
      integer, parameter :: nmlonp1 = nmlon+1 ! mlon+1 
      integer, parameter :: nmlatp1 = nmlat+1 ! mlat+1
      integer, parameter :: iulog   = 10

      real :: ylatm(0:nmlat) ! magnetic latitudes  (deg)
      real :: ylonm(0:nmlon) ! magnetic longitudes (deg)
      real :: dlonm          ! delon lon grid spacing
      real :: dlatm	     ! delat lat grid spacing

!---------------------------------------------------------------------- 
! array on magnetic grid:    
!---------------------------------------------------------------------- 
      real :: potent(0:nmlon,0:nmlat) ! electric potential   [V]  
      real :: ed1   (0:nmlon,0:nmlat) ! zonal electric field Ed1  [V/m] 
      real :: ed2   (0:nmlon,0:nmlat) ! meridional electric field Ed2/sin I_m  [V/m]  
       
      real :: date   ! iyear+iday+ut
      real :: day    ! iday+ut

      logical, parameter :: iutav=.false.  ! .true.  means UT-averaging 
                                           ! .false. means no UT-averaging
!     real, parameter ::  v_sw = 400.      ! solar wind velocity [km/s]
      real, parameter ::  v_sw = 450.      ! solar wind velocity [km/s]

!---------------------------------------------------------------------- 
! boundary for Weimer
!---------------------------------------------------------------------- 
      real, parameter :: bnd_wei = 44. ! colat. [deg]
      integer         :: nmlat_wei
      
!---------------------------------------------------------------------- 
! flag for choosing factors for empirical low latitude model      
!---------------------------------------------------------------------- 
      integer, parameter ::  iseasav = 0  ! flag for season 

!---------------------------------------------------------------------- 
! constants:
!---------------------------------------------------------------------- 
      real,parameter :: r_e  =  6.371e6 ! radius_earth [m] (same as for apex.F90)
      real,parameter :: h_r  = 130.0e3  ! reference height [m] (same as for apex.F90)
      real,parameter :: dy2yr= 365.24   ! day per avg. year used in Weimer
      real,parameter :: dy2mo= 30.6001  ! day per avg. month used in Weimer
      real,parameter :: pi=3.141592653

      real :: rtd                ! radians -> deg
      real :: dtr                ! deg -> radians
      real :: sqr2           
      real :: hr2rd              ! pi/12 hrs
      real :: dy2rd              ! 2*pi/365.24  average year
      real :: deg2mlt            ! for mlon to deg
      real :: mlt2deg            ! for mlt to mlon
      real :: sinIm_mag(0:nmlat) ! sinIm

      integer :: jmin, jmax   ! latitude index for interpolation of 
                              ! northward e-field ed2 at mag. equator

!---------------------------------------------------------------------- 
!  for spherical harmonics
!---------------------------------------------------------------------- 
      integer, parameter ::  nm = 19, mm = 18, mmp = mm + 1	  

      real :: r(0:nm,0:mm)      ! R_n^m
      real :: pmopmmo(0:mm)     ! sqrt(1+1/2m)

!---------------------------------------------------------------------- 
!  index for factors f_m(mlt),f_l(UT),f_-k(d)
!---------------------------------------------------------------------- 
      integer, parameter :: ni = 1091  ! for n=12 m=-18:18
      integer :: imax                                         ! max number of index
      integer,dimension(0:ni) :: kf,lf, mf, nf, jf
      real :: ft(1:3,0:2)  ! used for f_-k(season,k)

      real ::  a_klnm(0:ni)        !  A_klm
      real ::  a_lf  (0:ni)        ! A_klmn^lf for minimum  
      real ::  a_hf  (0:ni)        ! A_klmn^hf for maximum

!---------------------------------------------------------------------- 
! high_latitude boundary
!---------------------------------------------------------------------- 
      real,parameter ::ef_max  = 0.015   ! max e-field for high latitude boundary location [V/m]
      real,parameter ::lat_sft = 54.     ! shift of highlat_bnd to 54 deg
      integer            :: ilat_sft     ! index of shift for high latitude boundary
      integer, parameter :: nmax_sin = 2 ! max. wave number to be represented
      logical, parameter :: debug =.false.

!nm20121003:subroutines are separated into sub_efield.f.

      end module efield
