!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_adj_S_a
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
      public :: adj_S_a

      contains
!-----------------------------------------------------------------------
      subroutine adj_S_a
!------------------------------------------------------------------
! adjust S_a -> S_aM   eqn.8-11 Scherliess draft
!------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      implicit none

!-----------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: i
      real :: x2, y2, a90, a180, S_aM

      x2 = 90.*90.
      y2 = (90. - 65.)
      y2 = y2*y2
      a90  = atan2(y2,x2)
      y2 = (180. - 65.)
      y2 = y2*y2
      a180 = atan2(y2,x2)
!     y2 = (S_a-65.)
      y2 = (f107d - 65.)
      y2 = y2*y2
      S_aM = (atan2(y2,x2) - a90)/(a180 - a90) 
      S_aM = 90.*(1. + S_aM)
c     if(debug) write(iulog,*) 'f107d=',f107d,' S_aM =',S_aM
c     if(debug) write(iulog,*) 'By=',by

!-----------------------------------------------------------------
! inter/extrapolate to S_a (f107d)
!----------------------------------------------------------------
      do i = 0,ni                       ! eqn.8 Scherliess draft
        a_klnm(i) = S_aM*(a_hf(i)-a_lf(i))/90.+
     &2.*a_lf(i)- a_hf(i)
! for testing like in original code
!        a_klnm(i)=S_a*(a_hf(i)-a_lf(i))/90.+2.*a_lf(i)-a_hf(i)
!        a_klnm(i)=f107d*(a_hf(i)-a_lf(i))/90.+2.*a_lf(i)-a_hf(i)
      end do

      end subroutine adj_S_a
!-----------------------------------------------------------------------
      end module module_adj_S_a
