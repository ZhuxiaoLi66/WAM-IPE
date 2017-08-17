!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_get_efield
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
      public ::        get_efield     ! interface routine


      contains



      subroutine get_efield
!-----------------------------------------------------------------------
! Purpose: calculates the global electric potential field on the
!          geomagnetic grid (MLT in deg) and derives the electric field 
!
! Method:
!
! Author: A. Maute Dec 2003  am 12/17/03    
!-----------------------------------------------------------------------

!     use time_manager,   only : get_curr_calday, get_curr_date
!     use mo_solar_parms, only : get_solar_parms
!     use mag_parms,      only : get_mag_parms
!     use cam_control_mod, only: magfield_fix_year
!     use spmd_utils,      only: masterproc
!nm20121003
      USE efield !,ONLY:
      USE module_DerivPotential ,ONLY: DerivPotential
      USE module_GlobalElPotential ,ONLY:GlobalElPotential
      USE module_adj_S_a ,ONLY:adj_S_a
      USE module_input_parameters,ONLY:sw_debug

      integer :: idum1, idum2, tod ! time of day [s] 
      real kp

!-----------------------------------------------------------------------
! get current calendar day of year & date components 
! valid at end of current timestep
!-----------------------------------------------------------------------
!     iday = get_curr_calday()                   ! day of year
!     call get_curr_date (iyear,imo,iday_m,tod)! year, time of day [sec]
!     iyear = magfield_fix_year
!     iyear = 1995

!     if( iyear < 1900 ) then
!       write(iulog,"(/,'>>> get_efield: year < 1900 not possible: 
!    &year=',i5)") iyear
!       call endrun
!     end if

      tod=ut*3600.
!     ut = tod/3600.                   ! UT of day [sec]

!-----------------------------------------------------------------------
! get solar parms
!-----------------------------------------------------------------------
!     call get_solar_parms( f107_s = f107d )
!-----------------------------------------------------------------------
! get mag parms
!-----------------------------------------------------------------------
!     call get_mag_parms( by = by, bz = bz )
!     by=0.
!     bz=.433726-kp*(0.0849*kp+0.0810363)+
!    &f107d*(0.0079373-0.00219316*kp)
!     print*,kp,bz
!     bz=-2.5
      if ( sw_debug ) 
     & print*,'By=',by,' Bz=',bz,' F107d=',f107d,' UT[hr]',ut
!#ifdef EFIELD_DIAGS
!      if( masterproc ) then
!         write(iulog,*) 'get_efield: f107d,by,bz = ', f107d,by,bz 
!      end if
!#endif
!-----------------------------------------------------------------------
! ajust S_a
!-----------------------------------------------------------------------
      call adj_S_a
!-----------------------------------------------------------------------
! calculate global electric potential    
!-----------------------------------------------------------------------
      call GlobalElPotential

!-----------------------------------------------------------------------
! calculate derivative of global electric potential 
!-----------------------------------------------------------------------
      call DerivPotential

      end subroutine get_efield

      end module module_get_efield
