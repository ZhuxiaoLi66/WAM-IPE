!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_highlat_adjust
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
      public :: highlat_adjust

      contains
!-----------------------------------------------------------------------
      subroutine highlat_adjust( pot_highlats, pot_highlat, 
     &pot_midlat, ihlat_bnd )
!------------------------------------------------------------------
! Purpose: Adjust mid/low latitude electric potential and high latitude
!          potential such that there are continous across the mid to high 
!          latitude boundary
!
! Method:
! 1. integrate Phi_low/mid(phi,bnd) along the boundary mid to high latitude
! 2. integrate Phi_high(phi,bnd) along the boundary mid to high latitude
! 3. adjust Phi_high by delta =
!    Int_phi Phi_high(phi,bnd) d phi/360. - Int_phi Phi_low/mid(phi,bnd) d phi/360.
!
! Author: A. Maute Nov 2003  am 11/21/03
!------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
!------------------------------------------------------------------
!	... dummy arguments
!------------------------------------------------------------------
      integer, intent(in)     :: ihlat_bnd(0:nmlon) ! boundary mid to high latitude
      real, intent(in)    :: pot_midlat(0:nmlon,0:nmlat) ! low/mid latitude potentail
      real, intent(inout) :: pot_highlat(0:nmlon,0:nmlat)! high_lat potential
      real, intent(inout) :: pot_highlats(0:nmlon,0:nmlat)! high_lat potential! smoothed high_lat potential

!------------------------------------------------------------------
! local:     
!------------------------------------------------------------------
      integer  :: bnd, ilon, ilat, ilatS, ibnd60, ibnd_hl
      real :: pot60, pot_hl, del

!-------------------------------------------------------------------
! 1. integrate Phi_low/mid(phi,bnd) along the boundary mid to high latitude
! 2. integrate Phi_high(phi,bnd) along the boundary mid to high latitude
!-------------------------------------------------------------------
      pot60  = 0.
      pot_hl = 0.
      do ilon = 1,nmlon  ! long.           ! bnd -> eq to pole 0:90
    	ibnd60  = nmlat - ihlat_bnd(ilon)   ! 0:180 pole to pole
    	ibnd_hl = ihlat_bnd(ilon)         ! colatitude
        pot60   = pot60 + pot_midlat(ilon,ibnd60)
        pot_hl  = pot_hl + pot_highlats(ilon,ibnd_hl)
      end do
      pot60  = pot60/(nmlon)
      pot_hl = pot_hl/(nmlon)
      
c     if (debug) write(iulog,*) 'Mid-Latitude Boundary Potential =',
c    &pot60
c     if (debug) write(iulog,*) 'High-Latitude Boundary Potential=',
c    &pot_hl

!-------------------------------------------------------------------
! 3. adjust Phi_high by delta =
!    Int_phi Phi_high(phi,bnd) d phi/360. - Int_phi Phi_low/mid(phi,bnd) d phi/360.
!-------------------------------------------------------------------
      del = pot_hl - pot60

!$omp parallel do private(ilat,ilon,ilats)
      do ilat = 0,nmlat_wei      ! colatitude
        ilats = nmlat - ilat
        do ilon = 0,nmlon
	  pot_highlat(ilon,ilat)   = pot_highlat(ilon,ilat)   - del
	  pot_highlat(ilon,ilats)  = pot_highlat(ilon,ilats)  - del
	  pot_highlats(ilon,ilat)  = pot_highlats(ilon,ilat)  - del
	  pot_highlats(ilon,ilats) = pot_highlats(ilon,ilats) - del
        end do
      end do

      end subroutine highlat_adjust
!-----------------------------------------------------------------------
      end module module_highlat_adjust
