!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_interp_poten
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
      public :: interp_poten

      contains
!-----------------------------------------------------------------------
      subroutine interp_poten( pot_highlats, pot_highlat, pot_midlat, 
     &ihlat_bnd, itrans_width ) 
!-------------------------------------------------------------------
! Purpose: construct a smooth global electric potential field 
!
! Method: construct one global potential field
! 1. low/mid latitude: |lam| < bnd-trans_width
!   Phi(phi,lam) = Phi_low(phi,lam)
! 2. high latitude: |lam| > bnd+trans_width
!   Phi(phi,lam) = Phi_hl(phi,lam)
! 3. transition zone: bnd-trans_width <= lam <= bnd+trans_width 
! a. interpolate between high and low/midlatitude potential
!   Phi*(phi,lam) = 1/15*[ 5/(2*trans_width) * {Phi_low(phi,bnd-trans_width)*
!   [-lam+bnd+trans_width] + Phi_hl(phi,bnd+trans_width)*
!   [lam-bnd+trans_width]} + 10/(2*trans_width) {Phi_low(phi,lam)*
!   [-lam+bnd+trans_width] + Phi_hl(phi,lam)*
!   [lam-bnd+trans_width]}]
! b.  Interpolate between just calculated Potential and the high latitude
!    potential in a 3 degree zone poleward of the boundary:
!    bnd+trans_width < lam <= bnd+trans_width+ 3 deg 
!   Phi(phi,lam) = 1/3 { [3-(lam-bnd-trans_width)]* Phi*(phi,lam) +
!   [lam-bnd-trans_width)]* Phi_hl*(phi,lam) }
!
! Author: A. Maute Nov 2003  am 11/21/03      
!------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
!------------------------------------------------------------------
!	... dummy arguments
!------------------------------------------------------------------
      integer, intent(in)  :: ihlat_bnd(0:nmlon)
      integer, intent(in)  :: itrans_width(0:nmlon)
      real, intent(in) :: pot_highlats(0:nmlon,0:nmlat)
      real, intent(in) :: pot_highlat(0:nmlon,0:nmlat)
      real, intent(in) :: pot_midlat(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      real, parameter :: fac = 1./3.
      integer  :: ilon, ilat
      integer  :: ibnd, tw, hb1, hb2, lat_ind
      integer  :: j1, j2
      real :: a, b, lat, b1, b2
      real :: wrk1, wrk2
!nm20160401: imf by debug by tomoko
      integer  :: min_ilat

!$omp parallel do private(ilat,ilon,ibnd,tw)
      do ilon = 0,nmlon
        ibnd = ihlat_bnd(ilon)     ! high latitude boundary index
	tw   = itrans_width(ilon)  ! width of transition zone (index)
!-------------------------------------------------------------------
! 1. low/mid latitude: |lam| < bnd-trans_width
!   Phi(phi,lam) = Phi_low(phi,lam)
!-------------------------------------------------------------------
        do ilat = 0,nmlath-(ibnd+tw+1)
          potent(ilon,nmlath+ilat) = pot_midlat(ilon,nmlath+ilat)
          potent(ilon,nmlath-ilat) = pot_midlat(ilon,nmlath+ilat)
        end do
!------------------------------------------------------------------
! 2. high latitude: |lam| > bnd+trans_width
!   Phi(phi,lam) = Phi_hl(phi,lam)
!------------------------------------------------------------------
        do ilat = 0,ibnd-tw-1
          potent(ilon,ilat)       = pot_highlats(ilon,nmlat-ilat)
          potent(ilon,nmlat-ilat) = pot_highlats(ilon,nmlat-ilat)
        end do
      end do
!------------------------------------------------------------------
! 3. transition zone: bnd-trans_width <= lam <= bnd+trans_width 
!------------------------------------------------------------------
! a. interpolate between high and low/midlatitude potential
! update only southern hemisphere (northern hemisphere is copied
! after smoothing)
!------------------------------------------------------------------
!!$omp parallel do private(ilat,ilon,ibnd,tw,a,b,b1,b2,hb1,hb2,lat_ind,
!    &j1,j2,wrk1,wrk2)
      do ilon = 0,nmlon
        ibnd = ihlat_bnd(ilon)          ! high latitude boundary index
	tw   = itrans_width(ilon)       ! width of transition zone (index)
        a    = 1./(2.*tw)
	b1   = (nmlath - ibnd + tw)*a
	b2   = (nmlath - ibnd - tw)*a
	hb1  = nmlath - (ibnd + tw)
	j1   = nmlath - hb1
	hb2  = nmlath - (ibnd - tw)
	j2   = nmlath - hb2
!nm20160401 imf by debug by tomoko
        if (j2 < 0) j2 = 0              ! Tomoko's fix - j2 >= 0
	wrk1 = pot_midlat(ilon,j1)
	wrk2 = pot_highlats(ilon,j2)
!        write(iulog,*) 'pot_all ',ilon,hb1,hb2,nmlath -ibnd,tw
!nm20160401 imf by debug by tomoko
        min_ilat = ibnd-tw
        if (min_ilat < 0) min_ilat = 0  ! Tomoko's fix
        do ilat = min_ilat,ibnd+tw
!	do ilat = ibnd-tw,ibnd+tw
	  lat_ind = nmlath - ilat
          potent(ilon,ilat) =  
     &    fac*((wrk1 + 2.*pot_midlat(ilon,ilat))*(b1 - a*lat_ind)  
     &	  + (wrk2 + 2.*pot_highlats(ilon,ilat))*(a*lat_ind - b2))
          potent(ilon,nmlat-ilat) = potent(ilon,ilat)
        end do
!------------------------------------------------------------------
! b.  Interpolate between just calculated Potential and the high latitude
!    potential in a 3 degree zone poleward of the boundary
!------------------------------------------------------------------
	do ilat = hb2+1,nmlath
	  a = max( 3./dlatm - (ilat - hb2 - 1),0. )
	  b = 3./dlatm - a
          potent(ilon,nmlath-ilat) = (a*potent(ilon,nmlath-ilat)   
     &    + b*pot_highlat(ilon,nmlath-ilat))/3.*dlatm
          potent(ilon,nmlath+ilat) = potent(ilon,nmlath-ilat)
        end do
      end do      

      end subroutine interp_poten
!-----------------------------------------------------------------------
      end module module_interp_poten
