!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_highlat_getbnd
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
      public :: highlat_getbnd

      contains
!-----------------------------------------------------------------------
      subroutine highlat_getbnd( ihlat_bnd ) 
!------------------------------------------------------------------
! Purpose: calculate the height latitude bounday index ihl_bnd
!
! Method:
! 1. calculate E field from weimar model
!    boundary is set where the total electric field exceeds
!    0.015 V/m (corresp. approx. to 300 m/s)
! 2. moved halfways to 54 deg not necessarily equatorwards as in the
!    original comment from L. Scherliess- or?
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      USE module_GECMP,ONLY:gecmp
!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(out) :: ihlat_bnd(0:nmlon)

!------------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: ilon, ilat, ilat_sft_rvs
      real :: mlat, mlt, es, ez, e_tot

      ilat_sft_rvs = nmlath - ilat_sft          ! pole =0, equ=90
!$omp parallel do private(ilat,ilon,mlt,mlat,es,ez,e_tot)
      do ilon = 0,nmlon                         ! long.
	ihlat_bnd(ilon) = 0
        mlt  = ylonm(ilon)*deg2mlt              ! mag.local time ?
        do ilat = nmlat_wei+1,0,-1              ! lat. loop moving torwards pole
	  mlat = 90. - ylatm(ilat)           ! mag. latitude pole = 90 equator = 0
          call gecmp( mlat, mlt, es, ez )	! get electric field
          e_tot = sqrt( es**2 + ez**2 )
          if( abs(e_tot) >= ef_max ) then                        ! e-filed > limit -> boundary
            ihlat_bnd(ilon) = ilat - (ilat - ilat_sft_rvs)/2     ! shift boundary to lat_sft (54deg)
            exit
          end if
        end do
      end do     

!     write(iulog,"('highlat_getbnd: ihlat_bnd=',/,(12i6))") ihlat_bnd

      end subroutine highlat_getbnd
!-----------------------------------------------------------------------
      end module module_highlat_getbnd
