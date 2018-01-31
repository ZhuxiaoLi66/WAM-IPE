!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_GlobalElPotential
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
      public :: GlobalElPotential

      contains

      subroutine GlobalElPotential
!-----------------------------------------------------------------------
! Purpose: calculates the global electric potential field on the
!          geomagnetic grid (MLT in deg) 
!
! Method: rewritten code from Luedger Scherliess (11/20/99 LS)
!     routine to calculate the global electric potential in magnetic
!     Apex coordinates (Latitude and MLT).
!     High Latitude Model is Weimer 1996.
!     Midlatitude model is Scherliess 1999.
!     Interpolation in a transition region at about 60 degree 
!     magnetic apex lat
!
! Author: A. Maute Dec 2003  am 12/17/03 
!-----------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      USE module_interp_poten ,ONLY: interp_poten
      USE module_efield_mid ,ONLY:efield_mid
      USE module_prep_weimer ,ONLY:prep_weimer
      USE module_pot_latsmo ,ONLY: pot_latsmo
      USE module_pot_latsmo2 ,ONLY: pot_latsmo2
      USE module_pot_lonsmo ,ONLY: pot_lonsmo
      USE module_highlat_getbnd ,ONLY: highlat_getbnd
      USE module_bnd_sinus ,ONLY:      bnd_sinus
      USE module_highlat_adjust ,ONLY: highlat_adjust
      USE module_EpotVal ,ONLY: EpotVal
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: ilon, ilat, idlat
      integer  :: ihlat_bnd(0:nmlon)     ! high latitude boundary
      integer  :: itrans_width(0:nmlon)  ! width of transition zone
      real :: mlt, mlon, mlat, mlat_90, pot
      real :: pot_midlat(0:nmlon,0:nmlat) ! potential from L. Scherliess model
      real :: pot_highlat(0:nmlon,0:nmlat) ! potential from Weimer model
      real :: pot_highlats(0:nmlon,0:nmlat)! smoothed potential from Weimer model

!-----------------------------------------------------------------------
! Externals
!-----------------------------------------------------------------------
!nm20121012      real,external :: EpotVal        ! in wei96.f

!-----------------------------------------------------------------------
! convert to date and day	
!-----------------------------------------------------------------------
      day  = iday + ut/24.
      date = iyear + day/dy2yr

!-----------------------------------------------------------------------
! low/midlatitude electric potential - empirical model Scherliess 1999  
!-----------------------------------------------------------------------
!$omp parallel do private(ilat, ilon, mlat, pot)
      do ilat = 0,nmlath                        ! Calculate only for one magn. hemisphere
	mlat = ylatm(ilat)                      ! mag. latitude
        do ilon = 0,nmlon	                ! lon. loop
          call efield_mid( mlat, ylonm(ilon), pot )
	  pot_midlat(ilon,ilat+nmlath) = pot	! SH/NH symmetry 
	  pot_midlat(ilon,nmlath-ilat) = pot
        end do
      end do

!-----------------------------------------------------------------------
! hight latitude potential from Weimer model
! at the poles Weimer potential is not longitudinal dependent
!-----------------------------------------------------------------------
      call prep_weimer    ! calculate IMF angle & magnitude, tilt

!$omp parallel do private(ilat, ilon, mlat_90, pot)
      do ilat = 0,nmlat_wei  ! Calculate only for one magn. hemisphere
        mlat_90 = 90. - ylatm(ilat)  ! mag. latitude
        do ilon = 0,nmlon
    	  pot  = 1000.*EpotVal( mlat_90, ylonm(ilon)*deg2mlt ) ! calculate potential (kv -> v)
!-----------------------------------------------------------------------
! NH/SH symmetry
!-----------------------------------------------------------------------
    	  pot_highlat(ilon,ilat)        = pot
    	  pot_highlat(ilon,nmlat-ilat)  = pot
    	  pot_highlats(ilon,ilat)       = pot
    	  pot_highlats(ilon,nmlat-ilat) = pot
        end do
      end do     

!-----------------------------------------------------------------------
! weighted smoothing of high latitude potential
!-----------------------------------------------------------------------
      idlat = 2              ! smooth over -2:2 = 5 grid points
      call pot_latsmo( pot_highlats, idlat )
!-----------------------------------------------------------------------
! calculate the height latitude bounday ihl_bnd
! 1. calculate E field from weimar model
!    boundary is set where the total electric field exceeds
!    0.015 V/m (corresp. approx. to 300 m/s)
! 2. moved halfways to 54 deg 
! output : index 0-pole nmlath-equator
!-----------------------------------------------------------------------
      call highlat_getbnd( ihlat_bnd )
!-----------------------------------------------------------------------
! 3. adjust high latitude boundary sinusoidally
!    calculate width of transition zone
!-----------------------------------------------------------------------
      call bnd_sinus( ihlat_bnd, itrans_width ) 
!-----------------------------------------------------------------------
! 4. ajust high latitude potential to low latitude potential      
!-----------------------------------------------------------------------
      call highlat_adjust( pot_highlats, pot_highlat, pot_midlat, 
     &ihlat_bnd )
!-----------------------------------------------------------------------
! interpolation of high and low/midlatitude potential in the
! transition zone and put it into global potent array
!-----------------------------------------------------------------------
      call interp_poten( pot_highlats, pot_highlat, pot_midlat, 
     &ihlat_bnd, itrans_width) 
!-----------------------------------------------------------------------
! potential weighted smoothing in latitude
!-----------------------------------------------------------------------
      idlat = 2                 ! smooth over -2:2 = 5 grid points
      call pot_latsmo2( potent, idlat )
!-----------------------------------------------------------------------
! potential smoothing in longitude
!-----------------------------------------------------------------------
      idlat = nmlon/48          ! smooth over -idlat:idlat grid points
      call pot_lonsmo( potent, idlat )

!-----------------------------------------------------------------------
! output
!-----------------------------------------------------------------------
! output ( change later to netcdf file)
!      do ilat=0,nmlat
!       do ilon=0,nmlon
!         write(iulog,'(4(x,f12.5))') ylatm(ilat),ylonm(ilon), &
!           potent(ilon,ilat),potent(ilon,nmlat-ilat)
!         write(iulog,'(4(x,f12.5))') ylatm(ilat),ylonm(ilon), &
!           potent(ilon,ilat),potent(ilon,nmlat-ilat)
!	write(iulog,'(f10.3)') potent(ilon,ilat)
!       end do
!      end do

      end subroutine GlobalElPotential
      end module module_GlobalElPotential
