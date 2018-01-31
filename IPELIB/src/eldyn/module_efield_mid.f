!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_efield_mid
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
      public ::efield_mid 

      contains
!-----------------------------------------------------------------------
      subroutine efield_mid( mlat, mlon, pot )
!------------------------------------------------------------------
! Purpose: calculate the electric potential for low and 
!      midlatitudes from an empirical model (Scherliess 1999)
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      USE module_ff ,ONLY:ff
      USE module_pnm ,ONLY:pnm
      USE module_set_fkflfs ,ONLY:set_fkflfs
!------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      real, intent(in)  :: mlat, mlon
      real, intent(out) :: pot               ! electric potential (V)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: i, mp, np, nn
      real :: mod_mlat, ct, x
      real :: fk(0:2)      	    ! f_-k(day) 
      real :: fl(-2:2)          ! f_l(ut)  
      real :: fs(2)	            ! f_s(f10.7) 
      real :: f(-18:18)
      real :: p(0:nm,0:mm)      ! P_n^m	 spherical harmonics

      pot = 0. ! initialize                                        

      mod_mlat = mlat
      if( abs(mlat) <= 0.5 ) then
         mod_mlat = 0.5                     ! avoid geomag.equator
      end if

!------------------------------------------------------------------
! set f_-k, f_l, f_s depending on seasonal flag
!------------------------------------------------------------------
      call set_fkflfs( fk, fl, fs ) 
      
!------------------------------------------------------------------
! spherical harmonics 
!------------------------------------------------------------------
      ct = cos( (90. - mod_mlat)*dtr )  ! magnetic colatitude 
      call pnm( ct, p )	                   ! calculate P_n^m
      call ff( mlon, 18, f )               ! calculate f_m (phi) why 18 if N=12                              

      do i = 0,imax  
        mp  = mf(i)                                                      
        np  = nf(i)
        nn  = abs(mp)                      !   P_n^m = P_n^-m  
        x   = a_klnm(i)* fl(lf(i)) * fk(kf(i)) * fs(jf(i))
	pot = pot + x*f(mp)*p(np,nn) 
      end do 
      
      end subroutine efield_mid   
!-----------------------------------------------------------------------
      end module module_efield_mid
