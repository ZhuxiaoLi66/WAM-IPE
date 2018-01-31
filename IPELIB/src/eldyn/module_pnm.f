!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_pnm
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
      public :: pnm

      contains
                                                                      

      subroutine pnm( ct, p )
!----------------------------------------------------------------      
! Purpose: normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
! Method:
!   P_m^m    = sqrt(1+1/2m)*si*P_m-1^m-1                  m>0
!   P_n^m    = [cos*P_n-1^m - R_n-1^m*P_n-2^m ]/R_n^m     n>m>=0
!   dP/d phi = n*cos*P_n^m/sin-(2*n+1)*R_n^m*P_n-1^m/sin  n>=m>=0
!   R_n^m    = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!--------------------------------------------------------------------                                                                   
!nm20121003
      USE efield !,ONLY:
      implicit none

!-----------------------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------------------
      real, intent(inout) :: ct ! cos(colat)                 
      real, intent(inout) :: p(0:nm,0:mm)

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n, np
      real :: pm2, st

!      ct = min( ct,.99 )		! cos(colat)
      st = sqrt( 1. - ct*ct ) 	! sin(colat)

      p(0,0) = 1.  
      do mp = 1,mmp  ! m+1=1,mm+1
        m = mp - 1
	if( m >= 1 ) then
           p(m,m) = pmopmmo(m)*p(m-1,m-1)*st 			
	end if
	pm2 = 0.                                                                  
	do n = mp,nm                    ! n=m+1,N
	   np     = n + 1
	   p(n,m) = (ct*p(n-1,m) - r(n-1,m)*pm2)/r(n,m)
	   pm2    = p(n-1,m)
        end do
      end do

      end subroutine pnm                                                                         



      end module module_pnm
