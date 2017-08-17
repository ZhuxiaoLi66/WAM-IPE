!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_prep_pnm
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: prep_pnm

      contains

                                                                      



      subroutine prep_pnm
!-----------------------------------------------------------------      
! Purpose: constant factors for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:
!   PmoPmmo(m) = sqrt(1+1/2m)
!   R_n^m      = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!-----------------------------------------------------------------     
!nm20121003
      USE efield !,ONLY:
      implicit none                

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n
      real :: xms, xns, den

      do mp = 1, mmp            ! m+1 = 1,mm+1                                     
	m = mp - 1                                               
	xms = m*m                                                
	if( mp /= 1 ) then
           pmopmmo(m) = sqrt( 1. + .5/M )
	end if
	do n = m,nm      ! n = m,N                                     
	  xns    = n*n                                       
	  den    = max(4.*xns - 1.,1.)
	  r(n,m) = sqrt( (xns  - xms)/den )
	end do                 
      end do 

      end subroutine prep_pnm                                                                         



      end module module_prep_pnm

