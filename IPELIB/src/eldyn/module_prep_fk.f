!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_prep_fk
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: prep_fk

      contains
!-----------------------------------------------------------------------
      subroutine prep_fk
!-------------------------------------------------------------------
! Purpose: set up constants factors for f_-k(day) used for empirical model
!     to calculate the electric potential
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!-------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:

      ft(1,0) = .75*sqrt( 6.e0 )/pi			
      ft(1,1) = 2.e0*ft(1,0)					      
      ft(1,2) = 1.e0						      
      ft(2,0) = ft(1,0) 					      
      ft(2,1) = -ft(1,1)					      
      ft(2,2) = 1.e0						      
      ft(3,0) = ft(2,1) 					      
      ft(3,1) = 0.						      
      ft(3,2) = 1.e0							   

      end subroutine prep_fk
!-----------------------------------------------------------------------
      end module module_prep_fk
