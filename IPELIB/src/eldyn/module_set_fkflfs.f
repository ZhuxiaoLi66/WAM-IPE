!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_set_fkflfs
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: set_fkflfs

      contains
!-----------------------------------------------------------------------
      subroutine set_fkflfs( fk, fl, fs )
!------------------------------------------------------------------
! Purpose:  set f_-k(day) depending on seasonal flag used for empirical model
!     to calculate the electric potential
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!-----------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      USE module_ff ,ONLY:ff
!-----------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------
      real, intent(out) ::  
     &	fk(0:2),  	                ! f_-k(day) 
     &	fl(-2:2), 	                ! f_l(ut)  
     &	fs(2)		                ! f_s(f10.7) 
!------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: lp
      real :: ang
      real :: lon_ut

!------------------------------------------------------------------
! f_-k(day) 
! use factors for iseasav == 0 - Scherliess had iseasav as an input parameter
!------------------------------------------------------------------
      lp = iseasav
      if( iseasav == 0 ) then
        ang   = (day + 9.)*dy2rd
        fk(0) = sqr2*cos( 2.*ang )
        fk(1) = sqr2*cos( ang )
        fk(2) = 1.
      else if( iseasav >= 1 .and. iseasav <= 3 ) then
        fk(0) = ft(lp,0)
        fk(1) = ft(lp,1)
        fk(2) = ft(lp,2)
      else if( iseasav == 4 ) then
        fk(0) =0.
        fk(1) =0.
        fk(2) =1.
      end if

!-----------------------------------------------------------------
! f_l(ut) 
!-----------------------------------------------------------------
      lon_ut = 15.*ut        ! 15.*mlt - xmlon + 69. 
      call ff( lon_ut, 2, fl )                                                 
      if( iutav ) then  	! UT-averaging
     
	ang   = fl(0)
        fl(:) = 0.
        fl(0) = ang
	
      end if

!-----------------------------------------------------------------
! f_s(f10.7)  only fs(1) used  	
!-----------------------------------------------------------------
      fs(1) = 1.
!     fs(2) = S_a			  
      fs(2) = f107d			  

      end subroutine set_fkflfs
!-----------------------------------------------------------------------
      end module module_set_fkflfs
