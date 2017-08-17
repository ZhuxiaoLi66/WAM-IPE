!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_ADJUST
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: ADJUST

      contains
!-----------------------------------------------------------------------
	SUBROUTINE ADJUST(ANGLE)
!
!-----------------------------------------------------------------------
!	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
        real angle
!
!-----------------------------------------------------------------------
!
 10	CONTINUE
	IF(ANGLE.LT.0.)THEN
	  ANGLE=ANGLE+360.
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360.)THEN
	  ANGLE=ANGLE-360.
	  GOTO 20
	ENDIF
	RETURN
	END SUBROUTINE ADJUST

!-----------------------------------------------------------------------
      end module module_ADJUST
