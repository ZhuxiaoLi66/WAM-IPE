!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_prep_weimer
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: prep_weimer

      contains
!-----------------------------------------------------------------------
      subroutine prep_weimer
!-----------------------------------------------------------------
! Purpose:  for Weimer model calculate IMF angle, IMF magnitude
!  tilt of earth
!
! Method: using functions and subroutines from Weimer Model 1996
!     output:  angle, &  ! IMF angle
!     	       bt,    &  ! IMF magnitude
!     	       tilt      ! tilt of earth
!
! Author: A. Maute Nov 2003  am 11/20/03
!-----------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      USE module_SetModel ,ONLY: SetModel
      USE module_GET_TILT ,ONLY: GET_TILT
      USE module_ADJUST ,ONLY: ADJUST
!-----------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------
      real ::  
     &  angle,  ! IMF angle
     &  bt,    ! IMF magnitude
     &  tilt       ! tilt of earth

!-----------------------------------------------------------------
! function declarations
!-----------------------------------------------------------------
!nm20121012      real, external :: get_tilt	 ! in wei96.f

      if( by == 0. .and. bz == 0.) then
         angle = 0.
      else
         angle = atan2( by,bz )
      end if
      
      angle = angle*rtd
      call adjust( angle )
      bt = sqrt( by*by + bz*bz )
!-------------------------------------------------------------------
! use month and day of month - calculated with average no.of days per month
! as in Weimer
!-------------------------------------------------------------------
c     if(debug) write(iulog,*) 'prep_weimer: day->day of month',
c    &iday,imo,iday_m,ut
      tilt = get_tilt( iyear, imo, iday_m, ut )

c      if(debug) then
c       write(iulog,"(/,'efield prep_weimer:')")
c       write(iulog,*)  '  Bz   =',bz
c       write(iulog,*)  '  By   =',by
c       write(iulog,*)  '  Bt   =',bt
c       write(iulog,*)  '  angle=',angle
c       write(iulog,*)  '  VSW  =',v_sw
c       write(iulog,*)  '  tilt =',tilt
c      end if

      call SetModel( angle, bt, tilt, v_sw )

      end subroutine prep_weimer
!-----------------------------------------------------------------------
      end module module_prep_weimer
