!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_LEGENDRE
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
      public :: LEGENDRE

      contains
!-----------------------------------------------------------------------
	SUBROUTINE LEGENDRE(x,lmax,mmax,Plm)
!
!-----------------------------------------------------------------------
! compute Associate Legendre Function P_l^m(x)
! for all l up to lmax and all m up to mmax.
! returns results in array Plm
! if X is out of range ( abs(x)>1 ) then value is returned as if x=1.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
c       use cam_logfile,  only : iulog

        implicit none 
!
!------------------------------Arguments--------------------------------
!
        integer lmax, mmax
	real x, Plm(0:20,0:20)
!
!---------------------------Local variables-----------------------------
!
        integer m, lm2, l, iulog
        real xx, fact
        iulog=14
!
!-----------------------------------------------------------------------
!
	  DO l=0,20
	    DO m=0,20
		Plm(l,m)=0.
	    ENDDO
	  ENDDO
	xx=MIN(x,1.)
	xx=MAX(xx,-1.)
c         IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
c        write(iulog,*)'Bad arguments to Legendre'
c        RETURN
c        ENDIF
! First calculate all Pl0 for l=0 to l
	Plm(0,0)=1.
	IF(lmax.GT.0)Plm(1,0)=xx
	IF (lmax .GT. 1 )THEN
	  DO L=2,lmax
	    Plm(L,0)=( (2.*L-1)*xx*Plm(L-1,0) - 
     &(L-1)*Plm(L-2,0) )/L
	  ENDDO
	ENDIF
	IF (mmax .EQ. 0 )RETURN
	fact=SQRT( (1.-xx)*(1.+xx) )
	DO M=1,mmax
	  DO L=m,lmax
	    lm2=MAX(L-2,0)
	    Plm(L,M)=Plm(lm2,M) - ( 2*L-1)*fact*Plm(L-1,M-1)
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE LEGENDRE

!-----------------------------------------------------------------------
      end module module_LEGENDRE
