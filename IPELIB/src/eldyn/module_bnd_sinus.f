!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_bnd_sinus
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
      public :: bnd_sinus

      contains
!-----------------------------------------------------------------------
      subroutine bnd_sinus( ihlat_bnd, itrans_width )  
!------------------------------------------------------------------
! Purpose: 
!   1. adjust high latitude boundary (ihlat_bnd) sinusoidally
!   2. width of transition zone from midlatitude potential to high latitude
!      potential (itrans_width)
!
! Method:
! 1.adjust boundary sinusoidally
!   max. wave number to be represented nmax_sin
!   RHS(mi) = Sum_phi Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi)*hlat_bnd(phi) 
!   U(mi,mk)   = Sum_phi Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi) *
!                Sum_(mk=-nmax_sin)^_(mk=nmax_sin) f_mk(phi)
!   single values decomposition of U
!   solving U*LSG = RHS 
!   calculating hlat_bnd:
!   hlat_bnd = Sum_(mi=-nmax_sin)^_(mi=nmax_sin) f_mi(phi)*LSG(mi)
!
! 2. width of transition zone from midlatitude potential to high latitude
!    potential
!    trans_width(phi)=8.-2.*cos(phi) 
!
! Author: A. Maute Nov 2003  am 11/20/03
!------------------------------------------------------------------

c     use sv_decomp, only : svdcmp, svbksb
!nm20121003
      USE efield !,ONLY:  
      USE module_svdcmp ,ONLY: svdcmp 
      USE module_ff ,ONLY:ff   
      USE module_svbksb ,ONLY: svbksb
!----------------------------------------------------------------------------                                                                   
!	... dummy arguments
!----------------------------------------------------------------------------                                                                   
      integer, intent(inout) :: ihlat_bnd(0:nmlon)    ! loaction of boundary
      integer, intent(out)   :: itrans_width(0:nmlon) ! width of transition zone

!-----------------------------------------------------------------
! local variables
!-----------------------------------------------------------------
      integer, parameter :: nmax_a = 2*nmax_sin+1 ! absolute array length
      integer, parameter :: ishf   = nmax_sin+1
      integer  :: ilon, i, i1, j, bnd
      real :: sum, mlon
      real :: rhs(nmax_a)
      real :: lsg(nmax_a)
      real :: u(nmax_a,nmax_a)
      real :: v(nmax_a,nmax_a)
      real :: w(nmax_a,nmax_a)
      real :: f(-nmax_sin:nmax_sin,0:nmlon)

!------------------------------------------------------------------
!    Sinusoidal Boundary calculation
!------------------------------------------------------------------
      rhs(:) = 0.
      lsg(:) = 0.
      u(:,:) = 0.
      v(:,:) = 0.
      w(:,:) = 0.

      do ilon = 0,nmlon                  ! long.
        bnd  = nmlath - ihlat_bnd(ilon) ! switch from pole=0 to pole =90
        call ff( ylonm(ilon), nmax_sin, f(-nmax_sin,ilon) )
        do i = -nmax_sin,nmax_sin
	  i1 = i + ishf
          rhs(i1) = rhs(i1) + f(i,ilon) * bnd
!	  write(iulog,*) 'rhs ',ilon,i1,bnd,f(i,ilon),rhs(i+ishf)
          do j = -nmax_sin,nmax_sin 
            u(i1,j+ishf) = u(i1,j+ishf) + f(i,ilon)*f(j,ilon)
!	    write(iulog,*) 'u ',ilon,i1,j+ishf,u(i+ishf,j+ishf)
          end do
        end do
      end do

!     if (debug) write(iulog,*) ' Single Value Decomposition'
      call svdcmp( u, nmax_a, nmax_a, nmax_a, nmax_a, w, v )

!     if (debug) write(iulog,*) ' Solving'
      call svbksb( u, w, v, nmax_a, nmax_a, nmax_a, nmax_a, rhs, lsg )
!      
      do ilon = 0,nmlon  ! long.
!       sum = 0.
	sum = dot_product( lsg(-nmax_sin+ishf:nmax_sin+ishf),
     &f(-nmax_sin:nmax_sin,ilon) )
!       do i = -nmax_sin,nmax_sin
!         sum = sum + lsg(i+ishf)*f(i,ilon)  
!       end do
        ihlat_bnd(ilon)    = nmlath - int( sum + .5 )                                ! closest point
        itrans_width(ilon) = 
     &int( 8. - 2.*cos( ylonm(ilon)*dtr ) + .5 )/dlatm  ! 6 to 10 deg.
      end do      
!     write(iulog,"('bnd_sinus: ihlat_bnd=',/,(12i6))") ihlat_bnd
!     write(iulog,"('bnd_sinus: itrans_width=',/,(12i6))") itrans_width

      end subroutine bnd_sinus
!-----------------------------------------------------------------------
      end module module_bnd_sinus
