!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_DerivPotential
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: DerivPotential

      contains
!-----------------------------------------------------------------------
      subroutine DerivPotential
!-----------------------------------------------------------------
! Purpose: calulates the electric field [V/m] from the electric potential
!
! Method:  Richmond [1995] eqn 5.9-5.10
! ed1(:,:) = Ed1 = - 1/[R cos lam_m] d PHI/d phi_m
! ed2(:,:) = Ed2 = 1/R d PHI/d lam_m /sin I_m
! R = R_e + h_r we assume a reference height of 130 km which is also
! used in the TIEGCM code
!
! Author: A. Maute Dec 2003  am 12/16/03
!-----------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:

      integer  :: i, j, ip1f, ip2f, ip3f
      real :: coslm, r_t, fac, wrk
      real :: wrk1d(0:nmlon)

      r_t = r_e + h_r  ! earth radius + reference height [m]
!-----------------------------------------------------------------
! ed2= Ed2 is the equatorward/downward component of the electric field, at all 
! geomagnetic grid points (central differencing)
!-----------------------------------------------------------------
      fac = .5/(dlatm*dtr*r_t)
!$omp parallel do private(j, i, wrk )
      do j = 1,nmlath-1		! southern hemisphere
        wrk = fac/sinIm_mag(j)
        do i = 0,nmlon
          ed2(i,j) = (potent(i,j+1) - potent(i,j-1))*wrk
        end do
      end do

!$omp parallel do private(j, i, wrk )
      do j = nmlath+1,nmlat-1	! northern hemisphere
        wrk = fac/sinIm_mag(j)
        do i = 0,nmlon
          ed2(i,j) = (potent(i,j+1) - potent(i,j-1))*wrk
        end do
      end do

!-----------------------------------------------------------------------
! Interpolate of ed2 between between -12 <= lam_m <= 12 degrees:
!-----------------------------------------------------------------------
      wrk1d(:) = ed2(:,jmax) - ed2(:,jmin)
      do j = jmin+1,jmax-1
        fac = (ylatm(j) - ylatm(jmin))/(ylatm(jmax) - ylatm(jmin))
        do i = 0,nmlon
	    ed2(i,j) = ed2(i,jmin) + fac*wrk1d(i)
	end do
      end do

!-----------------------------------------------------------------------
! ed1= Ed1 is the zonal component of the electric field, at all 
! geomagnetic grid points (central differencing)
!-----------------------------------------------------------------------
      fac = .5/(dlonm*dtr*r_t)
!$omp parallel do private(j, i, wrk, coslm )
      do j = 1,nmlat-1
        coslm = ylatm(j) - 90.
        coslm = cos( coslm*dtr )
        wrk = fac/coslm
        do i = 1,nmlon-1
          ed1(i,j) = -(potent(i+1,j) - potent(i-1,j))*wrk
        end do
	i = 0
	ed1(i,j)  = -(potent(i+1,j) - potent(nmlon-1,j))*wrk
	ed1(nmlon,j) = ed1(i,j)
      end do

!-----------------------------------------------------------------------
! Poles:
!-----------------------------------------------------------------------
      do i = 0,nmlon
        ip1f = i + nmlon/4
        if( ip1f > nmlon ) then
           ip1f = ip1f - nmlon
        end if
        ip2f = i + nmlon/2
        if( ip2f > nmlon ) then
           ip2f = ip2f - nmlon
        end if
        ip3f = i + 3*nmlon/4
        if( ip3f > nmlon ) then
           ip3f = ip3f - nmlon
        end if
        ed1(i,0)=.25*(ed1(i,1)-ed1(ip2f,1)+ed2(ip1f,1)-ed2(ip3f,1))
        ed1(i,nmlat) = .25*(ed1(i,nmlat-1) - ed1(ip2f,nmlat-1)  
     &        + ed2(ip1f,nmlat-1) - ed2(ip3f,nmlat-1))
        ed2(i,0)=.25*(ed2(i,1)-ed2(ip2f,1)-ed1(ip1f,1)+ed1(ip3f,1))
        ed2(i,nmlat) = .25*(ed2(i,nmlat-1) - ed2(ip2f,nmlat-1)  
     &        - ed1(ip1f,nmlat-1) + ed1(ip3f,nmlat-1))
      end do

      end subroutine DerivPotential
!-----------------------------------------------------------------------
      end module module_DerivPotential
