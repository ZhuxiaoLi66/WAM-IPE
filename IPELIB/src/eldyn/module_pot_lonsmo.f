!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_pot_lonsmo
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: pot_lonsmo

      contains
!-----------------------------------------------------------------------
      subroutine pot_lonsmo( pot, idlon ) 
!-------------------------------------------------------------------
! Purpose: smoothing in longitude of potential
!
! Method:  weighted smoothing in longitude
!          assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(in)     :: idlon
      real, intent(inout) :: pot(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: ilon, ilat, id, iabs
      real :: wgt, del
      real :: w(-idlon:idlon)
      real :: pot_smo(0:nmlath) ! temp array for smooth. potential
      real :: tmp(-idlon:nmlon+idlon) ! temp array for smooth. potential

!-------------------------------------------------------------------
! weighting factors (regular grid spacing) 
!-------------------------------------------------------------------
      wgt = 0.
      do id = -idlon,idlon
        del   = abs(id)*dlonm	! delta lon_m
        w(id) = 1./(del + 1.)
        wgt   = wgt + w(id)
      end do
        wgt = 1./wgt

!-------------------------------------------------------------------
! averaging     
!-------------------------------------------------------------------
!     do ilon = 0,nmlon
!       do ilat = 0,nmlath
!         pot_smo(ilat) = 0.
!         do id = -idlon,idlon	                  !  org. was degree now grid points
!           iabs = ilon + id
!           if( iabs > nmlon ) then
!              iabs = iabs - nmlon ! test if wrap around
!           end if
!           if( iabs < 0 ) then
!              iabs = iabs + nmlon ! test if wrap around
!           end if
!           pot_smo(ilat) = pot_smo(ilat) + w(id)*pot(iabs,ilat)
!         end do
!         pot_smo(ilat)  = pot_smo(ilat)*wgt
!         pot(ilon,ilat) = pot_smo(ilat)       ! copy back into pot 
!         pot(ilon,nmlat-ilat) = pot_smo(ilat) ! copy back into pot    
!       end do   
!     end do

!$omp parallel do private(ilat,ilon,tmp)
      do ilat = 0,nmlath
          tmp(0:nmlon)             = pot(0:nmlon,ilat)
          tmp(-idlon:-1)           = pot(nmlon-idlon:nmlon-1,ilat)
          tmp(nmlon+1:nmlon+idlon) = pot(1:idlon,ilat)
          do ilon = 0,nmlon
      pot(ilon,ilat) = dot_product( tmp(ilon-idlon:ilon+idlon),w )*wgt
      pot(ilon,nmlat-ilat) = pot(ilon,ilat)
          end do   
      end do
      
      end subroutine pot_lonsmo
!-----------------------------------------------------------------------
      end module module_pot_lonsmo
