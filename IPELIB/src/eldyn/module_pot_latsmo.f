!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_pot_latsmo
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: pot_latsmo

      contains
!-----------------------------------------------------------------------
      subroutine pot_latsmo( pot, idlat )  ! pots == pot_highlats
!--------------------------------------------------------------------
! Purpose: smoothing in latitude of  potential
!
! Method: weighted smoothing in latitude 
! assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(in)     :: idlat
      real, intent(inout) :: pot(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!------------------------------------------------------------------
      integer  :: ilon, ilat, id
      real :: wgt, del
      real :: w(-idlat:idlat)
!     real :: pot_smo(0:nmlat) ! temp array for smooth. potential
      real :: pot_smo(0:nmlon,0:nmlat_wei) ! temp array for smooth. potential

!------------------------------------------------------------------
! weighting factors (regular grid spacing) 
!------------------------------------------------------------------
      wgt = 0. 
      do id = -idlat,idlat
	del   = abs(id)*dlatm	! delta lat_m
	w(id) = 1./(del + 1.)
        wgt   = wgt + w(id)
      end do
      wgt = 1./wgt

!     do ilon = 0,nmlon
!        do ilat = idlat,nmlat_wei-idlat
!       do ilat = idlat,nmlat-idlat
!         pot_smo(ilat) = 0.
!         do id = -idlat,idlat	!  org. was degree now grid points
!           pot_smo(ilat) = pot_smo(ilat) + w(id)*pot(ilon,ilat+id)
!           write(iulog,"('pot_latsmo: ilon=',i3,' ilat=',i3,' id=',i3,' pot(ilon,ilat+id)=',e12.4)") ilon,ilat,id,pot(ilon,ilat+id)
!         end do
!         pot_smo(ilat)       = pot_smo(ilat)*wgt
!         pot_smo(nmlat-ilat) = pot_smo(ilat)
!       end do
!      pot(ilon,idlat:nmlat-idlat) =  &        ! copy back into pot
!         pot_smo(idlat:nmlat-idlat)
!        pot(ilon,idlat:nmlat_wei-idlat)       = pot_smo(idlat:nmlat_wei-idlat)
!       pot(ilon,nmlat-nmlat_wei+idlat:nmlat) = pot_smo(nmlat-nmlat_wei+idlat:nmlat)
!        pot(ilon,nmlat-nmlat_wei+idlat:nmlat-idlat) = pot_smo(nmlat-nmlat_wei+idlat:nmlat-idlat)
!     end do

!$omp parallel do private(ilat)
      do ilat = idlat,nmlat_wei-idlat
         pot_smo(:,ilat) = matmul( pot(:,ilat-idlat:ilat+idlat),w )*wgt
      end do

      do ilat = idlat,nmlat_wei-idlat
         pot(:,ilat)       = pot_smo(:,ilat)
         pot(:,nmlat-ilat) = pot_smo(:,ilat)
      end do

      end subroutine pot_latsmo
!-----------------------------------------------------------------------
      end module module_pot_latsmo
