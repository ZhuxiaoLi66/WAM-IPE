! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
!NOTE: dbg20120125: only temporary used to switch on the transport only during the daytime...
      subroutine activate_perp_transport (utime,mp)
        USE module_precision
        USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_grid_3d,JMIN_IN,JMAX_IS,ISL,IBM,IGR,IQ,IGCOLAT,IGLON
        USE module_input_parameters,ONLY: sw_perp_transport
        USE module_physical_constants,ONLY: pi
        IMPLICIT NONE
        INTEGER (KIND=int_prec),  INTENT(IN)  :: utime    !universal time [sec]
        INTEGER (KIND=int_prec),  INTENT(IN)  :: mp
!---local variables---
        INTEGER (KIND=int_prec) :: lpj,midpoint
        REAL (KIND=real_prec) :: ltime
        REAL (KIND=real_prec),PARAMETER :: ltime_min= 9.00
        REAL (KIND=real_prec),PARAMETER :: ltime_max=17.00
!
!jicamarca
      lpj=130
! calculate LT at midpoint
      midpoint = JMIN_IN(lpj) + ( JMAX_IS(lpj) - JMIN_IN(lpj) )/2
!nm20110909: calculating LT should be made FUNCTION!!!
      ltime = REAL(utime,real_prec)/3600.0 + (plasma_grid_3d(midpoint,mp,IGLON)*180.0/pi)/15.0
      IF ( ltime > 24.0 )  ltime = MOD(ltime, 24.0)
       
      IF ( ltime>=ltime_min .AND. ltime<=ltime_max ) THEN
        sw_perp_transport(mp)=1
        print *, 'activating sw_perp_transport: mp=',mp,' LT=',ltime
      END IF

      END subroutine activate_perp_transport
