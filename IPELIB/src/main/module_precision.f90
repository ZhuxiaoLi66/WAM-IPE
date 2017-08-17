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
      MODULE module_precision
      IMPLICIT NONE
      INTEGER, PARAMETER, PUBLIC :: real_prec4 = SELECTED_REAL_KIND(6,37) !=4 
      INTEGER, PARAMETER, PUBLIC :: real_prec8 = SELECTED_REAL_KIND(15,300) !=8
      INTEGER, PARAMETER, PUBLIC :: real_prec  = 4 !SELECTED_REAL_KIND(15,300) !=8
      INTEGER, PARAMETER, PUBLIC :: int_prec  = SELECTED_INT_KIND(9) !=4
      END MODULE module_precision
