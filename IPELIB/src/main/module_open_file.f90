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
      MODULE module_open_file
      USE module_precision
      USE module_IPE_dimension   ,ONLY: ISPEC,ISPEV
      USE module_input_parameters,ONLY: mype
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_io.f90!

      PRIVATE
      PUBLIC :: open_file

      CONTAINS
!---------------------------
      SUBROUTINE open_file ( filename_dum, UNIT_dum, FORM_dum, STATUS_dum )  
      USE module_precision
      IMPLICIT NONE
      CHARACTER (LEN=*)      , INTENT(IN) :: filename_dum
      INTEGER (KIND=int_prec), INTENT(IN) :: UNIT_dum
      CHARACTER (LEN=*)      , INTENT(IN) :: FORM_dum
      CHARACTER (LEN=*)      , INTENT(IN) :: STATUS_dum

      OPEN(UNIT=UNIT_dum,FILE=TRIM(filename_dum),STATUS=STATUS_dum,FORM=FORM_dum)

      END SUBROUTINE open_file
END MODULE module_open_file
