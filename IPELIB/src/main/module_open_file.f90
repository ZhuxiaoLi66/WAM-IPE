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
      USE module_IO,ONLY: LUN_LOG
      IMPLICIT NONE
      CHARACTER (LEN=*)      , INTENT(IN) :: filename_dum
      INTEGER (KIND=int_prec), INTENT(IN) :: UNIT_dum
      CHARACTER (LEN=*)      , INTENT(IN) :: FORM_dum
      CHARACTER (LEN=*)      , INTENT(IN) :: STATUS_dum
!---local
      LOGICAL                             :: flag
      INTEGER (KIND=int_prec)             :: istat

      IF (UNIT_dum>=1) THEN
        print *,'SUBROUTINE open_file, unit number = ',UNIT_dum
      ELSE
        print *,'!STOP! unit number not provided!!!'
        STOP
      END IF
      INQUIRE ( UNIT=UNIT_dum, OPENED=flag ) 
      WRITE( UNIT=LUN_LOG, FMT=*) flag
      IF ( .NOT. flag ) THEN
        WRITE( UNIT=LUN_LOG, FMT=*) 'opening file:',filename_dum
        OPEN(UNIT=UNIT_dum,FILE=TRIM(filename_dum),STATUS=STATUS_dum,FORM=FORM_dum,IOSTAT=istat)
        IF ( istat /= 0 ) THEN
          WRITE( UNIT=LUN_LOG, FMT=*)'ERROR OPENING FILE',filename_dum
          STOP
        END IF
      END IF      

      END SUBROUTINE open_file
END MODULE module_open_file
