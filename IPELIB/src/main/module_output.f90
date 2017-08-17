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
      MODULE module_output
      USE module_precision
      USE module_IPE_dimension   ,ONLY: ISPEC,ISPEV
      USE module_input_parameters,ONLY: mype
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_io.f90!

      PRIVATE
      PUBLIC :: output

      CONTAINS
!---------------------------
        SUBROUTINE output ( utime )
        USE module_IO,ONLY: PRUNIT
        IMPLICIT NONE
        include "gptl.inc"
        integer ret
!------------------------
        INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]

        if(mype==0) then
          WRITE(UNIT=PRUNIT,FMT="('uts=',i7,i7)") utime,prunit
        endif

        END SUBROUTINE output
!---------------------------
END MODULE module_output
