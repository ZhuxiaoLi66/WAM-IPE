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
      MODULE module_close_files
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPEV
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_io.f90!

      PRIVATE
      PUBLIC :: close_files

      CONTAINS
!---------------------------
        SUBROUTINE close_files ( )
        USE module_IO,ONLY: PRUNIT,LUN_UT,lun_min1,lun_max1,LUN_PLASMA1 &
, lun_ipe_grid_neut_params_ut &
, lun_ipe_grid_neut_wind &
, lun_ipe_grid_neut_temp &
, lun_ipe_grid_neut_O_den &
, lun_ipe_grid_neut_N2_den &
, lun_ipe_grid_neut_O2_den
USE module_input_parameters,ONLY:sw_output_wind
        IMPLICIT NONE
        INTEGER (KIND=int_prec) :: i
!
        CLOSE(UNIT=PRUNIT)
        CLOSE(UNIT=LUN_UT)

DO i=lun_min1,lun_max1
        CLOSE( UNIT=LUN_PLASMA1(i) )
END DO

!nm20141001
        IF ( sw_output_wind ) THEN
          CLOSE(UNIT=lun_ipe_grid_neut_params_ut)
          CLOSE(UNIT=lun_ipe_grid_neut_wind)
          CLOSE(UNIT=lun_ipe_grid_neut_temp)
          CLOSE(UNIT=lun_ipe_grid_neut_O_den)
          CLOSE(UNIT=lun_ipe_grid_neut_N2_den)
          CLOSE(UNIT=lun_ipe_grid_neut_O2_den)
        END IF !( sw_output_wind ) THEN

        END SUBROUTINE close_files
!---------------------------
END MODULE module_close_files
