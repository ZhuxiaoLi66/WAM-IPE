!note:20120207: v36: used only activating the perp.transport gradually only during daytime...
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
!
PROGRAM  test_plasma

 USE module_precision
 USE module_input_parameters
 USE module_FIELD_LINE_GRID_MKS
 USE module_init_plasma_grid
 USE module_NEUTRAL_MKS 
 USE module_sub_PLASMA
!SMS$IGNORE BEGIN
 USE module_init_ELDYN
 USE module_sub_ELDYN
!SMS$IGNORE END
 USE module_open_output_files
 USE module_output
 USE module_close_files
 USE module_IPE_dimension
!SMS$IGNORE END
   IMPLICIT NONE
   INCLUDE "gptl.inc"
!SMS$INSERT   include "mpif.h"


   INTEGER(KIND=int_prec)           :: utime_driver ! Universal Time [sec]
   INTEGER(KIND=int_prec),parameter :: luntmp=300   !
   INTEGER(KIND=int_prec)           :: istat,mp,ret ! 
   INTEGER(KIND=int_prec)           :: utime_sub_loop, sub_time_loop, main_time_loop, n_time_steps
   CHARACTER(8)                     :: iterChar


     CALL gptlprocess_namelist ('GPTLnamelist', 77, ret) 
     ret = gptlinitialize ()
     ret = gptlstart ('Total')


!SMS$INSERT parallelBuild=.true.
! set up input parameters
!SMS$INSERT         MPI_COMM_IPE = MPI_COMM_WORLD
     CALL read_input_parameters ( )

! open Input/Output files
!SMS$SERIAL BEGIN
     CALL open_output_files ( )
!SMS$SERIAL END

! set up plasma grids by reading file
     CALL init_plasma_grid ( )


     IF ( sw_output_plasma_grid ) THEN
       CALL output_plasma_grid ( )
     END IF


! initialise the flux tubes from previous runs
     IF ( HPEQ_flip==0.0 ) THEN
       CALL io_plasma_bin ( 2, start_time )
     END IF

    ! Write the initial conditions to file
    WRITE( iterChar, '(I8.8)' )start_time
    CALL io_plasma_bin ( 1, start_time, 'iter_'//iterChar )


! initialization of electrodynamic module:
! read in E-field

     IF ( sw_perp_transport>=1 ) THEN
       CALL init_eldyn ( )
     ENDIF


     utime_driver = start_time
     n_time_steps = (stop_time - start_time)/(time_step)

     DO main_time_loop = 1, n_time_steps

       CALL neutral ( utime_driver, start_time )

       utime_sub_loop = utime_driver
       DO sub_time_loop = 1, time_step/solar_forcing_time_step

         CALL eldyn ( utime_sub_loop )

         CALL plasma ( utime_sub_loop )

         utime_sub_loop = utime_sub_loop + solar_forcing_time_step

       ENDDO

       utime_driver = utime_driver + time_step

       IF( MOD(utime_driver,ip_freq_output)==0)THEN
          WRITE( iterChar, '(I8.8)' )utime_driver
          CALL io_plasma_bin ( 1, utime_driver, 'iter_'//iterChar )
       ENDIF
       CALL output ( utime_driver )


     END DO


    ! Deallocate arrays
     CALL allocate_arrays ( 1 )

     ! close all open files
     CALL close_files ( )


END PROGRAM  test_plasma
