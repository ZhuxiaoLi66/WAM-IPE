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
!SMS$IGNORE BEGIN
!#ifdef TESTING
! USE module_MDI
 USE module_eldyn
!#endif
!SMS$IGNORE END
   IMPLICIT NONE
   INCLUDE "gptl.inc"
!SMS$INSERT   include "mpif.h"


   INTEGER(KIND=int_prec)           :: utime_driver ! Universal Time [sec]
   INTEGER(KIND=int_prec),parameter :: luntmp=300   !
   INTEGER(KIND=int_prec)           :: istat,mp,ret ! 
   INTEGER(KIND=int_prec)           :: iterate
   CHARACTER(8)                     :: iterChar


!SMS$IGNORE BEGIN
!#ifdef TESTING
!     CALL mdi % Build( )
!#endif
!SMS$IGNORE END
     CALL gptlprocess_namelist ('GPTLnamelist', 77, ret) 
     ret = gptlinitialize ()
     ret = gptlstart ('Total')


!SMS$INSERT parallelBuild=.true.
! set up input parameters
     ret = gptlstart ('read_input')
!SMS$INSERT         MPI_COMM_IPE = MPI_COMM_WORLD
     CALL read_input_parameters ( )
     ret = gptlstop  ('read_input')

! open Input/Output files
     ret = gptlstart ('open_output_files')
!SMS$SERIAL BEGIN
     CALL open_output_files ( )
!SMS$SERIAL END
     ret = gptlstop  ('open_output_files')

! set up plasma grids by reading file
     ret = gptlstart ('init_plasma_grid')
     CALL init_plasma_grid ( )
     ret = gptlstop  ('init_plasma_grid')

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-1")

     IF ( sw_output_plasma_grid ) THEN
       ret = gptlstart ('output_plasma_grid')
       PRINT *, 'sub-init_p: output plasma_grid'
       CALL output_plasma_grid ( )
       ret = gptlstop  ('output_plasma_grid')
     END IF

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-2")

! initialise the flux tubes from previous runs
     IF ( HPEQ_flip==0.0 ) THEN
       PRINT *,'before CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
       ret = gptlstart ('io_plasma_bin')
       CALL io_plasma_bin ( 2, start_time )
       ret = gptlstop  ('io_plasma_bin')
       PRINT *,'after CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
     END IF


! initialization of electrodynamic module:
! read in E-field

     ret = gptlstart ('init_eldyn')
     IF ( sw_perp_transport>=1 ) THEN
       CALL init_eldyn ( )
     ENDIF
     ret = gptlstop  ('init_eldyn')

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-4")

     ret = gptlstart ('time_loop')

     iterate = 0

!dbg20170916
print*,'driver:start_time=',start_time, stop_time, time_step
     DO utime_driver = start_time, stop_time, time_step
       iterate = iterate + 1

       PRINT*,'driver:utime_driver=',utime_driver

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-5")


!SMS$IGNORE BEGIN
!#ifdef TESTING

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "eldyn", &
!                          "ed1_90 before eldyn", &
!                          0, &
!                          SIZE(ed1_90), &
!                          PACK(ed1_90,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "eldyn", &
!                          "ed2_90 before eldyn", &
!                          0, &
!                          SIZE(ed2_90), &
!                          PACK(ed2_90,.TRUE.) )   
!
!#endif
!SMS$IGNORE END
       ret = gptlstart ('eldyn')
       IF ( sw_perp_transport>=1.AND. MOD( (utime_driver-start_time),ip_freq_eldyn)==0 ) THEN
!dbg20170916
print*,'driver:before eldyn, utime_driver',utime_driver
         CALL eldyn ( utime_driver )
       ENDIF
       ret = gptlstop  ('eldyn')

!SMS$IGNORE BEGIN
!#ifdef TESTING
!       CALL mdi % Update( "driver_ipe.f90", &
!                          "eldyn", &
!                          "ed1_90 after eldyn", &
!                          0, &
!                          SIZE(ed1_90), &
!                          PACK(ed1_90,.TRUE.) )   
!
!       CALL mdi % Update( "driver_ipe.f90", &
!                          "eldyn", &
!                          "ed2_90 after eldyn", &
!                          0, &
!                          SIZE(ed2_90), &
!                          PACK(ed2_90,.TRUE.) )   
!#endif
!SMS$IGNORE end

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-6")

!SMS$IGNORE BEGIN
!#ifdef TESTING
!        
!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "ON_m3 before neutral", &
!                          0, &
!                          SIZE(ON_m3), &
!                          PACK(ON_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "HN_m3 before neutral", &
!                          0, &
!                          SIZE(HN_m3), &
!                          PACK(HN_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "N2N_m3 before neutral", &
!                          0, &
!                          SIZE(N2N_m3), &
!                          PACK(N2N_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "O2N_m3 before neutral", &
!                          0, &
!                          SIZE(O2N_m3), &
!                          PACK(O2N_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "N4S_m3 before neutral", &
!                          0, &
!                          SIZE(N4S_m3), &
!                          PACK(N4S_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "TN_k before neutral", &
!                          0, &
!                          SIZE(TN_k), &
!                          PACK(TN_k,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "TINF_k before neutral", &
!                          0, &
!                          SIZE(TINF_k), &
!                          PACK(TINF_k,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "Un_ms1 before neutral", &
!                          0, &
!                          SIZE(Un_ms1), &
!                          PACK(Un_ms1,.TRUE.) )   

!#endif
!SMS$IGNORE END
        ! update neutral 3D structure: 
        ! use MSIS/HWM to get the values in the flux tube grid
       IF ( MOD( (utime_driver-start_time),ip_freq_msis)==0 ) THEN 

!SMS$IGNORE BEGIN
!#ifdef DEBUG
         PRINT *,'CALL MSIS',utime_driver,start_time, &
                  ip_freq_msis,(utime_driver-start_time), &
                  MOD( (utime_driver-start_time),ip_freq_msis)
!#endif
!SMS$IGNORE END
         ret = gptlstart ('neutral')
         CALL neutral ( utime_driver )
         ret = gptlstop  ('neutral')

       END IF

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-7")

!SMS$IGNORE BEGIN
!#ifdef TESTING
        
!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "ON_m3 after neutral", &
!                          0, &
!                          SIZE(ON_m3), &
!                          PACK(ON_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "HN_m3 after neutral", &
!                          0, &
!                          SIZE(HN_m3), &
!                          PACK(HN_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "N2N_m3 after neutral", &
!                          0, &
!                          SIZE(N2N_m3), &
!                          PACK(N2N_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "O2N_m3 after neutral", &
!                          0, &
!                          SIZE(O2N_m3), &
!                          PACK(O2N_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "N4S_m3 after neutral", &
!                          0, &
!                          SIZE(N4S_m3), &
!                          PACK(N4S_m3,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "TN_k after neutral", &
!                          0, &
!                          SIZE(TN_k), &
!                          PACK(TN_k,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "TINF_k after neutral", &
!                          0, &
!                          SIZE(TINF_k), &
!                          PACK(TINF_k,.TRUE.) )   

!       CALL mdi % Update( "driver_ipe.f90", &
!                          "neutral", &
!                          "Un_ms1 after neutral", &
!                          0, &
!                          SIZE(Un_ms1), &
!                          PACK(Un_ms1,.TRUE.) )   

!#endif
!SMS$IGNORE END

!SMS$IGNORE BEGIN
!#ifdef TESTING
!       CALL mdi % Update( "driver_ipe.f90", &
!                          "plasma", &
!                          "plasma_3d begin plasma", &
!                          0, &
!                          SIZE(plasma_3d), &
!                          PACK(plasma_3d,.TRUE.) )   
!#endif
!SMS$IGNORE END

! update plasma
        ret = gptlstart ('plasma')
!ghgm - a dummy timestamp (13 characters) needs to be here
! because we use timestamps in the fully coupleid WAM-IPE
! Obviously needs a better solution.....
        CALL plasma ( utime_driver  )
        ret = gptlstop  ('plasma')

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-8")

!SMS$IGNORE BEGIN
!#ifdef TESTING
!       CALL mdi % Update( "driver_ipe.f90", &
!                          "plasma", &
!                          "plasma_3d end plasma", &
!                          0, &
!                          SIZE(plasma_3d), &
!                          PACK(plasma_3d,.TRUE.) )   
!#endif
!SMS$IGNORE END

       IF( MOD(utime_driver,ip_freq_output)==0)THEN
          WRITE( iterChar, '(I8.8)' )utime_driver
          CALL io_plasma_bin ( 1, utime_driver, 'iter_'//iterChar )
       ENDIF
       ret = gptlstart ('output')
       CALL output ( utime_driver )
       ret = gptlstop  ('output')

!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-9")

!SMS$IGNORE BEGIN
!#ifdef TESTING
!       CALL mdi % Write_ModelDataInstances( "ipe" )
!       CALL mdi % CalculateStorageCost(  )
!#endif
!SMS$IGNORE END
     END DO

     ret = gptlstop  ('time_loop')

    ! Deallocate arrays
     ret = gptlstart ('allocate_arrays1')
     CALL allocate_arrays ( 1 )
     ret = gptlstop  ('allocate_arrays1')

     ! close all open files
     ret = gptlstart ('close_files')
     CALL close_files ( )
     ret = gptlstop  ('close_files')


     ret = gptlstop  ('Total')
!SMS$IGNORE BEGIN
!#ifdef TESTING
!     CALL mdi % Trash( )
!#endif
!SMS$IGNORE END
     CALL stop


END PROGRAM  test_plasma
