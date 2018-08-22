PROGRAM Eldyn_Converter

USE IPE_Precision
USE IPE_Model_Class

IMPLICIT NONE


  TYPE( IPE_Model ) :: ipe
  LOGICAL           :: init_success
  REAL(prec)        :: t0_ipe
  INTEGER           :: i

    CALL ipe % Build( init_success )

    IF( init_success )THEN

      DO i = 1, ipe % parameters % n_model_updates 

        t0_ipe = ipe % parameters % start_time + REAL(i-1,prec)*ipe % parameters % file_output_frequency

        CALL ipe % time_tracker % Update( t0_ipe )
        

        ! Here is where the electrodynamics "Update" routine is called.
        ! You can feel free to comment it out, and replace it with calls
        ! to your "Read", "Regrid", and "Merge" routines.
        CALL ipe % eldyn % Update( ipe % grid, &
                                   ipe % forcing, &
                                   ipe % time_tracker )


        CALL ipe % Write_NetCDF_IPE( "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".nc" ) 
        
        ! Interpolate to the Geographic Fixed Height Grid and write to file
        CALL ipe % Write_Geographic_NetCDF_IPE( "IPE_State.geo."//ipe % time_tracker % DateStamp( )//".nc" ) 

      ENDDO

    ENDIF



CONTAINS


! In this CONTAINS region ( before the END PROGRAM statement and after
! "CONTAINS" ), you can add whatever subroutines/functions you want and can call
! them within this program. Think of this area as your scratch space for quickly
! prototyping routines that accomplish what  you want. Once you're happy with
! the routine, then we can work on pushing it into the
! IPE_Electrodynamics_Class.F90 module.

END PROGRAM Eldyn_Converter
