PROGRAM IPE_Driver

USE IPE_Precision
USE IPE_Model_Class

IMPLICIT NONE


  TYPE( IPE_Model ) :: ipe
  LOGICAL           :: init_success
  REAL(prec)        :: t0_ipe, t1_ipe
  REAL(prec)        :: t2, t1
  INTEGER           :: i

    CALL ipe % Build( init_success )

    IF( init_success )THEN

      DO i = 1, ipe % parameters % n_model_updates 

        CALL CPU_TIME(t1)
        t0_ipe = ipe % parameters % start_time + REAL(i-1,prec)*ipe % parameters % file_output_frequency
        t1_ipe = ipe % parameters % start_time + REAL(i,prec)*ipe % parameters % file_output_frequency

        CALL ipe % Update( t0_ipe, t1_ipe )
        CALL CPU_TIME(t2)

!       CALL ipe % Write_NetCDF_IPE( "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".nc" ) 
        
        ! Interpolate to the Geographic Fixed Height Grid and write to file
!       CALL ipe % Write_Geographic_NetCDF_IPE( "IPE_State.geo."//ipe % time_tracker % DateStamp( )//".nc" ) 

        PRINT*, 'Total Update Time :', t2-t1

      ENDDO

    ENDIF

END PROGRAM IPE_Driver
