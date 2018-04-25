MODULE IPE_Time_Class

USE IPE_Precision

IMPLICIT NONE

  TYPE IPE_Time
    REAL(prec) :: utime
    INTEGER    :: year, month, day, hour, minute

    CONTAINS

    PROCEDURE :: Build => Build_IPE_Time
    PROCEDURE :: Set_Date => Set_Date_IPE_Time
    PROCEDURE :: Calculate_Hour_and_Minute => Calculate_Hour_and_Minute_IPE_Time
    PROCEDURE :: DateStamp => DateStamp_IPE_Time

  END TYPE IPE_Time



CONTAINS

  SUBROUTINE Build_IPE_Time( time_tracker, year, day_of_year, start_time )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    INTEGER, INTENT(in)              :: year, day_of_year
    REAL(prec), INTENT(in)           :: start_time

      time_tracker % utime = start_time 

      CALL time_tracker % Set_Date(  year, 0, day_of_year )

      CALL time_tracker % Calculate_Hour_and_Minute( )

  END SUBROUTINE Build_IPE_Time

  SUBROUTINE Set_Date_IPE_Time( time_tracker, year, month, day )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    INTEGER, INTENT(in)              :: year, month, day

      time_tracker % year  = year
      time_tracker % month = month
      time_tracker % day   = day


  END SUBROUTINE Set_Date_IPE_Time

  SUBROUTINE Calculate_Hour_and_Minute_IPE_Time( time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker

      time_tracker % hour = INT( time_tracker % utime/3600.0_prec )
      time_tracker % minute = INT( (time_tracker % utime - REAL(time_tracker % hour,prec)*3600.0_prec)/60.0_prec )

  END SUBROUTINE Calculate_Hour_and_Minute_IPE_Time
 
  FUNCTION DateStamp_IPE_Time( time_tracker ) RESULT( timestamp ) 
    CLASS( IPE_Time ) :: time_tracker
    CHARACTER(12)     :: timestamp

      WRITE( timestamp, '(I4,I2.2,I2.2,I2.2,I2.2)' ) time_tracker % year, &
                                                     time_tracker % month, &
                                                     time_tracker % day, &
                                                     time_tracker % hour, &
                                                     time_tracker % minute
                                                     
                                                   

  END FUNCTION DateStamp_IPE_Time

END MODULE IPE_Time_Class
