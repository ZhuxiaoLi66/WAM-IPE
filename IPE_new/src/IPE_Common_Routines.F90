MODULE IPE_Common_Routines

USE IPE_Precision

USE netcdf

IMPLICIT NONE

CONTAINS

  LOGICAL FUNCTION Almost_Equal( a, b ) 
    IMPLICIT NONE
    REAL(prec) :: a, b


      IF( a == 0.0_prec .OR. b == 0.0_prec )THEN
   
        IF( ABS(a-b) <= EPSILON(1.0_prec) )THEN

          Almost_Equal = .TRUE.

        ELSE

          Almost_Equal = .FALSE.
 
        ENDIF

      ELSE

        IF( (ABS(a-b) <= EPSILON(1.0_prec)*ABS(a)) .OR. (ABS(a-b) <= EPSILON(1.0_prec)*ABS(b)) )THEN
          
          Almost_Equal = .TRUE.

        ELSE

          Almost_Equal = .FALSE.

        ENDIF

      ENDIF
    
  END FUNCTION Almost_Equal

  INTEGER FUNCTION NewUnit(thisunit)
    IMPLICIT NONE
    INTEGER, INTENT(out), optional :: thisunit
    ! Local
    INTEGER, PARAMETER :: unitMin=100, unitMax=1000
    LOGICAL :: isopened
    INTEGER :: iUnit
 
      newunit=-1
 
      DO iUnit=unitMin, unitMax
 
         INQUIRE(UNIT=iUnit,opened=isopened)
 
         if( .not. isopened )then
            NewUnit = iUnit
            EXIT
         ENDif
 
      ENDDO
 
      IF( PRESENT(thisunit) ) thisunit = NewUnit
  
  END FUNCTION NewUnit

  SUBROUTINE Check(status)
    IMPLICIT NONE
    INTEGER, INTENT (in) :: status
    
     IF(status /= nf90_noerr) THEN 
       PRINT *, trim(nf90_strerror(status))
       STOP "NetCDF Error, Stopped"
     ENDIF
  END SUBROUTINE Check 

END MODULE IPE_Common_Routines
