MODULE IPE_Common_Routines

USE IPE_Precision


IMPLICIT NONE

CONTAINS

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


END MODULE IPE_Common_Routines
