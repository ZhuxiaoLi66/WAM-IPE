MODULE IPE_Model_Class

USE IPE_Precision
USE IPE_Model_Parameters_Class
USE IPE_Grid_Class
USE IPE_Neutrals_Class
!USE IPE_Forcing_Class
!USE IPE_Plasma_Class
!USE IPE_Electrodynamics_Class


IMPLICIT NONE


  ! The IPE_Model serves as a wrapper for all of the underlying attributes.
  ! This class should be used to orchestrate model setup, updating, and
  ! breakdown. At this level, we define the API for interacting with the
  ! deeper attributes within IPE.

  TYPE IPE_Model

    TYPE( IPE_Model_Parameters ) :: parameters
    TYPE( IPE_Grid )             :: grid
    !TYPE( IPE_Forcing )          :: forcing
    TYPE( IPE_Neutrals )         :: neutrals
    !TYPE( IPE_Plasma )           :: plasma
    !TYPE( IPE_Eldyn )            :: eldyn

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model
      PROCEDURE :: Trash => Trash_IPE_Model

!      PROCEDURE :: Write_NetCDF_IPE
!      PROCEDURE :: Read_NetCDF_IPE

!      PROCEDURE :: Write_Geographic_NetCDF_IPE
   
  END TYPE IPE_Model

CONTAINS

  SUBROUTINE Build_IPE_Model( ipe )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(out) :: ipe
    ! Local 
    LOGICAL :: init_success

      CALL ipe % parameters % Build( init_success )

      IF( init_success )THEN

        CALL ipe % grid % Read_IPE_Grid_NetCDF( ipe % parameters % netcdf_grid_file  )

#ifdef COUPLED_TO_WAM
        CALL ipe % neutrals % Build( nFluxtube       = ipe % grid % nFluxTube, &
                                     NLP             = ipe % grid % NLP, &
                                     NMP             = ipe % grid % NMP, &
                                     nCouplingFields = 7 )
#else
        CALL ipe % neutrals % Build( nFluxtube       = ipe % grid % nFluxTube, &
                                     NLP             = ipe % grid % NLP, &
                                     NMP             = ipe % grid % NMP )
#endif
                                      
      ELSE

        STOP
        
      ENDIF

  END SUBROUTINE Build_IPE_Model
!
  SUBROUTINE Trash_IPE_Model( ipe )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe


      CALL ipe % grid % Trash( )
      CALL ipe % neutrals % Trash( )

  END SUBROUTINE Trash_IPE_Model


END MODULE IPE_Model_Class
