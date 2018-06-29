MODULE IPE_Electrodynamics_Class


USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Grid

IMPLICIT NONE

  TYPE IPE_Electrodynamics
    INTEGER    :: nFluxTube, NLP, NMP
    INTEGER    :: NP_dynamo, MP_dynamo
    REAL(prec), ALLOCATABLE :: electric_potential(:,:) 
    REAL(prec), ALLOCATABLE :: electric_field(:,:,:) 
    REAL(prec), ALLOCATABLE :: v_ExB_geographic(:,:,:)  ! "zonal" and "meridional" direction on the geographic grid
    REAL(prec), ALLOCATABLE :: v_ExB_apex(:,:,:) ! "zonal" and "meridional" direction ( VEXBth, VEXBe ) on the apex grid

    ! Inputs for the potential solver (on the dynamo grid)
    REAL(prec), ALLOCATABLE :: hall_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: pedersen_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: b_parallel_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: neutral_apex_velocity(:,:,:) ! Components are apex directions

  END TYPE IPE_Electrodynamics

CONTAINS


  SUBROUTINE Build_IPE_Electrodynamics( eldyn, nFluxTube, NLP, NMP, NP_dynamo, MP_dynamo ) 
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(out) :: eldyn
    INTEGER, INTENT(in)                       :: nFluxTube
    INTEGER, INTENT(in)                       :: NLP
    INTEGER, INTENT(in)                       :: NMP
    INTEGER, INTENT(in)                       :: NP_dynamo
    INTEGER, INTENT(in)                       :: MP_dynamo


      eldyn % nFluxTube = nFluxTube
      eldyn % NLP       = NLP
      eldyn % NMP       = NMP
      eldyn % NP_dynamo = NP_dynamo
      eldyn % MP_dynamo = MP_dynamo

      ALLOCATE( eldyn % electric_potential(1:NLP,1:NMP), &
                eldyn % electric_field(1:3,1:NLP,1:NMP), &
                eldyn % v_ExB_geographic(1:3,1:NLP,1:NMP), &
                eldyn % v_ExB_apex(1:3,1:NLP,1:NMP), &
                eldyn % hall_conductivity(1:NP_dynamo,1:MP_dynamo), &
                eldyn % pedersen_conductivity(1:NP_dynamo,1:MP_dynamo), &
                eldyn % b_parallel_conductivity(1:NP_dynamo,1:MP_dynamo), &
                eldyn % neutral_apex_velocity(1:3,1:NP_dynamo,1:MP_dynamo) )

      eldyn % electric_potential      = 0.0_prec
      eldyn % electric_field          = 0.0_prec
      eldyn % v_ExB_geographic        = 0.0_prec
      eldyn % v_ExB_apex              = 0.0_prec
      eldyn % hall_conductiviy        = 0.0_prec
      eldyn % pedersen_conductiviy    = 0.0_prec
      eldyn % b_parallel_conductivity = 0.0_prec
      eldyn % neutral_apex_velocity   = 0.0_prec



  END SUBROUTINE Build_IPE_Electrodynamics

  SUBROUTINE Trash_IPE_Electrodynamics( eldyn ) 
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn

      DEALLOCATE( eldyn % electric_potential, &
                  eldyn % electric_field, &
                  eldyn % v_ExB_geographic, &
                  eldyn % v_ExB_apex, &
                  eldyn % hall_conductivity, &
                  eldyn % pedersen_conductivity, &
                  eldyn % b_parallel_conductivity, &
                  eldyn % neutral_apex_velocity )
 
  END SUBROUTINE Trash_IPE_Electrodynamics


END MODULE IPE_Electrodynamics_Class
