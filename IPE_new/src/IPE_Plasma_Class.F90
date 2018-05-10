MODULE IPE_Plasma_Class

USE IPE_Model_Precision

IMPLICIT NONE

  TYPE IPE_Plasma
    INTEGER :: nFluxTube, NLP, NMP

    REAL(prec), ALLOCATABLE :: ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperatures(:,:,:,:)

    REAL(prec), ALLOCATABLE :: electron_densitiy(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: hall_conductivity(:,:,:) 
    REAL(prec), ALLOCATABLE :: pedersen_conductivity(:,:,:) 
    REAL(prec), ALLOCATABLE :: b_parallel_conductivity(:,:,:) 

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Plasma
      PROCEDURE :: Trash => Trash_IPE_Plasma

  END TYPE IPE_Plasma


INTEGER, PARAMETER, PRIVATE    :: n_ion_species = 9
REAL(prec), PARAMETER, PRIVATE :: safe_density_minimum = 10.0_prec**(-4)
REAL(prec), PARAMETER, PRIVATE :: safe_temperature_minimum = 100.0_prec


CONTAINS

  SUBROUTINE Build_IPE_Plasma( plasma, nFluxTube, NLP, NMP )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(in) :: plasma
    INTEGER, INTENT(in)             :: nFluxTube
    INTEGER, INTENT(in)             :: NLP
    INTEGER, INTENT(in)             :: NMP

      plasma % nFluxTube = nFluxTube
      plasma % NLP       = NLP
      plasma % NMP       = NMP

      ALLOCATE( plasma % ion_densities(1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_velocities(1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_temperatures(1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_density(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_velocity(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_temperature(1:nFluxTube,1:NLP,1:NMP), &
                plasma % hall_conductivity(1:nFluxTube,1:NLP,1:NMP), &
                plasma % pedersen_conductivity(1:nFluxTube,1:NLP,1:NMP), &
                plasma % b_parallel_conductivity(1:nFluxTube,1:NLP,1:NMP) )
             
      plasma % ion_densities           = safe_density_minimum 
      plasma % ion_velocities          = 0.0_prec
      plasma % ion_temperatures        = safe_temperature_minimum
      plasma % electron_density        = safe_density_minimum
      plasma % electron_velocity       = 0.0_prec
      plasma % electron_temperature    = safe_temperature_minimum
      plasma % hall_conductivity       = 0.0_prec
      plasma % pedersen_conductivity   = 0.0_prec
      plasma % b_parallel_conductivity = 0.0_prec


  END SUBROUTINE Build_IPE_Plasma
!
  SUBROUTINE Trash_IPE_Plasma( plasma )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma

    
    DEALLOCATE( ion_densities, &
                ion_velocities, &
                ion_temperatures, &
                electron_densitiy, &
                electron_velocity, &
                electron_temperature, &
                hall_conductivity, & 
                pedersen_conductivity, & 
                b_parallel_conductivity )

  END SUBROUTINE Trash_IPE_Plasma



END MODULE IPE_Plasma_Class
