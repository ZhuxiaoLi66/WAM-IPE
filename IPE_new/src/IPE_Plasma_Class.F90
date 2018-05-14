MODULE IPE_Plasma_Class

USE IPE_Precision
USE IPE_Grid_Class
USE IPE_Neutrals_Class

IMPLICIT NONE

  TYPE IPE_Plasma
    INTEGER :: nFluxTube, NLP, NMP

    REAL(prec), ALLOCATABLE :: ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperatures(:,:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density(:,:,:)
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
REAL(prec), PARAMETER          :: DTMIN_FLIP = 10.0_prec


CONTAINS

  SUBROUTINE Build_IPE_Plasma( plasma, nFluxTube, NLP, NMP )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(out) :: plasma
    INTEGER, INTENT(in)              :: nFluxTube
    INTEGER, INTENT(in)              :: NLP
    INTEGER, INTENT(in)              :: NMP

      plasma % nFluxTube = nFluxTube
      plasma % NLP       = NLP
      plasma % NMP       = NMP

      ALLOCATE( plasma % ion_densities(1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_velocities(1:3,1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_temperatures(1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_density(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_velocity(1:3,1:nFluxTube,1:NLP,1:NMP), &
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

    
    DEALLOCATE( plasma % ion_densities, &
                plasma % ion_velocities, &
                plasma % ion_temperatures, &
                plasma % electron_density, &
                plasma % electron_velocity, &
                plasma % electron_temperature, &
                plasma % hall_conductivity, & 
                plasma % pedersen_conductivity, & 
                plasma % b_parallel_conductivity )

  END SUBROUTINE Trash_IPE_Plasma
!
  SUBROUTINE Update_IPE_Plasma( plasma, grid, neutrals, utime, flip_time_step )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    REAL(prec), INTENT(in)             :: utime
    REAL(prec), INTENT(in)             :: flip_time_step
    ! Local
    INTEGER :: i, lp, mp
    REAL(prec) :: ZX(1:grid % nFluxTube) 
    REAL(prec) :: SLX(1:grid % nFluxTube) 
    REAL(prec) :: GLX(1:grid % nFluxTube) 
    REAL(prec) :: BMX(1:grid % nFluxTube) 
    REAL(prec) :: GRX(1:grid % nFluxTube) 
    REAL(prec) :: OX(1:grid % nFluxTube) 
    REAL(prec) :: HX(1:grid % nFluxTube) 
    REAL(prec) :: N2X(1:grid % nFluxTube) 
    REAL(prec) :: O2X(1:grid % nFluxTube) 
    REAL(prec) :: HEX(1:grid % nFluxTube) 
    REAL(prec) :: N4SX(1:grid % nFluxTube) 
    REAL(prec) :: TNX(1:grid % nFluxTube) 
    REAL(prec) :: TINFX(1:grid % nFluxTube) 
    REAL(prec) :: UNX(1:grid % nFluxTube) 
    

  
      DO mp = 1, plasma % NMP    
        DO lp = 1, plasma % NLP    

          ! Copy over the grid information (for now)
          ZX(1:grid % flux_tube_max(lp))  = grid % altitude(1:grid % flux_tube_max(lp),lp)/1000.0_prec !convert from m to km
          !PCO = Pvalue(lp)  !Pvalue is a single value !*** Need PValue attribute for the grid
          SLX(1:grid % flux_tube_max(lp)) = grid % foot_point_distance(1:grid % flux_tube_max(lp),lp,mp)
          GLX(1:grid % flux_tube_max(lp)) = 0.5_prec*pi - grid % magnetic_colatitude(1:grid % flux_tube_max(lp),lp)  ! magnetic latitude [radians]
          BMX(1:grid % flux_tube_max(lp)) = grid % magnetic_field_strength(1:grid % flux_tube_max(lp),lp,mp)   !Tesla
          GRX(1:grid % flux_tube_max(lp)) = grid % grx(1:grid % flux_tube_max(lp),lp,mp)

          ! Copy over neutrals 
          OX(1:grid % flux_tube_max(lp))    = neutrals % oxygen(1:grid % flux_tube_max(lp),lp,mp) !(m-3)
          HX(1:grid % flux_tube_max(lp))    = neutrals % helium(1:grid % flux_tube_max(lp),lp,mp)
          N2X(1:grid % flux_tube_max(lp))   = neutrals % molecular_nitrogen(1:grid % flux_tube_max(lp),lp,mp)
          O2X(1:grid % flux_tube_max(lp))   = neutrals % molecular_oxygen(1:grid % flux_tube_max(lp),lp,mp)
          HEX(1:grid % flux_tube_max(lp))   = neutrals % helium(1:grid % flux_tube_max(lp),lp,mp)
          N4SX(1:grid % flux_tube_max(lp))  = neutrals % nitrogen(1:grid % flux_tube_max(lp),lp,mp)
          TNX(1:grid % flux_tube_max(lp))   = neutrals % temperature(1:grid % flux_tube_max(lp),lp,mp)
          TINFX(1:grid % flux_tube_max(lp)) = neutrals % temperature_inf(1:grid % flux_tube_max(lp),lp,mp)
          UNX(1:grid % flux_tube_max(lp))   = -neutrals % velocity_apex(3,1:grid % flux_tube_max(lp),lp,mp)


 !>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<!

      ENDDO
    ENDDO


  END SUBROUTINE Update_IPE_Plasma



END MODULE IPE_Plasma_Class
