MODULE IPE_Plasma_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Grid_Class
USE IPE_Neutrals_Class
USE IPE_Forcing_Class
USE IPE_Time_Class

IMPLICIT NONE

  TYPE IPE_Plasma
    INTEGER :: nFluxTube, NLP, NMP

    REAL(prec), ALLOCATABLE :: ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: ion_densities_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities_old(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperature_old(:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density_old(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature_old(:,:,:)

    REAL(prec), ALLOCATABLE :: ionization_rates(:,:,:,:)

    ! Interpolated Fields
    REAL(prec), ALLOCATABLE :: geo_ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ion_velocities(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ion_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ionization_rates(:,:,:,:) 
    

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Plasma
      PROCEDURE :: Trash => Trash_IPE_Plasma

      PROCEDURE :: Update => Update_IPE_Plasma
      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Plasma
      PROCEDURE :: Read_Legacy_Input => Read_Legacy_Input_IPE_Plasma

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: Buffer_Old_State
      PROCEDURE, PRIVATE :: Calculate_Pole_Values 
      PROCEDURE, PRIVATE :: Cross_Flux_Tube_Transport
      PROCEDURE, PRIVATE :: Auroral_Precipitation
      PROCEDURE, PRIVATE :: FLIP_Wrapper     

  END TYPE IPE_Plasma

  INTEGER, PARAMETER, PRIVATE :: NMP_reduce_factor = 1
  INTEGER, PARAMETER, PRIVATE :: NLP_reduce_factor = 1
  INTEGER, PARAMETER , PRIVATE   :: perp_transport_max_lp = 151
  INTEGER, PARAMETER, PRIVATE    :: n_ion_species = 9
  REAL(prec), PARAMETER, PRIVATE :: safe_density_minimum = 10.0_prec**(-4)
  REAL(prec), PARAMETER, PRIVATE :: safe_temperature_minimum = 100.0_prec
  
  ! Flip Parameters
  REAL(dp), PARAMETER, PRIVATE :: DTMIN    = 1.0D1
  REAL(dp), PARAMETER, PRIVATE :: FPAS     = 0.0D0 
  REAL(dp), PARAMETER, PRIVATE :: HEPRAT   = 9.0D-2 
  REAL(dp), PARAMETER, PRIVATE :: COLFACX  = 1.7D0 
  REAL(dp), PARAMETER, PRIVATE :: HPEQ     = 0.0D0 
  ! IHEPLS,INPLS turn on diffusive solutions if > 0. no solution if 0, chemical equilibrium if < 0
  INTEGER, PARAMETER, PRIVATE  :: IHEPLS   = 1 
  INTEGER, PARAMETER, PRIVATE  :: INPLS    = 1 
  INTEGER, PARAMETER, PRIVATE  :: INNO     = -1
  integer :: istop


! Temporary Timers
  REAL(prec) :: flip_time_avg, aur_precip_time_avg, transport_time_avg, buffer_time_avg
  REAL(prec) :: flip_count, aur_precip_count, transport_count, buffer_count

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
                plasma % ion_temperature(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_density(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_velocity(1:3,1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_temperature(1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_densities_old(1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_velocities_old(1:3,1:n_ion_species,1:nFluxTube,1:NLP,1:NMP), &
                plasma % ion_temperature_old(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_density_old(1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_velocity_old(1:3,1:nFluxTube,1:NLP,1:NMP), &
                plasma % electron_temperature_old(1:nFluxTube,1:NLP,1:NMP), &
                plasma % ionization_rates(1:4,1:nFluxTube,1:NLP,1:NMP) )
             
      plasma % ion_densities           = safe_density_minimum 
      plasma % ion_velocities          = 0.0_prec
      plasma % ion_temperature         = safe_temperature_minimum
      plasma % electron_density        = safe_density_minimum
      plasma % electron_velocity       = 0.0_prec
      plasma % electron_temperature    = safe_temperature_minimum
      plasma % ionization_rates        = 0.0_prec

      ALLOCATE( plasma % geo_ion_densities(1:n_ion_species,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ion_velocities(1:3,1:n_ion_species,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ion_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_density(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_velocity(1:3,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ionization_rates(1:4,1:nlon_geo,1:nlat_geo,1:nheights_geo) )

      plasma % geo_ion_densities           = 0.0_prec
      plasma % geo_ion_velocities          = 0.0_prec
      plasma % geo_ion_temperature         = 0.0_prec
      plasma % geo_electron_density        = 0.0_prec
      plasma % geo_electron_temperature    = 0.0_prec
      plasma % geo_electron_velocity       = 0.0_prec
      plasma % geo_ionization_rates        = 0.0_prec


      flip_time_avg = 0.0_prec
      aur_precip_time_avg = 0.0_prec
      transport_time_avg = 0.0_prec
      buffer_time_avg = 0.0_prec
      flip_count = 0.0_prec
      aur_precip_count = 0.0_prec
      transport_count = 0.0_prec
      buffer_count = 0.0_prec


  END SUBROUTINE Build_IPE_Plasma
!
  SUBROUTINE Trash_IPE_Plasma( plasma )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma

    
    DEALLOCATE( plasma % ion_densities, &
                plasma % ion_velocities, &
                plasma % ion_temperature, &
                plasma % electron_density, &
                plasma % electron_velocity, &
                plasma % electron_temperature, &
                plasma % ion_densities_old, &
                plasma % ion_velocities_old, &
                plasma % ion_temperature_old, &
                plasma % electron_density_old, &
                plasma % electron_velocity_old, &
                plasma % electron_temperature_old, &
                plasma % ionization_rates, &
                plasma % geo_ion_velocities, &
                plasma % geo_ion_temperature, &
                plasma % geo_electron_density, &
                plasma % geo_electron_velocity, &
                plasma % geo_electron_temperature, &
                plasma % geo_ionization_rates )

     PRINT*, 'Buffer ', buffer_time_avg/buffer_count
     PRINT*, 'Transport ', transport_time_avg/transport_count
     PRINT*, 'Auroral_Precipitation ', aur_precip_time_avg/aur_precip_count
     PRINT*, 'FLIP_Wrapper ', flip_time_avg/flip_count

  END SUBROUTINE Trash_IPE_Plasma
!
  SUBROUTINE Update_IPE_Plasma( plasma, grid, neutrals, forcing, time_tracker, v_ExB, time_step )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP,1:grid % NMP)
    REAL(prec), INTENT(in)             :: time_step
    ! Local
    INTEGER    :: time_loop, n_perp_transport_steps
    INTEGER    :: i, lp, mp, j
    REAL(prec) :: t2, t1

         
      CALL CPU_TIME( t1 )
      CALL plasma % Buffer_Old_State( )
      CALL CPU_TIME( t2 )
      buffer_time_avg = buffer_time_avg + t2 - t1
      buffer_count    = buffer_count + 1


      CALL CPU_TIME( t1 )
      CALL plasma % Cross_Flux_Tube_Transport( grid, v_ExB, time_step )
      CALL CPU_TIME( t2 )
      transport_time_avg = transport_time_avg + t2 - t1
      transport_count    = transport_count + 1


 
      CALL CPU_TIME( t1 )
      CALL plasma % Auroral_Precipitation( grid, &
                                           neutrals, &
                                           forcing, &
                                           time_tracker )
      CALL CPU_TIME( t2 )
      aur_precip_time_avg = aur_precip_time_avg + t2 - t1
      aur_precip_count    = aur_precip_count + 1

      CALL CPU_TIME( t1 )
      !CALL plasma % FLIP_Wrapper( grid, & 
      !                            neutrals, &
      !                            forcing, &
      !                            time_tracker, &
      !                            time_step )
      CALL CPU_TIME( t2 )
      flip_time_avg = flip_time_avg + t2 - t1
      flip_count    = flip_count + 1


  END SUBROUTINE Update_IPE_Plasma
!
  SUBROUTINE Buffer_Old_State( plasma )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma

      plasma % ion_densities_old   = plasma % ion_densities
      plasma % ion_temperature_old = plasma % ion_temperature
      plasma % ion_velocities_old  = plasma % ion_velocities

      plasma % electron_density_old   = plasma % electron_density
      plasma % electron_temperature_old = plasma % electron_temperature

  END SUBROUTINE Buffer_Old_State
!
  SUBROUTINE Calculate_Pole_Values( plasma, grid, ion_densities_pole_value, ion_temperature_pole_value, ion_velocities_pole_value, electron_density_pole_value, electron_temperature_pole_value ) 
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(in) :: plasma
    TYPE( IPE_Grid ), INTENT(in)    :: grid
    REAL(prec), INTENT(out)         :: ion_densities_pole_value(1:n_ion_species, 1:plasma % nFluxTube)
    REAL(prec), INTENT(out)         :: ion_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec), INTENT(out)         :: ion_velocities_pole_value(1:3,1:n_ion_species, 1:plasma % nFluxTube)
    REAL(prec), INTENT(out)         :: electron_density_pole_value(1:plasma % nFluxTube)
    REAL(prec), INTENT(out)         :: electron_temperature_pole_value(1:plasma % nFluxTube)
    ! Local
    INTEGER    :: i, lp, mp, j

      ion_densities_pole_value      = 0.0_prec
      ion_temperature_pole_value    = 0.0_prec
      ion_velocities_pole_value     = 0.0_prec

      electron_density_pole_value     = 0.0_prec
      electron_temperature_pole_value = 0.0_prec

      DO mp = 1, plasma % NMP , NMP_reduce_factor
        DO i = 1, grid % flux_tube_max(1)
          DO j = 1, n_ion_species
        
            ion_densities_pole_value(j,i)      = ion_densities_pole_value(j,i) + plasma % ion_densities(j,i,1,mp)
            ion_velocities_pole_value(1:3,j,i) = ion_velocities_pole_value(1:3,j,i) + plasma % ion_velocities(1:3,j,i,1,mp)

          ENDDO

          ion_temperature_pole_value(i)       = ion_temperature_pole_value(i) + plasma % ion_temperature(i,1,mp)
          electron_density_pole_value(i)      = electron_density_pole_value(i) + plasma % electron_density(i,1,mp)
          electron_temperature_pole_value(i)  = electron_temperature_pole_value(i) + plasma % electron_temperature(i,1,mp)

        ENDDO
      ENDDO

      ion_densities_pole_value        = ion_densities_pole_value/REAL( plasma % NMP )
      ion_temperature_pole_value      = ion_temperature_pole_value/REAL( plasma % NMP )
      ion_velocities_pole_value       = ion_velocities_pole_value/REAL( plasma % NMP )
      electron_density_pole_value     = electron_density_pole_value/REAL( plasma % NMP )
      electron_temperature_pole_value = electron_temperature_pole_value/REAL( plasma % NMP )

  END SUBROUTINE Calculate_Pole_Values
      
  SUBROUTINE Cross_Flux_Tube_Transport( plasma, grid, v_ExB, time_step )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP, 1:grid % NMP)
    REAL(prec), INTENT(in)             :: time_step

    ! Local
    REAL(prec) :: ion_densities_pole_value(1:n_ion_species, 1:plasma % nFluxTube)
    REAL(prec) :: ion_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: ion_velocities_pole_value(1:3,1:n_ion_species, 1:plasma % nFluxTube)
    REAL(prec) :: electron_density_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: electron_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: colat_90km(1:grid % NLP)
    REAL(prec) :: phi_t0 !magnetic longitude,phi[rad] at t0(previous time step)
    REAL(prec) :: theta_t0 !magnetic latitude,theta[rad] at t0
    REAL(prec) :: q_value
    REAL(prec) :: phi_i(1:2)
    REAL(prec) :: lp_comp_weight(1:2)
    REAL(prec) :: mp_comp_weight(1:2)
    REAL(prec) :: i_comp_weight(1:2)
    REAL(prec) :: q_int(1:2)
    REAL(prec) :: ion_densities_int(1:n_ion_species)
    REAL(prec) :: ion_temperature_int
    REAL(prec) :: ion_velocities_int(1:3,1:n_ion_species)
    REAL(prec) :: B_int, max_phi, ksi_fac
    INTEGER    :: lp_t0(1:2)
    INTEGER    :: mp_t0(1:2)
    INTEGER    :: lp_min, mp_min, isouth, inorth, ii, ispecial
    INTEGER    :: mp, lp, i, lpx, mpx, jth
    INTEGER    :: i_min(1:2)

      CALL plasma % Calculate_Pole_Values( grid, &
                                           ion_densities_pole_value,  &
                                           ion_temperature_pole_value, &
                                           ion_velocities_pole_value, &
                                           electron_density_pole_value, &
                                           electron_temperature_pole_value ) 

      colat_90km(1:grid % NLP) = grid % magnetic_colatitude(1,1:grid % NLP)

      max_phi = MAXVAL( grid % magnetic_longitude )

      DO mp = 1, grid % NMP , NMP_reduce_factor
        DO lp = 1, perp_transport_max_lp
  
          theta_t0 = colat_90km(lp) - v_ExB(1,lp,mp)*time_step
          phi_t0   = grid % magnetic_longitude(mp) - v_ExB(2,lp,mp)*time_step

          ! If a Lagrangian trajectory crosses the equator, we clip the colatitude
          ! so that the point resides at the equator.
          IF( theta_t0 > colat_90km( grid % NLP ) )THEN ! NLP ==> Equator
            theta_t0 = colat_90km( grid % NLP )
          ENDIF
  
          ! Apply periodicity in the zonal direction
          IF ( phi_t0 >= max_phi ) THEN
  
            phi_t0 = phi_t0 - pi*2.0_prec
  
          ELSE IF ( phi_t0 < 0.0_prec    ) THEN
  
            phi_t0 = phi_t0 + pi*2.0_prec
  
          ENDIF
  
          ! We are trying to find the lp indices for the flux tubes nearest to
          ! theta_t0. Output is lp_t0(1:2)

          lp_min = grid % NLP

          DO lpx = 1, grid % NLP

            IF( theta_t0 <= colat_90km(lpx) )THEN
              lp_min = lpx
              EXIT
            ENDIF

          ENDDO

          lp_t0(1) = lp_min-1
          lp_t0(2) = lp_min

          mp_min = grid % NMP

          DO mpx = 1, grid % NMP , NMP_reduce_factor

            IF( phi_t0 <= grid % magnetic_longitude(mpx) )THEN
              mp_min = mpx
              EXIT
            ENDIF

          ENDDO

          IF( mp_min == 1 )THEN

            mp_t0(1) = grid % NMP
            mp_t0(2) = 1

            phi_i(1) = grid % magnetic_longitude(mp_t0(1)) - 2.0_prec*pi
            phi_i(2) = grid % magnetic_longitude(mp_t0(2))

          ELSE

            mp_t0(1) = mp_min-1
            mp_t0(2) = mp_min

            phi_i(1) = grid % magnetic_longitude(mp_t0(1))
            phi_i(2) = grid % magnetic_longitude(mp_t0(2))
             
          ENDIF

          mp_comp_weight(1) =  ( phi_t0 - phi_i(2) )/( phi_i(1)-phi_i(2) )
          mp_comp_weight(2) = -( phi_t0 - phi_i(1) )/( phi_i(1)-phi_i(2) )

          IF( lp_min == 1 )THEN 

            ! For this lp, mp, we need to loop over the flux tube points. 
            ! For each flux tube point, we need to obtain the "q-value"
            ! On the four points nearest the source point (phi_t0, theta_t0) 
            ! we need to find the fluxe tube points nearest to the q-value
            ! to do along flux tube interpolation
            DO i = 1, grid % flux_tube_max(lp)
  
              q_value = grid % q_factor(i,lp,mp)
  
              ion_densities_int   = 0.0_prec
              ion_temperature_int = 0.0_prec
              ion_velocities_int  = 0.0_prec

              ispecial = 1
              isouth   = grid % flux_tube_max(1)
              inorth   = isouth-1

              DO ii = 1, grid % flux_tube_max(1)
                IF(  grid % q_factor(ii, 1, mp) < q_value )THEN
                  isouth   = ii
                  inorth   = ii-1
                  ispecial = 0
                  EXIT
                ENDIF
              ENDDO

              IF( isouth == 1 )THEN
                ispecial = 2
              ENDIF

              IF( ispecial == 0 )THEN

                q_int(1) = grid % q_factor(isouth, 1, mp)
                q_int(2) = grid % q_factor(inorth, 1, mp)
           
                i_comp_weight(1) = ( q_value - q_int(2) )/( q_int(1) - q_int(2) )
                i_comp_weight(2) = -( q_value - q_int(1) )/( q_int(1) - q_int(2) )

 
                ion_densities_int(1:n_ion_species) = ion_densities_int(1:n_ion_species) + &
                                      ( ion_densities_pole_value(1:n_ion_species, isouth)*i_comp_weight(1) + &
                                        ion_densities_pole_value(1:n_ion_species, inorth)*i_comp_weight(2) )

                ion_temperature_int = ion_temperature_int + &
                                      ( ion_temperature_pole_value(isouth)*i_comp_weight(1) + &
                                        ion_temperature_pole_value(inorth)*i_comp_weight(2) )

                ion_velocities_int(1:3,1:n_ion_species) = ion_velocities_int(1:3,1:n_ion_species) + &
                                      ( ion_velocities_pole_value(1:3,1:n_ion_species, isouth)*i_comp_weight(1) + &
                                        ion_velocities_pole_value(1:3,1:n_ion_species, inorth)*i_comp_weight(2) )

              ELSE

                ion_densities_int(1:n_ion_species) = ion_densities_int(1:n_ion_species) + ion_densities_pole_value(1:n_ion_species, isouth)

                ion_temperature_int = ion_temperature_int + ion_temperature_pole_value(isouth)

                ion_velocities_int(1:3,1:n_ion_species) = ion_velocities_int(1:3,1:n_ion_species) + ion_velocities_pole_value(1:3,1:n_ion_species, isouth)

              ENDIF

              plasma % ion_densities(1:n_ion_species,i,lp,mp) = ion_densities_int(1:n_ion_species)
              plasma % ion_temperature(i,lp,mp) = ion_temperature_int
              plasma % ion_velocities(1:3,1:n_ion_species,i,lp,mp) = ion_velocities_int(1:3,1:n_ion_species)

            ENDDO

  
          ELSEIF( lp_min == grid % NLP )THEN
  
 

            lp_comp_weight(1) =  ( theta_t0 - colat_90km(lp_t0(2)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )
            lp_comp_weight(2) = -( theta_t0 - colat_90km(lp_t0(1)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )
            ! In this case, information is propagating from beyond the NLPth
            ! flux tube, which just around the magnetic equator. In this case
            ! use the last solution in the lp direction for the lp interpolation

            DO i = 1, grid % flux_tube_max(lp)
  
              q_value = grid % q_factor(i,lp,mp)
  
              ion_densities_int   = 0.0_prec
              ion_temperature_int = 0.0_prec
              ion_velocities_int  = 0.0_prec
              B_int               = 0.0_prec
  
              DO mpx = 1, 2

               ispecial = 1
               isouth   = grid % flux_tube_max(grid % NLP)
               inorth   = isouth-1

               DO ii = 1, grid % flux_tube_max(grid % NLP)
                 IF(  grid % q_factor(ii, grid % NLP, mp_t0(mpx)) < q_value )THEN
                   isouth   = ii
                   inorth   = ii-1
                   ispecial = 0
                   EXIT
                 ENDIF
               ENDDO

               IF( isouth == 1 )THEN
                 ispecial = 2
               ENDIF

               IF( ispecial == 0 )THEN

                 q_int(1) = grid % q_factor(isouth, grid % NLP, mp_t0(mpx))
                 q_int(2) = grid % q_factor(inorth, grid % NLP, mp_t0(mpx))
             
                 i_comp_weight(1) = ( q_value - q_int(2) )/( q_int(1) - q_int(2) )
                 i_comp_weight(2) = -( q_value - q_int(1) )/( q_int(1) - q_int(2) )
  
                 B_int = B_int + ( grid % magnetic_field_strength(isouth, grid % NLP, mp_t0(mpx))*i_comp_weight(1) + &
                                   grid % magnetic_field_strength(inorth, grid % NLP, mp_t0(mpx))*i_comp_weight(2) )*&
                                 mp_comp_weight(mpx)

                 ion_densities_int(1:n_ion_species) = ion_densities_int(1:n_ion_species) + &
                                       ( plasma % ion_densities_old(1:n_ion_species, isouth, grid % NLP, mp_t0(mpx))*i_comp_weight(1) + &
                                         plasma % ion_densities_old(1:n_ion_species, inorth, grid % NLP, mp_t0(mpx))*i_comp_weight(2) )*&
                                       mp_comp_weight(mpx)

                 ion_temperature_int = ion_temperature_int + &
                                       ( plasma % ion_temperature_old(isouth, grid % NLP, mp_t0(mpx))*i_comp_weight(1) + &
                                         plasma % ion_temperature_old(inorth, grid % NLP, mp_t0(mpx))*i_comp_weight(2) )*&
                                       mp_comp_weight(mpx)

                 ion_velocities_int(1:3,1:n_ion_species) = ion_velocities_int(1:3,1:n_ion_species) + &
                                       ( plasma % ion_velocities_old(1:3,1:n_ion_species, isouth, grid % NLP, mp_t0(mpx))*i_comp_weight(1) + &
                                         plasma % ion_velocities_old(1:3,1:n_ion_species, inorth, grid % NLP, mp_t0(mpx))*i_comp_weight(2) )*&
                                       mp_comp_weight(mpx)
 
               ELSE

                 B_int = B_int + grid % magnetic_field_strength(isouth, grid % NLP, mp_t0(mpx))*mp_comp_weight(mpx)

                 ion_densities_int(1:n_ion_species) = ion_densities_int(1:n_ion_species) + &
                                                      plasma % ion_densities_old(1:n_ion_species, isouth, grid % NLP, mp_t0(mpx))*&
                                                      mp_comp_weight(mpx)

                 ion_temperature_int = ion_temperature_int + &
                                       plasma % ion_temperature_old(isouth, grid % NLP, mp_t0(mpx))*&
                                       mp_comp_weight(mpx)

                 ion_velocities_int(1:3,1:n_ion_species) = ion_velocities_int(1:3,1:n_ion_species) + &
                                                           plasma % ion_velocities_old(1:3,1:n_ion_species, isouth, grid % NLP, mp_t0(mpx))*&
                                                           mp_comp_weight(mpx)

               ENDIF
  
              ENDDO

              IF( lp <= 17 )THEN
                ksi_fac = 1.0_prec
              ELSE
                ksi_fac = grid % magnetic_field_strength(i,lp,mp)/B_int
              ENDIF

              plasma % ion_densities(1:n_ion_species,i,lp,mp) = ion_densities_int(1:n_ion_species)*( ksi_fac )**2
              plasma % ion_temperature(i,lp,mp) = ion_temperature_int*( ksi_fac**(4.0_prec/3.0_prec) )
              plasma % ion_velocities(1:3,1:n_ion_species,i,lp,mp) = ion_velocities_int(1:3,1:n_ion_species)

            ENDDO
  
          ELSE
      
            lp_comp_weight(1) =  ( theta_t0 - colat_90km(lp_t0(2)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )
            lp_comp_weight(2) = -( theta_t0 - colat_90km(lp_t0(1)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )

            DO i = 1, grid % flux_tube_max(lp)
  
              q_value = grid % q_factor(i,lp,mp)
  
              ion_densities_int   = 0.0_prec
              ion_temperature_int = 0.0_prec
              ion_velocities_int  = 0.0_prec
              B_int               = 0.0_prec
  
              DO mpx = 1, 2
                DO lpx = 1, 2
 
                  ispecial = 1
                  isouth   = grid % flux_tube_max(lp_t0(lpx))
                  inorth   = isouth-1

                  DO ii = 1, grid % flux_tube_max(lp_t0(lpx))
                    IF(  grid % q_factor(ii, lp_t0(lpx), mp_t0(mpx)) < q_value )THEN
                      isouth   = ii
                      inorth   = ii-1
                      ispecial = 0
                      EXIT
                    ENDIF
                  ENDDO

                  IF( isouth == 1 )THEN
                    ispecial = 2
                  ENDIF

                  IF( ispecial == 0 )THEN

                    q_int(1) = grid % q_factor(isouth, lp_t0(lpx), mp_t0(mpx))
                    q_int(2) = grid % q_factor(inorth, lp_t0(lpx), mp_t0(mpx))
             
                    i_comp_weight(1) = ( q_value - q_int(2) )/( q_int(1) - q_int(2) )
                    i_comp_weight(2) = -( q_value - q_int(1) )/( q_int(1) - q_int(2) )

                    B_int = B_int + ( grid % magnetic_field_strength(isouth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(1) + &
                                     grid % magnetic_field_strength(inorth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(2) )*&
                                   lp_comp_weight(lpx)*mp_comp_weight(mpx)
  
                    ion_densities_int(1:n_ion_species) = ion_densities_int(1:n_ion_species) + &
                                          ( plasma % ion_densities_old(1:n_ion_species, isouth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(1) + &
                                            plasma % ion_densities_old(1:n_ion_species, inorth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(2) )*&
                                          lp_comp_weight(lpx)*mp_comp_weight(mpx)

                    ion_temperature_int = ion_temperature_int + &
                                          ( plasma % ion_temperature_old(isouth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(1) + &
                                            plasma % ion_temperature_old(inorth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(2) )*&
                                          lp_comp_weight(lpx)*mp_comp_weight(mpx)

                    ion_velocities_int(1:3,1:n_ion_species) = ion_velocities_int(1:3,1:n_ion_species) + &
                                          ( plasma % ion_velocities_old(1:3,1:n_ion_species, isouth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(1) + &
                                            plasma % ion_velocities_old(1:3,1:n_ion_species, inorth, lp_t0(lpx), mp_t0(mpx))*i_comp_weight(2) )*&
                                          lp_comp_weight(lpx)*mp_comp_weight(mpx)
 
                  ELSE

                    B_int = B_int + grid % magnetic_field_strength(isouth, lp_t0(lpx), mp_t0(mpx))*&
                                    lp_comp_weight(lpx)*mp_comp_weight(mpx)

                    ion_densities_int(1:n_ion_species) = ion_densities_int(1:n_ion_species) + &
                                                         plasma % ion_densities_old(1:n_ion_species, isouth, lp_t0(lpx), mp_t0(mpx))*&
                                                         lp_comp_weight(lpx)*mp_comp_weight(mpx)

                    ion_temperature_int = ion_temperature_int + &
                                          plasma % ion_temperature_old(isouth, lp_t0(lpx), mp_t0(mpx))*&
                                          lp_comp_weight(lpx)*mp_comp_weight(mpx)

                    ion_velocities_int(1:3,1:n_ion_species) = ion_velocities_int(1:3,1:n_ion_species) + &
                                                              plasma % ion_velocities_old(1:3,1:n_ion_species, isouth, lp_t0(lpx), mp_t0(mpx))*&
                                                              lp_comp_weight(lpx)*mp_comp_weight(mpx)

                  ENDIF
  
                ENDDO
              ENDDO

              IF( lp <= 17 )THEN
                ksi_fac = 1.0_prec
              ELSE
                ksi_fac = grid % magnetic_field_strength(i,lp,mp)/B_int
              ENDIF

              plasma % ion_densities(1:n_ion_species,i,lp,mp) = ion_densities_int(1:n_ion_species)*( ksi_fac**2 )
              plasma % ion_temperature(i,lp,mp) = ion_temperature_int*( ksi_fac**(4.0_prec/3.0_prec) )
              plasma % ion_velocities(1:3,1:n_ion_species,i,lp,mp) = ion_velocities_int(1:3,1:n_ion_species)

            ENDDO
  

          ENDIF  
 
        ENDDO
      ENDDO
     
  END SUBROUTINE Cross_Flux_Tube_Transport

  SUBROUTINE Auroral_Precipitation( plasma, grid, neutrals, forcing, time_tracker )
  ! Previously : tiros_ionize_ipe 
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    ! Local
    INTEGER, PARAMETER :: jmaxwell = 6
    INTEGER    :: i, lp, mp, m, j, iband, l
    INTEGER    :: i1, i2, j1, j2, k, kk
    INTEGER    :: tiros_activity_level 
    REAL(prec) :: ratio(1:21)
    REAL(prec) :: bz, gscon, amu, e0, mlt, dfac
    REAL(prec) :: ch, chi, diff, mo, mo2, mn2, alpha 
    REAL(prec) :: eflux, essa, dl_lower,dl_upper,qiont_lower,qiont_upper 
    REAL(prec) :: ri, rj, th, THMagd, GW
    REAL(prec) :: pres(plasma % nFluxTube)
    REAL(prec) :: den(plasma % nFluxTube), grav(plasma % nFluxTube)
    REAL(prec) :: ntot(plasma % nFluxTube), meanmass(plasma % nFluxTube)
    REAL(prec) :: rno, rang, pr, ratioz, rlamz, mh, mhe, q
    REAL(prec) :: TE11(21),TE15(21),width_maxwell
    REAL(prec) :: ratio_ch,en_maxwell(jmaxwell),dl(jmaxwell)
    REAL(prec) :: qion_maxwell(plasma % nFluxTube)
    REAL(prec) :: en(1:15), width(1:15), RLAM(1:21), ionchr(1:21)

       en(1:15) = (/ 0.37_prec, 0.6_prec, 0.92_prec, 1.37_prec, 2.01_prec, &
                     2.91_prec, 4.19_prec, 6.0_prec, 8.56_prec, 12.18_prec, &
                     17.3_prec, 24.49_prec, 36.66_prec, 54.77_prec, 81.82_prec /)

       width(1:15) = (/ 0.158_prec, 0.315_prec, 0.315_prec, 0.63_prec, 0.631_prec, &
                        1.261_prec, 1.26_prec, 2.522_prec, 2.522_prec, 5.043_prec, &
                        5.043_prec, 10.0_prec, 14.81_prec, 22.13_prec, 33.06_prec /)

       rlam(1:21) = (/ 1.49_prec, 1.52_prec, 1.51_prec, 1.48_prec, 1.43_prec, &
                                  1.37_prec, 1.30_prec, 1.22_prec, 1.12_prec, 1.01_prec, &
                                  0.895_prec, 0.785_prec, 0.650_prec, 0.540_prec, 0.415_prec, &
                                  0.320_prec, 0.225_prec, 0.14_prec, 0.08_prec, 0.04_prec, 0.0_prec /)

       ionchr(1:21) = (/ 0.378_prec, 0.458_prec, 0.616_prec, 0.773_prec, 0.913_prec, &
                                     1.088_prec, 1.403_prec, 1.718_prec, 2.033_prec, 2.349_prec, &
                                     2.979_prec, 3.610_prec, 4.250_prec, 4.780_prec, 6.130_prec, &
                                     7.392_prec, 8.653_prec, 9.914_prec, 12.436_prec, 14.957_prec, 17.479_prec /)

      tiros_activity_level = forcing % nhemi_power_index( forcing % current_index )
      GW                   = forcing % nhemi_power( forcing % current_index )
      bz    = 1.38_prec*10.0_prec**(-23)
      gscon = 8.314_prec*10.0_prec**(3)
      mo    = 16.0_prec
      mo2   = 32.0_prec
      mn2   = 28.0_prec
      mh    = 1.0_prec
      mhe   = 4.0_prec
      amu   = 1.661_prec*10.0_prec**(-27)
      E0    = 0.035_prec
      WIDTH_maxwell = 0.050_prec

      DO m = 1, 21
        ratio(m) = REAL( (m-1), prec )*0.05_prec
      ENDDO

      DO iband=1,21

        te15(iband)=0.0_prec
        te11(iband)=0.0_prec

        DO m = 1, 15
          te15(iband) = te15(iband) + forcing % djspectra(m,iband)*en(m)*width(m)*1.6_prec*10.0_prec**(-6)
        ENDDO
        DO m = 1,11
          te11(iband) = te15(iband) + forcing % djspectra(m,iband)*en(m)*width(m)*1.6_prec*10.0_prec**(-6)
        ENDDO
       
      ENDDO

      DO j = 1, jmaxwell
        en_maxwell(j) = REAL( j, prec )*0.05_prec - 0.025_prec
      ENDDO

      DO mp = 1, grid % NMP , NMP_reduce_factor
        DO lp = 1, grid % NLP , NLP_reduce_factor

          DO i=1,grid % flux_tube_max(lp)
            plasma % ionization_rates(1,i,lp,mp) = 0.0_prec
            plasma % ionization_rates(2,i,lp,mp) = 0.0_prec
            plasma % ionization_rates(3,i,lp,mp) = 0.0_prec
            plasma % ionization_rates(4,i,lp,mp) = 0.0_prec
            qion_maxwell(i)                      = 0.0_prec
          ENDDO

          ! convert magnetic latitude from radians to degrees
          ! Use the geomagnetic latitude of the foot point
          thmagd = ( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )*180.0_prec/pi
          th = abs(thmagd) - 50.0_prec

          IF ( abs(thmagd) > 50.0_prec )  THEN

            mlt = time_tracker % utime/3600.0_prec + &
                  grid % longitude(1,lp,mp)*180.0_prec/pi/15.0_prec

            essa = (mlt + 12.0_prec)*15.0_prec

            IF ( essa >= 360.0_prec ) THEN
              essa = essa - 360.0_prec
            ELSEIF ( essa < 0.0_prec ) THEN
              essa = essa + 360.0_prec
            ENDIF

            dfac = 1.0_prec
            IF( tiros_activity_level > 9 .AND. GW > 96.0_prec ) THEN
              dfac = GW/96.0_prec
            ENDIF

            l = tiros_activity_level - 2
            IF ( l < 1 ) THEN
              l = 1
            ELSEIF ( l > 7 ) THEN
              l = 7
            ENDIF

            ri = essa/18.0_prec + 11.0_prec
            i1 = ri
            ri = ri - i1
            IF ( i1 > 20 ) THEN
              i1 = i1 - 20
            ENDIF

            i2 = i1 + 1
            IF ( i2 > 20 ) THEN
              i2 = i2 - 20
            ENDIF

            rj = th/2.0_prec + 1.0_prec
            j1 = rj
            rj = rj - j1
            j2 = j1 + 1

            eflux = rj*ri*forcing % emaps(j2,i2,l) + &
                    (1.0_prec-rj)*ri*forcing % emaps(j1,i2,l) + &
                    rj*(1.0_prec-ri)*forcing % emaps(j2,i1,l) + &
                    (1.0_prec-rj)*(1.0_prec-ri)*forcing % emaps(j1,i1,l)
            eflux = 10.0_prec**(eflux)/1000.0_prec

            ch = rj*ri*forcing % cmaps(j2,i2,l) + &
                 (1.0_prec-rj)*ri*forcing % cmaps(j1,i2,l) + &
                 rj*(1.0_prec-ri)*forcing % cmaps(j2,i1,l) + &
                 (1.0_prec-rj)*(1.0_prec-ri)*forcing % cmaps(j1,i1,l)


            IF ( ch < 0.378_prec ) THEN
              ch = 0.379_prec
            ENDIF

            DO kk = 2 , 21
              IF ( ch <= ionchr(kk) ) THEN
                k = kk - 1
                EXIT
              ENDIF
            ENDDO

            kk = k+1
            chi = ch - ionchr(k)
            diff = ionchr(kk) - ionchr(k)
            ratio_ch = chi/diff

            DO i = 1 , grid % flux_tube_max(lp)

              IF ( grid % altitude(i,lp) <= 10.0_prec**(6) )THEN

                grav(i) = -grid % grx(i,lp,mp)

                ntot(i) = ( neutrals % oxygen(i,lp,mp)  +&
                            neutrals % molecular_oxygen(i,lp,mp) +&
                            neutrals % molecular_nitrogen(i,lp,mp) +&
                            neutrals % hydrogen(i,lp,mp)  +&
                            neutrals % helium(i,lp,mp) )

                pres(i) = ntot(i)*bz*neutrals % temperature(i,lp,mp)

                meanmass(i) = ( neutrals % oxygen(i,lp,mp)*mo   +&
                                neutrals % molecular_oxygen(i,lp,mp)*mo2 +&
                                neutrals % molecular_nitrogen(i,lp,mp)*mn2 +&
                                neutrals % hydrogen(i,lp,mp)*mh   +&
                                neutrals % helium(i,lp,mp)*mhe )/ntot(i)

                den(i) = pres(i)*meanmass(i)/(gscon*neutrals % temperature(i,lp,mp))

                DO l = 1 , 15

                  dl_lower = forcing % djspectra(l,k)
                  dl_upper = forcing % djspectra(l,kk)
                  rang = 4.57_prec*10.0_prec**(-5)*en(l)**1.75_prec
                  pr = rang*grav(i)
                  ratioz = pres(i)/pr

                  IF ( ratioz > 1.0_prec ) THEN

                    rlamz = 0.0_prec

                  ELSE

                    DO m = 2 , 21 !! Used to loop from 1,21
                      IF ( ratioz <= ratio(m) ) THEN
                        rlamz = rlam(m-1) + (ratioz-ratio(m-1))*&
                                (rlam(m)-rlam(m-1))/(ratio(m)-ratio(m-1))
                        EXIT
                      ENDIF
                    ENDDO

                  ENDIF

                  qiont_lower = den(i)*en(l)*rlamz*dl_lower*width(l)*1.0_prec*10.0_prec**7/rang/e0
                  qiont_upper = den(i)*en(l)*rlamz*dl_upper*width(l)*1.0_prec*10.0_prec**7/rang/e0

                  plasma % ionization_rates(1,i,lp,mp) = plasma % ionization_rates(1,i,lp,mp) + &
                                                         (ratio_ch*qiont_upper+(1.0_prec-ratio_ch)*qiont_lower)*&
                                                         eflux/(ratio_ch*te11(kk)+(1.0_prec-ratio_ch)*te11(k))

                ENDDO

                alpha = ch/2.0_prec
                rno = eflux*6.24_prec*10.0_prec**(12)/2.0_prec/alpha**3

                DO l = 1 , JMAXWELL

                  dl(l) = rno*en_maxwell(l)*EXP(-en_maxwell(l)/alpha)
                  rang = 4.57_prec*10.0_prec**(-5)*en_maxwell(l)**1.75_prec
                  pr = rang*grav(i)
                  ratioz = pres(i)/pr

                  IF ( ratioz > 1.0_prec ) THEN
                    rlamz = 0.0_prec
                  ELSE

                    DO m = 2 , 21 !! Used to loop from 1, 21
                      IF ( ratioz <= ratio(m) ) THEN
                        rlamz = rlam(m-1) + (ratioz-ratio(m-1))*(rlam(m)-rlam(m-1))/(ratio(m)-ratio(m-1))
                        EXIT
                      ENDIF
                    ENDDO

                  ENDIF

                  qion_maxwell(i) = qion_maxwell(i) + den(i)*en_maxwell(l)*rlamz*dl(l)*width_maxwell/rang/e0

                ENDDO

                plasma % ionization_rates(1,i,lp,mp) = ( plasma % ionization_rates(1,i,lp,mp) + qion_maxwell(i) )*dfac

                q = plasma % ionization_rates(1,i,lp,mp)/( 0.92_prec*neutrals % molecular_nitrogen(i,lp,mp) +&
                                                           1.50_prec*neutrals % molecular_oxygen(i,lp,mp)  +&
                                                           0.56_prec*neutrals % oxygen(i,lp,mp))

                plasma % ionization_rates(2,i,lp,mp) = ( 0.50_prec*neutrals % molecular_oxygen(i,lp,mp) +& ! O+ Rate
                                                         0.56_prec*neutrals % oxygen(i,lp,mp) )*q

                plasma % ionization_rates(3,i,lp,mp) = neutrals % molecular_oxygen(i,lp,mp)*q             ! O2+ Rate

                plasma % ionization_rates(4,i,lp,mp) = 0.92_prec*neutrals % molecular_nitrogen(i,lp,mp)*q ! N2+ Rate

              ENDIF

            ENDDO
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE Auroral_Precipitation

  SUBROUTINE FLIP_Wrapper( plasma, grid, neutrals, forcing, time_tracker, flip_time_step )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    REAL(prec), INTENT(in)             :: flip_time_step
    ! Local
    INTEGER  :: i, lp, mp
    INTEGER  :: JMINX, JMAXX
    INTEGER  :: EFLAG(11,11)
    REAL(dp) :: PCO, UTHR
    REAL(dp) :: ZX(1:grid % nFluxTube) 
    REAL(dp) :: SLX(1:grid % nFluxTube) 
    REAL(dp) :: GLX(1:grid % nFluxTube) 
    REAL(dp) :: BMX(1:grid % nFluxTube) 
    REAL(dp) :: GRX(1:grid % nFluxTube) 
    REAL(dp) :: OX(1:grid % nFluxTube) 
    REAL(dp) :: HX(1:grid % nFluxTube) 
    REAL(dp) :: N2X(1:grid % nFluxTube) 
    REAL(dp) :: O2X(1:grid % nFluxTube) 
    REAL(dp) :: HEX(1:grid % nFluxTube) 
    REAL(dp) :: N4SX(1:grid % nFluxTube) 
    REAL(dp) :: TNX(1:grid % nFluxTube) 
    REAL(dp) :: TINFX(1:grid % nFluxTube) 
    REAL(dp) :: UNX(1:grid % nFluxTube) 
    REAL(dp) :: XIONNX(1:9,1:grid % nFluxTube)
    REAL(dp) :: XIONVX(1:9,1:grid % nFluxTube) 
    REAL(dp) :: TE_TIX(1:3,1:grid % nFluxTube)
    REAL(dp) :: EHTX(1:3,1:grid % nFluxTube)
    REAL(dp) :: AUR_PROD(1:3,1:grid % nFluxTube)
    REAL(dp) :: NNOX(1:grid % nFluxTube) 
    REAL(dp) :: NHEAT(1:grid % nFluxTube) 
    REAL(dp) :: SZA(1:grid % nFluxTube) 
    REAL(sp) :: F107D, F107A

      F107D = forcing % f107( forcing % current_index )   
      F107A = forcing % f107_81day_avg( forcing % current_index )   
      UTHR  = time_tracker % hour

      print *, 'GHGM TIME ', UTHR

      print *, 'GHGM ', plasma % NMP , plasma % NLP
      DO mp = 1, plasma % NMP , NMP_reduce_factor
          print *, 'GHGM mp ', mp
        DO lp = 1, plasma % NLP , NLP_reduce_factor
          print *, 'GHGM lp ', lp
!     DO mp = 1, 2    
!       DO lp = 1, 2    

          ! Copy over the grid information (for now)
          ZX(1:grid % flux_tube_max(lp))  = grid % altitude(1:grid % flux_tube_max(lp),lp)/1000.0_prec !convert from m to km
          PCO                             = grid % p_value(lp)  !Pvalue is a single value !*** Need PValue attribute for the grid
          SLX(1:grid % flux_tube_max(lp)) = grid % foot_point_distance(1:grid % flux_tube_max(lp),lp,mp)
          GLX(1:grid % flux_tube_max(lp)) = 0.5_prec*pi - grid % magnetic_colatitude(1:grid % flux_tube_max(lp),lp)  ! magnetic latitude [radians]
          BMX(1:grid % flux_tube_max(lp)) = grid % magnetic_field_strength(1:grid % flux_tube_max(lp),lp,mp)   !Tesla
          GRX(1:grid % flux_tube_max(lp)) = grid % grx(1:grid % flux_tube_max(lp),lp,mp)

          ! Copy over neutrals 
          OX(1:grid % flux_tube_max(lp))    = neutrals % oxygen(1:grid % flux_tube_max(lp),lp,mp) !(m-3)
          HX(1:grid % flux_tube_max(lp))    = neutrals % hydrogen(1:grid % flux_tube_max(lp),lp,mp)
          N2X(1:grid % flux_tube_max(lp))   = neutrals % molecular_nitrogen(1:grid % flux_tube_max(lp),lp,mp)
          O2X(1:grid % flux_tube_max(lp))   = neutrals % molecular_oxygen(1:grid % flux_tube_max(lp),lp,mp)
          HEX(1:grid % flux_tube_max(lp))   = neutrals % helium(1:grid % flux_tube_max(lp),lp,mp)
          N4SX(1:grid % flux_tube_max(lp))  = neutrals % nitrogen(1:grid % flux_tube_max(lp),lp,mp)
          TNX(1:grid % flux_tube_max(lp))   = neutrals % temperature(1:grid % flux_tube_max(lp),lp,mp)
          TINFX(1:grid % flux_tube_max(lp)) = neutrals % temperature_inf(1:grid % flux_tube_max(lp),lp,mp)
          UNX(1:grid % flux_tube_max(lp))   = -neutrals % velocity_apex(3,1:grid % flux_tube_max(lp),lp,mp)

          SZA(1:grid % flux_tube_max(lp)) = Solar_Zenith_Angle( time_tracker % utime, &
                                                                time_tracker % day_of_year, &
                                                                grid % magnetic_colatitude(1:grid % flux_tube_max(lp),lp), &
                                                                grid % longitude(1:grid % flux_tube_max(lp),lp,mp), &
                                                                grid % flux_tube_max(lp) )
          
          AUR_PROD(1:3,1:grid % flux_tube_max(lp)) = plasma % ionization_rates(2:4,1:grid % flux_tube_max(lp),lp,mp)
 
          EFLAG(1:11,1:11) = 0

          DO i=1, grid % nFluxTube

            ! Ion Densities
            XIONNX(1:9,i) = plasma % ion_densities(1:9,i,lp,mp)
            ! Along Flux Tube Ion Velocities
            XIONVX(1:9,i) = plasma % ion_velocities(3,1:9,i,lp,mp)

            ! Ion Temperatures
            TE_TIX(1,i) = plasma % ion_temperature(i,lp,mp)
            TE_TIX(2,i) = plasma % ion_temperature(i,lp,mp)
            ! Electron Temperature
            TE_TIX(3,i) = plasma % electron_temperature(i,lp,mp)

            EHTX(1:3,i) = 0.0_dp
            NNOX(i)     = 0.0_dp
            NHEAT(i)    = 0.0_dp
       
          ENDDO

          JMINX = 1
          JMAXX = grid % flux_tube_max(lp)

if(mp.eq.1.and.(lp.eq.168.or.lp.eq.169)) then
print *, ' GHGM LP ',lp
print *, ' GHGM1 ', JMINX
print *, ' GHGM2 ', JMAXX
print *, ' GHGM2.5 ', grid % flux_tube_max(lp)
print *, ' GHGM3 ', ZX(1:JMAXX)
print *, ' GHGM4 ', PCO
print *, ' GHGM5 ', SLX(1:JMAXX)
print *, ' GHGM6 ', GLX(1:JMAXX)
print *, ' GHGM7 ', BMX(1:JMAXX)
print *, ' GHGM8 ', GRX(1:JMAXX)
print *, ' GHGM9 ', OX(1:JMAXX)
print *, ' GHGM10 ', HX(1:JMAXX)
print *, ' GHGM11 ', N2X(1:JMAXX)
print *, ' GHGM12 ', O2X(1:JMAXX)
print *, ' GHGM13 ', HEX(1:JMAXX)
print *, ' GHGM14 ', N4SX(1:JMAXX)
print *, ' GHGM15 ', INNO
print *, ' GHGM16 ', NNOX(1:JMAXX)
print *, ' GHGM17 ', TNX(1:JMAXX)
print *, ' GHGM18 ', TINFX(1:JMAXX)
print *, ' GHGM19 ', UNX(1:JMAXX)
print *, ' GHGM20 ', flip_time_step
print *, ' GHGM21 ', DTMIN
print *, ' GHGM22 ', F107D
print *, ' GHGM23 ', F107A
print *, ' GHGM24 ', SZA(1:JMAXX)
print *, ' GHGM25 ', FPAS
print *, ' GHGM26 ', HPEQ
print *, ' GHGM27 ', HEPRAT
print *, ' GHGM28 ', COLFACX
print *, ' GHGM29 ', IHEPLS
print *, ' GHGM30 ', INPLS
print *, ' GHGM31 ', UTHR
print *, ' GHGM32 ', EHTX(1:3,1:JMAXX)
print *, ' GHGM33 ', AUR_PROD(1:3,1:JMAXX)
print *, ' GHGM34 ', TE_TIX(1:3,1:JMAXX)
print *, ' GHGM35 ', XIONNX(1:9,1:JMAXX)
print *, ' GHGM36 ', XIONVX(1:9,1:JMAXX)
print *, ' GHGM37 ', NHEAT(1:JMAXX)
print *, ' GHGM38 ', EFLAG
endif
!if(mp.eq.1.and.lp.eq.169) stop
 
          CALL CTIPINT( JMINX, & !.. index of the first point on the field line
                        JMAXX, & !.. index of the last point on the field line
                        grid % flux_tube_max(lp), & !.. CTIPe array dimension, must equal to FLDIM
                        ZX(1:JMAXX), & !.. array, altitude (km)
                        PCO, & !.. p coordinate (L-shell)
                        SLX(1:JMAXX), & !.. array, distance of point from northern hemisphere (meter)
                        GLX(1:JMAXX), & !.. array, magnetic latitude (radians)
                        BMX(1:JMAXX), & !.. array, magnetic field strength, (Tesla)
                        GRX(1:JMAXX), & !.. array, gravity, m2 s-1
                        OX(1:JMAXX), & !.. array, O density (m-3)
                        HX(1:JMAXX), & !.. array, H density (m-3)
                        N2X(1:JMAXX), & !.. array, N2 density (cm-3)
                        O2X(1:JMAXX), & !.. array, O2 density (cm-3)
                        HEX(1:JMAXX), & !.. array, He density (cm-3)
                        N4SX(1:JMAXX), & !.. array, N(4S) density (cm-3)
                        INNO, & !.. switch to turn on FLIP NO calculation IF <0
                        NNOX(1:JMAXX), & !.. array, NO density (cm-3)
                        TNX(1:JMAXX), & !.. array, Neutral temperature (K)
                        TINFX(1:JMAXX), & !.. array, Exospheric Neutral temperature (K)
                        UNX(1:JMAXX), & !.. array, Neutral wind (m/s), field aligned component, positive SOUTHward
                        flip_time_step, & !.. CTIPe time step (secs)
                        DTMIN, & !.. Minimum time step allowed (>=10 secs?)
                        F107D, & !.. Daily F10.7
                        F107A, & !.. 81 day average F10.7
                        SZA(1:JMAXX), & !.. Solar Zenith angle (radians)
                        FPAS, & !.. Pitch angle scattering fraction
                        HPEQ, & !.. Sets initial equatorial H+ density. See declaration below
                        HEPRAT, & !.. Intial He+/H+ ratio (.01 to 1.0)
                        COLFACX, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
                        IHEPLS, & !.. switches He+ dIFfusive solution on IF > 0
                        INPLS, & !.. switches N+ dIFfusive solution on IF > 0
                        UTHR, &  !.. Universal time in hours 
                        EHTX(1:3,1:JMAXX), & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
                        AUR_PROD(1:3,1:JMAXX), & ! IN 2D array, ionization rates for O+, O2+, and N2+
                        TE_TIX(1:3,1:JMAXX), & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
                        XIONNX(1:9,1:JMAXX), &
                        XIONVX(1:9,1:JMAXX), & !.. IN/OUT: 2D array, Storage for ion densities and velocities
                        NHEAT(1:JMAXX), & !.. OUT: array, Neutral heating rate (eV/cm^3/s)
                        EFLAG ) !.. OUT: 2D array, Error Flags


          DO i=1, grid % nFluxTube

            ! Ion Densities
            plasma % ion_densities(1:9,i,lp,mp) = XIONNX(1:9,i)
            ! Along Flux Tube Ion Velocities
            plasma % ion_velocities(3,1:9,i,lp,mp) = XIONVX(1:9,i)

            ! Ion Temperatures
            plasma % ion_temperature(i,lp,mp) = TE_TIX(1,i)
           ! plasma % ion_temperature(i,lp,mp) = TE_TIX(2,i)
            ! Electron Temperature
            plasma % electron_temperature(i,lp,mp) = TE_TIX(3,i) 

!         if (mp.eq.1.and.lp.eq.1) then
!         write(6,566) i , zx(i), XIONNX(1,i),XIONNX(2,i)                                              
!566 format(' GHGM ions ',i4,3e12.4)
!         endif

          ENDDO

      ENDDO
    ENDDO

    istop = 0
    if(istop.eq.1) stop


  END SUBROUTINE FLIP_Wrapper
!
  FUNCTION Solar_Zenith_Angle( utime, day, colatitude, longitude, flux_tube_max ) RESULT( sza )
    IMPLICIT NONE
    REAL(prec) :: utime
    INTEGER    :: day
    INTEGER    :: flux_tube_max
    REAL(prec) :: colatitude(1:flux_tube_max) 
    REAL(prec) :: longitude(1:flux_tube_max) 
    REAL(dp)   :: sza(1:flux_tube_max)
    ! Local
    REAL(dp) :: ty, sda, ssa, rlt, rlat, utime12, cos_sza
    INTEGER  :: i


        DO i=1, flux_tube_max

          ty = ( REAL(day, prec) + 15.5_prec )*12.0_prec/365.0_prec
          IF ( ty > 12.0_prec )THEN 
            ty = ty - 12.0_prec
          ENDIF

          sda = ATAN(0.434_prec*SIN(pi/6.0_prec*(ty-3.17_prec)))   

          IF ( utime >= 43200.0_prec ) THEN
            utime12 = utime - 43200.0_prec
          ELSE
            utime12 = utime + 43200.0_prec
          ENDIF

          ssa = longitude(i)*180.0_prec/pi + utime12/240.0_prec
          rlt = 180.0_prec + ssa
          IF ( rlt > 360.0_prec ) THEN
            rlt = rlt - 360.0_prec
          ENDIF
          rlt = rlt*pi/180.0_prec

          rlat = 0.5_prec*pi - colatitude(i) 
          cos_sza = -COS(rlat)*COS(sda)*COS(rlt)+SIN(rlat)*SIN(sda)

          sza(i)  = ACOS( cos_sza )

       END DO

  END FUNCTION Solar_Zenith_Angle

  SUBROUTINE Interpolate_to_GeographicGrid_IPE_Plasma( plasma, grid )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)         :: grid
    ! Local
    INTEGER :: i

      DO i = 1, n_ion_species

        CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_densities(i,:,:,:), plasma % geo_ion_densities(i,:,:,:) )
        CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_velocities(1,i,:,:,:), plasma % geo_ion_velocities(1,i,:,:,:) )
        CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_velocities(2,i,:,:,:), plasma % geo_ion_velocities(2,i,:,:,:) )
        CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_velocities(3,i,:,:,:), plasma % geo_ion_velocities(3,i,:,:,:) )

     ENDDO

     CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_temperature, plasma % geo_ion_temperature )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_density, plasma % geo_electron_density )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_temperature, plasma % geo_electron_temperature )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_velocity(1,:,:,:), plasma % geo_electron_velocity(1,:,:,:) )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_velocity(2,:,:,:), plasma % geo_electron_velocity(2,:,:,:) )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_velocity(3,:,:,:), plasma % geo_electron_velocity(3,:,:,:) )

     CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(1,:,:,:), plasma % geo_ionization_rates(1,:,:,:) )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(2,:,:,:), plasma % geo_ionization_rates(2,:,:,:) )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(3,:,:,:), plasma % geo_ionization_rates(3,:,:,:) )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(4,:,:,:), plasma % geo_ionization_rates(4,:,:,:) )


  END SUBROUTINE Interpolate_to_GeographicGrid_IPE_Plasma

  SUBROUTINE Read_Legacy_Input_IPE_Plasma( plasma, grid, filename )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    CHARACTER(100), INTENT(in)         :: filename
    ! Local
    INTEGER               :: ispec, i, lp, mp, ii, fUnit
    REAL(sp), ALLOCATABLE :: dumm(:,:,:) 

      ALLOCATE( dumm(1:grid % npts2d, 1:grid % NMP, 1:12) )      
      
      dumm = 0.0_sp

      OPEN( UNIT   = NewUnit( fUnit ), &
            FILE   = TRIM( filename ), &
            FORM   = 'UNFORMATTED', &
            STATUS = 'OLD' )


      DO ispec=1,12
        READ( fUnit ) dumm(:,:,ispec)
      ENDDO

      CLOSE( fUnit )

      DO mp = 1, grid % NMP , NMP_reduce_factor

        ii = 0

        DO lp = 1, grid % NLP , NLP_reduce_factor
          DO i = 1, grid % flux_tube_max(lp)

            ii = ii + 1
!            PRINT*, i, lp, mp, ii, dumm(ii,mp,11)
            plasma % ion_densities(1:9,i,lp,mp)    = dumm(ii,mp,1:9)
            plasma % electron_temperature(i,lp,mp) = dumm(ii,mp,10)
            plasma % ion_temperature(i,lp,mp)      = dumm(ii,mp,11)

          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE( dumm )

  END SUBROUTINE Read_Legacy_Input_IPE_Plasma

END MODULE IPE_Plasma_Class

