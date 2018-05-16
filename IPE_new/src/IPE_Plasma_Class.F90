MODULE IPE_Plasma_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Grid_Class
USE IPE_Neutrals_Class

IMPLICIT NONE

  TYPE IPE_Plasma
    INTEGER :: nFluxTube, NLP, NMP

    REAL(prec), ALLOCATABLE :: ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: hall_conductivity(:,:,:) 
    REAL(prec), ALLOCATABLE :: pedersen_conductivity(:,:,:) 
    REAL(prec), ALLOCATABLE :: b_parallel_conductivity(:,:,:) 

    ! Interpolated Fields
    REAL(prec), ALLOCATABLE :: geo_ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ion_velocities(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ion_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_hall_conductivity(:,:,:) 
    REAL(prec), ALLOCATABLE :: geo_pedersen_conductivity(:,:,:) 
    REAL(prec), ALLOCATABLE :: geo_b_parallel_conductivity(:,:,:) 
    

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Plasma
      PROCEDURE :: Trash => Trash_IPE_Plasma

      PROCEDURE :: Update => Update_IPE_Plasma
      PROCEDURE :: FLIP_Wrapper     

      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Plasma

  END TYPE IPE_Plasma


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
                plasma % hall_conductivity(1:nFluxTube,1:NLP,1:NMP), &
                plasma % pedersen_conductivity(1:nFluxTube,1:NLP,1:NMP), &
                plasma % b_parallel_conductivity(1:nFluxTube,1:NLP,1:NMP) )
             
      plasma % ion_densities           = safe_density_minimum 
      plasma % ion_velocities          = 0.0_prec
      plasma % ion_temperature         = safe_temperature_minimum
      plasma % electron_density        = safe_density_minimum
      plasma % electron_velocity       = 0.0_prec
      plasma % electron_temperature    = safe_temperature_minimum
      plasma % hall_conductivity       = 0.0_prec
      plasma % pedersen_conductivity   = 0.0_prec
      plasma % b_parallel_conductivity = 0.0_prec

      ALLOCATE( plasma % geo_ion_densities(1:n_ion_species,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ion_velocities(1:3,1:n_ion_species,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ion_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_density(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_velocity(1:3,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_hall_conductivity(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_pedersen_conductivity(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_b_parallel_conductivity(1:nlon_geo,1:nlat_geo,1:nheights_geo)  )

      plasma % geo_ion_densities           = 0.0_prec
      plasma % geo_ion_velocities          = 0.0_prec
      plasma % geo_ion_temperature         = 0.0_prec
      plasma % geo_electron_density        = 0.0_prec
      plasma % geo_electron_temperature    = 0.0_prec
      plasma % geo_electron_velocity       = 0.0_prec
      plasma % geo_hall_conductivity       = 0.0_prec
      plasma % geo_pedersen_conductivity   = 0.0_prec
      plasma % geo_b_parallel_conductivity = 0.0_prec

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
                plasma % hall_conductivity, & 
                plasma % pedersen_conductivity, & 
                plasma % b_parallel_conductivity, &
                plasma % geo_ion_velocities, &
                plasma % geo_ion_temperature, &
                plasma % geo_electron_density, &
                plasma % geo_electron_velocity, &
                plasma % geo_electron_temperature, &
                plasma % geo_hall_conductivity, & 
                plasma % geo_pedersen_conductivity, & 
                plasma % geo_b_parallel_conductivity )

  END SUBROUTINE Trash_IPE_Plasma
!
  SUBROUTINE Update_IPE_Plasma( plasma, grid, neutrals, utime, year, day, flip_time_step, f107d, f107a, sun_longitude )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    REAL(prec), INTENT(in)             :: utime
    INTEGER, INTENT(in)                :: year
    INTEGER, INTENT(in)                :: day
    REAL(prec), INTENT(in)             :: flip_time_step
    REAL(prec), INTENT(in)             :: f107d
    REAL(prec), INTENT(in)             :: f107a
    REAL(prec), INTENT(in)             :: sun_longitude
    ! Local

         
      !CALL plasma % Cross_Flux_Tube_Transport( )
 
      CALL plasma % FLIP_Wrapper( grid, & 
                                  neutrals, &
                                  utime, year, day, &
                                  flip_time_step, &
                                  f107d, f107a, &
                                  sun_longitude )


  END SUBROUTINE Update_IPE_Plasma
!
  SUBROUTINE FLIP_Wrapper( plasma, grid, neutrals, utime, year, day, flip_time_step, f107d, f107a, sun_longitude )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    REAL(prec), INTENT(in)             :: utime
    INTEGER, INTENT(in)                :: year
    INTEGER, INTENT(in)                :: day
    REAL(prec), INTENT(in)             :: flip_time_step
    REAL(prec), INTENT(in)             :: f107d
    REAL(prec), INTENT(in)             :: f107a
    REAL(prec), INTENT(in)             :: sun_longitude
    ! Local
    INTEGER  :: i, lp, mp
    INTEGER  :: JMINX, JMAXX
    INTEGER  :: EFLAG(11,11)
    REAL(dp) :: PCO, ltime, mlt
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
    REAL(dp) :: NNOX(1:grid % nFluxTube) 
    REAL(dp) :: NHEAT(1:grid % nFluxTube) 
    REAL(dp) :: SZA(1:grid % nFluxTube) 
    
      DO mp = 1, plasma % NMP    
        DO lp = 1, plasma % NLP    

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

          SZA(1:grid % flux_tube_max(lp)) = Solar_Zenith_Angle( utime, &
                                                                day, &
                                                                grid % magnetic_colatitude(1:grid % flux_tube_max(lp),lp), &
                                                                grid % longitude(1:grid % flux_tube_max(lp),lp,mp), &
                                                                grid % flux_tube_max(lp) )

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

          ltime = utime/3600.0_prec + grid % longitude(i,lp,mp)*180.0_prec/pi/15.0_prec
          IF ( ltime > 24.0_prec ) THEN
            ltime = MOD(ltime, 24.0_prec)
          ENDIF

          mlt = grid % magnetic_longitude(mp)*180.0_prec/pi/15.0_prec - &
                sun_longitude*12.0_prec/pi+12.0_prec 

          JMINX = 1
          JMAXX = grid % flux_tube_max(lp)

          CALL CTIPINT( JMINX, & !.. index of the first point on the field line
                        JMAXX, & !.. index of the last point on the field line
                        grid % flux_tube_max(lp), & !.. CTIPe array dimension, must equal to FLDIM
                        ZX, & !.. array, altitude (km)
                        PCO, & !.. p coordinate (L-shell)
                        SLX, & !.. array, distance of point from northern hemisphere (meter)
                        GLX, & !.. array, magnetic latitude (radians)
                        BMX, & !.. array, magnetic field strength, (Tesla)
                        GRX, & !.. array, gravity, m2 s-1
                        OX, & !.. array, O density (m-3)
                        HX, & !.. array, H density (m-3)
                        N2X, & !.. array, N2 density (cm-3)
                        O2X, & !.. array, O2 density (cm-3)
                        HEX, & !.. array, He density (cm-3)
                        N4SX, & !.. array, N(4S) density (cm-3)
                        INNO, & !.. switch to turn on FLIP NO calculation IF <0
                        NNOX, & !.. array, NO density (cm-3)
                        TNX, & !.. array, Neutral temperature (K)
                        TINFX, & !.. array, Exospheric Neutral temperature (K)
                        UNX, & !.. array, Neutral wind (m/s), field aligned component, positive SOUTHward
                        flip_time_step, & !.. CTIPe time step (secs)
                        DTMIN, & !.. Minimum time step allowed (>=10 secs?)
                        F107D, & !.. Daily F10.7
                        F107A, & !.. 81 day average F10.7
                        SZA, & !.. Solar Zenith angle (radians)
                        FPAS, & !.. Pitch angle scattering fraction
                        HPEQ, & !.. Sets initial equatorial H+ density. See declaration below
                        HEPRAT, & !.. Intial He+/H+ ratio (.01 to 1.0)
                        COLFACX, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
                        IHEPLS, & !.. switches He+ dIFfusive solution on IF > 0
                        INPLS, & !.. switches N+ dIFfusive solution on IF > 0
                        EHTX, & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
                        TE_TIX, & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
                        XIONNX,XIONVX, & !.. IN/OUT: 2D array, Storage for ion densities and velocities
                        NHEAT, & !.. OUT: array, Neutral heating rate (eV/cm^3/s)
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

          ENDDO

      ENDDO
    ENDDO


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

     CALL grid % Interpolate_to_Geographic_Grid( plasma % hall_conductivity, plasma % geo_hall_conductivity )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % pedersen_conductivity, plasma % geo_pedersen_conductivity )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % b_parallel_conductivity, plasma % geo_b_parallel_conductivity )

  END SUBROUTINE Interpolate_to_GeographicGrid_IPE_Plasma

END MODULE IPE_Plasma_Class
