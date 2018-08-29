MODULE IPE_Model_Class

USE IPE_Precision
USE IPE_Model_Parameters_Class
USE IPE_Time_Class
USE IPE_Grid_Class
USE IPE_Neutrals_Class
USE IPE_Forcing_Class
USE IPE_Plasma_Class
USE IPE_Electrodynamics_Class

USE netcdf

USE omp_lib

IMPLICIT NONE


  ! The IPE_Model serves as a wrapper for all of the underlying attributes.
  ! This class should be used to orchestrate model setup, updating, and
  ! breakdown. At this level, we define the API for interacting with the
  ! deeper attributes within IPE.

  TYPE IPE_Model

    TYPE( IPE_Time )             :: time_tracker
    TYPE( IPE_Model_Parameters ) :: parameters
    TYPE( IPE_Grid )             :: grid
    TYPE( IPE_Forcing )          :: forcing
    TYPE( IPE_Neutrals )         :: neutrals
    TYPE( IPE_Plasma )           :: plasma
    TYPE( IPE_Electrodynamics )  :: eldyn

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model
      PROCEDURE :: Trash => Trash_IPE_Model

      PROCEDURE :: Update => Update_IPE_Model
      PROCEDURE :: Write_NetCDF_IPE
      PROCEDURE :: Read_NetCDF_IPE
      PROCEDURE :: Write_Geographic_NetCDF_IPE

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: Geographic_Interpolation
   
  END TYPE IPE_Model


  REAL(prec), PARAMETER :: fillValue = -999999.9_prec

CONTAINS

  SUBROUTINE Build_IPE_Model( ipe, init_success )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(out)    :: ipe
    LOGICAL, INTENT(out)               :: init_success
    ! Local 
    CHARACTER(200) :: init_file
    LOGICAL        :: fileExists

      CALL ipe % parameters % Build( init_success )

      ! Initialize the clock
      CALL ipe % time_tracker % Build( ipe % parameters % initial_timestamp )

      IF( init_success )THEN

        ! ////// Forcing ////// !

        CALL ipe % forcing % Build( ipe % parameters % solar_forcing_time_step, &
                                    ipe % parameters % f107_kp_size )

        IF( ipe % parameters % use_f107_kp_file )THEN
          
          CALL ipe % forcing % Read_F107KP_IPE_Forcing( ipe % parameters % f107_kp_skip_size, &
                                                        ipe % parameters % f107_kp_file )

        ENDIF

        CALL ipe % forcing % Read_Tiros_IPE_Forcing( )


        ! ////// grid ////// !

        CALL ipe % grid % Read_IPE_Grid_NetCDF( ipe % parameters % netcdf_grid_file  )
        ipe % grid % npts2d = ipe % parameters % npts2d


        ! ////// neutrals ////// !

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


        ! ////// electric field /////// !

        CALL ipe % eldyn % Build( nFluxTube = ipe % grid % nFluxTube,&
                                  NLP       = ipe % grid % NLP, &
                                  NMP       = ipe % grid % NMP )


        ! ////// plasma ////// !

        CALL ipe % plasma % Build( nFluxTube = ipe % grid % nFluxTube, &
                                   NLP       = ipe % grid % NLP, &
                                   NMP       = ipe % grid % NMP )
                                      
      ELSE

        init_success = .FALSE.
        RETURN
        
      ENDIF

      init_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".nc" 
      INQUIRE( FILE = TRIM(init_file), EXIST = fileExists )
  
      IF( fileExists )THEN
        PRINT*, ' IPE_Model_Class : Build : Reading '//init_file
        CALL ipe % Read_NetCDF_IPE( init_file )
        PRINT*, ' IPE_Model_Class : Build : done    '//init_file
      ELSE
        PRINT*, ' IPE_Model_Class : Build : File not found : '//init_file
      ENDIF

  END SUBROUTINE Build_IPE_Model
!
  SUBROUTINE Trash_IPE_Model( ipe )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe

      CALL ipe % forcing  % Trash( )
      CALL ipe % grid     % Trash( )
      CALL ipe % neutrals % Trash( )
      CALL ipe % plasma   % Trash( )
      CALL ipe % eldyn    % Trash( )

  END SUBROUTINE Trash_IPE_Model
!
  SUBROUTINE Update_IPE_Model( ipe, t0, t1 )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe
    REAL(prec), INTENT(in)            :: t0, t1
    ! Local
    REAL(prec) :: AP(1:7)
    REAL(prec) :: wt1, wt2


      CALL ipe % time_tracker % Update( t0 )
      
      ! Need to add a call to update the index for capturing AP
      AP = ipe % forcing % GetAP( t0 )

      print *, 'GHGM updating neutrals'
      CALL CPU_TIME( wt1 )
      CALL ipe % neutrals % Update( ipe % grid, &
                                    ipe % time_tracker % utime, &
                                    ipe % time_tracker % year, &
                                    ipe % time_tracker % day_of_year,  & 
                                    ipe % forcing % f107( ipe % forcing % current_index ) , &
                                    ipe % forcing % f107_81day_avg( ipe % forcing % current_index ), &
                                    AP )
      CALL CPU_TIME( wt2 )
      PRINT*, 'Neutral : ', wt2-wt1

      print *, 'GHGM updating electrodynamics'
      CALL CPU_TIME( wt1 )
      CALL ipe % eldyn % Update( ipe % grid, &
                                 ipe % forcing, &
                                 ipe % time_tracker )
      CALL CPU_TIME( wt2 )
      PRINT*, 'Eldyn : ', wt2-wt1
   
      print *, 'GHGM updating plasma'
      CALL CPU_TIME( wt1 )
      CALL ipe % plasma % Update( ipe % grid, &
                                  ipe % neutrals, &
                                  ipe % forcing, &
                                  ipe % time_tracker, &
                                  ipe % eldyn % v_ExB_apex, &
                                  ipe % parameters % solar_forcing_time_step )

      CALL CPU_TIME( wt2 )
      PRINT*, 'Plasma : ', wt2-wt1

      ! Update the timer
      CALL ipe % time_tracker % Update( t1 )

  END SUBROUTINE Update_IPE_Model
!
  SUBROUTINE Geographic_Interpolation( ipe )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe


     IF( ipe % parameters % write_geographic_neutrals )THEN
       CALL ipe % neutrals % Interpolate_to_GeographicGrid( ipe % grid )
     ENDIF

     IF( ipe % parameters % write_geographic_eldyn )THEN
       CALL ipe % eldyn % Interpolate_to_GeographicGrid( ipe % grid )
     ENDIF

     CALL ipe % plasma % Interpolate_to_GeographicGrid( ipe % grid )

  END SUBROUTINE Geographic_Interpolation

  SUBROUTINE Write_NetCDF_IPE( ipe, filename )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(in) :: ipe
    CHARACTER(*), INTENT(in)       :: filename
    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid, time_varid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid
    INTEGER :: op_varid, hp_varid, hep_varid, np_varid, n2p_varid, o2p_varid, nop_varid
    INTEGER :: op2d_varid, op2p_varid, ion_temp_varid, phi_varid, mhd_phi_varid, hc_varid, pc_varid, bc_varid
    INTEGER :: opv_varid(1:3), hpv_varid(1:3), hepv_varid(1:3), npv_varid(1:3), n2pv_varid(1:3), o2pv_varid(1:3), nopv_varid(1:3)
    INTEGER :: ion_vel_op_varid, ion_vel_hp_varid, ion_vel_hep_varid
    INTEGER :: recStart(1:4), recCount(1:4)


      recStart = (/ 1, 1, 1, 1 /)
      recCount = (/ ipe % grid % nFluxTube, ipe % grid % NLP, ipe % grid % NMP, 1 /)


      time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
      CALL Check( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )

      CALL Check( nf90_def_dim( ncid, "s", ipe % grid % nFluxTube, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "lp", ipe % grid % NLP, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "mp", ipe % grid % NMP, y_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )


      CALL Check( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
      CALL Check( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )

      IF( ipe % parameters % write_apex_neutrals )THEN

        ! Neutrals
        CALL Check( nf90_def_var( ncid, "helium", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , helium_varid ) )
        CALL Check( nf90_put_att( ncid, helium_varid, "long_name", "Neutral Helium Density" ) )
        CALL Check( nf90_put_att( ncid, helium_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "oxygen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , oxygen_varid ) )
        CALL Check( nf90_put_att( ncid, oxygen_varid, "long_name", "Neutral Oxygen Density" ) )
        CALL Check( nf90_put_att( ncid, oxygen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "molecular_oxygen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , molecular_oxygen_varid ) )
        CALL Check( nf90_put_att( ncid, molecular_oxygen_varid, "long_name", "Neutral Molecular Oxygen Density" ) )
        CALL Check( nf90_put_att( ncid, molecular_oxygen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "molecular_nitrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , molecular_nitrogen_varid ) )
        CALL Check( nf90_put_att( ncid, molecular_nitrogen_varid, "long_name", "Neutral Molecular Nitrogen Density" ) )
        CALL Check( nf90_put_att( ncid, molecular_nitrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "nitrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , nitrogen_varid ) )
        CALL Check( nf90_put_att( ncid, nitrogen_varid, "long_name", "Neutral Nitrogen Density" ) )
        CALL Check( nf90_put_att( ncid, nitrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "hydrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , hydrogen_varid ) )
        CALL Check( nf90_put_att( ncid, hydrogen_varid, "long_name", "Neutral Hydrogen Density" ) )
        CALL Check( nf90_put_att( ncid, hydrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "temperature", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , temperature_varid ) )
        CALL Check( nf90_put_att( ncid, temperature_varid, "long_name", "Thermosphere Temperature" ) )
        CALL Check( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

        CALL Check( nf90_def_var( ncid, "u", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , u_varid ) )
        CALL Check( nf90_put_att( ncid, u_varid, "long_name", "Zonal Velocity" ) )
        CALL Check( nf90_put_att( ncid, u_varid, "units", "m s^{-1}" ) )

        CALL Check( nf90_def_var( ncid, "v", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , v_varid ) )
        CALL Check( nf90_put_att( ncid, v_varid, "long_name", "Meridional Velocity" ) )
        CALL Check( nf90_put_att( ncid, v_varid, "units", "m s^{-1}" ) )

        CALL Check( nf90_def_var( ncid, "w", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , w_varid ) )
        CALL Check( nf90_put_att( ncid, w_varid, "long_name", "Radial Velocity" ) )
        CALL Check( nf90_put_att( ncid, w_varid, "units", "m s^{-1}" ) )

      ENDIF

      IF( ipe % parameters % write_apex_eldyn )THEN

        CALL Check( nf90_def_var( ncid, "phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , phi_varid ) )
        CALL Check( nf90_put_att( ncid, phi_varid, "long_name", "Electric Potential" ) )
        CALL Check( nf90_put_att( ncid, phi_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "mhd_phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) ,mhd_phi_varid ) )
        CALL Check( nf90_put_att( ncid, mhd_phi_varid, "long_name", "Electric Potential - MHD Component" ) )
        CALL Check( nf90_put_att( ncid, mhd_phi_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "hc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hc_varid ) )
        CALL Check( nf90_put_att( ncid, hc_varid, "long_name", "Hall Conductivity" ) )
        CALL Check( nf90_put_att( ncid, hc_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "pc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , pc_varid ) )
        CALL Check( nf90_put_att( ncid, pc_varid, "long_name", "Pedersen Conductivity" ) )
        CALL Check( nf90_put_att( ncid, pc_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "bc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , bc_varid ) )
        CALL Check( nf90_put_att( ncid, bc_varid, "long_name", "Magnetic Field Aligned Conductivity" ) )
        CALL Check( nf90_put_att( ncid, bc_varid, "units", "[Unknown]" ) )

      ENDIF

      ! Plasma
      CALL Check( nf90_def_var( ncid, "O+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , op_varid ) )
      CALL Check( nf90_put_att( ncid, op_varid, "long_name", "Atomic oxygen ion number density (ground state)" ) )
      CALL Check( nf90_put_att( ncid, op_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "H+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , hp_varid ) )
      CALL Check( nf90_put_att( ncid, hp_varid, "long_name", "Hydrogen ion number density" ) )
      CALL Check( nf90_put_att( ncid, hp_varid, "units", " m^{-3}" ) )
      
      CALL Check( nf90_def_var( ncid, "He+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , hep_varid ) )
      CALL Check( nf90_put_att( ncid, hep_varid, "long_name", "Helium ion number density" ) )
      CALL Check( nf90_put_att( ncid, hep_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "N+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , np_varid ) )
      CALL Check( nf90_put_att( ncid, np_varid, "long_name", "Nitrogen ion number density" ) )
      CALL Check( nf90_put_att( ncid, np_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "NO+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , nop_varid ) )
      CALL Check( nf90_put_att( ncid, nop_varid, "long_name", "Nitrosonium ion number density" ) )
      CALL Check( nf90_put_att( ncid, nop_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "O2+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , o2p_varid ) )
      CALL Check( nf90_put_att( ncid, o2p_varid, "long_name", "Molecular Oxygen ion number density" ) )
      CALL Check( nf90_put_att( ncid, o2p_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "N2+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , n2p_varid ) )
      CALL Check( nf90_put_att( ncid, n2p_varid, "long_name", "Molecular Nitrogen ion number density" ) )
      CALL Check( nf90_put_att( ncid, n2p_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "O+(2D)", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , op2d_varid ) )
      CALL Check( nf90_put_att( ncid, op2d_varid, "long_name", "Atomic oxygen ion number density (first excited state)" ) )
      CALL Check( nf90_put_att( ncid, op2d_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "O+(2P)", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , op2p_varid ) )
      CALL Check( nf90_put_att( ncid, op2p_varid, "long_name", "Atomic oxygen ion number density (second excited state)" ) )
      CALL Check( nf90_put_att( ncid, op2p_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "ion_temp", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_temp_varid ) )
      CALL Check( nf90_put_att( ncid, ion_temp_varid, "long_name", "Ion temperature" ) )
      CALL Check( nf90_put_att( ncid, ion_temp_varid, "units", "K" ) )

      CALL Check( nf90_def_var( ncid, "ion_vel_op", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_vel_op_varid ) )
      CALL Check( nf90_put_att( ncid, ion_vel_op_varid, "long_name", "O+ Velocity (B Parallel)" ) )
      CALL Check( nf90_put_att( ncid, ion_vel_op_varid, "units", "?/s" ) )

      CALL Check( nf90_def_var( ncid, "ion_vel_hp", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_vel_hp_varid ) )
      CALL Check( nf90_put_att( ncid, ion_vel_hp_varid, "long_name", "H+ Velocity (B Parallel)" ) )
      CALL Check( nf90_put_att( ncid, ion_vel_hp_varid, "units", "?/s" ) )

      CALL Check( nf90_def_var( ncid, "ion_vel_hep", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_vel_hep_varid ) )
      CALL Check( nf90_put_att( ncid, ion_vel_hep_varid, "long_name", "He+ Velocity (B Parallel)" ) )
      CALL Check( nf90_put_att( ncid, ion_vel_hep_varid, "units", "?/s" ) )


      CALL Check( nf90_enddef(ncid) )
      
      CALL Check( nf90_put_var( ncid, time_varid, time ) )
      IF( ipe % parameters % write_apex_neutrals )THEN

        CALL Check( nf90_put_var( ncid, helium_varid, ipe % neutrals % helium, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, oxygen_varid, ipe % neutrals % oxygen, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, molecular_oxygen_varid, ipe % neutrals % molecular_oxygen, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, molecular_nitrogen_varid, ipe % neutrals % molecular_nitrogen, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, nitrogen_varid, ipe % neutrals % nitrogen, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, hydrogen_varid, ipe % neutrals % hydrogen, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, temperature_varid, ipe % neutrals % temperature, recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, u_varid, ipe % neutrals % velocity_geographic(1,:,:,:), recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, v_varid, ipe % neutrals % velocity_geographic(2,:,:,:), recStart, recCount ) )
        CALL Check( nf90_put_var( ncid, w_varid, ipe % neutrals % velocity_geographic(3,:,:,:), recStart, recCount ) )

      ENDIF

      IF( ipe % parameters % write_apex_eldyn )THEN
        CALL Check( nf90_put_var( ncid, phi_varid, ipe % eldyn % electric_potential ) )
        CALL Check( nf90_put_var( ncid, mhd_phi_varid, ipe % eldyn % mhd_electric_potential ) )
        CALL Check( nf90_put_var( ncid, hc_varid, ipe % eldyn % hall_conductivity ) )
        CALL Check( nf90_put_var( ncid, pc_varid, ipe % eldyn % pedersen_conductivity ) )
        CALL Check( nf90_put_var( ncid, bc_varid, ipe % eldyn % b_parallel_conductivity ) )
      ENDIF

      CALL Check( nf90_put_var( ncid, op_varid, ipe % plasma % ion_densities(1,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, hp_varid, ipe % plasma % ion_densities(2,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, hep_varid, ipe % plasma % ion_densities(3,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, np_varid, ipe % plasma % ion_densities(4,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, nop_varid, ipe % plasma % ion_densities(5,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, o2p_varid, ipe % plasma % ion_densities(6,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, n2p_varid, ipe % plasma % ion_densities(7,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, op2d_varid, ipe % plasma % ion_densities(8,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, op2p_varid, ipe % plasma % ion_densities(9,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, ion_temp_varid, ipe % plasma % ion_temperature, recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, ion_vel_op_varid, ipe % plasma % ion_velocities(1,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, ion_vel_hp_varid, ipe % plasma % ion_velocities(2,:,:,:), recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, ion_vel_hep_varid, ipe % plasma % ion_velocities(3,:,:,:), recStart, recCount ) )

      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Write_NetCDF_IPE
!
  SUBROUTINE Read_NetCDF_IPE( ipe, filename )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*), INTENT(in)          :: filename
    ! Local
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: dimid, varid
    INTEGER :: nFluxtube, NLP, NMP
    INTEGER :: istop
    CHARACTER(NF90_MAX_NAME) :: nameHolder


      
      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))


      ! Obtain the dimensions of the IPE_Grid
      CALL Check( nf90_inq_dimid( ncid, "s", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, nFluxTube ) )

      CALL Check( nf90_inq_dimid( ncid, "lp", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, NLP ) )

      CALL Check( nf90_inq_dimid( ncid, "mp", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, NMP ) )

      ! It is assumed that the IPE_Model has already been constructed with
      ! appropriate NMP, NLP, nFluxTube.

      IF( ipe % parameters % write_apex_neutrals )THEN

        CALL Check( nf90_inq_varid( ncid, "helium", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % helium ) )

        CALL Check( nf90_inq_varid( ncid, "oxygen", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % oxygen ) )

        CALL Check( nf90_inq_varid( ncid, "molecular_oxygen", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % molecular_oxygen ) )

        CALL Check( nf90_inq_varid( ncid, "molecular_nitrogen", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % molecular_nitrogen ) )

        CALL Check( nf90_inq_varid( ncid, "nitrogen", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % nitrogen ) )

        CALL Check( nf90_inq_varid( ncid, "hydrogen", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % hydrogen ) )

        CALL Check( nf90_inq_varid( ncid, "temperature", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % temperature ) )

        CALL Check( nf90_inq_varid( ncid, "u", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % velocity_geographic(1,:,:,:) ) )

        CALL Check( nf90_inq_varid( ncid, "v", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % velocity_geographic(2,:,:,:) ) )

        CALL Check( nf90_inq_varid( ncid, "w", varid ) )
        CALL Check( nf90_get_var( ncid, varid, ipe % neutrals % velocity_geographic(3,:,:,:) ) )

      ENDIF


      CALL Check( nf90_inq_varid( ncid, "O+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "H+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "He+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(3,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "N+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(4,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "NO+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(5,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "O2+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(6,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "N2+", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(7,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "O+(2D)", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(8,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "O+(2P)", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_densities(9,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "ion_temp", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_temperature ) )

      CALL Check( nf90_inq_varid( ncid, "ion_vel_op", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_velocities(1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "ion_vel_hp", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_velocities(2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "ion_vel_hep", varid ) )
      CALL Check( nf90_get_var( ncid, varid, ipe % plasma % ion_velocities(3,:,:,:) ) )

      CALL Check( nf90_close( ncid ) )


  END SUBROUTINE Read_NetCDF_IPE

  SUBROUTINE Write_Geographic_NetCDF_IPE( ipe, filename )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*), INTENT(in)          :: filename
    ! Local 
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid
    INTEGER :: op_varid, hp_varid, hep_varid, np_varid, n2p_varid, o2p_varid, nop_varid
    INTEGER :: op2d_varid, op2p_varid, phi_varid, mhd_phi_varid, exbu_varid, exbv_varid
    INTEGER :: ion_temp_varid, e_varid, hc_varid, pc_varid, bc_varid
    INTEGER :: ion_rate_varid, O_rate_varid, O2_rate_varid, N2_rate_varid


      CALL ipe % Geographic_Interpolation( ) 

      ! Time from reference time is calculated here
      time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
      CALL Check( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )

      CALL Check( nf90_def_dim( ncid, "Z", nheights_geo, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

      CALL Check( nf90_def_var( ncid, "Z", NF90_PREC, z_dimid, z_varid ) )
      CALL Check( nf90_put_att( ncid, z_varid, "long_name", "Altitude" ) )
      CALL Check( nf90_put_att( ncid, z_varid, "units", "km" ) )
      CALL Check( nf90_put_att( ncid, z_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, z_varid, "missing_value", fillValue) )


      CALL Check( nf90_def_var( ncid, "latitude", NF90_PREC, y_dimid, y_varid ) )
      CALL Check( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
      CALL Check( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )
      CALL Check( nf90_put_att( ncid, y_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, y_varid, "missing_value", fillValue) )


      CALL Check( nf90_def_var( ncid, "longitude", NF90_PREC, x_dimid, x_varid ) )
      CALL Check( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
      CALL Check( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )
      CALL Check( nf90_put_att( ncid, x_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, x_varid, "missing_value", fillValue) )

      CALL Check( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
      CALL Check( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )


      IF( ipe % parameters % write_geographic_neutrals )THEN

        CALL Check( nf90_def_var( ncid, "helium", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , helium_varid ) )
        CALL Check( nf90_put_att( ncid, helium_varid, "long_name", "Neutral Helium Density" ) )
        CALL Check( nf90_put_att( ncid, helium_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "oxygen", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , oxygen_varid ) )
        CALL Check( nf90_put_att( ncid, oxygen_varid, "long_name", "Neutral Oxygen Density" ) )
        CALL Check( nf90_put_att( ncid, oxygen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "molecular_oxygen", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , molecular_oxygen_varid ) )
        CALL Check( nf90_put_att( ncid, molecular_oxygen_varid, "long_name", "Neutral Molecular Oxygen Density" ) )
        CALL Check( nf90_put_att( ncid, molecular_oxygen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "molecular_nitrogen", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , molecular_nitrogen_varid ) )
        CALL Check( nf90_put_att( ncid, molecular_nitrogen_varid, "long_name", "Neutral Molecular Nitrogen Density" ) )
        CALL Check( nf90_put_att( ncid, molecular_nitrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "nitrogen", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , nitrogen_varid ) )
        CALL Check( nf90_put_att( ncid, nitrogen_varid, "long_name", "Neutral Nitrogen Density" ) )
        CALL Check( nf90_put_att( ncid, nitrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "hydrogen", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hydrogen_varid ) )
        CALL Check( nf90_put_att( ncid, hydrogen_varid, "long_name", "Neutral Hydrogen Density" ) )
        CALL Check( nf90_put_att( ncid, hydrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "temperature", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , temperature_varid ) )
        CALL Check( nf90_put_att( ncid, temperature_varid, "long_name", "Thermosphere Temperature" ) )
        CALL Check( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

        CALL Check( nf90_def_var( ncid, "u", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , u_varid ) )
        CALL Check( nf90_put_att( ncid, u_varid, "long_name", "Zonal Velocity" ) )
        CALL Check( nf90_put_att( ncid, u_varid, "units", "m s^{-1}" ) )

        CALL Check( nf90_def_var( ncid, "v", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , v_varid ) )
        CALL Check( nf90_put_att( ncid, v_varid, "long_name", "Meridional Velocity" ) )
        CALL Check( nf90_put_att( ncid, v_varid, "units", "m s^{-1}" ) )

        CALL Check( nf90_def_var( ncid, "w", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , w_varid ) )
        CALL Check( nf90_put_att( ncid, w_varid, "long_name", "Radial Velocity" ) )
        CALL Check( nf90_put_att( ncid, w_varid, "units", "m s^{-1}" ) )

      ENDIF

      IF( ipe % parameters % write_geographic_eldyn )THEN

        CALL Check( nf90_def_var( ncid, "phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , phi_varid ) )
        CALL Check( nf90_put_att( ncid, phi_varid, "long_name", "Electric Potential" ) )
        CALL Check( nf90_put_att( ncid, phi_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "mhd_phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , mhd_phi_varid ) )
        CALL Check( nf90_put_att( ncid, mhd_phi_varid, "long_name", "Electric Potential - MHD Component" ) )
        CALL Check( nf90_put_att( ncid, mhd_phi_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "exb_u", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , exbu_varid ) )
        CALL Check( nf90_put_att( ncid, exbu_varid, "long_name", "Zonal component of ExB drift velocity" ) )
        CALL Check( nf90_put_att( ncid, exbu_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "exb_v", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , exbv_varid ) )
        CALL Check( nf90_put_att( ncid, exbv_varid, "long_name", "Meridional component of ExB drift velocity" ) )
        CALL Check( nf90_put_att( ncid, exbv_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "hc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hc_varid ) )
        CALL Check( nf90_put_att( ncid, hc_varid, "long_name", "Hall Conductivity" ) )
        CALL Check( nf90_put_att( ncid, hc_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "pc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , pc_varid ) )
        CALL Check( nf90_put_att( ncid, pc_varid, "long_name", "Pedersen Conductivity" ) )
        CALL Check( nf90_put_att( ncid, pc_varid, "units", "[Unknown]" ) )

        CALL Check( nf90_def_var( ncid, "bc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , bc_varid ) )
        CALL Check( nf90_put_att( ncid, bc_varid, "long_name", "Magnetic Field Aligned Conductivity" ) )
        CALL Check( nf90_put_att( ncid, bc_varid, "units", "[Unknown]" ) )

      ENDIF

      ! Plasma
      CALL Check( nf90_def_var( ncid, "O+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op_varid ) )
      CALL Check( nf90_put_att( ncid, op_varid, "long_name", "Atomic oxygen ion number density (ground state)" ) )
      CALL Check( nf90_put_att( ncid, op_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "H+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hp_varid ) )
      CALL Check( nf90_put_att( ncid, hp_varid, "long_name", "Hydrogen ion number density" ) )
      CALL Check( nf90_put_att( ncid, hp_varid, "units", " m^{-3}" ) )
      
      CALL Check( nf90_def_var( ncid, "He+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hep_varid ) )
      CALL Check( nf90_put_att( ncid, hep_varid, "long_name", "Helium ion number density" ) )
      CALL Check( nf90_put_att( ncid, hep_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "N+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , np_varid ) )
      CALL Check( nf90_put_att( ncid, np_varid, "long_name", "Nitrogen ion number density" ) )
      CALL Check( nf90_put_att( ncid, np_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "NO+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , nop_varid ) )
      CALL Check( nf90_put_att( ncid, nop_varid, "long_name", "Nitrosonium ion number density" ) )
      CALL Check( nf90_put_att( ncid, nop_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "O2+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , o2p_varid ) )
      CALL Check( nf90_put_att( ncid, o2p_varid, "long_name", "Molecular Oxygen ion number density" ) )
      CALL Check( nf90_put_att( ncid, o2p_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "N2+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , n2p_varid ) )
      CALL Check( nf90_put_att( ncid, n2p_varid, "long_name", "Molecular Nitrogen ion number density" ) )
      CALL Check( nf90_put_att( ncid, n2p_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "O+(2D)", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op2d_varid ) )
      CALL Check( nf90_put_att( ncid, op2d_varid, "long_name", "Atomic oxygen ion number density (first excited state)" ) )
      CALL Check( nf90_put_att( ncid, op2d_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "O+(2P)", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op2p_varid ) )
      CALL Check( nf90_put_att( ncid, op2p_varid, "long_name", "Atomic oxygen ion number density  (second excited state)" ) )
      CALL Check( nf90_put_att( ncid, op2p_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "ion_temp", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , ion_temp_varid ) )
      CALL Check( nf90_put_att( ncid, ion_temp_varid, "long_name", "Ion temperature" ) )
      CALL Check( nf90_put_att( ncid, ion_temp_varid, "units", "K" ) )

      CALL Check( nf90_def_var( ncid, "e", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , e_varid ) )
      CALL Check( nf90_put_att( ncid, e_varid, "long_name", "Electron number density" ) )
      CALL Check( nf90_put_att( ncid, e_varid, "units", " m^{-3}" ) )

      CALL Check( nf90_def_var( ncid, "aur_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , ion_rate_varid ) )
      CALL Check( nf90_put_att( ncid, ion_rate_varid, "long_name", "Total Ionization Rate from Auroral Precipitation" ) )
      CALL Check( nf90_put_att( ncid, ion_rate_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "O+_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , O_rate_varid ) )
      CALL Check( nf90_put_att( ncid, O_rate_varid, "long_name", "Oxygen Ionization Rate from Auroral Precipitation" ) )
      CALL Check( nf90_put_att( ncid, O_rate_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "O2+_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , O2_rate_varid ) )
      CALL Check( nf90_put_att( ncid, O2_rate_varid, "long_name", "Molecular Oxygen Ionization Rate from Auroral Precipitation" ) )
      CALL Check( nf90_put_att( ncid, O2_rate_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "N2+_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , N2_rate_varid ) )
      CALL Check( nf90_put_att( ncid, N2_rate_varid, "long_name", "Molecular Nitrogen Ionization Rate from Auroral Precipitation" ) )
      CALL Check( nf90_put_att( ncid, N2_rate_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_enddef(ncid) )
      
      CALL Check( nf90_put_var( ncid, time_varid, time ) )
      IF( ipe % parameters % write_geographic_neutrals )THEN
        CALL Check( nf90_put_var( ncid, x_varid, ipe % grid % longitude_geo ) )
        CALL Check( nf90_put_var( ncid, y_varid, ipe % grid % latitude_geo ) )
        CALL Check( nf90_put_var( ncid, z_varid, ipe % grid % altitude_geo ) )
        CALL Check( nf90_put_var( ncid, helium_varid, ipe % neutrals % geo_helium ) )
        CALL Check( nf90_put_var( ncid, oxygen_varid, ipe % neutrals % geo_oxygen ) )
        CALL Check( nf90_put_var( ncid, molecular_oxygen_varid, ipe % neutrals % geo_molecular_oxygen ) )
        CALL Check( nf90_put_var( ncid, molecular_nitrogen_varid, ipe % neutrals % geo_molecular_nitrogen ) )
        CALL Check( nf90_put_var( ncid, nitrogen_varid, ipe % neutrals % geo_nitrogen ) )
        CALL Check( nf90_put_var( ncid, hydrogen_varid, ipe % neutrals % geo_hydrogen ) )
        CALL Check( nf90_put_var( ncid, temperature_varid, ipe % neutrals % geo_temperature ) )
        CALL Check( nf90_put_var( ncid, u_varid, ipe % neutrals % geo_velocity(1,:,:,:) ) )
        CALL Check( nf90_put_var( ncid, v_varid, ipe % neutrals % geo_velocity(2,:,:,:) ) )
        CALL Check( nf90_put_var( ncid, w_varid, ipe % neutrals % geo_velocity(3,:,:,:) ) )

      ENDIF

      IF( ipe % parameters % write_geographic_eldyn )THEN

        CALL Check( nf90_put_var( ncid, phi_varid, ipe % eldyn % geo_electric_potential ) )
        CALL Check( nf90_put_var( ncid, mhd_phi_varid, ipe % eldyn % geo_mhd_electric_potential ) )
        CALL Check( nf90_put_var( ncid, exbu_varid, ipe % eldyn % geo_v_ExB_geographic(1,:,:) ) )
        CALL Check( nf90_put_var( ncid, exbv_varid, ipe % eldyn % geo_v_ExB_geographic(2,:,:) ) )
        CALL Check( nf90_put_var( ncid, hc_varid, ipe % eldyn % geo_hall_conductivity ) )
        CALL Check( nf90_put_var( ncid, pc_varid, ipe % eldyn % geo_pedersen_conductivity ) )
        CALL Check( nf90_put_var( ncid, bc_varid, ipe % eldyn % geo_b_parallel_conductivity ) )

      ENDIF

      CALL Check( nf90_put_var( ncid, op_varid,  ipe % plasma % geo_ion_densities(1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, hp_varid,  ipe % plasma % geo_ion_densities(2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, hep_varid, ipe % plasma % geo_ion_densities(3,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, np_varid,  ipe % plasma % geo_ion_densities(4,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, nop_varid, ipe % plasma % geo_ion_densities(5,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, o2p_varid, ipe % plasma % geo_ion_densities(6,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, n2p_varid, ipe % plasma % geo_ion_densities(7,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, op2d_varid,  ipe % plasma % geo_ion_densities(8,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, op2p_varid,  ipe % plasma % geo_ion_densities(9,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, ion_temp_varid, ipe % plasma % geo_ion_temperature ) )
      CALL Check( nf90_put_var( ncid, e_varid, ipe % plasma % geo_electron_density ) )
      CALL Check( nf90_put_var( ncid, ion_rate_varid, ipe % plasma % geo_ionization_rates(1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, O_rate_varid, ipe % plasma % geo_ionization_rates(2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, O2_rate_varid, ipe % plasma % geo_ionization_rates(3,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, N2_rate_varid, ipe % plasma % geo_ionization_rates(4,:,:,:) ) )

      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Write_Geographic_NetCDF_IPE

END MODULE IPE_Model_Class
