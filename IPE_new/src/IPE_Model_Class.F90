MODULE IPE_Model_Class

USE IPE_Precision
USE IPE_Model_Parameters_Class
USE IPE_Grid_Class
USE IPE_Neutrals_Class
USE IPE_Forcing_Class
!USE IPE_Plasma_Class
!USE IPE_Electrodynamics_Class

USE netcdf

IMPLICIT NONE


  ! The IPE_Model serves as a wrapper for all of the underlying attributes.
  ! This class should be used to orchestrate model setup, updating, and
  ! breakdown. At this level, we define the API for interacting with the
  ! deeper attributes within IPE.

  TYPE IPE_Model

    REAL(prec)                   :: simulation_time
    TYPE( IPE_Model_Parameters ) :: parameters
    TYPE( IPE_Grid )             :: grid
    TYPE( IPE_Forcing )          :: forcing
    TYPE( IPE_Neutrals )         :: neutrals
    !TYPE( IPE_Plasma )           :: plasma
    !TYPE( IPE_Eldyn )            :: eldyn

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model
      PROCEDURE :: Trash => Trash_IPE_Model

      PROCEDURE :: Geographic_Interpolation
      PROCEDURE :: Write_NetCDF_IPE
!      PROCEDURE :: Read_NetCDF_IPE

      PROCEDURE :: Write_Geographic_NetCDF_IPE
   
  END TYPE IPE_Model


  REAL(prec), PARAMETER :: fillValue = -999999.9_prec

CONTAINS

  SUBROUTINE Build_IPE_Model( ipe, init_success )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(out) :: ipe
    LOGICAL, INTENT(out)            :: init_success
    ! Local 

      CALL ipe % parameters % Build( init_success )

      IF( init_success )THEN

        CALL ipe % forcing % Build( ipe % parameters % solar_forcing_time_step, &
                                    ipe % parameters % f107_kp_size )
                                    
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

        init_success = .FALSE.
        RETURN
        
      ENDIF

  END SUBROUTINE Build_IPE_Model
!
  SUBROUTINE Trash_IPE_Model( ipe )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe

      CALL ipe % forcing % Trash( )
      CALL ipe % grid % Trash( )
      CALL ipe % neutrals % Trash( )

  END SUBROUTINE Trash_IPE_Model
!
  SUBROUTINE Update_IPE_Model( ipe, t0, t1 )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe
    REAL(prec), INTENT(in)            :: t0, t1
    ! Local
    INTEGER    :: i, im3hr, im6hr, im9hr, im12hr, im36hr
    REAL(prec) :: AP(1:7)

    
      ! This call to update the neutrals is "diagnostic" and returns neutral
      ! parameter at the time t0
      i       = ipe % forcing % current_index
      im3hr   = MAX( 1, ipe % forcing % current_index-INT(3.0_prec*3600.0_prec/ipe % forcing % dt ) )
      im6hr   = MAX( 1, ipe % forcing % current_index-INT(6.0_prec*3600.0_prec/ipe % forcing % dt ) )
      im9hr   = MAX( 1, ipe % forcing % current_index-INT(9.0_prec*3600.0_prec/ipe % forcing % dt ) )
      im12hr  = MAX( 1, ipe % forcing % current_index-INT(12.0_prec*3600.0_prec/ipe % forcing % dt ) )
      im36hr  = MAX( 1, ipe % forcing % current_index-INT(36.0_prec*3600.0_prec/ipe % forcing % dt ) )

      !AP(1) = apa_eld(lpi)
      !AP(2) =  ap_eld(lpi)
      !AP(3) =  ap_eld(lpi-INT( 3.*60*60/input_params_interval))
      !AP(4) =  ap_eld(lpi-INT( 6.*60*60/input_params_interval))
      !AP(5) =  ap_eld(lpi-INT( 9.*60*60/input_params_interval))
      !AP(6) = apa_eld(lpi-INT(12.*60*60/input_params_interval))
      !AP(7) = apa_eld(lpi-INT(36.*60*60/input_params_interval))

      AP(1) = ipe % forcing % ap_1day_avg(i)
      AP(2) = ipe % forcing % ap(i)
      AP(3) = ipe % forcing % ap(im3hr)
      AP(4) = ipe % forcing % ap(im6hr)
      AP(5) = ipe % forcing % ap(im9hr)
      AP(6) = ipe % forcing % ap_1day_avg(im12hr)
      AP(7) = ipe % forcing % ap_1day_avg(im36hr)

      CALL ipe % neutrals % Update( ipe % grid, &
                                    ipe % simulation_time, &
                                    2000, 76,  & ! Hard code year 2000 and day 76 for now
                                    ipe % forcing % f107(i) , &
                                    ipe % forcing % f107_81day_avg(i), &
                                    AP )


  END SUBROUTINE Update_IPE_Model
!
  SUBROUTINE Geographic_Interpolation( ipe )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(inout) :: ipe


     IF( ipe % parameters % write_geographic_neutrals )THEN
       CALL ipe % neutrals % Interpolate_to_GeographicGrid( ipe % grid )
     ENDIF

  END SUBROUTINE Geographic_Interpolation

  SUBROUTINE Write_NetCDF_IPE( ipe, filename )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(in) :: ipe
    CHARACTER(*), INTENT(in)       :: filename
    ! Local
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, z_dimid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid


      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))

      CALL Check( nf90_def_dim( ncid, "s", ipe % grid % nFluxTube, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "lp", ipe % grid % NLP, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "mp", ipe % grid % NMP, y_dimid ) ) 

      IF( ipe % parameters % write_apex_neutrals )THEN

        CALL Check( nf90_def_var( ncid, "helium", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , helium_varid ) )
        CALL Check( nf90_put_att( ncid, helium_varid, "long_name", "Neutral Helium Density" ) )
        CALL Check( nf90_put_att( ncid, helium_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "oxygen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , oxygen_varid ) )
        CALL Check( nf90_put_att( ncid, oxygen_varid, "long_name", "Neutral Oxygen Density" ) )
        CALL Check( nf90_put_att( ncid, oxygen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "molecular_oxygen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , molecular_oxygen_varid ) )
        CALL Check( nf90_put_att( ncid, molecular_oxygen_varid, "long_name", "Neutral Molecular Oxygen Density" ) )
        CALL Check( nf90_put_att( ncid, molecular_oxygen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "molecular_nitrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , molecular_nitrogen_varid ) )
        CALL Check( nf90_put_att( ncid, molecular_nitrogen_varid, "long_name", "Neutral Molecular Nitrogen Density" ) )
        CALL Check( nf90_put_att( ncid, molecular_nitrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "nitrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , nitrogen_varid ) )
        CALL Check( nf90_put_att( ncid, nitrogen_varid, "long_name", "Neutral Nitrogen Density" ) )
        CALL Check( nf90_put_att( ncid, nitrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "hydrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , hydrogen_varid ) )
        CALL Check( nf90_put_att( ncid, hydrogen_varid, "long_name", "Neutral Hydrogen Density" ) )
        CALL Check( nf90_put_att( ncid, hydrogen_varid, "units", "kg m^{-3}" ) )

        CALL Check( nf90_def_var( ncid, "temperature", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , temperature_varid ) )
        CALL Check( nf90_put_att( ncid, temperature_varid, "long_name", "Thermosphere Temperature" ) )
        CALL Check( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

        CALL Check( nf90_def_var( ncid, "u", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , u_varid ) )
        CALL Check( nf90_put_att( ncid, u_varid, "long_name", "Zonal Velocity" ) )
        CALL Check( nf90_put_att( ncid, u_varid, "units", "m s^{-1}" ) )

        CALL Check( nf90_def_var( ncid, "v", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , v_varid ) )
        CALL Check( nf90_put_att( ncid, v_varid, "long_name", "Meridional Velocity" ) )
        CALL Check( nf90_put_att( ncid, v_varid, "units", "m s^{-1}" ) )

        CALL Check( nf90_def_var( ncid, "w", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , w_varid ) )
        CALL Check( nf90_put_att( ncid, w_varid, "long_name", "Radial Velocity" ) )
        CALL Check( nf90_put_att( ncid, w_varid, "units", "m s^{-1}" ) )

      ENDIF

      CALL Check( nf90_enddef(ncid) )
      
      IF( ipe % parameters % write_apex_neutrals )THEN

        CALL Check( nf90_put_var( ncid, helium_varid, ipe % neutrals % helium ) )
        CALL Check( nf90_put_var( ncid, oxygen_varid, ipe % neutrals % oxygen ) )
        CALL Check( nf90_put_var( ncid, molecular_oxygen_varid, ipe % neutrals % molecular_oxygen ) )
        CALL Check( nf90_put_var( ncid, molecular_nitrogen_varid, ipe % neutrals % molecular_nitrogen ) )
        CALL Check( nf90_put_var( ncid, nitrogen_varid, ipe % neutrals % nitrogen ) )
        CALL Check( nf90_put_var( ncid, hydrogen_varid, ipe % neutrals % hydrogen ) )
        CALL Check( nf90_put_var( ncid, temperature_varid, ipe % neutrals % temperature ) )
        CALL Check( nf90_put_var( ncid, u_varid, ipe % neutrals % velocity_geographic(1,:,:,:) ) )
        CALL Check( nf90_put_var( ncid, v_varid, ipe % neutrals % velocity_geographic(2,:,:,:) ) )
        CALL Check( nf90_put_var( ncid, w_varid, ipe % neutrals % velocity_geographic(3,:,:,:) ) )

      ENDIF



      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Write_NetCDF_IPE
!
!  SUBROUTINE Read_IPE_Model_NetCDF( grid, filename )
!    IMPLICIT NONE
!    CLASS( IPE_Grid ), INTENT(inout) :: grid
!    CHARACTER(*), INTENT(in)         :: filename
!    ! Local
!    INTEGER :: NF90_PREC
!    INTEGER :: ncid
!    INTEGER :: dimid, varid
!    INTEGER :: nFluxtube, NLP, NMP
!    CHARACTER(NF90_MAX_NAME) :: nameHolder
!
!
!      
!      IF( prec == sp )THEN
!        NF90_PREC = NF90_FLOAT
!      ELSE      
!        NF90_PREC = NF90_DOUBLE
!      ENDIF
!
!      CALL Check( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))
!  END SUBROUTINE Read_IPE_Model_NetCDF
!
  SUBROUTINE Write_Geographic_NetCDF_IPE( ipe, filename )
    IMPLICIT NONE
    CLASS( IPE_Model ), INTENT(in) :: ipe
    CHARACTER(*), INTENT(in)       :: filename
    ! Local
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid


      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))

      CALL Check( nf90_def_dim( ncid, "z", nheights_geo, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

      CALL Check( nf90_def_var( ncid, "z", NF90_PREC, z_dimid, z_varid ) )
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
      CALL Check( nf90_put_att( ncid, time_varid, "long_name", "Internal IPE Time" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "units", "s" ) )
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

      CALL Check( nf90_enddef(ncid) )
      
      IF( ipe % parameters % write_geographic_neutrals )THEN

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

      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Write_Geographic_NetCDF_IPE

END MODULE IPE_Model_Class
