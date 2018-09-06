MODULE IPE_Electrodynamics_Class


USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines
USE IPE_Grid_Class
USE IPE_Forcing_Class
USE IPE_Time_Class

USE efield_ipe

USE netcdf
!
!! For the TIEGCM Wrapper
!USE module_init_cons
!USE cons_module
!USE magfield_module
!USE dynamo_module
!USE module_sub_dynamo
!USE module_highlat
!USE module_init_heelis




IMPLICIT NONE

  TYPE IPE_Electrodynamics
    INTEGER    :: nFluxTube, NLP, NMP
    REAL(prec), ALLOCATABLE :: electric_potential(:,:)   
    REAL(prec), ALLOCATABLE :: mhd_electric_potential(:,:)   
    REAL(prec), ALLOCATABLE :: electric_field(:,:,:) 
    REAL(prec), ALLOCATABLE :: v_ExB_geographic(:,:,:)  ! "zonal" and "meridional" direction on the geographic grid
    REAL(prec), ALLOCATABLE :: v_ExB_apex(:,:,:) ! "zonal" and "meridional" direction ( VEXBth, VEXBe ) on the apex grid

    REAL(prec), ALLOCATABLE, PRIVATE :: lat_interp_weights(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude
    REAL(prec), ALLOCATABLE, PRIVATE :: lon_interp_weights(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude
    INTEGER, ALLOCATABLE, PRIVATE    :: lat_interp_index(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude
    INTEGER, ALLOCATABLE, PRIVATE    :: lon_interp_index(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude

    ! Inputs for the potential solver (on the dynamo grid)
    REAL(prec), ALLOCATABLE :: hall_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: pedersen_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: b_parallel_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: neutral_apex_velocity(:,:,:) ! Components are apex directions


    ! Geographic interpolated attributes
    REAL(prec), ALLOCATABLE :: geo_electric_potential(:,:)
    REAL(prec), ALLOCATABLE :: geo_mhd_electric_potential(:,:)
    REAL(prec), ALLOCATABLE :: geo_v_ExB_geographic(:,:,:)    ! ExB transport velocity with geographic components
    REAL(prec), ALLOCATABLE :: geo_hall_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: geo_pedersen_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: geo_b_parallel_conductivity(:,:)
    

    CONTAINS 

      PROCEDURE :: Build => Build_IPE_Electrodynamics
      PROCEDURE :: Trash => Trash_IPE_Electrodynamics

      PROCEDURE :: Update => Update_IPE_Electrodynamics
      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Electrodynamics

      PROCEDURE :: Read_Geospace_Potential

      PROCEDURE :: Read_MHD_Potential
      PROCEDURE :: Write_MHD_Potential
      PROCEDURE :: Interpolate_Geospace_to_MHDpotential

      PROCEDURE, PRIVATE :: Empirical_E_Field_Wrapper
      PROCEDURE, PRIVATE :: Regrid_Potential
      PROCEDURE, PRIVATE :: Calculate_Potential_Gradient
      PROCEDURE, PRIVATE :: Calculate_ExB_Velocity

  END TYPE IPE_Electrodynamics

  LOGICAL, PRIVATE :: setup_efield_empirical

  REAL(prec), PRIVATE :: theta90_rad(0:nmlat)
  REAL(prec), PRIVATE, ALLOCATABLE :: geospace_latitude(:), geospace_longitude(:), geospace_potential(:,:)

  INTEGER, PRIVATE :: n_lat_geospace, n_lon_geospace

CONTAINS


  SUBROUTINE Build_IPE_Electrodynamics( eldyn, nFluxTube, NLP, NMP ) 
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(out) :: eldyn
    INTEGER, INTENT(in)                       :: nFluxTube
    INTEGER, INTENT(in)                       :: NLP
    INTEGER, INTENT(in)                       :: NMP


      eldyn % nFluxTube = nFluxTube
      eldyn % NLP       = NLP
      eldyn % NMP       = NMP

      ALLOCATE( eldyn % electric_potential(1:NLP,1:NMP), &
                eldyn % mhd_electric_potential(1:NLP,1:NMP), &
                eldyn % electric_field(1:3,1:NLP,1:NMP), &
                eldyn % v_ExB_geographic(1:3,1:NLP,1:NMP), &
                eldyn % v_ExB_apex(1:3,1:NLP,1:NMP), &
                eldyn % hall_conductivity(1:NLP,1:NMP), &
                eldyn % pedersen_conductivity(1:NLP,1:NMP), &
                eldyn % b_parallel_conductivity(1:NLP,1:NMP), &
                eldyn % neutral_apex_velocity(1:3,1:NLP,1:NMP), &
                eldyn % lat_interp_weights(1:2,1:NLP), &
                eldyn % lon_interp_weights(1:2,1:NMP), &
                eldyn % lat_interp_index(1:2,1:NLP), &
                eldyn % lon_interp_index(1:2,1:NMP) )

      eldyn % electric_potential      = 0.0_prec
      eldyn % mhd_electric_potential  = 0.0_prec
      eldyn % electric_field          = 0.0_prec
      eldyn % v_ExB_geographic        = 0.0_prec
      eldyn % v_ExB_apex              = 0.0_prec
      eldyn % hall_conductivity       = 0.0_prec
      eldyn % pedersen_conductivity   = 0.0_prec
      eldyn % b_parallel_conductivity = 0.0_prec
      eldyn % neutral_apex_velocity   = 0.0_prec
      eldyn % lat_interp_weights      = 0.0_prec
      eldyn % lon_interp_weights      = 0.0_prec
      eldyn % lat_interp_index        = 0
      eldyn % lon_interp_index        = 0


      ALLOCATE( eldyn % geo_electric_potential(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_mhd_electric_potential(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_v_ExB_geographic(1:3,1:nlon_geo,1:nlat_geo), &
                eldyn % geo_hall_conductivity(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_pedersen_conductivity(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_b_parallel_conductivity(1:nlon_geo,1:nlat_geo) )

      ! When building the Electrodynamics data structure, we set this
      ! module-private switch to true to ensure that the appropriate 
      ! initialization is executed for the tiegcm model.
      !tie_gcm_init = .TRUE.

      setup_efield_empirical = .TRUE.

  END SUBROUTINE Build_IPE_Electrodynamics

  SUBROUTINE Trash_IPE_Electrodynamics( eldyn ) 
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn

      DEALLOCATE( eldyn % electric_potential, &
                  eldyn % mhd_electric_potential, &
                  eldyn % electric_field, &
                  eldyn % v_ExB_geographic, &
                  eldyn % v_ExB_apex, &
                  eldyn % hall_conductivity, &
                  eldyn % pedersen_conductivity, &
                  eldyn % b_parallel_conductivity, &
                  eldyn % neutral_apex_velocity, &
                  eldyn % lat_interp_weights, &
                  eldyn % lon_interp_weights, & 
                  eldyn % lat_interp_index, & 
                  eldyn % lon_interp_index )

      DEALLOCATE( eldyn % geo_electric_potential, &
                  eldyn % geo_v_ExB_geographic, &
                  eldyn % geo_hall_conductivity, &
                  eldyn % geo_pedersen_conductivity, &
                  eldyn % geo_b_parallel_conductivity )
 
      IF( ALLOCATED( geospace_latitude ) )  DEALLOCATE( geospace_latitude )
      IF( ALLOCATED( geospace_longitude ) ) DEALLOCATE( geospace_longitude )
      IF( ALLOCATED( geospace_potential ) ) DEALLOCATE( geospace_potential )

  END SUBROUTINE Trash_IPE_Electrodynamics

  SUBROUTINE Read_Geospace_Potential( eldyn, filename )
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    CHARACTER(*), INTENT(in)                    :: filename
    ! Local
    INTEGER :: ncid
    INTEGER :: dimid, varid
    INTEGER :: nFluxtube, NLP, NMP
    CHARACTER(NF90_MAX_NAME) :: nameHolder


      CALL Check( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))

      ! Obtain the dimensions of the Geospace grid
      CALL Check( nf90_inq_dimid( ncid, "lon", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, n_lon_geospace ) )

      CALL Check( nf90_inq_dimid( ncid, "lat", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, n_lat_geospace ) )

      IF( .NOT. ALLOCATED( geospace_latitude ) ) ALLOCATE( geospace_latitude(1:n_lat_geospace) )
      IF( .NOT. ALLOCATED( geospace_longitude ) ) ALLOCATE( geospace_longitude(1:n_lat_geospace) )
      IF( .NOT. ALLOCATED( geospace_potential ) ) ALLOCATE( geospace_potential(1:n_lon_geospace,1:n_lat_geospace) )

      CALL Check( nf90_inq_varid( ncid, "lat", varid ) )
      CALL Check( nf90_get_var( ncid, varid, geospace_latitude ) )

      CALL Check( nf90_inq_varid( ncid, "lon", varid ) )
      CALL Check( nf90_get_var( ncid, varid, geospace_longitude ) )

      CALL Check( nf90_inq_varid( ncid, "potential", varid ) )
      CALL Check( nf90_get_var( ncid, varid, geospace_potential ) )

      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Read_Geospace_Potential

   SUBROUTINE Interpolate_Geospace_to_MHDpotential( eldyn, grid, time_tracker)
     IMPLICIT NONE
     CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
     TYPE( IPE_Grid ), INTENT(in)             :: grid
     TYPE( IPE_Time ), INTENT(in)             :: time_tracker
     ! Local
     INTEGER :: j,latidx
     REAL(prec) :: theta110_rad,geospace_latitude_90_rad(1:n_lat_geospace)
     REAL(prec) :: colat_local(1:n_lat_geospace)
     REAL(prec) :: potential_local(1:n_lon_geospace,1:n_lat_geospace)
    
      DO j = 1, n_lat_geospace
        theta110_rad   = ( 90.0_prec - geospace_latitude(j) ) * dtr
        geospace_latitude_90_rad(j) = ASIN(SIN(theta110_rad)*SQRT((earth_radius+90000.0_prec)/(earth_radius+110000.0_prec)))
        IF ( theta110_rad > pi*0.50_prec ) geospace_latitude_90_rad(j) = pi - geospace_latitude_90_rad(j)
      ENDDO
 
      DO j = 1, n_lat_geospace
        latidx = n_lat_geospace-j+1
        colat_local(j)= geospace_latitude_90_rad(latidx)*rtd
        potential_local(:,j)=geospace_potential(:,latidx)
      END DO

      CALL eldyn % Regrid_Potential( grid, time_tracker, potential_local, geospace_longitude, colat_local, 1, n_lon_geospace, n_lat_geospace )

      eldyn % mhd_electric_potential= eldyn % electric_potential
 
   END SUBROUTINE Interpolate_Geospace_to_MHDpotential
 
  SUBROUTINE Write_MHD_Potential( eldyn, grid, time_tracker, filename )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(in) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)             :: grid
    TYPE( IPE_Time ), INTENT(in)             :: time_tracker
    CHARACTER(*), INTENT(in)                 :: filename
    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, time_dimid, time_varid
    INTEGER :: mhd_phi_varid
    INTEGER :: recStart(1:3), recCount(1:3)


      recStart = (/ 1, 1, 1 /)
      recCount = (/ grid % NLP, grid % NMP, 1 /)


      time = time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))

      CALL Check( nf90_def_dim( ncid, "lp", grid % NLP, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "mp", grid % NMP, y_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

      CALL Check( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
      CALL Check( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "units", "minutes" ) )

      CALL Check( nf90_def_var( ncid, "mhd_phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) ,mhd_phi_varid ) )
      CALL Check( nf90_put_att( ncid, mhd_phi_varid, "long_name", "Electric Potential - MHD Component" ) )
      CALL Check( nf90_put_att( ncid, mhd_phi_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_enddef(ncid) )
      
      CALL Check( nf90_put_var( ncid, time_varid, time ) )
      CALL Check( nf90_put_var( ncid, mhd_phi_varid, eldyn % mhd_electric_potential ) )

      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Write_MHD_Potential

  SUBROUTINE Read_MHD_Potential( eldyn, filename )
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    CHARACTER(*), INTENT(in)                    :: filename
    ! Local
    INTEGER :: ncid
    INTEGER :: dimid, varid
    CHARACTER(NF90_MAX_NAME) :: nameHolder


      CALL Check( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))

      CALL Check( nf90_inq_varid( ncid, "mhd_phi", varid ) )
      CALL Check( nf90_get_var( ncid, varid, eldyn % mhd_electric_potential) )

      CALL Check( nf90_close( ncid ) )


  END SUBROUTINE Read_MHD_Potential

  SUBROUTINE Update_IPE_Electrodynamics( eldyn, grid, forcing, time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_Forcing ), INTENT(in)             :: forcing
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    ! Local
    INTEGER :: lp, mp

      CALL eldyn % Empirical_E_Field_Wrapper( grid, forcing, time_tracker )


!      IF( geospace )THEN

!        CALL  eldyn % Read_Geospace_Potential( filename )
!        CALL  eldyn % Regrid_Geospace(  )
!#ifdef DEBUG
!        CALL eldyn % Write_MHD_Potential( )
!#endif
!        CALL eldyn % Merge_Geospace_Potential

!      ELSEIF( openggcm )THEN

!        CALL  eldyn % Read_OpenGGCM_Potential( filename )
!        CALL  eldyn % Regrid_OpenGGCM(  )
!#ifdef DEBUG
!        CALL eldyn % Write_MHD_Potential( )
!#endif
!        CALL eldyn % Merge_OpenGGCM_Potential

!      ENDIF


      ! Calculate the potential gradient in IPE coordinates.
      CALL eldyn % Calculate_Potential_Gradient( grid )

      ! Calculate ExB drift velocity
      CALL eldyn % Calculate_ExB_Velocity( grid ) 


  END SUBROUTINE Update_IPE_Electrodynamics
 
  SUBROUTINE Empirical_E_Field_Wrapper( eldyn, grid, forcing, time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_Forcing ), INTENT(in)             :: forcing
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    ! Local
    INTEGER :: i, j, year
    REAL(prec) :: theta130_rad, utime, lat
    REAL(prec) :: potent_local(0:nmlon,0:nmlat)
    REAL(prec) :: mlt_local(0:nmlon)
    REAL(prec) :: colat_local(0:nmlat)
     

       IF( setup_efield_empirical )THEN
         CALL efield_init
         setup_efield_empirical = .FALSE.

         ! Maps latitude from 130km to 90km along flux tube.
         DO j=0,nmlat

           theta130_rad   = ( 180.0_prec - ylatm(j) ) * dtr
           theta90_rad(j) = ASIN(SIN( theta130_rad )*SQRT((earth_radius+90000.0_prec)/(earth_radius+130000.0_prec)))
           IF ( theta130_rad > pi*0.50_prec ) theta90_rad(j) = pi-theta90_rad(j)

         END DO 

       ENDIF

       ! iday, iday_m, ut, f107d, bt, angle, v_sw, bz are all variables
       ! declared in efield_ipe.f
       CALL time_tracker % Get_Date( year, imo, iday_m )
       iday = time_tracker % day_of_year                 

       CALL time_tracker % Get_UTime( utime )
       ut=utime/3600.0_prec

       f107d = forcing % f107( forcing % current_index )

       bt = forcing % solarwind_Bt( forcing % current_index )
       angle = forcing % solarwind_angle( forcing % current_index )
       v_sw  = forcing % solarwind_velocity( forcing % current_index )
       bz    = forcing % solarwind_Bz( forcing % current_index )
        
       CALL  get_efield

       ! Interpolate the potential to the IPE grid
       potent_local = potent
       mlt_local    = ylonm
       colat_local  = theta90_rad*180.0_prec/pi
       CALL eldyn % Regrid_Potential( grid, time_tracker, potent_local, mlt_local, colat_local, 0, nmlon, nmlat )
        

  END SUBROUTINE Empirical_E_Field_Wrapper

  SUBROUTINE Calculate_ExB_Velocity( eldyn, grid )
    ! Calculates ExB drift velocity according to Eqs(4.18-4.19) in
    ! A.D. Richmond (1995)
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn 
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    ! Local
    INTEGER    :: mp, lp

      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP
  
          eldyn % v_ExB_apex(1,lp,mp) = eldyn % electric_field(2,lp,mp)/grid % apex_be3(lp,mp) 
          eldyn % v_ExB_apex(2,lp,mp) = -eldyn % electric_field(1,lp,mp)/grid % apex_be3(lp,mp)

        ENDDO
      ENDDO
    
  END SUBROUTINE Calculate_ExB_Velocity

  SUBROUTINE Calculate_Potential_Gradient( eldyn, grid )
    ! Uses 2nd order centered differencing (on the apex grid) to calculate the
    ! electric field components from the potential attribute from the
    ! IPE_Electrodynamics class.
    ! The electric field gradient is calculated using Sections 3 and 4 of A.D.
    ! Richmond (1995)
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn 
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    ! Local
    INTEGER :: mp, lp
    REAL(prec) :: r, d_lp, d_mp, coslam, sinim


      r = earth_radius + 90000.0_prec

      ! lp-component of e-field : Eq(4.8)
      DO mp = 1, grid % NMP

        lp = 1
        coslam = cos( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )
        d_lp = r*( grid % magnetic_colatitude(1,lp+1) - grid % magnetic_colatitude(1,lp) )*coslam

        eldyn % electric_field(1,lp,mp) = -( eldyn % electric_potential(lp+1,mp) - &  
                                            eldyn % electric_potential(lp,mp) )/d_lp

        DO lp = 2, grid % NLP-1

          coslam = cos( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )
          d_lp = r*( grid % magnetic_colatitude(1,lp+1) - grid % magnetic_colatitude(1,lp-1) )*coslam

          eldyn % electric_field(1,lp,mp) = -( eldyn % electric_potential(lp+1,mp) - &  
                                              eldyn % electric_potential(lp-1,mp) )/d_lp

        ENDDO

        lp = grid % NLP
        coslam = cos( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )
        d_lp = r*( grid % magnetic_colatitude(1,lp) - grid % magnetic_colatitude(1,lp-1) )*coslam

        eldyn % electric_field(1,lp,mp) = -( eldyn % electric_potential(lp,mp) - &  
                                            eldyn % electric_potential(lp-1,mp) )/d_lp

      ENDDO

    
      ! mp-component of e-field
      ! Do periodic boundary conditions ( @ mp == 1 )
      mp = 1
      DO lp = 1, grid % NLP

        coslam = cos( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )
        sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
        d_mp = sinim*r*( grid % magnetic_longitude(mp+1) - grid % magnetic_longitude( grid % NMP ) + 2.0_prec*pi )

        eldyn % electric_field(2,lp,mp) = ( eldyn % electric_potential(lp,mp+1) - &  
                                            eldyn % electric_potential(lp,grid % NMP) )/d_mp
      ENDDO
      
      DO mp = 2, grid % NMP-1
        DO lp = 1, grid % NLP

          coslam = cos( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )
          sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
          d_mp = sinim*r*( grid % magnetic_longitude(mp+1) - grid % magnetic_longitude(mp-1) )

          eldyn % electric_field(2,lp,mp) = ( eldyn % electric_potential(lp,mp+1) - &  
                                              eldyn % electric_potential(lp,mp-1) )/d_mp

        ENDDO
      ENDDO

      ! Do periodic boundary conditions ( @ mp == grid % NMP )
      mp = grid % NMP
      DO lp = 1, grid % NLP

        coslam = cos( 0.5_prec*pi - grid % magnetic_colatitude(1,lp) )
        sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
        d_mp = sinim*r*( grid % magnetic_longitude(1) - grid % magnetic_longitude(mp-1) + 2.0_prec*pi )

        eldyn % electric_field(2,lp,mp) = ( eldyn % electric_potential(lp,1) - &  
                                            eldyn % electric_potential(lp,mp-1) )/d_mp
      ENDDO

  END SUBROUTINE Calculate_Potential_Gradient

  SUBROUTINE Regrid_Potential( eldyn, grid, time_tracker, potential, mlt, colat, start_index, nlon, nlat )
  ! This subroutine regrids electric potential from a structured grid
  ! on (magnetic local time, magnetic colatitude) to IPE's grid.
  !
  ! Input :
  !
  !   grid
  !
  !   time_tracker
  !
  !   potential
  !
  !   mlt - magnetic local time [ deg ]
  !
  !   colat - magnetic colatitude [ deg ]
  !
  !   start_index - starting index for the potential, mlt, and colat arrays
  !
  !   nlon - last index in the longitude direction
  !
  !   nlat - last index in the latitude direction
  !
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn 
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    INTEGER, INTENT(in)                         :: start_index, nlon, nlat
    REAL(prec), INTENT(in)                      :: potential(start_index:nlon,start_index:nlat)
    REAL(prec), INTENT(in)                      :: mlt(start_index:nlon)
    REAL(prec), INTENT(in)                      :: colat(start_index:nlat)
    ! Local
    INTEGER :: mp, lp, i, j, i1, i2, j1, j2
    INTEGER :: ilat(1:grid % NLP), jlon(1:grid % NMP)
    REAL(prec) :: lat, lon, dlon, difflon, mindiff, wsum
    REAL(prec) :: lat_weight(1:2), lon_weight(1:2)
    REAL(prec) :: mlon90_rad(start_index:nlon)
    REAL(prec) :: mlat90_rad(start_index:nlon)
 

      mlon90_rad = MLT_to_MagneticLongitude( mlt, 1999, time_tracker % day_of_year, time_tracker % utime, start_index, nlon )
      mlat90_rad = colat*pi/180.0_prec
      

      ! Search for nearest grid points in the magnetic longitude/latitude grid
      DO mp = 1, grid % NMP
         lon = grid % magnetic_longitude(mp)        

         mindiff = 1000.0_prec ! This is just an initial guess for the minimum search algorithm
                               ! It is purposefully a large number that will
                               ! likely not be the minimum
         DO j = start_index, nlon

            difflon = ABS(mlon90_rad(j) - lon)
            IF( difflon < mindiff )THEN
              jlon(mp) = j
              mindiff = difflon
            ENDIF

         ENDDO


      ENDDO

      DO lp = 1, grid % NLP

        lat = grid % magnetic_colatitude(1,lp)        
        ! colatitude decreases with increasing lp
        ilat(lp) = nlat
        DO i = start_index, nlat

          ! Need to pass in colatitude through the call stack (theta90_rad ->
          ! colatitude )
          IF( mlat90_rad(i) < lat )THEN
            ilat(lp) = i
            EXIT
          ENDIF

        ENDDO

      ENDDO

      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP

          lat = grid % magnetic_colatitude(1,lp)        
          lon = grid % magnetic_longitude(mp)        

          IF( ilat(lp) == start_index )THEN

            i1 = ilat(lp)
            i2 = ilat(lp)
            lat_weight(1) = 1.0_prec
            lat_weight(2) = 0.0_prec 

          ELSE

            i1 = ilat(lp)-1
            i2 = ilat(lp)
            lat_weight(1) =  ( lat - mlat90_rad(i2) )/( mlat90_rad(i1) - mlat90_rad(i2) )
            lat_weight(2) = -( lat - mlat90_rad(i1) )/( mlat90_rad(i1) - mlat90_rad(i2) )

          ENDIF



         IF( jlon(mp) == start_index )THEN

           IF( mlon90_rad(jlon(mp)) > lon )THEN

             j1 = nlon
             j2 = 0

             dlon = mlon90_rad(j1) - mlon90_rad(j2)
             IF( dlon > 0.0_prec )THEN
               mlon90_rad(j1) = mlon90_rad(j1) - 2.0_prec*pi
               dlon = mlon90_rad(j1) - mlon90_rad(j2)
             ENDIF
 
             lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
             lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

           ELSE

             j1 = 0
             j2 = 1

             dlon = mlon90_rad(j1) - mlon90_rad(j2)
             IF( dlon > 0.0_prec )THEN
               mlon90_rad(j1) = mlon90_rad(j1) - 2.0_prec*pi
               dlon = mlon90_rad(j1) - mlon90_rad(j2)
             ENDIF
 
             lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
             lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

            ENDIF

         ELSEIF( jlon(mp) == nlon )THEN

           IF( mlon90_rad(jlon(mp)) > lon )THEN

             j1 = nlon-1
             j2 = nlon

             dlon = mlon90_rad(j1) - mlon90_rad(j2)
             IF( dlon > 0.0_prec )THEN
               mlon90_rad(j1) = mlon90_rad(j1) - 2.0_prec*pi
               dlon = mlon90_rad(j1) - mlon90_rad(j2)
             ENDIF
 
             lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
             lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

           ELSE

             j1 = nlon
             j2 = 0

             dlon = mlon90_rad(j1) - mlon90_rad(j2)
             IF( dlon > 0.0_prec )THEN
               mlon90_rad(j1) = mlon90_rad(j1) - 2.0_prec*pi
               dlon = mlon90_rad(j1) - mlon90_rad(j2)
             ENDIF
 
             lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
             lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

            ENDIF

         ELSE

           IF( mlon90_rad(jlon(mp)) > lon )THEN

             j1 = jlon(mp)-1
             j2 = jlon(mp)

             dlon = mlon90_rad(j1) - mlon90_rad(j2)
             IF( dlon > 0.0_prec )THEN
               mlon90_rad(j1) = mlon90_rad(j1) - 2.0_prec*pi
               dlon = mlon90_rad(j1) - mlon90_rad(j2)
             ENDIF
 
             lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
             lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

           ELSE

             j1 = jlon(mp)
             j2 = jlon(mp)+1

             dlon = mlon90_rad(j1) - mlon90_rad(j2)
             IF( dlon > 0.0_prec )THEN
               mlon90_rad(j1) = mlon90_rad(j1) - 2.0_prec*pi
               dlon = mlon90_rad(j1) - mlon90_rad(j2)
             ENDIF
 
             lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon 
             lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

            ENDIF

         ENDIF


         eldyn % electric_potential(lp,mp) = ( potential(j1,i1)*lat_weight(1)*lon_weight(1) +&
                                               potential(j2,i1)*lat_weight(1)*lon_weight(2) +&
                                               potential(j1,i2)*lat_weight(2)*lon_weight(1) +&
                                               potential(j2,i2)*lat_weight(2)*lon_weight(2) )


       ENDDO
     ENDDO

         
         
  END SUBROUTINE Regrid_Potential

  FUNCTION MLT_to_MagneticLongitude( mlt, year, day_of_year, utime, start_index, nlon ) RESULT( mag_longitude )
    INTEGER    :: start_index, nlon
    REAL(prec) :: mlt(start_index:nlon)
    INTEGER    :: year, day_of_year
    REAL(prec) :: utime
    REAL(prec) :: mag_longitude(start_index:nlon)
    ! Local
    INTEGER    :: i
    REAL(prec) :: sunlons

      ! Map magnetic local time to magnetic longitude
      CALL sunloc( year, day_of_year, utime, sunlons ) 
      DO i=start_index,nlon

        mag_longitude(i)=(mlt(i)-180.0_prec)*pi/180.0_prec+sunlons
        IF( mag_longitude(i) < 0.0_prec   ) mag_longitude(i)=mag_longitude(i)+pi*2.0
        IF( mag_longitude(i) >= pi*2.0_prec ) mag_longitude(i)=mag_longitude(i)-pi*2.0

      END DO

  END FUNCTION MLT_to_MagneticLongitude


  SUBROUTINE sunloc(iyr,iday,secs,sunlons)
    integer, intent(in) :: iyr, iday ! day of year
    REAL(prec), intent(in)   ::  secs    ! ut in seconds
    REAL(prec), INTENT(out) :: sunlons
    integer :: ihr,imn
    real :: sec,date,vp,xmlon ! apex magnetic longitude
    real ::  sbsllat    ! geographic latitude of subsolar point (degrees)
    real ::  sbsllon    ! geographic longitude of subsolar point (degrees)
    real ::  colat      ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    real ::  elon        ! East longitude of geomagnetic dipole north pole (deg)
    
    ihr = int(secs/3600.)
    imn = int((secs - float(ihr)*3600.)/60.)
    sec = secs - float(ihr)*3600. - float(imn)*60.
    
    !  calculate subsol point: given universal time
    !          input: iyr,iday,ihr,imn,sec
    !          output: sbsllat,sbsllon 
    !                  
    call subsol(iyr,iday,ihr,imn,sec ,sbsllat,sbsllon)
    
    date = float(iyr) + float(iday)/365. + float(ihr)/24./365. + &
           float(imn)/60./24./365.+ sec/60./60./24./365.
  
    call cofrm(date)
    call dypol(colat,elon,vp)
    
    ! calculate geomagn. diploe longitude
    !        input: aloni,sbsllat,sbsllon,colat,elon
    !        output: xmlon  
    call solgmlon(sbsllat,sbsllon,colat,elon,xmlon) 
    sunlons = xmlon*dtr

  END SUBROUTINE sunloc

  SUBROUTINE Interpolate_to_GeographicGrid_IPE_Electrodynamics( eldyn, grid )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)         :: grid

      CALL grid % Interpolate_2D_to_Geographic_Grid( eldyn % electric_potential, eldyn % geo_electric_potential )
      CALL grid % Interpolate_2D_to_Geographic_Grid( eldyn % mhd_electric_potential, eldyn % geo_mhd_electric_potential )

  END SUBROUTINE Interpolate_to_GeographicGrid_IPE_Electrodynamics

!!  SUBROUTINE TIEGCM_Wrapper( eldyn )
!!
!!
!!      IF( tie_gcm_init )THEN
!!
!!        CALL init_cons             
!!        CALL init_heelis           
!!
!!        tie_gcm_init = .FALSE.
!!
!!      ENDIF
!!      
!!      CALL sunloc(jyr,jday,jsecs)   ! every timestep
!!
!!! read in integrals      
!!!      if(input_type == 'NETCDF') then
!!!        call readin_netcdf
!!!      elseif(input_type == 'ASCII') then
!!!        call readin_ascii
!!!      else
!!!        write(6,*) 'Did not recognize input_type= ',
!!!      ENDif
!!       
!!
!!      call highlat
!!
!!      call dynamo
!!
!!  END SUBROUTINE TIEGCM_Wrapper
!
!! SUBROUTINE interface_field_line_integrals ( lp_plas,mp,utime,  &
!!                  sigma_ped_3d,sigma_hall_3d,Ue1_3d,Ue2_3d,Ne_3d )
!!
!!  USE module_precision
!!  USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_grid_3d,apexD,
!!JMIN_IN,JMAX_IS,Be3, apexDscalar,MaxFluxTube, ISL,IBM,east,north,up, &
!!                                       TN_k,ON_m3,O2N_m3,N2N_m3,Un_ms1,
!!plasma_3d
!!  USE module_input_PARAMETERs,ONLY:sw_3DJ,mype
!!  USE module_IPE_dimension,ONLY: NPTS2D,NMP
!!  USE module_calculate_field_line_integrals,ONLY:calculate_field_line_integrals
!!  USE module_save2fli_array,ONLY:save2fli_array
!!  USE cons_module               ,ONLY: lp_dyn_eq
!!  IMPLICIT NONE
!!  !---argument
!!  INTEGER (KIND=int_prec), INTENT(in)  ::  mp
!!  INTEGER (KIND=int_prec), INTENT(in)  ::  lp_plas
!!  INTEGER (KIND=int_prec), INTENT(IN)  ::  utime !universal time [sec]
!!  !---local
!!  INTEGER (KIND=int_prec)  ::  in,in1d    ! North footpoint of a flux tube
!!  INTEGER (KIND=int_prec)  ::  is,is1d    ! South footpoint of a flux tube
!!  INTEGER (KIND=int_prec)  ::  midpoint,midpoint1d
!!  REAL (KIND=REAL(prec)_prec)    ::  ds(MaxFluxTube)            !ds 1D [???
!!
!!  !nm20130830: Be3 is constant along a field line!
!!  REAL (KIND=REAL(prec)_prec)  ::  Apex_BE3 !B_e3 of reference above (=
!!Bmag/D),[Tesla] (4.13)
!!  REAL (KIND=REAL(prec)_prec)  ::  Apex_D(MaxFluxTube)
!!  REAL (KIND=REAL(prec)_prec)  ::  Apex_d1d1(MaxFluxTube)
!!  REAL (KIND=REAL(prec)_prec)  ::  Apex_d1d2(MaxFluxTube)
!!  REAL (KIND=REAL(prec)_prec)  ::  Apex_d2d2(MaxFluxTube)
!!  REAL (KIND=REAL(prec)_prec)  ::  Apex_BMAG(MaxFluxTube)
!!
!!  REAL (KIND=REAL(prec)_prec)  ::  ni_oplus_1d(MaxFluxTube)  !O+ densities [m-6]
!!  REAL (KIND=REAL(prec)_prec)  ::  ni_hplus_1d(MaxFluxTube)  !H+ densities [m-6]
!!  REAL (KIND=REAL(prec)_prec)  ::  ti_oplus_1d(MaxFluxTube)  !O+ temperatures [K]
!!  REAL (KIND=REAL(prec)_prec)  ::  no_plus(MaxFluxTube)  !NO+ density [m-6]
!!  REAL (KIND=REAL(prec)_prec)  ::  o2_plus(MaxFluxTube)  !O2+ density [m-6]
!!  REAL (KIND=REAL(prec)_prec)  ::  n2_plus(MaxFluxTube)  !N2+ density [m-6]
!!  REAL (KIND=REAL(prec)_prec)  ::  n_plus(MaxFluxTube)  !N+ density [m-6]
!!  REAL (KIND=REAL(prec)_prec)  ::  tn(MaxFluxTube)         !Tn [K]
!!  REAL (KIND=REAL(prec)_prec)  ::  o(MaxFluxTube)          !atomic oxygen density []
!!  REAL (KIND=REAL(prec)_prec)  ::  o2(MaxFluxTube)         !molecular oxygen density
!![]
!!  REAL (KIND=REAL(prec)_prec)  ::  n2(MaxFluxTube)         !molecular nitrogen density
!![]
!!
!!  REAL (KIND=REAL(prec)_prec)  ::  Ue(MaxFluxTube,2)        !Un parallel to
!!e1(1)/e2(2) [m/s]  (5.6)
!!  !nm20130830:  REAL (KIND=REAL(prec)_prec)  ::  Ue2(MaxFluxTube)        !Un parallel
!!  !to e2 [m/s]
!!
!!  ! Output Arguments: 
!!
!!  REAL (KIND=REAL(prec)_prec) ::  sigma_ped_1d(MaxFluxTube)  !pedersen conductivity
!![mho/m]
!!  REAL (KIND=REAL(prec)_prec) ::  sigma_hall_1d(MaxFluxTube) !hall conductivity
!!  REAL (KIND=REAL(prec)_prec) ::  sigma_ped_3d(MaxFluxTube,47,NMP)  !pedersen
!!conductivity [mho/m]
!!  REAL (KIND=REAL(prec)_prec) ::  sigma_hall_3d(MaxFluxTube,47,NMP)  !hall
!!conductivity [mho/m]
!!  REAL (KIND=REAL(prec)_prec) ::  Ue1_3d(MaxFluxTube,47,NMP) 
!!  REAL (KIND=REAL(prec)_prec) ::  Ue2_3d(MaxFluxTube,47,NMP) 
!!  REAL (KIND=REAL(prec)_prec) ::  electron_density_out(MaxFluxTube) !hall conductivity
!!  REAL (KIND=REAL(prec)_prec) ::  Ne_3d(MaxFluxTube,47,NMP) 
!!
!!  REAL (KIND=REAL(prec)_prec) ::    sigma_phph_dsi_1d(2)      !(5.13) divided by |sin
!!I_m |
!!  REAL (KIND=REAL(prec)_prec) ::    sigma_lmlm_msi_1d(2)      !(5.14) multiplied by |
!!sin I_m |
!!  REAL (KIND=REAL(prec)_prec) ::          sigma_h_1d(2)      !(5.17)
!!  REAL (KIND=REAL(prec)_prec) ::          sigma_c_1d(2)      !(5.18)
!!  REAL (KIND=REAL(prec)_prec) ::         Kdmph_dsi_1d(2)      !(5.19) divided by |sin
!!I_m |
!!  REAL (KIND=REAL(prec)_prec) ::          Kdmlm_1d(2)   !(5.20) plus or minus ????
!!! when sw_3DJ=1
!!  REAL (KIND=REAL(prec)_prec) ::       Je_1d(MaxFluxTube,2) !1D slice of the
!!3DJ(mp,lp)
!!!------other local variables----------------------
!!  INTEGER (KIND=int_prec) ::  CTIPDIM
!!  INTEGER (KIND=int_prec) ::  NPTS
!!  INTEGER (KIND=int_prec) ::  i,i1d,jth,jth1,jth2,i_set,lp_dyn
!!  REAL (KIND=REAL(prec)_prec) :: dotprod,d1xd2(3)
!!!---
!!  REAL (KIND=REAL(prec)_prec) ::  sigma_ped(NPTS2D,NMP)  !pedersen conductivity
!![mho/m]
!!  REAL (KIND=REAL(prec)_prec) ::  sigma_hall(NPTS2D,NMP) !hall conductivity
!!
!!  lp_dyn_loop: DO lp_dyn=1,lp_dyn_eq !from SH toward eq
!!    IN = JMIN_IN(lp_plas)
!!    IS = JMAX_IS(lp_plas)
!!    midpoint = IN + ( IS - IN )/2
!!    CTIPDIM = IS - IN + 1
!!
!!    i_loop: DO i=in,is
!!      i1d=i-in+1
!!      IF ( i==is ) THEN
!!        ds(i1d) = plasma_grid_3d(IS,lp_plas,mp,ISL) -
!!plasma_grid_3d(IS-1,lp_plas,mp,ISL)
!!      ELSE
!!        ds(i1d) = plasma_grid_3d(i+1,lp_plas,mp,ISL) -
!!plasma_grid_3d(i,lp_plas,mp,ISL)
!!      END IF
!!      !DOT_PRODUCT( D3(1:3,i,mp), Vn_ms1(1:3,i) ) / SQRT(  DOT_PRODUCT(
!!      !D3(1:3,i,mp), D3(1:3,i,mp) )  )     
!!      ! Un(1)=Ue1=d1*U: positive east, Eq(5.6) 
!!      ! Un(2)=Ue2=d2*U: positive down/equatorward, Eq(5.6) 
!!      ! un(3)=Ue3=d3*U: positive parallel to a field line, Eq(5.6) 
!!      DO jth=1,2
!!        Ue(i1d,jth) = Un_ms1(i,lp_plas,mp,jth)
!!      END DO
!!      DO jth=1,3
!!        IF (jth<=2) THEN
!!          jth1=jth
!!          jth2=jth
!!        ELSE IF (jth==3) THEN
!!          jth1=1
!!          jth2=2
!!        END IF
!!        dotprod = apexD(i,lp_plas,mp,east, jth1) * apexD(i,lp_plas,mp,east,
!!jth2)  &
!!             &  + apexD(i,lp_plas,mp,north,jth1) *
!!apexD(i,lp_plas,mp,north,jth2)  &
!!             &  + apexD(i,lp_plas,mp,up,   jth1) * apexD(i,lp_plas,mp,up
!!,jth2)
!!        IF ( jth==1 ) THEN
!!          apex_d1d1(i1d) = dotprod
!!        ELSE IF ( jth==2 ) THEN 
!!          apex_d2d2(i1d) = dotprod
!!        ELSE IF ( jth==3 ) THEN 
!!          apex_d1d2(i1d) = dotprod
!!        END IF
!!      END DO
!!
!!      apex_D(i1d) = apexDscalar(i,lp_plas,mp)
!!      apex_BMAG(i1d) = plasma_grid_3d(i,lp_plas,mp,IBM)
!!
!!      ni_oplus_1d(i1d)=plasma_3d(i,lp_plas,mp,1)
!!      ni_hplus_1d(i1d)=plasma_3d(i,lp_plas,mp,2)
!!            
!!      n_plus( i1d)=plasma_3d(i,lp_plas,mp,4)
!!      no_plus(i1d)=plasma_3d(i,lp_plas,mp,5)
!!      o2_plus(i1d)=plasma_3d(i,lp_plas,mp,6)
!!      n2_plus(i1d)=plasma_3d(i,lp_plas,mp,7)
!!      ti_oplus_1d(i1d)=plasma_3d(i,lp_plas,mp,11)
!!            
!!      tn(i1d)=TN_k(  i,lp_plas,mp) ![K]
!!      o(i1d)=ON_m3( i,lp_plas,mp) !m-3
!!      o2(i1d)=O2N_m3(i,lp_plas,mp)
!!      n2(i1d)=N2N_m3(i,lp_plas,mp)
!!            
!!    END DO i_loop
!!
!!    !nm20130830: Be3 is constant along magnetic field lines!
!!    Apex_BE3      =Be3(lp_plas,mp) 
!!    !nm20130830:  Apex_BE3(CTIPDIM)=Be3(2,mp,lp_plas) !SH
!!    in1d = JMIN_IN(lp_plas) - JMIN_IN(lp_plas) + 1
!!    is1d = JMAX_IS(lp_plas) - JMIN_IN(lp_plas) + 1
!!    midpoint1d = ( JMAX_IS(lp_plas) - JMIN_IN(lp_plas) )/2 + 1
!!    NPTS = MaxFluxTube
!!         
!!    CALL calculate_field_line_integrals (NPTS,in1d,is1d,ds,midpoint1d, &
!!      !nm20111107                                     Apex_d1,Apex_d2, &
!!                                                             Apex_BE3, &
!!                       apex_D,apex_d1d1,apex_d1d2,apex_d2d2,apex_BMAG, &
!!                  ni_oplus_1d,ni_hplus_1d,ti_oplus_1d,no_plus,o2_plus, &
!!                                                       n2_plus,n_plus, &
!!                                                           tn,o,o2,n2, &
!!      !nm20111107                              U_zonal,U_merid,U_vert, &
!!                                                                   Ue, &
!!                                                 electron_density_out, &
!!                                           sigma_ped_1d,sigma_hall_1d, & !output
!!                                  sigma_phph_dsi_1d,sigma_lmlm_msi_1d, &
!!                                                sigma_h_1d,sigma_c_1d, &
!!                                                Kdmph_dsi_1d,Kdmlm_1d, &
!!                                                      Je_1d,mp,lp_plas )
!!    ! save to 3D array                                         
!!    call save2fli_array (lp_plas,mp, &
!!                         sigma_phph_dsi_1d,sigma_lmlm_msi_1d, &
!!                         sigma_h_1d,sigma_c_1d, &
!!                         Kdmph_dsi_1d,Kdmlm_1d )
!!    ! Tzu-Wei:save conductivities and Ue
!!    do i_set=in1d,is1d
!!      Ue1_3d(i_set,lp_dyn,mp) = Ue(i_set,1)
!!      Ue2_3d(i_set,lp_dyn,mp) = Ue(i_set,2)
!!      sigma_ped_3d(i_set,lp_dyn,mp) = sigma_ped_1d(i_set)
!!      sigma_hall_3d(i_set,lp_dyn,mp) = sigma_hall_1d(i_set)
!!      Ne_3d(i_set,lp_dyn,mp) = electron_density_out(i_set)
!!    ENDdo
!!  END DO lp_dyn_loop
!!
!!  END SUBROUTINE Interface_Field_Line_Integrals
END MODULE IPE_Electrodynamics_Class

!-----------------------------------------------------------------------
!      SUBROUTINE GET_EFIELD90km ( utime_local )
!      USE module_precision
!      USE efield_ipe,ONLY: nmlat,ylatm,dlonm,potent,ylonm,nmlon
!      USE module_physical_constants,ONLY:pi,rtd,dtr,earth_radius,zero
!      USE module_eldyn,ONLY:j0,j1,theta90_rad,ed1_90,ed2_90,coslam_m
!      USE module_IPE_dimension,ONLY: NMP,NLP
!      USE module_FIELD_LINE_GRID_MKS
!      USE module_input_parameters
!      USE module_magfield,ONLY:sunlons
!      USE module_sunloc,ONLY:sunloc
!      IMPLICIT NONE
!      INTEGER (KIND=int_prec),INTENT(IN)   :: utime_local !universal time [sec]
!      REAL(KIND=real_prec) :: utsecs
!      INTEGER :: iyr
!      REAL(KIND=real_prec) :: theta130_rad
!      REAL(KIND=real_prec),PARAMETER :: ht130 = 130.0E+03     !height in meter
!      INTEGER (KIND=int_prec) :: j,i,mp,lp
!      INTEGER (KIND=int_prec) :: IN,IS
!      REAL(KIND=real_prec) :: mlon90_deg !deg
!      REAL(KIND=real_prec) :: d_phi_m, d_lam_m !in radian
!      REAL(KIND=real_prec) :: r !in meter
!      INTEGER (KIND=int_prec) :: jj0,jj1
!      INTEGER(KIND=int_prec) :: i0,i1,ihem
!      REAL(KIND=real_prec) :: pot_i0,pot_i1, pot_j0,pot_j1
!      REAL(KIND=real_prec),DIMENSION(0:nmlon) :: mlon130_rad
!      REAL(KIND=real_prec) :: mlon130_0
!      REAL(KIND=real_prec) :: cos2Lambda_m,sinLambda_m(2),sinI_m(2)
!      INTEGER (KIND=int_prec ) :: midpoint, ipts
!      REAL    (KIND=real_prec) :: v_e(2)   !1:ed2/be3 (4.18) ;2: -ed1/be3 (4.19)
!      REAL    (KIND=real_prec) :: vexbgeo(east:up) !EXB in gegraphic frame
!
!
!      END SUBROUTINE GET_EFIELD90km
!!
!!20110919: not used for now,so commented out. needs more debugging
!!     FUNCTION LINEAR_INTERPOLATION(X0,Y0,X1,Y1,XX)
!!     USE module_precision
!!     REAL(KIND=real_prec),INTENT(IN) :: X0,Y0  !(X0,Y0)
!!     REAL(KIND=real_prec),INTENT(IN) :: X1,Y1  !(X1,Y1)
!!     REAL(KIND=real_prec),INTENT(IN) :: XX     !X coordinate of the intepolated value YY
!!     REAL(KIND=real_prec) :: LINEAR_INTERPOLATION !the intepolated value YY at (XX)
!!     LINEAR_INTERPOLATION = Y0 + (XX - X0) * (Y1 - Y0) / (X1 - X0) 
!!     END  FUNCTION LINEAR_INTERPOLATION
!
!      SUBROUTINE get_theta1_at_ht1 (ht0,theta0, ht1,theta1)
!
!      USE module_precision
!      USE module_physical_constants,ONLY:earth_radius,pi
!      IMPLICIT NONE
!      REAL(KIND=real_prec),INTENT(IN) :: ht0 !m / km
!      REAL(KIND=real_prec),INTENT(IN) :: theta0 ![rad]
!      REAL(KIND=real_prec),INTENT(IN) :: ht1 !m / km
!      REAL(KIND=real_prec),INTENT(OUT) :: theta1
!!     ---local
!      REAL(KIND=real_prec) :: sintheta0
!      REAL(KIND=real_prec) :: sintheta1
!!     find out grid point at ht1=90km (theta1,phi1) of theta0 using dipole assumption
!      sintheta0 = SIN( theta0 )
!      sintheta1 = sintheta0*SQRT((earth_radius+ht1)/(earth_radius+ht0))
!      theta1    = ASIN(sintheta1)
!!     SH:example: mlat=30 comlat=90-30=60, mlatSH=-30,comlatSH=90+(90-60)=180-60
!      IF ( theta0 > pi*0.50 ) theta1 = pi-theta1
!
!      END SUBROUTINE get_theta1_at_ht1
