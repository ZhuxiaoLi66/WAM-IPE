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

      PROCEDURE, PRIVATE :: Empirical_E_Field_Wrapper
!      PROCEDURE, PRIVATE :: Regrid_Empirical_Potential

  END TYPE IPE_Electrodynamics

  LOGICAL, PRIVATE :: setup_efield_empirical

  REAL(prec), PRIVATE :: theta90_rad(0:nmlat)
  REAL(prec), PRIVATE, ALLOCATABLE :: geospace_latitude(:), geospace_longitude(:), geospace_potential(:,:)

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
    INTEGER :: n_lon_geospace, n_lat_geospace


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
 
      PRINT*, MINVAL( geospace_potential )
      PRINT*, MAXVAL( geospace_potential )

    ! read input ( maybe into local variable )

  END SUBROUTINE Read_Geospace_Potential
 
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
!        CALL eldyn % Merge_Potential

!      ELSEIF( openggcm )THEN

!      .....

!      ENDIF

!      CALL eldyn % Get_Efield_From_Potential( )
!      CALL eldyn % Get_ExB( )
        

  END SUBROUTINE Update_IPE_Electrodynamics
 
  SUBROUTINE Empirical_E_Field_Wrapper( eldyn, grid, forcing, time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_Forcing ), INTENT(in)             :: forcing
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    ! Local
    INTEGER :: i, j, year
    REAL(prec) :: mlon90_rad(0:nmlon)
    REAL(prec) :: theta130_rad, utime, lat, ylonm_local(0:nmlon)
     

        IF( setup_efield_empirical )THEN
          CALL efield_init
          setup_efield_empirical = .FALSE.

          ! Maps latitude from 130km to 90km along flux tube.
          DO j=0,nmlat

            theta130_rad   = ( 90.0_prec - ylatm(j) ) * dtr
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
        !CALL eldyn % Regrid_Empirical_Potential( grid )
        
        ! Calculate the potential gradient in IPE coordinates.


  END SUBROUTINE Empirical_E_Field_Wrapper

!  SUBROUTINE Regrid_Empirical_Potential( eldyn, grid )
!    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn 
!    TYPE( IPE_Grid ), INTENT(in)                :: grid
!    ! Local
!    INTEGER :: i, j, ii, jj, i_e, j_e
! 
!      Ed1_90=0.0_prec
!      Ed2_90=0.0_prec
!      theta90_rad = 0.0_prec
!      coslam_m    = 0.0_prec
!
!      DO lp=1,NLP
!
!        coslam_m(1,lp) = COS( pi*0.5_prec - grid % magnetic_colatitude(1,lp) )
!        coslam_m(2,lp) = COS( pi*0.5_prec - grid % magnetic_colatitude(grid % flux_tube_max(lp),lp) )
!
!        IF ( grid % magnetic_colatitude(1,lp) > theta90_rad(nmlat/2) ) THEN
!          j0(1,lp)=nmlat/2+1   !1:NH
!          j1(1,lp)=nmlat/2
!          j0(2,lp)=nmlat/2     !2:SH
!          j1(2,lp)=nmlat/2-1
!
!        ELSE 
!
!          DO j=nmlat,nmlat/2,-1
!
!            IF( theta90_rad(j) <= grid % magnetic_colatitude(1,lp) .AND. &
!                grid % magnetic_colatitude(IN,lp) <= theta90_rad(j-1) ) THEN
!
!              j0(1,lp)=j  !1:NH
!              j1(1,lp)=j-1
!              j0(2,lp)=nmlat-(j-1)  !2:SH
!              j1(2,lp)=nmlat-j
!
!              EXIT
!
!            ENDIF
!
!          ENDDO         
!
!        ENDIF
!
!      ENDDO 
!
!      d_phi_m = dlonm * dtr !constant
!      r = earth_radius + 90000.0_prec ![m]
!
!      DO mp=1,NMP
!        DO i=0,nmlon
!
!          i0=i
!          i1=i+1
!          IF ( i1>nmlon ) i1=i+1-nmlon
!          mlon130_0=mlon130_rad(i0)
!          IF ( mlon130_rad(i0)>mlon130_rad(i1) ) then
!            mlon130_0=mlon130_rad(i0)-pi*2.0  
!          ENDIF
!          IF ( mlon_rad(mp)>=mlon130_0.AND.                             &
!     &         mlon_rad(mp)<=mlon130_rad(i1) ) THEN
!            EXIT mlon130_loop1
!          ELSE
!          END IF
!        END DO mlon130_loop1 !: DO i=0,nmlon
!
!        mlat_loop90km1: DO lp=1,NLP
!          IN = JMIN_IN(lp)
!          IS = JMAX_IS(lp)
!
!!         computipng ed1_90(lp,mp)
!!         FUNC-linearinterpolation does not work! need more debugging. temporary use the average of the two potentials.
!!         linear interpolation of the potent at plasma_grid_3d(IN,mp) in mlat
!!         NH
!!         d          potent_i0 = LINEAR_INTERPOLATION(theta90_rad(j0(1,lp))
!!         d     &   ,potent(i0(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!!         d     &   ,potent(i0(mp),j1(1,lp)),plasma_grid_3d(IN,mp)%GL)
!!         d          potent_i1 = LINEAR_INTERPOLATION(theta90_rad(j0(1,lp))
!!         d     &   ,potent(i1(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!!         d     &   ,potent(i1(mp),j0(1,lp)),plasma_grid_3d(IN,mp)%GL)
!!         ii0=i0
!!         ii1=i1
!
!
!
!
!
!
!          jj0=j0(1,lp)              !1:NH
!          jj1=j1(1,lp)
!          pot_i1=( potent(i1,jj0)+potent(i1,jj1) )*0.50 
!          pot_i0=( potent(i0,jj0)+potent(i0,jj1) )*0.50
!
!          ed1_90(1,lp,mp)=-1.0/r/coslam_m(1,lp)                         &
!     &                         *(pot_i1-pot_i0)/d_phi_m
!!         SH
!!         d          potent_i0=LINEAR_INTERPOLATION(theta90_rad(j0(1,lp)) 
!!         d     &   ,potent(i0(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!!         d     &   ,potent(i0(mp),j1(1,lp)),plasma_grid_3d(IN,mp)%GL)
!!         d          potent_i1=LINEAR_INTERPOLATION(theta90_rad(j0(1,lp))
!!         d     &   ,potent(i1(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!!         d     &   ,potent(i1(mp),j0(1,lp)),plasma_grid_3d(IN,mp)%GL)
!!         ii0=i0
!!         ii1=i1
!          jj0=j0(2,lp)              !2:SH
!          jj1=j1(2,lp)              !2:SH
!          pot_i1=( potent(i1,jj0)+potent(i1,jj1) )*0.50
!          pot_i0=( potent(i0,jj0)+potent(i0,jj1) )*0.50
!          ed1_90(2,lp,mp)=-1.0/r/coslam_m(2,lp)*(pot_i1-pot_i0)/d_phi_m
!!         computing ed2_90(lp,mp) continues
!!         calculate sinIm !eq(3.7)
!          cos2Lambda_m = coslam_m(2,lp) * coslam_m(2,lp)      ! 0<cos2<1
!          sinLambda_m(1)  = + SQRT( 1.0 - cos2Lambda_m )  !>0 ---NH 
!          sinLambda_m(2)  = - SQRT( 1.0 - cos2Lambda_m )  !<0 ---SH 
!          sinI_m(1:2)= 2.0*sinLambda_m(1:2)/SQRT(4.0-3.0*cos2Lambda_m)
!!         NH
!          ihem=1
!          jj0=j0(ihem,lp)              !1:NH
!          jj1=j1(ihem,lp)              !1:NH
!          d_lam_m = theta90_rad( jj1 ) - theta90_rad( jj0 )
!          pot_j1=( potent(i0,jj1)+potent(i1,jj1) )*0.50
!          pot_j0=( potent(i0,jj0)+potent(i1,jj0) )*0.50
!          ed2_90(1,lp,mp)=+1.0/r/sinI_m(ihem)*(pot_j1-pot_j0)/d_lam_m   &
!     &*(-1.) !dbg20140224
!!         dbg20111108     &*(-1.)*sinI_m(ihem)     !E_m_lambda (5.10)
!!         SH
!          ihem=2
!          jj0=j0(ihem,lp)              !2:SH
!          jj1=j1(ihem,lp)              !2:SH
!          d_lam_m = theta90_rad( jj1 ) - theta90_rad( jj0 )
!          pot_j1=( potent(i0,jj1)+potent(i1,jj1) )*0.50
!          pot_j0=( potent(i0,jj0)+potent(i1,jj0) )*0.50
!          ed2_90(2,lp,mp)=+1.0/r/sinI_m(ihem)*(pot_j1-pot_j0) /d_lam_m  &
!     &*(-1.) !dbg20140224
!!         dbg20111108     &*(-1.)*sinI_m(ihem)  !E_m_lambda (5.10)
!
!          IF(sw_exb_up<=1.and.sw_perp_transport>=1) then
!             if(lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans) THEN 
!
!!initialization
!                VEXBup(lp,mp)=zero
!                VEXBe(lp,mp)=zero
!                VEXBth(lp,mp)=zero
!
!!           (0) self consistent electrodynamics comming soon...
!!           (1) WACCM E empirical model
!!           Ed1/2[V/m] at ( phi_t1(mp), theta_t1(lp) ), Be3[T]
!!           note: Ed1_90, Ed2_90, and Be3 are constant along magnetic field lines!!! 
!            midpoint = JMIN_IN(lp) + (JMAX_IS(lp) - JMIN_IN(lp))/2
!
!!nm20130830: Ed1/2_90 should be constant along magnetic field lines!!!
!            v_e(1) =   Ed2_90(1,lp,mp) / Be3(lp,mp) !(4.18) +mag-east(d1?) 
!            v_e(2) = - Ed1_90(1,lp,mp) / Be3(lp,mp) !(4.19) +down/equatorward(d2?)
!
!! EXB in geographic frame at (midpoint,lp,mp)
!            vexbgeo(east )=(v_e(1)*apexE(midpoint,lp,mp,east,1))        &
!     &                    +(v_e(2)*apexE(midpoint,lp,mp,east,2))
!            vexbgeo(north)=(v_e(1)*apexE(midpoint,lp,mp,north,1))       &
!     &                    +(v_e(2)*apexE(midpoint,lp,mp,north,2))
!            vexbgeo(up   )=(v_e(1)*apexE(midpoint,lp,mp,up   ,1))       &
!     &                    +(v_e(2)*apexE(midpoint,lp,mp,up   ,2))
!! EXB  magnetic exact upward at APEX (midpoint) 
!!           VEXBup(lp,mp) = v_e(2) * (-1.0) !convert from down to UPward
!            VEXBup(lp,mp) = vexbgeo(up)
!
!               ipts = JMIN_IN(lp) ! if ihem=1 NH
!!               ipts = JMAX_IS(lp) ! if ihem=2 SH
!! EXB in geographic frame at 90km NH (IN,lp,mp)
!            vexbgeo(east )=(v_e(1)*apexE(ipts,lp,mp,east,1))            &
!     &                    +(v_e(2)*apexE(ipts,lp,mp,east,2))
!            vexbgeo(north)=(v_e(1)*apexE(ipts,lp,mp,north,1))           &
!     &                    +(v_e(2)*apexE(ipts,lp,mp,north,2))
!            vexbgeo(up   )=(v_e(1)*apexE(ipts,lp,mp,up   ,1))           &
!     &                    +(v_e(2)*apexE(ipts,lp,mp,up   ,2))
!! vperp at 90km NH +poleward
!!VEXBth: horizontal component, positive EQUATORward: is calculated only for th-method
!!dbg20150319: it might be better to estimate sinI at JMIN_IN
!!              VEXBth(lp,mp) = v_e(2) * sinI_m(1) !if ihem=1 NH
!               VEXBth(lp,mp) = v_e(2) / sinI_m(1) !if ihem=1 NH
!
!!temporary solution to test the code
!               VEXBe(lp,mp) = 0.0
!
!
!              if ( sw_perp_transport==1 )   VEXBe(lp,mp)=zero
!
!            else                   ! (lp<lpmin_perp_trans.or.lp>lpmax_perp_trans) THEN 
!               VEXBup(lp,mp)=zero
!               VEXBe(lp,mp)=zero
!               VEXBth(lp,mp)=zero 
!            end if              !(lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans) THEN 
!          END IF !( sw_exb_up<=1.and. ... ) 
!
!
!        END DO mlat_loop90km1 !: DO lp=1,NLP
!      END DO mlon_loop90km0     !: DO mp=1,nmp
!!dbg20160408 sms debug 
!!SMS$PARALLEL END
!
!
!         IF ( MOD( (utime_local-start_time),ip_freq_output)==0 ) THEN
!!SMS$SERIAL(<vexbup,vexbe,vexbth,IN>:default=ignore) BEGIN
!            write(unit=2011,FMT='(20E12.4)') vexbup
!            write(unit=2012,FMT='(20E12.4)') vexbe
!            write(unit=2013,FMT='(E12.4)  ') sunlons(1)
!            write(unit=2014,FMT='(20E12.4)') vexbth
!!SMS$SERIAL END
!         END IF                 ! ( MOD( (utime_local-start_time),ip_freq_output)==0 ) THEN
!
!
!  END SUBROUTINE Regrid_Empirical_Potential

  FUNCTION MLT_to_MagneticLongitude( mlt, year, day_of_year, utime, nlon ) RESULT( mag_longitude )
    INTEGER    :: nlon
    REAL(prec) :: mlt(0:nlon)
    INTEGER    :: year, day_of_year
    REAL(prec) :: utime
    REAL(prec) :: mag_longitude(0:nlon)
    ! Local
    INTEGER    :: i
    REAL(prec) :: sunlons

      ! Map magnetic local time to magnetic longitude
      CALL sunloc( year, day_of_year, utime, sunlons ) 
      DO i=0,nmlon

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
