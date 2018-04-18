MODULE IPE_Grid_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines

USE netcdf

IMPLICIT NONE

  TYPE IPE_Grid
    INTEGER :: nFluxTube, NLP, NMP, NPTS2D

    REAL(prec), ALLOCATABLE :: altitude(:,:)
    REAL(prec), ALLOCATABLE :: latitude(:,:,:)
    REAL(prec), ALLOCATABLE :: longitude(:,:,:)
    REAL(prec), ALLOCATABLE :: grx(:,:,:)    ! GR index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: foot_point_distance(:,:,:)     ! ISL index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: magnetic_field_strength(:,:,:) ! IBM index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: magnetic_colatitude(:,:) ! plasma_grid_GL
    REAL(prec), ALLOCATABLE :: r_meter(:,:) ! rmeter2D

    REAL(prec), ALLOCATABLE :: q_factor(:,:,:) ! Q index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: l_magnitude(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: apex_d_vectors(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: d1xd2_magnitude(:,:,:)
    REAL(prec), ALLOCATABLE :: apex_e_vectors(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: apex_be3(:,:)

    INTEGER, ALLOCATABLE :: flux_tube_midpoint(:)
    INTEGER, ALLOCATABLE :: flux_tube_max(:)
    INTEGER, ALLOCATABLE :: southern_top_index(:)
    INTEGER, ALLOCATABLE :: northern_top_index(:)

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Grid
      PROCEDURE :: Trash => Trash_IPE_Grid

      ! -- Legacy I/O -- !
      PROCEDURE :: Read_IPE_Grid
      PROCEDURE :: Calculate_Grid_Attributes_From_Legacy_Input
      ! ---------------- !
 
      PROCEDURE :: Write_IPE_Grid_NetCDF
      PROCEDURE :: Read_IPE_Grid_NetCDF

      ! Functions
      PROCEDURE :: SinI
      
  END TYPE IPE_Grid 

CONTAINS

  SUBROUTINE Build_IPE_Grid( grid, nFluxTube, NLP, NMP, NPTS2D )
    IMPLICIT NONE
    CLASS( IPE_Grid ), INTENT(out) :: grid
    INTEGER, INTENT(in)            :: nFluxTube
    INTEGER, INTENT(in)            :: NLP
    INTEGER, INTENT(in)            :: NMP
    INTEGER, INTENT(in)            :: NPTS2D

      grid % nFluxTube = nFluxTube
      grid % NLP       = NLP
      grid % NMP       = NMP
      grid % NPTS2D    = NPTS2D
   
      ALLOCATE( grid % altitude(1:nFluxTube,1:NLP), &
                grid % latitude(1:nFluxTube,1:NLP,1:NMP), &
                grid % longitude(1:nFluxTube,1:NLP,1:NMP), &
                grid % grx(1:nFluxTube,1:NLP,1:NMP), &
                grid % foot_point_distance(1:nFluxTube,1:NLP,1:NMP), &
                grid % magnetic_field_strength(1:nFluxTube,1:NLP,1:NMP), &
                grid % magnetic_colatitude(1:nFluxTube,1:NLP), &
                grid % r_meter(1:nFluxTube,1:NLP), &
                grid % q_factor(1:nFluxTube,1:NLP,1:NMP), &
                grid % l_magnitude(1:3,1:2,1:nFluxTube,1:NLP,1:NMP), &
                grid % apex_e_vectors(1:3,1:2,1:nFluxTube,1:NLP,1:NMP), &
                grid % apex_d_vectors(1:3,1:3,1:nFluxTube,1:NLP,1:NMP), &
                grid % d1xd2_magnitude(1:nFluxTube,1:NLP,1:NMP), &
                grid % apex_be3(1:NLP,1:NMP), &
                grid % flux_tube_midpoint(1:NLP), &
                grid % flux_tube_max(1:NLP), &
                grid % southern_top_index(1:NLP), &
                grid % northern_top_index(1:NLP) )

      grid % altitude                = 0.0_prec
      grid % latitude                = 0.0_prec
      grid % longitude               = 0.0_prec
      grid % grx    = 0.0_prec
      grid % foot_point_distance     = 0.0_prec
      grid % magnetic_field_strength = 0.0_prec
      grid % magnetic_colatitude     = 0.0_prec
      grid % r_meter                 = 0.0_prec
      grid % q_factor                = 0.0_prec
      grid % l_magnitude             = 0.0_prec
      grid % apex_e_vectors          = 0.0_prec
      grid % apex_d_vectors          = 0.0_prec
      grid % d1xd2_magnitude         = 0.0_prec
      grid % apex_be3                = 0.0_prec

      grid % flux_tube_midpoint = 0
      grid % flux_tube_max      = 0
      grid % southern_top_index = 0
      grid % northern_top_index = 0


  END SUBROUTINE Build_IPE_Grid
!
  SUBROUTINE Trash_IPE_Grid( grid )
    IMPLICIT NONE
    CLASS( IPE_Grid ), INTENT(inout) :: grid
   
      DEALLOCATE( grid % altitude, &
                  grid % latitude, &
                  grid % longitude, &
                  grid % grx, &
                  grid % foot_point_distance, &
                  grid % magnetic_field_strength, &
                  grid % magnetic_colatitude, &
                  grid % r_meter, &
                  grid % q_factor, &
                  grid % l_magnitude, &
                  grid % apex_e_vectors, &
                  grid % apex_d_vectors, &
                  grid % d1xd2_magnitude, &
                  grid % apex_be3, &
                  grid % flux_tube_midpoint, &
                  grid % flux_tube_max, &
                  grid % southern_top_index, &
                  grid % northern_top_index )

  END SUBROUTINE Trash_IPE_Grid
!
  SUBROUTINE Read_IPE_Grid( grid, filename )
    IMPLICIT NONE
    CLASS( IPE_Grid ), INTENT(inout) :: grid
    CHARACTER(*), INTENT(in)         :: filename
    ! Local
    INTEGER    :: i, lp, mp
    INTEGER    :: fUnit, ii, NMP, NLP, NPTS2D, MaxFluxTube
    INTEGER    :: jmin(1:grid % NMP,1:grid % NLP)
    INTEGER    :: jmax(1:grid % NMP,1:grid % NLP)
    REAL(prec) ::  Be3_all1(1:grid % NMP, 1:grid % NLP)
    REAL(prec) ::  Be3_all2(1:grid % NMP, 1:grid % NLP)
    REAL(prec), ALLOCATABLE ::  dum0(:,:) 
    REAL(prec), ALLOCATABLE ::  dum1(:,:)
    REAL(prec), ALLOCATABLE ::  dum2(:,:)
    REAL(prec), ALLOCATABLE ::  dum3(:,:)
    REAL(prec), ALLOCATABLE ::  dum4(:,:,:)
    REAL(prec), ALLOCATABLE ::  dum5(:,:,:)
    REAL(prec), ALLOCATABLE ::  dum6(:,:,:)

 
      ! Place all of these dummies on the heap so that we don't wind up
      ! with a stack overflow. For typical values of NPTS2D, NLP, NMP,
      ! the "dum*" arrays will easily take up 100's of MB.

      ALLOCATE( dum0(1:grid % NPTS2D,1:grid % NMP), & 
                dum1(1:grid % NPTS2D,1:grid % NMP), &
                dum2(1:grid % NPTS2D,1:grid % NMP), &
                dum3(1:grid % NPTS2D,1:grid % NMP), &
                dum4(1:3,1:grid % NPTS2D,1:grid % NMP), &
                dum5(1:3,1:grid % NPTS2D,1:grid % NMP), &
                dum6(1:3,1:grid % NPTS2D,1:grid % NMP) )

      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = TRIM(filename), &
            STATUS = 'OLD', &
            FORM   = 'FORMATTED', &
            ACTION = 'READ' ) 

      READ( fUnit, * ) jmin, jmax

      DO lp = 1, grid % NLP

        grid % flux_tube_max(lp) = jmax(1,lp) - jmin(1,lp) + 1
         
      ENDDO

      MaxFluxTube = maxval(jmax(1,:)-jmin(1,:)) + 1

      ! Just in case we set up the grid initially with the improper number of
      ! flux tube points, we'll reset the grid here.
      IF( MaxFluxTube /= grid % nFluxTube )THEN

        NLP = grid % NLP
        NMP = grid % NMP
        NPTS2D = grid % NPTS2D
        
        CALL grid % Trash( )
        CALL grid % Build( MaxFluxTube, NLP, NMP, NPTS2D )
   
      ENDIF

     ! grid % flux_tube_max = grid % flux_tube_max - jmin(1,:) + 1
      DO lp = 1, grid % NLP

        grid % flux_tube_midpoint(lp) = 1 + ( grid % flux_tube_max(lp) - 1)/2

      ENDDO

      READ( fUnit, * ) dum0, dum1, dum2, dum3

      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          ii = jmin(1,lp) + (i-1)
          grid % r_meter(i,lp) =  dum0(ii,1)
          grid % altitude(i,lp) = grid % r_meter(i,lp) - earth_radius 
       
        ENDDO
      ENDDO     

      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP

          DO i = 1, grid % flux_tube_max(lp)

            ii = jmin(1,lp) + (i-1)
            grid % latitude(i,lp,mp)  = dum1(ii,mp)
            grid % longitude(i,lp,mp) = dum2(ii,mp)
            grid % q_factor(i,lp,mp)  = dum3(ii,mp)
    
          ENDDO

        ENDDO
      ENDDO


      READ( fUnit, * ) dum0

      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          ii = jmin(1,lp) + (i-1)
          grid % magnetic_colatitude(i,lp) = dum0(ii,1)
    
        ENDDO
      ENDDO

      READ( fUnit, * ) dum0, dum1
      
      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)
          
            ii = jmin(1,lp) + (i-1)
            grid % foot_point_distance(i,lp,mp)     = dum0(ii,mp)
            grid % magnetic_field_strength(i,lp,mp) = dum1(ii,mp)

          ENDDO
        ENDDO
      ENDDO


      READ( fUnit, * ) dum4, dum5, dum6

      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)
          
            ii = jmin(1,lp) + (i-1)
            grid % apex_d_vectors(1,1,i,lp,mp) = dum4(1,ii,mp)
            grid % apex_d_vectors(2,1,i,lp,mp) = dum4(2,ii,mp)
            grid % apex_d_vectors(3,1,i,lp,mp) = dum4(3,ii,mp)

            grid % apex_d_vectors(1,2,i,lp,mp) = dum5(1,ii,mp)
            grid % apex_d_vectors(2,2,i,lp,mp) = dum5(2,ii,mp)
            grid % apex_d_vectors(3,2,i,lp,mp) = dum5(3,ii,mp)

            grid % apex_d_vectors(1,3,i,lp,mp) = dum6(1,ii,mp)
            grid % apex_d_vectors(2,3,i,lp,mp) = dum6(2,ii,mp)
            grid % apex_d_vectors(3,3,i,lp,mp) = dum6(3,ii,mp)

          ENDDO
        ENDDO
      ENDDO

     
      READ( fUnit, * ) dum4, dum5

      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)
          
            ii = jmin(1,lp) + (i-1)
            grid % apex_e_vectors(1,1,i,lp,mp) = dum4(1,ii,mp)
            grid % apex_e_vectors(2,1,i,lp,mp) = dum4(2,ii,mp)
            grid % apex_e_vectors(3,1,i,lp,mp) = dum4(3,ii,mp)

            grid % apex_e_vectors(1,2,i,lp,mp) = dum5(1,ii,mp)
            grid % apex_e_vectors(2,2,i,lp,mp) = dum5(2,ii,mp)
            grid % apex_e_vectors(3,2,i,lp,mp) = dum5(3,ii,mp)

          ENDDO
        ENDDO
      ENDDO

      READ( fUnit, * ) Be3_all1, Be3_all2

      DO mp=1, grid % NMP
        DO lp=1, grid % NLP

          grid % apex_be3(lp,mp) = Be3_all1(mp,lp)

        ENDDO
      ENDDO

      CLOSE( fUnit )

      DEALLOCATE( dum0,&
                  dum1,&
                  dum2,&
                  dum3,&
                  dum4,&
                  dum5,&
                  dum6 )

  END SUBROUTINE Read_IPE_Grid
!
  SUBROUTINE Write_IPE_Grid_NetCDF( grid, filename )
    IMPLICIT NONE
    CLASS( IPE_Grid ), INTENT(in) :: grid
    CHARACTER(*), INTENT(in)      :: filename
    ! Local
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, z_dimid
    INTEGER :: altitude_varid, latitude_varid, longitude_varid, gravity_varid
    INTEGER :: fpd_varid, B_varid, colat_varid, r_varid, q_varid
    INTEGER :: l11_varid, l21_varid, l31_varid
    INTEGER :: l12_varid, l22_varid, l32_varid
    INTEGER :: apexe11_varid, apexe21_varid, apexe31_varid
    INTEGER :: apexe12_varid, apexe22_varid, apexe32_varid
    INTEGER :: apexd11_varid, apexd21_varid, apexd31_varid
    INTEGER :: apexd12_varid, apexd22_varid, apexd32_varid
    INTEGER :: apexd13_varid, apexd23_varid, apexd33_varid, d1xd2_varid
    INTEGER :: be3_varid, midpoint_varid, max_varid, south_varid, north_varid


      
      IF( prec == sp )THEN
        NF90_PREC = NF90_FLOAT
      ELSE      
        NF90_PREC = NF90_DOUBLE
      ENDIF

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))

      CALL Check( nf90_def_dim( ncid, "s", grid % nFluxTube, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "lp", grid % NLP, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "mp", grid % NMP, y_dimid ) ) 

      CALL Check( nf90_def_var( ncid, "altitude", NF90_PREC, (/ z_dimid, x_dimid /) , altitude_varid ) )
      CALL Check( nf90_put_att( ncid, altitude_varid, "long_name", "Radial distance above spherical earth" ) )
      CALL Check( nf90_put_att( ncid, altitude_varid, "units", "m" ) )

      CALL Check( nf90_def_var( ncid, "latitude", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , latitude_varid ) )
      CALL Check( nf90_put_att( ncid, latitude_varid, "long_name", "Geographic Latitude" ) )
      CALL Check( nf90_put_att( ncid, latitude_varid, "units", "radians" ) )

      CALL Check( nf90_def_var( ncid, "longitude", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , longitude_varid ) )
      CALL Check( nf90_put_att( ncid, longitude_varid, "long_name", "Geographic Longitude" ) )
      CALL Check( nf90_put_att( ncid, longitude_varid, "units", "radians" ) )

      CALL Check( nf90_def_var( ncid, "GRX", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , gravity_varid ) )
      CALL Check( nf90_put_att( ncid, gravity_varid, "long_name", "Gravity (used in FLIP)" ) )
      CALL Check( nf90_put_att( ncid, gravity_varid, "units", "m^{2} s^{-1}" ) )

      CALL Check( nf90_def_var( ncid, "foot_point_distance", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , fpd_varid ) )
      CALL Check( nf90_put_att( ncid, fpd_varid, "long_name", "Distance along flux tube" ) )
      CALL Check( nf90_put_att( ncid, fpd_varid, "units", "m" ) )

      CALL Check( nf90_def_var( ncid, "B", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , B_varid ) )
      CALL Check( nf90_put_att( ncid, B_varid, "long_name", "Magnetic Field Strength" ) )
      CALL Check( nf90_put_att( ncid, B_varid, "units", "Micro-Tesla" ) )

      CALL Check( nf90_def_var( ncid, "m_colat", NF90_PREC, (/ z_dimid, x_dimid /) , colat_varid ) )
      CALL Check( nf90_put_att( ncid, colat_varid, "long_name", "Magnetic Colatitude" ) )
      CALL Check( nf90_put_att( ncid, colat_varid, "units", "Radians" ) )

      CALL Check( nf90_def_var( ncid, "r_meter", NF90_PREC, (/ z_dimid, x_dimid /) , r_varid ) )
      CALL Check( nf90_put_att( ncid, r_varid, "long_name", "Radial distance from earth center" ) )
      CALL Check( nf90_put_att( ncid, r_varid, "units", "meters" ) )

      CALL Check( nf90_def_var( ncid, "q_factor", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , q_varid ) )
      CALL Check( nf90_put_att( ncid, q_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, q_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "l_magnitude_11", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , l11_varid ) )
      CALL Check( nf90_put_att( ncid, l11_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, l11_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "l_magnitude_21", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , l21_varid ) )
      CALL Check( nf90_put_att( ncid, l21_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, l21_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "l_magnitude_31", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , l31_varid ) )
      CALL Check( nf90_put_att( ncid, l31_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, l31_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "l_magnitude_12", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , l12_varid ) )
      CALL Check( nf90_put_att( ncid, l12_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, l12_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "l_magnitude_22", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , l22_varid ) )
      CALL Check( nf90_put_att( ncid, l22_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, l22_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "l_magnitude_32", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , l32_varid ) )
      CALL Check( nf90_put_att( ncid, l32_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, l32_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_e_11", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexe11_varid ) )
      CALL Check( nf90_put_att( ncid, apexe11_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexe11_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_e_21", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexe21_varid ) )
      CALL Check( nf90_put_att( ncid, apexe21_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexe21_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_e_31", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexe31_varid ) )
      CALL Check( nf90_put_att( ncid, apexe31_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexe31_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_e_12", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexe12_varid ) )
      CALL Check( nf90_put_att( ncid, apexe12_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexe12_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_e_22", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexe22_varid ) )
      CALL Check( nf90_put_att( ncid, apexe22_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexe22_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_e_32", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexe32_varid ) )
      CALL Check( nf90_put_att( ncid, apexe32_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexe32_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_11", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd11_varid ) )
      CALL Check( nf90_put_att( ncid, apexd11_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd11_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_21", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd21_varid ) )
      CALL Check( nf90_put_att( ncid, apexd21_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd21_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_31", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd31_varid ) )
      CALL Check( nf90_put_att( ncid, apexd31_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd31_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_12", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd12_varid ) )
      CALL Check( nf90_put_att( ncid, apexd12_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd12_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_22", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd22_varid ) )
      CALL Check( nf90_put_att( ncid, apexd22_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd22_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_32", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd32_varid ) )
      CALL Check( nf90_put_att( ncid, apexd32_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd32_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_13", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd13_varid ) )
      CALL Check( nf90_put_att( ncid, apexd13_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd13_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_23", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd23_varid ) )
      CALL Check( nf90_put_att( ncid, apexd23_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd23_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_d_33", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , apexd33_varid ) )
      CALL Check( nf90_put_att( ncid, apexd33_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, apexd33_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "d1xd2_mag", NF90_PREC, (/ z_dimid, x_dimid, y_dimid /) , d1xd2_varid ) )
      CALL Check( nf90_put_att( ncid, d1xd2_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, d1xd2_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "apex_be3", NF90_PREC, (/ x_dimid, y_dimid /) , be3_varid ) )
      CALL Check( nf90_put_att( ncid, be3_varid, "long_name", "Unknown" ) )
      CALL Check( nf90_put_att( ncid, be3_varid, "units", "[Unknown]" ) )

      CALL Check( nf90_def_var( ncid, "tube_midpoint", NF90_INT, (/ x_dimid /) , midpoint_varid ) )
      CALL Check( nf90_put_att( ncid, midpoint_varid, "long_name", "Flux Tube Midpoints" ) )
      CALL Check( nf90_put_att( ncid, midpoint_varid, "units", "Index" ) )

      CALL Check( nf90_def_var( ncid, "tube_max", NF90_INT, (/ x_dimid /) , max_varid ) )
      CALL Check( nf90_put_att( ncid, max_varid, "long_name", "Index of maximum flux tube height" ) )
      CALL Check( nf90_put_att( ncid, max_varid, "units", "Index" ) )

      CALL Check( nf90_def_var( ncid, "southern_top", NF90_INT, (/ x_dimid /) , south_varid ) )
      CALL Check( nf90_put_att( ncid, south_varid, "long_name", "Index of southern hemisphere top of flux tube" ) )
      CALL Check( nf90_put_att( ncid, south_varid, "units", "Index" ) )

      CALL Check( nf90_def_var( ncid, "northern_top", NF90_INT, (/ x_dimid /) , north_varid ) )
      CALL Check( nf90_put_att( ncid, north_varid, "long_name", "Index of southern hemisphere top of flux tube" ) )
      CALL Check( nf90_put_att( ncid, north_varid, "units", "Index" ) )

      CALL Check( nf90_enddef(ncid) )
      
      CALL Check( nf90_put_var( ncid, altitude_varid, grid % altitude ) )
      CALL Check( nf90_put_var( ncid, latitude_varid, grid % latitude ) )
      CALL Check( nf90_put_var( ncid, longitude_varid, grid % longitude ) )
      CALL Check( nf90_put_var( ncid, gravity_varid, grid % grx ) )
      CALL Check( nf90_put_var( ncid, fpd_varid, grid % foot_point_distance ) )
      CALL Check( nf90_put_var( ncid, B_varid, grid % magnetic_field_strength ) )
      CALL Check( nf90_put_var( ncid, colat_varid, grid % magnetic_colatitude ) )
      CALL Check( nf90_put_var( ncid, r_varid, grid % r_meter ) )
      CALL Check( nf90_put_var( ncid, q_varid, grid % q_factor ) )

      CALL Check( nf90_put_var( ncid, l11_varid, grid % l_magnitude(1,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, l21_varid, grid % l_magnitude(2,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, l31_varid, grid % l_magnitude(3,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, l12_varid, grid % l_magnitude(1,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, l22_varid, grid % l_magnitude(2,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, l32_varid, grid % l_magnitude(3,2,:,:,:) ) )

      CALL Check( nf90_put_var( ncid, apexe11_varid, grid % apex_e_vectors(1,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexe21_varid, grid % apex_e_vectors(2,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexe31_varid, grid % apex_e_vectors(3,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexe12_varid, grid % apex_e_vectors(1,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexe22_varid, grid % apex_e_vectors(2,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexe32_varid, grid % apex_e_vectors(3,2,:,:,:) ) )

      CALL Check( nf90_put_var( ncid, apexd11_varid, grid % apex_d_vectors(1,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd21_varid, grid % apex_d_vectors(2,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd31_varid, grid % apex_d_vectors(3,1,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd12_varid, grid % apex_d_vectors(1,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd22_varid, grid % apex_d_vectors(2,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd32_varid, grid % apex_d_vectors(3,2,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd13_varid, grid % apex_d_vectors(1,3,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd23_varid, grid % apex_d_vectors(2,3,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, apexd33_varid, grid % apex_d_vectors(3,3,:,:,:) ) )
      CALL Check( nf90_put_var( ncid, d1xd2_varid, grid % d1xd2_magnitude ) )

      CALL Check( nf90_put_var( ncid, be3_varid, grid % apex_be3 ) )
      CALL Check( nf90_put_var( ncid, midpoint_varid, grid % flux_tube_midpoint ) )
      CALL Check( nf90_put_var( ncid, max_varid, grid % flux_tube_max ) ) 
      CALL Check( nf90_put_var( ncid, south_varid, grid % southern_top_index ) ) 
      CALL Check( nf90_put_var( ncid, north_varid, grid % northern_top_index ) ) 

      CALL Check( nf90_close( ncid ) )


  END SUBROUTINE Write_IPE_Grid_NetCDF
!
  SUBROUTINE Read_IPE_Grid_NetCDF( grid, filename )
    IMPLICIT NONE
    CLASS( IPE_Grid ), INTENT(inout) :: grid
    CHARACTER(*), INTENT(in)         :: filename
    ! Local
    INTEGER :: NF90_PREC
    INTEGER :: ncid
    INTEGER :: dimid, varid
    INTEGER :: nFluxtube, NLP, NMP
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

      ! Allocate space for the grid before proceeding
      CALL grid % Build( nFluxTube, NLP, NMP, 1 )

      CALL Check( nf90_inq_varid( ncid, "altitude", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % altitude ) )

      CALL Check( nf90_inq_varid( ncid, "latitude", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % latitude ) )

      CALL Check( nf90_inq_varid( ncid, "longitude", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % longitude ) )

      CALL Check( nf90_inq_varid( ncid, "GRX", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % grx ) )

      CALL Check( nf90_inq_varid( ncid, "foot_point_distance", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % foot_point_distance ) )

      CALL Check( nf90_inq_varid( ncid, "B", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % magnetic_field_strength ) )

      CALL Check( nf90_inq_varid( ncid, "m_colat", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % magnetic_colatitude ) )

      CALL Check( nf90_inq_varid( ncid, "r_meter", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % r_meter ) )

      CALL Check( nf90_inq_varid( ncid, "q_factor", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % q_factor ) )

      CALL Check( nf90_inq_varid( ncid, "l_magnitude_11", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % l_magnitude(1,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "l_magnitude_21", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % l_magnitude(2,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "l_magnitude_31", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % l_magnitude(3,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "l_magnitude_12", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % l_magnitude(1,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "l_magnitude_22", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % l_magnitude(2,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "l_magnitude_32", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % l_magnitude(3,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_e_11", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_e_vectors(1,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_e_21", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_e_vectors(2,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_e_31", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_e_vectors(3,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_e_12", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_e_vectors(1,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_e_22", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_e_vectors(2,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_e_32", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_e_vectors(3,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_11", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(1,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_21", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(2,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_31", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(3,1,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_12", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(1,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_22", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(2,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_32", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(3,2,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_13", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(1,3,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_23", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(2,3,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "apex_d_33", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_d_vectors(3,3,:,:,:) ) )

      CALL Check( nf90_inq_varid( ncid, "d1xd2_mag", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % d1xd2_magnitude ) )

      CALL Check( nf90_inq_varid( ncid, "apex_be3", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % apex_be3 ) )

      CALL Check( nf90_inq_varid( ncid, "tube_midpoint", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % flux_tube_midpoint ) )
      
      CALL Check( nf90_inq_varid( ncid, "tube_max", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % flux_tube_max ) )

      CALL Check( nf90_inq_varid( ncid, "southern_top", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % southern_top_index ) )

      CALL Check( nf90_inq_varid( ncid, "northern_top", varid ) )
      CALL Check( nf90_get_var( ncid, varid, grid % northern_top_index ) )

      CALL Check( nf90_close( ncid ) )


  END SUBROUTINE Read_IPE_Grid_NetCDF

  SUBROUTINE Calculate_Grid_Attributes_From_Legacy_Input( grid )
    IMPLICIT NONE
    CLASS( IPE_Grid ) :: grid
    ! Local
    INTEGER    :: i, lp, mp, ii 
    REAL(prec) :: a(1:3), b(1:3), c(1:3)
    REAL(prec) :: bhat(1:3), ufac

      DO lp = 1, grid % NLP

        ! Calculate the northern top index
        DO i= 1, grid % flux_tube_midpoint(lp)
          
          grid % northern_top_index(lp) = i

          IF ( grid % altitude(i,lp) <= mesh_height_max .AND. grid % altitude(i+1,lp) > mesh_height_max )THEN
            EXIT
          ENDIF

        ENDDO 
         
        ! Calculate the southern top index
        DO i = grid % flux_tube_max(lp), grid % flux_tube_midpoint(lp), -1

          grid % southern_top_index(lp) = i

          IF ( grid % altitude(i,lp) <= mesh_height_max .AND. grid % altitude(i+1,lp) > mesh_height_max )THEN
            EXIT
          ENDIF

        ENDDO

      ENDDO

      DO mp = 1, grid % NMP
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)

            IF( Almost_Equal( grid % apex_d_vectors(1,1,i,lp,mp), 0.0_prec ) .AND. &
                Almost_Equal( grid % apex_d_vectors(1,2,i,lp,mp), 0.0_prec ) )THEN

               ii = i-1

               IF( i /= grid % flux_tube_midpoint(lp) ) then

                 PRINT*, ' '
                 PRINT*, '  Module IPE_Grid_Class : S/R Calculate_Grid_Attributes_From_Legacy_Input'
                 PRINT*, '    WARNING : Invalid apex_d_vectors'
                 PRINT*, ' '

               ENDIF

            ELSE

              ii = i

            ENDIF

            ! calculate D from eq 3.15: | d1 X d2 |
            a(1) = grid % apex_d_vectors(1,1,ii,lp,mp)
            a(2) = grid % apex_d_vectors(2,1,ii,lp,mp)
            a(3) = grid % apex_d_vectors(3,1,ii,lp,mp)

            b(1) = grid % apex_d_vectors(1,2,ii,lp,mp)
            b(2) = grid % apex_d_vectors(2,2,ii,lp,mp)
            b(3) = grid % apex_d_vectors(3,2,ii,lp,mp)

          
            c(1) = a(2)*b(3) - a(3)*b(2) 
            c(2) = a(3)*b(1) - a(1)*b(3) 
            c(3) = a(1)*b(2) - a(2)*b(1) 

            grid % d1xd2_magnitude(i,lp,mp) = SQRT ( c(1)**2 + c(2)**2 + c(3)**2 )

            bhat(1) = grid % apex_d_vectors(1,3,ii,lp,mp)*grid % d1xd2_magnitude(i,lp,mp)
            bhat(2) = grid % apex_d_vectors(2,3,ii,lp,mp)*grid % d1xd2_magnitude(i,lp,mp)
            bhat(3) = grid % apex_d_vectors(3,3,ii,lp,mp)*grid % d1xd2_magnitude(i,lp,mp)
 
            ufac  = SQRT(bhat(1)**2 + bhat(2)**2)

            !(1) magnetic eastward exactly horizontal
            !    l_e = bhat x k /|bhat x k|
            IF( ufac > 0.0_prec ) THEN
              grid % l_magnitude(1,1,i,lp,mp) =  bhat(2)/ufac
              grid % l_magnitude(2,1,i,lp,mp) = -bhat(1)/ufac
              grid % l_magnitude(3,1,i,lp,mp) =  0.0_prec
            ENDIF
      
            ! l_u = l_e x bhat
            grid % l_magnitude(i,lp,mp,1,2) = grid % l_magnitude(2,1,i,lp,mp)*bhat(3) - grid % l_magnitude(3,1,i,lp,mp)*bhat(2)
            grid % l_magnitude(i,lp,mp,2,2) = grid % l_magnitude(3,1,i,lp,mp)*bhat(1) - grid % l_magnitude(1,1,i,lp,mp)*bhat(3)
            grid % l_magnitude(i,lp,mp,3,2) = grid % l_magnitude(1,1,i,lp,mp)*bhat(2) - grid % l_magnitude(2,1,i,lp,mp)*bhat(1)


                                     
            grid % grx(i,lp,mp) = -( G0*earth_radius**2 )/( grid % r_meter(i,lp)**2 )*grid % sinI(i,lp,mp)

          ENDDO
        ENDDO
      ENDDO
        
  END SUBROUTINE Calculate_Grid_Attributes_From_Legacy_Input

  REAL(prec) FUNCTION SinI( grid, i, lp, mp ) 
    IMPLICIT NONE
    CLASS( IPE_Grid ) :: grid
    INTEGER           :: i, lp, mp

      SinI = grid % apex_d_vectors(3,3,i,lp,mp)/SQRT( grid % apex_d_vectors(1,3,i,lp,mp)**2 + &
                                                      grid % apex_d_vectors(2,3,i,lp,mp)**2 + &
                                                      grid % apex_d_vectors(3,3,i,lp,mp)**2  )

  END FUNCTION SinI

END MODULE IPE_Grid_Class

