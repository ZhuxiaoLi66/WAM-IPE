MODULE IPE_Electrodynamics_Class


USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Grid_Class
USE IPE_Forcing_Class
USE IPE_Time_Class

USE efield_ipe
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
    REAL(prec), ALLOCATABLE :: geo_v_ExB_geographic(:,:,:)    ! ExB transport velocity with geographic components
    REAL(prec), ALLOCATABLE :: geo_hall_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: geo_pedersen_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: geo_b_parallel_conductivity(:,:)
    

    CONTAINS 

      PROCEDURE :: Build => Build_IPE_Electrodynamics
      PROCEDURE :: Trash => Trash_IPE_Electrodynamics

      PROCEDURE :: Update => Update_IPE_Electrodynamics
      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Electrodynamics

      PROCEDURE, PRIVATE :: Empirical_E_Field_Wrapper
      PROCEDURE, PRIVATE :: Interpolate_Potential

  END TYPE IPE_Electrodynamics

  LOGICAL, PRIVATE :: setup_efield_empirical

  REAL(prec), PRIVATE :: theta90_rad(0:nmlat)

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
 
  END SUBROUTINE Trash_IPE_Electrodynamics

  SUBROUTINE Update_IPE_Electrodynamics( eldyn, grid, forcing, time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_Forcing ), INTENT(in)             :: forcing
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    ! Local
    INTEGER :: lp, mp

      CALL eldyn % Empirical_E_Field_Wrapper( grid, forcing, time_tracker )

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
    REAL(prec) :: theta130_rad, utime, sunlons, lat
     

        IF( setup_efield_empirical )THEN
          CALL efield_init
          setup_efield_empirical = .FALSE.

          ! Maps latitude from 130km to 90km along flux tube.
          DO j=0,nmlat

            theta130_rad   = ( 90.0_prec - ylatm(j) ) * dtr
            theta90_rad(j) = ASIN(SIN( theta130_rad )*SQRT((earth_radius+90000.0_prec)/(earth_radius+130000.0_prec)))
            IF ( theta130_rad > pi*0.50_prec ) theta90_rad(j) = pi-theta90_rad(j)

          END DO 

          ! Use CommonRoutines HatFunction to build interpolation weights from
          ! theta90_rad to grid % magnetic_latitude
          DO i = 1, grid % NLP
            CALL Hat_Weights( -theta90_rad, -grid % magnetic_colatitude(1,i), &
                              eldyn % lat_interp_weights(1:2,i), &
                              eldyn % lat_interp_index(1:2,i), &
                              nmlat )
          ENDDO

        ENDIF

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


        ! Map magnetic local time to magnetic longitude
        CALL sunloc( year, time_tracker % day_of_year, utime, sunlons ) !iyr,iday,secs)        
        DO i=0,nmlon

          mlon90_rad(i)=(ylonm(i)-180.0_prec)*pi/180.0_prec+sunlons
          IF( mlon90_rad(i) < 0.0_prec   ) mlon90_rad(i)=mlon90_rad(i)+pi*2.0
          IF( mlon90_rad(i) >= pi*2.0_prec ) mlon90_rad(i)=mlon90_rad(i)-pi*2.0

        END DO

        DO i = 1, grid % NMP
          CALL Hat_Weights( mlon90_rad, grid % magnetic_longitude(i), &
                            eldyn % lon_interp_weights(1:2,i), &
                            eldyn % lon_interp_index(1:2,i), &
                            nmlon )
        ENDDO


        ! Interpolate the potential to the IPE grid
        CALL eldyn % Interpolate_Potential( grid )
        
        ! Calculate the potential gradient in IPE coordinates.


  END SUBROUTINE Empirical_E_Field_Wrapper

  SUBROUTINE Interpolate_Potential( eldyn, grid )
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn 
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    ! Local
    INTEGER :: i, j, ii, jj, i_e, j_e
 
    DO j = 1, grid % NMP
      DO i = 1, grid % NLP


        eldyn % electric_potential(i,j) = 0.0_prec
        DO jj = 1, 2
          DO ii = 1, 2

            i_e = eldyn % lat_interp_index(ii,i) 
            j_e = eldyn % lon_interp_index(jj,j) 

            eldyn % electric_potential(i,j) = eldyn % electric_potential(i,j) + &
                                              potent(i_e,j_e)*&
                                              eldyn % lat_interp_weights(ii,i)*&
                                              eldyn % lon_interp_weights(jj,j)

          ENDDO
        ENDDO

      ENDDO
    ENDDO
          
    ! should interpolate from (mlon90_rad, theta90_rad) to ( grid %
    ! magnetic_longitude, grid % magnetic_latitude )
    !
    ! magnetic_latitude = pi/2 - magnetic_colatitude

  END SUBROUTINE Interpolate_Potential

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
