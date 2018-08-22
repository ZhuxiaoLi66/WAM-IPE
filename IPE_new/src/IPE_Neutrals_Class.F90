MODULE IPE_Neutrals_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Grid_Class

! MSIS
USE physics_msis ! gtd7

IMPLICIT NONE

  TYPE IPE_Neutrals

    INTEGER :: nFluxTube, NLP, NMP, nCouplingFields
  
    REAL(prec), ALLOCATABLE :: helium(:,:,:)  ! he
    REAL(prec), ALLOCATABLE :: oxygen(:,:,:)  ! on
    REAL(prec), ALLOCATABLE :: molecular_oxygen(:,:,:)  ! o2n
    REAL(prec), ALLOCATABLE :: molecular_nitrogen(:,:,:) ! n2n
    REAL(prec), ALLOCATABLE :: nitrogen(:,:,:) ! n4s 
    REAL(prec), ALLOCATABLE :: hydrogen(:,:,:) ! hn 
    REAL(prec), ALLOCATABLE :: temperature(:,:,:) 
    REAL(prec), ALLOCATABLE :: temperature_inf(:,:,:) 
    REAL(prec), ALLOCATABLE :: velocity_geographic(:,:,:,:) 
    REAL(prec), ALLOCATABLE :: velocity_apex(:,:,:,:) 
#ifdef COUPLED
    REAL(prec), ALLOCATABLE :: neutral_to_ipe_fields(:,:,:,:)
#endif
    
    ! Interpolated fields
    REAL(prec), ALLOCATABLE :: geo_helium(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_oxygen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_molecular_oxygen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_molecular_nitrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_nitrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_hydrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_velocity(:,:,:,:)


    CONTAINS

      PROCEDURE :: Build => Build_IPE_Neutrals
      PROCEDURE :: Trash => Trash_IPE_Neutrals
    
      PROCEDURE :: Update => Update_IPE_Neutrals
      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Neutrals


      ! PRIVATE Routines
#ifdef COUPLED
      PROCEDURE, PRIVATE :: Copy_WAM_to_IPE_CouplingFields
#endif
      PROCEDURE, PRIVATE :: Geographic_to_Apex_Velocity


  END TYPE IPE_Neutrals


REAL(prec), PARAMETER, PRIVATE :: min_density = 10.0_prec**(-12)
CHARACTER(250), PARAMETER :: hwm_path ='./'

CONTAINS

  SUBROUTINE Build_IPE_Neutrals( neutrals, nFluxTube, NLP, NMP, nCouplingFields )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    INTEGER, INTENT(in)                  :: nFluxTube
    INTEGER, INTENT(in)                  :: NLP
    INTEGER, INTENT(in)                  :: NMP
    INTEGER, OPTIONAL, INTENT(in)        :: nCouplingFields
  

      neutrals % nFluxTube = nFluxTube
      neutrals % NLP       = NLP
      neutrals % NMP       = NMP

      IF( PRESENT( nCouplingFields ) )THEN
        neutrals % nCouplingFields = nCouplingFields
      ENDIF

      ALLOCATE( neutrals % helium(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % oxygen(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % molecular_oxygen(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % molecular_nitrogen(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % nitrogen(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % hydrogen(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % temperature(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % temperature_inf(1:nFluxTube,1:NLP,1:NMP), &
                neutrals % velocity_geographic(1:3,1:nFluxTube,1:NLP,1:NMP), &
                neutrals % velocity_apex(1:3,nFluxTube,1:NLP,1:NMP) )

     neutrals % helium              = 0.0_prec
     neutrals % oxygen              = 0.0_prec
     neutrals % molecular_oxygen    = 0.0_prec
     neutrals % molecular_nitrogen  = 0.0_prec
     neutrals % nitrogen            = 0.0_prec
     neutrals % hydrogen            = 0.0_prec
     neutrals % temperature         = 0.0_prec
     neutrals % temperature_inf     = 0.0_prec
     neutrals % velocity_geographic = 0.0_prec
     neutrals % velocity_apex       = 0.0_prec

      ALLOCATE( neutrals % geo_helium(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_oxygen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_molecular_oxygen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_molecular_nitrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_nitrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_hydrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                neutrals % geo_velocity(1:3,1:nlon_geo,1:nlat_geo,1:nheights_geo) )

     neutrals % geo_helium             = 0.0_prec
     neutrals % geo_oxygen             = 0.0_prec
     neutrals % geo_molecular_oxygen   = 0.0_prec
     neutrals % geo_molecular_nitrogen = 0.0_prec
     neutrals % geo_nitrogen           = 0.0_prec
     neutrals % geo_hydrogen           = 0.0_prec
     neutrals % geo_temperature        = 0.0_prec
     neutrals % geo_velocity           = 0.0_prec

#ifdef COUPLED
     ALLOCATE( neutrals % couplingFields(1:nFluxTube,1:NLP,1:NMP,1:nCouplingFields) )
     neutrals % couplingFields = 0.0_prec
#endif


  END SUBROUTINE Build_IPE_Neutrals
!
  SUBROUTINE Trash_IPE_Neutrals( neutrals )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals


      DEALLOCATE( neutrals % helium, &
                  neutrals % oxygen, &
                  neutrals % molecular_oxygen, &
                  neutrals % molecular_nitrogen, &
                  neutrals % nitrogen, &
                  neutrals % hydrogen, &
                  neutrals % temperature, &
                  neutrals % temperature_inf, &
                  neutrals % velocity_geographic, &
                  neutrals % velocity_apex )

      DEALLOCATE( neutrals % geo_helium, &
                  neutrals % geo_oxygen, &
                  neutrals % geo_molecular_oxygen, &
                  neutrals % geo_molecular_nitrogen, &
                  neutrals % geo_nitrogen, &
                  neutrals % geo_hydrogen, &
                  neutrals % geo_temperature, &
                  neutrals % geo_velocity )

#ifdef COUPLED
     DEALLOCATE( neutrals % couplingFields )
#endif
   
  END SUBROUTINE Trash_IPE_Neutrals
!
  SUBROUTINE Update_IPE_Neutrals( neutrals, grid, utime, year, day, f107d, f107a, ap )
  !
  ! Usage :
  ! 
  !   CALL neutrals % Update( grid, utime, year, day, f107d, f107a, ap )
  ! ================================================================================================== !
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ), INTENT(in)         :: grid
    REAL(prec), INTENT(in)               :: utime
    INTEGER, INTENT(in)                  :: year
    INTEGER, INTENT(in)                  :: day
    REAL(prec), INTENT(in)               :: f107d
    REAL(prec), INTENT(in)               :: f107a
    REAL(prec), INTENT(in)               :: ap(1:7)
    ! Local
    INTEGER        :: i, lp, mp
    INTEGER        :: iyd
    integer :: istop , iht, ilat, ilon
    REAL(prec)     :: slt, lat, lon, alt
    REAL(4)     :: geo_longitude, geo_latitude
    REAL(prec)     :: densities(1:9), temperatures(1:2)
    REAL(4)        :: w(1:2), ap_msis(1:2)

      iyd = 99000 + day
      ap_msis(1:2) = ap(1:2) 

    print *, 'GHGM Starting neutrals'
    do ilon = 1 , 36
      geo_longitude = (ilon - 1) * 10.
      slt = utime/3600.0_prec + geo_longitude/15.0_prec
      do ilat = 1 , 19
        geo_latitude = ((ilat - 1) * 10.) - 90.0
        do iht = 1 , 72
          alt = 90.0 + (iht - 1) * 10.
!         print *, geo_longitude, geo_latitude, alt
          call hwm14( iyd, &           ! Input, year and day as yyddd
                        REAL(utime,4), &         ! Input, universal time ( sec )
                        REAL(alt,4), &           ! Input, altitude ( km )
                        REAL(geo_latitude,4), &           ! Input, geodetic latitude ( degrees )
                        REAL(geo_longitude,4), &           ! Input, geodetic longitude ( degrees )
                        REAL(slt,4), &           ! Input, local apparent solar time ( hrs )[ not used ]
                        REAL(f107a,4), &         ! Input, 3 month average of f10.7 flux [ not used ]
                        REAL(f107d,4), &         ! Input, daily average of f10.7 flux for the previous day [ not used ]
                        ap_msis, &           ! Input, magnetic index ( daily ), current 3hr ap index 
                        hwm_path, &      
                        w(1:2) )         ! Ouput, neutral wind velocity zonal and meridional components
          call gtd7( iyd, &                       ! Input, year and day as yyddd
                       utime, &                     ! Input, universal time ( sec )
                       alt, &                       ! Input, altitude ( km )
                       lat, &                       ! Input, geodetic latitude ( degrees )
                       lon, &                       ! Input, geodetic longitude ( degrees )
                       slt, &                       ! Input, local apparent solar time ( hrs )
                       f107a, &                     ! Input, 3 month average of f10.7 flux
                       f107d, &                     ! Input, daily average of f10.7 flux for the previous day
                       ap(1:7), &                   ! Input, magnetic index ( daily ), current, 3,6,9hrs prior 3hr ap index, 12-33 hr prior ap average, 36-57 hr prior ap average
                       48, &                        ! Mass number ( see src/msis/physics_msis.f90 for more details )
                       densities(1:9), &            ! Ouput, neutral densities in cubic meters
                       temperatures(1:2) )          ! Output, exospheric temperature and temperature at altitude
        enddo
      enddo
    enddo
    print *, 'GHGM Done neutrals'
    istop = 1
    if(istop.eq.1) stop




























      DO mp = 1, neutrals % NMP    
        print *, 'GHGM neutrals mp ', mp
        DO lp = 1, neutrals % NLP    
          print *, 'GHGM neutrals   lp       ', lp
          print *, 'GHGM flux tube size ', grid % flux_tube_max(lp)
          DO i = 1, grid % flux_tube_max(lp)

            lat = 90.0_prec-grid % latitude(i,lp,mp)*180.0_prec/pi
            lon = grid % longitude(i,lp,mp)*180.0_prec/pi
            alt = grid % altitude(i,lp)*0.001_prec
            slt = utime/3600.0_prec + lon/15.0_prec

            ap_msis(1:2) = ap(1:2) 

            ! If the model is not coupled to an external neutral wind model that fills the
            ! neutral wind velocity, then the neutral wind velocity is obtained via the gsw5
            ! routine, defined in src/msis/hwm14.f90
! GHGM only msis up to 800km) .......
            if(alt.le.800.0) then
            print *, 'GHGM neutrals       i           ', i, alt
            call hwm14( iyd, &           ! Input, year and day as yyddd
                        REAL(utime,4), &         ! Input, universal time ( sec )
                        REAL(alt,4), &           ! Input, altitude ( km )
                        REAL(lat,4), &           ! Input, geodetic latitude ( degrees )
                        REAL(lon,4), &           ! Input, geodetic longitude ( degrees )
                        REAL(slt,4), &           ! Input, local apparent solar time ( hrs )[ not used ]
                        REAL(f107a,4), &         ! Input, 3 month average of f10.7 flux [ not used ]
                        REAL(f107d,4), &         ! Input, daily average of f10.7 flux for the previous day [ not used ]
                        ap_msis, &           ! Input, magnetic index ( daily ), current 3hr ap index 
                        hwm_path, &      
                        w(1:2) )         ! Ouput, neutral wind velocity zonal and meridional components
            
            neutrals % velocity_geographic(1,i,lp,mp) = w(2) ! zonal velocity 
            neutrals % velocity_geographic(2,i,lp,mp) = w(1) ! meridional velocity 

            ! Neutral densities are calculated from the gtd7 routine, even if coupling is present
            ! We do this because we don't anticipate models of the thermosphere to output all of the
            ! necessary neutral parameters
            call gtd7( iyd, &                       ! Input, year and day as yyddd
                       utime, &                     ! Input, universal time ( sec )
                       alt, &                       ! Input, altitude ( km )
                       lat, &                       ! Input, geodetic latitude ( degrees )
                       lon, &                       ! Input, geodetic longitude ( degrees )
                       slt, &                       ! Input, local apparent solar time ( hrs )
                       f107a, &                     ! Input, 3 month average of f10.7 flux
                       f107d, &                     ! Input, daily average of f10.7 flux for the previous day
                       ap(1:7), &                   ! Input, magnetic index ( daily ), current, 3,6,9hrs prior 3hr ap index, 12-33 hr prior ap average, 36-57 hr prior ap average
                       48, &                        ! Mass number ( see src/msis/physics_msis.f90 for more details )
                       densities(1:9), &            ! Ouput, neutral densities in cubic meters
                       temperatures(1:2) )          ! Output, exospheric temperature and temperature at altitude

            ! We multiply my 10^6 to convert from cubic meters to cubic centimeters
            neutrals % helium(i,lp,mp)             = densities(1)*10.0_prec**6
            neutrals % oxygen(i,lp,mp)             = densities(2)*10.0_prec**6 
            neutrals % molecular_nitrogen(i,lp,mp) = densities(3)*10.0_prec**6
            neutrals % molecular_oxygen(i,lp,mp)   = densities(4)*10.0_prec**6
            neutrals % hydrogen(i,lp,mp)           = densities(7)*10.0_prec**6
            neutrals % nitrogen(i,lp,mp)           = densities(8)*10.0_prec**6

            neutrals % temperature_inf(i,lp,mp)    = temperatures(1)
            neutrals % temperature(i,lp,mp)        = temperatures(2)
            endif  ! GHGM if(alt.le.800.0) then
     

          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

#ifdef COUPLED

#ifdef COUPLED_TO_WAM

      CALL neutrals % Copy_WAM_to_IPE_CouplingFields( grid )

#endif

#endif

    CALL neutrals % Geographic_to_Apex_Velocity( grid )

  END SUBROUTINE Update_IPE_Neutrals
!
#ifdef COUPLED
  SUBROUTINE Copy_WAM_to_IPE_CouplingFields( neutrals, grid )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ), INTENT(in)         :: grid
    ! Local 
    INTEGER    :: i, lp, mp
    INTEGER    :: istep, ihtop, midpoints, ipts, ipts1, ihem
    REAL(prec) :: scale_height_O, scale_height_O2, scale_height_N2, dht

      !---------------------------------------------------------------------------
      ! The WamField array has indices of:
      !           1:Tn; 2:Vn_geographic_east; 3:Vn_geographic_north; 4:Vn_geographic_up; 
      !           5:O_density; 6:O2_density; 7:N2_density
      !
      ! NH and SH denote Northern and Southern Hemispheres respectively
      !---------------------------------------------------------------------------
      DO mp = 1, neutrals % NMP    
        DO lp = 1, neutrals % NLP    

          ! Northern hemisphere, altitude < 800 km
          DO i = 1, grid % northern_top_index(lp,mp)

            neutrals % temperature(i,lp,mp)           = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,1)
            neutrals % velocity_geographic(1,i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,2)
            neutrals % velocity_geographic(2,i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,3)
            neutrals % velocity_geographic(3,i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,4)
            neutrals % oxygen(i,lp,mp)                = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,5) 
            neutrals % molecular_oxygen(i,lp,mp)      = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,6)
            neutrals % molecular_nitrogen(i,lp,mp)    = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,7)

          ENDDO

          ! Southern hemisphere, altitude < 800 km
          DO i = 1, grid % northern_top_index(lp,mp)

            neutrals % temperature(i,lp,mp)           = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,1)
            neutrals % velocity_geographic(1,i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,2)
            neutrals % velocity_geographic(2,i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,3)
            neutrals % velocity_geographic(3,i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,4)
            neutrals % oxygen(i,lp,mp)                = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,5) 
            neutrals % molecular_oxygen(i,lp,mp)      = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,6)
            neutrals % molecular_nitrogen(i,lp,mp)    = neutrals % neutral_to_ipe_coupling_field(i,lp,mp,7)

          ENDDO

          ! At altitudes higher than 800 km, the neutral temperature is extrapolated through prolongation

          ! Northern hemisphere, altitude > 800 km
          DO i = grid % northern_top_index(lp,mp)+1, grid % flux_tube_midpoint(lp,mp) 

            neutrals % temperature(i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(grid % northern_top_index(lp,mp),lp,mp,1)

          ENDDO

          ! Southern hemisphere, altitude > 800 km
          DO i = grid % flux_tube_midpoint(lp,mp) + 1, grid % southern_top_index(lp,mp)-1

            neutrals % temperature(i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(grid % southern_top_index(lp,mp),lp,mp,1)

          ENDDO

          DO i = 1, grid % flux_tube_midpoint(lp,mp)
         
            neutrals % temperature_inf(i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(grid % northern_top_index(lp,mp),lp,mp,1)
         
          ENDDO


          DO i = grid % flux_tube_midpoint(lp,mp)+1, grid % flux_tube_max(lp,mp)
         
            neutrals % temperature_inf(i,lp,mp) = neutrals % neutral_to_ipe_coupling_field(grid % southern_top_index(lp,mp),lp,mp,1)
         
          ENDDO

          DO ihem=1, 2

            IF( ihem==1 )THEN 

              istep = 1
              ihTop = grid % northern_top_index(lp,mp)
              midPoints = grid % flux_tube_midpoint(lp,mp)

            ELSEIF( ihem==2 )THEN

              istep = -1
              ihTop = grid % southern_top_index(lp,mp)
              midPoints = grid % flux_tube_midpoint(lp,mp)+1            

            ENDIF
 
            DO ipts=ihTop+iStep, midPoints, iStep 
                     
              ! At altitudes higher than 800 km, neutral velocity components are extrapolated through prolongation 
              neutrals % velocity_geographic(1,ipts,lp,mp) =  neutrals % neutral_to_ipe_coupling_field(ihTop,lp,mp,2) 
              neutrals % velocity_geographic(2,ipts,lp,mp) =  neutrals % neutral_to_ipe_coupling_field(ihTop,lp,mp,3) 
              neutrals % velocity_geographic(3,ipts,lp,mp) =  neutrals % neutral_to_ipe_coupling_field(ihTop,lp,mp,4) 

              !---------------------------------------------------------------------------
              ! O, O2, and N2 densities above 800km decrease exponentially with height
              !---------------------------------------------------------------------------

              ipts1=ipts-iStep 

              dist = earth_radius/(grid % earth_radius + grid % altitude(ipts,lp))
              dist1 = earth_radius/(grid % earth_radius + grid % altitude(ipts1,lp))

              Hk_O  = GSCON * neutrals % temperature(ipts ,lp,mp) / (massn_kg(1)*G0*dist*dist)
              Hk1_O = GSCON * neutrals % temperature(ipts1,lp,mp) / (massn_kg(1)*G0*dist1*dist1)
              Hk_O2  = GSCON * neutrals % temperature(ipts ,lp,mp) / (massn_kg(2)*G0*dist*dist)
              Hk1_O2 = GSCON * neutrals % temperature(ipts1,lp,mp) / (massn_kg(2)*G0*dist1*dist1)
              Hk_N2  = GSCON * neutrals % temperature(ipts ,lp,mp) / (massn_kg(3)*G0*dist*dist)
              Hk1_N2 = GSCON * neutrals % temperature(ipts1,lp,mp) / (massn_kg(3)*G0*dist1*dist1)

              scale_height_O  = (Hk_O+Hk1_O)*0.5_prec
              scale_height_O2 = (Hk_O2+Hk1_O2)*0.5_prec
              scale_height_N2 = (Hk_N2+Hk1_N2)*0.5_prec

              dht = grid % altitude(ipts1,lp) - grid % altitude(ipts,lp)

              neutrals % oxygen(ipts,lp,mp)             =  neutrals % oxygen(ipts1,lp,mp) * exp(dht/scale_height_O) 
              neutrals % molecular_oxygen(ipts,lp,mp)   =  neutrals % molecular_oxygen(ipts1,lp,mp) * exp(dht/scale_height_O2)
              neutrals % molecular_nitrogen(ipts,lp,mp) =  neutrals % molecular_nitrogen(ipts1,lp,mp) * exp(dht/scale_height_N2) 
                        
              !---------------------------------------------------------------------------
              ! Finally, Make sure these densities do not go to zero at high altitudes....
              !---------------------------------------------------------------------------
              
              IF( neutrals % oxygen(ipts,lp,mp) <= min_density) neutrals % oxygen(ipts,lp,mp) = min_density
              IF( neutrals % molecular_oxygen(ipts,lp,mp) <= min_density) neutrals % molecular_oxygen(ipts,lp,mp) = min_density
              IF( neutrals % molecular_nitrogen(ipts,lp,mp) <=min_density) neutrals % molecular_nitrogen(ipts,lp,mp) = min_density

            ENDDO

          ENDDO
                        

        ENDDO
      ENDDO

  END SUBROUTINE Copy_WAM_to_IPE_CouplingFields
#endif
!
  SUBROUTINE Geographic_to_Apex_Velocity( neutrals, grid )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ), INTENT(in)         :: grid
    ! Local 
    INTEGER    :: i, lp, mp, idir
    REAL(prec) :: dotprod


      DO mp = 1, neutrals % NMP
        DO lp = 1, neutrals % NLP
          DO i= 1, neutrals % nFluxTube

            DO idir = 1, 3

               dotprod = grid % apex_d_vectors(1,idir,i,lp,mp)**2 + &
                         grid % apex_d_vectors(2,idir,i,lp,mp)**2 + &
                         grid % apex_d_vectors(3,idir,i,lp,mp)**2

               IF ( dotprod > 0.0_prec ) THEN
               
                 neutrals % velocity_apex(idir,i,lp,mp) = ( grid % apex_d_vectors(1,idir,i,lp,mp)*neutrals % velocity_geographic(1,i,lp,mp) + &
                                                            grid % apex_d_vectors(2,idir,i,lp,mp)*neutrals % velocity_geographic(2,i,lp,mp) + &
                                                            grid % apex_d_vectors(3,idir,i,lp,mp)*neutrals % velocity_geographic(3,i,lp,mp) )/&
                                                            SQRT( dotprod )

               ELSE

                  neutrals % velocity_apex(idir,i,lp,mp) = 0.0_prec

               END IF

             ENDDO

             ! dbg20110131 : the midpoint values become NaN otherwise because
             ! of inappropriate D1/3 values...
             IF ( lp >= 1 .AND. lp <= 6 )THEN
               IF( i == grid % flux_tube_midpoint(lp) )THEN
                 neutrals % velocity_apex(1:3,i,lp,mp) = neutrals % velocity_apex(1:3,i-1,lp,mp) 
               ENDIF
             ENDIF

          END DO
        END DO 
      END DO   

  END SUBROUTINE Geographic_to_Apex_Velocity
!
  SUBROUTINE Interpolate_to_GeographicGrid_IPE_Neutrals( neutrals, grid )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ), INTENT(in)         :: grid

      CALL grid % Interpolate_to_Geographic_Grid( neutrals % helium, neutrals % geo_helium )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % oxygen, neutrals % geo_oxygen )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % molecular_oxygen, neutrals % geo_molecular_oxygen )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % molecular_nitrogen, neutrals % geo_molecular_nitrogen )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % nitrogen, neutrals % geo_nitrogen )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % hydrogen, neutrals % geo_hydrogen )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % temperature, neutrals % geo_temperature )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % velocity_geographic(1,:,:,:), neutrals % geo_velocity(1,:,:,:) )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % velocity_geographic(2,:,:,:), neutrals % geo_velocity(2,:,:,:) )
      CALL grid % Interpolate_to_Geographic_Grid( neutrals % velocity_geographic(3,:,:,:), neutrals % geo_velocity(3,:,:,:) )

  END SUBROUTINE Interpolate_to_GeographicGrid_IPE_Neutrals



  SUBROUTINE INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE ( &
        geo_grid_longitudes_degrees,geo_grid_latitudes_degrees, &
        O_density,O2_density,N2_density, &
        tiegcm_op_density_fixed_ht , tiegcm_nop_density_fixed_ht , tiegcm_o2p_density_fixed_ht , &   ! test with tiegcm ion inputs.....
        NO_density, &
        N4S_density,N2D_density, &
        VX,VY,WVZ,TTS,telec, &
        IN,IS,IWRite2,TN_plasma_input_3d,O_plasma_input_3d, &
        O2_plasma_input_3d,N2_plasma_input_3d, &
        tiegcm_OP_input_3d,tiegcm_NOP_input_3d,tiegcm_O2P_input_3d, &
        NO_plasma_input_3d,N4S_plasma_input_3d,N2D_plasma_input_3d, &
        GLAt_3d, &
        GLOnd_3d, &
        PZ_3d, &
        um_plasma_input_3d,uz_plasma_input_3d,uv_plasma_input_3d, &
        te_plasma_input_3d, &
        ilon1_3d_fixed_ht,ilon2_3d_fixed_ht, &
        ilat1_3d_fixed_ht,ilat2_3d_fixed_ht,ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht,ihu_3d_fixed_ht, &
        sw_1st_call_int_fixed_ht, &
        GIP_switches)





! this calculates a neutral background for the
! plasmasphere code by interpolating values from the
! fixed height interface.....

  IMPLICIT NONE
    INTEGER :: N_heights
!nm032007:,N_Latitudes,N_longitudes
    PARAMETER (N_heights=interface_hts &
!nm032007: ,N_Latitudes=91,N_longitudes=20
     &)

  LOGICAL :: sw_1st_call_int_fixed_ht
  LOGICAL :: GIP_switches(20)
  LOGICAL :: sw_External_model_provides_NO_N4S_densities

  REAL(kind=8) :: dNO1 , dNO11 , dNO12 , dNO2 , dNO21 , dNO22 , &
  dNOl11 , dNOl12 , dNOl21 , dNOl22 , dNOu11 , dNOu12, &
  dNOu21 , dNOu22
  REAL(kind=8) :: dN4S1 , dN4S11 , dN4S12 , dN4S2 , dN4S21 , dN4S22 , &
  dN4Sl11 , dN4Sl12 , dN4Sl21 , dN4Sl22 , dN4Su11 , dN4Su12, &
  dN4Su21 , dN4Su22
  REAL(kind=8) :: dN2D1 , dN2D11 , dN2D12 , dN2D2 , dN2D21 , dN2D22 , &
  dN2Dl11 , dN2Dl12 , dN2Dl21 , dN2Dl22 , dN2Du11 , dN2Du12, &
  dN2Du21 , dN2Du22

  REAL(kind=8) :: dnn1 , dnn11 , dnn12 , dnn2 , dnn21 , dnn22 , &
  dnnl11 , dnnl12 , dnnl21 , dnnl22 , dnnu11 , dnnu12
  REAL(kind=8) :: &
  dnnu21 , dnnu22 , fach,  &
  faclat , faclon , &
  glond2

  REAL(kind=8) :: ol11 , ol12 , ol21 , ol22 ,  ou11 , ou12 , ou21 , ou22 , o11 , o12 , o21 , o22 , do1 , do2
  REAL(kind=8) :: ool11 , ool12 , ool21 , ool22 ,  oou11 , oou12 , oou21 , oou22 , oo11 , oo12 , oo21 , oo22  , doo1 , doo2
  REAL(kind=8) :: tnl11 , tnl12 , tnl21 , tnl22 ,  tnu11 , tnu12 , tnu21 , tnu22 , tn11 , tn12 , tn21 , tn22  , tn1 , tn2
  REAL(kind=8) :: tel11 , tel12 , tel21 , tel22 ,  teu11 , teu12 , teu21 , teu22 , te11 , te12 , te21 , te22 , te1 , te2
  REAL(kind=8) :: uml11 , uml12 , uml21 , uml22 ,  umu11 , umu12 , umu21 , umu22 , um11 , um12 , um21 , um22 , um1 , um2
  REAL(kind=8) :: uzl11 , uzl12 , uzl21 , uzl22 ,  uzu11 , uzu12 , uzu21 , uzu22 , uz11 , uz12 , uz21 , uz22 , uz1 , uz2
  REAL(kind=8) :: uvl11 , uvl12 , uvl21 , uvl22 ,  uvu11 , uvu12 , uvu21 , uvu22 , uv11 , uv12 , uv21 , uv22 , uv1 , uv2
  REAL(kind=8) :: topl11 , topl12 , topl21 , topl22 ,  topu11 , topu12 , topu21 , topu22 , top11 , top12 , top21 , top22  , top1 , top2
  REAL(kind=8) :: tnopl11 , tnopl12 , tnopl21 , tnopl22 ,  tnopu11 , tnopu12 , tnopu21 , tnopu22 , tnop11 , tnop12 , tnop21 , tnop22 , tnop1 , tnop2
  REAL(kind=8) :: to2pl11 , to2pl12 , to2pl21 , to2pl22 ,  to2pu11 , to2pu12 , to2pu21 , to2pu22 , to2p11 , to2p12 , to2p21 , to2p22 , to2p1 , to2p2

  INTEGER :: i , ih , ihl , ihu, ii , ilat , ilat1 , ilat2 , ilon , ilon1 , ilon2 , iprob
  INTEGER :: ispecial , IWRite , l , m , n ,  ifault , itube , iwrite1 , iwrite2 , istop

  REAL(kind=8) :: fixed_hts_in_km(N_heights)
  REAL(kind=8) :: O_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: O2_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: N2_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: tiegcm_op_density_fixed_ht(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: tiegcm_nop_density_fixed_ht(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: tiegcm_o2p_density_fixed_ht(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: NO_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: N4S_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: N2D_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: VX(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: VY(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: WVZ(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: TTS(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: telec(N_heights,N_Latitudes,N_longitudes)

  REAL(kind=8) :: pzh(N_heights)

  integer :: in(nmp,nlp),is(nmp,nlp),mp,lp
  REAL(kind=8) tn_plasma_input_3d(npts,nmp), &
  um_plasma_input_3d(npts,nmp),uz_plasma_input_3d(npts,nmp),uv_plasma_input_3d(npts,nmp), &
  o_plasma_input_3d(npts,nmp),o2_plasma_input_3d(npts,nmp),n2_plasma_input_3d(npts,nmp), &
  te_plasma_input_3d(npts,nmp),te_dum(npts), &
  NO_plasma_input_3d(npts,nmp), &
  N4S_plasma_input_3d(npts,nmp), N2D_plasma_input_3d(npts,nmp), &
  tiegcm_OP_input_3d(npts,nmp),tiegcm_NOP_input_3d(npts,nmp), tiegcm_O2P_input_3d(npts,nmp)

  REAL(kind=8) ::  glat_3d(npts,nmp),glond_3d(npts,nmp),pz_3d(npts,nmp)
  REAL(kind=8) ::  NO(npts)
  REAL(kind=8) ::  N4S(npts)
  REAL(kind=8) ::  N2D(npts)

  REAL(kind=8) ::  uz(NPTS)
  REAL(kind=8) ::  um(NPTS)
  REAL(kind=8) ::  uv(NPTS)

  REAL(kind=8) ::  TN(NPTS) , O(NPTS) , O2(NPTS) , N2(NPTS) , GLAt(NPTS) , &
  PZ(NPTS) , GLOnd(NPTS)
  REAL(kind=8) ::  tiegcm_OP_on_tube(npts)
  REAL(kind=8) ::  tiegcm_NOP_on_tube(npts)
  REAL(kind=8) ::  tiegcm_O2P_on_tube(npts)

  INTEGER ::   ilon1_3d_fixed_ht(npts,nmp),ilon2_3d_fixed_ht(npts,nmp)
  INTEGER ::   ilat1_3d_fixed_ht(npts,nmp),ilat2_3d_fixed_ht(npts,nmp)
  INTEGER ::   ispecial_3d_fixed_ht(npts,nmp)
  INTEGER ::   ihl_3d_fixed_ht(npts,nmp),ihu_3d_fixed_ht(npts,nmp)
  REAL(kind=8) pz_1000(npts)
  REAL(kind=8) geo_grid_longitudes_degrees(N_longitudes)
  REAL(kind=8) geo_grid_latitudes_degrees(N_Latitudes)
  REAL(kind=8) small_power,small_number
  small_power = -20.
  small_number = 1.d-20

  sw_External_model_provides_NO_N4S_densities = GIP_switches(5)

!   do l = 1 , N_longitudes
!     geo_grid_longitudes_degrees(l) = (float(l-1))*18.
!   enddo

!   do m = 1 , N_latitudes
!     geo_grid_latitudes_degrees(m) = (float(m - 46))*2.
!   enddo

  !write(6,*) '************************************'
  !write(6,*) 'longs ',geo_grid_longitudes_degrees
  !write(6,*) '************************************'
  !write(6,*) 'lats ',geo_grid_latitudes_degrees
  !write(6,*) '************************************'

!g
  iwrite1 = 0
  iwrite=0
  istop = 0
  if (istop == 1) stop
  do i = 1, N_heights
!   fixed_hts_in_km(i) = (float(i-1)*5.) + 90.
    fixed_hts_in_km(i) = fixed_heights_km(i)
    pzh(i) = fixed_hts_in_km(i)
  enddo
!g
!g  Big loop over all flux tubes (nmp and nlp).....
!g
!       write(6,*) '***** calling interface to plasma ***** '
  do mp = 1 , nmp
      !write(6,*) 'interface to plasma fixed_ht',mp,sw_1st_call_int_fixed_ht
      do lp = 1 , nlp
      !write(6,*) '     lp ',lp
      !g
      !g  calculate the 1D geographic inputs......
      !g
          do i=in(mp,lp),is(mp,lp)
              glat(i) = glat_3d(i,mp)
              glond(i) = glond_3d(i,mp)
              pz(i) = pz_3d(i,mp)
              pz_1000(i) = pz(i)*1000.
          enddo
      !g

      ! loop over all points on the tube...

          DO 900 i = IN(mp,lp) , IS(mp,lp)
              !write(6,*) 'points ',i
              glond2 = GLOnd(i)
              IF ( glond2 > 360. ) glond2 = glond2 - 360.
              IF ( glond2 < 0. ) glond2 = glond2 + 360.

              if(sw_1st_call_int_fixed_ht) then

                  ispecial = 0
                  DO 200 ilon = 1 , N_longitudes
                      IF ( geo_grid_longitudes_degrees(ilon) > glond2 ) THEN
                          ilon2 = ilon
                          ilon1 = ilon - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , glond2 , &
                          geo_grid_longitudes_degrees(ilon2) , geo_grid_longitudes_degrees(ilon1)
                          99001 FORMAT ('this 1 ',i3,3(2x,f10.1))
                          GOTO 250
                      ENDIF
                  200 ENDDO
                  ilon2 = 1
                  ilon1 = N_longitudes
                  ispecial = 1

                  250 DO 300 ilat = 1 , N_Latitudes
                      IF ( geo_grid_latitudes_degrees(ilat) > GLAt(i) ) THEN
                          ilat2 = ilat
                          ilat1 = ilat - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , GLAt(i) , &
                          geo_grid_latitudes_degrees(ilat2) , geo_grid_latitudes_degrees(ilat1)
                          GOTO 350
                      ENDIF
                  300 ENDDO
                  350 CONTINUE
              !g
              !g thus the required field line point lies
              !g within the square with latitudes ilat2,ilat1
              !g and longitudes ilon2,ilon1....at these four
              !g points there are two heights which are above
              !g and below the point.
              !g

                  DO 400 ih = 1 , N_heights
                      pzh(ih) = fixed_hts_in_km(ih)
                      IF ( pzh(ih) > PZ(i) ) THEN
                          ihu = ih
                          ihl = ih - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , PZ(i) , &
                          pzh(ihu) , pzh(ihl)
                          GOTO 450
                      ENDIF
                  400 ENDDO
                  ihu = N_heights
                  ihl = N_heights - 1
                  IF ( IWRite == 1 ) WRITE(6,99001) i , PZ(i) , pzh(ihu) , &
                  pzh(ihl)
                  450 IF ( ihu == 1 ) THEN
                      WRITE(6,*) 'ihu ' , ihu , pzh(ihu) , pzh(ihl)
                      DO 460 ii = IN(mp,lp) , IS(mp,lp)
                          WRITE(6,*) ii , PZ(ii)
                      460 ENDDO
                  ENDIF

              ! thus point lies in a box surrounded by the 8 points...

              ! ilon1 ilat1 ihl
              ! ilon1 ilat1 ihu
              ! ilon1 ilat2 ihl
              ! ilon1 ilat2 ihu
              ! ilon2 ilat1 ihl
              ! ilon2 ilat1 ihu
              ! ilon2 ilat2 ihl
              ! ilon2 ilat2 ihu

                  ilon1_3d_fixed_ht(i,mp) = ilon1
                  ilon2_3d_fixed_ht(i,mp) = ilon2
                  ilat1_3d_fixed_ht(i,mp) = ilat1
                  ilat2_3d_fixed_ht(i,mp) = ilat2
                  ispecial_3d_fixed_ht(i,mp) = ispecial
                  ihl_3d_fixed_ht(i,mp) = ihl
                  ihu_3d_fixed_ht(i,mp) = ihu

              else

                  ilon1 = ilon1_3d_fixed_ht(i,mp)
                  ilon2 = ilon2_3d_fixed_ht(i,mp)
                  ilat1 = ilat1_3d_fixed_ht(i,mp)
                  ilat2 = ilat2_3d_fixed_ht(i,mp)
                  ispecial = ispecial_3d_fixed_ht(i,mp)
                  ihl = ihl_3d_fixed_ht(i,mp)
                  ihu = ihu_3d_fixed_ht(i,mp)
!                   write(6,*) 'arggh ',ilon1,ilon2,ilat1,ilat2,ispecial,ihl,ihu
              endif

          ! neutral temperature on the eight surrounding points......

              tnu11 = TTS(ihu,ilat1,ilon1)
              tnl11 = TTS(ihl,ilat1,ilon1)
              tnu12 = TTS(ihu,ilat1,ilon2)
              tnl12 = TTS(ihl,ilat1,ilon2)
              tnu21 = TTS(ihu,ilat2,ilon1)
              tnl21 = TTS(ihl,ilat2,ilon1)
              tnu22 = TTS(ihu,ilat2,ilon2)
              tnl22 = TTS(ihl,ilat2,ilon2)

          ! electron temperature on the eight surrounding points......

!             teu11 = Telec(ihu,ilat1,ilon1)
!             tel11 = Telec(ihl,ilat1,ilon1)
!             teu12 = Telec(ihu,ilat1,ilon2)
!             tel12 = Telec(ihl,ilat1,ilon2)
!             teu21 = Telec(ihu,ilat2,ilon1)
!             tel21 = Telec(ihl,ilat2,ilon1)
!             teu22 = Telec(ihu,ilat2,ilon2)
!             tel22 = Telec(ihl,ilat2,ilon2)

          ! zonal wind on the eight surrounding points......

              uzu11 = VY(ihu,ilat1,ilon1)
              uzl11 = VY(ihl,ilat1,ilon1)
              uzu12 = VY(ihu,ilat1,ilon2)
              uzl12 = VY(ihl,ilat1,ilon2)
              uzu21 = VY(ihu,ilat2,ilon1)
              uzl21 = VY(ihl,ilat2,ilon1)
              uzu22 = VY(ihu,ilat2,ilon2)
              uzl22 = VY(ihl,ilat2,ilon2)

          ! meridional wind on the eight surrounding points......

              umu11 = VX(ihu,ilat1,ilon1)
              uml11 = VX(ihl,ilat1,ilon1)
              umu12 = VX(ihu,ilat1,ilon2)
              uml12 = VX(ihl,ilat1,ilon2)
              umu21 = VX(ihu,ilat2,ilon1)
              uml21 = VX(ihl,ilat2,ilon1)
              umu22 = VX(ihu,ilat2,ilon2)
              uml22 = VX(ihl,ilat2,ilon2)

          ! vertical wind on the eight surrounding points......

              uvu11 = WVZ(ihu,ilat1,ilon1)
              uvl11 = WVZ(ihl,ilat1,ilon1)
              uvu12 = WVZ(ihu,ilat1,ilon2)
              uvl12 = WVZ(ihl,ilat1,ilon2)
              uvu21 = WVZ(ihu,ilat2,ilon1)
              uvl21 = WVZ(ihl,ilat2,ilon1)
              uvu22 = WVZ(ihu,ilat2,ilon2)
              uvl22 = WVZ(ihl,ilat2,ilon2)

          ! atomic oxygen density on the eight surrounding points......

              ou11 = log10(O_density(ihu,ilat1,ilon1))
              ol11 = log10(O_density(ihl,ilat1,ilon1))
              ou12 = log10(O_density(ihu,ilat1,ilon2))
              ol12 = log10(O_density(ihl,ilat1,ilon2))
              ou21 = log10(O_density(ihu,ilat2,ilon1))
              ol21 = log10(O_density(ihl,ilat2,ilon1))
              ou22 = log10(O_density(ihu,ilat2,ilon2))
              ol22 = log10(O_density(ihl,ilat2,ilon2))

          ! molecular oxygen density on the eight surrounding points......

              oou11 = log10(O2_density(ihu,ilat1,ilon1))
              ool11 = log10(O2_density(ihl,ilat1,ilon1))
              oou12 = log10(O2_density(ihu,ilat1,ilon2))
              ool12 = log10(O2_density(ihl,ilat1,ilon2))
              oou21 = log10(O2_density(ihu,ilat2,ilon1))
              ool21 = log10(O2_density(ihl,ilat2,ilon1))
              oou22 = log10(O2_density(ihu,ilat2,ilon2))
              ool22 = log10(O2_density(ihl,ilat2,ilon2))

          ! molecular nitrogen density on the eight surrounding points......

              dnnu11 = log10(N2_density(ihu,ilat1,ilon1))
              dnnl11 = log10(N2_density(ihl,ilat1,ilon1))
              dnnu12 = log10(N2_density(ihu,ilat1,ilon2))
              dnnl12 = log10(N2_density(ihl,ilat1,ilon2))
              dnnu21 = log10(N2_density(ihu,ilat2,ilon1))
              dnnl21 = log10(N2_density(ihl,ilat2,ilon1))
              dnnu22 = log10(N2_density(ihu,ilat2,ilon2))
              dnnl22 = log10(N2_density(ihl,ilat2,ilon2))

          ! tiegcm O+ density on the eight surrounding points......

              topu11 = tiegcm_op_density_fixed_ht(ihu,ilat1,ilon1)
              topl11 = tiegcm_op_density_fixed_ht(ihl,ilat1,ilon1)
              topu12 = tiegcm_op_density_fixed_ht(ihu,ilat1,ilon2)
              topl12 = tiegcm_op_density_fixed_ht(ihl,ilat1,ilon2)
              topu21 = tiegcm_op_density_fixed_ht(ihu,ilat2,ilon1)
              topl21 = tiegcm_op_density_fixed_ht(ihl,ilat2,ilon1)
              topu22 = tiegcm_op_density_fixed_ht(ihu,ilat2,ilon2)
              topl22 = tiegcm_op_density_fixed_ht(ihl,ilat2,ilon2)

          ! tiegcm NO+ density on the eight surrounding points......

              tnopu11 = tiegcm_nop_density_fixed_ht(ihu,ilat1,ilon1)
              tnopl11 = tiegcm_nop_density_fixed_ht(ihl,ilat1,ilon1)
              tnopu12 = tiegcm_nop_density_fixed_ht(ihu,ilat1,ilon2)
              tnopl12 = tiegcm_nop_density_fixed_ht(ihl,ilat1,ilon2)
              tnopu21 = tiegcm_nop_density_fixed_ht(ihu,ilat2,ilon1)
              tnopl21 = tiegcm_nop_density_fixed_ht(ihl,ilat2,ilon1)
              tnopu22 = tiegcm_nop_density_fixed_ht(ihu,ilat2,ilon2)
              tnopl22 = tiegcm_nop_density_fixed_ht(ihl,ilat2,ilon2)

          ! tiegcm O2+ density on the eight surrounding points......

              to2pu11 = tiegcm_o2p_density_fixed_ht(ihu,ilat1,ilon1)
              to2pl11 = tiegcm_o2p_density_fixed_ht(ihl,ilat1,ilon1)
              to2pu12 = tiegcm_o2p_density_fixed_ht(ihu,ilat1,ilon2)
              to2pl12 = tiegcm_o2p_density_fixed_ht(ihl,ilat1,ilon2)
              to2pu21 = tiegcm_o2p_density_fixed_ht(ihu,ilat2,ilon1)
              to2pl21 = tiegcm_o2p_density_fixed_ht(ihl,ilat2,ilon1)
              to2pu22 = tiegcm_o2p_density_fixed_ht(ihu,ilat2,ilon2)
              to2pl22 = tiegcm_o2p_density_fixed_ht(ihl,ilat2,ilon2)

          if(sw_External_model_provides_NO_N4S_densities) then
          ! NO density on the eight surrounding points......

              dNOu11 = log10(NO_density(ihu,ilat1,ilon1))
              dNOl11 = log10(NO_density(ihl,ilat1,ilon1))
              dNOu12 = log10(NO_density(ihu,ilat1,ilon2))
              dNOl12 = log10(NO_density(ihl,ilat1,ilon2))
              dNOu21 = log10(NO_density(ihu,ilat2,ilon1))
              dNOl21 = log10(NO_density(ihl,ilat2,ilon1))
              dNOu22 = log10(NO_density(ihu,ilat2,ilon2))
              dNOl22 = log10(NO_density(ihl,ilat2,ilon2))

          ! N4S density on the eight surrounding points......

              dN4Su11 = log10(N4S_density(ihu,ilat1,ilon1))
              dN4Sl11 = log10(N4S_density(ihl,ilat1,ilon1))
              dN4Su12 = log10(N4S_density(ihu,ilat1,ilon2))
              dN4Sl12 = log10(N4S_density(ihl,ilat1,ilon2))
              dN4Su21 = log10(N4S_density(ihu,ilat2,ilon1))
              dN4Sl21 = log10(N4S_density(ihl,ilat2,ilon1))
              dN4Su22 = log10(N4S_density(ihu,ilat2,ilon2))
              dN4Sl22 = log10(N4S_density(ihl,ilat2,ilon2))

          ! N2D density on the eight surrounding points......

              dN2Du11 = log10(N2D_density(ihu,ilat1,ilon1))
              dN2Dl11 = log10(N2D_density(ihl,ilat1,ilon1))
              dN2Du12 = log10(N2D_density(ihu,ilat1,ilon2))
              dN2Dl12 = log10(N2D_density(ihl,ilat1,ilon2))
              dN2Du21 = log10(N2D_density(ihu,ilat2,ilon1))
              dN2Dl21 = log10(N2D_density(ihl,ilat2,ilon1))
              dN2Du22 = log10(N2D_density(ihu,ilat2,ilon2))
              dN2Dl22 = log10(N2D_density(ihl,ilat2,ilon2))

          endif

          ! now the 8 point interpolation.........

              fach = (PZ(i)-pzh(ihl))/(pzh(ihu)-pzh(ihl))
             if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'ZZZZZZZ ',ihl,ihu
!               write(6,*) 'XXXXXXX ',pz(i),pzh(ihl),pzh(ihu)
             endif

              IF ( fach > 1. ) THEN
                  tn11 = tnu11
                  uz11 = uzu11
                  um11 = umu11
                  uv11 = 0.
                  tn12 = tnu12
                  uz12 = uzu12
                  um12 = umu12
                  uv12 = 0.
                  tn21 = tnu21
                  uz21 = uzu21
                  um21 = umu21
                  uv21 = 0.
                  tn22 = tnu22
                  uz22 = uzu22
                  um22 = umu22
                  uv22 = 0.
!                 te11 = teu11+pz(i)-pzh(ihu)
!                 te12 = teu12+pz(i)-pzh(ihu)
!                 te21 = teu21+pz(i)-pzh(ihu)
!                 te22 = teu22+pz(i)-pzh(ihu)
              ELSE
                  tn11 = ((tnu11-tnl11)*fach) + tnl11
                  uz11 = ((uzu11-uzl11)*fach) + uzl11
                  um11 = ((umu11-uml11)*fach) + uml11
                  uv11 = ((uvu11-uvl11)*fach) + uvl11
                  tn12 = ((tnu12-tnl12)*fach) + tnl12
                  uz12 = ((uzu12-uzl12)*fach) + uzl12
                  um12 = ((umu12-uml12)*fach) + uml12
                  uv12 = ((uvu12-uvl12)*fach) + uvl12
                  tn21 = ((tnu21-tnl21)*fach) + tnl21
                  uz21 = ((uzu21-uzl21)*fach) + uzl21
                  um21 = ((umu21-uml21)*fach) + uml21
                  uv21 = ((uvu21-uvl21)*fach) + uvl21
                  tn22 = ((tnu22-tnl22)*fach) + tnl22
                  uz22 = ((uzu22-uzl22)*fach) + uzl22
                  um22 = ((umu22-uml22)*fach) + uml22
                  uv22 = ((uvu22-uvl22)*fach) + uvl22
!                 te11 = ((teu11-tel11)*fach) + tel11
!                 te12 = ((teu12-tel12)*fach) + tel12
!                 te21 = ((teu21-tel21)*fach) + tel21
!                 te22 = ((teu22-tel22)*fach) + tel22
              ENDIF

              o11 = (((ou11-ol11)*fach)+ol11)
              o12 = (((ou12-ol12)*fach)+ol12)
              o21 = (((ou21-ol21)*fach)+ol21)
              o22 = (((ou22-ol22)*fach)+ol22)
              oo11 =(((oou11-ool11)*fach)+ool11)
              oo12 = (((oou12-ool12)*fach)+ool12)
              oo21 = (((oou21-ool21)*fach)+ool21)
              oo22 = (((oou22-ool22)*fach)+ool22)
              dnn11 =(((dnnu11-dnnl11)*fach)+dnnl11)
              dnn12 = (((dnnu12-dnnl12)*fach)+dnnl12)
              dnn21 = (((dnnu21-dnnl21)*fach)+dnnl21)
              dnn22 = (((dnnu22-dnnl22)*fach)+dnnl22)

              top11 = (((topu11-topl11)*fach)+topl11)
              top12 = (((topu12-topl12)*fach)+topl12)
              top21 = (((topu21-topl21)*fach)+topl21)
              top22 = (((topu22-topl22)*fach)+topl22)
              tnop11 = (((tnopu11-tnopl11)*fach)+tnopl11)
              tnop12 = (((tnopu12-tnopl12)*fach)+tnopl12)
              tnop21 = (((tnopu21-tnopl21)*fach)+tnopl21)
              tnop22 = (((tnopu22-tnopl22)*fach)+tnopl22)
              to2p11 = (((to2pu11-to2pl11)*fach)+to2pl11)
              to2p12 = (((to2pu12-to2pl12)*fach)+to2pl12)
              to2p21 = (((to2pu21-to2pl21)*fach)+to2pl21)
              to2p22 = (((to2pu22-to2pl22)*fach)+to2pl22)

          if(sw_External_model_provides_NO_N4S_densities) then
              dNO11 =(((dNOu11-dNOl11)*fach)+dNOl11)
              dNO12 = (((dNOu12-dNOl12)*fach)+dNOl12)
              dNO21 = (((dNOu21-dNOl21)*fach)+dNOl21)
              dNO22 = (((dNOu22-dNOl22)*fach)+dNOl22)
              dN4S11 =(((dN4Su11-dN4Sl11)*fach)+dN4Sl11)
              dN4S12 = (((dN4Su12-dN4Sl12)*fach)+dN4Sl12)
              dN4S21 = (((dN4Su21-dN4Sl21)*fach)+dN4Sl21)
              dN4S22 = (((dN4Su22-dN4Sl22)*fach)+dN4Sl22)
              dN2D11 =(((dN2Du11-dN2Dl11)*fach)+dN2Dl11)
              dN2D12 = (((dN2Du12-dN2Dl12)*fach)+dN2Dl12)
              dN2D21 = (((dN2Du21-dN2Dl21)*fach)+dN2Dl21)
              dN2D22 = (((dN2Du22-dN2Dl22)*fach)+dN2Dl22)
           endif
          !g
!              if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'YYYYYY',ou11,ol11,fach
!               write(6,*) 'yabs 1',o11,fach
!              endif
              if(o11 > small_power) then
                  o11=10**o11
              else
                  o11=small_number
              endif
              if(o12 > small_power) then
                  o12=10**o12
              else
                  o12=small_number
              endif
              if(o21 > small_power) then
                  o21=10**o21
              else
                  o21=small_number
              endif
              if(o22 > small_power) then
                  o22=10**o22
              else
                  o22=small_number
              endif
!              if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'yabs ',o11
!              endif
          !g
              if(oo11 > small_power) then
                  oo11=10**oo11
              else
                  oo11=small_number
              endif
              if(oo12 > small_power) then
                  oo12=10**oo12
              else
                  oo12=small_number
              endif
              if(oo21 > small_power) then
                  oo21=10**oo21
              else
                  oo21=small_number
              endif
              if(oo22 > small_power) then
                  oo22=10**oo22
              else
                  oo22=small_number
              endif
          !g
              if(dnn11 > small_power) then
                  dnn11=10**dnn11
              else
                  dnn11=small_number
              endif
              if(dnn12 > small_power) then
                  dnn12=10**dnn12
              else
                  dnn12=small_number
              endif
              if(dnn21 > small_power) then
                  dnn21=10**dnn21
              else
                  dnn21=small_number
              endif
!if ( sw_debug ) print *, 'BEFORE dnn22',mp,lp,dnn22,small_power
              if(dnn22 > small_power) then
                  dnn22=10**dnn22
!if ( sw_debug ) print *, 'AFTER check dnn22',mp,lp,dnn22,small_number
              else
                  dnn22=small_number
              endif
          !g
          if(sw_External_model_provides_NO_N4S_densities) then
              if(dNO11 > small_power) then
                  dNO11=10**dNO11
              else
                  dNO11=small_number
              endif
              if(dNO12 > small_power) then
                  dNO12=10**dNO12
              else
                  dNO12=small_number
              endif
              if(dNO21 > small_power) then
                  dNO21=10**dNO21
              else
                  dNO21=small_number
              endif
              if(dNO22 > small_power) then
                  dNO22=10**dNO22
              else
                  dNO22=small_number
              endif
          !g
              if(dN4S11 > small_power) then
                  dN4S11=10**dN4S11
              else
                  dN4S11=small_number
              endif
              if(dN4S12 > small_power) then
                  dN4S12=10**dN4S12
              else
                  dN4S12=small_number
              endif
              if(dN4S21 > small_power) then
                  dN4S21=10**dN4S21
              else
                  dN4S21=small_number
              endif
              if(dN4S22 > small_power) then
                  dN4S22=10**dN4S22
              else
                  dN4S22=small_number
              endif
          !g
              if(dN2D11 > small_power) then
                  dN2D11=10**dN2D11
              else
                  dN2D11=small_number
              endif
              if(dN2D12 > small_power) then
                  dN2D12=10**dN2D12
              else
                  dN2D12=small_number
              endif
              if(dN2D21 > small_power) then
                  dN2D21=10**dN2D21
              else
                  dN2D21=small_number
              endif
              if(dN2D22 > small_power) then
                  dN2D22=10**dN2D22
              else
                  dN2D22=small_number
              endif

           endif
          !g

              IF ( ispecial == 0 ) THEN
                  faclon = (glond2-geo_grid_longitudes_degrees(ilon1)) / &
                           (geo_grid_longitudes_degrees(ilon2)-geo_grid_longitudes_degrees(ilon1))
              ELSEIF ( ispecial == 1 ) THEN
                  faclon = (glond2-geo_grid_longitudes_degrees(ilon1)) &
                  /(geo_grid_longitudes_degrees(ilon2)+360.-geo_grid_longitudes_degrees(ilon1))
              ENDIF
              tn2 = ((tn22-tn21)*faclon) + tn21
              tn1 = ((tn12-tn11)*faclon) + tn11
!             te2 = ((te22-te21)*faclon) + te21
!             te1 = ((te12-te11)*faclon) + te11
              uz2 = ((uz22-uz21)*faclon) + uz21
              uz1 = ((uz12-uz11)*faclon) + uz11
              um2 = ((um22-um21)*faclon) + um21
              um1 = ((um12-um11)*faclon) + um11
              uv2 = ((uv22-uv21)*faclon) + uv21
              uv1 = ((uv12-uv11)*faclon) + uv11
              do2 = ((o22-o21)*faclon) + o21
              do1 = ((o12-o11)*faclon) + o11
              doo2 = ((oo22-oo21)*faclon) + oo21
              doo1 = ((oo12-oo11)*faclon) + oo11
              dnn2 = ((dnn22-dnn21)*faclon) + dnn21
              dnn1 = ((dnn12-dnn11)*faclon) + dnn11

              top2 = ((top22-top21)*faclon) + top21
              top1 = ((top12-top11)*faclon) + top11
              tnop2 = ((tnop22-tnop21)*faclon) + tnop21
              tnop1 = ((tnop12-tnop11)*faclon) + tnop11
              to2p2 = ((to2p22-to2p21)*faclon) + to2p21
              to2p1 = ((to2p12-to2p11)*faclon) + to2p11

          if(sw_External_model_provides_NO_N4S_densities) then
              dNO2 = ((dNO22-dNO21)*faclon) + dNO21
              dNO1 = ((dNO12-dNO11)*faclon) + dNO11
              dN4S2 = ((dN4S22-dN4S21)*faclon) + dN4S21
              dN4S1 = ((dN4S12-dN4S11)*faclon) + dN4S11
              dN2D2 = ((dN2D22-dN2D21)*faclon) + dN2D21
              dN2D1 = ((dN2D12-dN2D11)*faclon) + dN2D11
          endif

              faclat = (GLAt(i)-geo_grid_latitudes_degrees(ilat1)) / &
                       (geo_grid_latitudes_degrees(ilat2)-geo_grid_latitudes_degrees(ilat1))
          ! write(6,*) 'faclat  ',faclat
              TN(i) = ((tn2-tn1)*faclat) + tn1
!             Te_dum(i) = ((te2-te1)*faclat) + te1
              !if(te_dum(i) > 5000.) then
              ! write(6,*) 'TE_DUM ',i,mp,lp,te_dum(i)
              !   te_dum(i) = 5000.
              !endif
              uz(i) = ((uz2-uz1)*faclat) + uz1
              um(i) = ((um2-um1)*faclat) + um1
              uv(i) = ((uv2-uv1)*faclat) + uv1
              O(i) = ((do2-do1)*faclat) + do1
              if(o(i) < small_number) o(i)=small_number
              O2(i) = ((doo2-doo1)*faclat) + doo1
              if(o2(i) < small_number) o2(i)=small_number
              N2(i) = ((dnn2-dnn1)*faclat) + dnn1
              if(n2(i) < small_number) n2(i)=small_number
          if(sw_External_model_provides_NO_N4S_densities) then
              NO(i) = ((dNO2-dNO1)*faclat) + dNO1
              if(NO(i) < small_number) NO(i)=small_number
              N4S(i) = ((dN4S2-dN4S1)*faclat) + dN4S1
              if(N4S(i) < small_number) N4S(i)=small_number
              N2D(i) = ((dN2D2-dN2D1)*faclat) + dN2D1
              if(N2D(i) < small_number) N2D(i)=small_number
          endif

              tiegcm_OP_on_tube(i) = ((top2-top1)*faclat) + top1
              tiegcm_NOP_on_tube(i) = ((tnop2-tnop1)*faclat) + tnop1
              tiegcm_O2P_on_tube(i) = ((to2p2-to2p1)*faclat) + to2p1
          !g
          !g THe above does not let any of the densities get lower
          !g than 'small_number'.  This can be a problem at low solar activity
          !g at the top of the larger flux-tubes.......
          !g
          ! write(199,*) i,mp,lp,faclon,faclat,fach
          900 ENDDO
          do i=in(mp,lp),is(mp,lp)
              TN_plasma_input_3d(i,mp) = tn(i)
              TN_plasma_input_3d(i,mp) = tn(i)
              if( TN_plasma_input_3d(i,mp) .lt. 0.) then
                  write(6,*) 'TN_plasma_input_3d lt 0)',i,mp,lp,tn(i)
                  TN_plasma_input_3d(i,mp) = 50.
              endif
              O_plasma_input_3d(i,mp) = o(i)
              O2_plasma_input_3d(i,mp) = o2(i)
              N2_plasma_input_3d(i,mp) = n2(i)
          if(sw_External_model_provides_NO_N4S_densities) then
              NO_plasma_input_3d(i,mp) = NO(i)
              N4S_plasma_input_3d(i,mp) = N4S(i)
              N2D_plasma_input_3d(i,mp) = N2D(i)
          endif
              um_plasma_input_3d(i,mp) = um(i)
              uz_plasma_input_3d(i,mp) = uz(i)
              uv_plasma_input_3d(i,mp) = uv(i)
!             te_plasma_input_3d(i,mp) = te_dum(i)

              tiegcm_OP_input_3d(i,mp) = tiegcm_OP_on_tube(i)
              tiegcm_NOP_input_3d(i,mp) = tiegcm_NOP_on_tube(i)
              tiegcm_O2P_input_3d(i,mp) = tiegcm_O2P_on_tube(i)

           enddo
      !g
      !g  end of the big flux tubes loop....
      !g
      enddo   !do lp = 1 , nlp
   enddo      !do mp = 1 , nmp
!       write(6,*) '***** done with interface to plasma ***** '

  sw_1st_call_int_fixed_ht = .FALSE.
  RETURN


  end SUBROUTINE INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE


END MODULE IPE_Neutrals_Class
