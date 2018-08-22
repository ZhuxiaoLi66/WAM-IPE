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


END MODULE IPE_Neutrals_Class
