MODULE IPE_Forcing_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines

IMPLICIT NONE

  TYPE IPE_Forcing

    ! forcing % n_time_levels = f107_kp_data_size    = f107_kp_size
    INTEGER                 :: n_time_levels
    REAL(prec)              :: dt
    REAL(prec), ALLOCATABLE :: time(:)
    REAL(prec)              :: current_time
    INTEGER                 :: current_index
    !
    REAL(prec), ALLOCATABLE :: f107(:)
    INTEGER, ALLOCATABLE    :: f107_flag(:)
    REAL(prec), ALLOCATABLE :: f107_81day_avg(:)
    REAL(prec), ALLOCATABLE :: kp(:)
    INTEGER, ALLOCATABLE    :: kp_flag(:)
    REAL(prec), ALLOCATABLE :: kp_1day_avg(:)
    REAL(prec), ALLOCATABLE :: ap(:)
    REAL(prec), ALLOCATABLE :: ap_1day_avg(:)
    REAL(prec), ALLOCATABLE :: nhemi_power(:)
    INTEGER, ALLOCATABLE    :: nhemi_power_index(:)
    REAL(prec), ALLOCATABLE :: shemi_power(:)
    INTEGER, ALLOCATABLE    :: shemi_power_index(:)
    ! Solar wind drivers
    REAL(prec), ALLOCATABLE :: solarwind_dBdt(:)
    REAL(prec), ALLOCATABLE :: solarwind_angle(:)
    REAL(prec), ALLOCATABLE :: solarwind_velocity(:)
    REAL(prec), ALLOCATABLE :: solarwind_Bz(:)

    ! Tiros
    REAL(prec), ALLOCATABLE :: emaps(:,:,:)
    REAL(prec), ALLOCATABLE :: cmaps(:,:,:)
    REAL(prec), ALLOCATABLE :: djspectra(:,:)

    !
    REAL(prec) :: sun_longitude

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Forcing
      PROCEDURE :: Trash => Trash_IPE_Forcing


      PROCEDURE :: Read_F107KP_IPE_Forcing
      PROCEDURE :: Read_Tiros_IPE_Forcing
 
      PROCEDURE :: Estimate_AP_from_KP

  END TYPE IPE_Forcing

CONTAINS

  SUBROUTINE Build_IPE_Forcing( forcing, dt, n_time_levels )
    IMPLICIT NONE
    CLASS( IPE_Forcing ), INTENT(out) :: forcing
    REAL(prec)                        :: dt
    INTEGER, INTENT(in)               :: n_time_levels
    

      forcing % n_time_levels = n_time_levels 
      forcing % dt            = dt
      forcing % current_time  = 0.0_prec
      forcing % current_index = 1

      ALLOCATE( forcing % time(1:n_time_levels), &
                forcing % f107(1:n_time_levels), &
                forcing % f107_flag(1:n_time_levels), &
                forcing % f107_81day_avg(1:n_time_levels), &
                forcing % kp(1:n_time_levels), &
                forcing % kp_flag(1:n_time_levels), &
                forcing % kp_1day_avg(1:n_time_levels), &
                forcing % ap(1:n_time_levels), &
                forcing % ap_1day_avg(1:n_time_levels), &
                forcing % nhemi_power(1:n_time_levels), &
                forcing % nhemi_power_index(1:n_time_levels), &
                forcing % shemi_power(1:n_time_levels), &
                forcing % shemi_power_index(1:n_time_levels), &
                forcing % solarwind_dBdt(1:n_time_levels), &
                forcing % solarwind_angle(1:n_time_levels), &
                forcing % solarwind_velocity(1:n_time_levels), &
                forcing % solarwind_Bz(1:n_time_levels), &
                forcing % emaps(1:maps_ipe_size(1),1:maps_ipe_size(2),1:maps_ipe_size(3)), &
                forcing % cmaps(1:maps_ipe_size(1),1:maps_ipe_size(2),1:maps_ipe_size(3)), &
                forcing % djspectra(1:n_flux_ipe,1:n_bands_ipe)  )

      ! Default settings
      forcing % f107              = 120.0_prec
      forcing % f107_flag         = 1
      forcing % f107_81day_avg    = 120.0_prec
      forcing % kp                = 3.0_prec
      forcing % kp_flag           = 3
      forcing % kp_1day_avg       = 3.0_prec
      forcing % ap                = 0.0_prec
      forcing % ap_1day_avg       = 0.0_prec
      forcing % nhemi_power       = 21.0_prec
      forcing % nhemi_power_index = 6
      forcing % shemi_power       = 21.0_prec
      forcing % shemi_power_index = 6

      forcing % solarwind_dBdt     = 1.0_prec
      forcing % solarwind_angle    = 0.0_prec
      forcing % solarwind_velocity = 400.0_prec
      forcing % solarwind_Bz       = 1.0_prec

      forcing % emaps     = 0.0_prec
      forcing % cmaps     = 0.0_prec
      forcing % djspectra = 0.0_prec

      CALL forcing % Estimate_AP_from_KP( )

  END SUBROUTINE Build_IPE_Forcing
!
  SUBROUTINE Trash_IPE_Forcing( forcing )
    IMPLICIT NONE
    CLASS( IPE_Forcing ), INTENT(inout) :: forcing

    DEALLOCATE( forcing % time, &
                forcing % f107, &
                forcing % f107_flag, &
                forcing % f107_81day_avg, &
                forcing % kp, &
                forcing % kp_flag, &
                forcing % kp_1day_avg, &
                forcing % ap, &
                forcing % ap_1day_avg, &
                forcing % nhemi_power, &
                forcing % nhemi_power_index, &
                forcing % shemi_power, &
                forcing % shemi_power_index, &
                forcing % solarwind_dBdt, &
                forcing % solarwind_angle, &
                forcing % solarwind_velocity, &
                forcing % solarwind_Bz )


  END SUBROUTINE Trash_IPE_Forcing
!
  SUBROUTINE Read_F107KP_IPE_Forcing( forcing, skip_size, filename )
    IMPLICIT NONE
    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    INTEGER, INTENT(in)                 :: skip_size 
    CHARACTER(*), INTENT(in)            :: filename
    ! Local
    INTEGER :: fUnit
    INTEGER :: i, read_in_size
    CHARACTER(20) :: date_work
    REAL(prec)    :: f107_work, kp_work, f107d_work, kpa_work
    REAL(prec)    :: hp_work, swbt_work, swang_work, swvel_work, swbz_work
    INTEGER       :: f107_flag_work, kp_flag_work, hpi_work
    
      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = TRIM(filename), &
            FORM   = 'FORMATTED', &
            ACTION = 'READ', &
            STATUS = 'OLD' )

      ! Skip over the header information
      DO i = 1, 5
        READ(fUnit, *)    
      END DO

      DO i = 1, skip_size

        READ(fUnit, *) date_work, f107_work, kp_work, f107_flag_work, &
                       kp_flag_work, f107d_work, kpa_work, hp_work, hpi_work, &
                       swbt_work, swang_work, swvel_work, swbz_work

      END DO

      read_in_size = forcing % n_time_levels - skip_size

      DO i = 1, MIN(read_in_size, forcing % n_time_levels)

        READ(fUnit, *) date_work, &
                       forcing % f107(i), &
                       forcing % kp(i), &
                       forcing % f107_flag(i), &
                       forcing % kp_flag(i), &
                       forcing % f107_81day_avg(i), &
                       forcing % kp_1day_avg(i), &
                       forcing % nhemi_power(i), &
                       forcing % nhemi_power_index(i), &
                       forcing % shemi_power(i), &
                       forcing % shemi_power_index(i), &
                       forcing % solarwind_dBdt(i), &
                       forcing % solarwind_angle(i), &
                       forcing % solarwind_velocity(i), &
                       forcing % solarwind_Bz(i)

      END DO

      CLOSE(fUnit)

      DO i = read_in_size + 1, forcing % n_time_levels

          forcing % f107(i)               = forcing % f107(read_in_size)
          forcing % f107_81day_avg(i)     = forcing % f107_81day_avg(read_in_size)
          forcing % kp(i)                 = forcing % kp(read_in_size)
          forcing % kp_1day_avg(i)        = forcing % kp_1day_avg(read_in_size)
          forcing % nhemi_power(i)        = forcing % nhemi_power(read_in_size)
          forcing % nhemi_power_index(i)  = forcing % nhemi_power_index(read_in_size)
          forcing % shemi_power(i)        = forcing % shemi_power(read_in_size)
          forcing % shemi_power_index(i)  = forcing % shemi_power_index(read_in_size)
          forcing % solarwind_dbdt(i)     = forcing % solarwind_dbdt(read_in_size)
          forcing % solarwind_angle(i)    = forcing % solarwind_angle(read_in_size)
          forcing % solarwind_velocity(i) = forcing % solarwind_velocity(read_in_size)
          forcing % solarwind_Bz(i)       = forcing % solarwind_Bz(read_in_size)

      END DO

      CALL forcing % Estimate_AP_from_KP( )
 
  END SUBROUTINE Read_F107KP_IPE_Forcing
!
  SUBROUTINE Read_Tiros_IPE_Forcing( forcing )
    IMPLICIT NONE
    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    ! Local
    INTEGER :: fUnit, iBand

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = './ionprof', &
            FORM = 'FORMATTED', &
            STATUS = 'OLD', &
            ACTION = 'READ' )

      READ (fUnit,'(1x,6E13.6)') forcing % emaps
      READ (fUnit,'(1x,6E13.6)') forcing % cmaps

      CLOSE(fUnit)

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = './tiros_spectra', &
            FORM = 'FORMATTED', &
            STATUS = 'OLD', &
            ACTION = 'READ' )

      READ(fUnit,*)
      READ(fUnit,*)
      READ(fUnit,*)

      DO iBand = 1, n_bands_ipe

         READ(fUnit,*)
         READ(fUnit,*)
         READ(fUnit,'(1X,5E10.4)') forcing % djspectra(1:n_flux_ipe,iBand)
         READ(fUnit,*)

      ENDDO

      CLOSE(fUnit)

  END SUBROUTINE Read_Tiros_IPE_Forcing
!
  SUBROUTINE Estimate_AP_from_KP( forcing )
    IMPLICIT NONE
    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    ! Local
    REAL(prec) :: lookup, remainder
    INTEGER    :: i
    INTEGER, PARAMETER  :: table(1:29) = (/  0,   2,   3, & ! 0-0.67
                                             4,   5,   6, & ! 1-1.67
                                             7,   9,  12, & ! 2-2.67
                                             15,  18,  22, & ! 3-3.67
                                             27,  32,  39, & ! 4-4.67
                                             48,  56,  67, & ! 5-5.67
                                             80,  94, 111, & ! 6-6.67
                                             132, 154, 179, & ! 7-7.67
                                             207, 236, 300, & ! 8-8.67
                                             400, 999/)       ! 9-dumm
        DO i=1, forcing % n_time_levels

          lookup    = forcing % kp(i)*3.0_prec + 1.0_prec
          remainder = lookup - INT(lookup)
          forcing % ap(i) = (1.0_prec - remainder) * table(INT(lookup)) + remainder *table(INT(lookup)+1)

          lookup    = forcing % kp_1day_avg(i)*3.0_prec + 1.0_prec
          remainder = lookup - INT(lookup)
          forcing % ap_1day_avg(i) = (1.0_prec - remainder) * table(INT(lookup)) + remainder *table(INT(lookup)+1)

        END DO

  END SUBROUTINE Estimate_AP_from_KP
!
END MODULE IPE_Forcing_Class
