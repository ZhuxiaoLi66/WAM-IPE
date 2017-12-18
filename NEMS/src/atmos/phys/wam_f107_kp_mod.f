      MODULE wam_f107_kp_mod
!
! VAY 2017-02-21 Orchestraion of OLD/NEW SWPC-drivers
!
      IMPLICIT none

      INTEGER                     :: f107_kp_size, f107_kp_interval
      INTEGER                     :: f107_kp_skip_size
      INTEGER                     :: f107_kp_read_in_size
      INTEGER                     :: kdt_interval
!
      INTEGER                     :: f107_kp_data_size
!
! mark them as f107_wy, kp_wy tp avoid mixture with real f107/kp
! add fixed sets for          f107_fix, f107d_fix, kp_fix
!
      REAL, POINTER, DIMENSION(:) :: f107_wy, kp_wy, kpa_wy, f107d_wy, hp_wy
      REAL, POINTER, DIMENSION(:) :: hpi_wy, swbz_wy, swvel_wy, swbt_wy, swang_wy
      REAL                        :: f107_fix, f107d_fix, kp_fix
      REAL                        :: swpcf107_fix, swpcf107d_fix, swpckp_fix
      REAL                        :: interpolate_weight
      CONTAINS

      SUBROUTINE read_wam_f107_kp_txt

! Subprogram:  read_wam_f107_kp_txt   read-in the inputted f10.7 and kp data. 
!   Prgmmr: Weiyu Yang          Date: 2015-10-19
!
! !revision history log:
!
!  13Apr2017   Houjun Wang, enable handling cases when f107=0 or too small

      CHARACTER*20 :: issuedate, realdate(f107_kp_size)
      CHARACTER*20 :: realdate_work
      INTEGER      :: i, j
      INTEGER      :: f107_flag_work, kp_flag_work
      REAL         :: f107_81d_avg
      REAL         :: f107_work, kp_work, f107d_work, kpa_work, hp_work, hpi_work
      REAL         :: swbz_work, swvel_work, swbt_work, swang_work

! Flags:   0=Forecast, 1=Estimated, 2=Observed
      INTEGER, DIMENSION(f107_kp_size) :: f107_flag, kp_flag 

! Skip the observation data before the forecast starting.
!--------------------------------------------------------
      OPEN(79, FILE='wam_input_f107_kp.txt', FORM='formatted')
      REWIND 79
      READ(79, 1000) issuedate
      DO i = 1, 4
          READ(79, *)    
      END DO

      DO i = 1, f107_kp_skip_size
          READ(79, *) realdate_work, f107_work, kp_work,           &
               f107_flag_work, kp_flag_work, f107d_work, kpa_work, &
               hp_work, hpi_work, swbt_work, swang_work, swvel_work, swbz_work
      END DO

! f107_kp_size is the forecast run length.
!----------------------------------------      
      f107_kp_read_in_size = f107_kp_data_size - f107_kp_skip_size
!
      DO i = 1, MIN(f107_kp_read_in_size, f107_kp_size)
!
!old      DO i = 1, f107_kp_read_in_size
          READ(79, *) realdate(i), f107_wy(i), kp_wy(i),             &
               f107_flag(i), kp_flag(i), f107d_wy(i), kpa_wy(i),     &
               hp_wy(i), hpi_wy(i), swbt_wy(i), swang_wy(i),         &
               swvel_wy(i), swbz_wy(i)
      END DO
      CLOSE(79)
!
! If run time longer than the f107 data, use the latest data to run
! continuously
!
      DO i = f107_kp_read_in_size + 1, f107_kp_size
          f107_wy (i) = f107_wy (f107_kp_read_in_size)
          kp_wy   (i) = kp_wy   (f107_kp_read_in_size)
          hp_wy   (i) = hp_wy   (f107_kp_read_in_size)
          f107d_wy(i) = f107d_wy(f107_kp_read_in_size)
          hpi_wy  (i) = hpi_wy  (f107_kp_read_in_size)
          kpa_wy  (i) = kpa_wy  (f107_kp_read_in_size)
          swbt_wy (i) = swbt_wy (f107_kp_read_in_size)
          swang_wy(i) = swang_wy(f107_kp_read_in_size)
          swvel_wy(i) = swvel_wy(f107_kp_read_in_size)
          swbz_wy (i) = swbz_wy (f107_kp_read_in_size)
      END DO
!For testing.
!------------
!      f107=70.0
!      kp=2.0

! HW 15May2016
! In case of f107 == 0 or too small, 
! use the last-readin 'correct' value or f107_81d_avg 
      do i = 1, f107_kp_size
      if (f107_wy(i) >= 70.0) then
         f107_work = f107_wy(i)
      endif
      enddo

      if (f107_work < 70.0) f107_work = f107_81d_avg

      do i = 1, f107_kp_size
      if (f107_wy(i) < 70.0) then
         f107_wy(i) = f107_work
      endif
      enddo

1000  FORMAT(20x, a20)
1001  FORMAT(20x, f3.0)
      PRINT*, 'issuedate=', issuedate
      PRINT*, 'f107_81d_avg=', f107_81d_avg
      DO i = 1, f107_kp_read_in_size
          PRINT*, i, f107_wy(i), kp_wy(i), f107_flag(i), kp_flag(i)
      END DO
      END SUBROUTINE read_wam_f107_kp_txt

!==========================================================
! Below two service subs to disable "read_wam_f107_kp_txt"
! during model tune-ups
!==========================================================
      SUBROUTINE fix_spweather_data
!=======================================================================
!VAY 2016: This is temporal "substitue" for "read_wam_f107_kp_txt"
!    with fixed Kp and F107 data to work with long-term WAM run
! TO DO "advance_solar" KP-F107 drivers with WAM calendar
!=======================================================================
      swpcf107_fix = 100.
      swpckp_fix   = 1.
      swpcf107d_fix = swpcf107_fix 

      f107_fix = 100.
      kp_fix   = 1.
      f107d_fix = f107_fix
      END SUBROUTINE fix_spweather_data
!
      SUBROUTINE read_spweather_real_data
!=======================================================================
!VAY 2016: This is temporal "substitue" for "read_wam_f107_kp_txt"
!    with fixed Kp and F107 data to work with long-term WAM run
! TO DO "advance_solar" KP-F107 drivers with WAM calendar
!=======================================================================
      f107_fix = 100.
      kp_fix   = 1.
      f107d_fix = f107_fix

      END SUBROUTINE read_spweather_real_data
!

      END MODULE wam_f107_kp_mod
