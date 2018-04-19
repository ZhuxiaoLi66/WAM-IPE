MODULE IPE_Constants_Dictionary

  USE IPE_Precision

  IMPLICIT NONE

  REAL(prec), PARAMETER :: earth_radius = 6.3712_prec*10.0_prec**6 !.. Earth radius [meter]
  REAL(prec), PARAMETER :: pi = 3.1415926536_prec
  REAL(prec), PARAMETER :: rtd = 180.0_prec/pi !radian-->degree
  REAL(prec), PARAMETER :: dtr = pi/180.0_prec !degree-->radian
  REAL(prec), PARAMETER :: G0 = 9.80665_prec        !.. strength of the Earth's gravity, nominal average value at the Earth's surface (standard gravity) [m s-2]
  REAL(prec), PARAMETER :: mass_kg(1:6)=(/ 16.0_prec, 1.0_prec, 28.0_prec, 32.0_prec, 4.0_prec, 14.0_prec /)
  REAL(prec), PARAMETER :: massn_kg(1:3)=(/ 16.0_prec, 32.0_prec, 28.0_prec /) ! Mass for major neutral species: 1:O; 2:O2; 3:N2
  REAL(prec), PARAMETER :: AMU = 1.66_prec*10.0_prec**(-27)     ! Atomic Mass Unit [kg]  
  REAL(prec), PARAMETER :: GSCON = 8.314_prec*10.0_prec**3  ! universal gas constant, as in tucan/ctipe


  REAL(prec), PARAMETER :: mesh_height_max = 782.0_prec*10.0_prec**(3)

  ! Parameters for controlling the fixed height grid interpolation
  INTEGER, PARAMETER    :: nheights_geo=183
  INTEGER, PARAMETER    :: nlon_geo=90 
  INTEGER, PARAMETER    :: nlat_geo=91
  REAL(prec), PARAMETER :: fillValue_geo = -999999.9999_prec

END MODULE IPE_Constants_Dictionary
