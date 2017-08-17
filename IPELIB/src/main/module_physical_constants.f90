MODULE module_physical_constants
      USE module_precision
      IMPLICIT NONE

      REAL(KIND=real_prec), PARAMETER, PUBLIC :: earth_radius = 6.3712E+06 !.. Earth radius [meter]
      REAL(KIND=real_prec), PARAMETER, PUBLIC :: pi = 3.1415926536
      REAL(KIND=real_prec), PARAMETER, PUBLIC :: rtd = 180.0/pi !radian-->degree
      REAL(KIND=real_prec), PARAMETER, PUBLIC :: dtr = pi/180.0 !degree-->radian

      REAL(KIND=real_prec), PARAMETER, PUBLIC :: G0 = 9.80665        !.. strength of the Earth's gravity, nominal average value at the Earth's surface (standard gravity) [m s-2]

!note: these parameters must be somewhere in flip! I should remove the duplicate parameters in the future and contain them all in a single module!!!
      REAL(KIND=real_prec),DIMENSION(6), PARAMETER, PUBLIC :: mass_kg=(/16.,1.,28.,32.,4.,14./)
!MASS for major neutral species: 1:O; 2:O2; 3:N2
      REAL(KIND=real_prec),DIMENSION(3), PARAMETER, PUBLIC :: massn_kg=(/16.,32.,28./)
      REAL(KIND=real_prec), PARAMETER, PUBLIC :: AMU = 1.66E-27     ! Atomic Mass Unit [kg]  

!---numbers with the precision (these can be included here because the numbers are universal...)----
      REAL (KIND=real_prec), PARAMETER, PUBLIC :: zero = 0.0_real_prec

!universal gass constant, as in tucan/ctipe [unit J / mol K???]
      REAL (KIND=real_prec), PARAMETER, PUBLIC :: GSCON = 8.314E+03

END MODULE module_physical_constants
