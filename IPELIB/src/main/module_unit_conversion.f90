MODULE module_unit_conversion
      USE module_precision
      IMPLICIT NONE

      REAL(KIND=real_prec), PARAMETER, PUBLIC :: M_TO_CM   = 1.0E+02
      REAL(KIND=real_prec), PARAMETER, PUBLIC :: M3_TO_CM3 = 1.0E-06
      REAL(KIND=real_prec), PARAMETER, PUBLIC :: M_TO_KM   = 1.0E-03


      REAL(KIND=real_prec), PARAMETER, PUBLIC :: eV2J = 1.6021892E-19  ! 1[eV] = 1.6021892E-19 [J]

END MODULE module_unit_conversion
