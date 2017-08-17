! calculate P value (L-shell) using a dipole assumption
!NOTE: in FLIP, PCO is only used in setting up the rough plasmasphere H+ initial profiles (See PROFIN). It does not have to be accurate.
SUBROUTINE Get_Pvalue_Dipole ( R_meter, Theta_rad, Pvalue_out )   
      USE module_precision
      USE module_physical_constants,ONLY: earth_radius
      IMPLICIT NONE
      REAL(KIND=real_prec), INTENT(IN)  ::  R_meter  != earth_radius + height [meter]
      REAL(KIND=real_prec), INTENT(IN)  ::  Theta_rad  != magnetic co-latitude[rad] at the northern foot point(JMIN)
      REAL(KIND=real_prec), INTENT(OUT) ::  Pvalue_out
!---
      REAL(KIND=real_prec) :: sinTheta  !sin(theta)
!---

      sinTheta = SIN( Theta_rad )
      Pvalue_out = R_meter / ( earth_radius * sinTheta * sinTheta )

END SUBROUTINE Get_Pvalue_Dipole   
