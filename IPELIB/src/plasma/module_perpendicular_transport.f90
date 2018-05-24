
MODULE module_perpendicular_transport
  PRIVATE
  PUBLIC :: perpendicular_transport
CONTAINS

  SUBROUTINE perpendicular_transport ( utime_local, mp,lp )
    USE module_precision
    USE module_input_parameters,ONLY: mype, sw_convection_footpoint_0_or_apex_1, sw_perp_transport
    USE module_find_neighbor_grid_TH, ONLY: find_neighbor_grid_TH
    USE module_find_neighbor_grid_R, ONLY: find_neighbor_grid_R
    USE module_stepback_mag_TH, ONLY: stepback_mag_TH
    USE module_stepback_mag_R, ONLY: stepback_mag_R
    IMPLICIT NONE
!--- INPUT ---
    INTEGER (KIND=int_prec), INTENT(IN) :: utime_local !universal time [sec]
    INTEGER (KIND=int_prec),INTENT(IN) :: mp
    INTEGER (KIND=int_prec),INTENT(IN) :: lp
!---

    REAL(KIND=REAL_prec8) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
    REAL(KIND=REAL_prec8) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
    INTEGER (KIND=int_prec),DIMENSION(2,2) :: mp_t0,lp_t0 !1st dim:ihem;2nd dim:i0/i1
    REAL(KIND=REAL_prec8) :: r0_apex ![meter]

!---

    ! calculate where the flux tube is coming from (semi-lagulangian issue)
    IF ( sw_convection_footpoint_0_or_apex_1==0 ) THEN

      CALL stepback_mag_TH (utime_local, mp,lp, phi_t0 , theta_t0, r0_apex )

      CALL find_neighbor_grid_TH ( mp,lp, phi_t0, theta_t0, r0_apex, mp_t0,lp_t0 )
      sw_perp_transport=2

    ELSE IF ( sw_convection_footpoint_0_or_apex_1==1 ) THEN

      CALL stepback_mag_R ( utime_local, mp,lp, phi_t0 , theta_t0, r0_apex )

      CALL find_neighbor_grid_R ( mp,lp, phi_t0, theta_t0, r0_apex, mp_t0,lp_t0 )
      sw_perp_transport=1

    ENDIF


    CALL interpolate_flux_tube ( mp,lp, phi_t0,theta_t0, r0_apex, mp_t0,lp_t0, utime_local)

  END SUBROUTINE perpendicular_transport

END MODULE module_perpendicular_transport
