
MODULE module_perpendicular_transport
  PRIVATE
  PUBLIC :: perpendicular_transport
CONTAINS

  SUBROUTINE perpendicular_transport ( utime_local, mp,lp )
    USE module_precision
    USE module_input_parameters,ONLY: mype, sw_perp_transport
    USE module_find_neighbor_grid, ONLY: find_neighbor_grid
    USE module_stepback, ONLY: stepback
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

      CALL stepback ( utime_local, mp,lp, phi_t0 , theta_t0, r0_apex )

      CALL find_neighbor_grid ( mp,lp, phi_t0, theta_t0, r0_apex, mp_t0,lp_t0 )
      sw_perp_transport=1

    CALL interpolate_flux_tube ( mp,lp, phi_t0,theta_t0, r0_apex, mp_t0,lp_t0, utime_local)

  END SUBROUTINE perpendicular_transport

END MODULE module_perpendicular_transport
