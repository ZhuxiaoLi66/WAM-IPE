!nm20130201: separated SUBROUTINE find_neighbor_grid_th/R into module subroutine.
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
      MODULE module_perpendicular_transport
        PRIVATE
        PUBLIC :: perpendicular_transport
      CONTAINS

      SUBROUTINE perpendicular_transport ( utime, mp,lp )
      USE module_precision
      USE module_input_parameters,ONLY: sw_debug,mype, sw_th_or_r, sw_perp_transport
      USE module_find_neighbor_grid_TH, ONLY: find_neighbor_grid_TH
      USE module_find_neighbor_grid_R, ONLY: find_neighbor_grid_R
      USE module_stepback_mag_TH, ONLY: stepback_mag_TH
      USE module_stepback_mag_R, ONLY: stepback_mag_R
      IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
!---

      REAL(KIND=real_prec8) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
      REAL(KIND=real_prec8) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      INTEGER (KIND=int_prec),DIMENSION(2,2) :: mp_t0,lp_t0 !1st dim:ihem;2nd dim:i0/i1
      REAL(KIND=real_prec8) :: r0_apex ![meter]

!---

      if(sw_debug) then
         print*,'sub-perpendicular transport:utime',utime,mp,lp
!         print "(' mp=',I4,' lp=',I4,3F10.4)", mp,lp,mlon_rad(mp),plasma_grid_3d(IN,mp)%GL, plasma_grid_3d(IS,mp)%GL       
      end if

! calculate where the flux tube is coming from (semi-lagulangian issue)
      IF ( sw_th_or_R==0 ) THEN
         CALL stepback_mag_TH (utime, mp,lp, phi_t0 , theta_t0, r0_apex )
         if(sw_debug) print *,'stepback_mag TH finished!'

      ELSE IF ( sw_th_or_R==1 ) THEN
         if(sw_debug) print*,'calling sub-stepback_mag_R:utime',utime,mp,lp
         CALL stepback_mag_R ( utime, mp,lp, phi_t0 , theta_t0, r0_apex )
         if(sw_debug) print *,'stepback_mag R finished!'
      END IF !( sw_th_or_R==1 ) THEN


      IF ( sw_th_or_R==0 ) THEN
         CALL find_neighbor_grid_TH ( mp,lp, phi_t0, theta_t0, r0_apex &
              &, mp_t0,lp_t0 )
         if(sw_debug) print *,'find_neighbor_grid_TH finished!'
         sw_perp_transport=2

      ELSE IF ( sw_th_or_R==1 ) THEN
         CALL find_neighbor_grid_R ( mp,lp, phi_t0, theta_t0, r0_apex &
              &, mp_t0,lp_t0 )
         if(sw_debug) print *,'find_neighbor_grid R finished!'
         sw_perp_transport=1
      END IF !( sw_th_or_R==1 ) THEN


! prepare all the parameters along the flux tube by interpolation, in addition to the adiabatic term, compressional term
      CALL interpolate_flux_tube ( mp,lp, phi_t0,theta_t0, r0_apex &
           &, mp_t0,lp_t0 ,&
           & utime) !dbg20141209
      if(sw_debug) print *,'interpolate_flux_tube finished!' 

      END SUBROUTINE perpendicular_transport
      END MODULE module_perpendicular_transport
