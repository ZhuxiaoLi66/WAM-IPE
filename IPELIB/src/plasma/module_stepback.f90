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
!nm20130201: separated into the new module SUBROUTINE!
!nm20130201:      SUBROUTINE stepback_mag (mp,lp &


MODULE module_stepback
  PRIVATE
  PUBLIC :: stepback
CONTAINS
  SUBROUTINE stepback (utime,mp,lp,phi_t0,theta_t0,r0_apex)
    USE module_precision
    USE module_IPE_dimension,ONLY: NLP
    USE module_FIELD_LINE_GRID_MKS,ONLY: mlon_rad,plasma_grid_Z,JMIN_IN,JMAX_IS,ht90,plasma_grid_mag_colat,plasma_grid_3d,east,north,up,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,VEXBup,minAltitude,maxAltitude, VEXBe,VEXBth
    USE module_physical_constants,ONLY: earth_radius,rtd,pi
    USE module_input_parameters,ONLY: time_step,sw_exb_up,start_time,lpmin_perp_trans, perp_transport_time_step
    IMPLICIT NONE
! INPUT
    INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
    INTEGER (KIND=int_prec), INTENT(IN) :: mp    !mag-lon index
    INTEGER (KIND=int_prec), INTENT(IN) :: lp    !mag-lat index
! OUTPUT
    REAL(KIND=REAL_prec8), INTENT(OUT) :: phi_t0(2)   !magnetic longitude,phi[rad] at T0
    REAL(KIND=REAL_prec8), INTENT(OUT) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
    REAL(KIND=REAL_prec8), INTENT(OUT) :: r0_apex     ![meter]
! local
    REAL   (KIND=REAL_prec) :: phi_t1     !magnetic longitude,phi at T1
    REAL   (KIND=REAL_prec) :: theta_t1(2)!magnetic latitude,theta at T1
    INTEGER(KIND=int_prec ) :: midpoint
    REAL   (KIND=REAL_prec) :: r,r_apex,sin2theta,sintheta,theta !meter
    REAL   (KIND=REAL_prec) :: coslambda_m
    INTEGER(KIND=int_prec ) :: ihem                              !1:NH; 2:SH
    REAL   (KIND=REAL_prec) :: GLON_deg, LT_SEC

    phi_t1 = mlon_rad(mp)
    theta_t1(1) = plasma_grid_mag_colat( JMIN_IN(lp),lp ) !NH
    theta_t1(2) = plasma_grid_mag_colat( JMAX_IS(lp),lp ) !SH

    r = earth_radius + ht90 ![m]
    midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
    r_apex = earth_radius + plasma_grid_Z(midpoint,lp) ![m]

!note: for the moment, Ed1/B is calculated only in NH, assuming that the flux tube is moving with the same velocity between N/SH.
    which_hemisphere: DO ihem=1,1 !ihem_max

      IF ( sw_exb_up<=1 ) THEN
!         (1) WACCM E empirical model moved to get_efield90km

!         dbg20120301:temp solution: make sure flux tube DOes not go beyond the sim region...
        IF ( lp==1.or.lp==NLP ) THEN
          VEXBup(lp,mp) = 0.0
        ENDIF

      ELSE IF ( sw_exb_up==2 ) THEN

!         (2) GIP empirical model

      ELSE IF ( sw_exb_up==3 ) THEN

!         (3) SUPIM empirical model:
!             note: becomes zero at R=4000km

        GLON_deg = plasma_grid_3d(midpoint,lp,mp,IGLON)*180./pi
        LT_SEC = utime + GLON_deg/15.*3600.
        IF ( LT_SEC>=86400.)  LT_SEC=LT_SEC-86400.
        IF ( LT_SEC<     0.)  LT_SEC=LT_SEC+86400.
        CALL supim_EXBV(utime,lp,LT_SEC,r_apex,GLON_deg,VEXBup(lp,mp))

      ELSE IF ( sw_exb_up==4 ) THEN

!         (4) zero for debug purpose
        VEXBup(lp,mp) = 0.0   !dbg20111101:v8

      ELSE IF ( sw_exb_up==5 ) THEN
!         (5) read in from a file
        IF ( mp==1.and.lp==lpmin_perp_trans.and.MOD( (utime-start_time),900 )==0 )    CALL read_vexb ( utime,lp,mp )

      ENDIF !ELSE IF ( sw_exb_up==3 ) THEN

       theta_t0(ihem) = theta_t1(ihem) - ( VEXBth(lp,mp) * REAL(perp_transport_time_step) ) / r

        coslambda_m  = COS ( pi*0.50 - theta_t0(ihem) )
        r0_apex = ( earth_radius + ht90 ) /  coslambda_m /  coslambda_m

           phi_t0(ihem)   = phi_t1 - ( VEXBe(lp,mp) * perp_transport_time_step ) / r
      IF ( phi_t0(ihem)>=pi*2.0 ) THEN
        phi_t0(ihem) = phi_t0(ihem) - pi*2.0
      ELSE IF ( phi_t0(ihem)< 0.0    ) THEN
        phi_t0(ihem) = phi_t0(ihem) + pi*2.0
      ENDIF

    ENDDO      which_hemisphere !: DO ihem=1,ihem_max

    ihem=1 !only
    sin2theta = r/( earth_radius+plasma_grid_Z(midpoint,lp) )
    sintheta = SQRT( sin2theta )
    theta_t1(ihem)    = ASIN ( sintheta )

  END SUBROUTINE stepback
END MODULE module_stepback
