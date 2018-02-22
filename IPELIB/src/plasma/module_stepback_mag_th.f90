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
!nm20140528: IF statement sw_r_or_th is added

MODULE module_stepback_mag_TH


  IMPLICIT NONE

  INTEGER ,private :: minlp,minmp,maxlp,maxmp
  REAL, private :: mindIF=+10.
  REAL, private :: maxdIF=-10.
  PRIVATE
  PUBLIC :: stepback_mag_TH

CONTAINS

  SUBROUTINE stepback_mag_TH (utime_local,mp,lp,phi_t0,theta_t0,r0_apex)
    USE module_precision
    USE module_IPE_dimension,ONLY: NLP
    USE module_FIELD_LINE_GRID_MKS,ONLY: mlon_rad,plasma_grid_Z,JMIN_IN,JMAX_IS,ht90,plasma_grid_GL,plasma_grid_3d,east,north,up,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,VEXBup,minAltitude,maxAltitude, VEXBe,VEXBth
    USE module_physical_constants,ONLY: earth_radius,rtd,pi,zero
    USE module_input_parameters,ONLY: time_step,sw_exb_up,sw_debug,start_time,lpmin_perp_trans,fac_exb_up, sw_perp_transport,sw_th_or_R,STOP_time,mype, perp_transport_time_step
    IMPLICIT NONE
! INPUT
    INTEGER (KIND=int_prec), INTENT(IN) :: utime_local !universal time [sec]
    INTEGER (KIND=int_prec), INTENT(IN) :: mp    !mag-lon index
    INTEGER (KIND=int_prec), INTENT(IN) :: lp    !mag-lat index
! OUTPUT
    REAL(KIND=REAL_prec8), INTENT(OUT) :: phi_t0(2)   !magnetic longitude,phi[rad] at T0
    REAL(KIND=REAL_prec8), INTENT(OUT) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
    REAL(KIND=REAL_prec8), INTENT(OUT) :: r0_apex     ![meter]
! local
    REAL   (KIND=REAL_prec8) :: phi_t1     !magnetic longitude,phi at T1
    REAL   (KIND=REAL_prec) :: theta_t1(2)!magnetic latitude,theta at T1
    INTEGER(KIND=int_prec ) :: midpoint
    REAL   (KIND=REAL_prec) :: r,r_apex,sin2theta,sintheta,theta !meter
    REAL   (KIND=REAL_prec) :: coslambda_m !dbg20140527
    INTEGER(KIND=int_prec ) :: ihem                              !1:NH; 2:SH
    REAL   (KIND=REAL_prec) :: GLON_deg, LT_SEC
    REAL   (8) :: r90,rph !meter
    REAL(8) :: sinLambda_m, cos2Lambda_m, sinIm !eq(3.7)
!

    phi_t1 = mlon_rad(mp)
    theta_t1(1) = plasma_grid_GL( JMIN_IN(lp),lp ) !NH
    theta_t1(2) = plasma_grid_GL( JMAX_IS(lp),lp ) !SH

    r90 = earth_radius + ht90 ![m]
    midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
    r_apex = earth_radius + plasma_grid_Z(midpoint,lp) ![m]

    DO ihem=1,1 !ihem_max  !<< Do we need this loop ?

      IF ( sw_exb_up<=1 ) THEN

        VEXBup(lp,mp) = VEXBup(lp,mp) * fac_exb_up

        IF ( lp==1.or.lp==NLP ) THEN
          VEXBup(lp,mp) = zero
        ENDIF


      ELSE IF ( sw_exb_up==3 ) THEN


        GLON_deg = plasma_grid_3d(midpoint,lp,mp,IGLON)*180./pi
        LT_SEC = utime_local + GLON_deg/15.*3600.

        IF ( LT_SEC>=86400.)  LT_SEC=LT_SEC-86400.

        IF ( LT_SEC<     0.)  LT_SEC=LT_SEC+86400.

        CALL supim_EXBV(utime_local,lp,LT_SEC,r_apex,GLON_deg,VEXBup(lp,mp))

      ELSE IF ( sw_exb_up==4 ) THEN

        VEXBup(lp,mp) = zero

      ELSE IF ( sw_exb_up==6 ) THEN

        IF ( lp<=27 ) THEN

          VEXBth(lp,mp) = VEXBth(lp,mp) * 0.5
          VEXBe(lp,mp)  = VEXBe(lp,mp) * 0.5

        ENDIF

      ELSE IF ( sw_exb_up==7 ) THEN

        VEXBe(lp,mp) = zero

      ELSE IF ( sw_exb_up==5 ) THEN

        IF ( mp==1.and.lp==lpmin_perp_trans.and.MOD( (utime_local-start_time),900 )==0 )    CALL read_vexb ( utime_local,lp,mp )

      ENDIF


      IF ( sw_th_or_R==1 ) THEN !R method (gip)

        r0_apex = r_apex - VEXBup(lp,mp) * REAL(perp_transport_time_step)


        IF ( r0_apex<(minAltitude+earth_radius) ) THEN

          r0_apex = minAltitude+earth_radius

        ELSE IF ( r0_apex>(maxAltitude+earth_radius) ) THEN

          r0_apex = maxAltitude+earth_radius

        ENDIF

        coslambda_m = SQRT( r90 / r0_apex )

        IF ( ihem==1 ) THEN       !NH

          theta_t0(ihem) = pi*0.50 - ACOS ( coslambda_m )

        ELSE IF ( ihem==2 ) THEN  !SH

          theta_t0(ihem) = pi*0.50 + ACOS ( coslambda_m )

        ENDIF

        rph = r_apex

      ELSE IF ( sw_th_or_R==0 ) THEN !th method (ctipe/shawn) ! Usually evaluates here

        theta_t0(ihem) = theta_t1(ihem) - ( VEXBth(lp,mp) * REAL(perp_transport_time_step) ) / r90
        rph = r90

        coslambda_m  = COS ( pi*0.50 - theta_t0(ihem) )
        r0_apex = ( earth_radius + ht90 ) *  coslambda_m *  coslambda_m

      ENDIF


      IF ( sw_perp_transport <= 1 ) THEN

        phi_t0(ihem) = phi_t1

      ELSE ! Usually evaluates here

        IF ( SIN(theta_t1(ihem))>zero ) THEN

          phi_t0(ihem) = phi_t1 - ( VEXBe(lp,mp) * REAL(perp_transport_time_step) ) / ( rph * SIN(theta_t1(ihem)) )

        ELSE 

!SMS$IGNORE begin
          PRINT*,mype,'sub-step: !STOP! INVALID sin theta_t1',SIN(theta_t1(ihem)),theta_t1(ihem),mp,lp
!SMS$IGNORE end
          STOP

        ENDIF

        IF ( phi_t0(ihem)<MINVAL(mlon_rad) .or. MAXVAL(mlon_rad)<phi_t0(ihem) ) THEN

!SMS$IGNORE begin
          PRINT*,mype,'sub-stepback_TH: !STOP! time step',perp_transport_time_step,'must be reduced!',MINVAL(mlon_rad),' phi_t0=',phi_t0(ihem),MAXVAL(mlon_rad),phi_t1,VEXBe(lp,mp),mp,lp &
          &,( - ( VEXBe(lp,mp) * REAL(perp_transport_time_step) ) / ( rph * SIN(theta_t1(ihem)) ) )
!SMS$IGNORE end
          STOP

        ENDIF

      ENDIF !( sw_perp_transport <= 1 ) THEN


      IF ( theta_t0(ihem)<zero.OR.theta_t0(ihem)>=pi   ) THEN

!SMS$IGNORE begin
        PRINT*,utime_local,mype,'sub-step: flux tube crosses the pole: !STOP! INVALID theta_t0',theta_t0(ihem),phi_t0(ihem),lp,mp,ihem
!SMS$IGNORE end
        STOP

      ENDIF


    ENDDO


  END SUBROUTINE stepback_mag_TH

END MODULE module_stepback_mag_TH
