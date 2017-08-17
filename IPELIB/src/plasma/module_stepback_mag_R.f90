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
!nm20130201: separated into the new module subroutine!
!nm20130201:      SUBROUTINE stepback_mag (mp,lp &


      MODULE module_stepback_mag_R
        PRIVATE
        PUBLIC :: stepback_mag_R
      CONTAINS
      SUBROUTINE stepback_mag_R (utime,mp,lp,phi_t0,theta_t0,r0_apex)
      USE module_precision
      USE module_IPE_dimension,ONLY: NLP
      USE module_FIELD_LINE_GRID_MKS,ONLY: mlon_rad,plasma_grid_Z,JMIN_IN,JMAX_IS,ht90,plasma_grid_GL,plasma_grid_3d,east,north,up,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,VEXBup,minAltitude,maxAltitude, VEXBe
      USE module_physical_constants,ONLY: earth_radius,rtd,pi
      USE module_input_parameters,ONLY: time_step,sw_exb_up,sw_debug,start_time,lpmin_perp_trans
      IMPLICIT NONE
! INPUT
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec), INTENT(IN) :: mp    !mag-lon index
      INTEGER (KIND=int_prec), INTENT(IN) :: lp    !mag-lat index
! OUTPUT
      REAL(KIND=real_prec), INTENT(OUT) :: phi_t0(2)   !magnetic longitude,phi[rad] at T0
      REAL(KIND=real_prec), INTENT(OUT) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      REAL(KIND=real_prec), INTENT(OUT) :: r0_apex     ![meter]
! local
      REAL   (KIND=real_prec) :: phi_t1     !magnetic longitude,phi at T1
      REAL   (KIND=real_prec) :: theta_t1(2)!magnetic latitude,theta at T1
      INTEGER(KIND=int_prec ) :: midpoint
      REAL   (KIND=real_prec) :: r,r_apex,sin2theta,sintheta,theta !meter
      INTEGER(KIND=int_prec ) :: ihem                              !1:NH; 2:SH
      REAL   (KIND=real_prec) :: GLON_deg, LT_SEC

      phi_t1 = mlon_rad(mp)
      theta_t1(1) = plasma_grid_GL( JMIN_IN(lp),lp ) !NH
      theta_t1(2) = plasma_grid_GL( JMAX_IS(lp),lp ) !SH

      r = earth_radius + ht90 ![m]
      midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
      r_apex = earth_radius + plasma_grid_Z(midpoint,lp) ![m]
if(sw_debug) print *,'sub-StR:',lp,mp,r,r_apex,plasma_grid_Z(midpoint,lp)

!note: for the moment, Ed1/B is calculated only in NH, assuming that the flux tube is moving with the same velocity between N/SH.
      which_hemisphere: DO ihem=1,1 !ihem_max
        
        IF ( sw_exb_up<=1 ) THEN 
!         (1) WACCM E empirical model moved to get_efield90km

!         dbg20120301:temp solution: make sure flux tube does not go beyond the sim region...
          if ( lp==1.or.lp==NLP ) then
            VEXBup(lp,mp) = 0.0
          endif

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
          if ( mp==1.and.lp==lpmin_perp_trans.and.MOD( (utime-start_time),900 )==0 )    CALL read_vexb ( utime,lp,mp )

        END IF !ELSE IF ( sw_exb_up==3 ) THEN 

!if(sw_debug)&
! print *,'sub-StR:',ihem,lp,mp,'v_exb_apex[m/s]',VEXBup(lp,mp)  ,utime

        r0_apex = r_apex - VEXBup(lp,mp) * time_step
if(sw_debug)&
& print *,'sub-StR:',r_apex,' r0 apex[m/s]',r0_apex

if ( r0_apex<(minAltitude+earth_radius) ) then
   print *,'!r0_apex too small!',r0_apex,VEXBup(lp,mp),lp,mp, (minAltitude+earth_radius)
   r0_apex = minAltitude+earth_radius
else if ( r0_apex>(maxAltitude+earth_radius) ) then
   print *,'!r0_apex too big!',r0_apex, VEXBup(lp,mp),lp,mp, (maxAltitude+earth_radius)
   r0_apex = maxAltitude+earth_radius
end if
!dbg20120301:
        sin2theta = r/r0_apex
        sintheta = SQRT( sin2theta )
        theta    = ASIN ( sintheta )
if(sw_debug) print *,'SIN',sin2theta,sintheta,theta,' mlat R0[deg]', (90.-  theta*180./pi)

        if ( ihem==1 ) then
           theta_t0(ihem) = pi*0.50 - ACOS ( sintheta )
        else if ( ihem==2 ) then
           theta_t0(ihem) = pi*0.50 + ACOS ( sintheta )
        end if

!temporary solution...
!        phi_t0(ihem)   = phi_t1
!nm20130201
        phi_t0(ihem)   = phi_t1 - ( VEXBe(lp,mp) * time_step ) / r_apex
        IF ( phi_t0(ihem)>=pi*2.0 ) THEN
          phi_t0(ihem) = phi_t0(ihem) - pi*2.0
        ELSE IF ( phi_t0(ihem)< 0.0    ) THEN
          phi_t0(ihem) = phi_t0(ihem) + pi*2.0
        END IF

END DO      which_hemisphere !: DO ihem=1,ihem_max

ihem=1 !only
if(sw_debug) print *,lp,mp, 'Z(mp,lp)',plasma_grid_Z(midpoint,lp), (plasma_grid_Z(midpoint,lp)+earth_radius)
if(sw_debug) print *, 'mlatN(mp,lp)',90.-plasma_grid_GL( JMIN_IN(lp),lp )*180./pi !NH
sin2theta = r/( earth_radius+plasma_grid_Z(midpoint,lp) )
sintheta = SQRT( sin2theta )
theta_t1(ihem)    = ASIN ( sintheta )
if(sw_debug) print *,'DIPOLE',sin2theta,sintheta,theta_t1(ihem),' mlat NH[deg]', (90.-  theta_t1(ihem)*180./pi)


if(sw_debug) print "('sub-StR:T1: phi=',F12.6,' theta=',F12.6,' r=',E13.5)",phi_t1*rtd,      (90.-theta_t1(ihem)*rtd), r_apex
if(sw_debug) print "('sub-StR:T0: phi=',F12.6,' theta=',F12.6,' r=',E13.5)",phi_t0(ihem)*rtd,(90.-theta_t0(ihem)*rtd), r0_apex
if(sw_debug) print "(3E12.4)", (theta_t1(ihem)-theta_t0(ihem))*rtd, (r_apex-r0_apex), (VEXBup(lp,mp) * time_step)

      END SUBROUTINE stepback_mag_R
      END MODULE module_stepback_mag_R
