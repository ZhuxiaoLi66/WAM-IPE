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
!nm20140528: if statement sw_r_or_th is added 

      MODULE module_stepback_mag_R
IMPLICIT NONE
!dbg20140528: theta1/2 dif diagonastic
INTEGER ,private :: minlp,minmp,maxlp,maxmp
REAL, private :: mindif=+10.
REAL, private :: maxdif=-10.
        PRIVATE
        PUBLIC :: stepback_mag_R
      CONTAINS
      SUBROUTINE stepback_mag_R (utime,mp,lp,phi_t0,theta_t0,r0_apex)
      USE module_precision
      USE module_IPE_dimension,ONLY: NLP
      USE module_FIELD_LINE_GRID_MKS,ONLY: mlon_rad,plasma_grid_Z,JMIN_IN,JMAX_IS,ht90,plasma_grid_GL,plasma_grid_3d,east,north,up,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,VEXBup,minAltitude,maxAltitude, VEXBe,VEXBth
      USE module_physical_constants,ONLY: earth_radius,rtd,pi,zero
      USE module_input_parameters,ONLY: time_step,sw_exb_up,sw_debug,start_time,lpmin_perp_trans,fac_exb_up, sw_perp_transport,sw_th_or_R,stop_time,mype
      IMPLICIT NONE
! INPUT
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec), INTENT(IN) :: mp    !mag-lon index
      INTEGER (KIND=int_prec), INTENT(IN) :: lp    !mag-lat index
! OUTPUT
      REAL(KIND=real_prec8), INTENT(OUT) :: phi_t0(2)   !magnetic longitude,phi[rad] at T0
      REAL(KIND=real_prec8), INTENT(OUT) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      REAL(KIND=real_prec8), INTENT(OUT) :: r0_apex     ![meter]
! local
      REAL   (KIND=real_prec8) :: phi_t1     !magnetic longitude,phi at T1
      REAL   (KIND=real_prec) :: theta_t1(2)!magnetic latitude,theta at T1
      INTEGER(KIND=int_prec ) :: midpoint
      REAL   (KIND=real_prec) :: r,r_apex,sin2theta,sintheta,theta !meter
      REAL   (KIND=real_prec) :: coslambda_m !dbg20140527
      INTEGER(KIND=int_prec ) :: ihem                              !1:NH; 2:SH
      REAL   (KIND=real_prec) :: GLON_deg, LT_SEC
!nm20140528: include theta method
      REAL   (8) :: r90,rph !meter
      REAL(8) :: sinLambda_m, cos2Lambda_m, sinIm !eq(3.7) 
!
      phi_t1 = mlon_rad(mp)
      theta_t1(1) = plasma_grid_GL( JMIN_IN(lp),lp ) !NH
      theta_t1(2) = plasma_grid_GL( JMAX_IS(lp),lp ) !SH

      r90 = earth_radius + ht90 ![m]
      midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
      r_apex = earth_radius + plasma_grid_Z(midpoint,lp) ![m]

!note: for the moment, Ed1/B is calculated only in NH, assuming that the flux tube is moving with the same velocity between N/SH.
      which_hemisphere: DO ihem=1,1 !ihem_max
        
         IF ( sw_exb_up<=1 ) THEN 
!         (1) WACCM E empirical model moved to get_efield90km

!nm20141006: fac_exb_up is introduced for test
            VEXBup(lp,mp) = VEXBup(lp,mp) * fac_exb_up 

!         dbg20120301:temp solution: make sure flux tube does not go beyond the sim region...
            if ( lp==1.or.lp==NLP ) then
               VEXBup(lp,mp) = zero
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

            !(4) zero for debug purpose
            VEXBup(lp,mp) = zero   

         ELSE IF ( sw_exb_up==6 ) THEN 

            !dbg20141210          VEXBth(lp,mp) = zero !dbg20141125
            if ( lp<=27 ) then   !mlat=59.54
               VEXBth(lp,mp) = VEXBth(lp,mp) * 0.5 !dbg20141210
               VEXBe(lp,mp)  = VEXBe(lp,mp) * 0.5  !dbg20141210
            endif

         ELSE IF ( sw_exb_up==7 ) THEN 
            !dbg20141125: 
            VEXBe(lp,mp) = zero !dbg20141125

         ELSE IF ( sw_exb_up==5 ) THEN 
            !         (5) read in from a file
            if ( mp==1.and.lp==lpmin_perp_trans.and.MOD( (utime-start_time),900 )==0 )    CALL read_vexb ( utime,lp,mp )

         END IF !ELSE IF ( sw_exb_up==3 ) THEN 


!nm20140606: here, for polar transport, i will need to introduce internal time step for more accurate semilagrangian estimation.
!nm20140606: internal_time: do t=t0, t1
         if ( sw_th_or_R==1 ) then !R method (gip)
            r0_apex = r_apex - VEXBup(lp,mp) * REAL(time_step)


            if ( r0_apex<(minAltitude+earth_radius) ) then
               print *,'!r0_apex too small!',r0_apex,VEXBup(lp,mp),lp,mp, (minAltitude+earth_radius)
               r0_apex = minAltitude+earth_radius
            else if ( r0_apex>(maxAltitude+earth_radius) ) then
               print *,'!r0_apex too big!',r0_apex, VEXBup(lp,mp),lp,mp, (maxAltitude+earth_radius)
               r0_apex = maxAltitude+earth_radius
            end if
!
!dbg20140527:  define APEX latitude[rad]: eq(3.3) of the imaginary FT(phi0,theta0)
            coslambda_m = SQRT( r90 / r0_apex )

            if ( ihem==1 ) then       !NH
                theta_t0(ihem) = pi*0.50 - ACOS ( coslambda_m )
               !           theta_t0(ihem) = pi*0.50 - ACOS ( sintheta )
            else if ( ihem==2 ) then  !SH
                theta_t0(ihem) = pi*0.50 + ACOS ( coslambda_m )
               !           theta_t0(ihem) = pi*0.50 + ACOS ( sintheta )
           end if

           rph = r_apex

!nm20140528: include theta method
        else if ( sw_th_or_R==0 ) then !th method (ctipe/shawn)

           theta_t0(ihem) = theta_t1(ihem) - ( VEXBth(lp,mp) * REAL(time_step) ) / r90
           rph = r90

!nm20140630: r0_apex needs to correspond to theta_t0. not sure if they are absolutely necessary??? 
           coslambda_m  = COS ( pi*0.50 - theta_t0(ihem) )
           r0_apex = ( earth_radius + ht90 ) *  coslambda_m *  coslambda_m
           
        end if !( sw_th_or_R==0 ) then !th method (ctipe)


!nm20140528 note: VEXBe definition changes depending on sw_th_or_R in eldyn/get_efield90km.f
        if ( sw_perp_transport <= 1 ) then 
           phi_t0(ihem) = phi_t1
        else !if ( sw_perp_transport > 1 ) then  

           if ( SIN(theta_t1(ihem))>zero ) then
              phi_t0(ihem) = phi_t1 - ( VEXBe(lp,mp) * REAL(time_step) ) / ( rph * SIN(theta_t1(ihem)) )
           else 
!SMS$IGNORE begin             
              print*,mype,'sub-step: !STOP! INVALID sin theta_t1',SIN(theta_t1(ihem)),theta_t1(ihem),mp,lp
!SMS$IGNORE end             
              STOP
           end if

!nm20160419 make sure phi_t0 is within mlon_rad range
           if ( phi_t0(ihem)<MINVAL(mlon_rad) .or. MAXVAL(mlon_rad)<phi_t0(ihem) ) then
!SMS$IGNORE begin
              print*,mype,'sub-stepback_R: !STOP! time step',time_step,'must be reduced!',MINVAL(mlon_rad),' phi_t0=',phi_t0(ihem),MAXVAL(mlon_rad),phi_t1,VEXBe(lp,mp),mp,lp &
     &,( - ( VEXBe(lp,mp) * REAL(time_step) ) / ( rph * SIN(theta_t1(ihem)) ) )
!SMS$IGNORE end
              STOP
           end if
        end if !( sw_perp_transport <= 1 ) then 


!nm20140528: make sure that:  0<=theta<=pi, otherwise flux tube goes over the poles to the other side
        IF ( theta_t0(ihem)<zero.OR.theta_t0(ihem)>=pi   ) THEN

!nm20160420: current grid does not need to allow flux tubes to go over the poles, because del lat=1.8 (difference between the maximum mlat=88.2 and the pole) is much bigger than the mlat resolution, ~0.5deg. Thus error stop is set up here for the moment.
!SMS$IGNORE begin
           print*,utime,mype,'sub-step: flux tube crosses the pole: !STOP! INVALID theta_t0',theta_t0(ihem),phi_t0(ihem),lp,mp,ihem
!SMS$IGNORE end
           STOP

!t           IF ( theta_t0(ihem)<zero  ) THEN
!t              theta_t0(ihem) = theta_t0(ihem) * (-1.0)
!t           ELSE IF ( theta_t0(ihem)>=pi   ) THEN
!t              theta_t0(ihem) = pi*2.0 - theta_t0(ihem)
!t           END IF
!t           !+or- 180deg, other side
!t           IF ( phi_t0(ihem) < pi ) THEN
!t              phi_t0(ihem)   = phi_t0(ihem) + pi  
!t           ELSE !IF ( phi_t0(ihem) >= pi ) THEN
!t              phi_t0(ihem)   = phi_t0(ihem) - pi
!t           END IF
!t           !SMS$IGNORE begin
!t           print*,utime,mype,'sub-step: new theta_t0',theta_t0(ihem),phi_t0(ihem),lp,mp,ihem
!t           !SMS$IGNORE end
           
        END IF !( theta_t0(ihem)<zero.OR.theta_t0(ihem)>=pi   ) THEN         


!nm20140606: end do internal_time !: do t=t0, t1; for polar transport, i will need to introduce internal time step for more accurate semilagrangian estimation.


     END DO      which_hemisphere !: DO ihem=1,ihem_max


!ihem=1 !only
!dbg20150309: special debug output to debug pole transport to be read by plt/read_fort9001.pro
!dbg20150309   write(9001,*) mp,lp, phi_t0(ihem),phi_t1 ,theta_t0(ihem), theta_t1(1)

      END SUBROUTINE stepback_mag_R
      END MODULE module_stepback_mag_R
