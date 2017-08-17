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
      MODULE module_stepback_mag_th
        PRIVATE
        PUBLIC :: stepback_mag_th
      CONTAINS
      SUBROUTINE stepback_mag_th (mp,lp &
!d     &, phi_t1 &
!d     &, theta_t1 &
     &, phi_t0 , theta_t0 )
      USE module_precision
      USE module_eldyn,ONLY: Ed1_90,Ed2_90,coslam_m
      USE module_FIELD_LINE_GRID_MKS,ONLY: Be3,mlon_rad,plasma_grid_GL,JMIN_IN,JMAX_IS,ht90,plasma_grid_Z !,apexE, geographic_coords
      USE module_physical_constants,ONLY: earth_radius,rtd,pi
!      USE cons_module,ONLY: h0 !potential solver reference height[cm] =90km
      USE module_input_parameters,ONLY: time_step,sw_debug,parallelBuild
      IMPLICIT NONE
! INPUT
      INTEGER (KIND=int_prec), INTENT(IN) :: mp  !mag-lon index
      INTEGER (KIND=int_prec), INTENT(IN) :: lp  !mag-lat index
!      REAL(KIND=real_prec), INTENT(IN) :: DT  ! time step [second]      
      REAL(KIND=real_prec) :: phi_t1  !magnetic longitude,phi at T1
      REAL(KIND=real_prec) :: theta_t1(2)!magnetic latitude,theta at T1
! OUTPUT
      REAL(KIND=real_prec), INTENT(OUT) :: phi_t0(2)!magnetic longitude,phi[rad] at T0
      REAL(KIND=real_prec), INTENT(OUT) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
!---local variables---
      REAL(KIND=real_prec) :: v_e(2)  !1:ed2/be3 (4.18) ;2: -ed1/be3 (4.19)
      REAL(KIND=real_prec) :: sinLambda_m  !
      REAL(KIND=real_prec) :: cos2Lambda_m  ! 
      REAL(KIND=real_prec) :: sinIm !eq(3.7) 
! these are geographic
!      TYPE(geographic_coords) :: VEXB !EXB velocity [m s-1] eq(4.17)
      REAL(KIND=real_prec) :: v_e2magsouth !EXB converted to positive:mag-southward
      REAL(KIND=real_prec) :: r !meter
      INTEGER (KIND=int_prec) :: ihem !1:NH; 2:SH
!     INTEGER (KIND=int_prec) :: ift
      REAL(KIND=real_prec) :: r_apex,sintheta,sin2theta,r_apex_dipole,r0_apex
      INTEGER (KIND=int_prec) :: midpoint
!---------------
!dbg
if(parallelBuild) then
  print*,'stepback_mag th does not work in parallel'
  print*,'Stopping in stepback_mag th'
  STOP
endif

phi_t1 = mlon_rad(mp)
theta_t1(1) = plasma_grid_GL( JMIN_IN(lp),lp ) !NH
theta_t1(2) = plasma_grid_GL( JMAX_IS(lp),lp ) !SH


!NOTE: electric fields do not have to be equal at both ends of the flux tubes
!    : transport problem over the poles
      r = earth_radius + ht90 ![m]

      which_hemisphere: DO ihem=1,2 !ihem_max
!        IF ( ihem==1 ) THEN
!          ift = JMIN_IN(lp)
!        ELSE IF ( ihem==2 ) THEN
!          ift = JMAX_IS(lp)
!        END IF

! Ed1/2[V/m] at ( phi_t1(mp), theta_t1(lp) ), Be3[T]
!note: Ed1_90, Ed2_90, Be3 are constant along magnetic field lines!!! 
        v_e(1) =   Ed2_90(1,lp,mp) / Be3(lp,mp) !(4.18) +mag-east(d1?) 
        v_e(2) = - Ed1_90(1,lp,mp) / Be3(lp,mp) !(4.19) +down/equatorward(d2?)
        
!dbg
if(sw_debug) &
& print *,'sub-St:',ihem,'ve2[m/s]',v_e(2),'ed1[mV/m]', Ed1_90(1,lp,mp)*1.0E+3,' Be3[tesla]',Be3(lp,mp) 
!dbg ,' BM-N',plasma_grid_3d( JMIN_IN(lp,mp) ,mp)%BM,' BM-S',plasma_grid_3d( JMAX_IS(lp,mp) ,mp)%BM


! calculate ExB drift (4.17) at IN/IS foot point
!CAUTION! these are geographic
!VEXB%east  = (v_e(1) * apexE(1,ift,mp)%east)  + (v_e(2) * apexE(2,ift,mp)%east)
!VEXB%north = (v_e(1) * apexE(1,ift,mp)%north) + (v_e(2) * apexE(2,ift,mp)%north)
!VEXB%up    = (v_e(1) * apexE(1,ift,mp)%up)    + (v_e(2) * apexE(2,ift,mp)%up)
!print *,'VEXB[m/s]',VEXB !%east,VEXB%north,VEXB%up

! assume Lambda_m in Richmond 95 equals to theta_t1
        cos2Lambda_m = coslam_m(ihem,lp) * coslam_m(ihem,lp) ! 0<cos2<1
        IF ( ihem==1 ) THEN
          sinLambda_m  = + SQRT( 1.0 - cos2Lambda_m )  !>0 ---NH 
        ELSE IF ( ihem==2 ) THEN 
          sinLambda_m  = - SQRT( 1.0 - cos2Lambda_m )  !<0 ---SH 
        END IF
        sinIm = 2.0 * sinLambda_m / SQRT(4.0-3.0*cos2Lambda_m)
        v_e2magsouth = v_e(2) * sinIm
if(sw_debug)&
& print *,'sub-St:',ihem,'v_e2*sinI[m/s]',v_e2magsouth,' sinI',sinIm

        theta_t0(ihem) = theta_t1(ihem) - ( v_e2magsouth * time_step) / r
        phi_t0(ihem)   = phi_t1
!make sure that:  0*dtr<=theta<=180*dtr=pi, otherwise goes over to the other side of the poles
        IF ( theta_t0(ihem)<0.0  ) THEN
          theta_t0(ihem) = theta_t0(ihem) * (-1.0)
          phi_t0(ihem)   = phi_t0(ihem) + pi  !other side
        ELSE IF ( theta_t0(ihem)>pi   ) THEN
          theta_t0(ihem) = pi*2.0 - theta_t0(ihem)
          phi_t0(ihem)   = phi_t0(ihem) + pi  !other side
        END IF
         
        phi_t0(ihem)   = phi_t0(ihem) - ( v_e(1) * time_step ) / ( r * sinLambda_m )
!make sure that:  0*dtr<=phi<360*dtr=2pi
        IF ( phi_t0(ihem)>=pi*2.0 ) THEN 
          phi_t0(ihem) = phi_t0(ihem) - pi*2.0
        ELSE IF ( phi_t0(ihem)< 0.0    ) THEN 
          phi_t0(ihem) = phi_t0(ihem) + pi*2.0
        END IF

      END DO which_hemisphere !: DO ihem=1,2

!if(sw_debug) then
if(lp==130) then
!sintheta=SIN(theta)
!sin2theta = sintheta*sintheta 
!r0_apex = r/sin2theta 
ihem=1 !only
if(sw_debug) print *,mlon_rad(mp)
midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
r_apex = earth_radius + plasma_grid_Z(midpoint,lp) ![m]

sintheta=SIN( theta_t1(ihem) )
sin2theta = sintheta*sintheta 
r_apex_dipole = r/sin2theta 
if(sw_debug) print "('sub-St:T1: phi=',F12.6,' theta=',F12.6,' r=',2E13.5)",phi_t1*rtd,(90.-theta_t1(ihem)*rtd)     ,r_apex_dipole,r_apex

sintheta=SIN( theta_t0(ihem) )
sin2theta = sintheta*sintheta 
r0_apex = r/sin2theta 
if(sw_debug) print "('sub-St:T0: phi=',F12.6,' theta=',F12.6,' r=',E13.5)",phi_t0(ihem)*rtd,(90.-theta_t0(ihem)*rtd),r0_apex
if(sw_debug) print "(4E12.4)", (theta_t1(ihem)-theta_t0(ihem))*rtd, (r_apex_dipole-r0_apex), (( v_e2magsouth * time_step) / r) , ( v_e2magsouth * time_step)
end if

      END SUBROUTINE stepback_mag_th

!nm20130201: separated into the new module subroutine!
!      SUBROUTINE stepback_mag_R (utime,mp,lp,phi_t0,theta_t0,r0_apex)
      END MODULE module_stepback_mag_th
