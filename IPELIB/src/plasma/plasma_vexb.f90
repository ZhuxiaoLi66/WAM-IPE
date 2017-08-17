!20111110: v18: included zonal transport
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

      SUBROUTINE plasma_vexb (utime,mp,lp,midpoint,r_apex,ihem )         
      USE module_precision
      USE module_input_parameters,ONLY: sw_exb_up,sw_debug,sw_perp_transport
      USE module_plasma,ONLY: VEXB3D
      USE module_eldyn,ONLY: Ed1_90,Ed2_90
      USE module_FIELD_LINE_GRID_MKS,ONLY: Be3,apexE,plasma_grid_3d,geographic_coords,l_mag,JMIN_IN,JMAX_IS
      USE module_physical_constants,ONLY: pi
      IMPLICIT NONE
! INPUT
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec), INTENT(IN) :: mp  !mag-lon index
      INTEGER (KIND=int_prec), INTENT(IN) :: lp  !mag-lat index
      INTEGER (KIND=int_prec), INTENT(IN) :: midpoint,ihem 
      REAL (KIND=real_prec),   INTENT(IN) :: r_apex
! OUTPUT

!---local
      REAL(KIND=real_prec) :: v_e(2)   !1:ed2/be3 (4.18) ;2: -ed1/be3 (4.19)
      TYPE(geographic_coords) :: VEXBgeo
      REAL(KIND=real_prec) :: GLON_deg, LT_SEC
      INTEGER (KIND=int_prec) :: i,IN,IS
!---
      IN=JMIN_IN(lp)
      IS=JMAX_IS(lp)

      IF ( sw_exb_up<=1 ) THEN 

! (0) self consistent electrodynamics
!...comming soon...

! (1) WACCM E empirical model
! Ed1/2[V/m] at ( phi_t1(mp), theta_t1(lp) ), Be3[T]
!note: Ed1_90, Ed2_90, Be3 are constant along magnetic field lines!!! 
        v_e(1) =   Ed2_90(1,lp,mp) / Be3(lp,mp) !(4.18) +mag-east(d1?) 
        v_e(2) = - Ed1_90(1,lp,mp) / Be3(lp,mp) !(4.19) +down/equatorward(d2?)
if(sw_debug)&
& print *,'sub-vexb:',ihem,'ve2[m/s]',v_e(2),'ed1[mV/m]', Ed1_90(1,lp,mp)*1.0E+3,' Be3[tesla]',Be3(lp,mp) 


        i_loop: DO i=IN,IS
!ExB in geographic frame 
           VEXBgeo%east  = (v_e(1) * apexE(1,i,mp)%east )  + (v_e(2) * apexE(2,i,mp)%east)
           VEXBgeo%north = (v_e(1) * apexE(1,i,mp)%north)  + (v_e(2) * apexE(2,i,mp)%north)
           VEXBgeo%up    = (v_e(1) * apexE(1,i,mp)%up   )  + (v_e(2) * apexE(2,i,mp)%up)

!ExB in magnetic APEX frame
!(1) magnetic eastward exact horizontal
           IF ( sw_perp_transport>=2 ) THEN
             VEXB3D(1,i,mp) = VEXBgeo%east * l_mag(1,i,mp)%east  +  VEXBgeo%north * l_mag(1,i,mp)%north  +  VEXBgeo%up * l_mag(1,i,mp)%up
           ELSE !           IF ( sw_perp_transport==1/0 )
             VEXB3D(1,i,mp) = 0.00 
           END IF
!(2) upward
           VEXB3D(2,i,mp) = VEXBgeo%east * l_mag(2,i,mp)%east  +  VEXBgeo%north * l_mag(2,i,mp)%north  +  VEXBgeo%up * l_mag(2,i,mp)%up
        END DO i_loop !: DO i=IN,IS

if(sw_debug)&
& print *,'sub-vexb:' &
&, v_e(1),apexE(1,midpoint,mp)%up,v_e(2),apexE(2,midpoint,mp)%up

        ELSE IF ( sw_exb_up==2 ) THEN 

!(2) GIP empirical model

        ELSE IF ( sw_exb_up==3 ) THEN 

!(3)  SUPIM empirical model: 
!note: becomes zero at R=4000km

          GLON_deg = plasma_grid_3d(midpoint,mp)%GLON*180./pi
          LT_SEC = utime + GLON_deg/15.*3600.
          IF ( LT_SEC>=86400.)  LT_SEC = MOD ( LT_SEC, 86400. )
          IF ( LT_SEC<     0.)  LT_SEC = LT_SEC + 86400.
          CALL supim_EXBV(utime,lp,LT_SEC,r_apex,GLON_deg ,VEXB3D(2,midpoint,mp))

        ELSE IF ( sw_exb_up==4 ) THEN 

!(4) both are zero for the purpose of debugging the tranport
          VEXB3D(1:2,IN:IS,mp) = 0.0

        END IF !ELSE IF ( sw_exb_up==3 ) THEN 

if(sw_debug)&
& print *,'sub-vexb:',ihem,lp,mp,'v_exb_apex[m/s]: east',VEXB3D(1,midpoint,mp),' upward',VEXB3D(2,midpoint,mp)          
      END SUBROUTINE plasma_vexb    

