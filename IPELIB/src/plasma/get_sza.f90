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
!
! I should modify the code to use the method used in FLIP!!!
! calculate Solar Zenith Angle [radians]
      SUBROUTINE Get_SZA ( utime,mp,lp,SZA_rad )
        USE module_precision
        USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_grid_3d,JMIN_IN,JMAX_IS,ISL,IBM,IGR,IQ,IGCOLAT,IGLON
        USE module_physical_constants,ONLY: pi,zero
        USE module_input_parameters,ONLY: NDAY,sw_debug
        USE module_IPE_dimension,ONLY: IPDIM
        IMPLICIT NONE
        INTEGER (KIND=int_prec),  INTENT(IN)  :: utime    !universal time [sec]
        INTEGER (KIND=int_prec),  INTENT(IN)  :: mp,lp
        DOUBLE PRECISION,DIMENSION(IPDIM), INTENT(OUT) :: SZA_rad !solar zenith angle [radians]
!---local variables---
        INTEGER (KIND=int_prec) :: i, ii
        REAL (KIND=real_prec) :: rlat  !.. geographic latitude[radian]
        REAL (KIND=real_prec) :: ty
        REAL (KIND=real_prec) :: sda   !.. declination angle [radian]
        REAL (KIND=real_prec) :: utime12   !.. time counting from 12UT [sec]
        REAL (KIND=real_prec) :: ssa   !.. sub solar angle [deg]
        REAL (KIND=real_prec) :: rlt   !.. =ssa+180 [deg]
        REAL (KIND=real_prec) :: rlt_r !.. =ssa+180 [rad]
        INTEGER (KIND=int_prec) :: IN,IS
        REAL (KIND=real_prec) :: cos_sza !.. COS(SZA_rad)
!        INTEGER (KIND=int_prec) :: stat_alloc

!
        SZA_rad(:) = zero
        IN = JMIN_IN(lp)
        IS = JMAX_IS(lp)
!        IF (.NOT. ALLOCATED(SZA_rad) )  ALLOCATE ( SZA_rad(IN:IS) ,STAT=stat_alloc)
!if ( stat_alloc/=0 ) then
!  print *, ALLOCATED( sza_rad )
!  print *,"!STOP! ALLOCATION FAILD! in get_sza:",stat_alloc,mp,lp,in,is
!  STOP
!endif

        field_line_loop: DO i=IN,IS
          ii = i-IN+1 !make sure SZA_rad(1:~)
!!        theta = aa - (m-1.0)*bb/2.0  !geographic CO-latitude [deg]
          rlat = pi/2.0 - plasma_grid_3d(i,lp,mp,IGCOLAT)  ![radian]
          ty = (NDAY+15.5)*12.0/365.0
          IF ( ty>12.0 ) ty = ty - 12.0
!.. declination angle [radian]
!NOTE: it should be consistent with APEX coordinates:
! declination angle = ATAN2 ( BEAST, BNORTH ); where BEAST/NORTH: east/northward magnetic field components
          sda = ATAN(0.434*SIN(pi/6.0*(ty-3.17)))   

!!        phi = (FLOAT(l)-1.0)*18.0
! FLOAT(nn)*DTIME:  nn: # of min from 12UT; DTIME=60[sec]
          IF ( utime >= 12*3600 ) then
             utime12 = REAL(utime) -12.0*3600.0
          ELSE  !IF ( utime < 12*3600 ) then
             utime12 = REAL(utime) +12.0*3600.0
          ENDIF
! Sub Solar Angle: the angle at some point on the Earth to the east from a line running north/south through the sub-solar point...
          ssa = plasma_grid_3d(i,lp,mp,IGLON)*180.0/pi + utime12/240.0  
          rlt = 180.0 + ssa
          IF ( rlt>360.0 ) rlt = rlt - 360.0
          rlt_r = rlt*pi/180.0  !convert from [deg] to [rad]
          cos_sza = -COS(rlat)*COS(sda)*COS(rlt_r)+SIN(rlat)*SIN(sda)
!dbg20111001: need to check the value before taking ACOS 
          IF ( cos_sza<-1.0 ) THEN
print *,'(1) !INVALID cos_sza',cos_sza,ii,mp,lp,COS(rlat),COS(sda),COS(rlt_r),SIN(rlat),SIN(sda)
            cos_sza=-1.0
          ELSE IF ( cos_sza>1.0 ) THEN
print *,'(2) !INVALID cos_sza',cos_sza,ii,mp,lp,COS(rlat),COS(sda),COS(rlt_r),SIN(rlat),SIN(sda)
            cos_sza=+1.0
          END IF
          SZA_rad(ii) = ACOS( cos_sza )

!dbg20110815
IF ( sw_debug )  write(unit=1005,FMT=*) i,ii,'sza_rad',sza_rad(ii),rlat,sda,rlt_r,-COS(rlat),COS(sda),COS(rlt_r),SIN(rlat),SIN(sda)
       END DO field_line_loop !: DO i=IN,IS

IF ( sw_debug ) & 
      print "('Get_SZA [deg]    =',2F10.4)",SZA_rad(IN-IN+1)*180./pi,SZA_rad(IS-IN+1)*180./pi
 
      END SUBROUTINE Get_SZA
