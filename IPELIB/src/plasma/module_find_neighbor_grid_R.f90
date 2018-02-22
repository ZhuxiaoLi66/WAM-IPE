!nm20130201: separated from perpendicular_transport.f90:
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
MODULE module_find_neighbor_grid_R
  PRIVATE
  PUBLIC :: find_neighbor_grid_R
CONTAINS
!20111005
! using r0_apex as in GIP
  SUBROUTINE find_neighbor_grid_R ( mp,lp &
  &, phi_t0 , theta_t0 &
  &, r0_apex &
  &,  mp_t0 ,    lp_t0 )
    USE module_precision
    USE module_physical_constants,ONLY: rtd,earth_radius
    USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,JMIN_IN,JMAX_IS,mlon_rad,dlonm90km,plasma_grid_Z,minTheta,maxTheta,midpnt
    USE module_IPE_dimension,ONLY: NMP,NLP
    USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,lpHaloSize,mpHaloSize,MaxLpHaloUSEd,MaxMpHaloUSEd,mype,parallelBuild
    IMPLICIT NONE
!--- INPUT ---
    INTEGER (KIND=int_prec),INTENT(IN) :: mp
    INTEGER (KIND=int_prec),INTENT(IN) :: lp
    REAL(KIND=REAL_prec8),INTENT(IN) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
    REAL(KIND=REAL_prec8),INTENT(IN) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
    REAL(KIND=REAL_prec8), INTENT(IN) :: r0_apex ![meter]
!---local
    INTEGER (KIND=int_prec) :: ihem,ihem_max
    REAL(KIND=REAL_prec) :: Z_t0
    INTEGER (KIND=int_prec),DIMENSION(2,2), INTENT(OUT) :: mp_t0,lp_t0 !1st rank ihem=2 is not USEd
    INTEGER (KIND=int_prec) :: lp_min,l
    INTEGER (KIND=int_prec) :: lp1,lp2,midpoint1,midpoint2,mpx,mpp,mpm,lpx,lpp,lpm
    INTEGER (KIND=int_prec),PARAMETER :: missing_value=-9999
!---
!array initialization
    mp_t0 = missing_value
    lp_t0 = missing_value

!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with dIFferent ExB drIFt
!dbg20120509 IF ( sw_perp_transport(mp)==3 ) THEN
    IF ( sw_perp_transport==3 ) THEN
      ihem_max=2
    ELSE
      ihem_max=1
    ENDIF

!NH only for debug purpose
!assume for the moment that flux tubes in NH/SH are moving together....
!IF NH flux tube is moving dIFferently from SH flux, run the loop upto ihem=2, THEN flux tube interpolation should be DOne separately between NHv.s.SH
    which_hemisphere: DO ihem=1,1  !ihem_max
!!!dbg20120125: why mp=19/20 are not working???
!!!dbg20120125:  mlon_deg = phi_t0(ihem)*rtd
!!!dbg20120125:  mp_t0(ihem,1) = INT( (mlon_deg/dlonm90km) , int_prec )+1
!!!dbg20120125:  mp_t0(ihem,2) = mp_t0(ihem,1)+1
      mpx_loop: DO mpx=0,NMP
        IF(mpx+1 > mpHaloSize) THEN
!SMS$ignore begin
          PRINT*,'mpx+1 > mpHaloSize in find_neighbor_grid_R: mpx=',mpx,' mpHaloSize=',mpHaloSize,' mp=',mp,' lp=',lp
          PRINT*,'Increase the halo size or take smaller time steps.'
          PRINT*,'STOPping in find_neighbor_grid_R:mype=',mype
!SMS$ignore end
          STOP
        ENDIF
        MaxMpHaloUSEd = max(MaxMpHaloUSEd,mpx+1)
        mpp=mp+mpx
        mpm=mp-mpx
        IF(.not.parallelBuild) THEN
          IF(mpp > NMP) mpp= mpp-NMP
          IF(mpm <   1) mpm= NMP+mpm
        ENDIF
        IF(mlon_rad(mpp)<=phi_t0(ihem).AND.phi_t0(ihem)<mlon_rad(mpp+1)) THEN
          mp_t0(ihem,1) =mpp
          mp_t0(ihem,2) =mpp+1

!dbg20140205 debug zonal transport
          IF( mpp==nmp ) THEN
            PRINT *,'!dbg20140205(1) mp=',mp, mp_t0(ihem,1), mp_t0(ihem,2), phi_t0(ihem)*rtd,mpp,ihem
          ENDIF
          EXIT mpx_loop
        ENDIF
        IF(mpm>1.and.mlon_rad(mpm)<=phi_t0(ihem).AND.phi_t0(ihem)<mlon_rad(mpm-1)) THEN
          mp_t0(ihem,1) =mpm-1
          mp_t0(ihem,2) =mpm
          EXIT mpx_loop
        ENDIF
      ENDDO mpx_loop !: DO mpx=0,NMP

      IF(sw_debug) THEN
        IF (mp==1.OR.mp==nmp) PRINT *,'dbg20140205! sub-find_neighbor_grid_R:mp=',mp, mp_t0(ihem,1:2),phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,1:2))*rtd
      ENDIF
!dbg20140205: correction IF mp>nmp
      IF (mp==nmp.and.mp_t0(ihem,2)>nmp ) THEN
        PRINT *,'!dbg20140205! mp_t0',mp_t0(ihem,2),mp
        mp_t0(ihem,2)=mp_t0(ihem,2)-nmp
        PRINT *,'!dbg20140205! corrected mp_t0',mp_t0(ihem,2)
      ENDIF


!find  lp0_t0:NH
      IF (ihem==1) THEN

!check pole regions!
!not totally sure whether I should USE theta_t0 or r0_apex???
        IF ( theta_t0(ihem) < minTheta ) THEN
          lp_t0(ihem,1)=missing_value !-999
          lp_t0(ihem,2)=1
          PRINT *,'sub-Fi_R: mp',mp,' lp',lp,'needs special pole interpolation'
          RETURN
        ELSE IF ( theta_t0(ihem) > maxTheta ) THEN
          PRINT *,'sub-Fi_R: !STOP! invalid theta_t0',mp,lp,theta_t0(ihem),maxTheta
          STOP
        ELSE   !IF ( plasma_grid_GL( JMIN_IN(lp),lp ) <= theta_t0(ihem) ) THEN

!!!UNDERCONSTRUCTION!!!
          IF(sw_debug) PRINT *,'sub-Fi_R: check GL NH[deg]',(90.-plasma_grid_GL( JMIN_IN(lp),lp )*rtd)

!    lp_min =lp-5 !not sure IF 10 is enough???
!    IF (lp_min<=0 ) lp_min=1
!    IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN
!      lp_min=lp-5
!      IF (lp_min<=0 ) lp_min=1
!      IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN
!        PRINT *,'sub-Fi_R:NH !STOP! not sure IF this is working???'
!        STOP
!      ENDIF
!    ENDIF

          z_t0 = r0_apex - earth_radius

!d l=130
!d PRINT *,JMIN_IN(l),JMAX_IS(l), midpnt(l),z_t0

          lpx_loop: DO lpx=0,NLP-1  !nearest point-->EQ
            IF(lpx+1 > lpHaloSize) THEN
!SMS$ignore begin
              PRINT'("Searching for inner,outer flux tube: lpx+1>lpHaloSize: lpx=",i4,"lpHaloSize=",i2,"lp=",i4)',lpx,lpHaloSize,lp
              PRINT*,'Increase the halo size or take smaller time steps.'
              PRINT*,'STOPping in find_neighbor_grid_R: mype=',mype
!SMS$ignore end
              STOP
            ENDIF
            MaxLpHaloUSEd = max(MaxLpHaloUSEd,lpx+1)
            lpp=lp+lpx
            IF(lpp > NLP-1) lpp= lpp-NLP+1
            lpm=lp-lpx
            IF(lpm < 1) lpm= NLP-1+lpm
            IF(plasma_grid_Z(midpnt(lpp+1),lpp+1)<=Z_t0.AND.Z_t0<plasma_grid_Z(midpnt(lpp),lpp)) THEN
              lp_t0(ihem,1)=lpp   !1=outer flux tube
              lp_t0(ihem,2)=lpp+1 !2=inner flux tube
              EXIT lpx_loop
            ENDIF
            IF(plasma_grid_Z(midpnt(lpm),lpm)<=Z_t0.AND.Z_t0<plasma_grid_Z(midpnt(lpm-1),lpm-1)) THEN
              lp_t0(ihem,1)=lpm-1 !1=outer flux tube
              lp_t0(ihem,2)=lpm   !2=inner flux tube
              EXIT lpx_loop
            ENDIF
            IF (lpx==NLP-1) THEN
              PRINT*,'Could not find inner,outer flux tube',lpp,lpm,midpnt(lpp),midpnt(lpp+1),midpnt(lpm),midpnt(lpm+1)
              PRINT*,Z_t0,plasma_grid_Z(midpnt(lpp+1),lpp+1),plasma_grid_Z(midpnt(lpp),lpp),plasma_grid_Z(midpnt(lpm+1),lpm+1),plasma_grid_Z(midpnt(lpm),lpm)
              PRINT*,'STOPping in find_neighbor_grid_R'
              STOP
            ENDIF
          ENDDO lpx_loop !: DO lpx=0,NLP-1

!OUT  lp_t0(ihem,1)=l
!IN   lp_t0(ihem,2)=l+1
          lp1 = lp_t0(ihem,1)
          midpoint1 = midpnt(lp1)
          lp2 = lp_t0(ihem,2) !=l+1
          midpoint2 = midpnt(lp2)
          IF(sw_debug) PRINT *,'lp1=',lp1,midpoint1,' lp2=',lp2,midpoint2
          IF(sw_debug) PRINT *,'sub-Fi_R: check R',(earth_radius+plasma_grid_Z(midpoint1,lp1)),(earth_radius+Z_t0),(earth_radius+plasma_grid_Z(midpoint2,lp2))

          IF(sw_debug) PRINT *,'sub-Fi_R: mp1=',mp_t0(ihem,1),' mp2=',mp_t0(ihem,2)
          IF(sw_debug) PRINT *,'sub-Fi_R: mlon',mlon_rad(mp_t0(ihem,1))*rtd, phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,2))*rtd, mp_t0(ihem,1:2)
          IF(sw_debug) PRINT *,'sub-Fi_R: mlat',(90.-plasma_grid_GL( JMIN_IN(lp_t0(ihem,1)),lp_t0(ihem,1) )*rtd), (90.-theta_t0(ihem)*rtd),(90.- plasma_grid_GL( JMIN_IN(lp_t0(ihem,2)),lp_t0(ihem,2) )*rtd)


        ENDIF! ( plasma_grid_3d(IN,lp)%GL <= theta_t0(ihem) ) THEN
      ENDIF !(ihem==1) THEN
    ENDDO which_hemisphere!:  DO ihem=1,ihem_max

    IF(sw_debug) PRINT *,  'sub-find_neighbor_grid R finished!'
  END SUBROUTINE find_neighbor_grid_R
END MODULE module_find_neighbor_grid_R
