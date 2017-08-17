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
      MODULE module_find_neighbor_grid_th
        PRIVATE
        PUBLIC :: find_neighbor_grid_th
      CONTAINS
      SUBROUTINE find_neighbor_grid_th ( mp,lp &
&, phi_t0 , theta_t0 &
&,  mp_t0 ,    lp_t0 )
      USE module_precision
      USE module_physical_constants,ONLY: rtd
      USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,JMIN_IN,JMAX_IS,mlon_rad,dlonm90km,minTheta
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,lpHaloSize,mpHaloSize,MaxLpHaloUsed,MaxMpHaloUsed,mype

     IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec ),INTENT(IN) :: mp
      INTEGER (KIND=int_prec ),INTENT(IN) :: lp
      REAL    (KIND=real_prec),INTENT(IN) :: phi_t0(2)  !magnetic longitude,phi[rad] at T0(previous time step)
      REAL    (KIND=real_prec),INTENT(IN) :: theta_t0(2)!magnetic latitude,theta[rad] at T0
!---local
      INTEGER (KIND=int_prec ) :: ihem,ihem_max
      INTEGER (KIND=int_prec ),DIMENSION(2,2), INTENT(OUT) :: mp_t0,lp_t0
      INTEGER (KIND=int_prec ) :: mpx,mpp,mpm,lpx,lpp,lpm
!---

!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
!dbg20120509: IF ( sw_perp_transport(mp)==3 ) THEN 
IF ( sw_perp_transport==3 ) THEN 
  ihem_max=2
ELSE
  ihem_max=1
END IF

!NH only for debug purpose
!assume for the moment that flux tubes in NH/SH are moving together....
!if NH flux tube is moving differently from SH flux, run the loop upto ihem=2, then flux tube interpolation should be done separately between NHv.s.SH
which_hemisphere: DO ihem=1,1  !ihem_max
!!!dbg20120125: why mp=19/20 are not working??? 
!!!dbg20120125:  mlon_deg = phi_t0(ihem)*rtd
!!!dbg20120125:  mp_t0(ihem,1) = INT( (mlon_deg/dlonm90km) , int_prec )+1
!!!dbg20120125:  mp_t0(ihem,2) = mp_t0(ihem,1)+1
  mpx_loop: DO mpx=0,NMP
    if(mpx+1 > mpHaloSize) then
      print*,'mpx+1 > mpHaloSize in find_neighbor_grid',mpx,mpHaloSize,mp
      print*,'Increase the halo size or take smaller time steps.'
      print*,'Stopping in find_neighbor_grid'
      STOP
    endif
    MaxMpHaloUsed = max(MaxMpHaloUsed,mpx+1)
    mpp=mp+mpx
    if(mpp > NMP) mpp= mpp-NMP
    mpm=mp-mpx
    if(mpm < 1) mpm= NMP+mpm
    IF ( mlon_rad(mpp)<=phi_t0(ihem) .AND. phi_t0(ihem)<mlon_rad(mpp+1) ) THEN
      mp_t0(ihem,1) =mpp
      mp_t0(ihem,2) =mpp+1
      EXIT mpx_loop
    END IF
    IF ( mpm>1.and.mlon_rad(mpm)<=phi_t0(ihem) .AND. phi_t0(ihem)<mlon_rad(mpm-1) ) THEN
      mp_t0(ihem,1) =mpm-1
      mp_t0(ihem,2) =mpm
      EXIT mpx_loop
    END IF
  END DO mpx_loop !: DO mpx=0,NMP
!dbg20120125:
if(sw_debug) print *,'dbg20120125! sub-find_neighbor_grid_th:', mp_t0(ihem,1:2),phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,1:2))*rtd
!STOP

!find  lp0_t0:NH
IF (ihem==1) THEN

!check pole regions!
  IF ( theta_t0(ihem) < minTheta ) THEN 
    lp_t0(ihem,1)=-999
    lp_t0(ihem,2)=1
    print *,'sub-Fi_th: mp',mp,' lp',lp,'needs special pole interpolation'
    RETURN
  END IF! ( plasma_grid_3d(IN,lp)%GL <= theta_t0(ihem) ) THEN 

  lpx_loop: DO lpx=0,NLP-1  !nearest point-->EQ
    IF(lpx+1 > lpHaloSize) THEN
      print*,'lpx+1 > lpHaloSize in find_neighbor_grid',lpx,lpHaloSize,lp
      print*,'Increase the halo size or take smaller time steps.'
      print*,'Stopping in find_neighbor_grid'
      STOP
    ENDIF
    MaxLpHaloUsed = max(MaxLpHaloUsed,lpx+1)
    lpp=lp+lpx
    IF(lpp > NLP-1) lpp= lpp-NLP+1
    lpm=lp-lpx
    IF(lpm < 1) lpm= NLP-1+lpm
    IF(plasma_grid_GL(JMIN_IN(lpp),lpp)<=theta_t0(ihem).AND.theta_t0(ihem)<plasma_grid_GL(JMIN_IN(lpp+1),lpp+1)) THEN
      lp_t0(ihem,1)=lpp
      lp_t0(ihem,2)=lpp+1
      EXIT lpx_loop
    ENDIF
    IF(lpm>1.and.plasma_grid_GL(JMIN_IN(lpm),lpm)<=theta_t0(ihem).AND.theta_t0(ihem)<plasma_grid_GL(JMIN_IN(lpm-1),lpm-1)) THEN
      lp_t0(ihem,1)=lpm-1
      lp_t0(ihem,2)=lpm
      EXIT lpx_loop
    ENDIF
    IF (lpx==NLP-1) THEN
      print*,'Could not find lp',lpp,lpm,plasma_grid_GL(JMIN_IN(lpp),lpp),plasma_grid_GL(JMIN_IN(lpm),lpm),theta_t0(ihem)
      print*,'Stopping in find_neighbor_grid'
      STOP
    ENDIF
  ENDDO lpx_loop !: DO lpx=0,NLP-1

END IF !(ihem==1) THEN

if(sw_debug) print *,'sub-Fi_th: mlon', mlon_rad(mp_t0(ihem,1))*rtd, phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,2))*rtd, mp_t0(ihem,1:2)
if(sw_debug) print *,'sub-Fi_th: mlat',plasma_grid_GL( JMIN_IN(lp_t0(ihem,1)),lp_t0(ihem,1) )*rtd, theta_t0(ihem)*rtd, plasma_grid_GL( JMIN_IN(lp_t0(ihem,2)),lp_t0(ihem,2) )*rtd

END DO which_hemisphere!:  DO ihem=1,ihem_max


if(sw_debug) print *,  'sub-find_neighbor_grid_th finished!'
      END SUBROUTINE find_neighbor_grid_th
      END MODULE module_find_neighbor_grid_th
