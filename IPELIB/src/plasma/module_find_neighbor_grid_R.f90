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
      USE module_physical_constants,ONLY: rtd,earth_radius,pi
      USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,JMIN_IN,JMAX_IS,mlon_rad,dlonm90km,plasma_grid_Z,minTheta,maxTheta,midpnt
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,lpHaloSize,mpHaloSize,MaxLpHaloUsed,MaxMpHaloUsed,mype,parallelBuild
      use module_PLASMA,only:utime_save
     IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
      REAL(KIND=real_prec8),INTENT(IN) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
      REAL(KIND=real_prec8),INTENT(IN) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      REAL(KIND=real_prec8), INTENT(IN) :: r0_apex ![meter]
!---local
      INTEGER (KIND=int_prec) :: ihem,ihem_max
      REAL(KIND=real_prec) :: Z_t0
      INTEGER (KIND=int_prec),DIMENSION(2,2), INTENT(OUT) :: mp_t0,lp_t0 !1st rank ihem=2 is not used
      INTEGER (KIND=int_prec) :: lp_min,l
      INTEGER (KIND=int_prec) :: lp1,lp2,midpoint1,midpoint2,mpx,mpp,mpm,lpx,lpp,lpm
      INTEGER (KIND=int_prec),PARAMETER :: missing_value=-999       
!dbg20141212 debug parallel
      INTEGER (KIND=int_prec) :: mpt02,mpt01
      REAL(KIND=real_prec) :: mlon0
!---
!array initialization
      mp_t0 = missing_value
      lp_t0 = missing_value

!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
!dbg20120509 IF ( sw_perp_transport(mp)==3 ) THEN 
      IF ( sw_perp_transport==3 ) THEN 
         ihem_max=2
      ELSE
         ihem_max=1
      END IF

!NH only for debug purpose: assume for the moment that flux tubes in NH/SH are moving together....
!if NH flux tube is moving differently from SH flux, run the loop upto ihem=2, then flux tube interpolation should be done separately between NH vs.SH
      which_hemisphere: DO ihem=1,1  !ihem_max
         mpx_loop: DO mpx=0,NMP
            if(mpx+1 > mpHaloSize) then
!SMS$ignore begin
               print*,mype,'mpx+1 > mpHaloSize in find_neighbor_grid_R: mpx=',mpx,' mpHaloSize=',mpHaloSize,' mp=',mp,' lp=',lp, ' phi_t0=',phi_t0(ihem)*rtd, (90.-theta_t0(ihem)*rtd)
               print*,'Increase the halo size or take smaller time steps.'
               print*,'Stopping in find_neighbor_grid_R'
!SMS$ignore end
               STOP
            endif !(mpx+1 > mpHaloSize) then
            MaxMpHaloUsed = max(MaxMpHaloUsed,mpx+1)
            mpp=mp+mpx
            mpm=mp-mpx

!nm20160419 error trap(1)
            if ( mpp<(1-mpHaloSize) .or. (NMP+mpHaloSize)<(mpp+1) )then
!SMS$IGNORE begin
               print*,mype,utime_save,'!STOP! INVALID mpp=',mpp,mpx,lp,mp,lbound(mlon_rad),ubound(mlon_rad),phi_t0(ihem),mpHaloSize
!SMS$IGNORE end
               STOP
            endif

!(1)when plasma flows eastward
            IF( mlon_rad(mpp)<=phi_t0(ihem).AND.phi_t0(ihem)<mlon_rad(mpp+1) ) THEN
               mp_t0(ihem,1) = mpp
               mp_t0(ihem,2) = mpp+1 
               EXIT mpx_loop
            END IF !( mlon_rad(mpp)<=phi_t0

!nm20160419 error trap(2)
            if ( (mpm-1)<(1-mpHaloSize) .or. (NMP+mpHaloSize)<mpm )then
!SMS$IGNORE begin
               print*,mype,utime_save,'!STOP! INVALID mpm(1)=',mpm,mpx,lp,mp,lbound(mlon_rad),ubound(mlon_rad),phi_t0(ihem),mpHaloSize
!SMS$IGNORE end
               STOP
            endif

!(2)when plasma flows westward
            IF( mlon_rad(mpm-1)<=phi_t0(ihem).AND.phi_t0(ihem)<mlon_rad(mpm) ) THEN
               mp_t0(ihem,1) = mpm-1
               mp_t0(ihem,2) = mpm
               EXIT mpx_loop
            END IF !( mlon_rad(mpm-1)<=phi_t0(i
 
         END DO mpx_loop !: DO mpx=0,NMP


!find  lp0_t0:NH
IF (ihem==1) THEN

!check pole regions! not totally sure whether I should use theta_t0 or r0_apex???
IF ( theta_t0(ihem) < minTheta ) THEN 

   lp_t0(ihem,1)=missing_value !-999
   lp_t0(ihem,2)=1

!SMS$IGNORE begin
   if(sw_debug)print"('mype=',i3,'subFin:specialPole:mp=',i3,'lp=',i3)",mype,mp,lp
!SMS$IGNORE end
   RETURN

ELSE IF ( theta_t0(ihem) > maxTheta ) THEN
   print *,'sub-Fi_R: !STOP! invalid theta_t0',mp,lp,theta_t0(ihem),maxTheta
   STOP
ELSE   !IF ( plasma_grid_GL( JMIN_IN(lp),lp ) <= theta_t0(ihem) ) THEN 



!    lp_min =lp-5 !not sure if 10 is enough???
!    if (lp_min<=0 ) lp_min=1
!    IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN 
!      lp_min=lp-5
!      if (lp_min<=0 ) lp_min=1
!      IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN 
!        print *,'sub-Fi_R:NH !STOP! not sure if this is working???'
!        STOP
!      END IF
!    END IF

z_t0 = r0_apex - earth_radius



lpx_loop: DO lpx=0,NLP-1  !nearest point-->EQ
  IF(lpx+1 > lpHaloSize) THEN
!SMS$ignore begin
    print*,'Searching for inner,outer flux tube: lpx+1 > lpHaloSize',lpx,lpHaloSize,lp
    print*,'Increase the halo size or take smaller time steps.'
    print*,'Stopping in find_neighbor_grid_R'
!SMS$ignore end
    STOP
  ENDIF
  MaxLpHaloUsed = max(MaxLpHaloUsed,lpx+1)
  lpp=lp+lpx
  IF(lpp > NLP-1) lpp= lpp-NLP+1
  lpm=lp-lpx



  IF(lpm < 1) lpm= NLP-1+lpm

!nm20140627
!r_method
!t  IF(plasma_grid_Z(midpnt(lpp+1),lpp+1)<=Z_t0.AND.Z_t0<plasma_grid_Z(midpnt(lpp),lpp)) THEN
!th method
  IF(plasma_grid_GL(JMIN_IN(lpp),lpp)<=theta_t0(ihem).AND.theta_t0(ihem)<plasma_grid_GL(JMIN_IN(lpp+1),lpp+1)) THEN
    lp_t0(ihem,1)=lpp   !1=outer flux tube
    lp_t0(ihem,2)=lpp+1 !2=inner flux tube
    EXIT lpx_loop
  ENDIF
!r_method
!t  IF(plasma_grid_Z(midpnt(lpm),lpm)<=Z_t0.AND.Z_t0<plasma_grid_Z(midpnt(lpm-1),lpm-1)) THEN
!th method
  IF(plasma_grid_GL(JMIN_IN(lpm-1),lpm-1)<=theta_t0(ihem).AND.theta_t0(ihem)<plasma_grid_GL(JMIN_IN(lpm),lpm)) THEN
    lp_t0(ihem,1)=lpm-1 !1=outer flux tube
    lp_t0(ihem,2)=lpm   !2=inner flux tube
    EXIT lpx_loop
  ENDIF
  IF (lpx==NLP-1) THEN
!SMS$IGNORE begin
    print*,'Could not find inner,outer flux tube',lpp,lpm,midpnt(lpp),midpnt(lpp+1),midpnt(lpm),midpnt(lpm+1)
    print*,Z_t0,plasma_grid_Z(midpnt(lpp+1),lpp+1),plasma_grid_Z(midpnt(lpp),lpp),plasma_grid_Z(midpnt(lpm+1),lpm+1),plasma_grid_Z(midpnt(lpm),lpm)
    print*,'Stopping in find_neighbor_grid_R'
!SMS$IGNORE end
    STOP
  ENDIF
ENDDO lpx_loop !: DO lpx=0,NLP-1

!OUT  lp_t0(ihem,1)=l  
!IN   lp_t0(ihem,2)=l+1
lp1 = lp_t0(ihem,1)
midpoint1 = midpnt(lp1)
lp2 = lp_t0(ihem,2) !=l+1
midpoint2 = midpnt(lp2)



END IF! ( plasma_grid_3d(IN,lp)%GL <= theta_t0(ihem) ) THEN 
END IF !(ihem==1) THEN
END DO which_hemisphere!:  DO ihem=1,ihem_max


      END SUBROUTINE find_neighbor_grid_R
      END MODULE module_find_neighbor_grid_R
