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
CHARACTER(LEN=2), PARAMETER :: title_lat='60'    !.. latitude of the flux tube

INTEGER(KIND=int_prec), PARAMETER :: IPDIM = &   !.. Field line grid dimension
       !&   265 !mlat=10[deg] at 90km
       !&   339 !mlat=35[deg] ---modified110210
       !&  1357 !mlat=35[deg] at 90km
       !&  2141 !mlat=60[deg] at 90km
       !&  4501 !mlat=85[deg] at 90km highres
                !mlat=88[deg] ---modified20101116 
       !&  6415 !mlat=88[deg] at 90km
       !&  1069 !mlat=88[deg] at 90km
       !&  1473 !lp=41: mlat=38.85[deg] at 90km
       !&  1435 !lp=42: mlat=37.40[deg] at 90km
       !&  1393 !lp=43: mlat=35.81[deg] at 90km
       !&  1357 !lp=44: mlat=34.39[deg]
       !&  1317 !lp=45: mlat=32.82[deg]
       !&   495 !lp=46: mlat=31.43[deg]
       !&  1117 !mlat=88[deg] at 90km lowres
       &  1115 !mlat=88[deg] at 90km lowres

INTEGER(KIND=int_prec) :: NMP !# of points along all the flux tubes in one longitude plane
INTEGER(KIND=int_prec) :: NLP !# of flux tubes in longitude the entire grid

INTEGER(KIND=int_prec) :: NPTS2D
!sw_grid=0:NPTS2D=44438; 1:44514; 2:31287
!INTEGER(KIND=int_prec), PARAMETER :: NMP_all =80
!INTEGER(KIND=int_prec), PARAMETER :: NMP0  =1       !mp lower bound 4 allocating array
!INTEGER(KIND=int_prec), PARAMETER :: NMP1  =NMP_all !mp upper bound
!INTEGER(KIND=int_prec), PARAMETER :: NLP_all =170 
  !low res # of flux tubes in a meridional(height-latitutde) plane (the entire grid)
!INTEGER(KIND=int_prec), PARAMETER :: NLP =NLP_all
  !low res # of flux tubes in a meridional(height-latitutde) plane for allocation
INTEGER(KIND=int_prec), PARAMETER :: ISPEC = 9   !.. Species dimension for density
INTEGER(KIND=int_prec), PARAMETER :: ISPEV = 4   !.. Species dimension for V field line velocity
INTEGER(KIND=int_prec), PARAMETER :: ISPET = 2   !.. Species dimension for T temperature
INTEGER(KIND=int_prec), PARAMETER :: ISTOT =16   !..=ISPEC+3+ISEV,TOTAL Species
