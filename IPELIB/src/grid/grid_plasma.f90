!20120112: copied from vapor: 20110909: source/main/module_field_line_grid.3d.glb.f90 to read new Q and TD file, called from driver_grid.f90
!20110414: modified to run the low resolution version
!20110120: 
!          3D version of module_field_line_grid.f90
MODULE module_FIELD_LINE_GRID_MKS
      USE module_precision
!nm20120112:       USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
      IMPLICIT NONE

!nm20120112:
INTEGER(KIND=int_prec), PARAMETER :: sw_newQ=1

!---IPE_dimension---got rid of the dependance on the module IPE_dimension...
INTEGER(KIND=int_prec), PARAMETER :: NPTS2D = 44438
INTEGER(KIND=int_prec), PARAMETER :: NMP    =    80
INTEGER(KIND=int_prec), PARAMETER :: NLP    =   170 
!nm20120112:for module_IO
      CHARACTER (LEN=150) :: filename
      INTEGER (KIND=int_prec), PARAMETER :: LUN_pgrid=7
! --- PRIVATE ---
!
! --- PUBLIC ---
!... read in parameters
      INTEGER(KIND=int_prec), DIMENSION(NMP,NLP),TARGET,PUBLIC :: JMIN_IN,JMAX_IS  !.. first and last indices on field line grid
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC ::  GCOLAT    !.. geographic co-latitude [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC ::  GLON      !.. geographic longitude [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC ::  Qvalue
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC ::  GL_rad    !.. magnetic co-latitude Eq(6.1) [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC ::  SL_meter  !.. distance of point from northern hemisphere foot point [meter]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC ::  BM_nT     !.. magnetic field strength [nT]
!nm20110822:no longer used
!      REAL(KIND=real_prec),               ALLOCATABLE, PUBLIC ::  SZA_rad(:) !solar zenith angle [radians]
! components (east, north, up) of base vectors
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC ::  D1    !.. Eq(3.8) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC ::  D2    !.. Eq(3.9) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC ::  D3    !.. Eq(3.10) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC ::  E1    !.. Eq(3.11) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC ::  E2    !.. Eq(3.12) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(2,NMP,NLP), PUBLIC ::  Be3         ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [nT]

!
!------------
!...calculated parameters
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP),     PUBLIC ::  Z_meter !.. altitude [meter]
      REAL(KIND=real_prec), DIMENSION(       NMP,NLP), PUBLIC ::  Pvalue  !.. p coordinate (L-shell)
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP),     PUBLIC ::  GR_mks  !.. Gravity [m2 s-1]


      PRIVATE :: read_plasma_grid_global,test_grid_output
      PUBLIC :: init_plasma_grid


      CONTAINS
!---------------------------
! initialise plasma grids
        SUBROUTINE init_plasma_grid ( )
        USE module_physical_constants,ONLY: earth_radius, pi, G0
        IMPLICIT NONE

!---for module_input_parameters...
      LOGICAL, PARAMETER :: sw_debug=.TRUE.

        INTEGER (KIND=int_prec) :: i, mp,lp
        REAL (KIND=real_prec) :: sinI
        INTEGER (KIND=int_prec), parameter :: sw_sinI=0  !0:flip; 1:APEX
        INTEGER (KIND=int_prec) :: in,is
        REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  r_meter_all     !.. distance from the center of the Earth[meter]
!---

!if the new GLOBAL 3D version
        CALL read_plasma_grid_global ( r_meter_all )


        apex_longitude_loop: DO mp = 1,NMP

IF ( sw_debug )  print *,'mp=',mp        

! make sure to use the MKS units.
      Z_meter(1:NPTS2D,mp) = r_meter_all(1:NPTS2D,mp) - earth_radius ![meter]
print *,"Z_meter calculation completed"

!.. p coordinate (L-shell) is a single value along a flux tube 
!NOTE: in FLIP, PCO is only used in setting up the rough plasmasphere H+ initial profiles (See PROFIN). It does not have to be accurate.
          apex_latitude_height_loop:   DO lp = 1,NLP

            IN = JMIN_IN(mp,lp)
            IS = JMAX_IS(mp,lp)

!nm20120112:            CALL Get_Pvalue_Dipole ( r_meter_all(IN,mp), GL_rad(IN,mp), Pvalue(mp,lp) )

!debug write
IF ( sw_debug ) THEN
print "('lp=',i6,'  IN=',i6,'  IS=',i6,'  NPTS=',i6)", lp,IN,IS,(IS-IN+1)
print "('r [m]      =',2E12.4)", r_meter_all(in,mp),r_meter_all(is,mp)
print "('G-LAT [deg]=',2f10.4)",(90.-GCOLAT(in,mp)*180./pi),(90.-GCOLAT(is,mp)*180./pi)
print "('M-LAT [deg]=',2f10.4)",(90.-GL_rad(in,mp)*180./pi),(90.-GL_rad(is,mp)*180./pi)
print "('GLON  [deg]=',2f10.4)",(GLON(in,mp)*180./pi),(GLON(is,mp)*180./pi)
print "('Qvalue     =',2E12.4)", Qvalue(in,mp), Qvalue(is,mp)
print "('BM [nT]    =',2E12.4)", BM_nT(in,mp), BM_nT(is,mp)
print "('D1         =',2E12.4)", D1(1,in,mp), D1(1,is,mp)
print "('D2         =',2E12.4)", D2(1,in,mp), D2(1,is,mp)
print "('D3         =',6E12.4)", D3(1:3,in,mp), D3(1:3,is,mp)
print "('E1         =',2E12.4)", E1(1,in,mp), E1(1,is,mp)
print "('E2         =',2E12.4)", E2(1,in,mp), E2(1,is,mp)
print "('Be3 [nT] NH/SH  =',2E12.4)", Be3(1,mp,lp), Be3(2,mp,lp)

print "('SL [m]     =',4E13.5)", SL_meter(in:in+1,mp), SL_meter(is-1:is,mp)
print "('Z  [m]     =',4E13.5)",  Z_meter(in:in+1,mp),  Z_meter(is-1:is,mp)
!nm20120112:print "('Pvalue     =',F10.4)", Pvalue(mp,lp)
END IF !( sw_debug ) THEN


! assuming Newtonian gravity: G0 is gravity at the sea level (z=0) 
!NOTE: positive in NORTHern hemisphere; negative in SOUTHern hemisphere
         flux_tube: DO i=IN,IS
!nm20120112:           CALL Get_sinI ( sw_sinI, sinI, GL_rad(i,mp), D3(1:3,i,mp) ) 
!nm20120112:           GR_mks(i,mp)  =  G0 * ( earth_radius * earth_radius ) / ( r_meter_all(i,mp) * r_meter_all(i,mp) ) * sinI * (-1.0)

!IF ( sw_debug )  print "(4E12.4)", Z_meter(i,mp),sinI, (G0 * ( earth_radius * earth_radius ) / ( r_meter_all(i,mp) * r_meter_all(i,mp) )),  GR_mks(i,mp)
        
         END DO flux_tube

!nm20120112:
!IF ( sw_debug )  then
!  if ( sw_sinI==0 ) then
!    print *,'sinI: flip'
!  else if ( sw_sinI==1 ) then
!    print *, 'sinI: APEX'
!  endif
!  print "('GRavity[m2 s-1]=',4E12.4)",GR_mks(in:in+2,mp),GR_mks(is,mp)
!END IF

       END DO apex_latitude_height_loop   !: DO lp = 1,NLP
     END DO apex_longitude_loop         !: DO mp = 1,NMP 

        END SUBROUTINE init_plasma_grid
!---------------------------
!20110726: the new GLOBAL 3D version: NMP=80
! global grid with the low resolution version
! magnetic longitude used for the grid. from 0 to 355.5 with 4.5 degree interval

      SUBROUTINE read_plasma_grid_global ( r_meter_all )
!nm20120112:        USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
!nm20120112:        USE module_physical_constants,ONLY: earth_radius,pi
!nm20120112:        USE module_IO,ONLY: filename,LUN_pgrid
        IMPLICIT NONE


!
!        integer (kind=int_prec), parameter :: NMP=80
!        integer (kind=int_prec), parameter :: NLP=170
!        integer (kind=int_prec), parameter :: NPTS2D=44438 
!        real(kind=8) gr_2d(npts2,nmp)
!        real(kind=8) gcol_2d(npts2,nmp)
!    real(kind=8) glon_2d(npts2,nmp)
!    real(kind=8) q_coordinate_2d(npts2,nmp)
!    real(kind=8) bcol_2d(npts2,nmp)
!    real(kind=8) Apex_D1_2d(3,npts2,nmp)
!    real(kind=8) Apex_D2_2d(3,npts2,nmp)
!    real(kind=8) Apex_D3_2d(3,npts2,nmp)
!    real(kind=8) Apex_E1_2d(3,npts2,nmp)
!    real(kind=8) Apex_E2_2d(3,npts2,nmp)
!    real(kind=8) Apex_grdlbm2_2d(3,npts2,nmp)
!    real(kind=8) integral_ds_2d(npts2,nmp)
!    real(kind=8) apex_BMAG_2d(npts2,nmp)
!    real(kind=8)  Apex_BE3_N(nmp,nlp), Apex_BE3_S(nmp,nlp)

!-------------
!... read in parameters
      INTEGER(KIND=int_prec), DIMENSION(NMP,NLP) :: JMIN_IN_all,JMAX_IS_all  !.. first and last indices on field line grid
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  GCOLAT_all  !.. geographic co-latitude [rad]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  GLON_all    !.. geographic longitude [rad]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  Qvalue_all
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  GL_rad_all      !.. magnetic co-latitude Eq(6.1) [rad]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  SL_meter_all  !.. distance of point from northern hemisphere foot point [meter]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP) ::  BM_nT_all      !.. magnetic field strength [nT]
! components (east, north, up) of base vectors
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP) ::  D1_all    !.. Eq(3.8) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP) ::  D2_all    !.. Eq(3.9) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP) ::  D3_all    !.. Eq(3.10) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP) ::  E1_all    !.. Eq(3.11) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP) ::  E2_all    !.. Eq(3.12) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(2,NMP,NLP) ::  Be3_all         ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [nT]

!-------------local
        CHARACTER (LEN=11) :: FORM_dum
        CHARACTER (LEN=7)  :: STATUS_dum
        REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(OUT) ::  r_meter_all     !.. distance from the center of the Earth[meter]
        CHARACTER(LEN=*), PARAMETER :: filepath_pgrid= &
!     & '/mnt/lfs0/projects/idea/maruyama/sandbox/ipe/field_line_grid/20110419lowres_global/'!low res global grid: APEX
     & '/mnt/lfs0/projects/idea/maruyama/sandbox/ipe/field_line_grid/20111025dipole/'       !low res global grid: tilted dipole
        CHARACTER(LEN=*), PARAMETER :: filename_pgrid= &
!      & 'GIP_apex_coords_global_lowres'        !APEX: original
!     & 'GIP_apex_coords_global_lowresCORRECTED20110824'        !APEX: 20110824: corrected for lp=1-6
     & 'GIP_tilt_coords_global_lowres'        !tilted dipole: original


      INTEGER(KIND=int_prec) :: sw_test_grid=1  !1: ON testgrid; 0: OFF 
!---

      filename =filepath_pgrid//filename_pgrid
      FORM_dum ='formatted' 
      STATUS_dum ='old'
      CALL open_file ( filename, LUN_pgrid, FORM_dum, STATUS_dum ) 
 
      print *,"open file completed"
      READ (UNIT=LUN_pgrid, FMT=*) JMIN_IN_all, JMAX_IS_all  !IN_2d_3d , IS_2d_3d
      print *,"reading JMIN_IN etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) r_meter_all,GCOLAT_all, GLON_all, Qvalue_all !gr_2d, gcol_2d, glon_2d, q_coordinate_2d
      print *,"reading r_meter etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) GL_rad_all          !bcol_2d
      print *,"reading GL_rad etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) SL_meter_all, BM_nT_all !integral_ds_2d, apex_BMAG_2d
      print *,"reading SL_meter etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) D1_all, D2_all, D3_all      !Apex_D1_2d
      print *,"reading D1-3 etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) E1_all, E2_all          !Apex_E1_2d
      print *,"reading E1/2 etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) Be3_all(1,1:NMP,1:NLP),Be3_all(2,1:NMP,1:NLP) !Apex_BE3_N
      print *,"reading Be3 etc completed"
      CLOSE (UNIT=LUN_pgrid)
      print *,"global grid reading finished, file closed..."

!nm20120112: new Q
if ( sw_newQ==1 ) then
      filename ='/mnt/lfs0/projects/idea/maruyama/sandbox/ipe/field_line_grid/20111005newq/new_q'
      FORM_dum ='formatted' 
      STATUS_dum ='old'
      CALL open_file ( filename, LUN_pgrid, FORM_dum, STATUS_dum )  
      print *,"new Q open file completed"
      READ (UNIT=LUN_pgrid, FMT=*) Qvalue_all  !q_coordinate_2d
      print *,"reading new Qvalue completed"
      CLOSE (UNIT=LUN_pgrid)
      print *,"new Q grid reading finished, file closed..."
end if !( sw_newQ==1 ) then


!dbg20110811:
IF ( sw_test_grid==1 ) THEN

  CALL test_grid_output ( JMIN_IN_all,JMAX_IS_all,r_meter_all,SL_meter_all,BM_nT_all,GL_rad_all &
& ,GCOLAT_all,GLON_all,Qvalue_all,D1_all,D2_all,D3_all,E1_all,E2_all,Be3_all &
&,filepath_pgrid,filename_pgrid)
END IF 

!assign only necessary parts
JMIN_IN(1:NMP,1:NLP)=JMIN_IN_all(    1:NMP,1:NLP)
JMAX_IS(1:NMP,1:NLP)=JMAX_IS_all(    1:NMP,1:NLP)
Be3(1:2,1:NMP,1:NLP)=    Be3_all(1:2,1:NMP,1:NLP)

GCOLAT(  1:NPTS2D,1:NMP) =  GCOLAT_all(1:NPTS2D,1:NMP)
GLON(    1:NPTS2D,1:NMP) =    GLON_all(1:NPTS2D,1:NMP)
Qvalue(  1:NPTS2D,1:NMP) =  Qvalue_all(1:NPTS2D,1:NMP)
GL_rad(  1:NPTS2D,1:NMP) =  GL_rad_all(1:NPTS2D,1:NMP)
SL_meter(1:NPTS2D,1:NMP) =SL_meter_all(1:NPTS2D,1:NMP)
BM_nT (  1:NPTS2D,1:NMP) =   BM_nT_all(1:NPTS2D,1:NMP)

  D1(1:3,1:NPTS2D,1:NMP) =  D1_all(1:3,1:NPTS2D,1:NMP)
  D2(1:3,1:NPTS2D,1:NMP) =  D2_all(1:3,1:NPTS2D,1:NMP)
  D3(1:3,1:NPTS2D,1:NMP) =  D3_all(1:3,1:NPTS2D,1:NMP)
  E1(1:3,1:NPTS2D,1:NMP) =  E1_all(1:3,1:NPTS2D,1:NMP)
  E2(1:3,1:NPTS2D,1:NMP) =  E2_all(1:3,1:NPTS2D,1:NMP)

      END SUBROUTINE  read_plasma_grid_global
!---
!---
SUBROUTINE test_grid_output ( JMIN_IN_all,JMAX_IS_all,r_meter_all,SL_meter_all,BM_nT_all,GL_rad_all &
& ,GCOLAT_all,GLON_all,Qvalue_all,D1_all,D2_all,D3_all,E1_all,E2_all,Be3_all &
& ,filepath_pgrid,filename_pgrid)
      USE module_physical_constants,ONLY: earth_radius,pi
!nm20120112:      USE module_IO,ONLY: filename,LUN_pgrid
      IMPLICIT NONE



      INTEGER(KIND=int_prec), DIMENSION(NMP,NLP),INTENT(INOUT) :: JMIN_IN_all,JMAX_IS_all  !.. first and last indices on field line grid
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  r_meter_all     !.. distance from the center of the Earth[meter]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  SL_meter_all  !.. distance of point from northern hemisphere foot point [meter]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  BM_nT_all      !.. magnetic field strength [nT]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  GL_rad_all      !.. magnetic co-latitude Eq(6.1) [rad]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  GCOLAT_all  !.. geographic co-latitude [rad]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  GLON_all    !.. geographic longitude [rad]
      REAL(KIND=real_prec8), DIMENSION(NPTS2D,NMP), INTENT(INOUT) ::  Qvalue_all
!
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP), INTENT(INOUT) ::  D1_all    !.. Eq(3.8) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP), INTENT(INOUT) ::  D2_all    !.. Eq(3.9) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP), INTENT(INOUT) ::  D3_all    !.. Eq(3.10) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP), INTENT(INOUT) ::  E1_all    !.. Eq(3.11) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(3,NPTS2D,NMP), INTENT(INOUT) ::  E2_all    !.. Eq(3.12) Richmond 1995
      REAL(KIND=real_prec8), DIMENSION(2,NMP,NLP), INTENT(INOUT) ::  Be3_all 




        CHARACTER(LEN=*), INTENT(IN) :: filepath_pgrid
        CHARACTER(LEN=*), INTENT(IN) :: filename_pgrid
! local
      INTEGER(KIND=int_prec) :: CTIPDIM, ctipdim1, midpoint,JMINX,JMAXX,ii,j, k,IN0(12),IS0(12),jj(12),IPDIM5,i
      INTEGER(KIND=int_prec),PARAMETER :: lp_start=1
      INTEGER(KIND=int_prec),PARAMETER :: lp_stop =8
      INTEGER(KIND=int_prec) :: mp,lp,dgp,lp0,counter,i7,i0,i1,jth
      REAL(KIND=real_prec),PARAMETER ::invalid_value=-9.99999999E+10
        CHARACTER (LEN=11) :: FORM_dum
        CHARACTER (LEN=7)  :: STATUS_dum

        REAL(KIND=real_prec), DIMENSION(3,  NPTS2D) :: tmp
        REAL(KIND=real_prec), DIMENSION(5,3,NPTS2D) :: tmp1
      INTEGER(KIND=int_prec),PARAMETER :: sw_test_grid_correct=1 !0:OFF; 1:ON
!

print *,"sub-test_grid_output: check filepath=",filepath_pgrid
print *,"filename=",filename_pgrid

!dbg20110824: fix SL_meter
print *,'!dbg! check SL_meter'
dgp=10
!do mp=1,1
mp=1
do lp=lp_start,lp_stop
      IN0(lp)=JMIN_IN_all(mp,lp)
      IS0(lp)=JMAX_IS_all(mp,lp)
      print *,"mp=",mp, " lp=",lp,IN0(lp),IS0(lp)

print *,"NH",SL_meter_all(IN0(lp):IN0(lp)+dgp,mp)
print *,"SH",SL_meter_all(IS0(lp)-dgp:IS0(lp),mp)
end do !lp 




!end do !mp


!JMINX   = IN1 - IN1 + 1
!JMAXX   = IS1 - IN1 + 1
!CTIPDIM1 = IS1 - IN1 + 1
!
!!test1: SL
! midpoint = IN1 + (JMAXX/2)+1 -1
!!NH
!DO i=IN1,midpoint
! SL_meter(i-IN1+1) = SL_meter3d(i,mp)
!END DO
!!SH
!DO i=IS1,midpoint+1, -1
! SL_meter(i-IN1+1) = SL_meter3d(i,mp)
!END DO

write(UNIT=1000,FMT="('Z_km',10i13)") (lp, lp=lp_start,lp_stop)
write(UNIT=1001,FMT="('SL_meter',10i13)")  (lp, lp=lp_start,lp_stop)
write(UNIT=1002,FMT="('Bm_nT',10i13)")  (lp, lp=lp_start,lp_stop)
write(UNIT=1003,FMT="('GL_rad',10i13)") (lp, lp=lp_start,lp_stop) 
write(UNIT=1004,FMT="('JMIN JMAX',10i13)")  (lp, lp=lp_start,lp_stop)
write(UNIT=1005,FMT="('GCOLAT',10i13)") (lp, lp=lp_start,lp_stop) 

IPDIM5 = JMAX_IS_all(mp,lp_start) - JMIN_IN_all(mp,lp_start) +1  !
DO i=1,IPDIM5  !IN0(5),IS0(5)

  ii = i-1
  do  lp=lp_start,lp_stop
    jj(lp) = ii+JMIN_IN_all(mp,lp)
    if ( jj(lp)>IS0(lp) )  jj(lp)=IS0(lp)
  end do

write(UNIT=1000,FMT="(i8,8E13.4)") i, ((r_meter_all(jj(lp),mp)-earth_radius)*1.0E-3, lp=lp_start,lp_stop)
write(UNIT=1001,FMT="(i8,8E13.4)") i, (SL_meter_all(jj(lp),mp), lp=lp_start,lp_stop) 
write(UNIT=1002,FMT="(i8,8E13.4)") i, ( BM_nT_all(jj(lp),mp), lp=lp_start,lp_stop)

write(UNIT=1003,FMT="(i8,8F10.4)") i, ( (90.-GL_rad_all(jj(lp),mp)*180./pi), lp=lp_start,lp_stop)
write(UNIT=1005,FMT="(i8,8F10.4)") i, ( (90.-GCOLAT_all(jj(lp),mp)*180./pi), lp=lp_start,lp_stop)
END DO

write(UNIT=1004,FMT="(4i8)") (lp,JMIN_IN_all(mp,lp),JMAX_IS_all(mp,lp), (JMAX_IS_all(mp,lp)-JMIN_IN_all(mp,lp)+1), lp=lp_start,lp_stop)
!dbg

if ( sw_test_grid_correct==1 ) then 
! CORRECT 
!(1) Infinity... in SL_meter
!(2) value(lp=1,6)<---value(lp=7) for r_meter,SL,BM
!(3) ctipdim(lp=1,6)=1117 <--- ctipdim(lp=7)=1115 

lp0=7
IPDIM5 = IS0(lp0) - IN0(lp0) +1 
do mp=1,NMP
do lp=1,6

 CTIPDIM  = IS0(lp) - JMIN_IN_all(mp,lp) + 1
 midpoint = JMIN_IN_all(mp,lp) + (CTIPDIM/2)+1 -1
print *,mp,lp,CTIPDIM,midpoint 


counter=0
 do i=JMIN_IN_all(mp,lp),IS0(lp)
   counter = counter + 1

   if ( counter <= IPDIM5 ) then
     i7 = IN0(lp0)+ (i-JMIN_IN_all(mp,lp))
if ( mp==1 .and. lp==1 ) print *,'!CORRECT!', IN0(lp0), i7, IS0(lp0)
     r_meter_all( i,mp) =   r_meter_all(i7,mp)
     SL_meter_all(i,mp) =  SL_meter_all(i7,mp)
     BM_nT_all(   i,mp) =     BM_nT_all(i7,mp)
     Qvalue_all(  i,mp) =  Qvalue_all(i7,mp)

     if ( counter == IPDIM5 ) JMAX_IS_all(mp,lp)=i

   else !if ( counter > IPDIM5 )
if ( mp==1 ) print *,'assign invalid value',i,lp,mp
      r_meter_all(i,mp) = invalid_value
     SL_meter_all(i,mp) = invalid_value
     BM_nT_all(   i,mp) = invalid_value
     Qvalue_all(  i,mp) = invalid_value
   end if

!original SH
if ( i>=midpoint ) then
  tmp(1,  i-2) = GCOLAT_all( i,mp)
  tmp(2,  i-2) = GLON_all(   i,mp)
  tmp(3,  i-2) = GL_rad_all(   i,mp)

  tmp1(1, 1:3,i-2) = D1_all(1:3,i,mp)
  tmp1(2, 1:3,i-2) = D2_all(1:3,i,mp)
  tmp1(3, 1:3,i-2) = D3_all(1:3,i,mp)
  tmp1(4, 1:3,i-2) = E1_all(1:3,i,mp)
  tmp1(5, 1:3,i-2) = E2_all(1:3,i,mp)
end if



     
end do

!assign tmp values
i0=midpoint-1
i1=IS0(lp)-2  !original IS
GCOLAT_all(i0:i1,mp) = tmp(1, i0:i1)
GLON_all(  i0:i1,mp) = tmp(2, i0:i1)
GL_rad_all(i0:i1,mp) = tmp(3, i0:i1)

! the last two points
GCOLAT_all(i1+1:i1+2,mp) = invalid_value
GLON_all(  i1+1:i1+2,mp) = invalid_value
GL_rad_all(i1+1:i1+2,mp) = invalid_value

do jth=1,3
  D1_all(jth,i0:i1,mp) = tmp1(1, jth,i0:i1)
  D2_all(jth,i0:i1,mp) = tmp1(2, jth,i0:i1)
  D3_all(jth,i0:i1,mp) = tmp1(3, jth,i0:i1)
  E1_all(jth,i0:i1,mp) = tmp1(4, jth,i0:i1)
  E2_all(jth,i0:i1,mp) = tmp1(5, jth,i0:i1)

! the last two points
  D1_all(jth,i1+1:i1+2,mp) = invalid_value
  D2_all(jth,i1+1:i1+2,mp) = invalid_value
  D3_all(jth,i1+1:i1+2,mp) = invalid_value
  E1_all(jth,i1+1:i1+2,mp) = invalid_value
  E2_all(jth,i1+1:i1+2,mp) = invalid_value

end do !jth

end do  !lp
end do  !mp


! CORRECTED data
write(UNIT=2000,FMT="('Z_km',10i13)") (lp, lp=lp_start,lp_stop)
write(UNIT=2001,FMT="('SL_meter',10i13)")  (lp, lp=lp_start,lp_stop)
write(UNIT=2002,FMT="('Bm_nT',10i13)")  (lp, lp=lp_start,lp_stop)
write(UNIT=2003,FMT="('GL_rad',10i13)") (lp, lp=lp_start,lp_stop) 
write(UNIT=2004,FMT="('JMIN JMAX',10i13)")  (lp, lp=lp_start,lp_stop)
write(UNIT=2005,FMT="('GCOLAT',10i13)") (lp, lp=lp_start,lp_stop) 

mp=1
IPDIM5 = IS0(lp_start) - IN0(lp_start) +1  !
DO i=1,IPDIM5 !IN0(5),IS0(5)

  ii = i-1
  do  lp=lp_start,lp_stop
    jj(lp) = ii+JMIN_IN_all(mp,lp)
    if ( jj(lp)>IS0(lp) )  jj(lp)=IS0(lp)
  end do

write(UNIT=2000,FMT="(i8,8E13.4)") i, ((r_meter_all(jj(lp),mp)-earth_radius)*1.0E-3, lp=lp_start,lp_stop)
write(UNIT=2001,FMT="(i8,8E13.4)") i, (SL_meter_all(jj(lp),mp), lp=lp_start,lp_stop) 
write(UNIT=2002,FMT="(i8,8E13.4)") i, ( BM_nT_all(jj(lp),mp), lp=lp_start,lp_stop)

write(UNIT=2003,FMT="(i8,8F10.4)") i, ( (90.-GL_rad_all(jj(lp),mp)*180./pi), lp=lp_start,lp_stop)
write(UNIT=2005,FMT="(i8,8F10.4)") i, ( (90.-GCOLAT_all(jj(lp),mp)*180./pi), lp=lp_start,lp_stop)
END DO

write(UNIT=2004,FMT="(4i8)") (lp,JMIN_IN_all(mp,lp),JMAX_IS_all(mp,lp), (JMAX_IS_all(mp,lp)-JMIN_IN_all(mp,lp)+1), lp=lp_start,lp_stop)

!---
!output the CORRECTED grid to a file...
      if ( sw_newQ==0 ) then
        filename =filepath_pgrid//filename_pgrid//'CORRECTED20120112'
      else if ( sw_newQ==1 ) then 
        filename =filepath_pgrid//filename_pgrid//'CORRECTED20120112nQ'
      endif
print *,'sw_newQ',sw_newQ,filename

      FORM_dum ='formatted' 
      STATUS_dum ='UNKNOWN'
      CALL open_file ( filename, LUN_pgrid, FORM_dum, STATUS_dum ) 
 
      print *,"MODIFIED grid open file completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) JMIN_IN_all, JMAX_IS_all  !IN_2d_3d , IS_2d_3d
      print *,"WRITEing JMIN_IN etc completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) r_meter_all,GCOLAT_all, GLON_all, Qvalue_all !gr_2d, gcol_2d, glon_2d, q_coordinate_2d
      print *,"WRITEing r_meter etc completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) GL_rad_all          !bcol_2d
      print *,"WRITEing GL_rad etc completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) SL_meter_all, BM_nT_all !integral_ds_2d, apex_BMAG_2d
      print *,"WRITEing SL_meter etc completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) D1_all, D2_all, D3_all      !Apex_D1_2d
      print *,"WRITEing D1-3 etc completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) E1_all, E2_all          !Apex_E1_2d
      print *,"WRITEing E1/2 etc completed"
      WRITE (UNIT=LUN_pgrid, FMT=*) Be3_all(1,1:NMP,1:NLP),Be3_all(2,1:NMP,1:NLP) !Apex_BE3_N
      print *,"WRITEing Be3 etc completed"
      CLOSE(UNIT=LUN_pgrid)
      print *,"MODIFIED global grid WRITEing finished, file closed..."
end if !( sw_test_grid_correct ) then 

STOP

END SUBROUTINE test_grid_output

!---------------------------
END MODULE module_FIELD_LINE_GRID_MKS

!20120112:copied from module_io.3d.f90
!-------
        SUBROUTINE open_file ( filename_dum, UNIT_dum, FORM_dum, STATUS_dum )  
        USE module_precision
!20120112:        USE module_IO,ONLY: LUN_LOG
        IMPLICIT NONE
!from module_io
      INTEGER (KIND=int_prec), PARAMETER :: LUN_LOG=9
        CHARACTER (LEN=*), INTENT(IN) :: filename_dum
        INTEGER (KIND=int_prec), INTENT(IN) :: UNIT_dum
        CHARACTER (LEN=*), INTENT(IN) :: FORM_dum
        CHARACTER (LEN=*), INTENT(IN) :: STATUS_dum
!---local
        LOGICAL :: flag
        INTEGER (KIND=int_prec) :: istat

IF (UNIT_dum>=1) THEN
  print *,'unit number=',UNIT_dum
ELSE
  print *,'!STOP! unit number not provided!!!'
  STOP
END IF

        INQUIRE ( UNIT=UNIT_dum, OPENED=flag ) 
WRITE( UNIT=LUN_LOG, FMT=*) flag
        IF ( .NOT. flag ) THEN
WRITE( UNIT=LUN_LOG, FMT=*) 'opening file:',filename_dum
          OPEN(UNIT=UNIT_dum,FILE=TRIM(filename_dum),STATUS=STATUS_dum,FORM=FORM_dum,IOSTAT=istat)
          IF ( istat /= 0 ) THEN
            WRITE( UNIT=LUN_LOG, FMT=*)'ERROR OPENING FILE',filename_dum
            STOP
          END IF
        END IF      

        END SUBROUTINE open_file
!---
