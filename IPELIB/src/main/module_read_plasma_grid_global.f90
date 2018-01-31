      MODULE module_read_plasma_grid_global

      CONTAINS
!---------------------------
!20110726: the new GLOBAL 3D version: NMP=80
! global grid with the low resolution version
! magnetic longitude used for the grid. from 0 to 355.5 with 4.5 degree interval

      SUBROUTINE read_plasma_grid_global
        USE module_precision
        USE module_IPE_dimension,ONLY: NMP,NLP,ISTOT
        USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
        USE module_physical_constants,ONLY: earth_radius,pi,zero
        USE module_input_parameters,ONLY:read_input_parameters,sw_debug,sw_neutral_heating_flip,nprocs,mype
        USE module_IO,ONLY: filename,LUN_pgrid
        USE module_FIELD_LINE_GRID_MKS,ONLY: &
&           JMIN_IN_all,JMAX_IS_all          &
&,          JMIN_ING   ,JMAX_ISG             &
&,          JMIN_IN    ,JMAX_IS              &
&,          r_meter2D, plasma_grid_GL,plasma_grid_3d,apexD,apexE,Be3,plasma_grid_Z &
&,          ISL,IBM,IGR,IQ,IGCOLAT,IGLON,east,north,up                             &
&,          MaxFluxTube,minTheta,maxTheta,minAltitude,maxAltitude,midpnt,plasma_3d
        USE module_open_file,ONLY: open_file
        IMPLICIT NONE

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
      INTEGER(KIND=int_prec) lp,mp,stat_alloc,midpoint_min,midpoint_max

      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum0    !.. distance from the center of the Earth[meter]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum1    !.. geographic co-latitude [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum2    !.. geographic longitude [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum3
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  GL_rad_all    !.. magnetic co-latitude Eq(6.1) [rad]
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  SL_meter_all  !.. distance of point from northern hemisphere foot point [meter]
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  BM_T_all      !.. magnetic field strength [T]
! components (east, north, up) of base vectors
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  dum4    !.. Eq(3.8) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  dum5    !.. Eq(3.9) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  dum6    !.. Eq(3.10) Richmond 1995
!dbg20110927      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  E1_all    !.. Eq(3.11) Richmond 1995
!dbg20110927      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  E2_all    !.. Eq(3.12) Richmond 1995
!JFM  REAL(KIND=real_prec), DIMENSION(2,NMP,NLP) ::  Be3_all         ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [T]
!SMS$DISTRIBUTE(dh,NLP,NMP) BEGIN
      REAL(KIND=real_prec), DIMENSION(NMP,NLP) ::  Be3_all1,Be3_all2 ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [T]: "Ed1, Ed2, and Be3 are constant along magnetic field lines"
!SMS$DISTRIBUTE END

!-------------local
        CHARACTER (LEN=11) :: FORM_dum
        CHARACTER (LEN=7)  :: STATUS_dum
        CHARACTER(LEN=*), PARAMETER :: filepath_pgrid= &
!     & '../../field_line_grid/20110419lowres_global/'  !20110419 low res global grid
     &  './'
        CHARACTER(LEN=*), PARAMETER :: filename_pgrid= &
!     &  'GIP_apex_coords_global_lowresCORRECTED'        !20110824: corrected for lp=1-6
     &  'ipe_grid'        !20110824: corrected for lp=1-6

!dbg      INTEGER(KIND=int_prec) :: sw_test_grid=0  !1: ON testgrid; 0: OFF 
!---

        ALLOCATE ( JMIN_IN_all ( NMP,NLP  ) &
     &,            JMAX_IS_all ( NMP,NLP  ) &
     &,            JMIN_ING    (     NLP  ) &
     &,            JMAX_ISG    (     NLP  ) &
     &,            STAT=stat_alloc        )
 
      IF ( stat_alloc/=0 ) THEN
        print *,"!STOP! ALLOCATION of JMAX* and JMIN* FAILD!:",stat_alloc
        STOP
      END IF

!JMIN_IN_all and JMAX_IS_all are OUT variables to workaround an SMS bug
!SMS$SERIAL(<JMIN_ING,JMAX_ISG,MaxFluxTube,JMIN_IN_all,JMAX_IS_all,OUT> : default=ignore) BEGIN
      filename =filepath_pgrid//filename_pgrid
      FORM_dum ='formatted' 
      STATUS_dum ='old'
      CALL open_file ( filename, LUN_pgrid, FORM_dum, STATUS_dum ) 
      print *,"open file completed"
      READ (UNIT=LUN_pgrid, FMT=*) JMIN_IN_all, JMAX_IS_all  !IN_2d_3d , IS_2d_3d
      print *,"reading JMIN_IN etc completed"
      JMIN_ING = JMIN_IN_all(1,:)
      JMAX_ISG = JMAX_IS_all(1,:)
      MaxFluxTube = maxval(JMAX_ISG-JMIN_ING+1)
!SMS$SERIAL END

      DEALLOCATE ( JMIN_IN_all,JMAX_IS_all,STAT=stat_alloc )
      IF ( stat_alloc/=0 ) THEN
        print *,"!STOP! DEALLOCATION of JMIN_IN_all/JMAX_IS_all FAILD!",stat_alloc
        STOP
      END IF

      CALL allocate_arrays ( 0 )

      JMIN_IN = 1
      JMAX_IS = JMAX_ISG - JMIN_ING + 1
      write(88,*) NLP
      DO lp=1,NLP  !longest -->shortest flux tube
        midpnt(lp) = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
       write(88,*) lp,JMAX_IS(lp)-JMIN_IN(lp)
      END DO

!array initialization
!SMS$IGNORE BEGIN
      Be3            = zero
      plasma_grid_3d = zero
      plasma_grid_Z  = zero
      plasma_grid_GL = zero
      plasma_3d      = zero
      apexD          = zero
      apexE          = zero
      r_meter2D      = zero
!SMS$IGNORE END

!SMS$SERIAL(<r_meter2D,plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,OUT> : default=ignore) BEGIN
READ (UNIT=LUN_pgrid, FMT=*) dum0, dum1, dum2, dum3 !gr_2d, gcol_2d, glon_2d, q_coordinate_2d
do lp=1,NLP
  r_meter2D    (JMIN_IN(lp):JMAX_IS(lp),lp) = dum0(JMIN_ING(lp):JMAX_ISG(lp),1)                !r_meter
  plasma_grid_Z(JMIN_IN(lp):JMAX_IS(lp),lp) = dum0(JMIN_ING(lp):JMAX_ISG(lp),1) - earth_radius ![meter]
enddo
do mp=1,NMP
  do lp=1,NLP
    plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,mp,IGCOLAT) = dum1(JMIN_ING(lp):JMAX_ISG(lp),mp) !GCOLAT
    plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,mp,IGLON  ) = dum2(JMIN_ING(lp):JMAX_ISG(lp),mp) !GLON
    plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,mp,IQ     ) = dum3(JMIN_ING(lp):JMAX_ISG(lp),mp) !Q
  enddo
enddo
print *,"reading r_meter etc completed"

READ (UNIT=LUN_pgrid, FMT=*) dum0          !bcol_2d
do lp=1,NLP
  plasma_grid_GL(JMIN_IN(lp):JMAX_IS(lp),lp) = dum0(JMIN_ING(lp):JMAX_ISG(lp),1) !GL
enddo
print *,"reading GL_rad etc completed"
READ (UNIT=LUN_pgrid, FMT=*) dum0, dum1 !integral_ds_2d, apex_BMAG_2d
do lp=1,NLP
  plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,ISL) = dum0(JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !SL
  plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,IBM) = dum1(JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !BM
enddo
print *,"reading SL_meter etc completed"
!SMS$SERIAL END
!SMS$EXCHANGE(plasma_grid_3d)

minTheta=plasma_grid_GL(JMIN_IN(  1),  1)
maxTheta=plasma_grid_GL(JMIN_IN(NLP),NLP) 
!dbbg20120301: temporary solution in stepback_mag_R to keep flux tube within the sim region, instead of stopping
midpoint_min = JMIN_IN(NLP) + ( JMAX_IS(NLP) - JMIN_IN(NLP) )/2
midpoint_max = JMIN_IN(  1) + ( JMAX_IS(  1) - JMIN_IN(  1) )/2
minAltitude  = plasma_grid_Z(midpoint_min,NLP)
maxAltitude  = plasma_grid_Z(midpoint_max,  1)

!SMS$SERIAL(<apexD,OUT> : default=ignore) BEGIN
      READ (UNIT=LUN_pgrid, FMT=*) dum4, dum5, dum6      !Apex_D1_2d
!D2
!dbg20110923  apexD(1,1:NPTS2D,1:NMP)%east  =  dum4(1,1:NPTS2D,1:NMP) !D1
!dbg20110923  apexD(1,1:NPTS2D,1:NMP)%north =  dum4(2,1:NPTS2D,1:NMP)
!dbg20110923  apexD(1,1:NPTS2D,1:NMP)%up    =  dum4(3,1:NPTS2D,1:NMP)
!D2
!dbg20110923  apexD(2,1:NPTS2D,1:NMP)%east  =  dum5(1,1:NPTS2D,1:NMP) !D2
!dbg20110923  apexD(2,1:NPTS2D,1:NMP)%north =  dum5(2,1:NPTS2D,1:NMP)
!dbg20110923  apexD(2,1:NPTS2D,1:NMP)%up    =  dum5(3,1:NPTS2D,1:NMP)
!D3
do lp=1,NLP
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ,1) =  dum4(1,JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !D1
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north,1) =  dum4(2,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ,1) =  dum4(3,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)

  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ,2) =  dum5(1,JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !D2
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north,2) =  dum5(2,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ,2) =  dum5(3,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)

  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ,3) =  dum6(1,JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !D3
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north,3) =  dum6(2,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
  apexD(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ,3) =  dum6(3,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
enddo
      print *,"reading D1-3 etc completed"
!SMS$SERIAL END

!SMS$SERIAL(<apexE,OUT> : default=ignore) BEGIN
      READ (UNIT=LUN_pgrid, FMT=*) dum4, dum5          !Apex_E1_2d
!E1
do lp=1,NLP
  apexE(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ,1) =  dum4(1,JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !E1
  apexE(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north,1) =  dum4(2,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
  apexE(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ,1) =  dum4(3,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
!E2
  apexE(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ,2) =  dum5(1,JMIN_ING(lp):JMAX_ISG(lp),1:NMP) !E2
  apexE(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north,2) =  dum5(2,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
  apexE(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ,2) =  dum5(3,JMIN_ING(lp):JMAX_ISG(lp),1:NMP)
enddo
      print *,"reading E1/2 etc completed"
!SMS$SERIAL END

!JFM Be3_all1 and Be3_all2 are OUT variables to workaround an SMS limitation.
!SMS$SERIAL(<Be3,Be3_all1,Be3_all2,OUT> : default=ignore) BEGIN
!JFM  READ (UNIT=LUN_pgrid, FMT=*) Be3_all(1,1:NMP,1:NLP),Be3_all(2,1:NMP,1:NLP) !Apex_BE3_N
!JFM  READ (UNIT=LUN_pgrid, FMT=*) Be3_all(1,:,:),Be3_all(2,:,:) !Apex_BE3_N
      READ (UNIT=LUN_pgrid, FMT=*) Be3_all1,Be3_all2 !Apex_BE3_N
!JFM  Be3(1:2,1:NMP,1:NLP)=    Be3_all(1:2,1:NMP,1:NLP)

!nm20130830
!d print *, "!nm20130830: make sure Be3 is constant!? Be3_all1=", Be3_all1(1,130)," Be3_all2=",Be3_all2(1,130)

      do mp=1,NMP
        do lp=1,NLP
          Be3(lp,mp)=Be3_all1(mp,lp)
!nm20130830: Be3 is constant along a flux tube!
!nm20130830          Be3(2,lp,mp)=Be3_all2(mp,lp)
        enddo
      enddo
      print *,"reading Be3 completed"
      CLOSE(UNIT=LUN_pgrid)
      print *,"global grid reading finished, file closed..."
!SMS$SERIAL END

!dbg20110811:
!dbg IF ( sw_test_grid==1 ) THEN

!dbg  CALL test_grid_output ( JMIN_IN_all,JMAX_IS_all,r_meter_all,SL_meter_all,BM_T_all,GL_rad_all &
!dbg & ,GCOLAT_all,GLON_all,Qvalue_all,D1_all,D2_all,D3_all,E1_all,E2_all,Be3_all &
!dbg &,filepath_pgrid,filename_pgrid)
!dbg END IF 

      END SUBROUTINE  read_plasma_grid_global
!---
!20110927:test_grid_output deleted
!---------------------------
      END MODULE module_read_plasma_grid_global
