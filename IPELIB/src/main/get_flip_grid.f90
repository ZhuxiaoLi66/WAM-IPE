!20120223; copied from main.bk/module_field_line_grid.1d.f90
! used only for debug purpose
! CAUTION! apexD apexE Be3 cannot be calculated, 
! thus do not use when calculating electrodynamics!
! also parallel component of the neutral wind requires apexD in this version 
!
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
      SUBROUTINE get_FLIP_grid ( mp,lp )
      USE module_precision
      USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_3d,plasma_grid_GL,plasma_grid_Z,mlon_rad,ht90,JMIN_IN,JMAX_IS,Pvalue,apexD,east,north,up,ISL,IBM,IGR,IQ,IGCOLAT,IGLON
      USE module_input_parameters,ONLY:NYEAR,NDAY,sw_debug,iout  
!, PCO_flip,BLON_flip  !20120223
      USE module_physical_constants,ONLY: pi
      USE module_unit_conversion,ONLY:M_to_CM
IMPLICIT NONE
!---INPUT
      INTEGER (KIND=int_prec),INTENT(IN) :: mp,lp
!---local
      DOUBLE PRECISION  PCO,BLON,Z0,SCAL
      INTEGER JMAX,IDAY
      ! FLIP GRID
      INTEGER (KIND=int_prec), PARAMETER :: FLDIM0=9001     !.. Array dimensions

      DOUBLE PRECISION RE,VTOT
      DOUBLE PRECISION GLATD(FLDIM0),GLOND(FLDIM0)
      DOUBLE PRECISION COSDIP(FLDIM0),Q(FLDIM0),GL(FLDIM0)
      DOUBLE PRECISION Z(FLDIM0),SL(FLDIM0),BM(FLDIM0),GR(FLDIM0),R(FLDIM0)
      DOUBLE PRECISION SINDIP(FLDIM0)
      INTEGER (KIND=int_prec) :: midpoint,istp,i
      INTEGER (KIND=int_prec) :: IN,IS
      REAL (KIND=real_prec) :: dotprod

      JMAX=-300
      IDAY=NYEAR*1000 + NDAY      !.. Day of year in the form YYYYddd	 
      Z0=ht90*1.0E-3        ![km].. Lower boundary altitude
      SCAL=-1.0      !.. SCAL distributes the grid points along B
!      PCO=PCO_flip   !Pvalue       !.. L-shell or geographic latitude
PCO=Pvalue(lp)
!      BLON=BLON_flip ![deg?].. Magnetic if PCO is L-shell (or geographic) longitude
BLON=mlon_rad(mp)*180.0/pi

      !.. # of grid points on field line (set negative for default)
      print *,'FLDIM0=',FLDIM0,'JMAX',JMAX,'IDAY',IDAY,'PCO',PCO,'BLON',BLON,'Z0',Z0,'SCAL',SCAL

      !..  Call the Field line grid routine
      CALL FLIP_GRID(FLDIM0,JMAX,IDAY,PCO,BLON,Z0,SCAL,RE,GL,Q,COSDIP, &
     &   SINDIP,Z,SL,GLATD,GLOND,BM,GR,R,VTOT)


      print "('L NH',3E12.4)",SL(1),SL(2),SL(3)
      print "('SL SH',3E12.4)",SL(JMAX),SL(JMAX-1),SL(JMAX-2)
      
!assing the FLIP grid to the IPE arrays
! make sure that JMIN_IN starts from 1            
!      JMIN_IN(lp) = +1
     JMAX_IS(lp) = JMIN_IN(lp) + JMAX -1
     print *,'JMAX=',JMAX,' new JMAX_IS=',JMAX_IS(lp)

     IN = JMIN_IN(lp)
     IS = JMAX_IS(lp)


!      midpoint = 1 + (JMAX_IS/2)+1 -1 
      midpoint = IN + ( IS - IN )/2
      print *,'midpoint',midpoint,plasma_grid_Z(midpoint,lp)

! make sure to use the RIGHT SI UNIT!!!
plasma_grid_Z (in:is,lp)  =  Z(1:JMAX)* 1.0E+3    !convert from km to m
plasma_grid_3d(IN:IS,lp,mp,IGLON) =  GLOND(1:JMAX)*pi/180.0 !convert from [deg] to [rad]
plasma_grid_3d(in:is,lp,mp,IQ)  =  Q(1:JMAX)
plasma_grid_3d(in:is,lp,mp,ISL) = SL(1:JMAX)/M_to_CM  !cm-->m

      print "('SL_meter NH',3E12.4)",plasma_grid_3d(in:in+2,lp,mp,ISL)
      print "('SL_meter SH',3E12.4)",plasma_grid_3d(is-2:is,lp,mp,ISL)

plasma_grid_3d(in:is,lp,mp,IBM)  =  BM(1:JMAX)*1.0E-4    !from gauss to tesla
plasma_grid_GL(in:is,lp) = pi*0.5 - GL(1:JMAX)   ![rad]


plasma_grid_3d(in:is,lp,mp,IGCOLAT)  =  (90.0-GLATD(1:JMAX))*pi/180.0 !convert from lat[deg] to CO-LAT[rad]
plasma_grid_3d(in:is,lp,mp,IGR)  =  GR(1:JMAX)/M_to_CM  !cm2 --> m2

!debug write
!IF ( sw_debug ) THEN
print *,'FLIP GRID',iout(1),iout(2),pco,blon
print "('JMIN(IN)=',i6,'  JMAX(IS)=',i6,'  midpoint=',i5)", IN,IS, midpoint

print "('G-LAT [deg]=NH',21f10.4)",(90.-plasma_grid_3d(in:in+20,lp,mp,IGCOLAT)*180./pi)
print "('G-LAT [deg]=SH',21f10.4)",(90.-plasma_grid_3d(is-20:is,lp,mp,IGCOLAT)*180./pi)

print "('M-LAT [deg]=',2f10.4)",(90.-plasma_grid_GL(in,lp)*180./pi),(90.-plasma_grid_GL(is,lp)*180./pi)
print "('GLON  [deg]=',2f10.4)",(plasma_grid_3d(IN,lp,mp,IGLON)*180./pi),( plasma_grid_3d(IS,lp,mp,IGLON)*180./pi)
print "('Qvalue     =',2E12.4)", plasma_grid_3d(in,lp,mp,IQ ) , plasma_grid_3d(is,lp,mp,IQ )
print "('BM [nT]    =',2E12.4)", plasma_grid_3d(in,lp,mp,IBM) , plasma_grid_3d(is,lp,mp,IBM)

print "('SL [m]     =',4E13.5)", plasma_grid_3d(in:in+1,lp,mp,ISL), plasma_grid_3d(is-1:is,lp,mp,ISL)

print "('Z  [m]     =',4E13.5)",  plasma_grid_Z(in:in+1,lp),  plasma_grid_Z(is-1:is,lp)

!END IF !( sw_debug ) THEN


print "('SL [m]     =',4E13.5)", plasma_grid_Z(in,lp), plasma_grid_Z(is,lp)



!do i=1,jmax
!   istp=JMIN_IN(lp)+i-1

   apexD(in:is,lp,mp,north,3) = COSDIP(1:JMAX)
   apexD(in:is,lp,mp,up   ,3) = SINDIP(1:JMAX)

!end do !i=1,jmax

print *,'sub-get_FLIP_grid finished'
END SUBROUTINE get_FLIP_grid
