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
SUBROUTINE flux_tube_solver ( utime,mp,lp )
  USE module_precision
  USE module_IPE_dimension,ONLY: ISPEC,ISPEV,IPDIM
  USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,Pvalue,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,plasma_3d,ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3,TN_k,TINF_k,un_ms1,mlon_rad
  USE module_input_parameters,ONLY: time_step,F107D_new,F107_new,DTMIN_flip  &
  &, sw_INNO,FPAS_flip,HEPRAT_flip,COLFAC_flip,sw_IHEPLS,sw_INPLS,sw_debug,iout, start_time  &
  &, HPEQ_flip, sw_wind_flip, sw_depleted_flip, start_time_depleted &
  &, sw_output_fort167,mpfort167,lpfort167 &
  &, sw_neutral_heating_flip, ip_freq_output, parallelBuild,mype &
  &, sw_ctip_input,utime0lpi,lpi,input_params_begin,input_params_interval,solar_forcing_time_step
  USE module_physical_constants,ONLY: pi,zero
  USE module_IO,ONLY: PRUNIT,LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4
  USE module_unit_conversion,ONLY: M_TO_KM
  USE module_heating_rate,ONLY: get_neutral_heating_rate
  USE module_deplete_flux_tube,ONLY: deplete_flux_tube
  USE module_magfield,ONLY:sunlons
!
  IMPLICIT NONE
!------------------------
  INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
  INTEGER (KIND=int_prec), INTENT(IN) :: mp    !longitude
  INTEGER (KIND=int_prec), INTENT(IN) :: lp    !latitude
  REAL (KIND=REAL_prec) :: ltime !local time [hour]

!--- for CTIPINT
  INTEGER :: IN,IS !.. lcv + spatial grid indices
  INTEGER JMINX,JMAXX !.. lcv + spatial grid indices
  INTEGER CTIPDIM         !.. CTIPe array dimension, must equal to FLDIM
  DOUBLE PRECISION ::  PCO
  INTEGER ::  INNO            !.. switch to turn on FLIP NO calculation IF <0
  DOUBLE PRECISION, DIMENSION(IPDIM) ::  ZX,GLX,SLX,BMX,GRX,OX,HX,N2X,O2X,HEX,N4SX,NNOX &
  &, UNX &  ! field aligned component of neutral wind velocity [m s-1], positive SOUTHward
  &, TNX &
    !.. TINFX has to be an array for grazing incidence column densities
  &, TINFX & !.. Exospheric temperature [k]
  &, SZA_dum
  DOUBLE PRECISION :: DTMIN !.. Minimum time step allowed (>=10 secs?)
  DOUBLE PRECISION DT

  REAL ::  F107D_dum,F107A_dum         !.. daily and 81 day average F10.7      !..
  DOUBLE PRECISION ::  FPAS,HEPRAT,COLFACX
  DOUBLE PRECISION HPEQ

  INTEGER ::  IHEPLS,INPLS  !.. switches He+ and N+ dIFfusive solutions on
  DOUBLE PRECISION &
    !.. EHTX(3,J) = e heating rate, EHTX(1,J) = ion heating rate, EHTX(2,J) unUSEd
  &  EHTX(3,IPDIM) &
    !.. TE_TI(3,J) = Te, TE_TIX(2,J) = Ti = TE_TIX(2,J)
  & ,TE_TIX(3,IPDIM) &
  & ,XIONNX(ISPEC,IPDIM),XIONVX(ISPEC,IPDIM) &
  & ,NHEAT(IPDIM) &  !.. Neutral heating rate [eV/cm^3/s]
  & ,hrate_cgs(22,IPDIM)   !.. heating rates [eV/cm^3/s]

!nm20121020       REAL(KIND=REAL_prec), DIMENSION(7,MaxFluxTube,NLP,NMP) :: hrate_mks !.. each component of the Neutral heating rate (eV/kg/s)
!nm20121020      REAL(KIND=REAL_prec) :: min_hrate,max_hrate

  INTEGER EFLAG(11,11)    !.. error flags, check =0 on return from FLIP
  INTEGER :: PRUNIT_dum !.. Unit number to PRINT results
  INTEGER (KIND=int_prec) :: midpoint
!      INTEGER (KIND=int_prec) :: stat_alloc
  INTEGER (KIND=int_prec) :: ipts,i,ret
  INTEGER (KIND=int_prec),PARAMETER :: ip_freq_output_fort=900
  INTEGER (KIND=int_prec) :: jth !dbg20120501
  REAL :: mlt
!----------------------------------

  IN = JMIN_IN(lp)
  IS = JMAX_IS(lp)

! make sure that JMINX equals to 1
  JMINX   = IN - IN + 1
  JMAXX   = IS - IN + 1
  CTIPDIM = IS - IN + 1

  ZX(1:CTIPDIM)  = plasma_grid_Z(IN:IS,lp) * M_TO_KM !convert from m to km
  PCO = Pvalue(lp)  !Pvalue is a single value
  SLX(1:CTIPDIM) = plasma_grid_3d(IN:IS,lp,mp,ISL)
  GLX(1:CTIPDIM) = pi/2. - plasma_grid_GL(IN:IS,lp)  ! magnetic latitude [radians]
  BMX(1:CTIPDIM) = plasma_grid_3d(IN:IS,lp,mp,IBM)   !Tesla
  GRX(1:CTIPDIM) = plasma_grid_3d(IN:IS,lp,mp,IGR)

  OX(1:CTIPDIM) = ON_m3(IN:IS,lp,mp) !(m-3)
  HX(1:CTIPDIM) = HN_m3(IN:IS,lp,mp)
  N2X(1:CTIPDIM) = N2N_m3(IN:IS,lp,mp)
  O2X(1:CTIPDIM) = O2N_m3(IN:IS,lp,mp)
  HEX(1:CTIPDIM) = HE_m3(IN:IS,lp,mp)
  N4SX(1:CTIPDIM) = N4S_m3(IN:IS,lp,mp)

  INNO = sw_INNO

  TNX(1:CTIPDIM) = TN_k(IN:IS,lp,mp)
  TINFX(1:CTIPDIM) = TINF_k(IN:IS,lp,mp)

! FLIP expects positive SOUTHward along a field line
  UNX(1:CTIPDIM)  = (-1.) * Un_ms1(IN:IS,lp,mp,3)

!nm20160420: i am not totally sure about the FLIP time step???
!      DT        = REAL(time_step)
  DT        = REAL(solar_forcing_time_step) !nm20160420
  DTMIN     = DTMIN_flip
  IF ( sw_ctip_input ) THEN
    LPI = INT( ( utime - utime0LPI ) / REAL(input_params_interval) ) + 1 + input_params_begin
!t        IF(sw_debug)
    PRINT*,'sub-eld: LPI=',lpi
!t        IF(sw_debug)
    PRINT*,'sub-eld: utime',utime,'dt_m=',((utime-utime0LPI)/60.)
  else
    LPI=1
  ENDIF
  F107D_dum = F107_new(lpi)
  F107A_dum = F107d_new(lpi)

!nm20110822: moved from module_plasma
!! I need to get the new get_sza based on phil's method
!! calculate Solar Zenith Angle [radians]
  CALL Get_SZA ( utime,mp,lp, SZA_dum )
!      SZA_dum(JMINX:JMAXX)   = SZA_rad(IN:IS)
!nm20110822: no more allocatable arrays
!      IF ( ALLOCATED( SZA_rad  ) )  DEALLOCATE( SZA_rad, STAT=stat_alloc  )
!IF ( stat_alloc/=0 ) THEN
!  PRINT *, ALLOCATED( sza_rad )
!  PRINT *,"!STOP! SZA_rad DEALLOCATION FAILED! in flux_tube_solver:",stat_alloc,mp,lp,in,is,jminx,jmaxx
!  STOP
!ENDIF

  FPAS      = FPAS_flip

!nm110510: test the depleted flux tube
!nm20140729: moved to module_deplete_flux_tube.f90
  CALL deplete_flux_tube ( utime, mp,lp, HPEQ )
!      IF ( utime == start_time ) THEN
!        HPEQ      = HPEQ_flip
!      ELSE
!        HPEQ      = 0.0
!
!        IF ( sw_depleted_flip==1 .AND. utime == start_time_depleted ) THEN
!          HPEQ      = - 0.1
!        ENDIF
!      ENDIF


  HEPRAT    = HEPRAT_flip
  COLFACX   = COLFAC_flip

  IHEPLS    = sw_IHEPLS
  INPLS     = sw_INPLS

!! array initialization
  EFLAG(:,:)=0

!dbg20110802: 3D multiple-lp run
!IF (mp==1 .and.lp==10) THEN
!  sw_debug=.true.
!else
!  sw_debug=.false.
!ENDIF


  IF ( sw_debug ) THEN

    PRINT *,'sub-flux_tube_solver'
    PRINT "('mp=',i6,' lp=',i6,' JMINX=',I6,' JMAXX=',I6)", mp,lp,jminx,jmaxx
    PRINT "('CTIPDIM=',I6)", ctipdim
    PRINT "('Z [km]     =',2F10.4)", ZX(jminx),ZX(jmaxx)
    PRINT "('PCO        =',2F10.4)", PCO,Pvalue(lp)
    PRINT "('SLX [m]    =',2E12.4)", SLX(jminx), SLX(jmaxx)
    PRINT "('GLX [deg]  =',2F10.4)",(GLX(jminx)*180./pi),(GLX(jmaxx)*180./pi)
    PRINT "('BMX [Tesla]    =',2E12.4)", BMX(jminx), BMX(jmaxx)
    PRINT "('GRX[m2 s-1]=',2E12.4)",GRX(jminx),GRX(jmaxx)
!---neutral parameters
    PRINT "('LOG10 OX [m-3]     =',2F10.4)",LOG10(OX(jminx)),LOG10(OX(jmaxx))
    PRINT "('LOG10 HX [m-3]     =',2F10.4)",LOG10(HX(jminx)),LOG10(HX(jmaxx))
    PRINT "('LOG10 N2X [m-3]    =',2F10.4)",LOG10(N2X(jminx)),LOG10(N2X(jmaxx))
    PRINT "('LOG10 O2X [m-3]    =',2F10.4)",LOG10(O2X(jminx)),LOG10(O2X(jmaxx))
    PRINT "('LOG10 HEX [m-3]    =',2F10.4)",LOG10(HEX(jminx)),LOG10(HEX(jmaxx))
    PRINT "('LOG10 N4SX [m-3]   =',2F10.4)",LOG10(N4SX(jminx)),LOG10(N4SX(jmaxx))

    PRINT "('INNO =',I6)",INNO
    IF ( INNO>=0 )  & !when CTIPe calculates NO
    & PRINT "('LOG10 NNOX [m-3]   =',2F10.4)",LOG10(NNOX(jminx)),LOG10(NNOX(jmaxx))
    PRINT "('Tn [K]       =',2F10.4)",TNX(jminx),TNX(jmaxx)
    PRINT "('TINF [K]     =',2F10.4)",TINFX(jminx),TINFX(jmaxx)
    PRINT "('UNX [m s-1]  =',2F10.4)",UNX(jminx),UNX(jmaxx)

    PRINT "('DT [sec]     =',F10.4)",DT
    PRINT "('DTMIN [sec]  =',F10.4)",DTMIN
    PRINT "('F107D_dum    =',F10.4)",F107D_dum
    PRINT "('F107A_dum    =',F10.4)",F107A_dum
    PRINT "('SZA [deg]    =',2F10.4)",SZA_dum(jminx)*180./pi,SZA_dum(jmaxx)*180./pi
    PRINT "('FPAS         =',F10.4)",FPAS
    PRINT "('HPEQ         =',F10.4)",HPEQ
    PRINT "('HEPRAT       =',F10.4)",HEPRAT
    PRINT "('COLFACX      =',F10.4)",COLFACX
    PRINT "('IHEPLS       =',I6)",IHEPLS
    PRINT "('INPLS        =',I6)",INPLS

  ENDIF !( sw_debug ) THEN

!dbg20120125:      midpoint = JMINX + (CTIPDIM-1)/2
  midpoint = IN + (IS-IN)/2
  IF ( lp>=1 .AND. lp<=6 )  midpoint = midpoint - 1
!nm20110909: calculating LT should be made FUNCTION!!!
  ltime = REAL(utime)/3600.0 + (plasma_grid_3d(midpoint,lp,mp,IGLON)*180.0/pi)/15.0
  IF ( ltime > 24.0 )  ltime = MOD(ltime, 24.0)

  IF( sw_output_fort167.AND.mp==mpfort167.AND.lp==lpfort167 ) THEN
!sms$ignore begin
    WRITE(UNIT=LUN_FLIP1,FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)") mp,lp,LUN_FLIP1,REAL(UTIME)/3600., ltime
    WRITE(UNIT=LUN_FLIP2,FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)") mp,lp,LUN_FLIP2,REAL(UTIME)/3600., ltime
    WRITE(UNIT=LUN_FLIP3,FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)") mp,lp,LUN_FLIP3,REAL(UTIME)/3600., ltime
    WRITE(UNIT=LUN_FLIP4,FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)") mp,lp,LUN_FLIP4,REAL(UTIME)/3600., ltime
!sms$ignore end
  ENDIF !( sw_output_fort167
  IF(sw_debug) PRINT*,'sub-fl: UTs=',UTIME,' LThr=',ltime,' mp',mp,' lp',lp


!dbg20110131:
  IF ( sw_debug ) THEN
!sms$ignore begin
    WRITE(UNIT=PRUNIT,FMT="('mp=',i6,' lp=',i6,' UT=',F10.2,i7)") mp,lp,REAL(UTIME)/3600.,mype
!sms$ignore end
  ENDIF! ( sw_debug ) THEN
  DO ipts=1,CTIPDIM
!N&T from the previous time step are absolute necesary for the solver...
!dbg20120501
    DO jth=1,ISPEC
      XIONNX(jth,ipts) = plasma_3d(ipts,lp,mp,jth)
    ENDDO !jth
!te
    TE_TIX(3,ipts) = plasma_3d(ipts,lp,mp,ISPEC+1)
!ti
    DO jth=1,2
      TE_TIX(jth,ipts) = plasma_3d(ipts,lp,mp,jth+ISPEC+1)
    ENDDO !jth
!vi
    DO jth=1,ISPEC
      IF ( jth<=2 ) THEN
        XIONVX(jth,ipts) = plasma_3d(ipts,lp,mp,jth+ISPEC+3)
      ELSE
        XIONVX(jth,ipts) = zero
      ENDIF
    ENDDO !jth



    EHTX(  1:3        ,ipts)=zero  !dbg20110927


    IF ( INNO<0 ) THEN   !when flip calculates NO
      NNOX(              ipts)=zero !dbg20110927
    ELSE 
      PRINT *,'CTIPe calculates NO'
    ENDIF
    NHEAT(             ipts)=zero !dbg20110927
    hrate_cgs(1:22,    ipts)=zero !nm20121020

  ENDDO 

  IF ( sw_ihepls==0 ) XIONNX(3,:) = zero
  IF ( sw_inpls==0  ) XIONNX(4,:) = zero

  mlt = mlon_rad(mp)*180./pi/15.0D0-sunlons(1)*12.0D0/pi+12.0 !;[hr]

! call the flux tube solver (FLIP)
  CALL CTIPINT( &
  &             JMINX, & !.. index of the first point on the field line
  &             JMAXX, & !.. index of the last point on the field line
  &           CTIPDIM, & !.. CTIPe array dimension, must equal to FLDIM
  &                ZX, & !.. array, altitude (km)
  &               PCO, & !.. p coordinate (L-shell)
  &               SLX, & !.. array, distance of point from northern hemisphere (meter)
  &               GLX, & !.. array, magnetic latitude (radians)
  &               BMX, & !.. array, magnetic field strength, (Tesla)
  &               GRX, & !.. array, gravity, m2 s-1
  &                OX, & !.. array, O density (m-3)
  &                HX, & !.. array, H density (m-3)
  &               N2X, & !.. array, N2 density (cm-3)
  &               O2X, & !.. array, O2 density (cm-3)
  &               HEX, & !.. array, He density (cm-3)
  &              N4SX, & !.. array, N(4S) density (cm-3)
  &              INNO, & !.. switch to turn on FLIP NO calculation IF <0
  &              NNOX, & !.. array, NO density (cm-3)
  &               TNX, & !.. array, Neutral temperature (K)
  &             TINFX, & !.. array, Exospheric Neutral temperature (K)
  &               UNX, & !.. array, Neutral wind (m/s), field aligned component, positive SOUTHward
  &                DT, & !.. CTIPe time step (secs)
  &             DTMIN, & !.. Minimum time step allowed (>=10 secs?)
  &         F107D_dum, & !.. Daily F10.7
  &         F107A_dum, & !.. 81 day average F10.7
  &           SZA_dum, & !.. Solar Zenith angle (radians)
  &              FPAS, & !.. Pitch angle scattering fraction
  &              HPEQ, & !.. Sets initial equatorial H+ density. See declaration below
  &            HEPRAT, & !.. Intial He+/H+ ratio (.01 to 1.0)
  &           COLFACX, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
  &            IHEPLS, & !.. switches He+ dIFfusive solution on IF > 0
  &             INPLS, & !.. switches N+ dIFfusive solution on IF > 0
  &              EHTX, & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
  &            TE_TIX, & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
  &     XIONNX,XIONVX, & !.. IN/OUT: 2D array, Storage for ion densities and velocities
  &             NHEAT, & !.. OUT: array, Neutral heating rate (eV/cm^3/s)
  &             EFLAG, & !.. OUT: 2D array, Error Flags
  &                mp, &
  &                lp, &
  &                utime, & !dbg20141209
  &         hrate_cgs,mlt  ) !.. heating rates [eV/cm^3/s] !nm20121020





! output
  DO ipts=1,CTIPDIM

    DO jth=1,ISPEC
      plasma_3d(ipts+IN-1,lp,mp,jth) = XIONNX(jth,ipts)
    ENDDO !jth


!dbg20160421 debug
    IF( EFLAG(2,1)/=0 .and. 32<(ipts+IN-1) .and. (ipts+IN-1)<39 )THEN
!SMS$IGNORE begin
      PRINT"(i2,' XION=',e10.2,f7.0,f7.1,i4,' lp=',i3,' mp=',i2,' LT=',f6.1)",mype,plasma_3d(ipts+IN-1,lp,mp,1),(plasma_grid_Z(ipts+IN-1,lp)*1.e-3),((pi/2. - plasma_grid_GL(ipts+IN-1,lp))*180./pi),(ipts+IN-1),lp,mp,ltime
!SMS$IGNORE end
    ENDIF!(EFLAG




!te
    plasma_3d(ipts+IN-1,lp,mp,ISPEC+1) = TE_TIX(3,ipts)
!ti
    DO jth=1,2
      plasma_3d(ipts+IN-1,lp,mp,jth+ISPEC+1) = TE_TIX(jth,ipts)
    ENDDO !jth
!vi
    DO jth=1,ISPEV
      plasma_3d(ipts+IN-1,lp,mp,jth+ISPEC+3) = XIONVX(jth,ipts)
    ENDDO !jth


!nm20110404: save each component of heating rate for output
    IF ( sw_neutral_heating_flip==1 .AND. &
    &  MOD( (utime-start_time),ip_freq_output)==0) THEN
      IF(parallelBuild) THEN
        PRINT*,'sw_neutral_heating_flip=1 DOes not work in parallel'
        PRINT*,'STOPping in Neut_heating'
        STOP
      ENDIF

      CALL get_neutral_heating_rate ( hrate_cgs , lp,mp )

    ENDIF !( sw_neutral_heating_flip==1 ) THEN

  ENDDO       !DO ipts=1,CTIPDIM

  PRUNIT_dum = PRUNIT
  CALL WRITE_EFLAG(PRUNIT_dum, &  !.. Unit number to PRINT results
  &                      EFLAG, &  !.. Error flag array
  &                         mp, &
  &                         lp,utime,ltime )



END SUBROUTINE flux_tube_solver
