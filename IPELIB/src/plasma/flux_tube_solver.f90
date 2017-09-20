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
      USE module_input_parameters,ONLY: time_step,F107D,F107AV,DTMIN_flip  &
     &, sw_INNO,FPAS_flip,HEPRAT_flip,COLFAC_flip,sw_IHEPLS,sw_INPLS,sw_debug,iout, start_time  &
     &, HPEQ_flip, sw_wind_flip, sw_depleted_flip, start_time_depleted &
     &, sw_output_fort167,mpfort167,lpfort167 &
     &, sw_neutral_heating_flip, ip_freq_output, parallelBuild,mype,ip_freq_paraTrans
      USE module_PLASMA,ONLY: plasma_1d !dbg20120501
!dbg20110927      USE module_heating_rate,ONLY: NHEAT_cgs
      USE module_physical_constants,ONLY: pi,zero
      USE module_IO,ONLY: PRUNIT,LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4
      USE module_unit_conversion,ONLY: M_TO_KM
      USE module_heating_rate,ONLY: get_neutral_heating_rate
      USE module_deplete_flux_tube,ONLY: deplete_flux_tube
      USE module_magfield,ONLY:sunlons
!
      IMPLICIT NONE
      include "gptl.inc"
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec), INTENT(IN) :: mp    !longitude
      INTEGER (KIND=int_prec), INTENT(IN) :: lp    !latitude
      REAL (KIND=real_prec) :: ltime !local time [hour]

!--- for CTIPINT
      INTEGER :: IN,IS !.. lcv + spatial grid indices
      INTEGER JMINX,JMAXX !.. lcv + spatial grid indices
      INTEGER CTIPDIM         !.. CTIPe array dimension, must equal to FLDIM
      DOUBLE PRECISION ::  PCO
      INTEGER ::  INNO            !.. switch to turn on FLIP NO calculation if <0
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

      INTEGER ::  IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on 
      DOUBLE PRECISION &
      !.. EHTX(3,J) = e heating rate, EHTX(1,J) = ion heating rate, EHTX(2,J) unused
     &  EHTX(3,IPDIM) &
      !.. TE_TI(3,J) = Te, TE_TIX(2,J) = Ti = TE_TIX(2,J)
     & ,TE_TIX(3,IPDIM) &
     & ,XIONNX(ISPEC,IPDIM),XIONVX(ISPEC,IPDIM) &
     & ,NHEAT(IPDIM) &  !.. Neutral heating rate [eV/cm^3/s]
     & ,hrate_cgs(22,IPDIM)   !.. heating rates [eV/cm^3/s]

!nm20121020       REAL(KIND=real_prec), DIMENSION(7,MaxFluxTube,NLP,NMP) :: hrate_mks !.. each component of the Neutral heating rate (eV/kg/s) 
!nm20121020      REAL(KIND=real_prec) :: min_hrate,max_hrate

      INTEGER EFLAG(11,11)    !.. error flags, check =0 on return from FLIP
      INTEGER :: PRUNIT_dum !.. Unit number to print results
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
      DT        = REAL(ip_freq_paraTrans) !nm20160420
      DTMIN     = DTMIN_flip
      F107D_dum = F107D
      F107A_dum = F107AV

!nm20110822: moved from module_plasma
!! I need to get the new get_sza based on phil's method
!! calculate Solar Zenith Angle [radians]
      ret = gptlstart ('Get_SZA')
      CALL Get_SZA ( utime,mp,lp, SZA_dum )
      ret = gptlstop  ('Get_SZA')
!      SZA_dum(JMINX:JMAXX)   = SZA_rad(IN:IS)
!nm20110822: no more allocatable arrays
!      IF ( ALLOCATED( SZA_rad  ) )  DEALLOCATE( SZA_rad, STAT=stat_alloc  )
!if ( stat_alloc/=0 ) then
!  print *, ALLOCATED( sza_rad )
!  print *,"!STOP! SZA_rad DEALLOCATION FAILED! in flux_tube_solver:",stat_alloc,mp,lp,in,is,jminx,jmaxx
!  STOP
!endif

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
!if (mp==1 .and.lp==10) then
!  sw_debug=.true.
!else
!  sw_debug=.false.
!end if


IF ( sw_debug ) THEN

print *,'sub-flux_tube_solver'
print "('mp=',i6,' lp=',i6,' JMINX=',I6,' JMAXX=',I6)", mp,lp,jminx,jmaxx
print "('CTIPDIM=',I6)", ctipdim
print "('Z [km]     =',2F10.4)", ZX(jminx),ZX(jmaxx)
print "('PCO        =',2F10.4)", PCO,Pvalue(lp)
print "('SLX [m]    =',2E12.4)", SLX(jminx), SLX(jmaxx) 
print "('GLX [deg]  =',2F10.4)",(GLX(jminx)*180./pi),(GLX(jmaxx)*180./pi)
print "('BMX [Tesla]    =',2E12.4)", BMX(jminx), BMX(jmaxx)
print "('GRX[m2 s-1]=',2E12.4)",GRX(jminx),GRX(jmaxx)
!---neutral parameters
print "('LOG10 OX [m-3]     =',2F10.4)",LOG10(OX(jminx)),LOG10(OX(jmaxx))
print "('LOG10 HX [m-3]     =',2F10.4)",LOG10(HX(jminx)),LOG10(HX(jmaxx))
print "('LOG10 N2X [m-3]    =',2F10.4)",LOG10(N2X(jminx)),LOG10(N2X(jmaxx))
print "('LOG10 O2X [m-3]    =',2F10.4)",LOG10(O2X(jminx)),LOG10(O2X(jmaxx))
print "('LOG10 HEX [m-3]    =',2F10.4)",LOG10(HEX(jminx)),LOG10(HEX(jmaxx))
print "('LOG10 N4SX [m-3]   =',2F10.4)",LOG10(N4SX(jminx)),LOG10(N4SX(jmaxx))

print "('INNO =',I6)",INNO
IF ( INNO>=0 )  & !when CTIPe calculates NO
     & print "('LOG10 NNOX [m-3]   =',2F10.4)",LOG10(NNOX(jminx)),LOG10(NNOX(jmaxx))
print "('Tn [K]       =',2F10.4)",TNX(jminx),TNX(jmaxx)
print "('TINF [K]     =',2F10.4)",TINFX(jminx),TINFX(jmaxx)
print "('UNX [m s-1]  =',2F10.4)",UNX(jminx),UNX(jmaxx)

print "('DT [sec]     =',F10.4)",DT
print "('DTMIN [sec]  =',F10.4)",DTMIN
print "('F107D_dum    =',F10.4)",F107D_dum
print "('F107A_dum    =',F10.4)",F107A_dum
print "('SZA [deg]    =',2F10.4)",SZA_dum(jminx)*180./pi,SZA_dum(jmaxx)*180./pi
print "('FPAS         =',F10.4)",FPAS
print "('HPEQ         =',F10.4)",HPEQ
print "('HEPRAT       =',F10.4)",HEPRAT
print "('COLFACX      =',F10.4)",COLFACX
print "('IHEPLS       =',I6)",IHEPLS
print "('INPLS        =',I6)",INPLS

END IF !( sw_debug ) then

!dbg20120125:      midpoint = JMINX + (CTIPDIM-1)/2
      midpoint = IN + (IS-IN)/2
      IF ( lp>=1 .AND. lp<=6 )  midpoint = midpoint - 1
!nm20110909: calculating LT should be made FUNCTION!!!
      ltime = REAL(utime)/3600.0 + (plasma_grid_3d(midpoint,lp,mp,IGLON)*180.0/pi)/15.0
      IF ( ltime > 24.0 )  ltime = MOD(ltime, 24.0)

ret = gptlstart ('sw_output_fort167')
IF( sw_output_fort167.AND.mp==mpfort167.AND.lp==lpfort167 ) THEN
      WRITE(UNIT=LUN_FLIP1,FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)") mp,lp,LUN_FLIP1,REAL(UTIME)/3600., ltime
      WRITE(UNIT=LUN_FLIP2,FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)") mp,lp,LUN_FLIP2,REAL(UTIME)/3600., ltime
      WRITE(UNIT=LUN_FLIP3,FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)") mp,lp,LUN_FLIP3,REAL(UTIME)/3600., ltime
      WRITE(UNIT=LUN_FLIP4,FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)") mp,lp,LUN_FLIP4,REAL(UTIME)/3600., ltime

END IF !( sw_output_fort167
ret = gptlstop  ('sw_output_fort167')
if(sw_debug) print*,'sub-fl: UTs=',UTIME,' LThr=',ltime,' mp',mp,' lp',lp


!dbg20110131:
      IF ( sw_debug ) then
!sms$ignore begin
         WRITE(UNIT=PRUNIT,FMT="('mp=',i6,' lp=',i6,' UT=',F10.2,i7)") mp,lp,REAL(UTIME)/3600.,mype
!sms$ignore end
      END IF! ( sw_debug ) then
      ret = gptlstart ('flux_tube_solver_loop1')
      DO ipts=1,CTIPDIM
!N&T from the previous time step are absolute necesary for the solver...
!dbg20120501
        DO jth=1,ISPEC
          XIONNX(jth,ipts) = plasma_1d(jth,ipts)
        END DO !jth
!te
        TE_TIX(3,ipts) = plasma_1d(ISPEC+1,ipts)
!ti
        DO jth=1,2
          TE_TIX(jth,ipts) = plasma_1d(jth+ISPEC+1,ipts)
        END DO !jth
!vi
        DO jth=1,ISPEC
IF ( jth<=2 ) THEN
          XIONVX(jth,ipts) = plasma_1d(jth+ISPEC+3,ipts) 
ELSE
          XIONVX(jth,ipts) = zero
END IF
        END DO !jth


!### need to change these derived data type to a simpler arrays!!!
!nm20120412: need to restore the save V//(1:2) for parallel conv. effect
!nm20120412         XIONVX(1:ISPEC,ipts) = zero  !dbg20110927
!dbg20120501: add v// to perp transport
!dbg20120501         XIONVX(1,ipts)=plasma_1d(1,1,ipts) 
!dbg20120501         XIONVX(2,ipts)=plasma_1d(2,1,ipts) 
!dbg20120501         XIONVX(3:ISPEC,ipts) = zero



!20110930: inclusion for future neutral coupling:
! auroral electron heating-->EHTX(3)
! frictional heating for ions-->EHTX(1) 

         EHTX(  1:3        ,ipts)=zero  !dbg20110927


         IF ( INNO<0 ) THEN   !when flip calculates NO
           NNOX(              ipts)=zero !dbg20110927
         ELSE !when CTIPe calculates NO
           print *,'CTIPe calculates NO'
         END IF
         NHEAT(             ipts)=zero !dbg20110927
         hrate_cgs(1:22,    ipts)=zero !nm20121020
      END DO !ipts=
!nm20170328: correct setting for the He+/N+ switches
          if ( sw_ihepls==0 ) XIONNX(3,:) = zero
          if ( sw_inpls==0  ) XIONNX(4,:) = zero
      ret = gptlstop ('flux_tube_solver_loop1')

!tmp20151112: mlt for aurora in flip
      mlt = mlon_rad(mp)*180./pi/15.0D0-sunlons(1)*12.0D0/pi+12.0 !;[hr]

! call the flux tube solver (FLIP)
      ret = gptlstart ('CTIPINT')
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
     &              INNO, & !.. switch to turn on FLIP NO calculation if <0
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
     &            IHEPLS, & !.. switches He+ diffusive solution on if > 0
     &             INPLS, & !.. switches N+ diffusive solution on if > 0
     &              EHTX, & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
     &            TE_TIX, & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
     &     XIONNX,XIONVX, & !.. IN/OUT: 2D array, Storage for ion densities and velocities
     &             NHEAT, & !.. OUT: array, Neutral heating rate (eV/cm^3/s) 
     &             EFLAG, & !.. OUT: 2D array, Error Flags
     &                mp, &
     &                lp, &
     &                utime, & !dbg20141209
     &         hrate_cgs,mlt  ) !.. heating rates [eV/cm^3/s] !nm20121020
      ret = gptlstop  ('CTIPINT')






! output
      ret = gptlstart ('flux_tube_solver_loop2')
      DO ipts=1,CTIPDIM

         DO jth=1,ISPEC
            plasma_3d(ipts+IN-1,lp,mp,jth) = XIONNX(jth,ipts)
         END DO !jth


!dbg20160421 debug 
if( EFLAG(2,1)/=0 .and. 32<(ipts+IN-1) .and. (ipts+IN-1)<39 )then
!SMS$IGNORE begin
print"(i2,' XION=',e10.2,f7.0,f7.1,i4,' lp=',i3,' mp=',i2,' LT=',f6.1)",mype,plasma_3d(ipts+IN-1,lp,mp,1),(plasma_grid_Z(ipts+IN-1,lp)*1.e-3),((pi/2. - plasma_grid_GL(ipts+IN-1,lp))*180./pi),(ipts+IN-1),lp,mp,ltime
!SMS$IGNORE end
endif!(EFLAG




!te
         plasma_3d(ipts+IN-1,lp,mp,ISPEC+1) = TE_TIX(3,ipts)
!ti
         DO jth=1,2
            plasma_3d(ipts+IN-1,lp,mp,jth+ISPEC+1) = TE_TIX(jth,ipts)
         END DO !jth
!vi
         DO jth=1,ISPEV
           plasma_3d(ipts+IN-1,lp,mp,jth+ISPEC+3) = XIONVX(jth,ipts)
         END DO !jth


!nm20110404: save each component of heating rate for output
         IF ( sw_neutral_heating_flip==1 .AND. &
            &  MOD( (utime-start_time),ip_freq_output)==0) THEN
            if(parallelBuild) then
               print*,'sw_neutral_heating_flip=1 does not work in parallel'
               print*,'Stopping in Neut_heating'
               stop
            endif

!nm20121020            DO ij=1,MaxFluxTube
!nm20121020            j2d=ij+JMIN_IN(lp)-1
!nm20121020               DO jth=1,7
!nm20121020                  hrate_cgs_save(jth,ij)=hrate_cgs(jth,ij) !!(1) PGR neu_neu
!      hrate_cgs_save(2,j2d,mp_save)=hrate(2) !!PGR O1D
!      hrate_cgs_save(3,j2d,mp_save)=hrate(3) !!PGR ion_neu
!      hrate_cgs_save(4,j2d,mp_save)=hrate(4) !!PGR elec_ion
!      hrate_cgs_save(5,j2d,mp_save)=hrate(5) !!PGR SRO2dis
!      hrate_cgs_save(6,j2d,mp_save)=hrate(6) !!PGR UVN2dis
!      hrate_cgs_save(7,j2d,mp_save)=hrate(7) !!PGR 3bod CHECK dimension
!nm20121020               END DO !jth
!nm20121020            END DO !ij


!nm20121020            IF ( lp==NLP ) THEN
!20111118 output hrate_mks
!nm20121020               IF ( mp==1 ) write(5000,*) utime
               CALL get_neutral_heating_rate ( hrate_cgs , lp,mp )
!nm20121020: need to move the output to io_plasma_bin!!!
!nm20121020               DO jth=1,7
!nm20121020                  lun=5000+jth  
!nm20121020                  min_hrate =  huge(min_hrate)
!nm20121020                  max_hrate = -huge(max_hrate)
!nm20121020                  do lp0=1,NLP
!nm20121020                     min_hrate=min(min_hrate,MINVAL( &
!nm20121020                          & hrate_mks(jth,1:MaxFluxTube,lp0,mp)))
!nm20121020                     max_hrate=max(max_hrate,MAXVAL( &
!nm20121020                          & hrate_mks(jth,1:MaxFluxTube,lp0,mp)))
!nm20121020                  enddo
!nm20121020                  print *,'hrate=',lun,jth,mp,lp,min_hrate,max_hrate
!nm20121020                  write(lun,*) mp
!nm20121020                  do lp0=1,NLP
!nm20121020                     write(lun,*)hrate_mks(jth,1:MaxFluxTube,lp0,mp)
!nm20121020                  enddo
!nm20121020               END DO !jth=1,7
!nm20121020            END IF !( lp==NLP ) THEN

         END IF !( sw_neutral_heating_flip==1 ) THEN
!nm20121020:

      END DO       !DO ipts=1,CTIPDIM
      ret = gptlstop  ('flux_tube_solver_loop2')
!dbg20110927      NHEAT_cgs(IN:IS,lp,mp) =  NHEAT(        1:CTIPDIM) 

      ret = gptlstart ('WRITE_EFLAG')
      PRUNIT_dum = PRUNIT
      CALL WRITE_EFLAG(PRUNIT_dum, &  !.. Unit number to print results
     &                      EFLAG, &  !.. Error flag array
     &                         mp, &
     &                         lp,utime,ltime )
      ret = gptlstop  ('WRITE_EFLAG')


!nm20110822: no more allocatable arrays
!      IF ( ALLOCATED( ZX  ) )  DEALLOCATE(  ZX &
!     & ,STAT=stat_alloc         )
!      IF ( ALLOCATED( SLX  ) )  DEALLOCATE(  SLX ) !
!      IF ( ALLOCATED( GLX ) )  DEALLOCATE( GLX )
!      IF ( ALLOCATED( BMX  ) )  DEALLOCATE(  BMX ) !
!      IF ( ALLOCATED( GRX  ) )  DEALLOCATE(  GRX ) !
!      IF ( ALLOCATED( OX  ) )  DEALLOCATE(  OX ) !
!      IF ( ALLOCATED( HX  ) )  DEALLOCATE(  HX ) !
!      IF ( ALLOCATED( N2X  ) )  DEALLOCATE(  N2X ) !
!      IF ( ALLOCATED( O2X  ) )  DEALLOCATE(  O2X ) !
!      IF ( ALLOCATED( HEX  ) )  DEALLOCATE(  HEX ) !
!      IF ( ALLOCATED( N4SX  ) )  DEALLOCATE(  N4SX ) !
!      IF ( ALLOCATED( TNX  ) )  DEALLOCATE(  TNX ) !
!      IF ( ALLOCATED( TINFX  ) )  DEALLOCATE(  TINFX ) !
!      IF ( ALLOCATED( UNX ) )  DEALLOCATE( UNX )
!      IF ( ALLOCATED( SZA_dum  ) )  DEALLOCATE( SZA_dum  )
!      IF ( ALLOCATED( NNOX ) )  DEALLOCATE( NNOX )
!      IF ( ALLOCATED( EHTX ) )  DEALLOCATE( EHTX )
!      IF ( ALLOCATED( TE_TIX ) )  DEALLOCATE( TE_TIX )
!      IF ( ALLOCATED( XIONNX ) )  DEALLOCATE( XIONNX )
!      IF ( ALLOCATED( XIONVX ) )  DEALLOCATE( XIONVX )
!      IF ( ALLOCATED( NHEAT  ) )  DEALLOCATE( NHEAT  )
!if ( stat_alloc/=0 ) then
!  print *, ALLOCATED( ZX )
!  print *,"!STOP! ZX DEALLOCATION FAILD! in flux_tube_solver:",stat_alloc,mp,lp,in,is,jminx,jmaxx,ctipdim
!  STOP
!endif

      END SUBROUTINE flux_tube_solver
