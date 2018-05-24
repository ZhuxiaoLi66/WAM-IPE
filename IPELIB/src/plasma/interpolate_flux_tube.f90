!nm20130201: separated the Q interpolation part to module_Qinterpolate.f90
!dbg20120330: next version! should USE the apex routine to precisely estimate x(0),r(0), b0(0) for the imaginary flux tube...(one change at a time...)
!           :also density interpolation should be DOne with log???
!v17: previous versions were all WRONG!!! regarding how to implement ksi factor...
!v11: 20111101: ksi fac 1.0
!v15: 20111212: included Te/i transport
! DATE: 08 September, 2011
!
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
SUBROUTINE interpolate_flux_tube (mp,lp &
&, phi_t0,theta_t0 &
&, r0_apex &
&,mp_t0,lp_t0 ,&
& utime_local) !dbg20141209
  USE module_precision
!     plasma_grid_3d,plasma_grid_Z,plasma_grid_mag_colat,plasma_3d_old are all IN arrays
  USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,JMAX_IS,plasma_grid_3d,plasma_grid_Z,plasma_grid_mag_colat,ht90,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,plasma_3d,plasma_3d_old, mlon_rad,maxAltitude,minAltitude,minTheta,poleVal
  USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,sw_ksi,mype,lps,lpe,mps,mpe,nprocs,sw_ihepls,sw_inpls,sw_th_or_R
  USE module_IPE_dimension,ONLY: ISPEC,ISPET,IPDIM, ISTOT, NMP
  USE module_physical_constants,ONLY: earth_radius,pi,zero,rtd
  USE module_Qinterpolation,ONLY:Qinterpolation
  IMPLICIT NONE
!--- INPUT ---
  INTEGER (KIND=int_prec),INTENT(IN) :: mp
  INTEGER (KIND=int_prec),INTENT(IN) :: lp
  INTEGER (KIND=int_prec),INTENT(IN) :: utime_local !dbg20141209
  REAL(KIND=REAL_prec8),INTENT(IN) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
  REAL(KIND=REAL_prec8),INTENT(IN) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
  INTEGER (KIND=int_prec),DIMENSION(2,2), INTENT(IN) :: mp_t0,lp_t0 !1st dim:ihem, 2nd dim:ilp/imp
  REAL(KIND=REAL_prec8),INTENT(IN) :: r0_apex ![meter]
!---local
  INTEGER (KIND=int_prec) :: imp_max,ihem_max, ihem
  INTEGER (KIND=int_prec) :: ilp,lp0,imp,mp0
  INTEGER (KIND=int_prec) :: i,jth
  INTEGER (KIND=int_prec) :: i1d,midpoint ,kk


  REAL(KIND=REAL_prec8) :: factor, factor_ksi
!nm20140206      REAL(KIND=REAL_prec), DIMENSION(2) :: ksi_fac  !dim:ilp
  REAL(KIND=REAL_prec) :: ksi_fac
  REAL(KIND=REAL_prec8),DIMENSION(0:2) :: r
  REAL(KIND=REAL_prec8),DIMENSION(0:2) :: lambda_m,rapex,B0,x,y
!nm20170328            INTEGER (KIND=int_prec),PARAMETER :: TSP=3      !N(1:3) perp.transport
!nm20170328: N+ must be included in perp. transport
  INTEGER (KIND=int_prec),PARAMETER :: TSP=4      !N(1:4) perp.transport
  INTEGER (KIND=int_prec),PARAMETER :: iT=ISPEC+3   !add T(1:3)
!nm20130228      REAL(KIND=REAL_prec) :: n0(iT,2) !1d:species(N&T); 2dim:ilp

  INTEGER (KIND=int_prec),PARAMETER :: iB=iT+1 !add B
  INTEGER (KIND=int_prec),PARAMETER :: iR=iB+1 !add R
  REAL(KIND=REAL_prec8) :: Qint(iR,IPDIM,2,2)  !1d:species; ; 4d:imp; 5d:ilp
  REAL(KIND=REAL_prec8) :: Qint_dum(iR,IPDIM)  !1d:species;
  REAL(KIND=REAL_prec8),DIMENSION(ISTOT,IPDIM,2) :: plasma_2d !3d:imp
!dbg20140205 debug zonal transport
  REAL(KIND=REAL_prec8) :: mlon1,mlon2
  INTEGER (KIND=int_prec) :: mp1,mp2, mp3
!dbg20140528
  REAL(KIND=REAL_prec) :: tmp_fac1
  REAL(KIND=REAL_prec8) :: tmp_fac2
  REAL(KIND=REAL_prec) :: tmp_lam
!---



!array initialization: may not be necessary becaUSE they are local parameters...
  Qint(:,:,:,:)=zero

  IF ( sw_perp_transport==1 ) THEN !THETA only transport included
    imp_max=1
  ELSE IF ( sw_perp_transport>=2 ) THEN
    imp_max=2
  ENDIF

  !3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with dIFferent ExB drIFt
  IF ( sw_perp_transport==3 ) THEN
    ihem_max=2
  ELSE
    ihem_max=1
  ENDIF



  which_hemisphere: DO ihem=1,ihem_max

    IF ( lp_t0(ihem,1)>=1 ) THEN


      mp_t0_loop: DO imp=1,imp_max
        mp0 = mp_t0(ihem,imp)
        IF(nprocs == 1)THEN
          IF ( mp0<1   )  mp0 = mp0 + NMP
          IF ( mp0>NMP )  mp0 = mp0 - NMP
        ENDIF

!(1) Q interpolation: from Q_T0(mp_t0,lp_t0) --> Q_T1(mp,lp) for all 4 flux tubes
        lp_t0_loop: DO ilp=1,2 !outer/inner flux tubes

          ! (1.1) i0:  lp_t0(ihem,1)=l
          lp0 = lp_t0(ihem,ilp)
          CALL Qinterpolation (mp,lp &
          &, lp0, mp0 &
          &, iR, Qint_dum, TSP ) !nm20170328:

          Qint(1:iR,1:IPDIM,imp,ilp) = Qint_dum(1:iR,1:IPDIM)

        ENDDO lp_t0_loop !: DO ilp=1,2
      ENDDO mp_t0_loop!: DO imp=1,2




! (2) intepolate between the 4 flux tubes each point with the same Q value
! what KIND of interpolation is the most appropriate for the 4 points???
!note: 2 flux tubes only USEd IF imp_max=1/sw_perp_transport==1 ---THETA only transport included
!I should check quadratic interpolation CTIPe has using the third flux tube: tubes_quadratic_interpolate.f
!here in this interpolation should be DOne for both mp=j0,j1 simultaneously...

      rapex(0) = r0_apex
      IF ( sw_th_or_R==0 ) THEN
        lambda_m(0) = pi*0.50 - theta_t0(ihem)
      ELSE IF ( sw_th_or_R==1 ) THEN
        lambda_m(0) = ACOS(SQRT((earth_radius+ht90)/r0_apex  ))
      ENDIF !sw_th_or_R

!R of apex altitude for the two IN/OUT FTs
      DO ilp=1,2  !outer/inner flux tubes



        lp0 = lp_t0(ihem,ilp)
        midpoint = JMIN_IN(lp0) + ( JMAX_IS(lp0) - JMIN_IN(lp0) )/2
        rapex(ilp)=plasma_grid_Z(midpoint,lp0) + earth_radius



!note: this factor cannot work when lp<=6!!!
!nm20140630 i may not need this line at all???
        IF ( sw_th_or_R==1.and.lp>6 ) THEN
          lambda_m(ilp) = ACOS(SQRT((earth_radius+ht90)/rapex(ilp)))
!t PRINT "('!dbg20140627 lambda_m=',f8.6,' ilp=',i4,' rapex',f9.0)",lambda_m(ilp),ilp,rapex(ilp)
        ELSE
          !note: this factor cannot work when lp<=6!!!
          !becaUSE rapex DOes not mean anything for lp<=6
          lambda_m(ilp) = pi*0.5 - plasma_grid_mag_colat( JMIN_IN(lp0),lp0 )

        ENDIF !
      ENDDO  !DO ilp=1,2

!not sure which factor is more correct??? either r- or lambda (gip) base???
!note: tmp_fac1 is not available lp<=6 becaUSE rapex(1)=(2)

      IF ( sw_th_or_R==1.and.lp>6.and.rapex(1)/=rapex(2) ) THEN

        factor = ( rapex(0)-rapex(2) ) / ( rapex(1)-rapex(2) )

      ELSE IF ( lambda_m(1)/=lambda_m(2) ) THEN
        ! the values are only for NH
        factor = ( lambda_m(0) - lambda_m(2)) / (lambda_m(1) - lambda_m(2))
      ELSE
!SMS$ignore begin
        PRINT *,'!sub-interpolate_ft: STOP! INVALID factor:lambda', factor,lambda_m,mp,lp,mype
!SMS$ignore end
        STOP
      ENDIF





!error trap
      IF ( factor>1.0.OR.factor<0.0) THEN
!SMS$ignore begin
        PRINT"('sub-interpolate_flux_tube:INVALID factor=',e12.4,'mp=',i3,'lp=',i4,'rapex=',3e12.4,'lambda_m=',3e12.4,'mype=',i4)",factor,mp,lp,rapex(0:2),lambda_m(0:2),mype
!SMS$ignore end
        IF ( factor>1.0 )factor=1.0
        IF ( factor<0.0 )factor=0.0
      ENDIF



      mp_t0_loop1: DO imp=1,imp_max
        flux_tube_loopT1_fac: DO i=JMIN_IN(lp),JMAX_IS(lp) !9000
          i1d=i-JMIN_IN(lp)+1
          r(1:2)=Qint(iR,i1d,imp,1:2) !R


          r(0) = factor * ( r(1)-r(2) ) + r(2)


!weighting of X between Nin & Nout
! X can be either R or lambda (but only at IN/IS!!!)
          IF( r(1)/=r(2) ) THEN

            IF ( sw_th_or_R==1 ) THEN
              x(0:2) = r(0:2)
            else IF ( sw_th_or_R==0 ) THEN
              x(0:2) = lambda_m(0:2)
            ENDIF !sw_th_or_R

          ELSE IF( r(1)==r(2) ) THEN

            IF(sw_debug) THEN
!SMS$IGNORE BEGIN
              PRINT "('!R1=R2! r0',E14.6,' r1=',E14.6,' r2=',E14.6,' Qr1',E13.5,' Qr2=',E13.5,' i=',i6,' mp=',i3,' lp=',i4,' imp=',i3)",r(0:2),Qint(iR,i1d,imp,1:2),i,mp,lp,imp
              PRINT *, r(2),(earth_radius+ht90),(maxAltitude+earth_radius)
!SMS$IGNORE END
            ENDIF

            IF( r(2)==(earth_radius+ht90) ) THEN

!dbg20111006: somehow IS DOes not fit here???
              IF(sw_debug) THEN
!SMS$ignore begin
                PRINT *,'i=IN/S',i,i1d,' IN=',JMIN_IN(lp),' IS=',JMAX_IS(lp),lambda_m(0:2)*180./pi,mype
!SMS$ignore end
              ENDIF
              x(0:2) = lambda_m(0:2)

              IF(sw_debug) THEN
!SMS$IGNORE BEGIN
                PRINT *,'!dbg20120503: lambda',x(0:2)*180./pi
!SMS$IGNORE END
              ENDIF
            ELSE IF ( r(2)==(earth_radius+maxAltitude) ) THEN

              IF(sw_debug) THEN
!SMS$IGNORE BEGIN
                PRINT *,'midpoint lp<=6',i,i1d,' IN=',JMIN_IN(lp),' IS=',JMAX_IS(lp),lambda_m(0:2)*180./pi
!SMS$IGNORE END
              ENDIF
              x(0:2) = lambda_m(0:2)

            ELSE IF ( lp<7 ) THEN

              IF(sw_debug) THEN
!SMS$IGNORE BEGIN
                PRINT *,'lp<7',i,i1d,' IN=',JMIN_IN(lp),' IS=',JMAX_IS(lp),lambda_m(0:2)*180./pi
!SMS$IGNORE END
              ENDIF
              x(0:2) = lambda_m(0:2)

            ELSE
!dbg20140630: why becomes r1=r2 for lp>7???
!SMS$IGNORE BEGIN
              PRINT"('sub-Intrp:!STOP! INVALID R12!!!,i1d,imp,mp,mype',4i7)",i1d,imp,mp,mype
              PRINT *,'why lp>=7?',i,i1d,' IN=',JMIN_IN(lp),' IS=',JMAX_IS(lp),lambda_m(0:2)*180./pi
!SMS$IGNORE END
              x(0:2) = lambda_m(0:2)
!dbg20140630         STOP
            ENDIF
          ENDIF !   IF( r(1)/=r(2) ) THEN




! 1. interpolate Bfield intensity Bt0 at the imaginary FT(phi0,theta0) using dipole assumption
          IF(sw_debug) THEN
!SMS$ignore begin
            PRINT "('QintB=',3E12.4)",Qint(iB,i1d,imp,1:2),mype
!SMS$ignore end
          ENDIF
          B0(1:2)=Qint(iB,i1d,imp,1:2) * ( r(1:2)*r(1:2)*r(1:2) )/(r(0)*r(0)*r(0))
          B0(0)=( B0(1)+B0(2) )*0.50
          IF(sw_debug) THEN
!SMS$ignore begin
            PRINT "('B=',3E12.4)",B0(0),B0(1),B0(2),mype
!SMS$ignore end
          ENDIF
          !dbg20111101:v14
          !dbg20140627: need IF statement same as line 320!!!
          !       IF ( x(1)/=x(2) ) THEN
          B0(0) = ( (x(1)-x(0))*Qint(iB,i1d,imp,2) + (x(0)-x(2))*Qint(iB,i1d,imp,1) ) / ( x(1)-x(2) )
          IF(sw_debug) THEN
!SMS$ignore begin
            PRINT "('v14:B=',3E12.4)",B0(0) !,B0(1),B0(2),mype
!SMS$ignore end
          ENDIF

!dbg20141210: polar cap boundary: mlat(17)=71.69
!dbg20141210    IF ( sw_ksi==0 ) THEN
          IF ( lp<=17.or.sw_ksi==0 ) THEN
!nm20140206      ksi_fac(1:2) =1.000
            ksi_fac =1.000
          else IF ( sw_ksi==1 ) THEN
!nm20140206      ksi_fac(1:2) = B0(0) / Qint(iB,i1d,imp,1:2)
!SMS$ignore begin
            WRITE(6,*)'!STOP! INVALID option! sw_ksi=',sw_ksi
!SMS$ignore end
            STOP
!dbg20120330: new and CORRECT method to apply the ksi factor!
          else IF ( sw_ksi==2 ) THEN
!nm20140206      ksi_fac(1) = plasma_grid_3d(i,lp,mp,IBM) / B0(0)
            ksi_fac = plasma_grid_3d(i,lp,mp,IBM) / B0(0)
          ENDIF



!4. calculate N(phi0,theta0) with weighting of X between Nin & Nout
! X can be either R or lambda (but only at IN/IS)


          jth_loop4: DO jth=1,iT !TSP+3
            IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop4
!nm20170328: he+ ihepls<=0
            IF ( jth==3.and.sw_ihepls<=0 ) CYCLE jth_loop4
!nm20170328: n+ inpls<=0
            IF ( jth==4.and.sw_inpls<=0 ) CYCLE jth_loop4

            IF ( x(1)/=x(2) ) THEN
              IF(jth<=TSP)THEN      !for densities
!               plasma_2d(jth,i1d,imp) = (10**(( (x(1)-x(0))*DLOG10(Qint(jth,i1d,imp,2)) + (x(0)-x(2))*DLOG10(Qint(jth,i1d,imp,1)) ) / ( x(1)-x(2) )))*(ksi_fac*ksi_fac)
                plasma_2d(jth,i1d,imp) = ( (x(1)-x(0))*Qint(jth,i1d,imp,2) + (x(0)-x(2))*Qint(jth,i1d,imp,1) ) / ( x(1)-x(2) )*(ksi_fac*ksi_fac)
              else  !for temperatures
                plasma_2d(jth,i1d,imp) = ( (x(1)-x(0))*Qint(jth,i1d,imp,2) + (x(0)-x(2))*Qint(jth,i1d,imp,1) )*(ksi_fac**(4./3.)) / ( x(1)-x(2) )
              ENDIF

!error check
              IF (plasma_2d(jth,i1d,imp)<=zero) THEN

!SMS$ignore begin
                PRINT "(i3,i7,'subInt:!STOP! INVALID N/T: plasma_2d',E12.4,' Qint1',E12.4,' Qint2',E12.4,' jth',i2,' i1d=',i7,'i=',i7,' mp=',i3,' lp=',i3,' imp=',i2)" &
                &, mype,utime_local&
                &, plasma_2d(jth,i1d,imp)&
                &, Qint(jth,i1d,imp,1)&
                &, Qint(jth,i1d,imp,2)&
                &, jth,i1d,i,mp,lp,imp
                PRINT*,'check dIFX', (x(1)-x(0)), (x(0)-x(2)), ( x(1)-x(2) ),ksi_fac
                PRINT"('!check X!=',3E16.8)",x(1),x(0),x(2)
                PRINT"('!check B!=',3E12.4)",B0(1),B0(0),B0(2)
!SMS$ignore end
!               STOP
              ENDIF !IF (plasma_2d(jth,i1d,imp)<=zero) THEN

            ELSE !IF ( x(1)/=x(2) ) THEN
!SMS$ignore begin
              WRITE(6,*)'sub-Intrp:!STOP! INVALID x(1:2)',x(0:2),i1d,mp,lp
!SMS$ignore end
              STOP
            ENDIF !IF ( x(1)/=x(2) ) THEN

          ENDDO jth_loop4!jth=1,iT !=TSP+3

        ENDDO flux_tube_loopT1_fac !: DO i=in(lp),is(lp) !9000

      ENDDO mp_t0_loop1 !: DO imp=1,imp_max

!zonal interpolation
      flux_tube_loopT1_fac1: DO i=JMIN_IN(lp),JMAX_IS(lp)
        i1d=i-JMIN_IN(lp)+1

        jth_loop5: DO jth=1,iT
          IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop5
!nm20170328: he+ ihepls<=0
          IF ( jth==3.and.sw_ihepls<=0 ) CYCLE jth_loop5
!nm20170328: n+ inpls<=0
          IF ( jth==4.and.sw_inpls<=0 ) CYCLE jth_loop5
!---
          IF ( sw_perp_transport>=2 ) THEN

!nm20160419: replace with mlon1&2
            mlon1 = mlon_rad( mp_t0(ihem,1) )
            mlon2 = mlon_rad( mp_t0(ihem,2) )

            mp1 = mp_t0(ihem,1)
            mp2 = mp_t0(ihem,2)

!nm20160419: adjust mp0 for serial (for parallel, it is assumed mp0 is within the array bound)
            IF(nprocs==1)THEN
              IF ( mp1<1   )  mp1 = mp1 + NMP
              IF ( mp2<1   )  mp2 = mp2 + NMP
              IF ( mp2>NMP )  mp2 = mp2 - NMP
              IF ( mp1>NMP )  mp1 = mp1 - NMP
            ENDIF

            IF ( (mlon1-mlon2)/=zero) THEN

              IF (  mlon1 < mlon2 ) THEN

! B interpolation
                B0(0)             = ( (mlon1        - phi_t0(ihem) ) * plasma_grid_3d(i,lp,mp2,IBM)   &
                &                           + (  phi_t0(ihem) - mlon2        ) * plasma_grid_3d(i,lp,mp1,IBM)   &
                &                           ) / (mlon1 - mlon2)

                plasma_3d(i1d,lp,mp,jth) = ( (mlon1        - phi_t0(ihem) ) * plasma_2d(jth,i1d,2)   &
                &                           + (  phi_t0(ihem) - mlon2        ) * plasma_2d(jth,i1d,1)   &
                &                           ) / (mlon1 - mlon2)


                IF (jth==1.and. plasma_3d(i1d,lp,mp,jth)<=zero ) THEN
!SMS$IGNORE begin
                  PRINT *,mype,utime_local,"sub-interpFT:!STOP! INVALID plasma_3d L430",mlon1,phi_t0(ihem),plasma_2d(jth,i1d,2),mlon2,plasma_2d(jth,i1d,1),jth,i1d,ihem,lp,mp
                  PRINT *,(mlon1        - phi_t0(ihem) ),(phi_t0(ihem) - mlon2        ),(mlon1 - mlon2)
!SMS$IGNORE end
                  STOP
                ENDIF !jth==1


                ! calculate ksi_factor
                ksi_fac = plasma_grid_3d(i,lp,mp,IBM) / B0(0)  !is this correct???

                !apply ksi_factor to plasma_3d
                IF ( jth<=TSP ) THEN
                  factor_ksi = ksi_fac * ksi_fac
                ELSE !             IF ( jth>TSP ) THEN
                  factor_ksi = ksi_fac**(4./3.)
                ENDIF !             IF ( jth<=TSP ) THEN
                plasma_3d(i1d,lp,mp,jth) = plasma_3d(i1d,lp,mp,jth) * factor_ksi


                IF(jth==1.and.plasma_3d(i1d,lp,mp,jth)<=zero)THEN
!SMS$IGNORE begin
                  PRINT *,mype,utime_local,'sub-interpFT:!STOP! INVALID plasma_3d L451',mp,lp,plasma_3d(i1d,lp,mp,jth),factor_ksi,i1d,jth
!SMS$IGNORE end
                  STOP
                ENDIF !jth


              ELSE !    mlon1 >= mlon2
!
!SMS$ignore begin
                PRINT *, 'sub-interpFT:!STOP! INVALID mlon order! ihem=',ihem,' mp=',mp,' lp=',lp,'mp(1)',mp_t0(ihem,1),'mp(2)',mp_t0(ihem,2),mlon1*rtd,mlon2*rtd,mype
!SMS$ignore end
                STOP
              ENDIF
            ELSE    ! IF ( (mlon1-mlon2)==0.) THEN
!SMS$ignore begin
              PRINT *, 'sub-interpFT:!STOP! INVALID same mlon1&2!',ihem,mp,lp,mp_t0(ihem,1),mp_t0(ihem,2),mlon1*rtd,mlon2*rtd,mype
!SMS$ignore end
              STOP
              !       plasma_3d(jth,i1d) = plasma_2d(jth,i1d,1)
            ENDIF
          ELSE  !IF ( sw_perp_transport<2 ) THEN
            plasma_3d(i1d,lp,mp,jth) = plasma_2d(jth,i1d,1)
          ENDIF
!---



        ENDDO jth_loop5
      ENDDO flux_tube_loopT1_fac1


    ELSE IF ( lp_t0(ihem,1)==-9999 ) THEN !missing_value in module_find_nei...

!SMS$IGNORE begin
      IF(sw_debug)PRINT"('mype=',i3,'subInt:specialPole:mp=',i3,'lp=',i3,'mpt0=',2i3,'lpt0=',i4,i2)",mype,mp,lp,mp_t0(ihem,1:2),lp_t0(ihem,1:2)
!SMS$IGNORE end

!CAUTION! poleVal is calculated (in module_sub_plasma.f90) ONLY at mype=0 but is broadcast to all processors.
!nm20140630 temporary solution for pole
      flux_tube_loop: DO i=JMIN_IN(lp),JMAX_IS(lp)
        i1d=i-JMIN_IN(lp)+1

        jth_loop6: DO jth=1,iT
          IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop6
!nm20160420 special 3 point pole interpolation
          x(1)=zero
          x(0)=theta_t0(ihem)
          x(2)=minTheta

          DO imp=1,imp_max

!nm20180112 bug found!!!
            mp3 = mp_t0(ihem,imp)
            IF(nprocs==1)THEN
              IF ( mp3<1   ) THEN
                mp3 = mp3 + NMP
!SMS$IGNORE begin
                PRINT*,mype,'sub-interpFT:(1)corrected mp3!',mp3, mp_t0(ihem,imp),ihem,imp,lp
!SMS$IGNORE end
              ENDIF !mp3<
              IF ( mp3>NMP ) THEN
                mp3 = mp3 - NMP
!SMS$IGNORE begin
                PRINT*,mype,'sub-interpFT:(2)corrected mp3!',mp3, mp_t0(ihem,imp),ihem,imp,lp
!SMS$IGNORE end
              ENDIF !mp3>

            ENDIF !nprocs

            y(imp) = ( (x(1)-x(0))*plasma_3d_old(i,lp_t0(ihem,2),mp3,jth) + (x(0)-x(2))*poleVal(i,jth) ) / ( x(1)-x(2) )
          ENDDO !imp

          x(1)=mlon_rad( mp_t0(ihem,1) )
          x(0)=phi_t0(ihem)
          x(2)=mlon_rad( mp_t0(ihem,2) )
          plasma_3d(i1d,lp,mp,jth) = ( (x(1)-x(0))*y(2) + (x(0)-x(2))*y(1) ) / ( x(1)-x(2) )


          IF(jth==1.and.plasma_3d(i1d,lp,mp,jth)<=zero)THEN
!SMS$IGNORE begin
            PRINT*,mype,utime_local,'sub-interpFT:!STOP! INVALID plasma_3d L511',plasma_3d(i1d,lp,mp,jth),lp,mp,i,i1d
            PRINT*,mype,'sub-interpFT:mlon=',mlon_rad( mp_t0(ihem,1) ),phi_t0(ihem),mlon_rad( mp_t0(ihem,2) ),y,mp_t0(ihem,1),mp_t0(ihem,2),lp,mp,ihem
!SMS$IGNORE end
            STOP
          ENDIF !jth==1


        ENDDO jth_loop6
      ENDDO flux_tube_loop!: DO i=JMIN_IN(lp),JMAX_IS(lp)

    ELSE
!SMS$ignore begin
      PRINT *,'sub-Intrp:!STOP! INVALID lp_t0:',lp_t0,' mp',mp,' lp',lp,mype
!SMS$ignore end
      STOP
    ENDIF


  ENDDO which_hemisphere !: DO ihem=1,ihem_max


END SUBROUTINE interpolate_flux_tube
