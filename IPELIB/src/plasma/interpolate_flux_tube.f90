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
SUBROUTINE interpolate_flux_tube (mp,lp,phi_t0,theta_t0,r0_apex,mp_t0,lp_t0,utime_local)                              
  USE module_precision
  USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,JMAX_IS,plasma_grid_3d,plasma_grid_Z,plasma_grid_mag_colat, &
                                      ht90,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,plasma_3d,plasma_3d_old, mlon_rad, &
                                      maxAltitude,minAltitude,minTheta,poleVal
  USE module_input_parameters,ONLY:sw_perp_transport,mype,lps,lpe,mps,mpe,nprocs,sw_ihepls,sw_inpls,sw_convection_footpoint_0_or_apex_1
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
  REAL(KIND=REAL_prec) :: ksi_fac
  REAL(KIND=REAL_prec8),DIMENSION(0:2) :: r
  REAL(KIND=REAL_prec8),DIMENSION(0:2) :: lambda_m,r_apex,B0,x,y
  INTEGER (KIND=int_prec),PARAMETER :: TSP=4      !N(1:4) perp.transport
  INTEGER (KIND=int_prec),PARAMETER :: iT=ISPEC+3   !add T(1:3)

  INTEGER (KIND=int_prec),PARAMETER :: iB=iT+1 !add B
  INTEGER (KIND=int_prec),PARAMETER :: iR=iB+1 !add R
  REAL(KIND=REAL_prec8) :: Qint(iR,IPDIM,2,2)  !1d:species; ; 4d:imp; 5d:ilp
  REAL(KIND=REAL_prec8) :: Qint_dum(iR,IPDIM)  !1d:species;
  REAL(KIND=REAL_prec8),DIMENSION(ISTOT,IPDIM,2) :: plasma_2d !3d:imp
  REAL(KIND=REAL_prec8) :: mlon1,mlon2
  INTEGER (KIND=int_prec) :: mp1,mp2, mp3
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

      r_apex(0) = r0_apex
      IF ( sw_convection_footpoint_0_or_apex_1==0 ) THEN
        lambda_m(0) = pi*0.50 - theta_t0(ihem)
      ELSE IF ( sw_convection_footpoint_0_or_apex_1==1 ) THEN
        lambda_m(0) = ACOS(SQRT((earth_radius+ht90)/r0_apex  ))
      ENDIF !sw_convection_footpoint_0_or_apex_1

!R of apex altitude for the two IN/OUT FTs
      DO ilp=1,2  !outer/inner flux tubes



        lp0 = lp_t0(ihem,ilp)
        midpoint = JMIN_IN(lp0) + ( JMAX_IS(lp0) - JMIN_IN(lp0) )/2
        r_apex(ilp)=plasma_grid_Z(midpoint,lp0) + earth_radius



!note: this factor cannot work when lp<=6!!!
!nm20140630 i may not need this line at all???
        IF ( sw_convection_footpoint_0_or_apex_1==1.and.lp>6 ) THEN
          lambda_m(ilp) = ACOS(SQRT((earth_radius+ht90)/r_apex(ilp)))
!t PRINT "('!dbg20140627 lambda_m=',f8.6,' ilp=',i4,' r_apex',f9.0)",lambda_m(ilp),ilp,r_apex(ilp)
        ELSE
          !note: this factor cannot work when lp<=6!!!
          !becaUSE r_apex DOes not mean anything for lp<=6
          lambda_m(ilp) = pi*0.5 - plasma_grid_mag_colat( JMIN_IN(lp0),lp0 )

        ENDIF !
      ENDDO  !DO ilp=1,2

!not sure which factor is more correct??? either r- or lambda (gip) base???

      IF ( sw_convection_footpoint_0_or_apex_1==1.and.lp>6.and.r_apex(1)/=r_apex(2) ) THEN

        factor = ( r_apex(0)-r_apex(2) ) / ( r_apex(1)-r_apex(2) )

      ELSE IF ( lambda_m(1)/=lambda_m(2) ) THEN
        ! the values are only for NH
        factor = ( lambda_m(0) - lambda_m(2)) / (lambda_m(1) - lambda_m(2))
      ELSE
      ENDIF

        IF ( factor>1.0 )factor=1.0
        IF ( factor<0.0 )factor=0.0

      mp_t0_loop1: DO imp=1,imp_max
        flux_tube_loopT1_fac: DO i=JMIN_IN(lp),JMAX_IS(lp) !9000
          i1d=i-JMIN_IN(lp)+1
          r(1:2)=Qint(iR,i1d,imp,1:2) !R


          r(0) = factor * ( r(1)-r(2) ) + r(2)


!weighting of X between Nin & Nout
! X can be either R or lambda (but only at IN/IS!!!)
          IF( r(1)/=r(2) ) THEN

            IF ( sw_convection_footpoint_0_or_apex_1==1 ) THEN
              x(0:2) = r(0:2)
            else IF ( sw_convection_footpoint_0_or_apex_1==0 ) THEN
              x(0:2) = lambda_m(0:2)
            ENDIF !sw_convection_footpoint_0_or_apex_1

          ELSE IF( r(1)==r(2) ) THEN

            IF( r(2)==(earth_radius+ht90) ) THEN

              x(0:2) = lambda_m(0:2)

            ELSE IF ( r(2)==(earth_radius+maxAltitude) ) THEN

              x(0:2) = lambda_m(0:2)

            ELSE IF ( lp<7 ) THEN

              x(0:2) = lambda_m(0:2)

            ELSE
            ENDIF
          ENDIF !   IF( r(1)/=r(2) ) THEN




! 1. interpolate Bfield intensity Bt0 at the imaginary FT(phi0,theta0) using dipole assumption
          B0(1:2)=Qint(iB,i1d,imp,1:2) * ( r(1:2)*r(1:2)*r(1:2) )/(r(0)*r(0)*r(0))
          B0(0) = ( (x(1)-x(0))*Qint(iB,i1d,imp,2) + (x(0)-x(2))*Qint(iB,i1d,imp,1) ) / ( x(1)-x(2) )

!dbg20141210: polar cap boundary: mlat(17)=71.69
          IF ( lp<=17 ) THEN
            ksi_fac =1.000
          else
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
                plasma_2d(jth,i1d,imp) = ( (x(1)-x(0))*Qint(jth,i1d,imp,2) + (x(0)-x(2))*Qint(jth,i1d,imp,1) ) / ( x(1)-x(2) )*(ksi_fac*ksi_fac)
              else  !for temperatures
                plasma_2d(jth,i1d,imp) = ( (x(1)-x(0))*Qint(jth,i1d,imp,2) + (x(0)-x(2))*Qint(jth,i1d,imp,1) )*(ksi_fac**(4./3.)) / ( x(1)-x(2) )
              ENDIF

            ELSE !IF ( x(1)/=x(2) ) THEN
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

! B interpolation
                B0(0)             = ( (mlon1        - phi_t0(ihem) ) * plasma_grid_3d(i,lp,mp2,IBM)   &
                &                           + (  phi_t0(ihem) - mlon2        ) * plasma_grid_3d(i,lp,mp1,IBM)   &
                &                           ) / (mlon1 - mlon2)

                plasma_3d(i1d,lp,mp,jth) = ( (mlon1        - phi_t0(ihem) ) * plasma_2d(jth,i1d,2)   &
                &                           + (  phi_t0(ihem) - mlon2        ) * plasma_2d(jth,i1d,1)   &
                &                           ) / (mlon1 - mlon2)

                ! calculate ksi_factor
                ksi_fac = plasma_grid_3d(i,lp,mp,IBM) / B0(0)  !is this correct???

                !apply ksi_factor to plasma_3d
                IF ( jth<=TSP ) THEN
                  factor_ksi = ksi_fac * ksi_fac
                ELSE !             IF ( jth>TSP ) THEN
                  factor_ksi = ksi_fac**(4./3.)
                ENDIF !             IF ( jth<=TSP ) THEN
                plasma_3d(i1d,lp,mp,jth) = plasma_3d(i1d,lp,mp,jth) * factor_ksi
          ELSE  !IF ( sw_perp_transport<2 ) THEN
            plasma_3d(i1d,lp,mp,jth) = plasma_2d(jth,i1d,1)
          ENDIF
!---



        ENDDO jth_loop5
      ENDDO flux_tube_loopT1_fac1


    ELSE IF ( lp_t0(ihem,1)==-9999 ) THEN !missing_value in module_find_nei...

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

            mp3 = mp_t0(ihem,imp)
            IF(nprocs==1)THEN
              IF ( mp3<1   ) THEN
                mp3 = mp3 + NMP
              ENDIF !mp3<
              IF ( mp3>NMP ) THEN
                mp3 = mp3 - NMP
              ENDIF !mp3>

            ENDIF !nprocs

            y(imp) = ( (x(1)-x(0))*plasma_3d_old(i,lp_t0(ihem,2),mp3,jth) + (x(0)-x(2))*poleVal(i,jth) ) / ( x(1)-x(2) )
          ENDDO !imp

          x(1)=mlon_rad( mp_t0(ihem,1) )
          x(0)=phi_t0(ihem)
          x(2)=mlon_rad( mp_t0(ihem,2) )
          plasma_3d(i1d,lp,mp,jth) = ( (x(1)-x(0))*y(2) + (x(0)-x(2))*y(1) ) / ( x(1)-x(2) )

        ENDDO jth_loop6
      ENDDO flux_tube_loop!: DO i=JMIN_IN(lp),JMAX_IS(lp)

    ELSE
    ENDIF


  ENDDO which_hemisphere !: DO ihem=1,ihem_max


END SUBROUTINE interpolate_flux_tube
