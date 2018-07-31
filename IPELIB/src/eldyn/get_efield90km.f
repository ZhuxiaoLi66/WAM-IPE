!FUNC-linearinterpolation does not work! need more debugging. temporary use the average of the two potentials.
!another idea is to 2Dbilinear interpolation of potential onto the ipe grid, and then one can do the usual central differencing.
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
!note: potential in efield.f is the value at 130km
!     & h_r  = 130.0e3,    ! reference height [m] (same as for apex.F90)     
      SUBROUTINE GET_EFIELD90km ( utime_local )
      USE module_precision
      USE efield_ipe,ONLY: nmlat,ylatm,dlonm,potent,ylonm,nmlon
      USE module_physical_constants,ONLY:pi,rtd,dtr,earth_radius,zero
      USE module_eldyn,ONLY:j0,j1,theta90_rad,ed1_90,ed2_90,coslam_m
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_FIELD_LINE_GRID_MKS
      USE module_input_parameters
      USE module_magfield,ONLY:sunlons
      USE module_sunloc,ONLY:sunloc
      IMPLICIT NONE
!
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime_local !universal time [sec]
!
      REAL(KIND=real_prec) :: utsecs
      INTEGER :: iyr
      REAL(KIND=real_prec) :: theta130_rad
      REAL(KIND=real_prec),PARAMETER :: ht130 = 130.0E+03     !height in meter
      INTEGER (KIND=int_prec) :: j,i,mp,lp
      INTEGER (KIND=int_prec) :: IN,IS
!      REAL(KIND=real_prec) :: potent_i0,potent_i1
!      REAL(KIND=real_prec) :: LINEAR_INTERPOLATION !the intepolated value YY at (XX)
      REAL(KIND=real_prec) :: mlon90_deg !deg
      REAL(KIND=real_prec) :: d_phi_m, d_lam_m !in radian
      REAL(KIND=real_prec) :: r !in meter
      INTEGER (KIND=int_prec) :: jj0,jj1
      INTEGER(KIND=int_prec) :: i0,i1,ihem
      REAL(KIND=real_prec) :: pot_i0,pot_i1, pot_j0,pot_j1
      REAL(KIND=real_prec),DIMENSION(0:nmlon) :: mlon130_rad
      REAL(KIND=real_prec) :: mlon130_0
      REAL(KIND=real_prec) :: cos2Lambda_m,sinLambda_m(2),sinI_m(2)
      INTEGER (KIND=int_prec ) :: midpoint, ipts
      REAL    (KIND=real_prec) :: v_e(2)   !1:ed2/be3 (4.18) ;2: -ed1/be3 (4.19)
      REAL    (KIND=real_prec) :: vexbgeo(east:up) !EXB in gegraphic frame
!      REAL    (KIND=real_prec) :: vperp !EXB at 90km only for th method


! array initialization
      Ed1_90=zero
      Ed2_90=zero


      IF ( j0(1,lps)<0 ) THEN
! array initialization
        theta90_rad = zero
        coslam_m    = zero

!NOTE: ylatm-90: magnetic latitude [deg] --> ylatm:degrees from SP
! find out grid point at r0=90km(ht1) (theta90,phi90) of ylonm(ii),ylatm(jj) using dipole assumption
        mlat_loop130km0: DO j=0,nmlat

          theta130_rad = ( 90.0 - (ylatm(j)-90.0) ) * dtr
          CALL get_theta1_at_ht1(ht130,theta130_rad,ht90,theta90_rad(j))
        END DO mlat_loop130km0

!dbg20160408 sms debug: original parallel begin location
!SMS$PARALLEL(dh, lp, mp) BEGIN

! note that mlat90km is the constant in m-lon
        mlat_loop90km0: DO lp=1,NLP

! NH
!memo: mlat90_deg=(90.-plasma_grid_3d(IN,mp)%GL*rtd)
          IN = JMIN_IN(lp)
          coslam_m(1,lp)=COS(pi*0.5-plasma_grid_mag_colat(IN,lp))

          IS = JMAX_IS(lp)
          coslam_m(2,lp) = COS( pi*0.5-plasma_grid_mag_colat(IS,lp) )

!EQ: potential difference is constant below 4.4988 < APEX=130km
          IF ( plasma_grid_mag_colat(IN,lp)>theta90_rad(nmlat/2) ) THEN
            j0(1,lp)=nmlat/2+1   !1:NH
            j1(1,lp)=nmlat/2
            j0(2,lp)=nmlat/2     !2:SH
            j1(2,lp)=nmlat/2-1

          ELSE ! plasma_grid_3d(IN,mp)%GL>=4.49(R130 apex)

!           NH find the closest j of mlat90 from 130km
            mlat_loop130km1: DO j=nmlat,nmlat/2,-1 !NP(j=90)-->EQ(j=45)

              IF (theta90_rad(j)<=plasma_grid_mag_colat(IN,lp).AND.            &
     &            plasma_grid_mag_colat(IN,lp)<=theta90_rad(j-1) ) THEN
!               dbg
                j0(1,lp)=j  !1:NH
                j1(1,lp)=j-1
                j0(2,lp)=nmlat-(j-1)  !2:SH
                j1(2,lp)=nmlat-j

                EXIT mlat_loop130km1
              ELSE
              END IF

            END DO mlat_loop130km1 !: DO j=0,nmlat        

          END IF                   ! ( plasma_grid_3d(IN,mp)%GL>theta90_rad(nmlat/2) ) THEN

        END DO mlat_loop90km0!: DO lp=1,NLP
!dbg20160408 sms debug
!SMS$PARALLEL END

      END IF                    !( j0(1,lps)>0 ) THEN

      d_phi_m = dlonm * dtr !constant
      r = earth_radius + ht90 ![m]


!     prepare sunlons(1) before converting from MLT(ylonm) to mlon
!     note: NYEAR=2000: is above the MAX recommended for extrapolation!!!
      iyr=1999 
!     iday=97
      utsecs=REAL(utime_local, real_prec)
      CALL sunloc(iyr,NDAY,utsecs) !iyr,iday,secs)        
!     convert from MLT(ylonm)[rad] to mlon[deg]
      mlon130_loop0: DO i=0,nmlon
        mlon130_rad(i)=(ylonm(i)-180.)*pi/180.+sunlons(1)
!       make sure that 0 <=mlon130< 2*pi
        IF( mlon130_rad(i)< 0.0   ) mlon130_rad(i)=mlon130_rad(i)+pi*2.0
        IF( mlon130_rad(i)>=pi*2.0) mlon130_rad(i)=mlon130_rad(i)-pi*2.0
      END DO mlon130_loop0 !: DO i=0,nmlon


!dbg20160408 sms debug: 
!SMS$PARALLEL(dh, lp, mp) BEGIN
      mlon_loop90km0: DO mp=1,NMP

        mlon130_loop1: DO i=0,nmlon
          i0=i
          i1=i+1
          IF ( i1>nmlon ) i1=i+1-nmlon
          mlon130_0=mlon130_rad(i0)
          IF ( mlon130_rad(i0)>mlon130_rad(i1) ) then
            mlon130_0=mlon130_rad(i0)-pi*2.0  
          ENDIF
          IF ( mlon_rad(mp)>=mlon130_0.AND.                             &
     &         mlon_rad(mp)<=mlon130_rad(i1) ) THEN
            EXIT mlon130_loop1
          ELSE
          END IF
        END DO mlon130_loop1 !: DO i=0,nmlon

        mlat_loop90km1: DO lp=1,NLP
          IN = JMIN_IN(lp)
          IS = JMAX_IS(lp)

!         computipng ed1_90(lp,mp)
!         FUNC-linearinterpolation does not work! need more debugging. temporary use the average of the two potentials.
!         linear interpolation of the potent at plasma_grid_3d(IN,mp) in mlat
!         NH
!         d          potent_i0 = LINEAR_INTERPOLATION(theta90_rad(j0(1,lp))
!         d     &   ,potent(i0(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!         d     &   ,potent(i0(mp),j1(1,lp)),plasma_grid_3d(IN,mp)%GL)
!         d          potent_i1 = LINEAR_INTERPOLATION(theta90_rad(j0(1,lp))
!         d     &   ,potent(i1(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!         d     &   ,potent(i1(mp),j0(1,lp)),plasma_grid_3d(IN,mp)%GL)
!         ii0=i0
!         ii1=i1






          jj0=j0(1,lp)              !1:NH
          jj1=j1(1,lp)
          pot_i1=( potent(i1,jj0)+potent(i1,jj1) )*0.50 
          pot_i0=( potent(i0,jj0)+potent(i0,jj1) )*0.50

          ed1_90(1,lp,mp)=-1.0/r/coslam_m(1,lp)                         &
     &                         *(pot_i1-pot_i0)/d_phi_m
!         SH
!         d          potent_i0=LINEAR_INTERPOLATION(theta90_rad(j0(1,lp)) 
!         d     &   ,potent(i0(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!         d     &   ,potent(i0(mp),j1(1,lp)),plasma_grid_3d(IN,mp)%GL)
!         d          potent_i1=LINEAR_INTERPOLATION(theta90_rad(j0(1,lp))
!         d     &   ,potent(i1(mp),j0(1,lp)),theta90_rad(j1(1,lp))
!         d     &   ,potent(i1(mp),j0(1,lp)),plasma_grid_3d(IN,mp)%GL)
!         ii0=i0
!         ii1=i1
          jj0=j0(2,lp)              !2:SH
          jj1=j1(2,lp)              !2:SH
          pot_i1=( potent(i1,jj0)+potent(i1,jj1) )*0.50
          pot_i0=( potent(i0,jj0)+potent(i0,jj1) )*0.50
          ed1_90(2,lp,mp)=-1.0/r/coslam_m(2,lp)*(pot_i1-pot_i0)/d_phi_m
!         computing ed2_90(lp,mp) continues
!         calculate sinIm !eq(3.7)
          cos2Lambda_m = coslam_m(2,lp) * coslam_m(2,lp)      ! 0<cos2<1
          sinLambda_m(1)  = + SQRT( 1.0 - cos2Lambda_m )  !>0 ---NH 
          sinLambda_m(2)  = - SQRT( 1.0 - cos2Lambda_m )  !<0 ---SH 
          sinI_m(1:2)= 2.0*sinLambda_m(1:2)/SQRT(4.0-3.0*cos2Lambda_m)
!         NH
          ihem=1
          jj0=j0(ihem,lp)              !1:NH
          jj1=j1(ihem,lp)              !1:NH
          d_lam_m = theta90_rad( jj1 ) - theta90_rad( jj0 )
          pot_j1=( potent(i0,jj1)+potent(i1,jj1) )*0.50
          pot_j0=( potent(i0,jj0)+potent(i1,jj0) )*0.50
          ed2_90(1,lp,mp)=+1.0/r/sinI_m(ihem)*(pot_j1-pot_j0)/d_lam_m   &
     &*(-1.) !dbg20140224
!         dbg20111108     &*(-1.)*sinI_m(ihem)     !E_m_lambda (5.10)
!         SH
          ihem=2
          jj0=j0(ihem,lp)              !2:SH
          jj1=j1(ihem,lp)              !2:SH
          d_lam_m = theta90_rad( jj1 ) - theta90_rad( jj0 )
          pot_j1=( potent(i0,jj1)+potent(i1,jj1) )*0.50
          pot_j0=( potent(i0,jj0)+potent(i1,jj0) )*0.50
          ed2_90(2,lp,mp)=+1.0/r/sinI_m(ihem)*(pot_j1-pot_j0) /d_lam_m  &
     &*(-1.) !dbg20140224
!         dbg20111108     &*(-1.)*sinI_m(ihem)  !E_m_lambda (5.10)

          IF(sw_exb_up<=1.and.sw_perp_transport>=1) then
             if(lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans) THEN 

!initialization
                VEXBup(lp,mp)=zero
                VEXBe(lp,mp)=zero
                VEXBth(lp,mp)=zero

!           (0) self consistent electrodynamics comming soon...
!           (1) WACCM E empirical model
!           Ed1/2[V/m] at ( phi_t1(mp), theta_t1(lp) ), Be3[T]
!           note: Ed1_90, Ed2_90, and Be3 are constant along magnetic field lines!!! 
            midpoint = JMIN_IN(lp) + (JMAX_IS(lp) - JMIN_IN(lp))/2

!nm20130830: Ed1/2_90 should be constant along magnetic field lines!!!
            v_e(1) =   Ed2_90(1,lp,mp) / Be3(lp,mp) !(4.18) +mag-east(d1?) 
            v_e(2) = - Ed1_90(1,lp,mp) / Be3(lp,mp) !(4.19) +down/equatorward(d2?)

! EXB in geographic frame at (midpoint,lp,mp)
            vexbgeo(east )=(v_e(1)*apexE(midpoint,lp,mp,east,1))        &
     &                    +(v_e(2)*apexE(midpoint,lp,mp,east,2))
            vexbgeo(north)=(v_e(1)*apexE(midpoint,lp,mp,north,1))       &
     &                    +(v_e(2)*apexE(midpoint,lp,mp,north,2))
            vexbgeo(up   )=(v_e(1)*apexE(midpoint,lp,mp,up   ,1))       &
     &                    +(v_e(2)*apexE(midpoint,lp,mp,up   ,2))
! EXB  magnetic exact upward at APEX (midpoint) 
!           VEXBup(lp,mp) = v_e(2) * (-1.0) !convert from down to UPward
            VEXBup(lp,mp) = vexbgeo(up)

               ipts = JMIN_IN(lp) ! if ihem=1 NH
!               ipts = JMAX_IS(lp) ! if ihem=2 SH
! EXB in geographic frame at 90km NH (IN,lp,mp)
            vexbgeo(east )=(v_e(1)*apexE(ipts,lp,mp,east,1))            &
     &                    +(v_e(2)*apexE(ipts,lp,mp,east,2))
            vexbgeo(north)=(v_e(1)*apexE(ipts,lp,mp,north,1))           &
     &                    +(v_e(2)*apexE(ipts,lp,mp,north,2))
            vexbgeo(up   )=(v_e(1)*apexE(ipts,lp,mp,up   ,1))           &
     &                    +(v_e(2)*apexE(ipts,lp,mp,up   ,2))
! vperp at 90km NH +poleward
!VEXBth: horizontal component, positive EQUATORward: is calculated only for th-method
!dbg20150319: it might be better to estimate sinI at JMIN_IN
!              VEXBth(lp,mp) = v_e(2) * sinI_m(1) !if ihem=1 NH
               VEXBth(lp,mp) = v_e(2) / sinI_m(1) !if ihem=1 NH

!temporary solution to test the code
               VEXBe(lp,mp) = 0.0


              if ( sw_perp_transport==1 )   VEXBe(lp,mp)=zero

            else                   ! (lp<lpmin_perp_trans.or.lp>lpmax_perp_trans) THEN 
               VEXBup(lp,mp)=zero
               VEXBe(lp,mp)=zero
               VEXBth(lp,mp)=zero 
            end if              !(lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans) THEN 
          END IF !( sw_exb_up<=1.and. ... ) 


        END DO mlat_loop90km1 !: DO lp=1,NLP
      END DO mlon_loop90km0     !: DO mp=1,nmp
!dbg20160408 sms debug 
!SMS$PARALLEL END


         IF ( MOD( (utime_local-start_time),ip_freq_output)==0 ) THEN
!SMS$SERIAL(<vexbup,vexbe,vexbth,IN>:default=ignore) BEGIN
            write(unit=2011,FMT='(20E12.4)') vexbup
            write(unit=2012,FMT='(20E12.4)') vexbe
            write(unit=2013,FMT='(E12.4)  ') sunlons(1)
            write(unit=2014,FMT='(20E12.4)') vexbth
!SMS$SERIAL END
         END IF                 ! ( MOD( (utime_local-start_time),ip_freq_output)==0 ) THEN

      END SUBROUTINE GET_EFIELD90km
!
!20110919: not used for now,so commented out. needs more debugging
!     FUNCTION LINEAR_INTERPOLATION(X0,Y0,X1,Y1,XX)
!     USE module_precision
!     REAL(KIND=real_prec),INTENT(IN) :: X0,Y0  !(X0,Y0)
!     REAL(KIND=real_prec),INTENT(IN) :: X1,Y1  !(X1,Y1)
!     REAL(KIND=real_prec),INTENT(IN) :: XX     !X coordinate of the intepolated value YY
!     REAL(KIND=real_prec) :: LINEAR_INTERPOLATION !the intepolated value YY at (XX)
!     LINEAR_INTERPOLATION = Y0 + (XX - X0) * (Y1 - Y0) / (X1 - X0) 
!     END  FUNCTION LINEAR_INTERPOLATION

      SUBROUTINE get_theta1_at_ht1 (ht0,theta0, ht1,theta1)

      USE module_precision
      USE module_physical_constants,ONLY:earth_radius,pi
      IMPLICIT NONE
      REAL(KIND=real_prec),INTENT(IN) :: ht0 !m / km
      REAL(KIND=real_prec),INTENT(IN) :: theta0 ![rad]
      REAL(KIND=real_prec),INTENT(IN) :: ht1 !m / km
      REAL(KIND=real_prec),INTENT(OUT) :: theta1
!     ---local
      REAL(KIND=real_prec) :: sintheta0
      REAL(KIND=real_prec) :: sintheta1
!     find out grid point at ht1=90km (theta1,phi1) of theta0 using dipole assumption
      sintheta0 = SIN( theta0 )
      sintheta1 = sintheta0*SQRT((earth_radius+ht1)/(earth_radius+ht0))
      theta1    = ASIN(sintheta1)
!     SH:example: mlat=30 comlat=90-30=60, mlatSH=-30,comlatSH=90+(90-60)=180-60
      IF ( theta0 > pi*0.50 ) theta1 = pi-theta1

      END SUBROUTINE get_theta1_at_ht1
