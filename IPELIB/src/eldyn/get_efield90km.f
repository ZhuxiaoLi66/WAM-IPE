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
      SUBROUTINE GET_EFIELD90km ( utime )
      USE module_precision
      USE efield,ONLY: nmlat,ylatm,dlonm,potent,ylonm,nmlon
      USE module_physical_constants,ONLY:pi,rtd,dtr,earth_radius,zero
      USE module_eldyn,ONLY:j0,j1,theta90_rad,ed1_90,ed2_90,coslam_m
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,JMIN_IN,JMAX_IS &
     &,mlon_rad,ht90,Be3,apexE,VEXBup,east,north,up, VEXBe,l_mag, VEXBth
      USE module_input_parameters,ONLY: sw_debug,NYEAR,NDAY,sw_exb_up    &
     &,sw_perp_transport,lpmin_perp_trans,lpmax_perp_trans,mype          &
     &,lps,start_time,ip_freq_output                                     &
     &,sw_th_or_R
      USE module_magfield,ONLY:sunlons
      USE module_sunloc,ONLY:sunloc
      IMPLICIT NONE
!
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime !universal time [sec]
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
      REAL    (KIND=real_prec) :: vperp !EXB at 90km only for th method


! array initialization
      Ed1_90=zero
      Ed2_90=zero

!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print"('utime=',i7,' (1)mype=',i2,' j0s=',i5,' lps=',i3)",utime   &
!     &,mype,j0(1,lps),lps
!!SMS$IGNORE end

      IF ( j0(1,lps)<0 ) THEN
! array initialization
        theta90_rad = zero
        coslam_m    = zero

!NOTE: ylatm-90: magnetic latitude [deg] --> ylatm:degrees from SP
! find out grid point at r0=90km(ht1) (theta90,phi90) of ylonm(ii),ylatm(jj) using dipole assumption
        mlat_loop130km0: DO j=0,nmlat

          theta130_rad = ( 90.0 - (ylatm(j)-90.0) ) * dtr
          CALL get_theta1_at_ht1(ht130,theta130_rad,ht90,theta90_rad(j))
          if(sw_debug) then
            print *,'mlat90[deg]=',j,'ori=',(ylatm(j)-90.),             &
     &               (90.-theta90_rad(j)*rtd)
          endif
        END DO mlat_loop130km0

!dbg20160408 sms debug: original parallel begin location
!SMS$PARALLEL(dh, lp, mp) BEGIN

! note that mlat90km is the constant in m-lon
        mlat_loop90km0: DO lp=1,NLP

! NH
!memo: mlat90_deg=(90.-plasma_grid_3d(IN,mp)%GL*rtd)
          IN = JMIN_IN(lp)
          coslam_m(1,lp)=COS(pi*0.5-plasma_grid_GL(IN,lp))

          if(coslam_m(1,lp)<=0.or.coslam_m(1,lp)>=1.)then
!SMS$IGNORE BEGIN
            print*,'sub-get_e:NH!STOP! INVALID coslam!',lp,IN,mype
            print*, coslam_m(1,lp),(90.-plasma_grid_GL(IN,lp)*rtd),     &
     &              plasma_grid_GL(IN,lp) 
!SMS$IGNORE end
            STOP
          end if

          IS = JMAX_IS(lp)
          coslam_m(2,lp) = COS( pi*0.5-plasma_grid_GL(IS,lp) )

          if(coslam_m(2,lp)<=0.or.coslam_m(2,lp)>=1.)then
!SMS$IGNORE BEGIN
            print*,'sub-get_e:SH!STOP! INVALID coslam!',lp,IS,mype
            print*, coslam_m(2,lp),(90.-plasma_grid_GL(IS,lp)*rtd),     & 
     &              plasma_grid_GL(IS,lp) 
!SMS$IGNORE end
            STOP
          end if

!EQ: potential difference is constant below 4.4988 < APEX=130km
          IF ( plasma_grid_GL(IN,lp)>theta90_rad(nmlat/2) ) THEN
            if(sw_debug) then
              print *,'check EQ',lp,(90.-plasma_grid_GL(IN,lp)*rtd)     &
     &                             ,(90.-theta90_rad(nmlat/2) *rtd)
            endif
            j0(1,lp)=nmlat/2+1   !1:NH
            j1(1,lp)=nmlat/2
            j0(2,lp)=nmlat/2     !2:SH
            j1(2,lp)=nmlat/2-1

!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print"('(2)mype=',i2,' j0=',i3,' lp=',i3)",mype,j0(1,lp),lp       &
!     &,j1(1,lp)
!!SMS$IGNORE end
            
          ELSE ! plasma_grid_3d(IN,mp)%GL>=4.49(R130 apex)

!           NH find the closest j of mlat90 from 130km
            mlat_loop130km1: DO j=nmlat,nmlat/2,-1 !NP(j=90)-->EQ(j=45)
!           d print *,'BEFORE',lp,j,theta90_rad(j),plasma_grid_3d(IN,mp)%GL
!           d     &             ,theta90_rad(j-1)


              IF (theta90_rad(j)<=plasma_grid_GL(IN,lp).AND.            &
     &            plasma_grid_GL(IN,lp)<=theta90_rad(j-1) ) THEN
!               dbg
                if(sw_debug) then
                  print *,'AFTER',lp,j,theta90_rad(j)                   &
     &                           ,plasma_grid_GL(IN,lp),theta90_rad(j-1)
                endif
                j0(1,lp)=j  !1:NH
                j1(1,lp)=j-1
                j0(2,lp)=nmlat-(j-1)  !2:SH
                j1(2,lp)=nmlat-j
!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print"('(3)mype=',i2,' j0=',i3,' lp=',2i3)",mype,j0(1,lp),lp      &
!     &,j1(1,lp)
!!SMS$IGNORE end

                EXIT mlat_loop130km1
              ELSE
                IF ( j==nmlat/2 ) THEN
                  print *,'sub-get_e:!STOP! could not find j!',lp,IN,j  &
     &                   ,plasma_grid_GL(IN,lp),theta90_rad(j)
                  STOP
                END IF
              END IF

            END DO mlat_loop130km1 !: DO j=0,nmlat        

!           dbg20110919
            if(sw_debug)then
              print *,lp,' NH',(90.-theta90_rad(j0(1,lp))*rtd)          &
     &                        ,(90.-plasma_grid_GL(IN,lp)*rtd)          &
     &                        ,(90.-theta90_rad(j1(1,lp))*rtd)
              print *,lp,' SH',(90.-theta90_rad(j0(2,lp))*rtd)          &
     &                        ,(90.-plasma_grid_GL(IS,lp)*rtd)          &
     &                        ,(90.-theta90_rad(j1(2,lp))*rtd)
            end if !(sw_debug)then
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
      utsecs=REAL(utime, real_prec)
      CALL sunloc(iyr,NDAY,utsecs) !iyr,iday,secs)        
      if (sw_debug) print *,'sunlons(1)',sunlons(1),' iyr',iyr,utsecs
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

!       mlon_rad: from 0 to 355.5 with dlon_ipe=4.5 deg resolution
!       mlon90_deg = mlon_rad(mp)*rtd
!       dbg REAL((mp-1),real_prec) * dlonm90km 
!       dbg
!       d   print *,'mp',mp,' mlon90_deg=',mlon90_deg(mp),' dlonm=',dlonm
!       find i of mlon 130km
!       dlonm: delon lon grid spacing in degree
!       BUG?  i = INT( (mlon90_deg/dlonm) , int_prec )
!


        mlon130_loop1: DO i=0,nmlon
          i0=i
          i1=i+1
          IF ( i1>nmlon ) i1=i+1-nmlon
          mlon130_0=mlon130_rad(i0)
          IF ( mlon130_rad(i0)>mlon130_rad(i1) ) then
            mlon130_0=mlon130_rad(i0)-pi*2.0  
          ENDIF
          if(sw_debug) then
            print *,mp,i0,'mlon130(i0)=',mlon130_0,' mlon(mp)='         &
     &        ,mlon_rad(mp),' mlon130(i1)=',mlon130_rad(i1)
          endif
          IF ( mlon_rad(mp)>=mlon130_0.AND.                             &
     &         mlon_rad(mp)<=mlon130_rad(i1) ) THEN
            EXIT mlon130_loop1
          ELSE
            if ( i==nmlon ) then
              print *,'sub-get_e:(2) !STOP! could not find mlon'        &
     &               ,mp,mlon_rad(mp),i0,mlon130_rad(i0)
              STOP
            end if
          END IF
        END DO mlon130_loop1 !: DO i=0,nmlon
!       dbg
        if(sw_debug)then
          print *,'i0',i0,'i1',i1
          print *,'mlon130_rad(i0)=',mlon130_rad(i0),                   &
     &           ' mlon130_rad(i1)=',mlon130_rad(i1)
        end if

        mlat_loop90km1: DO lp=1,NLP
          IN = JMIN_IN(lp)
          IS = JMAX_IS(lp)
          if(mp==1.and.lp>150.and.lp<158) then
             print"('mp',i3,' lp',i3,' IN=',i5,' latN',F6.2             &
     &            ,' IS=',i5,' latS',F6.2)",mp,lp                       &
     &            ,IN,(90.-plasma_grid_GL(IN,lp)*rtd)                   &
     &            ,IS,(90.-plasma_grid_GL(IS,lp)*rtd)
          endif
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



!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print                                                             &
!     &"('(4)mype=',i2,' j0=',i3,' lp=',i3,' mp=',i2,' i0=',i3,'i1=',i3)"&
!     &,mype,j0(1,lp),lp,mp,i0,i1
!      print*,"(43)mype=",mype,j1(1,lp)
!!SMS$IGNORE end



          jj0=j0(1,lp)              !1:NH
          jj1=j1(1,lp)
          pot_i1=( potent(i1,jj0)+potent(i1,jj1) )*0.50 
          pot_i0=( potent(i0,jj0)+potent(i0,jj1) )*0.50
          if (r<=0..or.coslam_m(1,lp)==0..or.d_phi_m==0.)then
!SMS$IGNORE BEGIN
            print*,'sub-get_e:NH!STOP! INVALID',lp,ihem,mp,mype
            print*, r,coslam_m(1,lp),d_phi_m
!SMS$IGNORE END
            STOP
          endif
!!dbg20160408 sms new decomp debug
!      if (lp==1.and.mp==1) then
!!SMS$IGNORE BEGIN
!      print *,'(44)mype=',mype,mp,lp,pot_i1,pot_i0                      &
!     &,potent(i1,jj1),potent(i0,jj1),jj1                                &
!     &,coslam_m(1,lp)
!!SMS$IGNORE end
!      end if                    !mp==41
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
          if (r<=0..or.coslam_m(2,lp)==0..or.d_phi_m==0.)then
!SMS$IGNORE BEGIN
            print*,'sub-get_e:SH!STOP! INVALID',lp,ihem,mp,mype
            print*, r,coslam_m(2,lp),d_phi_m
!SMS$IGNORE END
            STOP
          endif
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
          if (d_lam_m==0.)then
            print *,'sub-get_ed2:NH!STOP! INVALID',lp,d_lam_m
            STOP
          endif
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
      if(sw_debug.AND.mp==1) print *,lp,mp,'!nm20130830: ed1_90NH='     &
     &,  ed1_90(1,lp,mp),' SH=',ed1_90(2,lp,mp) 

            v_e(1) =   Ed2_90(1,lp,mp) / Be3(lp,mp) !(4.18) +mag-east(d1?) 
            v_e(2) = - Ed1_90(1,lp,mp) / Be3(lp,mp) !(4.19) +down/equatorward(d2?)
!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print"('(5)mype=',i2,' mp=',i2,' lp=',i3)",mype,mp,lp
!      print *,mype,'ve2[m/s]=',v_e(2),'ed1[mV/m]=',Ed1_90(1,lp,mp)*1.E3 &
!     &,' Be3[tesla]=',Be3(lp,mp) 
!!SMS$IGNORE end

! EXB in geographic frame at (midpoint,lp,mp)
            vexbgeo(east )=(v_e(1)*apexE(midpoint,lp,mp,east,1))        &
     &                    +(v_e(2)*apexE(midpoint,lp,mp,east,2))
            vexbgeo(north)=(v_e(1)*apexE(midpoint,lp,mp,north,1))       &
     &                    +(v_e(2)*apexE(midpoint,lp,mp,north,2))
            vexbgeo(up   )=(v_e(1)*apexE(midpoint,lp,mp,up   ,1))       &
     &                    +(v_e(2)*apexE(midpoint,lp,mp,up   ,2))
            if(sw_debug) then
              print *,'sub-StR:',v_e(1),apexE(midpoint,lp,mp,up,1)      &
     &                          ,v_e(2),apexE(midpoint,lp,mp,up,2)
            endif
! EXB  magnetic exact upward at APEX (midpoint) 
            VEXBup(lp,mp)=                                              &
!nm20130201: the new, correct definition
     &       vexbgeo(east)  * l_mag(midpoint,lp,mp,east ,2)             &
     &     + vexbgeo(north) * l_mag(midpoint,lp,mp,north,2)             &
     &     + vexbgeo(up   ) * l_mag(midpoint,lp,mp,up   ,2)
!vexbgeo(up) !nm20140528
!     &                    (v_e(1)*apexE(midpoint,lp,mp,up,1))           &
!     &                   +(v_e(2)*apexE(midpoint,lp,mp,up,2))

!nm20140528: VEXBth at 90km required for th method, sinI_m for NH
            if ( sw_th_or_R==0 ) then !th method(ctipe/shawn)
               ipts = JMIN_IN(lp) ! if ihem=1 NH
!               ipts = JMAX_IS(lp) ! if ihem=2 SH
      if(sw_debug)                                                      &
     & write(9001,*)'!dbg20140616 sw_th_or_R=',sw_th_or_R,ipts,mp,lp
! EXB in geographic frame at 90km NH (IN,lp,mp)
            vexbgeo(east )=(v_e(1)*apexE(ipts,lp,mp,east,1))            &
     &                    +(v_e(2)*apexE(ipts,lp,mp,east,2))
            vexbgeo(north)=(v_e(1)*apexE(ipts,lp,mp,north,1))           &
     &                    +(v_e(2)*apexE(ipts,lp,mp,north,2))
            vexbgeo(up   )=(v_e(1)*apexE(ipts,lp,mp,up   ,1))           &
     &                    +(v_e(2)*apexE(ipts,lp,mp,up   ,2))
! vperp at 90km NH +poleward
!dbg20150319: astrid confirmed this!!!
               vperp =                                                  &
     &  v_e(2) * (-1.)  !converted from equator to poleward
!d     &       vexbgeo(east)  * l_mag(ipts,lp,mp,east ,2)                 &
!d     &     + vexbgeo(north) * l_mag(ipts,lp,mp,north,2)                 &
!d     &     + vexbgeo(up   ) * l_mag(ipts,lp,mp,up   ,2)

!nm20140616: kensyo: compare VEXBup vs. vperp
!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print"('(6)mype=',i2,' mp=',i2,' lp=',i3)",mype,mp,lp
!      write(6,"(i2,'vperp=',3f9.2)")mype,vperp,v_e(2)
!!SMS$IGNORE end

!VEXBth: horizontal component, positive EQUATORward: is calculated only for th-method
!dbg20150319: it might be better to estimate sinI at JMIN_IN
               VEXBth(lp,mp)= vperp * (-1.) * sinI_m(1) !if ihem=1 NH
!dbg20160408 sms
!!SMS$IGNORE BEGIN
!      print"('(7)mype=',i2,' mp=',i2,' lp=',i3)",mype,mp,lp
!      write(6,"(i2,' VEXBth=',3f10.3)")mype,VEXBth(lp,mp),vperp         &
!     &,sinI_m(1)!if ihem=1 NH
!!SMS$IGNORE end

!nm20150516(1):SED test
!      VEXBth(lp,mp)=0.0

            else if ( sw_th_or_R==1 ) then !R method(gip)
               ipts = midpoint
            end if !( sw_th_or_R==1
! exact magnetic eastward  
! R method: at APEX
! th method: at 90km JMIN_IN
!dbg20150319: astrid confirmed this!!!
            VEXBe(lp,mp)=                                               &
     &  v_e(1)  !+eastward

!d vexbgeo(east )*l_mag(ipts,lp,mp,east ,1)       &
!d     &                  +vexbgeo(north)*l_mag(ipts,lp,mp,north,1)       &
!d     &                  +vexbgeo(up  ) *l_mag(ipts,lp,mp,up   ,1)


!nm20150516(2):SED test
!      VEXBe(lp,mp)=0.0


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


!dbg20140228
!      if ( sw_debeg_eld ) then 
         IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN
!SMS$SERIAL(<vexbup,vexbe,vexbth,IN>:default=ignore) BEGIN
            write(unit=2011,FMT='(20E12.4)') vexbup
            write(unit=2012,FMT='(20E12.4)') vexbe
            write(unit=2013,FMT='(E12.4)  ') sunlons(1)
            write(unit=2014,FMT='(20E12.4)') vexbth
!SMS$SERIAL END
         END IF                 ! ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN
!      endif                     !( sw_debug_eld ) then 

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
