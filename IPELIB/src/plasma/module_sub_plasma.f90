!dbg20120501: add v// to perp transport
!20110911: note: openmp was tried on jet but did not work: only thread 0 was used not the other thread...although other threads did exist...needs more investigation...
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
      MODULE module_sub_PLASMA
      USE module_precision
      USE module_IPE_dimension      ,ONLY: ISPEC,ISPET,ISPEV,IPDIM,NLP,NMP,ISTOT
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d,VEXBup,plasma_3d_old
!sms$insert use module_calcPoleVal  ,ONLY: calcPoleVal
      IMPLICIT NONE
      include "gptl.inc"

!nm20121003module parameters are separated into module_plasma.f90

      PRIVATE
      PUBLIC :: plasma !dbg20120501 ,plasma_data_1d,plasma_data_1d4n

      CONTAINS
!---------------------------
      SUBROUTINE plasma ( utime )
      USE module_input_parameters,ONLY:mpstop,ip_freq_output,start_time,stop_time,&
&     sw_neutral_heating_flip,sw_perp_transport,lpmin_perp_trans,lpmax_perp_trans,&
&     sw_para_transport,sw_debug,sw_dbg_perp_trans,sw_exb_up,nprocs,mype,         &
&     HPEQ_flip,ip_freq_paraTrans,barriersOn
      USE module_physical_constants,ONLY:rtd,zero
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_3d,plasma_grid_GL,  &
&     plasma_grid_Z,JMAX_IS,hrate_mks3d,poleVal
      USE module_PLASMA,ONLY:utime_save,plasma_1d
      USE module_perpendicular_transport,ONLY:perpendicular_transport
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
!--- local variables ---
      INTEGER (KIND=int_prec)  :: mp
      INTEGER (KIND=int_prec)  :: lp
      INTEGER (KIND=int_prec)  :: i,j,midpoint, i1d,k,ret
      INTEGER (KIND=int_prec)  :: jth,status
!d      INTEGER :: lun_dbg=999
!t      REAL(KIND=real_prec) :: phi_t0   !magnetic longitude,phi at T0
!t      REAL(KIND=real_prec) :: theta_t0 !magnetic latitude,theta at T0

      if(barriersOn) then
        ret = gptlstart ('sub_plasma_barrier')
!sms$insert      call ppp_barrier(status)
        ret = gptlstop  ('sub_plasma_barrier')
      endif

! save ut so that other subroutines can refer to it
      utime_save=utime

      ret = gptlstart ('apex_lon_loop') !24772.857

!SMS$PARALLEL(dh, lp, mp) BEGIN

!!SMS$IGNORE begin
      plasma_3d_old = plasma_3d
!!SMS$IGNORE end

      if(nprocs==1) then ! Store special pole values for the serial case
        ret = gptlstart ('poleVal')
        do i=JMIN_IN(1),JMAX_IS(1)
          do jth=1,ISTOT
            poleVal(i,jth) = SUM( plasma_3d(i,1,1:NMP,jth) ) / REAL(NMP)
          end do
        end do
        ret = gptlstop ('poleVal')
      else
!sms$insert call calcPoleVal
      endif

!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-1")
      if(barriersOn) then
        ret = gptlstart ('exchange_barrier')
!sms$insert      call ppp_barrier(status)
        ret = gptlstop ('exchange_barrier')
      endif
      ret = gptlstart ('EXCHANGE')
!SMS$EXCHANGE(plasma_3d_old)
      ret = gptlstop  ('EXCHANGE')
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-2")

        IF ( sw_neutral_heating_flip==1 )  THEN
           hrate_mks3d(:,:,:,:)=zero!0.0_real_prec
        ENDIF

!     apex_longitude_loop: DO mp = mpstrt,mpstop,mpstep !1,NMP
      apex_longitude_loop: DO mp = 1,mpstop
!nm20121115        mp_save=mp
        if ( sw_debug )  WRITE (0,"('sub-p: mp=',I4)")mp
!d        n0_2dbg(:)=zero
!       apex_latitude_height_loop: DO lp = lpstrt,lpstop,lpstep
        apex_latitude_height_loop: DO lp = 1,NLP

!20111025: not sure if these lines work when ut=0 & HPEQ=0.5(initial profiles are prepared within flip) , or maybe it is ok if they are zero?
!save the values from the previous time step...
          DO i=JMIN_IN(lp),JMAX_IS(lp)
             i1d=i-JMIN_IN(lp)+1
             DO jth=1,ISTOT
                plasma_1d(jth,i1d) = plasma_3d(i,lp,mp,jth)


                if(jth==1.and.plasma_3d(i,lp,mp,jth)<=zero)then
!SMS$IGNORE begin
                   print*,utime,'!STOP! INVALID plasma3d:mype',mype,plasma_3d(i,lp,mp,jth),i,lp,mp,jth
!SMS$IGNORE end
                   STOP
                endif !jth

             END DO !jth
          END DO !i
          if ( sw_debug )  WRITE (0,"('sub-p: lp=',I4)")lp
!nm20130401: transport is not called when HPEQ_flip=0.5 as initial profiles do not exist!
        IF ( HPEQ_flip==0.5 .AND. utime==0 ) THEN

           print *,lp,mp,'utime=',utime,' plasma perp transport is not called when HPEQ_flip=0.5 & start_time=0 because initial profiles do not exist!'

        ELSE IF ( utime>0 ) THEN

           ret = gptlstart ('perp_transport')
           IF ( sw_perp_transport>=1 ) THEN
              IF ( lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans ) THEN

                 CALL perpendicular_transport ( utime,mp,lp )


              ELSE  !IF ( lp>lpmin_perp_trans ) THEN

                 if(utime==start_time.AND.mp==1) then
                    midpoint=JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
                    print "('NO PERP. TRANS: mp=',I3,' lp=',I4,' mlatNd',F8.3,' apht',F10.2)", mp,lp,(90.-plasma_grid_GL(JMIN_IN(lp),lp)*rtd),plasma_grid_Z(midpoint,lp)*1.0e-3
                 endif
              END IF
           END IF !( sw_perp_transport>=1 ) THEN
           ret = gptlstop ('perp_transport')

        END IF !( HPEQ_flip==0.5 .AND. utime==0 ) THEN

! call flux tube solver
        ret = gptlstart ('flux_tube_solver')
        IF ( sw_para_transport==1 .and. MOD( (utime-start_time),ip_freq_paraTrans)==0 ) THEN 
           CALL flux_tube_solver ( utime,mp,lp )
        ELSE !IF ( sw_para_transport==0 ) THEN 

           DO i=JMIN_IN(lp),JMAX_IS(lp)
              i1d=i-JMIN_IN(lp)+1
              DO jth=1,ISTOT
                 plasma_3d(i,lp,mp,jth) = plasma_1d(jth,i1d)

                 if(jth==1.and.plasma_1d(jth,i1d)<=zero)then
!SMS$IGNORE begin
                    print*,utime,'sub-plasma:!STOP! INVALID plasma_1d af-paraT:mype',mype,plasma_1d(jth,i1d),lp,mp,i1d,i,jth
!SMS$IGNORE end
                    STOP
                 endif !jth

              END DO !jth
           END DO !i
        END IF !( sw_para_transport==1 ) THEN           
        ret = gptlstop ('flux_tube_solver')

! calculate neutral heating rate: NHEAT_mks in [eV kg-1 s-1]
!20110729: temporarily commented out to save time...
!dbg20110927          IF ( sw_neutral_heating_flip==1 ) &
!dbg20110927     &       CALL get_neutral_heating_rate ( )

! calculate the field line integrals for the electrodynamic solver
!t        CALL calculate_field_line_integrals ( )

        END DO apex_latitude_height_loop !: DO lp = 1
      END DO apex_longitude_loop !: DO mp = 

!SMS$PARALLEL END
      ret = gptlstop ('apex_lon_loop')
      if(barriersOn) then
        ret = gptlstart ('barrierAfterApexLongitudeLoop')
!sms$insert      call ppp_barrier(status)
        ret = gptlstop  ('barrierAfterApexLongitudeLoop')
      endif
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-4")
!!!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-5")
!!!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-6")
!!!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-7")

      END SUBROUTINE plasma

      END MODULE module_sub_PLASMA
