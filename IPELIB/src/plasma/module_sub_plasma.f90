!dbg20120501: add v// to perp transport
!20110911: note: openmp was tried on jet but did not work: only thread 0 was USEd not the other thread...although other threads did exist...needs more investigation...
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
!sms$insert USE module_calcPoleVal  ,ONLY: calcPoleVal
  IMPLICIT NONE

!nm20121003module parameters are separated into module_plasma.f90

  PRIVATE
  PUBLIC :: plasma !dbg20120501 ,plasma_DATA_1d,plasma_DATA_1d4n

CONTAINS
!---------------------------
  SUBROUTINE plasma ( utime )
    USE module_input_parameters,ONLY:mpSTOP,ip_freq_output,start_time,STOP_time,&
    &     sw_neutral_heating_flip,sw_perp_transport,lpmin_perp_trans,lpmax_perp_trans,&
    &     sw_para_transport,sw_dbg_perp_trans,sw_exb_up,nprocs,mype,         &
    &     HPEQ_flip,barriersOn, solar_forcing_time_step, &
    &     perp_transport_time_step
    USE module_physical_constants,ONLY:rtd,zero
    USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_3d,plasma_grid_mag_colat,  &
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
    INTEGER :: utime_perp_transport, time_loop
    REAL :: t1, t2

!SMS$SERIAL BEGIN
    CALL CPU_TIME( t1 )
!SMS$SERIAL END
    

! save ut so that other SUBROUTINEs can refer to it
    utime_save=utime

    IF ( sw_neutral_heating_flip==1 )  THEN
      hrate_mks3d(:,:,:,:)=zero!0.0_REAL_prec
    ENDIF


!SMS$SERIAL (<plasma_3d,IN> :DEFAULT=IGNORE)BEGIN
    PRINT*, 'min/max plasma : ',MINVAL(plasma_3d), MAXVAL(plasma_3d), utime
!SMS$SERIAL END


    IF ( utime > 0 ) THEN

      utime_perp_transport = utime
      DO time_loop = 1, solar_forcing_time_step/perp_transport_time_step

        IF(nprocs==1) THEN ! Store special pole values for the serial CASE

          DO i=JMIN_IN(1),JMAX_IS(1)
            DO jth=1,ISTOT
              poleVal(i,jth) = SUM( plasma_3d(i,1,1:NMP,jth) ) / REAL(NMP)
            ENDDO
          ENDDO

        ELSE

!sms$insert call calcPoleVal

        ENDIF

!SMS$PARALLEL(dh, lp, mp) BEGIN

        plasma_3d_old = plasma_3d

!SMS$EXCHANGE(plasma_3d_old)

        DO mp = 1,mpSTOP
          DO lp = lpmin_perp_trans,lpmax_perp_trans

            CALL perpendicular_transport ( utime_perp_transport, mp, lp )

          ENDDO
        ENDDO
!SMS$PARALLEL END

        utime_perp_transport = utime_perp_transport + perp_transport_time_step

      ENDDO

    ENDIF


!SMS$PARALLEL(dh, lp, mp) BEGIN
    DO mp = 1,mpSTOP
      DO lp = lpmin_perp_trans,lpmax_perp_trans

        CALL flux_tube_solver ( utime,mp,lp )

      ENDDO
    ENDDO
!SMS$PARALLEL END

!SMS$SERIAL BEGIN
    CALL CPU_TIME( t2 )
    PRINT*, "PLASMA_TIME :", t2-t1
!SMS$SERIAL END

  END SUBROUTINE plasma

END MODULE module_sub_PLASMA
