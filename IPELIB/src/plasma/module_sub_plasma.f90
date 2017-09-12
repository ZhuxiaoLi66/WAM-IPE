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
      USE module_IPE_dimension,ONLY: ISPEC,ISPET,ISPEV,IPDIM,NLP,NMP,ISTOT
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d,VEXBup,plasma_3d_old
      IMPLICIT NONE
      include "gptl.inc"

      PRIVATE
      PUBLIC :: plasma 

      CONTAINS
!---------------------------
      SUBROUTINE plasma ( utime, timestamp_for_IPE )
      USE module_input_parameters,ONLY:mpstop,ip_freq_output,start_time,stop_time,&
&     sw_neutral_heating_flip,sw_perp_transport,lpmin_perp_trans,lpmax_perp_trans,sw_para_transport,sw_debug,        &
&     sw_dbg_perp_trans,sw_exb_up,parallelBuild,mype, &
& HPEQ_flip, ut_start_perp_trans, dumpFrequency
      USE module_physical_constants,ONLY:rtd,zero
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_3d,plasma_grid_GL,plasma_grid_Z,JMAX_IS,hrate_mks3d
      USE module_PLASMA,ONLY:utime_save,plasma_1d
      USE module_perpendicular_transport,ONLY:perpendicular_transport
      IMPLICIT NONE
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
!--- local variables ---
      INTEGER (KIND=int_prec) :: mp
      INTEGER (KIND=int_prec) :: lp
      INTEGER (KIND=int_prec) :: i,j,midpoint, i1d,k,ret  
      INTEGER (KIND=int_prec) :: jth  
      integer :: status

      character(len=12), INTENT(IN) :: timestamp_for_IPE


!---------------
! array initialization
!dbg20120313 note! needed to comment out temporarily for exb reading test:sw_exb=5
if ( sw_exb_up/=5 ) then
!JFM moved to allocate_arrays
!     VEXBup=zero
end if


! save ut so that other subroutines can refer to it
      utime_save=utime
      print *, 'GHGM IN PLASMA, UTIME : ', utime

      ret = gptlstart ('apex_lon_loop') !24772.857
!SMS$PARALLEL(dh, lp, mp) BEGIN
      plasma_3d_old = plasma_3d
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-1")
      ret = gptlstart ('exchange_barrier')
!sms$insert      call ppp_barrier(status)
      ret = gptlstop ('exchange_barrier')
      ret = gptlstart ('EXCHANGE')
!SMS$EXCHANGE(plasma_3d_old)
      ret = gptlstop  ('EXCHANGE')
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-2")
!     apex_longitude_loop: DO mp = mpstrt,mpstop,mpstep !1,NMP
      apex_longitude_loop: DO mp = 1,mpstop
!nm20121115        mp_save=mp
        IF ( sw_neutral_heating_flip==1 )  hrate_mks3d(:,:,mp,:)=zero
        if ( sw_debug )  WRITE (0,"('sub-p: mp=',I4)")mp
!d        n0_2dbg(:)=zero

!dbg20120412: sw_divvpar4t
!        DO i=1,NPTS2D
!          plasma_3d4n(i,mp)%V_ms1(1:ISPEV)=zero
!        END DO


!!!dbg20120125: only temporary used to switch on the transport only during the daytime...
!dbg20120509        IF ( sw_rw_sw_perp_trans.AND.sw_perp_transport(mp)==0 )  CALL activate_perp_transport (utime,mp)
!!!dbg20120125:
!       apex_latitude_height_loop: DO lp = lpstrt,lpstop,lpstep
        apex_latitude_height_loop: DO lp = 1,NLP
!nm20121115          lp_save=lp


!dbg20120228: debug how2validate the transport
if(sw_dbg_perp_trans.and.utime==start_time.and.lp==1)then
  if(parallelBuild) then
    print*,'sw_dbg_perp_trans=true does not work for parallel runs'
    print*,'Stopping module_sub_plasma'
    stop
  endif
!  DO j=1,NLP
!    DO i=JMIN_IN(j),JMAX_IS(j)
!      DO jth=1,ISTOT
!JFM     plasma_3d(jth,i,lp,mp)=100.0
!dbg20120501      plasma_3d(mp,j)%N_m3( 1:ISPEC,i)=100.0
!dbg20120501      plasma_3d(mp,j)%Te_k(         i)=100.0
!dbg20120501      plasma_3d(mp,j)%Ti_k( 1:ISPET,i)=100.0
!      END DO !jth
!    END DO !i
!  END DO !j
end if
!if(sw_dbg_perp_trans) print *, '1!dbg max o+',MAXVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) ),MINVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) )




!20111025: not sure if these lines work when ut=0 & HPEQ=0.5(initial profiles are prepared within flip) , or maybe it is ok if they are zero?
!save the values from the previous time step...
!     print*,'JFM,mp,lp,JMIN_IN,JMAX_IS',mp,lp,JMIN_IN(lp),JMAX_IS(lp)!1,  1,    1, 1115
                                                                      !1,  2, 1118, 2232
                                                                      !1,  3, 2235, 3349
                                                                      !1,170,44430,44438
          DO i=JMIN_IN(lp),JMAX_IS(lp)
             i1d=i-JMIN_IN(lp)+1
             DO jth=1,ISTOT
                plasma_1d(jth,i1d) = plasma_3d(i,lp,mp,jth)
!dbg20120501          DO i=1,IPDIM
!dbg20120501            n0_1d%N_m3( 1:ISPEC,i) = plasma_3d(mp,lp)%N_m3( 1:ISPEC,i)
!dbg20120501            n0_1d%Te_k(         i) = plasma_3d(mp,lp)%Te_k(         i)
!dbg20120501            n0_1d%Ti_k( 1:ISPET,i) = plasma_3d(mp,lp)%Ti_k( 1:ISPET,i)
             END DO !jth
          END DO !i
!dbg20120501
!dbg20120501          j=1
!dbg20120501          DO i=JMIN_IN(lp),JMAX_IS(lp)
!dbg20120501            i1d=i-JMIN_IN(lp)+1
!dbg20120501            DO k=1,2
!dbg20120501              plasma_1d(k,j,i1d)=plasma_3d4n(i,mp)%V_ms1(k)
!dbg20120501            END DO
!dbg20120501          END DO

!if(sw_dbg_perp_trans) print *, '2!dbg max o+',MAXVAL( n0_1d%N_m3( 1,1:IPDIM) ),MINVAL( n0_1d%N_m3( 1,1:IPDIM) )

          if ( sw_debug )  WRITE (0,"('sub-p: lp=',I4)")lp





!dbg20120509          IF ( sw_perp_transport(mp)>=1 ) THEN
!nm20130401: transport is not called when HPEQ_flip=0.5 as initial profiles do
!not exist!
        IF ( HPEQ_flip==0.5 .AND. utime==ut_start_perp_trans ) THEN

          print *,lp,mp,'utime=',utime,' plasma perp transport is not called when HPEQ_flip=0.5 & start_time=0 because initial profiles do not exist!'

        ELSE IF ( utime>0 ) THEN

          ret = gptlstart ('perp_transport')
          IF ( sw_perp_transport>=1 ) THEN
            IF ( lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans ) THEN
              CALL perpendicular_transport ( utime,mp,lp )




            ELSE  !IF ( lp>lpmin_perp_trans ) THEN
if(utime==start_time) then
midpoint=JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
print "('NO PERP. TRANS: mp=',I3,' lp=',I4,' mlatNd',F8.3,' apht',F10.2)", mp,lp,(90.-plasma_grid_GL(JMIN_IN(lp),lp)*rtd) &
& ,plasma_grid_Z(midpoint,lp)*1.0e-3
endif
            END IF
          END IF !( sw_perp_transport>=1 ) THEN
          ret = gptlstop ('perp_transport')

        END IF !( HPEQ_flip==0.5 .AND. utime==0 ) THEN

! update the boundary conditions if the top of the flux tube is open
!t        CALL update_flux_tube_boundary_condition ( )

!!nm20110822: moved to flux_tube_plasma  
!! calculate Solar Zenith Angle [radians]
!!          CALL Get_SZA ( utime,mp,lp ) 

! call flux tube solver
          ret = gptlstart ('flux_tube_solver')
          IF ( sw_para_transport==1 ) THEN 
            CALL flux_tube_solver ( utime,mp,lp )
          ELSE IF ( sw_para_transport==0 ) THEN 

!dbg20111101:v9: temporary ...
            DO i=JMIN_IN(lp),JMAX_IS(lp)
             i1d=i-JMIN_IN(lp)+1
              DO jth=1,ISTOT
                 plasma_3d(i,lp,mp,jth) = plasma_1d(jth,i1d)
!dbg20120501            DO i=1,IPDIM
!dbg20120501              plasma_3d(mp,lp)%N_m3( 1:ISPEC,i)=            n0_1d%N_m3( 1:ISPEC,i)
!dbg20120501              plasma_3d(mp,lp)%Te_k(         i)=            n0_1d%Te_k(         i)
!dbg20120501              plasma_3d(mp,lp)%Ti_k( 1:ISPET,i)=            n0_1d%Ti_k( 1:ISPET,i)
              END DO !jth
           END DO !i

!dbg20120501
!dbg20120501            j=1
!dbg20120501            DO i=JMIN_IN(lp),JMAX_IS(lp)
!dbg20120501            i1d=i-JMIN_IN(lp)+1
!dbg20120501            DO k=1,2
!dbg20120501              plasma_3d4n(i,mp)%V_ms1(k)=                   plasma_1d(k,j,i1d)
!dbg20120501            END DO
!dbg20120501          END DO            


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
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-4")


!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-5")

!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-6")

      
!sms$compare_var(plasma_3d,"module_sub_plasma.f90 - plasma_3d-7")

      END SUBROUTINE plasma
      END MODULE module_sub_PLASMA
