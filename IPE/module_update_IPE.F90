!nm20140919: subroutine to update IPE
! DATE: 19 September, 2014
!********************************************
!***      Copyright 2014 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Whole Atmosphere Model(WAM)-Ionosphere Plasmasphere Electrodynamics(IPE) model Interface
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!-------------------------------------------- 
!
module module_update_IPE
   implicit none	

   private

   public :: update_IPE
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine update_IPE ( clock, rc )
    use ESMF
    use nnt_types_module
    use module_decomp
    USE module_precision
    USE module_input_parameters,ONLY: sw_perp_transport,utime,start_time,time_step, solar_forcing_time_step, &
                                      ip_freq_msis,sw_debug,nTimeStep,mype,ip_freq_eldyn,ip_freq_plasma, &
                                      swEsmfTime,internalTimeLoopMax, ip_freq_output, time_step
    USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d
!SMS$IGNORE BEGIN
    USE module_sub_eldyn,ONLY: eldyn
!SMS$IGNORE END
    USE module_NEUTRAL_MKS,ONLY: neutral 
    USE module_sub_PLASMA,ONLY: plasma
    USE module_output,ONLY: output
    implicit none
!---
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
!---
!    INTEGER(KIND=int_prec)           :: utime !universal time [sec]
!nm20160723: gptl
    integer :: ret
!---
!nm20161003: esmf timing lib
    real(ESMF_KIND_R8) :: beg_time, end_time
!---
!nm20161225: internal time loop
    integer :: utime_sub_loop, sub_time_loop !internal time loop counter
!---
       type(ESMF_Time) :: currTime
       type(ESMF_Time) :: stopTime
       integer(ESMF_KIND_I4) :: yy
       integer(ESMF_KIND_I4) :: mm
       integer(ESMF_KIND_I4) :: dd
       integer(ESMF_KIND_I4) :: h
       integer(ESMF_KIND_I4) :: m
       integer(ESMF_KIND_I4) :: s
       character(len=8) :: fmt ! format descriptor
       character(len=8) :: fmt1 ! format descriptor
       character(len=8) :: yy_str
       character(len=8) :: mm_str
       character(len=8) :: dd_str
       character(len=8) :: h_str
       character(len=8) :: m_str
       character(len=8) :: s_str
       character(len=12) :: timestamp_for_IPE_output_files
       character(len=16) :: timestamp_for_logfile


       if(swEsmfTime) CALL ESMF_VMWtime(beg_time)

       CALL ESMF_LogWrite("sub-update_IPE start:", ESMF_LOGMSG_INFO, rc=rc)


       CALL ESMF_LogWrite("sub-update_IPE neutral:", ESMF_LOGMSG_INFO, rc=rc)
       CALL neutral ( utime, start_time )

       utime_sub_loop = utime
       DO sub_time_loop = 1, time_step/solar_forcing_time_step

         CALL ESMF_LogWrite("sub-update_IPE eldyn:", ESMF_LOGMSG_INFO, rc=rc)  
         CALL eldyn ( utime_sub_loop )

         CALL ESMF_LogWrite("sub-update_IPE plasma:", ESMF_LOGMSG_INFO, rc=rc) 
         CALL plasma ( utime_sub_loop )

         utime_sub_loop = utime_sub_loop + solar_forcing_time_step

       ENDDO
       utime = utime + time_step

       IF( swEsmfTime )THEN  

         CALL ESMF_VMWtime(end_time)
         IF( mype==0 )write(unit=9999,FMT=*)mype,nTimeStep,"module_updateIPE RunTime=",(end_time-beg_time)

       ENDIF !swEsmfTime



       CALL ESMF_ClockGet(clock, currTime=currTime, rc=rc)
       CALL ESMF_TimeGet(currtime,yy=yy,mm=mm,dd=dd,h=h,m=m,s=s,rc=rc)                                 
       fmt = '(i2.2)'
       fmt1 = '(i4.4)'
       write (yy_str,fmt1) yy
       write (mm_str,fmt) MM
       write (dd_str,fmt) DD
       write (h_str,fmt) H
       write (m_str,fmt) M
       write (s_str,fmt) S
       CALL ESMF_ClockGet(clock, stopTime=stopTime, rc=rc)
       CALL ESMF_TimeGet(stopTime,yy=yy,mm=mm,dd=dd,h=h,m=m,s=s,rc=rc)                                 
       fmt = '(i2.2)'
       fmt1 = '(i4.4)'
       write (yy_str,fmt1) yy
       write (mm_str,fmt) MM
       write (dd_str,fmt) DD
       write (h_str,fmt) H
       write (m_str,fmt) M
       write (s_str,fmt) S
       timestamp_for_IPE_output_files = trim(yy_str)//trim(mm_str)//trim(dd_str)//trim(h_str)//trim(m_str)                  
       timestamp_for_logfile = trim(yy_str)//'-'//trim(mm_str)//'-'//trim(dd_str)//' '//trim(h_str)//':'//trim(m_str)
       if(mype==0) print *, "IPE progress : "//timestamp_for_logfile


        if(swEsmfTime) then
          CALL ESMF_VMWtime(end_time)
          if(mype==0.or.mype==8)write(unit=9999,FMT=*)mype,nTimeStep,"plasma endT=",(end_time-beg_time)
        end if

        if( MOD( utime, ip_freq_output) == 0) THEN
          CALL io_plasma_bin( 1, utime, timestamp_for_IPE_output_files )
        endif


        CALL ESMF_LogWrite("sub-update_IPE output:", ESMF_LOGMSG_INFO, rc=rc)  
        CALL output ( utime )

        nTimeStep = nTimeStep + 1

        CALL ESMF_LogWrite("sub-update_IPE finished:", ESMF_LOGMSG_INFO, rc=rc)

        rc = ESMF_SUCCESS

  end subroutine update_IPE

end module module_update_IPE
