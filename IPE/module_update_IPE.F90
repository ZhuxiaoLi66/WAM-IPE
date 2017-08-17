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
    USE module_input_parameters,ONLY: sw_perp_transport,utime,start_time,time_step,ip_freq_msis,sw_debug,nTimeStep,mype,ip_freq_eldyn,ip_freq_plasma,swEsmfTime,internalTimeLoopMax
    USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d
!SMS$IGNORE BEGIN
    USE module_sub_eldyn,ONLY: eldyn
!SMS$IGNORE END
    USE module_NEUTRAL_MKS,ONLY: neutral 
    USE module_sub_PLASMA,ONLY: plasma
    USE module_output,ONLY: output
    implicit none
!g    include "gptl.inc"
!---SMS
      integer PPP__GlobalSize(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__GlobalUpper(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__GlobalLower(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__Permute(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__GlobalStart(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__GlobalStop(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__HaloLower(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__HaloUpper(PPP_MAX_RANK,PPP_MAX_VARS)
      integer PPP__DecompType(PPP_MAX_VARS)
      integer PPP__DataType(PPP_MAX_VARS)
      logical IAM_ROOT
      integer GlobalLower1(PPP_MAX_RANK,PPP_MAX_EXCHANGE_VARS)
      integer GlobalUpper2(PPP_MAX_RANK,PPP_MAX_EXCHANGE_VARS)
      integer Permute_Decomp3(PPP_MAX_RANK,PPP_MAX_EXCHANGE_VARS)
      integer GlobalStart4(PPP_MAX_RANK,PPP_MAX_EXCHANGE_VARS)
      integer GlobalStop5(PPP_MAX_RANK,PPP_MAX_EXCHANGE_VARS)
      character*32 VarName6
      logical SMS_DEBUGGING_ON
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
    integer :: ncnt !internal time loop counter
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
       character(len=13) :: timestamp_for_IPE_output_files


    call ESMF_LogWrite("sub-update_IPE start:", ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_LogWrite("sub-update_IPE ClockPrint:", ESMF_LOGMSG_INFO, rc=rc)
    call ESMF_ClockPrint(clock, options="currTime name startTime stopTime timeStep", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!nm20160315---copied from driver_ipe_sms.f90


!nm20161225 internal time loop
      internalTimeLoop: DO ncnt=1,internalTimeLoopMax

! GHGM moved from the bottom of the loop - we need to update the time first...


      if (IAM_ROOT()) then
        print"('Utime=',2i7,f11.4,' nTimeStep',2i2)",utime,(MOD(utime,86400)),(MOD(utime,86400)/3600.),nTimeStep,ncnt
        if(swEsmfTime) write(unit=9999,FMT=*)'update_IPE',nTimeStep,':ut=',utime
      endif
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-5")

      if (SMS_DEBUGGING_ON()) then
      GlobalLower1(:,1:1) = 1
      GlobalUpper2(:,1:1) = 1
      Permute_Decomp3(:,1:1) = 0
      GlobalStart4(:,1:1) = 1
      GlobalStop5(:,1:1) = 1

!-SMS          storing information about variable: plasma_3d

      GlobalLower1(1,1) = LBOUND(plasma_3d,1)
      GlobalLower1(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalLower1(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalLower1(4,1) = LBOUND(plasma_3d,4)
      GlobalUpper2(1,1) = UBOUND(plasma_3d,1)
      GlobalUpper2(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalUpper2(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalUpper2(4,1) = UBOUND(plasma_3d,4)
      Permute_Decomp3(2,1) = 1
      Permute_Decomp3(3,1) = 2
      GlobalStart4(1,1) = LBOUND(plasma_3d,1)
      GlobalStart4(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalStart4(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalStart4(4,1) = LBOUND(plasma_3d,4)
      GlobalStop5(1,1) = UBOUND(plasma_3d,1)
      GlobalStop5(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalStop5(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalStop5(4,1) = UBOUND(plasma_3d,4)

      VarName6= 'plasma_3d'
      ppp__str= 'driver_ipe.f90 - plasma_3d-5'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('module_update_IPE.F90:L109.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()


!---empirical E-field
!g      ret = gptlstart ('eldyn')
      IF ( sw_perp_transport>=1.AND. MOD( (utime-start_time),ip_freq_eldyn)==0 ) THEN
        call ESMF_LogWrite("sub-update_IPE eldyn:", ESMF_LOGMSG_INFO, rc=rc)  

        if(sw_debug.and.IAM_ROOT())print*,'call eldyn: utime=',utime,nTimeStep,ncnt
        CALL eldyn ( utime )
      ENDIF
!g      ret = gptlstop  ('eldyn')
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-6")

      if (SMS_DEBUGGING_ON()) then
      GlobalLower1(:,1:1) = 1
      GlobalUpper2(:,1:1) = 1
      Permute_Decomp3(:,1:1) = 0
      GlobalStart4(:,1:1) = 1
      GlobalStop5(:,1:1) = 1

!-SMS          storing information about variable: plasma_3d

      GlobalLower1(1,1) = LBOUND(plasma_3d,1)
      GlobalLower1(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalLower1(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalLower1(4,1) = LBOUND(plasma_3d,4)
      GlobalUpper2(1,1) = UBOUND(plasma_3d,1)
      GlobalUpper2(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalUpper2(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalUpper2(4,1) = UBOUND(plasma_3d,4)
      Permute_Decomp3(2,1) = 1
      Permute_Decomp3(3,1) = 2
      GlobalStart4(1,1) = LBOUND(plasma_3d,1)
      GlobalStart4(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalStart4(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalStart4(4,1) = LBOUND(plasma_3d,4)
      GlobalStop5(1,1) = UBOUND(plasma_3d,1)
      GlobalStop5(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalStop5(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalStop5(4,1) = UBOUND(plasma_3d,4)

      VarName6= 'plasma_3d'
      ppp__str= 'driver_ipe.f90 - plasma_3d-6'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('module_update_IPE.F90:L159.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()

!--- update neutral
        IF ( MOD( (utime-start_time),ip_freq_msis)==0 ) THEN 

          if(sw_debug.and.IAM_ROOT())print*,ncnt,'call neutral: utime=',utime,nTimeStep

!g          ret = gptlstart ('neutral')
          call ESMF_LogWrite("sub-update_IPE neutral:", ESMF_LOGMSG_INFO, rc=rc)

!nm20161003 esmf timing library
          if(swEsmfTime) CALL ESMF_VMWtime(beg_time)

          CALL neutral ( utime )

          if(swEsmfTime) then  
            CALL ESMF_VMWtime(end_time)
!          if(mype<=1)
            if(mype==0.or.mype==8)write(unit=9999,FMT=*)mype,nTimeStep,"neutral endT=",(end_time-beg_time)
          end if !swEsmfTime


!g          ret = gptlstop  ('neutral')
        END IF ! ( MOD( (utime-start_time),ip_freq_msis)==0 ) THEN
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-7")

      if (SMS_DEBUGGING_ON()) then
      GlobalLower1(:,1:1) = 1
      GlobalUpper2(:,1:1) = 1
      Permute_Decomp3(:,1:1) = 0
      GlobalStart4(:,1:1) = 1
      GlobalStop5(:,1:1) = 1

!-SMS          storing information about variable: plasma_3d

      GlobalLower1(1,1) = LBOUND(plasma_3d,1)
      GlobalLower1(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalLower1(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalLower1(4,1) = LBOUND(plasma_3d,4)
      GlobalUpper2(1,1) = UBOUND(plasma_3d,1)
      GlobalUpper2(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalUpper2(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalUpper2(4,1) = UBOUND(plasma_3d,4)
      Permute_Decomp3(2,1) = 1
      Permute_Decomp3(3,1) = 2
      GlobalStart4(1,1) = LBOUND(plasma_3d,1)
      GlobalStart4(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalStart4(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalStart4(4,1) = LBOUND(plasma_3d,4)
      GlobalStop5(1,1) = UBOUND(plasma_3d,1)
      GlobalStop5(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalStop5(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalStop5(4,1) = UBOUND(plasma_3d,4)

      VarName6= 'plasma_3d'
      ppp__str= 'driver_ipe.f90 - plasma_3d-7'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('module_update_IPE.F90:L210.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()


! update plasma
      IF ( MOD( (utime-start_time),ip_freq_plasma)==0 ) THEN 
!g        ret = gptlstart ('plasma')
        call ESMF_LogWrite("sub-update_IPE plasma:", ESMF_LOGMSG_INFO, rc=rc) 


!nm20161003 esmf timing library
        if(swEsmfTime) CALL ESMF_VMWtime(beg_time)

        if(sw_debug.and.IAM_ROOT())print*,ncnt,'call plasma: utime=',utime,nTimeStep

     call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
     call ESMF_TimeGet(currtime,yy=yy,mm=mm,dd=dd,h=h,m=m,s=s,rc=rc)                                 
     fmt = '(i2.2)'
     fmt1 = '(i4.4)'
     write (yy_str,fmt1) yy
     write (mm_str,fmt) MM
     write (dd_str,fmt) DD
     write (h_str,fmt) H
     write (m_str,fmt) M
     write (s_str,fmt) S
     print *, ' GHGM THISSSS BEGIN TIME ', trim(yy_str)//trim(mm_str)//trim(dd_str)//'T'//trim(h_str)//trim(m_str)     
     call ESMF_ClockGet(clock, stopTime=stopTime, rc=rc)
     call ESMF_TimeGet(stopTime,yy=yy,mm=mm,dd=dd,h=h,m=m,s=s,rc=rc)                                 
     fmt = '(i2.2)'
     fmt1 = '(i4.4)'
     write (yy_str,fmt1) yy
     write (mm_str,fmt) MM
     write (dd_str,fmt) DD
     write (h_str,fmt) H
     write (m_str,fmt) M
     write (s_str,fmt) S
     print *, ' GHGM THISSSS END TIME ', trim(yy_str)//trim(mm_str)//trim(dd_str)//'T'//trim(h_str)//trim(m_str)     
     timestamp_for_IPE_output_files = trim(yy_str)//trim(mm_str)//trim(dd_str)//'T'//trim(h_str)//trim(m_str)                  

        CALL plasma ( utime, timestamp_for_IPE_output_files )

        if(swEsmfTime) then
	  CALL ESMF_VMWtime(end_time)
          if(mype==0.or.mype==8)write(unit=9999,FMT=*)mype,nTimeStep,"plasma endT=",(end_time-beg_time)
        end if


!g        ret = gptlstop  ('plasma')
      END IF   !  ( MOD( (utime-start_time),ip_freq_msis)==0 ) THEN 
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-8")

      if (SMS_DEBUGGING_ON()) then
      GlobalLower1(:,1:1) = 1
      GlobalUpper2(:,1:1) = 1
      Permute_Decomp3(:,1:1) = 0
      GlobalStart4(:,1:1) = 1
      GlobalStop5(:,1:1) = 1

!-SMS          storing information about variable: plasma_3d

      GlobalLower1(1,1) = LBOUND(plasma_3d,1)
      GlobalLower1(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalLower1(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalLower1(4,1) = LBOUND(plasma_3d,4)
      GlobalUpper2(1,1) = UBOUND(plasma_3d,1)
      GlobalUpper2(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalUpper2(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalUpper2(4,1) = UBOUND(plasma_3d,4)
      Permute_Decomp3(2,1) = 1
      Permute_Decomp3(3,1) = 2
      GlobalStart4(1,1) = LBOUND(plasma_3d,1)
      GlobalStart4(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalStart4(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalStart4(4,1) = LBOUND(plasma_3d,4)
      GlobalStop5(1,1) = UBOUND(plasma_3d,1)
      GlobalStop5(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalStop5(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalStop5(4,1) = UBOUND(plasma_3d,4)

      VarName6= 'plasma_3d'
      ppp__str= 'driver_ipe.f90 - plasma_3d-8'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('module_update_IPE.F90:L257.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()


! output UTIME to FLIP_ERR
!g        ret = gptlstart ('output')
        call ESMF_LogWrite("sub-update_IPE output:", ESMF_LOGMSG_INFO, rc=rc)  
        CALL output ( utime )
!g        ret = gptlstop  ('output')
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-9")

      if (SMS_DEBUGGING_ON()) then
      GlobalLower1(:,1:1) = 1
      GlobalUpper2(:,1:1) = 1
      Permute_Decomp3(:,1:1) = 0
      GlobalStart4(:,1:1) = 1
      GlobalStop5(:,1:1) = 1

!-SMS          storing information about variable: plasma_3d

      GlobalLower1(1,1) = LBOUND(plasma_3d,1)
      GlobalLower1(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalLower1(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalLower1(4,1) = LBOUND(plasma_3d,4)
      GlobalUpper2(1,1) = UBOUND(plasma_3d,1)
      GlobalUpper2(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalUpper2(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalUpper2(4,1) = UBOUND(plasma_3d,4)
      Permute_Decomp3(2,1) = 1
      Permute_Decomp3(3,1) = 2
      GlobalStart4(1,1) = LBOUND(plasma_3d,1)
      GlobalStart4(2,1) = dh__LowBounds(1,dh__NestLevel)
      GlobalStart4(3,1) = dh__LowBounds(2,dh__NestLevel)
      GlobalStart4(4,1) = LBOUND(plasma_3d,4)
      GlobalStop5(1,1) = UBOUND(plasma_3d,1)
      GlobalStop5(2,1) = dh__UpperBounds(1,dh__NestLevel)
      GlobalStop5(3,1) = dh__UpperBounds(2,dh__NestLevel)
      GlobalStop5(4,1) = UBOUND(plasma_3d,4)

      VarName6= 'plasma_3d'
      ppp__str= 'driver_ipe.f90 - plasma_3d-9'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('module_update_IPE.G90:L305.0',' ', ppp__status,  &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()
!---copy end

!nm20160608: IPE internal time management
      utime = utime + time_step
      nTimeStep = nTimeStep + 1

!nm20161225 internal time loop
      end do internalTimeLoop !: ncnt=1,3

    call ESMF_LogWrite("sub-update_IPE finished:", ESMF_LOGMSG_INFO, rc=rc)

    rc = ESMF_SUCCESS

  end subroutine update_IPE
end module module_update_IPE
