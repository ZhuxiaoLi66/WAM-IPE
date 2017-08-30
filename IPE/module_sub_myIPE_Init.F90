!nm20140919: subroutine to initialize IPE
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
module module_sub_myIPE_Init
   implicit none

   private

   public :: myIPE_Init
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine myIPE_Init ( gridComp, importState, exportState, clock, rc )
  use ESMF
!nm20160307 SMS should be included
  use nnt_types_module
  use  module_decomp
!---
  USE module_myIPE_Init
  USE module_input_parameters, ONLY: read_input_parameters,utime,start_time,stop_time, &
                                     HPEQ_flip,sw_perp_transport,NYEAR,NDAY,sw_output_plasma_grid, &
                                     parallelBuild,mype, MPI_COMM_IPE, ut_start_perp_trans
  USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d
  USE module_open_output_files,ONLY: open_output_files
  USE module_init_plasma_grid, ONLY: init_plasma_grid
!SMS$IGNORE BEGIN
  USE module_init_ELDYN,       ONLY: init_eldyn
!SMS$IGNORE END
  USE module_IO,ONLY: PRUNIT
  USE module_open_file,ONLY: open_file
  implicit none
!g  include "gptl.inc"
  include "mpif.h"
!---
  type(ESMF_GridComp)  :: gridComp
  type(ESMF_State)     :: importState
  type(ESMF_State)     :: exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc
!---MPI communicator
       integer      :: ierr, size
       type(ESMF_VM):: vm_IPE !virtual machine  !?is this vm_local or vm_global???
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
!---clock
       type(ESMF_TimeInterval)   :: timeStep
       type(ESMF_Time)           :: startTime
       type(ESMF_Time)           :: stopTime
!       type(ESMF_TimeInterval)   :: runDuration
!       real(ESMF_KIND_R8)        :: runTimeStepCount
!       type(ESMF_Time)           :: refTime
       type(ESMF_Time)           :: currTime
!       type(ESMF_Time)           :: prevTime
!       type(ESMF_TimeInterval)   :: currSimTime
!       type(ESMF_TimeInterval)   :: prevSimTime
       type(ESMF_Calendar)       :: calendar
!       type(ESMF_CalKind_Flag)   :: calkindflag
!       integer                   :: timeZone
!       integer(ESMF_KIND_I8)     :: advanceCount
!       integer                   :: alarmCount
!       type(ESMF_Direction_Flag) :: direction
!t       character (len=*)         :: name
!---time
!       integer(ESMF_KIND_I4) :: yy
!       integer(ESMF_KIND_I8),   intent(out), optional :: yy_i8
       integer :: mm
       integer :: dd
!       integer(ESMF_KIND_I4),   intent(out), optional :: d
       integer(ESMF_KIND_I8) :: d_i8
       integer(ESMF_KIND_I4) :: h
       integer(ESMF_KIND_I4) :: m
       integer(ESMF_KIND_I4) :: s
       integer(ESMF_KIND_I4) :: yystop
       integer(ESMF_KIND_I4) :: mmstop
       integer(ESMF_KIND_I4) :: ddstop
       integer(ESMF_KIND_I4) :: hstop
       integer(ESMF_KIND_I4) :: mstop
       integer(ESMF_KIND_I4) :: sstop
       character(len=8) :: fmt ! format descriptor
       character(len=8) :: fmt1 ! format descriptor
       character(len=8) :: yy_str
       character(len=8) :: mm_str
       character(len=8) :: dd_str
       character(len=8) :: h_str
       character(len=8) :: m_str
       character(len=8) :: s_str
       character(len=8) :: yy_stop_str
       character(len=8) :: mm_stop_str
       character(len=8) :: dd_stop_str
       character(len=8) :: h_stop_str
       character(len=8) :: m_stop_str
       character(len=8) :: s_stop_str
       character(len=56) :: model_start_time
       character(len=56) :: model_stop_time
!       integer(ESMF_KIND_I8),   intent(out), optional :: s_i8
!       integer(ESMF_KIND_I4),   intent(out), optional :: ms
!       integer(ESMF_KIND_I4),   intent(out), optional :: us
!       integer(ESMF_KIND_I4),   intent(out), optional :: ns
!       real(ESMF_KIND_R8),      intent(out), optional :: d_r8
!       real(ESMF_KIND_R8),      intent(out), optional :: h_r8
!       real(ESMF_KIND_R8),      intent(out), optional :: m_r8
!       real(ESMF_KIND_R8),      intent(out), optional :: s_r8
!       real(ESMF_KIND_R8),      intent(out), optional :: ms_r8
!       real(ESMF_KIND_R8),      intent(out), optional :: us_r8
!       real(ESMF_KIND_R8),      intent(out), optional :: ns_r8
!       integer(ESMF_KIND_I4),   intent(out), optional :: sN
!       integer(ESMF_KIND_I8),   intent(out), optional :: sN_i8
!       integer(ESMF_KIND_I4),   intent(out), optional :: sD
!       integer(ESMF_KIND_I8),   intent(out), optional :: sD_i8
!       type(ESMF_Calendar),     intent(out), optional :: calendar
!       type(ESMF_CalKind_Flag), intent(out), optional :: calkindflag
!       integer,                 intent(out), optional :: timeZone ! not imp
!       character (len=*),       intent(out), optional :: timeString
!       character (len=*),       intent(out), optional :: timeStringISOFrac
!       integer,                 intent(out), optional :: dayOfWeek
!       type(ESMF_Time),         intent(out), optional :: midMonth
!       integer(ESMF_KIND_I4),   intent(out), optional :: dayOfYear
!       real(ESMF_KIND_R8),      intent(out), optional :: dayOfYear_r8
!       type(ESMF_TimeInterval), intent(out), optional :: dayOfYear_intvl
!---calendarget
!       type(ESMF_CalKind_Flag),intent(out), optional :: calkindflag
       integer :: daysPerMonth(12)
       integer :: monthsPerYear
       integer(ESMF_KIND_I4) :: secondsPerDay
!       integer(ESMF_KIND_I4) :: secondsPerYear
!       integer(ESMF_KIND_I4) :: daysPerYear
!       integer(ESMF_KIND_I4),  intent(out), optional :: daysPerYearDn
!       integer(ESMF_KIND_I4),  intent(out), optional :: daysPerYearDd
!       character (len=*),      intent(out), optional :: name
!nm20160723 gptl
      INTEGER :: ret
!---

  call ESMF_LogWrite("sub-init_IPE STARTING:", ESMF_LOGMSG_INFO, rc=rc)

!t  if (IAM_ROOT()) then
    print *,"sub-init_IPE: ESMF_ClockPrint"
!t  endif
  CALL ESMF_LogWrite("sub-init_IPE: ESMF_ClockPrint", ESMF_LOGMSG_INFO, rc=rc)
  call ESMF_ClockPrint(clock, options="currTime name startTime stopTime timeStep", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
      	  file=__FILE__)) &
          return  ! bail out
  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)


!t  if (IAM_ROOT()) then
    print *,"sub-init_IPE: ESMF_ClockGet"
!t  endif
  CALL ESMF_LogWrite("sub-init_IPE: ESMF_ClockGet", ESMF_LOGMSG_INFO, rc=rc)
  call ESMF_ClockGet(clock, &
         timeStep=timeStep, &
         startTime=startTime, stopTime=stopTime, &
         currTime=currTime, &
         calendar=calendar, &
         rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
      	  file=__FILE__)) &
          return  ! bail out

!t  if (IAM_ROOT()) then
    print *,"sub-init_IPE: ESMF_TimePrint"
!t  endif

  CALL ESMF_LogWrite("sub-init_IPE: ESMF_TimePrint", ESMF_LOGMSG_INFO, rc=rc)
  call ESMF_TimePrint(currTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
      	  file=__FILE__)) &
          return  ! bail out


!t  if (IAM_ROOT()) then
    print *,"sub-init_IPE: ESMF_TimeGet"
!t  endif
  CALL ESMF_LogWrite("sub-init_IPE: ESMF_TimeGet", ESMF_LOGMSG_INFO, rc=rc)

call ESMF_TimeGet(currtime, yy=NYEAR,mm=mm,dd=dd,d_i8=d_i8,h=h,m=m,s=s,rc=rc)

if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
      	  file=__FILE__)) &
          return  ! bail out
!t   if (IAM_ROOT()) then
     print *,'NYEAR=', NYEAR,' d_i8=',d_i8
     print *, "The clock's final current time is NYEAR=", NYEAR, "/MM=", MM, "/DD=", DD, &
               " H=", H, ":M=", M, ":S=", S
     print *, "GHGM START DATETIME ", NYEAR, MM, DD, H, M, S
     fmt = '(i2.2)'
     fmt1 = '(i4.4)'
     write (yy_str,fmt1) NYEAR
     write (mm_str,fmt) MM
     write (dd_str,fmt) DD
     write (h_str,fmt) H
     write (m_str,fmt) M
     write (s_str,fmt) S
     print *, trim(yy_str)//trim(mm_str)//trim(dd_str)//'T'//trim(h_str)//trim(m_str)
     model_start_time = &
&          'MODEL_START_TIME.'//trim(yy_str)//trim(mm_str)//trim(dd_str)//'T'//trim(h_str)//trim(m_str)
     open(unit=5995,FILE=model_start_time,status='unknown',form='formatted')
     write(5995,*) model_start_time
     close(5995)
call ESMF_TimeGet(stoptime, yy=yystop,mm=mmstop,dd=ddstop,h=hstop,m=mstop,s=sstop,rc=rc)
     write (yy_stop_str,fmt1) yystop
     write (mm_stop_str,fmt) mmstop
     write (dd_stop_str,fmt) ddstop
     write (h_stop_str,fmt) hstop
     write (m_stop_str,fmt) mstop
     write (s_stop_str,fmt) sstop
     print *, trim(yy_stop_str)//trim(mm_stop_str)//trim(dd_stop_str)//'T'//trim(h_stop_str)//trim(m_stop_str)
     model_stop_time = &
&          'MODEL_STOP_TIME.'//trim(yy_stop_str)//trim(mm_stop_str)//trim(dd_stop_str)//'T'//trim(h_stop_str)//trim(m_stop_str)
     open(unit=5994,FILE=model_stop_time,status='unknown',form='formatted')
     write(5994,*) model_stop_time
     close(5994)
!t   endif

!input to IPE (values in IPE.inp will no longer be used!)
  ! get NDAY
!t     if (IAM_ROOT()) then
       print *,"sub-init_IPE: ESMF_CalendarGet"
!t     endif
     CALL ESMF_LogWrite("sub-init_IPE: CalendarGet", ESMF_LOGMSG_INFO, rc=rc)
     call ESMF_CalendarGet(calendar, daysPerMonth=daysPerMonth, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
      	  file=__FILE__)) &
          return  ! bail out

     print *,"Days Per Month ",daysPerMonth
     NDAY = dd + SUM ( daysPerMonth(1:mm-1) )
     print *,'GHGM IPE INIT : NDAY = ', NDAY
     start_time = H*3600 + M*60 + S  !UT in seconds
     print *,'GHGM IPE INIT : start_time = ', start_time
 

! i would need to get corresponding kp & ap here...  


! Wrapping layer in which native arrays are extracted from model
  ! data structures, and referenced or copied into ESMF Arrays,
  ! Array Bundles, Fields, or Field Bundles.  References to these
  ! objects are then placed into import and export States.

  ! Scientific content of initialize routine goes here.

!--- obtain global vm information
      CALL ESMF_LogWrite("sub-init_IPE: ESMF_VMGetCurrent", ESMF_LOGMSG_INFO, rc=rc)
      call ESMF_VMGetCurrent(vm=vm_IPE, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
      	  file=__FILE__)) &
          return  ! bail out

!--- obtain MPI communicator:
      CALL ESMF_LogWrite("sub-init_IPE: ESMF_VMGet", ESMF_LOGMSG_INFO, rc=rc)
      call ESMF_VMGet(vm=vm_IPE, mpiCommunicator=MPI_COMM_IPE, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      	 line=__LINE__, &
         file=__FILE__)) &
     	 return  ! bail out

      call MPI_Barrier(MPI_COMM_IPE, ierr)
!t if (IAM_ROOT()) then
      print *,'NUOPC MPI communicator=',MPI_COMM_IPE
!t endif

! after obtaining comm from ESMF
      call MPI_Comm_size(MPI_COMM_IPE, size, ierr)
!t if (IAM_ROOT()) then
      print *, "Number of PETs in communicator: ", size 
!t endif

!nm20160315---start copy driver_ipe_sms.f90
  call sms_start(ppp__status)

  call nnt_chkstat('module_sub_myIPE_Init.F90:L272',' ',ppp__status,          &
     &   NNT_ABORT_ON_ERROR, ppp__status)


!g      call gptlprocess_namelist ('GPTLnamelist', 77, ret) 
!g      ret = gptlinitialize ()
!g      ret = gptlstart ('Total')

!!SMS$INSERT parallelBuild=.true.
 parallelBuild=.true.
! set up input parameters
!g      ret = gptlstart ('read_input')
      CALL ESMF_LogWrite("sub-initialize_IPE: read_input_parameters", ESMF_LOGMSG_INFO, rc=rc)
      CALL read_input_parameters ( )
      !ut_start_perp_trans = 0
      ut_start_perp_trans = start_time 

!g      ret = gptlstop  ('read_input')

! open Input/Output files
!g      ret = gptlstart ('open_output_files')
!--- unit=8
!filename ='FLIP_ERROR_FLAG_'//TRIM(string_tmp)//'.log'
      CALL open_file ( 'FLIP_ERR', PRUNIT,'formatted','unknown') !open by all processors
!!SMS$SERIAL BEGIN

!t      if (IAM_ROOT()) then
      CALL ESMF_LogWrite("sub-initialize_IPE: open_output_files", ESMF_LOGMSG_INFO, rc=rc)	
      CALL open_output_files ( )
!t      endif
!!SMS$SERIAL END
!g      ret = gptlstop  ('open_output_files')

! set up plasma grids by reading file
!g      ret = gptlstart ('init_plasma_grid')
      CALL ESMF_LogWrite("sub-initialize_IPE: init_plasma_grid", ESMF_LOGMSG_INFO, rc=rc)
      CALL init_plasma_grid ( )
!g      ret = gptlstop  ('init_plasma_grid')

!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-1")

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
      ppp__str= 'driver_ipe.f90 - plasma_3d-1'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('driver_ipe_sms.f90:133.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()
IF ( sw_output_plasma_grid ) THEN
!g  ret = gptlstart ('output_plasma_grid')

      if (IAM_ROOT()) then
  print *, 'sub-init_p: output plasma_grid'
      endif
  CALL ESMF_LogWrite("sub-initialize_IPE: output_plasma_grid", ESMF_LOGMSG_INFO, rc=rc)
  CALL output_plasma_grid ( )
!g  ret = gptlstop  ('output_plasma_grid')
END IF
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-2")

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
      ppp__str= 'driver_ipe.f90 - plasma_3d-2'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('driver_ipe_sms.f90:181.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()

! initialise the flux tubes from previous runs
      IF ( HPEQ_flip==0.0 ) THEN

        if (IAM_ROOT()) then
      	  print *,'before CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
        endif
!g        ret = gptlstart ('io_plasma_bin')
        CALL ESMF_LogWrite("sub-initialize_IPE: io_plasma_bin", ESMF_LOGMSG_INFO, rc=rc)
        CALL io_plasma_bin ( 2, start_time )
!g        ret = gptlstop  ('io_plasma_bin')
        if (IAM_ROOT()) then
          print *,'after CALL io_plasma_bin finished! READ: start_time=', start_time,' stop_time=',stop_time
        endif

      END IF
!nm20160608: IPE internal time management
      utime = start_time
 
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-3")

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
      ppp__str= 'driver_ipe.f90 - plasma_3d-3'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('driver_ipe_sms.f90:236.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()


! initialization of electrodynamic module:
! read in E-field
!g      ret = gptlstart ('init_eldyn')
      IF ( sw_perp_transport>=1 ) THEN
        CALL ESMF_LogWrite("sub-initialize_IPE: init_eldyn", ESMF_LOGMSG_INFO, rc=rc)
        CALL init_eldyn ( )
      ENDIF
!g      ret = gptlstop  ('init_eldyn')
!!sms$compare_var(plasma_3d,"driver_ipe.f90 - plasma_3d-4")

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
      ppp__str= 'driver_ipe.f90 - plasma_3d-4'
      call PPP_COMPARE_VAR(dh(dh__NestLevel),plasma_3d,NNT_REAL,        &
     &   GlobalUpper2(:,1),Permute_Decomp3(:,1),GlobalStart4(:,1),      &
     &   GlobalStop5(:,1),GlobalLower1(:,1),4,VarName6,ppp__str,        &
     &   ppp__status)
      call NNT_CHKSTAT('driver_ipe_sms.f90:285.0',' ', ppp__status,   &
     &   NNT_ABORT_ON_ERROR, ppp__status)
      endif    ! end of SMS_DEBUGGING_ON()

!nm20160315---end copy driver_ipe_sms.f90

  call ESMF_LogWrite("sub-initialize_IPE finished:", ESMF_LOGMSG_INFO, rc=rc)

  rc = ESMF_SUCCESS

  end subroutine myIPE_Init
end module module_sub_myIPE_Init
