















!-----------------------------------------------------------------------
!
      MODULE module_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the ATM component.
!-----------------------------------------------------------------------
!
!***  The ATM component lies in the hierarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!              /|!             / | !          ATM/OCN/ICE components
!          |    |
!          |    |
!          |    |
!          |    (MOM5, HYCOM, etc.)
!          |
!          CORE component (GSM, NMM, FIM, GEN, etc.)
!
!-----------------------------------------------------------------------
!  2011-05-11  Theurich & Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  2011-10/04  Yang  - Modified for using the ESMF 5.2.0r library.
!  2013-07     Theurich - NUOPC option to be compliant with ESMF 6.2.0 reference.
!-----------------------------------------------------------------------
!
      USE ESMF

      use NUOPC
      use NUOPC_Model, only: &
        model_routine_SS            => SetServices, &
        model_label_DataInitialize  => label_DataInitialize, &
        model_label_CheckImport     => label_CheckImport, &
        model_label_Advance         => label_Advance
      use module_CPLFIELDS

!
      USE module_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE            &
                                         ,WRAP_ATM_INTERNAL_STATE
!
      USE module_NMM_GRID_COMP,ONLY: NMM_REGISTER
      USE module_GFS_GRID_COMP,ONLY: GFS_REGISTER
      USE module_FIM_GRID_COMP,ONLY: FIM_REGISTER
      USE module_GEN_GRID_COMP,ONLY: GEN_REGISTER   ! For the "Generic Core" gridded component.
!
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!

!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE

  integer, parameter :: MAXNAMELEN = 128

  ! private internal state to keep instance data
  type InternalStateStruct
    integer(ESMF_KIND_I4), pointer :: numPerRow(:)
    real(ESMF_KIND_R8), pointer :: lons(:,:), lats(:,:)
    integer, pointer :: rowinds(:), indList(:)
    integer :: dims(2)
    integer :: myrows
    integer :: wamtotalnodes, localnodes
    integer :: PetNo, PetCnt
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type
!
      PUBLIC :: ATM_REGISTER
!
!-----------------------------------------------------------------------
!
      TYPE(ATM_INTERNAL_STATE),POINTER,SAVE :: ATM_INT_STATE
      TYPE(WRAP_ATM_INTERNAL_STATE)   ,SAVE :: WRAP
!
      TYPE(ESMF_Clock),SAVE :: CLOCK_ATM                                   !<-- The Clock of the ATM component
!
      character(len=160) :: nuopcMsg
      logical :: write_diagnostics = .false.
      logical :: profile_memory = .false.

!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_REGISTER(ATM_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  Register the Init, Run, and Finalize routines of 
!***  the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      INTEGER            ,INTENT(OUT)   :: RC_REG                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      ! the NUOPC model component will register the generic methods
      call NUOPC_CompDerive(ATM_GRID_COMP, model_routine_SS, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=137, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
        
      ! Provide InitializeP0 to switch from default IPDv00 to IPDv02
      call ESMF_GridCompSetEntryPoint(ATM_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        InitializeP0, phase=0, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=145, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      ! IPDv02p1: advertise Fields
      call NUOPC_CompSetEntryPoint(ATM_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        phaseLabelList=(/"IPDv02p1"/), userRoutine=InitializeIPDv02p1, &
        rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=154, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      ! IPDv02p3: realize connected Fields
      call NUOPC_CompSetEntryPoint(ATM_GRID_COMP, ESMF_METHOD_INITIALIZE, &
        phaseLabelList=(/"IPDv02p3"/), userRoutine=InitializeIPDv02p3, &
        rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=163, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      ! attach specializing method(s)
      call NUOPC_CompSpecialize(ATM_GRID_COMP, &
        specLabel=model_label_DataInitialize, specRoutine=ATM_DATAINIT, &
        rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=172, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      call NUOPC_CompSpecialize(ATM_GRID_COMP, &
        specLabel=model_label_Advance, specRoutine=ATM_ADVANCE, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=178, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      call ESMF_MethodRemove(ATM_GRID_COMP, model_label_CheckImport, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=183, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      call NUOPC_CompSpecialize(ATM_GRID_COMP, &
        specLabel=model_label_CheckImport, specRoutine=ATM_CHECKIMPORT, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=189, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      ! Overwrite generic NUOPC_Model Finalize method
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP, ESMF_METHOD_FINALIZE, &
        ATM_FINALIZE, rc=RC_REG)
      if (ESMF_LogFoundError(rcToCheck=RC_REG, msg=ESMF_LOGERR_PASSTHRU, &
        line=197, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_REGISTER succeeded'
      ELSE
        WRITE(0,*)' ATM_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_REGISTER
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    character(len=10)                         :: value
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv02 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv02"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=231, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out

    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=238, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
    write_diagnostics=(trim(value)=="true")

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=246, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
    profile_memory=(trim(value)/="false")
    
  end subroutine

  !-----------------------------------------------------------------------

  subroutine InitializeIPDv02p1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS

    ! advertise Fields

    ! importable fields:
    call NUOPC_Advertise(importState, StandardNames=ImportFieldsList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=268, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out

    ! exportable fields:
    call NUOPC_Advertise(exportState, StandardNames=ExportFieldsList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=275, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------

  subroutine InitializeIPDv02p3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer                         :: n
    type(ESMF_Grid)                 :: gridIn, gridOut
    type(ESMF_array)                :: array

    rc = ESMF_SUCCESS
    
    ! call into the actual NEMS/ATM initialize routine
    
    call ATM_INITIALIZE(gcomp, importState, exportState, clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=299, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
    
    ! use the regular Gaussian Grid that was setup during the NEMS/ATM init
    gridIn  = gauss2d ! for imported Fields
    gridOut = gauss2d ! for exported Fields
    if (atm_int_state%CORE == "nmm") then
    endif

    ! conditionally realize or remove Fields from States ...

    do n = 1,nImportFields
      call realizeConnectedInternCplField(importState, &
        field=importFields(n), standardName=trim(importFieldsList(n)), &
        grid=gridIn, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=320, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
    enddo

    do n = 1,nExportFields
      call realizeConnectedInternCplField(exportState, &
        field=exportFields(n), standardName=trim(exportFieldsList(n)), &
        grid=gridOut, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=330, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
    enddo

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    subroutine realizeConnectedInternCplField(state, field, standardName, grid, rc)
      type(ESMF_State)                :: state
      type(ESMF_Field)                :: field
      character(len=*)                :: standardName
      type(ESMF_Grid)                 :: grid
      integer, intent(out), optional  :: rc
      
      ! local variables
      character(len=80)               :: fieldName
      type(ESMF_ArraySpec)            :: arrayspec
      integer                         :: i

      if (present(rc)) rc = ESMF_SUCCESS
      
      fieldName = standardName  ! use standard name as field name

      !! Create fields using wam2dmesh if they are WAM fields 
      if (NUOPC_IsConnected(state, fieldName=fieldName)) then
         if (fieldName == "northward_wind_neutral" .or. &
             fieldName == "eastward_wind_neutral" .or. &
             fieldName == "upward_wind_neutral" .or. &
             fieldName == "temp_neutral" .or. &
             fieldName == "O_Density" .or. &
             fieldName == "O2_Density" .or. &
             fieldName == "N2_Density" .or. &
             fieldName == "height")  then
           call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=365, &
              file="module_ATM_GRID_COMP.F90")) &
              return  ! bail out
           field = ESMF_FieldCreate(wam2dmesh, arrayspec, &
              ungriddedLBound=(/1/), ungriddedUBound=(/wamlevels/),&
              name=fieldName, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=372, &
              file="module_ATM_GRID_COMP.F90")) &
              return  ! bail out
         else
           ! realize the connected Field pass back up for internal cpl fields
           field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=fieldName, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=379, &
             file="module_ATM_GRID_COMP.F90")) &
             return  ! bail out
         endif
         call NUOPC_Realize(state, field=field, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=385, &
           file="module_ATM_GRID_COMP.F90")) &
           return  ! bail out
      else
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/fieldName/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=392, &
          file="module_ATM_GRID_COMP.F90")) &
          return  ! bail out
      endif
    end subroutine

  end subroutine
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!

      SUBROUTINE ATM_INITIALIZE(ATM_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_EARTH                             &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  The Initialize step of the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The ATM import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The ATM export state
      TYPE(ESMF_Clock)                  :: CLOCK_EARTH                     !<-- The Clock of the EARTH component
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
      TYPE(ESMF_Config)       :: CF
      real(ESMF_KIND_R8)      :: medAtmCouplingIntervalSec
      type(ESMF_Clock)        :: atmClock
      type(ESMF_TimeInterval) :: atmStep, earthStep
      type(ESMF_Time)         :: currTime, stopTime
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_INIT = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
      call ESMF_ClockPrint(CLOCK_EARTH, options="currTime", &
        preString="entering ATM_INITIALIZE with CLOCK_EARTH current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_EARTH, options="startTime", &
        preString="entering ATM_INITIALIZE with CLOCK_EARTH start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_EARTH, options="stopTime", &
        preString="entering ATM_INITIALIZE with CLOCK_EARTH stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

!
!-----------------------------------------------------------------------
!***  Allocate the ATM component's internal state, point at it,
!***  and attach it to the ATM component.
!-----------------------------------------------------------------------
!
      ALLOCATE(ATM_INT_STATE,stat=RC)
      wrap%ATM_INT_STATE=>ATM_INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the ATM Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(ATM_GRID_COMP                  &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
!-----------------------------------------------------------------------
!***  Create the configure object for the ATM configure file which
!***  specifies the dynamic core.
!-----------------------------------------------------------------------
!
      CF=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Load the ATM configure file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigLoadFile(config=CF ,filename='atmos.configure' ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the configure object to the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the configure file to the ATM component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM component
                           ,config  =CF                                 &  !<-- The associated configure object
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
!-----------------------------------------------------------------------
!***  Setup a Clock instance that the ATM core is allowed to modify.
!***  This Clock is NOT under direct NUOPC control!
!-----------------------------------------------------------------------
!
      atm_int_state%CLOCK_ATM = ESMF_ClockCreate(CLOCK_EARTH, rc=RC_INIT)
      if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=521, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return

!-----------------------------------------------------------------------
!***  Extract the dynamic core name from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract dynamic core from the ATM configure file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ATM configure object
                                  ,value =atm_int_state%CORE            &  !<-- The dynamic core name
                                  ,label ='atm_model:'                  &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ATM subcomponent and its associated import/export
!***  states for the core name that was extracted.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      atm_int_state%CORE_GRID_COMP=ESMF_GridCompCreate(name=TRIM(atm_int_state%CORE)//' component' &
                                                      ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the configure object to the CORE component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     MESSAGE_CHECK="Attach the configure file to the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!     CALL ESMF_GridCompSet(gridcomp=atm_int_state%CORE_GRID_COMP       &  !<-- The ATM component
!                          ,config  =CF                                 &  !<-- The associated configure object
!                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the subcomponent's Init, Run, and Finalize subroutines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the CORE component's Init, Run, and Finalize steps"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      SELECT CASE(atm_int_state%CORE)
!
        CASE('nmm')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,NMM_REGISTER                   &
                                        ,rc=RC)
!
        CASE('gfs')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GFS_REGISTER                   &
                                        ,rc=RC)
!
        CASE('gsm')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GFS_REGISTER                   &
                                        ,rc=RC)
!
        CASE('fim')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,FIM_REGISTER                   &
                                        ,rc=RC)

        CASE('gen')
          CALL ESMF_GridCompSetServices (atm_int_state%CORE_GRID_COMP   &
                                        ,GEN_REGISTER                   &
                                        ,rc=RC)
        CASE DEFAULT
          write(0,*)' ATM_INITIALIZE requires unknown core: ',TRIM(atm_int_state%CORE)                      
!
      END SELECT
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Core component's import/export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the CORE import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      atm_int_state%CORE_IMP_STATE=ESMF_StateCreate(name   = "CORE Import"            &
                                                   ,stateintent = ESMF_STATEINTENT_IMPORT  &
                                                   ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the CORE export state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      atm_int_state%CORE_EXP_STATE=ESMF_StateCreate(name   = "CORE Export"            &
                                                   ,stateintent = ESMF_STATEINTENT_EXPORT  &
                                                   ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Nest the import/export states of the CORE component into the
!***  analgous states of the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK= "Add the CORE states into the ATMOS states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Initialize the CORE component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%CORE_GRID_COMP &
                                  ,importState=atm_int_state%CORE_IMP_STATE &
                                  ,exportState=atm_int_state%CORE_EXP_STATE &
                                  ,clock      =atm_int_state%CLOCK_ATM      &
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_INITIALIZE succeeded'
      ELSE
        WRITE(0,*)' ATM_INITIALIZE failed  RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
! - Under NUOPC, the EARTH driver clock is a separate instance from the 
! - ATM clock. However, the ATM clock may have been reset during ATM initialize
! - and therefore the EARTH driver clock must also be adjusted.
! - Affected: currTime, timeStep
      call ESMF_ClockGet(atm_int_state%CLOCK_ATM, currTime=currTime, &
        stopTime=stopTime, rc=RC_INIT)
      if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=706, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return
      call ESMF_ClockGet(CLOCK_EARTH, timeStep=earthStep, rc=RC_INIT)
      if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=708, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return

      if (earthStep>(stopTime-currTime)) earthStep=stopTime-currTime
      call ESMF_ClockSet(CLOCK_EARTH, currTime=currTime, &
        timeStep=earthStep, rc=RC_INIT)
      if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=713, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return

      ! Set ATM component clock as copy of EARTH clock.
      call NUOPC_CompSetClock(ATM_GRID_COMP, CLOCK_EARTH, rc=RC_INIT)
      if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=717, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return

      ! Read in the ATM coupling interval
      call ESMF_ConfigGetAttribute(CF, medAtmCouplingIntervalSec, &
        label="atm_coupling_interval_sec:", default=-1.0_ESMF_KIND_R8, &
        rc=RC_INIT)
      if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=723, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return

      if (medAtmCouplingIntervalSec>0._ESMF_KIND_R8) then
        ! The coupling time step was provided
        call ESMF_TimeIntervalSet(atmStep, s_r8=medAtmCouplingIntervalSec, &
          rc=RC_INIT)
        if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=729, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return
        call ESMF_GridCompGet(ATM_GRID_COMP, clock=atmClock, rc=RC_INIT)
        if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=731, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return
        call ESMF_ClockSet(atmClock, timestep=atmStep, rc=RC_INIT)
        if (ESMF_LogFoundError(RC_INIT, msg="Breaking out of subroutine", line=733, file="module_ATM_GRID_COMP.F90", rcToReturn=RC_INIT)) return
      endif

      call ESMF_ClockPrint(CLOCK_EARTH, options="currTime", &
        preString="leaving  ATM_INITIALIZE with CLOCK_EARTH current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_EARTH, options="startTime", &
        preString="leaving  ATM_INITIALIZE with CLOCK_EARTH start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_EARTH, options="stopTime", &
        preString="leaving  ATM_INITIALIZE with CLOCK_EARTH stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="currTime", &
        preString="leaving  ATM_INITIALIZE with CLOCK_ATM current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="startTime", &
        preString="leaving  ATM_INITIALIZE with CLOCK_ATM start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="stopTime", &
        preString="leaving  ATM_INITIALIZE with CLOCK_ATM stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
!
!-----------------------------------------------------------------------
!

      END SUBROUTINE ATM_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!

  subroutine ATM_DATAINIT(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_State)              :: exportState

    rc = ESMF_SUCCESS
    
    ! the ATM initializes export Fields that the MED initialize depends on

    ! query the Component for its exportState
    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=786, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
    
    ! -> set Updated Field Attribute to "true", indicating to the IPDv02p5
    ! generic code to set the timestamp for this Field
    
    call setAllFieldsUpdated(exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=795, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
      
    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=804, &
      file="module_ATM_GRID_COMP.F90")) &
      return  ! bail out
        
  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    subroutine setAllFieldsUpdated(state, rc)
      type(ESMF_State)                :: state
      integer, intent(out), optional  :: rc
      
      integer                         :: i, fieldCount
      character(len=80), allocatable  :: fieldNameList(:)
      type(ESMF_Field)                :: field
      type(ESMF_StateItem_Flag)       :: itemType
      real(ESMF_KIND_R8), pointer     :: fptr(:,:)

      if (present(rc)) rc = ESMF_SUCCESS
      
      call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=824, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      allocate(fieldNameList(fieldCount))
      call ESMF_StateGet(State, itemNameList=fieldNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=830, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      
      do i=1, fieldCount
        call ESMF_StateGet(state, itemName=fieldNameList(i), &
          itemType=itemType, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=838, &
          file="module_ATM_GRID_COMP.F90")) &
          return  ! bail out
        if (itemType /= ESMF_STATEITEM_NOTFOUND) then
          ! item exists -> set "Updated"
          call ESMF_StateGet(state, itemName=fieldNameList(i), field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=845, &
            file="module_ATM_GRID_COMP.F90")) &
            return  ! bail out
	  call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=850, &
            file="module_ATM_GRID_COMP.F90")) &
            return  ! bail out
          fptr=0.d0 ! zero out the entire field
          call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=856, &
            file="module_ATM_GRID_COMP.F90")) &
            return  ! bail out
        endif
      enddo
    end subroutine
    
  end subroutine
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_ADVANCE(ATM_GRID_COMP, rc)
!
!-----------------------------------------------------------------------
!***  Advance the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      type(ESMF_GridComp)   :: ATM_GRID_COMP
      integer, intent(out)  :: rc
!
!---------------------
!***  Local Variables
!---------------------
!
      type(ESMF_Clock)              :: clock
      type(ESMF_Time)               :: stopTime, currTime
      type(ESMF_State)              :: importState, exportState
      type(ESMF_Field)              :: field
      type(ESMF_StateItem_Flag)     :: itemType
      
      !TODO: move the slice counter into an internal state to be instance safe
      integer, save                 :: slice=1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
    if(profile_memory) call ESMF_VMLogMemInfo("Entering ATM ATM_ADVANCE ")
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Use the internal Clock set by NUOPC layer for ATM but update stopTime
!-----------------------------------------------------------------------

      ! Component internal Clock gets updated per NUOPC rules
      call ESMF_GridCompGet(ATM_GRID_COMP, clock=clock, rc=rc)
      if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=908, file="module_ATM_GRID_COMP.F90", rcToReturn=rc)) return
      
      ! The stopTime will be updated to be the next ATM-OCN coupling time
      call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
      if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=912, file="module_ATM_GRID_COMP.F90", rcToReturn=rc)) return
      
      ! Set the ATM-OCN coupling time to be stopTime in Clock that ATM core uses
      call ESMF_ClockSet(atm_int_state%CLOCK_ATM, currTime=currTime, &
        stopTime=stopTime, rc=rc)
      if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=917, file="module_ATM_GRID_COMP.F90", rcToReturn=rc)) return

      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="currTime", &
        preString="entering ATM_ADVANCE with CLOCK_ATM current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="startTime", &
        preString="entering ATM_ADVANCE with CLOCK_ATM start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="stopTime", &
        preString="entering ATM_ADVANCE with CLOCK_ATM stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

!-----------------------------------------------------------------------
!***  Execute the Run step of the selected dynamic core.
!-----------------------------------------------------------------------

      if(profile_memory) call ESMF_VMLogMemInfo("Entering ATM GridCompRun ")
      CALL ESMF_GridCompRun(gridcomp   =atm_int_state%CORE_GRID_COMP    &
                           ,importState=atm_int_state%CORE_IMP_STATE    &
                           ,exportState=atm_int_state%CORE_EXP_STATE    &
                           ,clock      =atm_int_state%CLOCK_ATM         &
                           ,rc         =rc)
      if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=942, file="module_ATM_GRID_COMP.F90", rcToReturn=rc)) return
      if(profile_memory) call ESMF_VMLogMemInfo("Leaving ATM GridCompRun ")

!-----------------------------------------------------------------------

      ! query the Component for its importState and exportState
      call ESMF_GridCompGet(ATM_GRID_COMP, exportState=exportState, &
        importState=importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=951, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out
      
      if(write_diagnostics) then
        ! for testing write all of the Fields in the importState to file
        call NUOPC_Write(importState, fileNamePrefix="field_atm_import_", &
          timeslice=slice, relaxedFlag=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=960, &
          file="module_ATM_GRID_COMP.F90")) &
          return  ! bail out
        ! for testing write all of the Fields in the exportState to file
        call NUOPC_Write(exportState, fileNamePrefix="field_atm_export_", &
          timeslice=slice, relaxedFlag=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=967, &
          file="module_ATM_GRID_COMP.F90")) &
          return  ! bail out
        ! advance the time slice counter
        slice = slice + 1
      endif

      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="currTime", &
        preString="leaving  ATM_ADVANCE with CLOCK_ATM current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="startTime", &
        preString="leaving  ATM_ADVANCE with CLOCK_ATM start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(atm_int_state%CLOCK_ATM, options="stopTime", &
        preString="leaving  ATM_ADVANCE with CLOCK_ATM stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
    if(profile_memory) call ESMF_VMLogMemInfo("Leaving ATM ATM_ADVANCE ")

!-----------------------------------------------------------------------

      END SUBROUTINE ATM_ADVANCE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

      SUBROUTINE ATM_CHECKIMPORT(ATM_GRID_COMP, rc)
!
!-----------------------------------------------------------------------
!***  Check the import state fields
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      type(ESMF_GridComp)   :: ATM_GRID_COMP
      integer, intent(out)  :: rc
!
!---------------------
!***  Local Variables
!---------------------
!
      integer :: n, nf
      type(ESMF_Clock)   :: clock
      type(ESMF_Time)    :: currTime, invalidTime
      type(ESMF_State)   :: importState
      logical            :: timeCheck1,timeCheck2
      type(ESMF_Field),pointer  :: fieldList(:)
      character(len=128) :: fldname

      ! query the Component for its clock
      call ESMF_GridCompGet(ATM_GRID_COMP, clock=clock, &
         importState=importState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1025, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      ! get the current time out of the clock
      call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1032, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      ! set up invalid time (by convention)
      call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, &
        h=00, m=00, s=00, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1040, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      nullify(fieldList)
      call NUOPC_GetStateMemberLists(importState, &
        fieldList=fieldList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=1048, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      !--------------------------------
      ! set the importFieldsValid flag
      ! have to be VERY CAREFUL here, the order of fieldList does NOT
      ! match the order of importFieldsList.
      ! associated(fieldList) will be false if there are no fields
      !--------------------------------

      importFieldsValid(:) = .true.
      if (associated(fieldList)) then
      do n = 1,size(fieldList)
        call ESMF_FieldGet(fieldList(n), name=fldname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1064, &
          file="module_ATM_GRID_COMP.F90")) &
          return  ! bail out
        nf = QueryFieldList(ImportFieldsList,fldname)
        timeCheck1 = NUOPC_IsAtTime(fieldList(n), invalidTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=1070, &
          file="module_ATM_GRID_COMP.F90")) &
          return  ! bail out
        if (timeCheck1) then
          importFieldsValid(nf) = .false.
        else
          timeCheck2 = NUOPC_IsAtTime(fieldList(n), currTime, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=1078, &
            file="module_ATM_GRID_COMP.F90")) &
            return  ! bail out
          if (.not.timeCheck2) then
            !TODO: introduce and use INCOMPATIBILITY return codes!!!!
            call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
              msg="NUOPC INCOMPATIBILITY DETECTED: "//&
              "Import Field not at current time", &
              line=1086, file="module_ATM_GRID_COMP.F90", &
              rcToReturn=rc)
              return  ! bail out
          endif
        endif
        write(MESSAGE_CHECK,'(A,2i4,l3)') &
          "ATM_CHECKIMPORT "//trim(fldname),n,nf,importFieldsValid(nf)
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
      enddo
      endif

!-----------------------------------------------------------------------

      END SUBROUTINE ATM_CHECKIMPORT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!

      SUBROUTINE ATM_FINALIZE(ATM_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_EARTH                               &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  Finalize the ATM component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: ATM_GRID_COMP                   !<-- The ATM component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The ATM import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The ATM import state
      TYPE(ESMF_Clock)                  :: CLOCK_EARTH                     !<-- The Clock of the EARTH component
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE                     !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Finalize step of the CORE component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%CORE_GRID_COMP &
                                ,importState=atm_int_state%CORE_IMP_STATE &
                                ,exportState=atm_int_state%CORE_EXP_STATE &
                                ,clock      =atm_int_state%CLOCK_ATM      &
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      call ESMF_ClockDestroy(atm_int_state%CLOCK_ATM, rc=RC_FINALIZE)
      if (ESMF_LogFoundError(rcToCheck=RC_FINALIZE, msg=ESMF_LOGERR_PASSTHRU, &
        line=1155, &
        file="module_ATM_GRID_COMP.F90")) &
        return  ! bail out

      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_FINALIZE succeeded'
      ELSE
        WRITE(0,*)' ATM_FINALIZE failed  RC_FINALIZE=',RC_FINALIZE
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE module_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
