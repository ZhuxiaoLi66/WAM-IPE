#include "./ESMFVersionDefine.h"

module module_MEDSpaceWeather

  use ESMF
  use ESMF_IO_NCPutGetMod

  use NUOPC
  use NUOPC_Mediator, &
    mediator_routine_SS             => SetServices, &
    mediator_routine_Run            => routine_Run, &
    mediator_label_DataInitialize   => label_DataInitialize, &
    mediator_label_Advance          => label_Advance, &
    mediator_label_CheckImport      => label_CheckImport, &
    mediator_label_TimestampExport  => label_TimestampExport, &
    mediator_label_SetRunClock      => label_SetRunClock, &
    mediator_label_Finalize         => label_Finalize

#define ESMF_NETCDF
#ifdef ESMF_NETCDF
  use netcdf
#endif

  implicit none
  
  private

  include "mpif.h"

!#define USE_CART3D_COORDSYS
!#define OUT_WEIGHT

  integer, parameter :: MAXNAMELEN = 128

  ! private internal state to keep instance data
  type InternalStateStruct
    integer :: wamdims(3)
    real(ESMF_KIND_R8), pointer :: wamhgt(:)
    integer :: myrows, startlevel, totallevels
    integer :: wamtotalnodes, localnodes, startnode
    type(ESMF_Mesh):: wam2dMesh, wamMesh, ipeMesh
    type(ESMF_RouteHandle) :: routehandle
    integer :: PetNo, PetCnt

    ! Fields to save WAM Cardinal directions
    type(ESMF_Field) :: wam_north, wam_east

    ! Unit vectors for cardinal direction to 3D Cartesian vector transformation
    type(ESMF_Field) :: wam_north_uvec, wam_east_uvec

    ! WAM wind vector expressed in 3D cartesian
    type(ESMF_Field) :: wam_wind_3Dcart_vec


    ! Unit vectors for 3D Cartesian vector to cardinal direction transformation
    type(ESMF_Field) :: ipe_north_uvec, ipe_east_uvec

    ! IPE wind vector expressed in 3D cartesian
    type(ESMF_Field) :: ipe_wind_3Dcart_vec
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type

  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    ! local variables
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(mediator, mediator_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv03
    call ESMF_GridCompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p4"/), userRoutine=InitializeP4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeP5, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)

     call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_DataInitialize, &
       specRoutine=DataInitialize, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_Advance, &
      specRoutine=MediatorAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_Finalize, &
      specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(mediator, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: mediator
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(mediator, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  subroutine InitializeAdvertise(mediator, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS

    ! exportable fields: WAM export fields
    call NUOPC_Advertise(importState, StandardNames=(/ &
      "northward_wind_neutral         ", &
      "eastward_wind_neutral          ", &
      "upward_wind_neutral            ", &
      "temp_neutral                   ", &
      "O_Density                      ", &
      "O2_Density                     ", &
      "N2_Density                     ", &
      "height                         "  &
      /), transferOfferGeomObject = "will provide", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! exportable fields: IPE import fields
    call NUOPC_Advertise(exportState, StandardNames=(/ &
      "northward_wind_neutral         ", &
      "eastward_wind_neutral          ", &
      "upward_wind_neutral            ", &
      "temp_neutral                   ", &
      "O_Density                      ", &
      "O2_Density                     ", &
      "N2_Density                     "  &
      /), TransferOfferGeomObject="cannot provide", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

        WRITE(0,*)' MEDIATOR complete advertizing fields '

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(mediator, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    type(ESMF_Mesh):: wam2dMesh

    rc = ESMF_SUCCESS

    ! calling into Peggy's code
    call initGrids(mediator, wam2dMesh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! create fields
    call realizeConnectedFields(importState, wam2dMesh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

        WRITE(0,*)' MEDIATOR complete realized fields '


  contains

    subroutine realizeConnectedFields(state, mesh, rc)
      ! TODO: this method may move into the NUOPC_ utility layer
      type(ESMF_State)                :: state
      type(ESMF_Mesh)                 :: mesh
      integer, intent(out), optional  :: rc

      ! local variables
      character(len=80), allocatable  :: fieldNameList(:)
      integer                         :: i, itemCount, k
      type(ESMF_ArraySpec)            :: arrayspec
      type(ESMF_Field)                :: field
      integer                         :: levels

      if (present(rc)) rc = ESMF_SUCCESS

      levels = 150 ! hardcode level, probably should get it from internal state

      call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList(itemCount))
      call ESMF_StateGet(state, itemNameList=fieldNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      k=1 ! initialize
      do i=1, itemCount
          ! create a Field with one undistributed dimension
          call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          field = ESMF_FieldCreate(mesh, arrayspec, &
	   ungriddedLBound=(/1/), ungriddedUBound=(/levels/), &
	   name=fieldNameList(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          ! realize the connected Field using the just created Field
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      enddo

    end subroutine realizeConnectedFields

  end subroutine initializeRealize
    
  subroutine InitializeP4(mediator, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    type(ESMF_VM)        :: vm
    type(ESMF_Mesh)      :: ipemesh, medmesh
    type(ESMF_Field)     :: field
    character(len=80), allocatable  :: fieldNameList(:)
    integer              :: i, itemCount, k
    type(ESMF_DistGrid)  :: ipedistgrid, meddistgrid
    integer              :: minIndices(1,1), maxIndices(1,1)
    integer              :: decount, PetNo, PetCnt
    type(InternalState)  :: is
    logical              :: freeflag

    rc = ESMF_SUCCESS

        WRITE(0,*)' MEDIATOR starts InitializeP4 '

  ! query component for its internal state
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(mediator, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    PetNo = is%wrap%PetNo
    PetCnt = is%wrap%PetCnt

    ! Get the IPE field from the IPE module, get the elemDistgrid and redistribute in over
    ! the number of processors used by the Mediator
    
      call ESMF_StateGet(exportState, itemCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList(itemCount))
      call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! Get the mesh from the first itme
      call ESMF_StateGet(exportState, field=field, itemName=fieldNameList(1), rc=rc)   
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, mesh=ipemesh, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_MeshGet(ipemesh, nodalDistgrid=ipedistgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_DistGridGet(ipedistgrid, deCount=decount, &
      	   minIndexPTile=minIndices, maxIndexPTile=maxIndices, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! If the original distribution (decount) use the same number of 
      ! PETs (PetCnt), not need to redistribute.  Otherwise, redistribute
      ! the elementDistgrid to use PetCnt processors
      if (PetCnt /= decount) then
         ! create a new distgrid evenly distribute the nodes over all the PEs.
         meddistgrid = ESMF_DistGridCreate((/minIndices(1,1)/), (/maxIndices(1,1)/), &
	 	       rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
            	file=__FILE__)) &
            	return  ! bail out
         medmesh = ESMF_MeshCreate(meddistgrid, meddistgrid, &
             	  rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out
         call ESMF_MeshDestroy(ipemesh, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out
        ! replace the field with new mesh     
        k=1 ! initialize
        do i=1, itemCount
          ! create a Field with one undistributed dimension
          call ESMF_StateGet(exportState, field=field, itemName=fieldNameList(i), rc=rc)   
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
        	file=__FILE__)) &
		return  ! bail out
            		
          call ESMF_FieldEmptySet(field, medmesh, meshloc=ESMF_MESHLOC_NODE,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
        	file=__FILE__)) &
		return  ! bail out
        enddo
     endif

     WRITE(0,*)' MEDIATOR done initializeP4'
      
  end subroutine initializeP4

  !----------------------------------------------------------------------------

  subroutine InitializeP5(mediator, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_VM)        :: vm
    character(len=80), allocatable  :: fieldNameList(:)
    integer              :: i, itemCount, k
    type(ESMF_Field)     :: field, wamfield, ipefield
    type(ESMF_ArraySpec) :: arrayspec
    real(ESMF_KIND_R8)   :: starttime, endtime, timesend(1), timereport(1)
    type(InternalState)  :: is
    real(ESMF_KIND_R8), pointer :: fptr(:)

    

    rc = ESMF_SUCCESS

     WRITE(0,*)' MEDIATOR start initializeP5'

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  ! query component for its internal state
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(mediator, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(exportState, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    allocate(fieldNameList(itemCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    k=1 ! initialize
    do i=1, itemCount
       ! create a Field with one undistributed dimension
       call ESMF_StateGet(exportState, field=field, itemName=fieldNameList(i), rc=rc)   
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
        	file=__FILE__)) &
		return  ! bail out
            		
       call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
        	file=__FILE__)) &
		return  ! bail out
       
      ! Intialize the fields to -999.0 (null value)
      call ESMF_FieldGet(field, farrayptr=fptr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
        	file=__FILE__)) &
		return  ! bail out
      fptr = -999.0
 
      if (i==1) then
         call ESMF_FieldGet(field, mesh=is%wrap%ipeMesh, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, &
           file=__FILE__)) &
           return  ! bail out
      endif
    enddo


    ! Create Fields to hold WAM cardinal vectors in preparation for transformation
    is%wrap%wam_north=ESMF_FieldCreate(is%wrap%wamMesh,typekind=ESMF_TYPEKIND_R8, &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    is%wrap%wam_east=ESMF_FieldCreate(is%wrap%wamMesh,typekind=ESMF_TYPEKIND_R8, &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    ! Create Fields to hold WAM unit vectors for vector transformation
    is%wrap%wam_north_uvec=ESMF_FieldCreate(is%wrap%wamMesh,typekind=ESMF_TYPEKIND_R8, &
         gridToFieldMap=(/2/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    is%wrap%wam_east_uvec=ESMF_FieldCreate(is%wrap%wamMesh,typekind=ESMF_TYPEKIND_R8, &
         gridToFieldMap=(/2/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    ! Set 3D Cartesian unit vectors for IPE
    call Set_Field_Cardinal_UVecs(is%wrap%wam_north_uvec, is%wrap%wam_east_uvec, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out


    ! Create Fields to hold ipe unit vectors for vector transformation
    is%wrap%ipe_north_uvec=ESMF_FieldCreate(is%wrap%ipeMesh,typekind=ESMF_TYPEKIND_R8, &
         gridToFieldMap=(/2/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    is%wrap%ipe_east_uvec=ESMF_FieldCreate(is%wrap%ipeMesh,typekind=ESMF_TYPEKIND_R8, &
         gridToFieldMap=(/2/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    ! Set 3D Cartesian unit vectors for IPE
    ! USE CART VERSION BECAUSE OF A BUG IN MESHREDIST() IN ESMF 7.0.0 and before, once in 7.1.0 USE NON-CART VERSION
    if (ESMF_VERSION_MAJOR > 7) then
       call Set_Field_Cardinal_UVecs(is%wrap%ipe_north_uvec, is%wrap%ipe_east_uvec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

    else if ((ESMF_VERSION_MAJOR == 7) .and. (ESMF_VERSION_MINOR > 0)) then
       call Set_Field_Cardinal_UVecs(is%wrap%ipe_north_uvec, is%wrap%ipe_east_uvec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    else
       call Set_Field_Cardinal_UVecsCart(is%wrap%ipe_north_uvec, is%wrap%ipe_east_uvec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    endif

    ! Create Fields to hold Cartesian versions of wind
    is%wrap%wam_wind_3Dcart_vec=ESMF_FieldCreate(is%wrap%wamMesh,typekind=ESMF_TYPEKIND_R8, &
         gridToFieldMap=(/2/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    is%wrap%ipe_wind_3Dcart_vec=ESMF_FieldCreate(is%wrap%ipeMesh,typekind=ESMF_TYPEKIND_R8, &
         gridToFieldMap=(/2/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

      ! Now we have both wammesh and ipemesh, call ESMF_FieldRegridStore() to create
      ! a routehandle
      ! Create src and dst fields and run RegridStore()
      call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out

      wamField = ESMF_FieldCreate(is%wrap%wamMesh, arrayspec, meshloc=ESMF_MESHLOC_NODE,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out
      ipeField = ESMF_FieldCreate(is%wrap%ipeMesh, arrayspec, meshloc=ESMF_MESHLOC_NODE,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out


#if 1
      call ESMF_VMBarrier(vm)
      call ESMF_VMWTime(starttime)
#endif
      call ESMF_FieldRegridStore(wamField, ipeField, &
       	 unmappedaction =ESMF_UNMAPPEDACTION_IGNORE, &
	 regridmethod = ESMF_REGRIDMETHOD_BILINEAR, &
	 polemethod = ESMF_POLEMETHOD_NONE, &
         lineType = ESMF_LINETYPE_GREAT_CIRCLE, &
	 routehandle = is%wrap%routehandle, &
	 rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out

#if 1
      call ESMF_VMWtime(endtime)
      timesend(1)=endtime-starttime
      call ESMF_VMReduce(vm, sendData=timesend, recvData=timereport, count=1, &
    	 reduceflag=ESMF_REDUCE_MAX, rootPet=0, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            	line=__LINE__, &
        	file=__FILE__)) &
        	return  ! bail out
      if (is%wrap%PetNo==0) then
         print *, 'Time to do RegridStore WAM->IPE is ', timereport(1)*1000, 'msec'
      endif
#endif
      call ESMF_FieldDestroy(wamfield)
      call ESMF_FieldDestroy(ipefield)

     WRITE(0,*)' MEDIATOR done initializeP5'

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(mediator, rc)
     type(ESMF_GridComp)  :: mediator
     integer, intent(out) :: rc

     ! indicate that data initialization is complete (breaking out of init-loop)
     call NUOPC_CompAttributeSet(mediator, &
       name="InitializeDataComplete", value="true", rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

   end subroutine

  !-----------------------------------------------------------------------------

  subroutine MediatorAdvance(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(mediator, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! HERE THE MEDIATOR does the mediation of Fields that come in on the
    ! importState with a timestamp consistent to the currTime of the 
    ! mediators Clock.
    
    ! The Mediator uses the data on the import Fields to update the data
    ! held by Fields in the exportState.
    
    ! After this routine returns the generic Mediator will correctly
    ! timestamp the export Fields and update the Mediator Clock to:
    !
    !       currTime -> currTime + timeStep
    !
    ! Where the timeStep is equal to the parent timeStep.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="-------->MED Advance() mediating for: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! calling into Peggy's code
    call RunRegrid(mediator, importstate, exportstate, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="----------------> model time step to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
     
  end subroutine

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Peggy's routines below....
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  subroutine initGrids(model, wam2dMesh, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_MESH) :: wam2dmesh
    integer, intent(out) :: rc
    

  type(ESMF_VM) :: vm
  type(InternalState)     :: is
  type(ESMF_MESH) :: wammesh
  real(ESMF_KIND_R8), pointer :: wamlon(:,:), wamlat(:,:), wamhgt(:)
  real(ESMF_KIND_R8), pointer :: ipelon(:,:), ipelat(:,:), ipehgt(:), ipedata(:)
  real(ESMF_KIND_R8), pointer :: hgtbuf(:,:,:), varbuf(:,:,:,:)
  real(ESMF_KIND_R8), pointer :: databuf(:)
  integer(ESMF_KIND_I4), pointer :: maxlevs(:)
  integer(ESMF_KIND_I4), pointer :: numPerRow(:), shuffleOrder(:)
  integer :: nc1, ncid
  integer :: varid, dimid1, dimid2, dimid3
  integer :: ndims, dimids(3)
  integer :: wamdims(3), ipedims(3)
  integer :: PetNo, PetCnt
  integer(ESMF_KIND_I4), pointer :: elementIds(:), elementTypes(:), elementConn(:)
  integer(ESMF_KIND_I4), pointer :: nodeIds(:), nodeOwners(:)
  real(ESMF_KIND_R8), pointer :: nodeCoords(:)
  integer(ESMF_KIND_I4), pointer :: southind(:), northind(:), totallats(:), gap(:)
  integer :: minheight, maxheight, halo, neighbor, remind
  integer(ESMF_KIND_I4), pointer :: totalheight(:)
  real(ESMF_KIND_R8) :: lon, lat, hgt, lon1, lat1, hgt1
  real(ESMF_KIND_R4) :: interval
  integer :: i,j, k, l, ii, jj, kk, count1, count3, count8, localcount, countup, save, base, base1
  logical :: even
  integer :: start, count, diff, lastlat, totalelements, totalnodes, localnodes, startid
  integer :: wamtotalnodes
  integer :: elmtcount, increment
  integer :: startlevel, next, ind, ind1, totalnodes2d, totallevels, myrows, trigs
  integer :: reminder, steps, startrow, startnode
  integer, pointer :: rowinds(:), petTable(:), baseind(:)
  integer(ESMF_KIND_I4), pointer :: elementCnt(:), nodeCnt(:), sendbuf(:), recvbuf(:)
  integer(ESMF_KIND_I4), allocatable :: indList(:)
  real(ESMF_KIND_R8), pointer :: conntbl(:), globalCoords(:,:), fptr2d(:,:), fptr1d(:)
  type(ESMF_Arrayspec) :: arrayspec
  type(ESMF_Array) :: array, array1, array2
  real(ESMF_KIND_R8) :: maxerror, minerror, totalerrors, deg2rad
  real(ESMF_KIND_R8) :: starttime, endtime, timesend(1), timereport(1)
  real(ESMF_KIND_R8) :: differr
  real(ESMF_KIND_R8), pointer :: varout(:), lonbuf(:),latbuf(:)
  real(ESMF_KIND_R8), pointer :: weights(:)
  integer(ESMF_KIND_I4), pointer :: indices(:,:)
  character(len=MAXNAMELEN) :: wamfilename
  character(len=MAXNAMELEN) :: filename
  integer :: wgtcount(1)
  integer :: count2(1), start2(1)
  integer, pointer :: allCounts(:), connectbase(:)
  real, parameter :: PI=3.1415927
  integer :: j1
  integer :: localrc, status
  real, parameter :: earthradius=6371.0  !in kilometers

  ! For output
  real(ESMF_KIND_R8), pointer :: lontbl(:), lattbl(:), hgttbl(:)
  integer :: lonid, latid, hgtid, vertid, elmtid, nodeid, numid, connid, timeid
  integer :: data1id, data2id, wgtid, wamid, ipeid
  integer :: globalTotal, globalTotalelmt, nodestartid, totalwgts
  type(ESMF_Distgrid) :: nodalDistgrid, distgrid

  rc = ESMF_SUCCESS

  !------------------------------------------------------------------------
  ! get global vm information
  !
  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  ! set up local pet info
  call ESMF_VMGet(vm, localPet=PetNo, petCount=PetCnt, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  !------------------------------------------------------------------------
  ! Allocate memory for the internal state and set it in the Component.
    allocate(is%wrap, stat=rc)
    if (ESMF_LogFoundAllocError(statusToCheck=rc, &
      msg="Allocation of the internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  !wamfilename = 'data/wam3dgridnew.nc'
  wamfilename = 'wam3dgridnew_20160427.nc'
  filename = 'wam2dmesh.nc'
  minheight = 90
  maxheight = 800
  deg2rad = PI/180.0

#ifdef ESMF_NETCDF
  !!-------------------------------------
  !! Create WAM mesh
  !!-------------------------------------
  ! We need to create a 2D mesh with distgrid only and use it to get the data
  ! from the DATAWAM
  ! we also need to create the 3D intermediate WAM mesh to be used to regrid
  ! with IPE grid
 
  ! Read in WAM grid from wam3dgrid.nc
  !!
  status = nf90_open(path= wamfilename, mode=nf90_nowrite, ncid=nc1)
  call CheckNCError(status, wamfilename)
  status = nf90_inq_varid(nc1,'lons', varid)
  call CheckNCError(status, 'lons')
  status = nf90_inquire_variable(nc1, varid, ndims=ndims, dimids = dimids)
  call CheckNCError(status, 'lons')
  status = nf90_inquire_dimension(nc1,dimids(1), len=wamdims(1))
  call CheckNCError(status, 'lons 1st dimension')
  status = nf90_inquire_dimension(nc1,dimids(2), len=wamdims(2))
  call CheckNCError(status, 'lons 2nd dimension')

  ! WAM dimension order:  lons, lats (192, 94)
  allocate(wamlon(wamdims(1), wamdims(2)), &
  	   wamlat(wamdims(1), wamdims(2)))
  status = nf90_get_var(nc1, varid, wamlon)
  call CheckNCError(status, 'lons')
  status = nf90_inq_varid(nc1,'lats', varid)
  call CheckNCError(status, 'lats')
  status = nf90_get_var(nc1, varid, wamlat)
  call CheckNCError(status, 'lats')

  ! intermediate height fields
  status = nf90_inq_varid(nc1,'height', varid)
  call CheckNCError(status, 'height')
  status = nf90_inquire_variable(nc1, varid, ndims=ndims, dimids = dimids)
  call CheckNCError(status, 'height')
  status = nf90_inquire_dimension(nc1,dimids(1), len=wamdims(3))
  call CheckNCError(status, 'height 1st dimension')

  allocate(wamhgt(wamdims(3)))
  status = nf90_get_var(nc1, varid, wamhgt)
  call CheckNCError(status, 'height')
  allocate(NumPerRow(wamdims(2)), ShuffleOrder(wamdims(2)))
  status = nf90_inq_varid(nc1,'NumPerRow', varid)
  call CheckNCError(status, 'NumPerRow')
  status = nf90_get_var(nc1, varid, NumPerRow)
  call CheckNCError(status, 'NumPerRow')
  status = nf90_inq_varid(nc1,'ShuffleOrder', varid)
  call CheckNCError(status, 'ShuffleOrder')
  status = nf90_get_var(nc1, varid, ShuffleOrder)
  call CheckNCError(status, 'ShuffleOrder')
  status= nf90_close(nc1)
  call CheckNCError(status, wamfilename)

  ! Use the shuffle order to create the 2D mesh first
  ! find the total number of nodes in each processor and create local index table
  localnodes=0
  myrows = wamdims(2)/PetCnt
#if 0
  if ((wamdims(2)-myrows*PetCnt) > PetNo) myrows = myrows+1
  allocate(rowinds(myrows))      !my local row index
  allocate(petTable(wamdims(2))) !the owner PET for each row
  next = 0
  ind1 = 0
  do i=1,wamdims(2)
    ind=ShuffleOrder(i)
    petTable(ind)=next
    if (next == PetNo) then
       ind1=ind1+1
       rowinds(ind1)=ind
       localnodes=localnodes + numPerRow(ind)
    endif
    next=next+1
    if (next == PetCnt) next=0
  enddo
  print *, PetNo, 'localnodes/myrows:', localnodes, myrows
#else
  reminder = wamdims(2)-myrows*PetCnt
  if (reminder > PetNo)  then
     myrows = myrows+1
     startrow = myrows*PetNo+1
  else
     startrow = myrows*PetNo + reminder+1
  endif
  allocate(rowinds(myrows))      !my local row index
  allocate(petTable(wamdims(2))) !the owner PET for each row
  do i=1,myrows
    rowinds(i)=ShuffleOrder(startrow+i-1)
    localnodes = localnodes + numPerRow(rowinds(i))
  enddo
  print *, PetNo, 'start rows:', startrow, localnodes, myrows
  ind1 = 1
  steps = wamdims(2)/PetCnt
  do next=0,PetCnt-1
    if (next<reminder) then
       do i=ind1, ind1+steps
          PetTable(ShuffleOrder(i))=next
       enddo
       ind1 = ind1+steps+1
    else
       do i=ind1, ind1+steps-1
          PetTable(ShuffleOrder(i))=next
       enddo
       ind1 = ind1+steps
    endif
  enddo
  if (ind1 /= wamdims(2)+1) then
     print *, 'MEDIATOR Wrong ind1 ', ind1, wamdims(2)
  endif
#endif
  
  ! sort rowinds
  call ESMF_UtilSort(rowinds, ESMF_SORTFLAG_ASCENDING, rc)

  print *, PetNo, 'sorted rowinds ', rowinds

  ! Save the lat/lon in the order the distgrid
  allocate(lonbuf(localnodes), latbuf(localnodes))

  ! Create a distgrid using a collapsed 1D index array based on the local row index
  allocate(indList(localnodes))
  k=1
  do i=1,myrows
    ind=rowinds(i)
    do j=1,numPerRow(ind)
       indList(k)=(ind-1)*wamdims(1)+j
       lonbuf(k)=wamlon(j,ind)
       latbuf(k)=wamlat(j,ind)
       k=k+1
    enddo
  enddo

  ! write out the index into a NetCDF file
  ! find the total nodes first
  totalnodes = sum(numPerRow)
  ! find out starting node id
  startnode = 1
  do i=1,startrow-1
    startnode = startnode+numPerRow(ShuffleOrder(i))
  enddo
  call ESMF_VMBarrier(vm)

  distgrid = ESMF_DistGridCreate(indList, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
  ! Create mesh using the distgrid as the nodaldistgrid,  no elemdistgrid available
  ! just use nodeldistgrid for both
  wam2dmesh = ESMF_MeshCreate(distgrid,distgrid,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
  
  ! Create the 3D mesh with fixed height
  ! find the lower height level where height > minheight
  do i=1,wamdims(3)
     if (wamhgt(i) > minheight) then
    	startlevel = i-1
	exit
     endif
  enddo
  totallevels = wamdims(3)-startlevel+1

  if (PetNo == 0) then
    ! create the output file
    status = nf90_create(filename, NF90_CLOBBER, ncid)
    call CheckNCError(status, filename)
    status = nf90_def_dim(ncid, 'nodes', totalnodes, dimid1)
    call CheckNCError(status, 'dimension nodes')
    status = nf90_def_dim(ncid, 'levels', totallevels, dimid2)
    call CheckNCError(status, 'dimension nodes')
    status = nf90_def_dim(ncid, 'vars', 7, dimid3)
    call CheckNCError(status, 'dimension nodes')
    status = nf90_def_var(ncid, 'lons', NF90_FLOAT, (/dimid1/), varid)
    call CheckNCError(status, 'variable lons')
    status = nf90_def_var(ncid, 'lats', NF90_FLOAT, (/dimid1/), varid)
    call CheckNCError(status, 'variable lats')
    status = nf90_def_var(ncid, 'wamdata', NF90_FLOAT, (/dimid1, dimid2, dimid3/), varid)
    call CheckNCError(status, 'variable wamdata')
    status = nf90_close(ncid)
    call CheckNCError(status, filename)
  endif

  ! write out the lat/lon from each processor
  do i=0, PetCnt-1
    if (PetNo == i) then
      status = nf90_open(filename, NF90_WRITE, ncid)
      call CheckNCError(status, filename)
      count2(1)=localnodes
      start2(1)=startnode
      status = nf90_inq_varId(ncid, 'lons', varid)
      call CheckNCError(status, 'lons')
      status = nf90_put_var(ncid, varid, lonbuf, & 
		    start2, count2)
      call CheckNCError(status, 'lons')
      status = nf90_inq_varId(ncid, 'lats', varid)
      call CheckNCError(status, 'lats')
      status = nf90_put_var(ncid, varid, latbuf, & 
		    start2, count2)
      call CheckNCError(status, 'lats')
      status = nf90_close(ncid)
      call CheckNCError(status, filename)
    endif
    call ESMF_VMBarrier(vm)
  enddo
  deallocate(lonbuf, latbuf)

  ! create the node table, find the total number of nodes in each processor, including the not-owned node
  totalnodes=0
  totalelements=0 
  allocate(baseind(myrows))
  do i=1,myrows
     ind=rowinds(i)
     baseind(i)=totalnodes
     ! Add the neighbor nodes
     ! If PetCnt==1, no need to add neighbor node
     ! If last row, no need to add neighbors
     ! If the neighbor is local, no need to add
     if (ind < wamdims(2)) then
       if (PetCnt>1) then
        if ((i < myrows .and. rowinds(i+1) /= ind+1) .or. i==myrows) then
          totalnodes=totalnodes+numPerRow(ind)+numPerRow(ind+1)
        else
          totalnodes=totalnodes+numPerRow(ind)
        endif
       endif
       if (numPerRow(ind) >= numPerRow(ind+1)) then
         totalelements = totalelements + numPerRow(ind)
       else
         totalelements = totalelements + numPerRow(ind+1)
       endif
     else
       totalnodes=totalnodes+numPerRow(ind)
       !Add extra elements at the top
       totalelements = totalelements+numPerRow(ind)-2
     endif

     ! Add extra elements for south pole
      if (ind==1) then 
        totalelements = totalelements+numPerRow(ind)-2
     endif
  enddo
  if (PetCnt == 1) then
     baseind(1)=0
     do i=2,wamdims(2)
       baseind(i)=baseind(i-1)+numPerRow(i-1)
     enddo
  endif

  totalnodes2d=totalnodes  ! totalnodes includes neighboring nodes and my own nodes
  totalnodes = totalnodes * totallevels 
  localnodes = localnodes * totallevels  ! localnodes are locally owned nodes
  totalelements = totalelements * (totallevels-1)
  allocate(nodeIds(totalnodes), nodeOwners(totalnodes), nodeCoords(totalnodes*3))

  ! Fill nodeIds, nodeOwners, and nodeCoords arrays, longitude first, latitude, then height
  count1=1
  localcount=1
  count3=1
  if (PetCnt > 1) then
  do k=1, totallevels
    do i=1,myrows
       ind=rowinds(i)
       do j=1,numPerRow(ind)
          ! Global id based on the 3D indices
       	  nodeIds(count1)= j+wamdims(1)*(ind-1)+wamdims(1)*wamdims(2)*(k-1)
          nodeOwners(count1)=PetNo
	  lon = wamlon(j,ind)
          lat  = wamlat(j,ind)
          hgt = wamhgt(startlevel+k-1)
          call convert2Cart(lon, lat, hgt, nodeCoords(count3:count3+2))
          count1=count1+1
          localcount=localcount+1
          count3=count3+3
       enddo
       ! if not the last row, add the neighbor row's nodes and the neighbor
       ! is not local
       if (ind < wamdims(2)) then
         if (i==myrows .or. (i < myrows .and. rowinds(i+1)/= ind+1)) then
       	  do j=1, numPerRow(ind+1)
            ! Global id based on the 3D indices
            nodeIds(count1)= j+wamdims(1)*ind+wamdims(1)*wamdims(2)*(k-1)
	    if (PetTable(ind+1) == PetNo) then
	       print *, PetNo, 'wrong neighbor ', count1, PetTable(ind+1)
            endif
            nodeOwners(count1)=PetTable(ind+1)
	    lon = wamlon(j,ind+1)
            lat  = wamlat(j,ind+1)
            hgt = wamhgt(k+startlevel-1)
            call convert2Cart(lon, lat, hgt, nodeCoords(count3:count3+2))
            count1=count1+1
            count3=count3+3
          enddo
         endif
	endif
     enddo
  enddo
  else ! PetCnt==1
  ! For sequential case, store the rows in its order, do not shuffle
  do k=1, totallevels
    do ind=1,myrows
       do j=1,numPerRow(ind)
          ! Global id based on the 3D indices
       	  nodeIds(count1)= j+wamdims(1)*(ind-1)+wamdims(1)*wamdims(2)*(k-1)
          nodeOwners(count1)=PetNo
	  lon = wamlon(j,ind)
          lat  = wamlat(j,ind)
          hgt = wamhgt(k+startlevel-1)
          call convert2Cart(lon, lat, hgt, nodeCoords(count3:count3+2))
          count1=count1+1
          localcount=localcount+1
          count3=count3+3
       enddo
     enddo
  enddo
  endif ! PetCnt > 1

  if (count1-1 /= totalnodes .or. localcount-1 /= localnodes) then
     print *, 'totalcount mismatch ', count1-1, totalnodes, localcount-1, localnodes
  endif

#ifdef USE_CART3D_COORDSYS
  wamMesh = ESMF_MeshCreate(3,3,coordSys=ESMF_COORDSYS_CART, rc=rc)
#else
  wamMesh = ESMF_MeshCreate(3,3,coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
#endif
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
  call ESMF_MeshAddNodes(wamMesh, nodeIds, nodeCoords, nodeOwners, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
  
  deallocate(wamlon, wamlat)
  deallocate(nodeIds, nodeCoords, nodeOwners)

  allocate(elementIds(totalelements), elementTypes(totalelements), &
           elementConn(totalelements*8))

  elementTypes(:)=ESMF_MESHELEMTYPE_HEX

  ! find out the starting global id of the local element
  allocate(elementCnt(PetCnt),sendbuf(1))
  sendbuf(1)=totalelements
  call ESMF_VMAllGather(vm, sendbuf, elementCnt, 1, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
  
  ! find the starting elementID
  startid=0
  do i=1,PetNo
    startid=startid+elementCnt(i)
  enddo
  globaltotalelmt = 0
  do i=1,PetCnt
    globaltotalelmt=globaltotalelmt+elementCnt(i)
  enddo
  deallocate(elementCnt, sendbuf)

  ! Build the local elementConn table using local node indices
  ! If PetCnt=1, no shuffle, the rows are in order
  count1=1
  count8=1
  do k=1, totallevels-1
    do i=1,myrows
       if (PetCnt > 1) then
           ind=rowinds(i)
       else
	   ind = i
       endif
       base = baseind(i)+totalnodes2d*(k-1)
       if (ind == wamdims(2)) then
#if 0
         ! create dummy elements by connecting every other nodes to form a triangle to cover the pole
         do j=1, numPerRow(ind)-2, 2
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+j+1
	   elementConn(count8+2)=base+j+2
	   elementConn(count8+3)=base+j+2
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+j+1
	   elementConn(count8+6)=base+totalnodes2d+j+2
	   elementConn(count8+7)=base+totalnodes2d+j+2
	   count1=count1+1
	   count8=count8+8
         enddo
         ! Last one, connect it back to the first node 
     	 elementIds(count1)=startid+count1
	 elementConn(count8)= base+j
	 elementConn(count8+1)=base+j+1
	 elementConn(count8+2)=base+1
	 elementConn(count8+3)=base+1
	 elementConn(count8+4)=base+totalnodes2d+j
	 elementConn(count8+5)=base+totalnodes2d+j+1
	 elementConn(count8+6)=base+totalnodes2d+1
	 elementConn(count8+7)=base+totalnodes2d+1
	 count1=count1+1
	 count8=count8+8
         cycle       
#else
         ! using the zigzag method to create triangles that covers the pole
         ! First half
         do j=1, numPerRow(ind)/2-1
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+j+1
	   elementConn(count8+2)=base+numPerRow(ind)-j
	   elementConn(count8+3)=base+numPerRow(ind)-j
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+j+1
	   elementConn(count8+6)=base+totalnodes2d+numPerRow(ind)-j
	   elementConn(count8+7)=base+totalnodes2d+numPerRow(ind)-j
	   count1=count1+1
	   count8=count8+8
         enddo
         ! second half
         do j=numPerRow(ind)/2+1, numPerRow(ind)-1
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+j+1
	   elementConn(count8+2)=base+numPerRow(ind)-j
	   elementConn(count8+3)=base+numPerRow(ind)-j
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+j+1
	   elementConn(count8+6)=base+totalnodes2d+numPerRow(ind)-j
	   elementConn(count8+7)=base+totalnodes2d+numPerRow(ind)-j
	   count1=count1+1
	   count8=count8+8
         enddo	   
         cycle
#endif
       endif

      ! South pole
      if (ind == 1) then
         ! using the zigzag method to create triangles that covers the pole
         ! First half
         do j=1, numPerRow(ind)/2-1
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+j+1
	   elementConn(count8+2)=base+numPerRow(ind)-j
	   elementConn(count8+3)=base+numPerRow(ind)-j
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+j+1
	   elementConn(count8+6)=base+totalnodes2d+numPerRow(ind)-j
	   elementConn(count8+7)=base+totalnodes2d+numPerRow(ind)-j
	   count1=count1+1
	   count8=count8+8
         enddo

	 ! Second half
         do j=numPerRow(ind)/2+1, numPerRow(ind)-1
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+j+1
	   elementConn(count8+2)=base+numPerRow(ind)-j
	   elementConn(count8+3)=base+numPerRow(ind)-j
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+j+1
	   elementConn(count8+6)=base+totalnodes2d+numPerRow(ind)-j
	   elementConn(count8+7)=base+totalnodes2d+numPerRow(ind)-j
	   count1=count1+1
	   count8=count8+8
         enddo	   
       endif

       ! the two adjacent rows have the same number of points, elements are cubes
       if (numPerRow(ind+1) == numPerRow(ind)) then
         do j=1,numPerRow(ind)-1
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+j+1
	   elementConn(count8+2)=base+numPerRow(ind)+j+1
	   elementConn(count8+3)=base+numPerRow(ind)+j
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+j+1
	   elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+j+1
	   elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+j
	   count1=count1+1
	   count8=count8+8
         enddo
         ! last one in the row, wrap around
       	 elementIds(count1)=startid+count1
	 elementConn(count8)= base+j
	 elementConn(count8+1)=base+1
	 elementConn(count8+2)=base+numPerRow(ind)+1
	 elementConn(count8+3)=base+numPerRow(ind)+j
	 elementConn(count8+4)=base+totalnodes2d+j
	 elementConn(count8+5)=base+totalnodes2d+1
	 elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+1
	 elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+j
	 count1=count1+1
	 count8=count8+8
       else
         ! the number of nodes are different, make prism elements
	 diff=numPerRow(ind)-numPerRow(ind+1)
	 if (diff > 0) then 
          ! make triangles with base at lower row
          ! triangles will be evenly distributed
	  interval = real(numPerRow(ind))/(diff+1)
	  jj=1
          trigs=1
          do j=1,numPerRow(ind)-1
 	     if (j > trigs*interval) then
!	     if (mod(j,increment)==0) then
               ! triangles - base at bottom
               trigs=trigs+1
               elementIds(count1)=startid+count1
	       elementConn(count8)= base+j
	       elementConn(count8+1)=base+j+1
	       elementConn(count8+2)=base+numPerRow(ind)+jj
	       elementConn(count8+3)=base+numPerRow(ind)+jj
	       elementConn(count8+4)=base+totalnodes2d+j
	       elementConn(count8+5)=base+totalnodes2d+j+1
	       elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+jj
	       elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+jj
             else
               elementIds(count1)=startid+count1
	       elementConn(count8)= base+j
	       elementConn(count8+1)=base+j+1
	       elementConn(count8+2)=base+numPerRow(ind)+jj+1
	       elementConn(count8+3)=base+numPerRow(ind)+jj
	       elementConn(count8+4)=base+totalnodes2d+j
	       elementConn(count8+5)=base+totalnodes2d+j+1
	       elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+jj+1
	       elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+jj
	       jj=jj+1
	     endif
	     count1=count1+1
	     count8=count8+8
           enddo
           ! last one in the row, wrap around
           elementIds(count1)=startid+count1
	   elementConn(count8)= base+j
	   elementConn(count8+1)=base+1
	   elementConn(count8+2)=base+numPerRow(ind)+1
	   elementConn(count8+3)=base+numPerRow(ind)+jj
	   elementConn(count8+4)=base+totalnodes2d+j
	   elementConn(count8+5)=base+totalnodes2d+1
	   elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+1
	   elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+jj
   	   count1=count1+1
	   count8=count8+8
	   if (k==1 .and. (jj /= numPerRow(ind+1))) then
	      print *, PetNo, 'Upper row index mismatch', ind, jj, numPerRow(ind+1)
           endif
        else  ! diff < 0 
          ! make triangles with base at upper row
          ! triangles will be evenly distributed
	  interval = real(numPerRow(ind+1))/(-1*diff+1)
	  jj=1
	  trigs=1
          do j=1,numPerRow(ind+1)-1
	     if (j > trigs*interval) then
	       trigs = trigs+1              
               ! triangles - base at bottom
               elementIds(count1)=startid+count1
	       elementConn(count8)= base+jj
	       elementConn(count8+1)=base+jj
	       elementConn(count8+2)=base+numPerRow(ind)+j+1
	       elementConn(count8+3)=base+numPerRow(ind)+j
	       elementConn(count8+4)=base+totalnodes2d+jj
	       elementConn(count8+5)=base+totalnodes2d+jj
	       elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+j+1
	       elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+j
             else
               elementIds(count1)=startid+count1
	       elementConn(count8)= base+jj
	       elementConn(count8+1)=base+jj+1
	       elementConn(count8+2)=base+numPerRow(ind)+j+1
	       elementConn(count8+3)=base+numPerRow(ind)+j
	       elementConn(count8+4)=base+totalnodes2d+jj
	       elementConn(count8+5)=base+totalnodes2d+jj+1
	       elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+j+1
	       elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+j
	       jj=jj+1
	     endif
	     count1=count1+1
 	     count8=count8+8
           enddo
           ! last one in the row, wrap around
       	   elementIds(count1)=startid+count1
	   elementConn(count8)= base+jj
	   elementConn(count8+1)=base+1
	   elementConn(count8+2)=base+numPerRow(ind)+1
	   elementConn(count8+3)=base+numPerRow(ind)+j
	   elementConn(count8+4)=base+totalnodes2d+jj
	   elementConn(count8+5)=base+totalnodes2d+1
	   elementConn(count8+6)=base+numPerRow(ind)+totalnodes2d+1
	   elementConn(count8+7)=base+numPerRow(ind)+totalnodes2d+j
 	   count1=count1+1
	   count8=count8+8
	   if (k==1 .and. (jj /= numPerRow(ind))) then
	      print *, PetNo, 'Lower row index mismatch', ind, jj, numPerRow(ind)
           endif
        endif
       endif         	            
    enddo
  enddo 

  if (count1-1 /= totalelements) then
     print *, 'total element mismatch ', count1-1, totalelements
  endif

  do i=1, totalelements*8
     if (elementConn(i) > totalnodes) then
          print *, PetNo, 'node id out of bound', i/8, elementConn(i)
     endif
  enddo  
  print *, PetNo, "Before MeshAddElements"
  call ESMF_MeshAddElements(wamMesh, elementIds, elementTypes, elementConn,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

   print *, PetNo, "After MeshAddElements"

  deallocate(NumPerRow, ShuffleOrder, rowinds, petTable)
  deallocate(indList, baseind)
  deallocate(elementIds, elementTypes, elementConn)

  ! Info passed to the run routine
  is%wrap%wamdims = wamdims
  is%wrap%wamhgt => wamhgt
  is%wrap%startlevel = startlevel
  is%wrap%totallevels = totallevels
  is%wrap%wamtotalnodes = totalnodes
  is%wrap%localnodes = localnodes
  is%wrap%startnode = startnode
  is%wrap%wam2dMesh = wam2dMesh
  is%wrap%wamMesh = wamMesh
  is%wrap%PetNo = PetNo
  is%wrap%PetCnt = PetCnt

  return 
#else
    call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, & 
                 msg="- ESMF_NETCDF not defined when lib was compiled") 
    return
#endif

end subroutine InitGrids

! Log function which handles negative numbers
function log_hneg(in)
  real(ESMF_KIND_R8) :: in, log_hneg
  real(ESMF_KIND_R8), parameter :: log_min=1.0E-10

  if (in >= log_min) then
     log_hneg=log(in)
  else 
     log_hneg=log(log_min)
  endif
end function log_hneg


! Run Routine
subroutine RunRegrid(model, importState, exportState, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    integer, intent(out) :: rc
    
  type(ESMF_VM) :: vm
  type(InternalState)     :: is
  type(ESMF_RouteHandle)  :: routehandle
  type(ESMF_Field)        :: hgtfield, datafield, ipefield, wamfield
  type(ESMF_Field) :: ipe_north, ipe_east
  real(ESMF_KIND_R8), pointer :: hgtbuf(:,:), varbuf(:,:)
  real(ESMF_KIND_R8), allocatable :: tempbuf(:)
  real(ESMF_KIND_R8), pointer :: databuf(:), dstdata(:), wamhgt(:)
  real(ESMF_KIND_R8), pointer :: wamdata(:,:)
  integer :: nc1, ncid
  integer :: varid, data2id
  integer :: ndims, dimids(3), wamdims(3)
  integer :: PetNo, PetCnt
  integer :: totalnodes
  real(ESMF_KIND_R8), pointer :: fptr2d(:,:), fptr1d(:), fieldarray(:)
  type(ESMF_Arrayspec) :: arrayspec
  type(ESMF_Array) :: array
  type(ESMF_DistGrid) :: distgrid
  real(ESMF_KIND_R8) :: maxerror, minerror, totalerrors, deg2rad
  real(ESMF_KIND_R8) :: starttime, endtime, timesend(1), timereport(1)
  real(ESMF_KIND_R8) :: differr
  real(ESMF_KIND_R8), pointer :: varout(:), lonbuf(:)
  character(len=MAXNAMELEN) :: filename
  real, parameter :: PI=3.1415927
  integer :: i, ii, j, jj, j1, k, kk, l, count1
  integer :: localrc, status
  integer :: itemCount, localnodes, startlevel, totallevels, inlevels
  integer :: ubnd(2), lbnd(2)
  character(len=80), allocatable :: fieldNameList(:)
  integer :: numNodes, numElmts
  integer, save      :: slice=1
  integer :: count3(3), start3(3)
  ! constants requred for extrapolation of the density fields
  real(ESMF_KIND_R8) :: H, R, g0, re, mass, dist
  real(ESMF_KIND_R8) :: hgt_prev, H_prev, data_prev
  real(ESMF_KIND_R8) :: hgt_curr, H_curr, data_curr
  real(ESMF_KIND_R8) :: H_avg
  integer :: extrap_start_level

  ! Debug
  type(ESMF_Field) :: tmp_north, tmp_east


  rc = ESMF_SUCCESS
  R = 8.3141
  g0 = 9.80665
  re = 6.3712e03
  filename = 'wam2dmesh.nc'

  ! The WAM source level above which extrapolation is used to compute the values in the fixed height grid
  extrap_start_level=149  ! must be >=1 and <=inlevels

  !------------------------------------------------------------------------
  ! get global vm information
  !
  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  ! set up local pet info
  call ESMF_VMGet(vm, localPet=PetNo, petCount=PetCnt, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
  !
  ! query component for its internal state
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  wamdims=is%wrap%wamdims
  wamhgt => is%wrap%wamhgt
  PetCnt = is%wrap%PetCnt
  PetNo  = is%wrap%PetNo  
  startlevel = is%wrap%startlevel
  totallevels = is%wrap%totallevels

  ! Get the data from DATAWAM import fields
  ! Do 1D linear interpolation of the variables in the z direction to the fixed height grid first,

  call ESMF_StateGet(importstate, itemCount=itemCount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
  allocate(fieldNameList(itemCount))
  call ESMF_StateGet(importstate, itemNameList=fieldNameList, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

  call ESMF_StateGet(importState, itemName="height", &
       field=hgtfield, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldGet(hgtfield, farrayPtr=hgtbuf, computationalLbound=lbnd, &
       computationalUbound=ubnd, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
       return  ! bail out
  ! the heights are in meters, convert to kms
  hgtbuf = hgtbuf/1000.0
  inlevels = ubnd(2)-lbnd(2)+1

  call ESMF_StateGet(importState, itemName="temp_neutral", &
       field=datafield, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
       return  ! bail out
  call ESMF_FieldGet(datafield, farrayPtr=varbuf, computationalLbound=lbnd, &
     	 computationalUbound=ubnd, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
  localnodes=ubnd(1)-lbnd(1)+1
  totalnodes = localnodes * totallevels

  ! Error check extrap start level
  if ((extrap_start_level>inlevels) .or. &
      (extrap_start_level<1)) then
    call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
         msg="- extrap_start_level out of range", &
         line=__LINE__, &
         file=__FILE__,  &
         rcToReturn=rc)
    return
  endif			       

  ! save the 149th level of temp field
  allocate(tempbuf(localnodes))
  do i=1,localnodes
        tempbuf(i)=varbuf(i,extrap_start_level)
  enddo

#if 0
   if (slice==2) then
        ! write out text file for heights
        write(filename, "(A7,I2)") 'height.', PetNo+32
        open(100, file=filename)
        write(100, '(8(1X, F6.2))') wamhgt
        do i=lbnd(1), ubnd(1)
          write(100, '(8(1X, F6.2))') hgtbuf(i,:)
        enddo
        close(100)
   endif
#endif
   jj = 0
   do j=1, itemCount
     if (trim(fieldNameList(j)) .ne. "height") then
        jj = jj+1
        call ESMF_StateGet(importstate, itemname=fieldNameList(j), &
	   	   field=datafield, rc=rc)
  	if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
	call ESMF_FieldGet(datafield, farrayPtr=varbuf, computationalLbound=lbnd, &
       	   	computationalUbound=ubnd, rc=rc)
  	if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out

	! interpolate dataptr(localnodes, inlevels) to wamdata(localnodes, totallevels) 
        ! alone the second dimension
        ! using hgtbuf(localnodes,inlevels) as the source height and wamhgt(totallevels)
        ! as the destination heights 
        ! kk is the source index in 2nd dimension, k is the destination index   
        ! note hgtbuf has the heigth of the original WAM grid (150), the wamdim(3) is the
        ! fixed height WAM grid extended to 800KM, which is > 150

        if (j==1) allocate(wamdata(localnodes, totallevels))
        ! At the first time step, the values from the importstate are all invalid
        if (slice==1) then 
	  wamdata(:,:)=1.0
        else
	  if (PetNo == 0) then
             print *, 'fill field ',trim(fieldNameList(j))
          endif
          if (trim(fieldNameList(j)) == "O_Density" .or. &
              trim(fieldNameList(j)) == "O2_Density" .or. &
              trim(fieldNameList(j)) == "N2_Density") then
   
            if (trim(fieldNameList(j)) == "O_Density") then
	       mass = 16
            elseif (trim(fieldNameList(j)) == "O2_Density") then
               mass = 32
            else
               mass = 28
            endif
            ! do log interpolation 
  	    do i=1,localnodes
              kk = 1  ! source ind
              do k=startlevel, wamdims(3)
                 do while (kk<=inlevels .and. hgtbuf(i,kk)<wamhgt(k))
	            kk=kk+1
                 enddo
	         if (kk>extrap_start_level) then
		    hgt_prev=hgtbuf(i,extrap_start_level)
		    dist=re/(re+hgt_prev)
                    H_prev=R*tempbuf(i)/(mass*g0*dist*dist)
		    data_prev=varbuf(i,extrap_start_level)

                    do l=k,wamdims(3) 
                       ! use the following equation:
                       ! val(cur) = val(prev)*exp(h(prev)-h(cur))/H_avg)
                       !   Where:
                       !   H_avg=0.5*(H(prev)+H(cur))
                       !   H(x)=R*T(esl)/(mass*g(x))
                       !   R=8.3141,  
                       !   mass=16,32,28 for O, O2 and N2, 
                       !   g(x)=g0(re/(re+h(x)))^2
                       !   g0=9.80665
                       !   re=6.3712e+03
                       !   h(prev)= height in km at the previous level
                       !   val(prev)= value at the previous level
                       !   h(cur)= height in km at the current level
                       !   T(esl)= temperature at the extrap start level (temp values are extrapolated
                       !            upwards by just copying the top level)

                       ! Calculate info for this level
                       hgt_curr=wamhgt(l)
		       dist=re/(re+hgt_curr)	
                       H_curr=R*tempbuf(i)/(mass*g0*dist*dist)

		       ! Extrapolate data to this level
		       H_avg=0.5*(H_prev+H_curr)
                       data_curr=data_prev*exp((hgt_prev-hgt_curr)/H_avg)

		       ! Set extrapolated data in output array
 	               wamdata(i,l-startlevel+1)=data_curr

		       ! Set info for next pass
   		       hgt_prev=hgt_curr
                       H_prev=H_curr
		       data_prev=data_curr
                    enddo
	            exit
                 endif
                 if (kk>1) then
                    wamdata(i,k-startlevel+1)=exp((log_hneg(varbuf(i,kk))*(wamhgt(k)-hgtbuf(i,kk-1))+ &
	                 log_hneg(varbuf(i,kk-1))*(hgtbuf(i,kk)-wamhgt(k)))/ &
		         (hgtbuf(i,kk)-hgtbuf(i,kk-1)))
                 else 
	            wamdata(i,k-startlevel+1)=varbuf(i,kk)
                 endif
              enddo
            enddo
          else  
	    ! do linear interpolation
  	    do i=1,localnodes
              kk = 1  ! source ind
              do k=startlevel, wamdims(3)
                 do while (kk<=inlevels .and. hgtbuf(i,kk)<wamhgt(k))
	            kk=kk+1
                 enddo
	         if (kk>extrap_start_level) then
                    do l=k,wamdims(3) ! use the value as the highest level in the source grid
                                      ! to fill the remaining levels in the destination
 	               wamdata(i,l-startlevel+1)=varbuf(i,extrap_start_level)
                    enddo
	            exit
                 endif
                 if (kk>1) then
                    wamdata(i,k-startlevel+1)=(varbuf(i,kk)*(wamhgt(k)-hgtbuf(i,kk-1))+ &
	                 varbuf(i,kk-1)*(hgtbuf(i,kk)-wamhgt(k)))/ &
		         (hgtbuf(i,kk)-hgtbuf(i,kk-1))
                 else 
	            wamdata(i,k-startlevel+1)=varbuf(i,kk)
                 endif
              enddo
            enddo
          endif
        endif

        ! write out the variables from each processor
        if (slice == 2) then
        do i=0, PetCnt-1
         if (PetNo == i) then
            status = nf90_open(filename, NF90_WRITE, ncid)
            call CheckNCError(status, filename)
            count3(1)=localnodes
            count3(2)=totallevels
            count3(3)=1
            start3(1)=is%wrap%startnode
            start3(2)=1
            start3(3)=jj
            status = nf90_inq_varId(ncid, 'wamdata', varid)
            call CheckNCError(status, trim(fieldNameList(j)))
            status = nf90_put_var(ncid, varid, wamdata, & 
		    start3, count3)
            call CheckNCError(status, trim(fieldNameList(j)))
            status = nf90_close(ncid)
            call CheckNCError(status, filename)
          endif
          call ESMF_VMBarrier(vm)
        enddo
        endif

#define NEW_WAY_VECTOR
#ifdef NEW_WAY_VECTOR
        ! If these are winds, then just save results to be handled later 
        if (trim(fieldNameList(j)) .eq. "northward_wind_neutral") then

           call Copy_data_into_Field(is%wrap%wam_north, wamdata, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__)) &
                return  ! bail out
        
           ! Move on to next field
           cycle
        endif

        if (trim(fieldNameList(j)) .eq. "eastward_wind_neutral")  then

           call Copy_data_into_Field(is%wrap%wam_east, wamdata, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__)) &
                return  ! bail out

           ! Move on to next field
           cycle
        endif
#endif

        ! Create Field around wam data
        wamfield=ESMF_FieldCreate(is%wrap%wammesh, reshape(wamdata, (/localnodes*totallevels/)), &
		ESMF_INDEX_DELOCAL, datacopyflag=ESMF_DATACOPY_VALUE, meshloc=ESMF_MESHLOC_NODE, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
   
        ! Find the field of the same name in the export state -- that is built on IPE mesh
        call ESMF_StateGet(exportstate, itemname=fieldNameList(j), &
     	  field=ipefield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out

        ! Regrid!!
        call ESMF_FieldRegrid(wamfield, ipefield, is%wrap%routehandle, &
	       zeroregion=ESMF_REGION_SELECT, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out

        call ESMF_FieldDestroy(wamfield)
      endif
   enddo
   deallocate(wamdata)


   ! Handle winds using vector transformation
#ifdef NEW_WAY_VECTOR

#ifdef ANALYTIC_TEST
   ! Test by setting vectors to analytic field
   call Set_Tst_Cardinal_Vecs(is%wrap%wam_north, is%wrap%wam_east, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

   !! Get ipe northwind field
   call ESMF_StateGet(exportstate, itemname="northward_wind_neutral", &
        field=ipe_north, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

   ! Get ipe eastwind field
   call ESMF_StateGet(exportstate, itemname="eastward_wind_neutral", &
        field=ipe_east, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

#define VECTOR_XFORM
#ifdef VECTOR_XFORM
   !! Convert to Cartesian vector
   call Cardinal_to_Cart3D(is%wrap%wam_north, is%wrap%wam_east, &
                           is%wrap%wam_north_uvec, &
                           is%wrap%wam_east_uvec, &
                           is%wrap%wam_wind_3Dcart_vec, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

   !! Regrid
   call ESMF_FieldRegrid(is%wrap%wam_wind_3Dcart_vec, &
                         is%wrap%ipe_wind_3Dcart_vec, &
                         is%wrap%routehandle, &
                         zeroregion=ESMF_REGION_SELECT, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


   !! Convert back to cardinal vectors
   call Cart3D_to_Cardinal(is%wrap%ipe_wind_3Dcart_vec, &
                           is%wrap%ipe_north_uvec, &
                           is%wrap%ipe_east_uvec, &
                           ipe_north, ipe_east, rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

#else
   !! Regrid north
   call ESMF_FieldRegrid(is%wrap%wam_north, &
                         ipe_north, &
                         is%wrap%routehandle, &
                         zeroregion=ESMF_REGION_SELECT, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

   !! Regrid east
   call ESMF_FieldRegrid(is%wrap%wam_east, &
                         ipe_east, &
                         is%wrap%routehandle, &
                         zeroregion=ESMF_REGION_SELECT, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

#ifdef ANALYTIC_TEST
   ! Test by comparing to analytic field
   call Cmp_Tst_Cardinal_Vecs(ipe_north, ipe_east, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

#endif

    ! advance the time slice counter
    slice = slice + 1

  return

end subroutine RunRegrid

!Finalize Routine
subroutine Finalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

  type(ESMF_VM) :: vm
  type(InternalState)     :: is

  ! query component for its internal state
  nullify(is%wrap)
  call ESMF_GridCompGetInternalState(model, is, rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
  

  ! Destroy ESMF objects
  call ESMF_FieldDestroy(is%wrap%wam_north, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldDestroy(is%wrap%wam_east, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldDestroy(is%wrap%wam_north_uvec, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldDestroy(is%wrap%wam_east_uvec, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldDestroy(is%wrap%wam_wind_3Dcart_vec, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out


  call ESMF_FieldDestroy(is%wrap%ipe_north_uvec, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldDestroy(is%wrap%ipe_east_uvec, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldDestroy(is%wrap%ipe_wind_3Dcart_vec, rc=rc)  
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_MeshDestroy(is%wrap%wam2dmesh, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_MeshDestroy(is%wrap%wammesh, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_MeshDestroy(is%wrap%ipemesh, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldRegridRelease(is%wrap%routehandle)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! deallocate
  deallocate(is%wrap%wamhgt)
  deallocate(is%wrap)
  
  print *, 'Complete MEDIATOR'
end subroutine Finalize

subroutine ErrorMsgAndAbort(localPet)
    integer ::  localPet
  
    if (localPet >= 0) then
      write(*,*) "ERROR: Problem on processor ",localPet,". Please see the PET*.ESMF_LogFile files for a traceback."
    else
      write(*,*) "ERROR: Please see the PET*.LogFile files for a traceback."
    endif
  
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
end subroutine ErrorMsgAndAbort

!------------------------------------------------------------------------------
!
!  check CDF file error code
!
#undef  ESMF_METHOD
#define ESMF_METHOD "CheckNCError"
subroutine CheckNCError (ncStatus, errmsg)

    integer,          intent(in)  :: ncStatus
    character(len=*), intent(in)  :: errmsg

    integer, parameter :: nf90_noerror = 0

#ifdef ESMF_NETCDF
    if ( ncStatus .ne. nf90_noerror) then
        print '("NetCDF Error: ", A, " : ", A)', &
    		trim(errmsg),trim(nf90_strerror(ncStatus))
        call ErrorMsgAndAbort(-1)
    end if
#else
    call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, & 
                 msg="- ESMF_NETCDF not defined when lib was compiled") 
    return
#endif

end subroutine CheckNCError

!------------------------------------------------------------------------------
!
!  Convert 3D Spherical to 3D Cartisian if USE_CART3D_COORDSYS is set,
!    otherwise, just normalize the z field
!
#undef  ESMF_METHOD
#define ESMF_METHOD "convet2Cart"
#ifdef USE_CART3D_COORDSYS
subroutine convert2Cart (lon, lat, hgt, coords, rc)
   real(ESMF_KIND_R8):: lon, lat, hgt
   real(ESMF_KIND_R8):: coords(3)
   integer, optional :: rc

   real(ESMF_KIND_R8) :: earthradius, nhgt
   integer :: localrc

   if (present(rc)) rc=ESMF_FAILURE
   earthradius = 6371.0
   nhgt = 1+hgt/earthradius

   call c_esmc_sphdeg_to_cart(lon, lat, &
               coords(1), coords(2), coords(3), &
               localrc)
   if (localrc /= ESMF_SUCCESS) return 

   coords(1)=nhgt*coords(1)
   coords(2)=nhgt*coords(2)
   coords(3)=nhgt*coords(3)

   if (present(rc)) rc=ESMF_SUCCESS

end subroutine convert2Cart
#else
subroutine convert2Cart (lon, lat, hgt, coords, rc)
   real(ESMF_KIND_R8):: lon, lat, hgt
   real(ESMF_KIND_R8):: coords(3)
   integer, optional :: rc

   real(ESMF_KIND_R8) :: earthradius, nhgt
   integer :: localrc

   if (present(rc)) rc=ESMF_FAILURE
   earthradius = 6371.0
   nhgt = 1+hgt/earthradius

   coords(1)=lon
   coords(2)=lat
   coords(3)=nhgt

   if (present(rc)) rc=ESMF_SUCCESS

end subroutine convert2Cart
#endif
#undef  ESMF_METHOD
#define ESMF_METHOD "convet2Sphdeg"
subroutine convert2Sphdeg (coord1, coord2, coord3, lon, lat, hgt)
   real(ESMF_KIND_R8):: coord1, coord2, coord3
   real(ESMF_KIND_R8):: lon, lat, hgt

   real(ESMF_KIND_R8) :: earthradius, nhgt, rad2deg
   real, parameter :: PI=3.1415927
   integer :: localrc

   earthradius = 6371.0
   rad2deg = 180.0/PI
   nhgt = sqrt(coord1*coord1+coord2*coord2+coord3*coord3)
   hgt = (nhgt-1)*earthradius
   lon = atan(coord2/coord1)*rad2deg
   if (coord1 < 0) lon = lon + 180.0
   if (coord1 > 0 .and. coord2 < 0) lon = 360.0 + lon
   lat = 90-acos(coord3/nhgt)*rad2deg

end subroutine convert2Sphdeg


! Fill a Field with 3D Cartesian unit vectors corresponding to cardinal directions
! This assumes that the incoming fields are built on the same Mesh, and have an undistributed dim of 3
subroutine Set_Field_Cardinal_UVecs(north_field, east_field, rc)
  type(ESMF_Field) :: north_field, east_field
  type(ESMF_Mesh) :: mesh
  type(ESMF_VM) :: vm
  integer :: localPet
  integer :: numNodes
  real(ESMF_KIND_R8), allocatable :: nodeCoords(:)
  real(ESMF_KIND_R8), pointer :: north_field_ptr(:,:)
  real(ESMF_KIND_R8), pointer :: east_field_ptr(:,:)
  real(ESMF_KIND_R8) :: lat, lon
  integer :: localDECount, i
  integer :: rc
  
  ! Error checking
  real(ESMF_KIND_R8) :: max_coord(2), g_max_coord(2)
  real(ESMF_KIND_R8) :: min_coord(2), g_min_coord(2)
    

  ! debug
  real(ESMF_KIND_R8) :: len

  ! Get mesh
  call ESMF_FieldGet(north_field, mesh=mesh, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! If there are no DEs on this processor, then leave
  if (localDECount .eq. 0) then
     return
  endif

  ! If there is more than 1 DE then complain, because we aren't handling 
  ! that case right now
  if (localDECount .gt. 1) then
     return
  endif


  ! Get Coordinates
  call ESMF_MeshGet(mesh, numOwnedNodes=numNodes, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

   ! Allocate space for coordinates                                                                                                                          
   allocate(nodeCoords(3*numNodes))

   ! Set interpolated function                                                                                                                              
   call ESMF_MeshGet(mesh, ownedNodeCoords=nodeCoords, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Get pointer to north field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(north_field, 0, north_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Error checking of Field Bounds
  if ((lbound(north_field_ptr,1) .ne. 1) .or. &
       (ubound(north_field_ptr,1) .ne. 3) .or. &
       (lbound(north_field_ptr,2) .ne. 1) .or. &
       (ubound(north_field_ptr,2) .ne. numNodes)) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="north Field bounds wrong", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif

  ! Get pointer to east field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(east_field, 0, east_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Error checking of Field Bounds
  if ((lbound(east_field_ptr, 1) .ne. 1) .or. &
       (ubound(east_field_ptr, 1) .ne. 3) .or. &
       (lbound(east_field_ptr, 2) .ne. 1) .or. &
       (ubound(east_field_ptr, 2) .ne. numNodes)) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="east Field bounds wrong", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif

  ! For error checking
  min_coord= HUGE(min_coord)
  max_coord=-HUGE(max_coord)

  ! Loop setting unit vectors
  do i=1,numNodes

     ! Get position on sphere
     lon=nodeCoords(3*(i-1)+1)
     lat=nodeCoords(3*(i-1)+2)

     ! Get min/max of coord for error checking     
     if (lon < min_coord(1)) min_coord(1)=lon
     if (lon > max_coord(1)) max_coord(1)=lon
     if (lat < min_coord(2)) min_coord(2)=lat
     if (lat > max_coord(2)) max_coord(2)=lat

     ! Convert to radians
     lon=lon*ESMF_COORDSYS_DEG2RAD
     lat=lat*ESMF_COORDSYS_DEG2RAD
     
     ! Set east vector
     east_field_ptr(1,i)=cos(lon)
     east_field_ptr(2,i)=sin(lon)
     east_field_ptr(3,i)=0.0

     ! Set north vector
     north_field_ptr(1,i)=-sin(lat)*sin(lon)
     north_field_ptr(2,i)= sin(lat)*cos(lon)
     north_field_ptr(3,i)= cos(lat)

  enddo


  ! Get current vm
  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! set up local pet info
  call ESMF_VMGet(vm, localPet=localPet, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out


  ! Compute global max
  call ESMF_VMAllReduce(vm, max_coord, g_max_coord, 2, &
       ESMF_REDUCE_MAX, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_VMAllReduce(vm, min_coord, g_min_coord, 2, &
       ESMF_REDUCE_MIN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! Report error
  if ((g_max_coord(1) - g_min_coord(1)) < 20.0) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
 msg="Longitude range of grid unexpectedly small (< 20 deg) possibly using radians or 7.1.0 snapshot before 16", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif

  if ((g_max_coord(2) - g_min_coord(2)) < 20.0) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
msg="Latitude range of grid unexpectedly small (< 20 deg) possibly using radians or 7.1.0 snapshot before 16", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif

  ! Get rid of coordinates
  deallocate(nodeCoords)

  ! return success
  rc=ESMF_SUCCESS

end subroutine Set_Field_Cardinal_UVecs


! This function does the same thing as the above. However, it starts from Cart coordinates in the Mesh. 
! This is due to a bug in redisting the mesh. When we switch to 7.1.0 you can get rid of this subroutine and
! use the above for both meshes. 
subroutine Set_Field_Cardinal_UVecsCart(north_field, east_field, rc)
  type(ESMF_Field) :: north_field, east_field
  type(ESMF_Mesh) :: mesh
  integer :: numNodes
  real(ESMF_KIND_R8), allocatable :: nodeCoords(:)
  real(ESMF_KIND_R8), pointer :: north_field_ptr(:,:)
  real(ESMF_KIND_R8), pointer :: east_field_ptr(:,:)
  real(ESMF_KIND_R8) :: lat, lon, r
  integer :: localDECount, i
  integer :: rc
  real(ESMF_KIND_R8) :: x,y,z  
  real(ESMF_KIND_R8),parameter :: half_pi=1.5707963267949_ESMF_KIND_R8
  real(ESMF_KIND_R8),parameter :: two_pi=6.28318530717959_ESMF_KIND_R8

  ! debug
  real(ESMF_KIND_R8) :: len

  ! Get mesh
  call ESMF_FieldGet(north_field, mesh=mesh, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! If there are no DEs on this processor, then leave
  if (localDECount .eq. 0) then
     return
  endif

  ! If there is more than 1 DE then complain, because we aren't handling 
  ! that case right now
  if (localDECount .gt. 1) then
     return
  endif


  ! Get Coordinates
  call ESMF_MeshGet(mesh, numOwnedNodes=numNodes, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

   ! Allocate space for coordinates                                                                                                                          
   allocate(nodeCoords(3*numNodes))

   ! Set interpolated function                                                                                                                              
   call ESMF_MeshGet(mesh, ownedNodeCoords=nodeCoords, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Get pointer to north field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(north_field, 0, north_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Error checking of Field Bounds
  if ((lbound(north_field_ptr,1) .ne. 1) .or. &
       (ubound(north_field_ptr,1) .ne. 3) .or. &
       (lbound(north_field_ptr,2) .ne. 1) .or. &
       (ubound(north_field_ptr,2) .ne. numNodes)) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="north Field bounds wrong", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif

  ! Get pointer to east field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(east_field, 0, east_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Error checking of Field Bounds
  if ((lbound(east_field_ptr, 1) .ne. 1) .or. &
       (ubound(east_field_ptr, 1) .ne. 3) .or. &
       (lbound(east_field_ptr, 2) .ne. 1) .or. &
       (ubound(east_field_ptr, 2) .ne. numNodes)) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="east Field bounds wrong", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif

  ! Loop setting unit vectors
  do i=1,numNodes

     ! Get position on sphere
     x=nodeCoords(3*(i-1)+1)
     y=nodeCoords(3*(i-1)+2)
     z=nodeCoords(3*(i-1)+3)

     ! convert to lon/lat/r
     r=sqrt(x*x+y*y+z*z)

     lon=atan2(y,x)
     if (lon < 0.0) lon = lon + two_pi

     lat=half_pi-acos(z/r)

     ! Set east vector
     east_field_ptr(1,i)=cos(lon)
     east_field_ptr(2,i)=sin(lon)
     east_field_ptr(3,i)=0.0

     ! Set north vector
     north_field_ptr(1,i)=-sin(lat)*sin(lon)
     north_field_ptr(2,i)= sin(lat)*cos(lon)
     north_field_ptr(3,i)= cos(lat)

  enddo

  ! Get rid of coordinates
  deallocate(nodeCoords)

  ! return success
  rc=ESMF_SUCCESS

end subroutine Set_Field_Cardinal_UVecsCart


! Convert Cardinal vectors to one 3D cartesian vector
subroutine Cardinal_to_Cart3D(north_field, east_field, &
                              north_uvec, east_uvec, &
                              cart_vec, rc)
  type(ESMF_Field) :: north_field, east_field
  type(ESMF_Field) :: north_uvec, east_uvec
  type(ESMF_Field) :: cart_vec
  integer :: rc
  real(ESMF_KIND_R8), pointer :: north_field_ptr(:)
  real(ESMF_KIND_R8), pointer :: east_field_ptr(:)
  real(ESMF_KIND_R8), pointer :: north_uvec_ptr(:,:)
  real(ESMF_KIND_R8), pointer :: east_uvec_ptr(:,:)
  real(ESMF_KIND_R8), pointer :: cart_vec_ptr(:,:)
  integer :: localDECount, lDE
  integer :: clbnd(1), cubnd(1), i


  ! Get localDECount
  ! (Asssumes that all the incoming Fields have the same local DE Count)
  call ESMF_FieldGet(north_field, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out


  ! Loop over local DEs processing data
  do lDE=0,localDECount-1

     ! Get pointer to north field array
     call ESMF_FieldGet(north_field, lDE, north_field_ptr, &
          computationalLBound=clbnd, computationalUBound=cubnd, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to east field array
     call ESMF_FieldGet(east_field, lDE, east_field_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to north unit vector array
     call ESMF_FieldGet(north_uvec, lDE, north_uvec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to north unit vector array
     call ESMF_FieldGet(east_uvec, lDE, east_uvec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to east unit vector array
     call ESMF_FieldGet(east_uvec, lDE, east_uvec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     


     ! Get pointer to east unit vector array
     call ESMF_FieldGet(cart_vec, lDE, cart_vec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Loop over points processing
     do i=clbnd(1), cubnd(1)


#if 0
        cart_vec_ptr(1,i)=east_uvec_ptr(1,i)
        cart_vec_ptr(2,i)=east_uvec_ptr(2,i)
        cart_vec_ptr(3,i)=east_uvec_ptr(3,i)
#endif

        cart_vec_ptr(1,i)=east_field_ptr(i)*east_uvec_ptr(1,i)+ &
                          north_field_ptr(i)*north_uvec_ptr(1,i)

        cart_vec_ptr(2,i)=east_field_ptr(i)*east_uvec_ptr(2,i)+ &
                          north_field_ptr(i)*north_uvec_ptr(2,i)

        cart_vec_ptr(3,i)=east_field_ptr(i)*east_uvec_ptr(3,i)+ &
                          north_field_ptr(i)*north_uvec_ptr(3,i)
     enddo
  enddo

  ! return success
  rc=ESMF_SUCCESS

end subroutine Cardinal_to_Cart3D



! Convert 3D cartesian vector to Cardinal vectors
subroutine Cart3D_to_Cardinal(cart_vec, &
                              north_uvec, east_uvec, &
                              north_field, east_field, rc)
  type(ESMF_Field) :: cart_vec
  type(ESMF_Field) :: north_uvec, east_uvec
  type(ESMF_Field) :: north_field, east_field
  integer :: rc
  real(ESMF_KIND_R8), pointer :: north_field_ptr(:)
  real(ESMF_KIND_R8), pointer :: east_field_ptr(:)
  real(ESMF_KIND_R8), pointer :: north_uvec_ptr(:,:)
  real(ESMF_KIND_R8), pointer :: east_uvec_ptr(:,:)
  real(ESMF_KIND_R8), pointer :: cart_vec_ptr(:,:)
  integer :: localDECount, lDE
  integer :: clbnd(1), cubnd(1), i


  ! Get localDECount
  ! (Asssumes that all the incoming Fields have the same local DE Count)
  call ESMF_FieldGet(north_field, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out


  ! Loop over local DEs processing data
  do lDE=0,localDECount-1

     ! Get pointer to north field array
     call ESMF_FieldGet(north_field, lDE, north_field_ptr, &
          computationalLBound=clbnd, computationalUBound=cubnd, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to east field array
     call ESMF_FieldGet(east_field, lDE, east_field_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to north unit vector array
     call ESMF_FieldGet(north_uvec, lDE, north_uvec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to north unit vector array
     call ESMF_FieldGet(east_uvec, lDE, east_uvec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Get pointer to east unit vector array
     call ESMF_FieldGet(east_uvec, lDE, east_uvec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     


     ! Get pointer to east unit vector array
     call ESMF_FieldGet(cart_vec, lDE, cart_vec_ptr, &
          rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out     

     ! Loop over points processing
     do i=clbnd(1), cubnd(1)

        east_field_ptr(i)=cart_vec_ptr(1,i)*east_uvec_ptr(1,i)+ &
                          cart_vec_ptr(2,i)*east_uvec_ptr(2,i)+ &
                          cart_vec_ptr(3,i)*east_uvec_ptr(3,i)

        north_field_ptr(i)=cart_vec_ptr(1,i)*north_uvec_ptr(1,i)+ &
                           cart_vec_ptr(2,i)*north_uvec_ptr(2,i)+ &
                           cart_vec_ptr(3,i)*north_uvec_ptr(3,i)

!        if (i .lt. 1000) then
!           write(*,*) i," east=",east_field_ptr(i)," north=",north_field_ptr(i)
!        endif

     enddo
  enddo

  ! return success
  rc=ESMF_SUCCESS

end subroutine Cart3D_to_Cardinal



! Convert 3D cartesian vector to Cardinal vectors
subroutine Copy_data_into_Field(field, data, rc)
  type(ESMF_Field) :: field
  real(ESMF_KIND_R8), pointer :: data(:,:)
  integer :: rc
  real(ESMF_KIND_R8), pointer :: field_ptr(:)
  integer :: localDECount, lDE
  integer :: clbnd(1), cubnd(1), i


  ! Get localDECount
  ! (Asssumes that all the incoming Fields have the same local DE Count)
  call ESMF_FieldGet(field, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! If there's no local DEs then leave
  if (localDECount == 0) return

  ! Should only be 1 localDE
  if (localDECount > 1) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="Only Fields with one or less localDEs supported", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif


  ! Get Pointer
  call ESMF_FieldGet(field, 0, field_ptr, &
       computationalLBound=clbnd, computationalUBound=cubnd, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out     

  ! Error check bounds
  if ((cubnd(1)-clbnd(1)+1) .ne. &
       size(data)) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="Size of data array is different than size of Field", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif


  ! Copy data
  field_ptr=reshape(data, (/size(data,1)*size(data,2)/))

  ! return success
  rc=ESMF_SUCCESS

end subroutine Copy_data_into_Field

! Convert 3D cartesian vector to Cardinal vectors
subroutine Write_Field(field, rc)
  type(ESMF_Field) :: field
  integer :: rc
  real(ESMF_KIND_R8), pointer :: field_ptr(:)
  integer :: localDECount, lDE
  integer :: clbnd(1), cubnd(1), i


  ! Get localDECount
  ! (Asssumes that all the incoming Fields have the same local DE Count)
  call ESMF_FieldGet(field, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! If there's no local DEs then leave
  if (localDECount == 0) return

  ! Should only be 1 localDE
  if (localDECount > 1) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, & 
          msg="Only Fields with one or less localDEs supported", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
     return     
  endif


  ! Get Pointer
  call ESMF_FieldGet(field, 0, field_ptr, &
       computationalLBound=clbnd, computationalUBound=cubnd, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out     

   do i=clbnd(1), cubnd(1)
      write(*,*) i,"north=",field_ptr(i)
   enddo

  ! return success
  rc=ESMF_SUCCESS

end subroutine Write_Field

subroutine Set_Tst_Cardinal_Vecs(north_field, east_field, rc)
  type(ESMF_Field) :: north_field, east_field
  type(ESMF_Mesh) :: mesh
  integer :: numNodes
  real(ESMF_KIND_R8), allocatable :: nodeCoords(:)
  real(ESMF_KIND_R8), pointer :: north_field_ptr(:)
  real(ESMF_KIND_R8), pointer :: east_field_ptr(:)
  real(ESMF_KIND_R8) :: lat, lon, x, y, z
  integer :: localDECount, i
  integer :: rc
  
  ! debug
  real(ESMF_KIND_R8) :: len

  ! Get mesh
  call ESMF_FieldGet(north_field, mesh=mesh, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! If there are no DEs on this processor, then leave
  if (localDECount .eq. 0) then
     return
  endif

  ! If there is more than 1 DE then complain, because we aren't handling 
  ! that case right now
  if (localDECount .gt. 1) then
     return
  endif


  ! Get Coordinates
  call ESMF_MeshGet(mesh, numOwnedNodes=numNodes, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

   ! Allocate space for coordinates                                                                                                                          
   allocate(nodeCoords(3*numNodes))

   ! Set interpolated function                                                                                                                              
   call ESMF_MeshGet(mesh, ownedNodeCoords=nodeCoords, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Get pointer to north field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(north_field, 0, north_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Get pointer to east field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(east_field, 0, east_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

  ! Loop setting unit vectors
  do i=1,numNodes

     ! Get position on sphere
     lon=nodeCoords(3*(i-1)+1)
     lat=nodeCoords(3*(i-1)+2)
     

     ! Get x,y,z coordinates
     call c_esmc_sphdeg_to_cart(lon, lat, &
               x, y, z, &
               rc)

     ! Set east vector
     east_field_ptr(i)=cos(lon*ESMF_COORDSYS_DEG2RAD)

     ! Set north vector
     north_field_ptr(i)=(-sin(lat*ESMF_COORDSYS_DEG2RAD)* &
                         sin(lon*ESMF_COORDSYS_DEG2RAD))

#if 0
     ! Set east vector
     east_field_ptr(i)=100.0

     ! Set north vector
     north_field_ptr(i)=200.0
#endif

  enddo

  ! Get rid of coordinates
  deallocate(nodeCoords)

  ! return success
  rc=ESMF_SUCCESS

end subroutine Set_Tst_Cardinal_Vecs

subroutine Cmp_Tst_Cardinal_Vecs(north_field, east_field, rc)
  type(ESMF_Field) :: north_field, east_field
  type(ESMF_Mesh) :: mesh
  type(ESMF_VM) :: vm
  integer :: numNodes, localPet
  real(ESMF_KIND_R8), allocatable :: nodeCoords(:)
  real(ESMF_KIND_R8), pointer :: north_field_ptr(:)
  real(ESMF_KIND_R8), pointer :: east_field_ptr(:)
  real(ESMF_KIND_R8) :: lat, lon, r
  integer :: localDECount, i
  integer :: rc
  real(ESMF_KIND_R8) :: x,y,z  

  real(ESMF_KIND_R8) :: relerr, err
  real(ESMF_KIND_R8) :: max_relerr(2),g_max_relerr(2)
  real(ESMF_KIND_R8) :: max_err(2), g_max_err(2)
  real(ESMF_KIND_R8) :: north_exact, east_exact
  real(ESMF_KIND_R8),parameter :: half_pi=1.5707963267949_ESMF_KIND_R8
  real(ESMF_KIND_R8),parameter :: two_pi=6.28318530717959_ESMF_KIND_R8
  

  ! debug
  real(ESMF_KIND_R8) :: max_lat, max_lon
  real(ESMF_KIND_R8) :: min_lat, min_lon

  ! Get mesh
  call ESMF_FieldGet(north_field, mesh=mesh, &
                     localDECount=localDECount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! If there are no DEs on this processor, then leave
  if (localDECount .eq. 0) then
     return
  endif


  ! If there is more than 1 DE then complain, because we aren't handling 
  ! that case right now
  if (localDECount .gt. 1) then
     return
  endif


  ! Get Coordinates
  call ESMF_MeshGet(mesh, numOwnedNodes=numNodes, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

   ! Allocate space for coordinates                                                                                                                          
   allocate(nodeCoords(3*numNodes))

   ! Set interpolated function                                                                                                                              
   call ESMF_MeshGet(mesh, ownedNodeCoords=nodeCoords, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Get pointer to north field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(north_field, 0, north_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Get pointer to east field array
  ! (Should only be 1 localDE)                                                                                                                                 
  call ESMF_FieldGet(east_field, 0, east_field_ptr, &
       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out


  ! Debug
  min_lon= 1000.00
  max_lon=-1000.00

  min_lat= 1000.00
  max_lat=-1000.00


  ! Init max
  max_relerr= -HUGE(max_relerr)
  max_err= -HUGE(max_err)


  ! Loop setting unit vectors
  do i=1,numNodes

     ! Get position on sphere
     if (ESMF_VERSION_MAJOR > 7) then
        lon=nodeCoords(3*(i-1)+1)*ESMF_COORDSYS_DEG2RAD
        lat=nodeCoords(3*(i-1)+2)*ESMF_COORDSYS_DEG2RAD
        r=nodeCoords(3*(i-1)+3)
     else if ((ESMF_VERSION_MAJOR == 7) .and. (ESMF_VERSION_MINOR > 0)) then
        lon=nodeCoords(3*(i-1)+1)*ESMF_COORDSYS_DEG2RAD
        lat=nodeCoords(3*(i-1)+2)*ESMF_COORDSYS_DEG2RAD
        r=nodeCoords(3*(i-1)+3)
     else 
        ! Get position on sphere
        x=nodeCoords(3*(i-1)+1)
        y=nodeCoords(3*(i-1)+2)
        z=nodeCoords(3*(i-1)+3)
        
        ! convert to lon/lat/r
        r=sqrt(x*x+y*y+z*z)
        
        lon=atan2(y,x)
        if (lon < 0.0) lon = lon + two_pi
        lat=half_pi-acos(z/r)
     endif


     ! Debug     
     if (lon < min_lon) min_lon=lon
     if (lon > max_lon) max_lon=lon
     if (lat < min_lat) min_lat=lat
     if (lat > max_lat) max_lat=lat


     ! calc exact values
     east_exact=cos(lon)
     north_exact=(-sin(lat)*sin(lon))

#if 0
     east_exact=100.0
     north_exact=200.0

     if (lat*ESMF_COORDSYS_RAD2DEG < -85.0) cycle
     if (lat*ESMF_COORDSYS_RAD2DEG >  85.0) cycle
#endif


#if 0
     if (relerr > 0.3) then
        write(*,*) i,"east rel error=",east_field_ptr(i),east_exact,relerr,&
             lon*ESMF_COORDSYS_RAD2DEG,lat*ESMF_COORDSYS_RAD2DEG
     endif

     if (err > 0.3) then
        write(*,*) i,"east error=",east_field_ptr(i),east_exact,err,&
             lon*ESMF_COORDSYS_RAD2DEG,lat*ESMF_COORDSYS_RAD2DEG
     endif
#endif

     ! Calc north error
     err=abs(north_field_ptr(i)-north_exact)
     if (north_exact .ne. 0.0) then
        relerr=abs(err/north_exact)
     else
        relerr=err
     endif

    if (err > max_err(1)) then
        max_err(1)=err
     endif

     if (relerr > max_relerr(1)) then
        max_relerr(1)=relerr
     endif

#if 0
     if (relerr > 0.3) then
        write(*,*) i,"north rel error=",north_field_ptr(i),north_exact,relerr, &
             lon*ESMF_COORDSYS_RAD2DEG,lat*ESMF_COORDSYS_RAD2DEG
     endif

     if (err > 0.3) then
        write(*,*) i,"north error=",north_field_ptr(i),north_exact,err, &
             lon*ESMF_COORDSYS_RAD2DEG,lat*ESMF_COORDSYS_RAD2DEG
     endif
#endif

     ! Calc east error
     err=abs(east_field_ptr(i)-east_exact)
     if (east_exact .ne. 0.0) then
        relerr=abs(err/east_exact)
     else
        relerr=err
     endif

     if (err > max_err(2)) then
        max_err(2)=err
     endif

     if (relerr > max_relerr(2)) then
        max_relerr(2)=relerr
     endif
  enddo

  ! Get current vm
  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! set up local pet info
  call ESMF_VMGet(vm, localPet=localPet, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! Compute global max
  call ESMF_VMReduce(vm, max_err, g_max_err, 2, &
       ESMF_REDUCE_MAX, 0, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_VMReduce(vm, max_relerr, g_max_relerr, 2, &
       ESMF_REDUCE_MAX, 0, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
  
  
  ! Only write out if Pet 0
  if (localPet==0) then
     write(*,*) "max relerr (north/east)=",g_max_relerr(1),g_max_relerr(2)
     write(*,*) "max err    (north/east)=",g_max_err(1),g_max_err(2)
  endif

  !write(*,*) "c min/max lon=",min_lon,max_lon
  !write(*,*) "c min/max lat=",min_lat,max_lat


  ! Get rid of coordinates
  deallocate(nodeCoords)

  ! return success
  rc=ESMF_SUCCESS

end subroutine Cmp_Tst_Cardinal_Vecs


end module
