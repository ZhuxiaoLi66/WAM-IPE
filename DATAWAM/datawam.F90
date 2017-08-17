module DATAWAM
 
  !-----------------------------------------------------------------------------
  ! DATAWAM: DATA WAM Component
  !
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

#ifdef ESMF_NETCDF
  use netcdf
#endif
  
  implicit none
  
  private

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
  
  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! Derive from the generic NUOPC model component
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! Provide InitializeP0 to switch from default IPDv00 to IPDv01
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p2"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS

    ! exportable fields: typical ATM export fields
    call NUOPC_Advertise(exportState, StandardNames=(/ &
      "northward_wind_neutral         ", &
      "eastward_wind_neutral          ", &
    !  "upward_wind_neutral            ", &
      "temp_neutral                   ", &
    !  "O_Density                      ", &
    !  "O2_Density                     ", &
    !  "N2_Density                     ", &
    !  "NO_Density                     ", &
    !  "N4S_Density                    ", &
    !  "N2D_Density                    ", &
    !  "H_Density                      ", &
    !  "He_Density                     ", &
      "height                         "  &
      /), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set Component name so it becomes identifiable in the output
    call ESMF_GridCompSet(gcomp, name="DATAWAM", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Mesh)   :: mesh
    integer           :: levels

    rc = ESMF_SUCCESS
    
    ! Hardcode the original WAM grid vertical levels for now.
    levels = 150
    ! Create a LocStream to represent the reduced Gassian grid from the original WAM model 
    ! with the same distribution as used in the WAM.  Only create the 2D layer and put
    ! the vertical layer as an undistributed dimension.  
    ! For the timebeing, read the grid information from a grid file
    ! call createWAMGrid(locStream, rc)
    call createWAMGrid(gcomp, mesh, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! realize connected Fields in the exportState and data initialize
    call realizeConnectedFields(exportState, mesh, &
      levels, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if 0
  ! Create 2D WAM as a location stream
  subroutine createWAMGrid(locStream, rc)

     type(ESMF_LocStream), intent(out) :: locStream
     integer, intent(out)              :: rc

     integer             :: PetNo, PetCnt
     type(ESMF_VM)       :: vm
     integer             :: nc1, varid
     integer             :: ndims, dimids(2), dims(2)
     real(ESMF_KIND_R8), pointer    :: dstlon(:,:), dstlat(:,:)
     real(ESMF_KIND_R8), pointer    :: lons(:), lats(:)
     integer(ESMF_KIND_I4), pointer :: NumPerRow(:), ShuffleOrder(:)
     integer             :: myrows
     integer(ESMF_KIND_I4), pointer :: rowinds(:), startind(:)
     integer             :: i, j, k, next, ind1
     integer(ESMF_KIND_I4), pointer :: indList(:)
     type(ESMF_DistGrid) :: distgrid
     character(len=ESMF_MAX_STRLEN) :: gridfilename
     integer             :: status
     type(InternalState)  :: is

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

#if 0
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
#endif

      gridfilename = 'data/wam3dgridnew.nc'

      ! Read in the WAM grid 
      ! This part of the code will be replaced by getting the information from the WAM 
      ! model directly

#ifdef ESMF_NETCDF

      status = nf90_open(path=gridfilename, mode=nf90_nowrite, ncid=nc1)
      call CheckNCError(status, gridfilename)
      status = nf90_inq_varid(nc1,'lons', varid)
      call CheckNCError(status, 'lons')
      status = nf90_inquire_variable(nc1, varid, ndims=ndims, dimids = dimids)
      call CheckNCError(status, 'lons')
      status = nf90_inquire_dimension(nc1,dimids(1), len=dims(1))
      call CheckNCError(status, 'lons 1st dimension')
      status = nf90_inquire_dimension(nc1,dimids(2), len=dims(2))
      call CheckNCError(status, 'lons 2nd dimension')

      ! WAM dimension order:  lons, lats (192, 94)
      allocate(dstlon(dims(1), dims(2)), &
  	   dstlat(dims(1), dims(2)))
      status = nf90_get_var(nc1, varid, dstlon)
      call CheckNCError(status, 'lons')
      status = nf90_inq_varid(nc1,'lats', varid)
      call CheckNCError(status, 'lats')
      status = nf90_get_var(nc1, varid, dstlat)
      call CheckNCError(status, 'lats')

      allocate(NumPerRow(dims(2)), ShuffleOrder(dims(2)))
      status = nf90_inq_varid(nc1,'NumPerRow', varid)
      call CheckNCError(status, 'NumPerRow')
      status = nf90_get_var(nc1, varid, NumPerRow)
      call CheckNCError(status, 'NumPerRow')
      status = nf90_inq_varid(nc1,'ShuffleOrder', varid)
      call CheckNCError(status, 'ShuffleOrder')
      status = nf90_get_var(nc1, varid, ShuffleOrder)
      call CheckNCError(status, 'ShuffleOrder')
      status= nf90_close(nc1)
      call CheckNCError(status, is%wrap%dstfilename)

      ! find the total number of nodes in each processor and create local index table
      localnodes=0
      myrows = dims(2)/PetCnt
      if ((dims(2)-myrows*PetCnt) > PetNo) myrows = myrows+1
      allocate(rowinds(myrows))
      allocate(startind(dims(2)))
      next = 0
      ind1 = 0
      do i=1,dims(2)
         ind=ShuffleOrder(i)
         if (next == PetNo) then
           ind1=ind1+1
           rowinds(ind1)=ind
           localnodes=localnodes + numPerRow(ind)
         endif
         next=next+1
         if (next == PetCnt) next=0
         if (i>1) then
            startind(i)=startind(i-1)+numPerRow(i-1)
         else
            startind(i)=1
         endif
      enddo

      ! Create a distgrid using a collapsed 1D index array based on the local row index
      allocate(indList(localnodes))
      k=1
      do i=1,myrows
        ind=rowinds(i)
        do j=1,numPerRow(ind)
           indList(k)=startInd(ind)+j
           k=k+1
        enddo
      enddo
      distgrid = ESMF_DistGridCreate(indList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      ! Create LocStream
      locSream = ESMF_LocStreamCreate(distgrid, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)  
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
 
      ! Add the Key::Lons and Key:Lats
      allocate(lons(localnodes), lats(localnodes))
      k=1
      do i=1,myrows
         ind=rowinds(i)
         do j=1,numPerRow(ind)
	    lons(k) = dstlon(j,ind)
            lats(k)  = dstlat(j,ind)
	    k=k+1
         enddo
      enddo

      call ESMF_LocStreamAddKey(locStream, "ESMF::Lons",lons, keyUnits="degree_east", &
         keyLongName="Longitude", rc=localrc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      call ESMF_LocStreamAddKey(locStream, "ESMF::Lats",lats, keyUnits="degree_north", &
       keyLongName="Latitude", rc=localrc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

      deallocate(dstlon, dstlat)
      deallocate(NumPerRow, ShuffleOrder)
      deallocate(rowinds, startind, indList)

#else
    call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, & 
                 msg="- ESMF_NETCDF not defined when lib was compiled") 
    return
#endif

    end subroutine createWAMgrid
#else
  ! Create 2D WAM as a ESMF_Mesh with only distgrid (no coordinates)
  subroutine createWAMGrid(gcomp, mesh, rc)

     type(ESMF_GridComp)  :: gcomp
     type(ESMF_Mesh), intent(out) :: mesh
     integer, intent(out)              :: rc

     integer             :: PetNo, PetCnt
     type(ESMF_VM)       :: vm
     integer             :: nc1, varid
     integer             :: ndims, dimids(2), dims(2)
     integer(ESMF_KIND_I4), pointer :: NumPerRow(:), ShuffleOrder(:)
     integer             :: myrows
     real(ESMF_KIND_R8), pointer    :: dstlon(:,:), dstlat(:,:)
     integer(ESMF_KIND_I4), pointer :: rowinds(:), startind(:)
     integer             :: i, j, k, next, ind1
     integer(ESMF_KIND_I4), pointer :: indList(:)
     type(ESMF_DistGrid) :: distgrid
     character(len=MAXNAMELEN) :: gridfilename
     integer             :: status
     integer             :: localnodes, ind
     type(InternalState)  :: is

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

     ! Allocate memory for the internal state and set it in the Component.
       allocate(is%wrap, stat=rc)
       if (ESMF_LogFoundAllocError(statusToCheck=rc, &
         msg="Allocation of the internal state memory failed.", &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
       call ESMF_GridCompSetInternalState(gcomp, is, rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

      gridfilename = 'data/wam3dgridnew.nc'

      ! Read in the WAM grid 
      ! This part of the code will be replaced by getting the information from the WAM 
      ! model directly

#ifdef ESMF_NETCDF

      status = nf90_open(path=gridfilename, mode=nf90_nowrite, ncid=nc1)
      call CheckNCError(status, gridfilename)
      status = nf90_inq_varid(nc1,'lons', varid)
      call CheckNCError(status, 'lons')
      status = nf90_inquire_variable(nc1, varid, ndims=ndims, dimids = dimids)
      call CheckNCError(status, 'lons')
      status = nf90_inquire_dimension(nc1,dimids(1), len=dims(1))
      call CheckNCError(status, 'lons 1st dimension')
      status = nf90_inquire_dimension(nc1,dimids(2), len=dims(2))
      call CheckNCError(status, 'lons 2nd dimension')

      ! WAM dimension order:  lons, lats (192, 94)
      allocate(dstlon(dims(1), dims(2)), &
  	   dstlat(dims(1), dims(2)))
      status = nf90_get_var(nc1, varid, dstlon)
      call CheckNCError(status, 'lons')
      status = nf90_inq_varid(nc1,'lats', varid)
      call CheckNCError(status, 'lats')
      status = nf90_get_var(nc1, varid, dstlat)
      call CheckNCError(status, 'lats')

      allocate(NumPerRow(dims(2)), ShuffleOrder(dims(2)))
      status = nf90_inq_varid(nc1,'NumPerRow', varid)
      call CheckNCError(status, 'NumPerRow')
      status = nf90_get_var(nc1, varid, NumPerRow)
      call CheckNCError(status, 'NumPerRow')
      status = nf90_inq_varid(nc1,'ShuffleOrder', varid)
      call CheckNCError(status, 'ShuffleOrder')
      status = nf90_get_var(nc1, varid, ShuffleOrder)
      call CheckNCError(status, 'ShuffleOrder')
      status= nf90_close(nc1)
      call CheckNCError(status, gridfilename)

      ! find the total number of nodes in each processor and create local index table
      localnodes=0
      myrows = dims(2)/PetCnt
      if ((dims(2)-myrows*PetCnt) > PetNo) myrows = myrows+1
      allocate(rowinds(myrows))
      allocate(startind(dims(2)))
      next = 0
      ind1 = 0
      do i=1,dims(2)
         ind=ShuffleOrder(i)
         if (next == PetNo) then
           ind1=ind1+1
           rowinds(ind1)=ind
           localnodes=localnodes + numPerRow(ind)
         endif
         next=next+1
         if (next == PetCnt) next=0
         !if (i>1) then
         !   startind(i)=startind(i-1)+numPerRow(i-1)
         !else
         !   startind(i)=1
         !endif
      enddo

      ! Create a distgrid using a collapsed 1D index array based on the local row index
      allocate(indList(localnodes))
      k=1
      do i=1,myrows
        ind=rowinds(i)
        do j=1,numPerRow(ind)
           indList(k)=dims(1)*(ind-1)+j
           k=k+1
        enddo
      enddo
      distgrid = ESMF_DistGridCreate(indList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      ! Create mesh using the distgrid as the nodaldistgrid,  no elemdistgrid available
      ! just use nodeldistgrid for both
      mesh = ESMF_MeshCreate(distgrid,distgrid,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      !deallocate(dstlon, dstlat)
      !deallocate(NumPerRow, ShuffleOrder)
      !deallocate(rowinds, startind, indList)

      is%wrap%dims = dims
      is%wrap%lons => dstlon
      is%wrap%lats => dstlat
      is%wrap%numPerRow => numPerRow
      is%wrap%rowinds => rowinds
      is%wrap%myrows = myrows
      is%wrap%indlist => indList
      is%wrap%PetNo = PetNo
      is%wrap%PetCnt = PetCnt

#else
    call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, & 
                 msg="- ESMF_NETCDF not defined when lib was compiled") 
    return
#endif
    end subroutine createWAMgrid
#endif

    subroutine realizeConnectedFields(state, mesh, levels, rc)
      ! TODO: this method may move into the NUOPC_ utility layer
      type(ESMF_State)                :: state
      type(ESMF_Mesh)                 :: mesh
      integer, intent(in)             :: levels
      integer, intent(out), optional  :: rc
      ! local variables
      character(len=80), allocatable  :: fieldNameList(:)
      integer                         :: i, itemCount, k
      type(ESMF_Field)                :: field
      type(ESMF_ArraySpec)            :: arrayspec

      if (present(rc)) rc = ESMF_SUCCESS
      
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
	   ungriddedLBound=(/1/), ungriddedUBound=(/levels/),&
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

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_State)                :: importState, exportState
    type(ESMF_Time)                 :: currTime, stopTime
    type(ESMF_TimeInterval)         :: timeStep
    character(len=80), allocatable  :: fieldNameList(:)
    integer                         :: i, itemCount
    real(ESMF_KIND_R8), pointer     :: farrayptr(:,:)
    type(ESMF_Field)                :: field

    integer, save                   :: slice=1
    
    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing DATAWAM from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! update data in Fields in the exportState
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
      
    ! Fill the height field first
    call ESMF_StateGet(exportState, itemName="height", field=field, &
        rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    call fillWAMHeight(gcomp, field, farrayptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    do i=1, itemCount
      call ESMF_StateGet(exportState, itemName=fieldNameList(i), field=field, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! Fill the data using analytical data
      ! will get the data from the WAM model once hook it up with the real WAM model
      if (fieldNameList(i) /= "height") then
        call fillWAMField(gcomp, field, farrayptr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    enddo

    call ESMF_ClockGet(clock, stopTime=stopTime, currTime=currTime, &
      timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! advance the time slice counter
    slice = slice + 1

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if 0 
  ! Code to use location stream and construct synthetic fields using the coordinates
  subroutine fillWAMField(field, fieldname, rc)
 
    type(ESMF_Field)     :: field
    character(len=*)     :: fieldname
    integer, intent(out) :: rc

    type(ESMF_LocStream) :: locstream
    real(ESMF_KIND_R8)   :: lons(:), lats(:)
    real(ESMF_KIND_R8)   :: fptr(:,:)
    interger             :: lbnd(2), ubnd(2)
    real(ESMF_KIND_R8)   :: x, highest, lowest, interval

    ! create analytical fields and will get the data from the WAM model eventually
    call ESMF_FieldGet(field, locstream=locstream, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
 
    call ESMF_FileGet(field, localDe=0, exclusiveLBound=lbnd, exclusiveUBound=ubnd,
    	 farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    call ESMF_LocStreamGetKey(locstream, localDE=0, keyName="ESMF:Lons", &
         farray=lons, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    call ESMF_LocStreamGetKey(locstream, localDE=0, keyName="ESMF:Lats", &
         farray=lats, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    if (fieldname /= "height") then
      do i=lbnd(1),ubnd(1)
         do j=lbnd(2), ubnd(2)
            fptr(i,j)=(151-j)*0.35*(1.0+cos(lats(i)*deg2rad)* &
		    cos(lats(i)*deg2rad)*cos(2*lons(i)*deg2rad))		
         enddo
      enddo
    else ! construct height fields, must be monotonically increase from 0 to 800
      call RANDOM_SEED(is%wrap%PetNo)
      levels = ubnds(2)-lbnds(1)+1
      do i=lbnd(1),ubnd(1)
         ! use random number generator to generate lowest and highest heights
         ! The lowest height will be between 10 and 20 KM
         call RANDOM_NUMBER(x)
         lowest = x*10+10;
         ! The highest height with be between 400 and 800 KM
         highest = x*400+400;
         intreval = (hightest-lowest)/(levels-1);
         do j=lbnd(2), ubnd(2)
           fptr(i,j)=lowest+interval*(j-lbnd(2))
         enddo
      enddo
    endif
  end subroutine fillWAMField

#else

  ! Create analytical fields for the 2D WAM built on a ESMF_Mesh
  subroutine fillWAMField(gcomp, field, hgtptr, rc)
    
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_Field)     :: field
    real(ESMF_KIND_R8), pointer :: hgtptr(:,:)
    integer, intent(out) :: rc

    type(ESMF_Array) :: array
    real(ESMF_KIND_R8), pointer   :: lons(:,:), lats(:,:)
    real(ESMF_KIND_R8), pointer   :: fptr(:,:)
    integer             :: lbnd(2), ubnd(2)
    real(ESMF_KIND_R8) :: lon, lat, deg2rad
    integer            :: xind, yind, i, j, ind, levels
    real, parameter :: PI=3.14159265
    type(InternalState)  :: is
    real(ESMF_KIND_R8)   :: x, highest, lowest, interval

    deg2rad = PI/180.0

    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lons => is%wrap%lons
    lats => is%wrap%lats
        
    call ESMF_FieldGet(field, localDe=0, exclusiveLBound=lbnd, &
    	 exclusiveUBound=ubnd, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! use is%wrap%indList to find the index in the original lat/lon array
    do i=lbnd(1),ubnd(1)
         ind=is%wrap%indList(i-lbnd(1)+1)
         yind=(ind/is%wrap%dims(1))+1
         xind = ind-(yind-1)*is%wrap%dims(1)
         lon=lons(xind,yind)
         lat=lats(xind,yind)
         do j=lbnd(2), ubnd(2)
            fptr(i,j)=1.0+hgtptr(i,j)*0.004+cos(lat*deg2rad)* &
		    cos(lat*deg2rad)*cos(2*lon*deg2rad)
	 enddo
    enddo
    rc = ESMF_SUCCESS
  end subroutine fillWAMField

  ! Create a random fields for the 2D WAM built on a ESMF_Mesh
  subroutine fillWAMHeight(gcomp, field, fptr, rc)

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_Field)     :: field
    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    integer, intent(out) :: rc

    type(InternalState)  :: is
    real(ESMF_KIND_R8)   :: x, highest, lowest, interval
    integer             :: lbnd(2), ubnd(2)
    integer             :: levels, i, j
    
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldGet(field, localDe=0, exclusiveLBound=lbnd, &
    	 exclusiveUBound=ubnd, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! construct height fields, must be monotonically increase from 0 to 800
    call RANDOM_SEED(is%wrap%PetNo)
    levels = ubnd(2)-lbnd(2)+1
    do i=lbnd(1),ubnd(1)
         ! use random number generator to generate lowest and highest heights
         ! The lowest height will be between 10 and 20 KM
         call RANDOM_NUMBER(x)
         lowest = x*10+10;
         ! The highest height with be between 500 and 800 KM
         highest = x*300+500;
         interval = (highest-lowest)/(levels-1);
         do j=lbnd(2), ubnd(2)
           fptr(i,j)=lowest+interval*(j-lbnd(2))
         enddo
     enddo
     rc = ESMF_SUCCESS

  end subroutine fillWAMHeight
#endif

  end subroutine ModelAdvance
  !-----------------------------------------------------------------------------

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
  
  deallocate(is%wrap%numPerRow, is%wrap%rowinds)
  deallocate(is%wrap%lons, is%wrap%lats)
  deallocate(is%wrap%indList)  
  deallocate(is%wrap)
  
  print *, 'Complete DATAWAM'
end subroutine Finalize

!------------------------------------------------------------------------------
!
!  check CDF file error code
!
#undef  ESMF_METHOD
#define ESMF_METHOD "CheckNCError"
subroutine CheckNCError (ncStatus, errmsg)

    integer,          intent(in)  :: ncStatus
    character(len=*), intent(in)  :: errmsg

    character(len=256) :: msg
    integer, parameter :: nf90_noerror = 0

#ifdef ESMF_NETCDF
    if ( ncStatus .ne. nf90_noerror) then
        write(msg, '("NetCDF Error: ", A, " : ", A)') &
    		trim(errmsg),trim(nf90_strerror(ncStatus))
        call ESMF_LogSetError(ESMF_FAILURE, &
	      msg=msg, &
	      line=__LINE__, &
              file=__FILE__)

        !bail out 
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
#else
    call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, & 
                 msg="- ESMF_NETCDF not defined when lib was compiled") 
    return
#endif

end subroutine CheckNCError

end module
