!nm20140919: originally copied from MODEL.F90 but modified for IPE
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
#define ESMF_CONTEXT  line=__LINE__,file=__FILE__,method=ESMF_METHOD

module ipeCap

  !-----------------------------------------------------------------------------
  ! IPE Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS           => SetServices,          &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance,        &
    model_label_Finalize       => label_Finalize
  use module_finalize_IPE,        only: finalize_IPE
  use module_sub_myIPE_Init,      only: myIPE_Init
  use module_update_IPE,          only: update_IPE
  use module_IPE_dimension,       only: NMP, NLP, NPTS2D
  use module_FIELD_LINE_GRID_MKS, only: JMIN_IN, JMAX_IS, &
    plasma_grid_3d, plasma_grid_Z, IGCOLAT, IGLON, &
    wamField
  use module_input_parameters,    only: lps, lpe, mps, mpe, &
    mesh_height_min, mesh_height_max, mesh_write, mesh_write_file, &
    mype, swESMFTime


  implicit none

  integer, parameter :: importFieldCount = 7
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "temp_neutral          ", &
      "eastward_wind_neutral ", &
      "northward_wind_neutral", &
      "upward_wind_neutral   ", &
      "O_Density             ", &
      "O2_Density            ", &
      "N2_Density            "  &
      /)
  
  integer,                   parameter :: iHemiStart  = 0  ! Use Northern (0) and
  integer,                   parameter :: iHemiEnd    = 1  ! Southern (1) hemisphere (for debug only)

  ! -- mesh data
  integer                              :: numLocalNodes
  integer, dimension(:),   allocatable :: jmin
  integer, dimension(:),   allocatable :: jSouth
  integer, dimension(:,:), allocatable :: numLineNodes

  private

  public :: SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::SetServices()"

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
!---
!nm20161003: esmf timing lib
    real(ESMF_KIND_R8) :: beg_time, end_time
!---
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv01
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, specRoutine=InitializeData, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    
    ! attach specializing method(s)
!nm20161003 esmf timing library
    if (swESMFTime) then
      call ESMF_VMWtime(beg_time, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
    end if
    ! -- advance method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    if (swESMFTime) then
      call ESMF_VMWtime(end_time, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
      write(unit=9999,FMT=*)mype,"ModelAdvance endT=",(end_time-beg_time)
    end if

    ! -- finalize method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
     specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

  end subroutine SetServices
  
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeP0()"

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv02p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    
  end subroutine InitializeP0

  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeAdvertise()"

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! import fields from WAM
    call NUOPC_Advertise(importState, StandardNames=importFieldNames, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

  end subroutine InitializeAdvertise
  
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeRealize()"

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field) :: field
    type(ESMF_Mesh)  :: mesh

    integer          :: item
    logical          :: isConnected
    real(ESMF_KIND_R8), dimension(:), pointer :: fptr

    rc = ESMF_SUCCESS
    
    ! check if all required fields are connected
    item = 0
    isConnected = .true.
    do while (isConnected .and. (item < importFieldCount))
      item = item + 1
      isConnected = NUOPC_IsConnected(importState, &
        fieldName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    if (item < importFieldCount) then
      call ESMF_LogSetError(ESMF_RC_NOT_FOUND, &
        msg="Not all required fields are connected", &
        ESMF_CONTEXT, &
        rcToReturn=rc)
      return
    end if

    ! initialize IPE
    call myIPE_Init (clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! create 3D IPE mesh
    call IPEMeshCreate(gcomp, mesh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! write IPE mesh to VTK file if requested
    if (mesh_write > 0) then
      call ESMF_MeshWrite(mesh, trim(mesh_write_file), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! realize connected Fields in the importState
    do item = 1, importFieldCount
      field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
        name=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      fptr = 0._ESMF_KIND_R8
      call NUOPC_Realize(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeRealize::IPEMeshCreate()"

    !---------------------------------------------------------------------------
    ! SUBROUTINE: IPEMeshCreate
    !
    ! DESCRIPTION:
    !> \brief
    !> Build a 3D mesh from magnetic QD coordinates distributed in rectangular
    !! tiles in the (magnetic longitude, magnetic latitude) index space (mp,lp)
    !
    !> \details
    !> This subroutine selects all available heights between
    !! \c mesh_height_min and \c mesh_height_max, read as km from the
    !! \c &ipecap namelist in IPE.inp, then internally converted to m.
    !! The mesh's nodes and elements are created and indexed sequentially from
    !! the Northern to the Southern hemisphere, with nodes on the magnetic
    !! equator belonging only to the Northern hemisphere.
    !! The hemispheres are connected as hexahedral cells are built along the
    !! magnetic equator.
    !! Information about the horizontal domain decomposition is retrieved
    !! from the IPE module:
    !! \code
    !!   module_input_parameters
    !! \endcode
    !! as bounds of the local tile in magnetic index space (lps:lpe,mps:mpe).
    !! For debugging purposes, each hemisphere can be built independently by
    !! setting the module \c iHemiStart and \c iHemiEnd parameters to
    !! 0 (Northern hemisphere) or 1 (Southern hemisphere).
    !---------------------------------------------------------------------------

    subroutine IPEMeshCreate(gcomp, mesh, rc)

      ! --- input/output variables
      type(ESMF_GridComp)            :: gcomp
      type(ESMF_Mesh)                :: mesh
      integer, optional, intent(out) :: rc

      ! -- local variables
      integer :: localrc
      integer :: localPet, petCount, pet
      integer :: lp, lpp, lpu, mp, mpp, mpu, kp, kpp, kpe
      integer :: lpStart, lpOffset, mpStart, mpOffset, mHalo
      integer :: kpStart, kpEnd, kpStep, kpOffset
      integer :: numNodes, numHemiNodes
      integer :: numElems, numLongElems
      integer :: nCount
      integer :: id, in, is, im, i, j
      integer :: iHemi, jHemi

      integer, dimension(1)                :: jtop
      integer, dimension(:),   allocatable :: numLongNodes
      integer, dimension(:),   allocatable :: numLineElems
      integer, dimension(:),   allocatable :: nodeIds, nodeOwners
      integer, dimension(:),   allocatable :: elemIds, elemTypes, elemConn
      integer, dimension(:,:), allocatable :: idNodeOffset
      integer, dimension(:,:), allocatable :: petMap

      integer(ESMF_KIND_I4), dimension(:), allocatable :: localBounds, globalBounds
      real   (ESMF_KIND_R8), dimension(:), allocatable :: nodeCoords

      type(ESMF_VM) :: vm

      ! -- local parameters
      integer,               parameter :: logLevel    = 0  ! set log level
      integer, dimension(8), parameter :: iConn       = (/ 0,0,0,1,1,1,1,0 /)  ! connection array for hexahedron face (i,i+1) -> (i+2,i+3)...
      real(ESMF_KIND_R8),    parameter :: rad2deg     = 57.29577951308232087721_ESMF_KIND_R8
      real(ESMF_KIND_R8),    parameter :: earthRadius = 6371008.8_ESMF_KIND_R8 ! IUGG Earth Mean Radius (Moritz, 2000)

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      ! -- check that input parameters are acceptable
      ! -- verify that min and max mesh height are provided in the right order
      if (mesh_height_min >= mesh_height_max) then
        call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
        msg="Minimim height must be smaller than maximum height", &
        ESMF_CONTEXT, &
        rcToReturn=rc)
        return
      end if

      ! -- identify neighboring PETs
      ! -- get VM for current gridded component
      call ESMF_GridCompGet(gcomp, vm=vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- retrieve PET information
      call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- gather information about domain decomposition
      allocate(localBounds(4), globalBounds(4*petCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      localBounds = (/ mps, mpe, lps, lpe /)
      globalBounds = 0
      call ESMF_VMAllGather(vm, localBounds, globalBounds, 4, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- debug
      if (logLevel > 0) then
        if (localPet == 0) then
          write(6,'("====== IPEMeshCreate: tile bounds ======"/&
                   &"localPET     mps     mpe     lps     lpe"/&
                   &40("-"))')
          do i = 1, size(globalBounds), 4
            write(6,'(5i8)') i/4, globalBounds(i:i+3)
          end do
          write(6,'(40("="))')
        end if
      end if

      ! -- set bounds for local DE
      mHalo = min(NMP-mpe+mps-1,1)   ! magnetic longitude: halo size
      mpu = mpe + mHalo              ! magnetic longitude: periodic dimension
      lpu = min(lpe + 1, NLP)        ! magnetic latitude:  non-periodic dimension

      ! -- determine PETs that own neighboring DEs
      allocate(petMap(0:lpu-lpe,0:mpu-mpe), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- initialize PET map array
      petMap = -1

      i = 0
      do pet = 0, petCount-1
        do mpp = mpe, mpu
          mp = mod(mpp-1, NMP) + 1    ! returns periodic longitude index
          if ((globalBounds(i+1) <= mp) .and. (mp <= globalBounds(i+2))) then
            do lp = lpe, lpu
              if ((globalBounds(i+3) <= lp) .and. (lp <= globalBounds(i+4))) then
                petMap(lp-lpe,mpp-mpe) = pet
              end if
            end do
          end if
        end do
        i = i + 4
      end do

      deallocate(localBounds, globalBounds, stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- check if PET map is complete
      if (any(petMap < 0)) then
        call ESMF_LogSetError(ESMF_RC_OBJ_BAD, &
        msg="Failure to determine tile to PET mapping", &
        ESMF_CONTEXT, &
        rcToReturn=rc)
        return
      end if

      ! -- select grid points between min and max height
      allocate(jmin(NLP), &
        jSouth(NLP), &
        numLineNodes(NLP,0:1), &
        numLongNodes(0:1), &
        stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      jmin = 0
      jSouth = 0
      numLineNodes = 0
      numLongNodes = 0
      do lp = 1, NLP
        in = jmin_in(lp)
        is = jmax_is(lp)
        im = (in + is)/2    ! highest point (apex) along field line
        if (plasma_grid_Z(im,lp) >= mesh_height_min) then
          jtop = minloc(plasma_grid_Z(in:im,lp), plasma_grid_Z(in:im,lp) >= mesh_height_min)
          jmin(lp) = jtop(1) + in - 1
          is = is - jtop(1) + 1
          in = jmin(lp)
          if (plasma_grid_Z(im,lp) > mesh_height_max) then
            jtop = maxloc(plasma_grid_Z(in:im,lp), plasma_grid_Z(in:im,lp) <= mesh_height_max)
            jSouth(lp) = is - jtop(1) + 1
            numLineNodes(lp,0) = jtop(1)
          else
            jSouth(lp) = im + 1
            numLineNodes(lp,0) = im - in + 1
          end if
          numLineNodes(lp,1) = is - jSouth(lp) + 1
        end if
        numLongNodes = numLongNodes + numLineNodes(lp,:)
      end do

      numNodes      = (mpu-mps+1) * sum(numLineNodes(lps:lpu,iHemiStart:iHemiEnd))
      numLocalNodes = (mpe-mps+1) * sum(numLineNodes(lps:lpe,iHemiStart:iHemiEnd))

      ! -- build mesh nodes
      allocate(nodeIds(numNodes), &
        nodeOwners(numNodes), &
        nodeCoords(3*numNodes), &
        stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      nodeIds = 0
      nodeOwners = 0
      nodeCoords = 0._ESMF_KIND_R8

      i = 0
      nCount = 0
      do iHemi = iHemiStart, iHemiEnd
        mpStart = (iHemi-iHemiStart)*NMP*numLongNodes(0)
        lpStart = sum(numLineNodes(:lps-1,iHemi))
        do mpp = mps, mpu
          mp = mod(mpp-1, NMP) + 1
          j = max(mpp-mpe, 0)
          mpOffset = (mp-1) * numLongNodes(iHemi) + mpStart
          lpOffset = lpStart
          do lp = lps, lpu
            kpStep   = numLineNodes(lp,iHemi) - 1
            kpStart  = iHemi * kpStep + 1
            kpEnd    = (1 - iHemi) * kpStep + 1
            kpStep   = 1 - 2 * iHemi
            kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
            pet = petMap(max(lp-lpe,0),j)
            id = 0
            do kpp = kpStart, kpEnd, kpStep
              id = id + 1
              nCount = nCount + 1
              nodeIds(nCount) = id + lpOffset + mpOffset
              nodeOwners(nCount) = pet
              kp = kpp + kpOffset
              nodeCoords(i + 1) = rad2deg * plasma_grid_3d(kp,lp,mp,IGLON)
              nodeCoords(i + 2) = 90._ESMF_KIND_R8 - rad2deg * plasma_grid_3d(kp,lp,mp,IGCOLAT)
              nodeCoords(i + 3) = (1._ESMF_KIND_R8 + plasma_grid_Z(kp,lp)/earthRadius)
              i = i + 3
            end do
            lpOffset = lpOffset + numLineNodes(lp,iHemi)
          end do
        end do
      end do

      ! -- free up memory
!     deallocate(jmin, jSouth, petMap, stat=localrc)
      deallocate(petMap, stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- create mesh object
      mesh = ESMF_MeshCreate(parametricDim=3, spatialDim=3,  &
        coordSys=ESMF_COORDSYS_SPH_DEG, &
        rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- add nodes to mesh
      call ESMF_MeshAddNodes(mesh, nodeIds=nodeIds, &
        nodeCoords=nodeCoords, nodeOwners=nodeOwners, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- free up memory
      deallocate(nodeIds, nodeCoords, nodeOwners, stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- create mesh elements
      ! -- allocate work arrays
      allocate(numLineElems(NLP-1), idNodeOffset(lps:lpu,iHemiStart:iHemiEnd), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      numLineElems = 0
      numLongElems = 0
      do lp = 1, NLP - 1
        numLineElems(lp) = max(maxval(numLineNodes(lp:lp+1,0)) - 1, 0)
        numLongElems = numLongElems + numLineElems(lp)
      end do

      numElems = (iHemiEnd-iHemiStart+1)* (mpe-mps+1) * sum(numLineElems(lps:lpu-1))

      allocate(elemIds(numElems), &
        elemTypes(numElems), &
        elemConn(8*numElems), &
        stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      elemIds   = 0
      elemConn  = 0
      elemTypes = ESMF_MESHELEMTYPE_HEX

      ! -- recompute for local nodes
      numLongNodes = sum(numLineNodes(lps:lpu,:), dim=1)
      numHemiNodes = (mpu-mps+1) * numLongNodes(0)

      idNodeOffset = 0
      do lp = lps, lpu - 1
        idNodeOffset(lp + 1,:) = idNodeOffset(lp,:) + numLineNodes(lp,:)
      end do

      j = 0
      nCount = 0
      lpStart = sum(numLineElems(:lps-1))
      do iHemi = iHemiStart, iHemiEnd
        mpOffset = (mps-1) * numLongElems + (iHemi-iHemiStart) * NMP * numLongElems
        do mp = mps, mpe
          im = mp-mps+1
          im = im * min(1,NMP-im)
          lpOffset = lpStart
          do lp = lps, lpu - 1
            kpOffset = lpOffset + mpOffset
            do kp = 1, numLineElems(lp)
              nCount = nCount + 1
              elemIds(nCount) = kp + kpOffset
              do i = 1, 7, 2
                j = j + 1
                kpp = kp + iConn(i)
                lpp = lp + iConn(i+1)
                kpe = min(kpp, numLineNodes(lpp,0))
                jHemi = iHemi * max(min(numLineNodes(lpp,iHemi)-kpe,0)+1,0) !  select hemisphere for connected elements
                id = min(kpp,numLineNodes(lpp,jHemi)) + idNodeOffset(lpp,jHemi) + jHemi*(iHemi-iHemiStart)*numHemiNodes
                elemConn(j  ) = id + (mp-mps)*numLongNodes(jHemi)
                elemConn(j+4) = id + im * numLongNodes(jHemi)
              end do
              j = j + 4
            end do
            lpOffset = lpOffset + numLineElems(lp)
          end do
          mpOffset = mpOffset + numLongElems
        end do
      end do

      ! -- free up memory
!     deallocate(numLineElems, numLineNodes, idNodeOffset, numLongNodes, stat=localrc)
      deallocate(numLineElems, idNodeOffset, numLongNodes, stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- add elements to mesh
      call ESMF_MeshAddElements(mesh, elementIds=elemIds,&
             elementTypes=elemTypes, elementConn=elemConn, &
             rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

      ! -- free up memory
      deallocate(elemIds, elemTypes, elemConn, stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT, &
        rcToReturn=rc)) return

    end subroutine IPEMeshCreate

  end subroutine InitializeRealize
  
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::InitializeData()"

  subroutine InitializeData(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_State)     :: importState

    rc = ESMF_SUCCESS
    
    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
        
  end subroutine InitializeData

  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::ModelAdvance()"

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)   :: clock
    type(ESMF_State)   :: importState
    type(ESMF_Field)   :: field
    type(ESMF_Mesh)    :: mesh
!---
!nm20161003: esmf timing lib
    real(ESMF_KIND_R8) :: beg_time, end_time
!---
    integer :: item
    integer :: iHemi
    integer :: kp, kpp, kpStart, kpEnd, kpStep, kpOffset
    integer :: lp, mp, mpp
    integer :: nCount, numOwnedNodes, spatialDim
    integer :: localrc
    logical :: checkFields
    character(len=ESMF_MAXSTR) :: errmsg, value
    real(ESMF_KIND_R8)         :: dataValue
    real(ESMF_KIND_R8), dimension(:), pointer :: dataPtr, ownedNodeCoords

    real(ESMF_KIND_R8), parameter :: BAD_VALUE = -999._ESMF_KIND_R8

    rc = ESMF_SUCCESS

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! Look for CheckFields attribute to see if we should be
    ! checking the fields for errors
    call ESMF_AttributeGet(gcomp, &
      name="CheckFields", value=value, defaultvalue="false", &
      convention="NUOPC", purpose="Instance", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    ! convert value to checkfields
    checkFields=.false.
    if (trim(value)=="true") checkFields=.true.

    ! HERE IPE ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.

    ! -- print start time
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing IPE from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out
    
    ! -- import data
    ! -- initialize field to bad value
    wamField = BAD_VALUE

    do item = 1, importFieldCount
      ! --- retrieve field
      call ESMF_StateGet(importState, field=field, &
        itemName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out

      ! --- get field data
      nullify(dataPtr)
      call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out

      ! -- check if data size is consistent with local mesh size
      if (size(dataPtr) /= numLocalNodes) then
        errmsg = ""
        write(errmsg, '("Field: ",a,": retrieved data size (",i0,&
          &") does not match expected data size (",i0,")")') &
          trim(importFieldNames(item)), size(dataPtr), numLocalNodes
        call ESMF_LogSetError(ESMF_RC_PTR_BAD, &
          msg=errmsg, &
          ESMF_CONTEXT, &
          rcToReturn=rc)
        return  ! bail out
      end if

      nCount = 0
      do iHemi = iHemiStart, iHemiEnd
        do mp = mps, mpe
          do lp = lps, lpe
            kpStep   = numLineNodes(lp,iHemi) - 1
            kpStart  = iHemi * kpStep + 1
            kpEnd    = (1 - iHemi) * kpStep + 1
            kpStep   = 1 - 2 * iHemi
            kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
            do kpp = kpStart, kpEnd, kpStep
              nCount = nCount + 1
              kp = kpp + kpOffset
              wamField(kp,lp,mp,item) = dataPtr(nCount)
            end do
          end do
        end do
      end do

    end do

    ! -- check values of imported fields, if requested
    if (checkFields) then

      do item = 1, importFieldCount
        ! --- retrieve field
        call ESMF_StateGet(importState, field=field, &
          itemName=trim(importFieldNames(item)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          ESMF_CONTEXT)) &
          return  ! bail out

        ! --- get field data
        nullify(dataPtr)
        call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          ESMF_CONTEXT)) &
          return  ! bail out

        if (.not.associated(ownedNodeCoords)) then
          ! -- retrieve node coordinates
          call ESMF_FieldGet(field, mesh=mesh, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            ESMF_CONTEXT)) &
            return  ! bail out

          call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedNodes=numOwnedNodes, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            ESMF_CONTEXT)) &
            return  ! bail out

          allocate(ownedNodeCoords(spatialDim * numOwnedNodes), stat=localrc)
          if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            ESMF_CONTEXT, &
            rcToReturn=rc)) return

          call ESMF_MeshGet(mesh, ownedNodeCoords=ownedNodeCoords, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            ESMF_CONTEXT)) &
            return  ! bail out

        end if

        nCount = 0
        do iHemi = iHemiStart, iHemiEnd
          do mp = mps, mpe
            do lp = lps, lpe
              kpStep   = numLineNodes(lp,iHemi) - 1
              kpStart  = iHemi * kpStep + 1
              kpEnd    = (1 - iHemi) * kpStep + 1
              kpStep   = 1 - 2 * iHemi
              kpOffset = (1 - iHemi) * jmin(lp) + iHemi * jSouth(lp) - 1
              do kpp = kpStart, kpEnd, kpStep
                kp = kpp + kpOffset
                nCount = nCount + 1
                dataValue = dataPtr(nCount)
                errmsg = ""
                if (dataValue == BAD_VALUE) then
                   ! -- field array contains unfilled points
                   errmsg = "BAD_VALUE"
                else if (abs(dataValue - BAD_VALUE) < 0.1_ESMF_KIND_R8) then
                   ! -- field array contains approximate unfilled points
                   errmsg = "approx BAD_VALUE"
                else if (.not.(dataValue == dataValue)) then
                   ! -- field array contains NaN
                   errmsg = "NaN"
                end if
                if (len_trim(errmsg) > 0) then
                  write(6,'("ERROR: ",a,": ",a," found (",g20.8,") at (mp,lp,kp) = (",2(i0,","),i0,")"' &
                    //'," coord: ",2(f12.6,","),g20.8)') &
                    trim(importFieldNames(item)), trim(errmsg), mp, lp, kp, &
                    ownedNodeCoords(nCount:nCount+spatialDim-1)
                endif
                errmsg = ""
                select case (trim(importFieldNames(item)))
                  case ("temp_neutral", "O_Density", "O2_Density", "N2_Density")
                    ! -- these fields must contain positive values
                    if (dataValue <= 0._ESMF_KIND_R8) errmsg = "value <= 0.0"
                end select
                if (len_trim(errmsg) > 0) then
                  write(6,'("ERROR: ",a,": ",a," found (",g20.8,") at (mp,lp,kp) = (",2(i0,","),i0,")"' &
                    //'," coord: ",2(f12.6,","),g20.8)') &
                    trim(importFieldNames(item)), trim(errmsg), mp, lp, kp, &
                    ownedNodeCoords(nCount:nCount+spatialDim-1)
                endif
              end do
            end do
          end do
        end do

      end do

      if (associated(ownedNodeCoords)) then
        deallocate(ownedNodeCoords, stat=localrc)
        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          ESMF_CONTEXT, &
          rcToReturn=rc)) return
      end if

    end if

    nullify(dataPtr)

    ! -- advance IPE model
!nm20161003 esmf timing library
    if (swESMFTime) then
      call ESMF_VMWtime(beg_time, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
    end if

    call update_IPE(clock, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

    if (swESMFTime) then
      call ESMF_VMWtime(end_time, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        ESMF_CONTEXT)) &
        return  ! bail out
      write(unit=9999,FMT=*)mype,"update_IPE endT=",(end_time-beg_time)
    end if

    ! -- print stop time
    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT)) &
      return  ! bail out

  end subroutine ModelAdvance
 
  !-----------------------------------------------------------------------------

#undef  ESMF_METHOD
#define ESMF_METHOD "IPECap::Finalize()"

  subroutine Finalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin

    ! -- free up memory
    deallocate(jmin, jSouth, numLineNodes, stat=localrc)
    if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      ESMF_CONTEXT, &
      rcToReturn=rc)) return

    ! -- finalize IPE model
    call finalize_IPE(gcomp, rc)

  end subroutine Finalize

end module ipeCap
