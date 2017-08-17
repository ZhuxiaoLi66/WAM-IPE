module DATAIPE

  !-----------------------------------------------------------------------------
  ! DATAIPE: Data IPE Component
  !
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance

#ifdef ESMF_NETCDF
  use netcdf
#endif

  implicit none
  
  private

  integer, parameter :: MAXNAMELEN = 128

  ! private internal state to keep instance data
  type InternalStateStruct
    integer, pointer :: maxlevels(:), totallats(:)
    integer :: ipedims(3)
    integer :: startlon, countlon, highestlevel
    integer :: localnodes
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

    ! import fields from WAM
    call NUOPC_Advertise(importState, StandardNames=(/ &
      "northward_wind_neutral         ", &
      "eastward_wind_neutral          ", &
    !  "upward_wind_neutral            ", &
      "temp_neutral                   " &
    !  "O_Density                      ", &
    !  "O2_Density                     ", &
    !  "N2_Density                     ", &
    !  "NO_Density                     ", &
    !  "N4S_Density                    ", &
    !  "N2D_Density                    ", &
    !  "H_Density                      ", &
    !  "He_Density                     ", &
      /), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set Component name so it becomes identifiable in the output
    call ESMF_GridCompSet(gcomp, name="DATAIPE", rc=rc)
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
    type(ESMF_Mesh)      :: mesh
    logical              :: freeflag
  
    rc = ESMF_SUCCESS
    
    ! create a 3D IPE Mesh using definition in ipe3dgrid.nc
    call createIPEGrid(gcomp, mesh, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! realize connected Fields in the importState
    call realizeConnectedFields(importState, mesh=mesh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine createIPEGrid(gcomp, mesh, rc)
      type(ESMF_GridComp) :: gcomp
      type(ESMF_Mesh) :: mesh
      integer :: rc
      
      type(ESMF_VM) :: vm
      integer :: PetNo, PetCnt
      character(len=MAXNAMELEN) :: ipefilename
      integer :: nc1, status
      integer :: varid, ndims
      integer :: dimids(3), ipedims(3)
      real(ESMF_KIND_R8), pointer    :: ipelon(:,:,:), ipelat(:,:,:), ipehgt(:,:)
      integer(ESMF_KIND_I4), pointer :: totalheight(:)
      integer(ESMF_KIND_I4), pointer :: elementIds(:), elementTypes(:), elementConn(:)
      integer(ESMF_KIND_I4), pointer :: nodeIds(:), nodeOwners(:)
      real(ESMF_KIND_R8), pointer    :: nodeCoords(:)
      integer(ESMF_KIND_I4), pointer :: southind(:), northind(:), totallats(:), gap(:)
      integer(ESMF_KIND_I4), pointer :: maxlevs(:)
      integer(ESMF_KIND_I4), pointer :: elementCnt(:), nodeCnt(:), sendbuf(:), recvbuf(:)
      real(ESMF_KIND_R8), pointer :: conntbl(:), globalCoords(:,:), fptr2d(:,:), fptr1d(:)
      integer, pointer :: allCounts(:), connectbase(:)
      integer :: count, remind, start, halo, neightbor, base, increment
      integer :: totalnodes, localnodes, totalelements, globalTotal, globalTotalelmt
      real(ESMF_KIND_R8) :: maxheight, lon, lat, hgt
      integer :: i, j, k, ii, jj, count1, count3, count8, localcount
      integer :: startid, nodestartid
      integer :: neighbor, elmtcount, countup, diff, lastlat, base1, save
      logical :: even
      type(InternalState)  :: is
 
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

      ipefilename = 'data/ipe3dgrid.nc'
      maxheight = 800.0

#ifdef ESMF_NETCDF
      ! Read in the ipe grid 
      status = nf90_open(path=ipefilename, mode=nf90_nowrite, ncid=nc1)
      call CheckNCError(status, ipefilename)
      status = nf90_inq_varid(nc1,'lons', varid)
      call CheckNCError(status, 'lons')
      status = nf90_inquire_variable(nc1, varid, ndims=ndims, dimids = dimids)
      call CheckNCError(status, 'lons')
      status = nf90_inquire_dimension(nc1,dimids(1), len=ipedims(1))
      call CheckNCError(status, 'lons 1st dimension')
      status = nf90_inquire_dimension(nc1,dimids(2), len=ipedims(2))
      call CheckNCError(status, 'lons 2nd dimension')
      status = nf90_inquire_dimension(nc1,dimids(3), len=ipedims(3))
      call CheckNCError(status, 'lons 3rd dimension')

      ! IPE dimension order:  nlevels, lats and lons (1115, 170, 80)
      allocate(ipelon(ipedims(1), ipedims(2), ipedims(3)), &
  	   ipelat(ipedims(1), ipedims(2), ipedims(3)), &
	   ipehgt(ipedims(1), ipedims(2)))

      allocate(maxlevs(ipedims(2)))
      status = nf90_get_var(nc1, varid, ipelon)
      call CheckNCError(status, 'lons')
      status = nf90_inq_varid(nc1,'lats', varid)
      call CheckNCError(status, 'lats')
      status = nf90_get_var(nc1, varid, ipelat)
      call CheckNCError(status, 'lats')
      status = nf90_inq_varid(nc1,'height', varid)
      call CheckNCError(status, 'height')
      status = nf90_get_var(nc1, varid, ipehgt)
      call CheckNCError(status, 'height')
      status = nf90_inq_varid(nc1,'maxLevels', varid)
      call CheckNCError(status, 'maxLevels')
      status = nf90_get_var(nc1, varid, maxlevs)
      call CheckNCError(status, 'maxLevels')
      status = nf90_close(nc1)
      call CheckNCError(status, ipefilename)

      ! Construct IPE Mesh based on block decomposition
      ! The decomposition is stored in a namelist called SMSnamelist, 
      ! the keyword is process_layout and it is a 2D decomposition of n,m
      ! where nxm has to be equal to the total number of processors and the
      ! latitude lines (170) will be decomposed to n blocks and longitude lines (80)
      ! to m blocks.  Usually we keep n=1 if total number of processors are < 80
  
      ! Let's hardcode the docomposition to 1xPetCnt if PetCnt < 80 for now
      ! 
      ! Also we only build the mesh with a subset of the height levels that
      ! intersect with the WAM grid.  Currently, we set it to 90KM to 800KM
      !
      ! Step 1: calculate number of longitude lines in the local processor
      if (ipedims(3) > PetCnt) then
        count = ipedims(3)/PetCnt
        remind = ipedims(3) - PetCnt*count
        if (PetNo < remind) then
           count = count+1
           start = count*PetNo+1
        else
           start = count*PetNo+remind+1
        endif
      else 
        call ESMF_LogSetError(rcToCheck=ESMF_FAILURE, &
	      msg="Total number of proessor exceed the longitude dimension", &
              line=__LINE__, file=__FILE__, &
	      rcToReturn=rc)
        return
      endif
      !print *, PetNo, ' start from ', start, ' with ', count, ' columns'

      !Find the total node count and cell count in the local PE
      ! hardcode the maxheight to 800km, find how many height levels are within the range
      ! for each latitute tube
      ! Need to separate it into south hemisphere and north hemisphere 
      allocate(totalheight(ipedims(2)),northind(ipedims(2)), southind(ipedims(2)))
      allocate(gap(ipedims(2)))
      do j=1,ipedims(2)
        do i=1,(maxlevs(j)-1)/2
          if (ipehgt(i,j)> maxheight) then
             northind(j)=i-1
	     exit
          endif
        enddo
        if (i==(maxlevs(j)/2+1)) northind(j)=i-1 
        do i=maxlevs(j),(maxlevs(j)+1)/2+1,-1
          if (ipehgt(i,j)> maxheight) then
             southind(j)=i+1
	     exit
          endif
        enddo
        if (i==(maxlevs(j)+1)/2) southind(j)=i+1 
        if (southind(j)==northind(j)) southind(j)=southind(j)+1
        totalheight(j) = northind(j)+maxlevs(j)-southind(j)+1
        gap(j) = southind(j)-northind(j)+1   
        !if (PetNo == 0) then
        !  print *, 'min/max levels for ', j, 'th field line are', northind(j), southind(j)
        !endif
      enddo

      ! The layers in between northind(1) and and southind(1) will not be included in 
      ! the mesh, will adjust the node index by taking out the gap layers to reduce the number
      ! of missing indices.
      maxheight=totalheight(1)
      ! Now construct the node arrays for Mesh creation
      ! First count how many nodes used in the local PE, each longitude has the same number of
      ! nodes, need to add one halo column
      totalnodes = 0
      localnodes = 0
      if (PetCnt > 1) then
         do j=1,ipedims(2)
            totalnodes = totalheight(j)*(count+1)+totalnodes  
            localnodes = totalheight(j)*count+localnodes
         enddo
      else 
         do j=1,ipedims(2)
            totalnodes = totalheight(j)*count+totalnodes  
         enddo
         localnodes=totalnodes
      endif
      allocate(nodeIds(totalnodes), nodeCoords(totalnodes*3),nodeOwners(totalnodes))

      ! find out the starting global id of the local node
      allocate(nodeCnt(PetCnt*2),sendbuf(2))
      sendbuf(1)=localnodes
      sendbuf(2)=count
      call ESMF_VMAllGather(vm, sendbuf, nodeCnt, 2, rc=rc)
      ! find the starting globalNodeID
      ! and global total nodes
      globalTotal = 0
      startid=0
      do i=1,PetNo
         startid=startid+nodeCnt(i*2-1)
      enddo
      do i=1,PetCnt
         globalTotal=globalTotal+nodeCnt(i*2-1)
      enddo
      deallocate(sendbuf)
      ! ipetotalnodes = globalTotal !??????? never used

      ! count number of latitudes and local elements at each height 
      totalelements=0
      allocate(totallats(northind(1)))
      do j=1,northind(1)
         do i=1,ipedims(2)
            if (northind(i)<j) exit
         enddo
         totallats(j)=i-1
         ! do not count the toppest level
         if (j<northind(1)) then
            totalelements=totalelements + count*(totallats(j)-1)
         endif
      enddo
      ! totalelements include one extra row connecting the north and south hemisphere
      totalelements = totalelements*2 + count*(northind(1)-1)

      ! Do the halo longitude
      if (PetNo == PetCnt-1) then
         halo=1
         neighbor=0
         base=1
         increment=nodeCnt(2)
      else
         halo=start+count
         neighbor=PetNo+1
         base=startid+localnodes+1
         increment=nodeCnt(PetNo*2+2)
      endif
      deallocate(nodeCnt)
      !print *, PetNo, ' halo, neighbor, base, increment ', halo, neighbor, base, increment  
      count1 = 1
      localcount = 1
      count3 = 1

      !----------------------------
      ! The local node will be arranged by longitude first, followed by latitude, height and
      ! and hemisphere. There will be fewer latitude columns at higher height. 
      ! North hemiphere first, then South hemiphere.  Two hemisphere should have the same number of nodes.
      ! The global ID is the sequential index of the 3D array coord(height,lat, lon).  But, the nodes at the
      ! South Hemisphere has adjusted heights by deducting the number of layers of the gap between southind(:)
      ! and northind(:)
      !-----------------------------
      ! North Hemisphere

      do i=1,northind(1)
         do j=1,totallats(i)
            do k=start,count+start-1
               !nodeIds(count1)=i+(j-1)*maxheight+(k-1)*j*maxheight
               nodeIds(count1)=localcount+startid
	       lon=ipelon(i,j,k)
               lat=ipelat(i,j,k)
               hgt=ipehgt(i,j)
               call convert2Cart(lon,lat,hgt,nodeCoords(count3:count3+2))
               nodeOwners(count1)=PetNo
               count1=count1+1
               localcount=localcount+1
               count3=count3+3
            enddo
            if (PetCnt > 1) then
               ! Add hallo
               !nodeIds(count1)=i+(j-1)*maxheight+(halo-1)*j*maxheight
               ! Need to find out the globalID of the neightboring node
               nodeIds(count1)=base
               base=base+increment
               lon=ipelon(i,j,halo)
               lat=ipelat(i,j,halo)
               hgt=ipehgt(i,j)
               call convert2Cart(lon,lat,hgt,nodeCoords(count3:count3+2))
               nodeOwners(count1)=neighbor
               count1=count1+1
               count3=count3+3
            endif
         enddo
      enddo
   ! South Hemisphere
  do ii=1,northind(1)
    do j=1,totallats(ii)
      !adjust the height index based on the maxlevel of the latitude
      i=maxlevs(j)-ii+1
      do k=start,count+start-1
         !nodeIds(count1)=(i-gap(j))+(j-1)*maxheight+(k-1)*j*maxheight
         nodeIds(count1)=localcount+startid
	 lon=ipelon(i,j,k)
         lat=ipelat(i,j,k)
         hgt=ipehgt(i,j)
         call convert2Cart(lon,lat,hgt,nodeCoords(count3:count3+2))
         nodeOwners(count1)=PetNo
         count1=count1+1
         localcount=localcount+1
         count3=count3+3
       enddo
       ! Add hallo
       if (PetCnt > 1) then
       !nodeIds(count1)=(i-gap(j))+(j-1)*maxheight+(halo-1)*j*maxheight
       nodeIds(count1)=base
       base=base+increment
       lon=ipelon(i,j,halo)
       lat=ipelat(i,j,halo)
       hgt=ipehgt(i,j)
       call convert2Cart(lon,lat,hgt,nodeCoords(count3:count3+2))
       nodeOwners(count1)=neighbor
       count1=count1+1
       count3=count3+3
       endif 
     enddo
  enddo
  nodestartid = startid+1
  !print *, PetNo, 'totalnodes, localnodes ', count1-1, localcount-1
  if (count1-1 /= totalnodes) then
     print *, 'count mismatch:', count1-1, totalnodes
  endif

  ! Create mesh and add node
#ifdef USE_CART3D_COORDSYS
  mesh = ESMF_MeshCreate(3,3,coordSys=ESMF_COORDSYS_CART, rc=rc)
#else
  mesh = ESMF_MeshCreate(3,3,coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
#endif
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  call ESMF_MeshAddNodes(mesh, nodeIds, nodeCoords, nodeOwners, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  ! Now construct Mesh elements
  ! Mesh are separated at the megnetic equator, so the mesh elements are two less than the 
  ! total latitude points at each layer
  ! First find out how many elements there are per processor
  ! Start from the bottom, calculate the elements in the north hemisphere, then double it
  ! The number of cells at each height is determined by how many grid points at that layer

  allocate(elementIds(totalelements), elementTypes(totalelements), &
  	   elementConn(totalelements*8))

  elementTypes(:)=ESMF_MESHELEMTYPE_HEX
  
  ! find out the starting global id of the local element
  allocate(elementCnt(PetCnt),sendbuf(1))
  sendbuf(1)=totalelements
  call ESMF_VMAllGather(vm, sendbuf, elementCnt, 1, rc=rc)
  ! find the starting elementID
  startid=0
  do i=1,PetNo
    startid=startid+elementCnt(i)
  enddo
  globalTotalelmt=0
  do i=1,PetCnt
    globalTotalelmt=globalTotalelmt+elementCnt(i)
  enddo
  deallocate(elementCnt, sendbuf)

  !print *, PetNo, 'Total Elements: ', totalelements, 'starting at', startid

  ! Build the local elementConn table using local node indices
  count1=1
  count8=1
  ! base is the starting local index of the node at level i
  base=0
  ! Do the bottom rectangle of each element first
  ! North hemisphere first

  !!!!!!!!!!!!!!!!!!!!!!
  ! For sequential case, make count one less since no neighbor column is needed
  if (PetCnt == 1) count=count-1
  !!!!!!!!!!!!!!!!!!!!!!
  allocate(connectbase(northind(1)))
  do i=1,northind(1)-1
    do j=1,totallats(i)-1
      do k=1,count
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*j+k
         elementConn(count8+1)=base+(count+1)*j+k+1
         elementConn(count8+2)=base+(count+1)*(j-1)+k+1
         elementConn(count8+3)=base+(count+1)*(j-1)+k
	 count1=count1+1
	 count8=count8+8
      enddo      
      if (PetCnt ==1) then
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*j+k
         elementConn(count8+1)=base+(count+1)*j+1
         elementConn(count8+2)=base+(count+1)*(j-1)+1
         elementConn(count8+3)=base+(count+1)*(j-1)+k
	 count1=count1+1
	 count8=count8+8
      endif
    enddo
    base=base+totallats(i)*(count+1)
    connectbase(i)=base-(count+1)
  enddo
  ! add the top level to base
  base = base+totallats(i)*(count+1)
  connectbase(i)=base-(count+1)
  ! South Hemisphere
  do i=1,northind(1)-1
    do j=1,totallats(i)-1
      do k=1,count
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*(j-1)+k
         elementConn(count8+1)=base+(count+1)*(j-1)+k+1
         elementConn(count8+2)=base+(count+1)*j+k+1
         elementConn(count8+3)=base+(count+1)*j+k
	 count1=count1+1
	 count8=count8+8
      enddo      
      ! If PetCnt==1, connect the last column back to the first column
      if (PetCnt == 1) then
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*(j-1)+k
         elementConn(count8+1)=base+(count+1)*(j-1)+1
         elementConn(count8+2)=base+(count+1)*j+1
         elementConn(count8+3)=base+(count+1)*j+k
	 count1=count1+1
	 count8=count8+8
      endif
    enddo
    ! Add one extra row of cells to connect the north and south hemisphere
    ! Need to figure out the local node id of the last row in North Hemisphere
    do k=1,count
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*(totallats(i)-1)+k
         elementConn(count8+1)=base+(count+1)*(totallats(i)-1)+k+1
         elementConn(count8+2)=connectbase(i)+k+1
         elementConn(count8+3)=connectbase(i)+k
	 count1=count1+1
	 count8=count8+8
    enddo      
    ! If PetCnt==1, connect the last column back to the first column
    if (PetCnt == 1) then
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*(totallats(i)-1)+k
         elementConn(count8+1)=base+(count+1)*(totallats(i)-1)+1
         elementConn(count8+2)=connectbase(i)+1
         elementConn(count8+3)=connectbase(i)+k
	 count1=count1+1
	 count8=count8+8
    endif

    base=base+totallats(i)*(count+1)
  enddo

  ! Check the total number of elements match with totalelements
  if (count1-1 /= totalelements) then
     print *, PetNo, ' total elements mismatch ', count1-1, totalelements
  endif
  ! Now building the top face (quad) of the cube, note if the top level has fewer latitudes, we need 
  ! to collape some the the top into on single line.
  count8=1
  base = 0
  ! first time for North hemisphere, second time for South hemisphere, should be identical
  if (PetCnt==1) then 
     elmtcount=count+1
  else
     elmtcount=count
  endif
  do ii=1,2  
    do i=1, northind(1)-2
      ! base is the starting local index to the elementConn table for the level i+1
      if (ii==1) then
         base = base+count*(totallats(i)-1)*8
      else
         ! For south hemisphere, we have one extra row to connect the south to north
         base = base+count*totallats(i)*8
      endif
      ! countup is the sequential index to elementConn at level i+1 (starting at 1)
      countup=1
      if (totallats(i)==totallats(i+1)) then
        do j=1,totallats(i)-1
          do k=1,elmtcount
            elementConn(count8+4)=elementConn(base+countup)
            elementConn(count8+5)=elementConn(base+countup+1)
            elementConn(count8+6)=elementConn(base+countup+2)
            elementConn(count8+7)=elementConn(base+countup+3)
            count8=count8+8
  	    countup=countup+8
          enddo
        enddo
	if (ii==2) then
          do k=1,elmtcount
            elementConn(count8+4)=elementConn(base+countup)
            elementConn(count8+5)=elementConn(base+countup+1)
            elementConn(count8+6)=elementConn(base+countup+2)
            elementConn(count8+7)=elementConn(base+countup+3)
            count8=count8+8
  	    countup=countup+8
          enddo
        endif
      else ! the top level has fewer rows than the bottom layer
        diff = totallats(i)-totallats(i+1)
        ! put the prism cells close the the last few latitudes (close to the equator)
        lastlat=totallats(i)-diff*2
        do j=1,lastlat-1
          do k=1,elmtcount
            elementConn(count8+4)=elementConn(base+countup)
            elementConn(count8+5)=elementConn(base+countup+1)
            elementConn(count8+6)=elementConn(base+countup+2)
            elementConn(count8+7)=elementConn(base+countup+3)
            count8=count8+8
	    countup=countup+8
          enddo
        enddo
        even=.TRUE.
        jj=lastlat
        do j=lastlat,totallats(i)-1
  	  if (even) then
            save = countup
            do k=1,elmtcount
              elementConn(count8+4)=elementConn(base+countup)
              elementConn(count8+5)=elementConn(base+countup+1)
              elementConn(count8+6)=elementConn(base+countup+1)
              elementConn(count8+7)=elementConn(base+countup)
              count8=count8+8
	      countup=countup+8
	    enddo
	    even = .FALSE.
	    countup = save            
          else
            do k=1,elmtcount
              elementConn(count8+4)=elementConn(base+countup)
              elementConn(count8+5)=elementConn(base+countup+1)
              elementConn(count8+6)=elementConn(base+countup+2)
              elementConn(count8+7)=elementConn(base+countup+3)
              count8=count8+8
 	      countup=countup+8
            enddo
            jj=jj+1
	    even = .TRUE.
	  endif
        enddo
	if (ii==2) then
          do k=1,elmtcount
            elementConn(count8+4)=elementConn(base+countup)
            elementConn(count8+5)=elementConn(base+countup+1)
            elementConn(count8+6)=elementConn(base+countup+2)
            elementConn(count8+7)=elementConn(base+countup+3)
            count8=count8+8
  	    countup=countup+8
          enddo
        endif
        !print *, PetNo, 'Level ', i, 'Top level elements ', countup/8
      endif 	    
   enddo

  ! Now the top row of the cube i=northind(1)-1
  ! First find out the startind nodeind of the top level 
  ! base1 is the starting index of the local node id at the top level 
  if (ii==1) base1=0
  do j=1,northind(1)-1
     base1=base1+totallats(j)*(count+1)
  enddo
  ! If the top two levels have the same number of latitude points, i=northind(1)-1
  if (totallats(i)==totallats(i+1)) then
      if (ii==1) then 
        do j=1,totallats(i)-1
          do k=1,count
             elementConn(count8+4)=base1+(count+1)*j+k
             elementConn(count8+5)=base1+(count+1)*j+k+1
             elementConn(count8+6)=base1+(count+1)*(j-1)+k+1
             elementConn(count8+7)=base1+(count+1)*(j-1)+k
             count8=count8+8
           enddo
	   if (PetCnt==1) then
             elementConn(count8+4)=base1+(count+1)*j+k
             elementConn(count8+5)=base1+(count+1)*j+1
             elementConn(count8+6)=base1+(count+1)*(j-1)+1
             elementConn(count8+7)=base1+(count+1)*(j-1)+k
             count8=count8+8
           endif
        enddo
	base1 = base1+totallats(i+1)*(count+1)
     else ! (ii==2)
        do j=1,totallats(i)-1
          do k=1,count
             elementConn(count8+4)=base1+(count+1)*(j-1)+k
             elementConn(count8+5)=base1+(count+1)*(j-1)+k+1
             elementConn(count8+6)=base1+(count+1)*j+k+1
             elementConn(count8+7)=base1+(count+1)*j+k
             count8=count8+8
           enddo
           if (PetCnt==1) then
             elementConn(count8+4)=base1+(count+1)*(j-1)+k
             elementConn(count8+5)=base1+(count+1)*(j-1)+1
             elementConn(count8+6)=base1+(count+1)*j+1
             elementConn(count8+7)=base1+(count+1)*j+k
             count8=count8+8
           endif
        enddo
        ! The extra row of cells
        do k=1,count
             elementConn(count8+4)=base1+(count+1)*(totallats(i)-1)+k
             elementConn(count8+5)=base1+(count+1)*(totallats(i)-1)+k+1
             elementConn(count8+6)=base1-(count+1)+k+1
             elementConn(count8+7)=base1-(count+1)+k
  	     count8=count8+8
        enddo      
        ! If PetCnt==1, connect the last column back to the first column
        if (PetCnt == 1) then
             elementConn(count8+4)=base1+(count+1)*(totallats(i)-1)+k
             elementConn(count8+5)=base1+(count+1)*(totallats(i)-1)+1
             elementConn(count8+6)=base1-(count+1)+1
             elementConn(count8+7)=base1-(count+1)+k
  	     count8=count8+8
        endif
      endif ! ii
    else ! totallats(i) /= totallats(i+1)
      diff = totallats(i)-totallats(i+1)
      ! put the prism cells close the the last few latitudes (close to the equator)
      lastlat=totallats(i)-diff*2
      ! the first 1 to lastlat rows will be cubes
      if (ii==1) then 
        do j=1,lastlat-1
          do k=1,count
             elementConn(count8+4)=base1+(count+1)*j+k
             elementConn(count8+5)=base1+(count+1)*j+k+1
             elementConn(count8+6)=base1+(count+1)*(j-1)+k+1
             elementConn(count8+7)=base1+(count+1)*(j-1)+k
             count8=count8+8
           enddo
	   if (PetCnt==1) then
             elementConn(count8+4)=base1+(count+1)*j+k
             elementConn(count8+5)=base1+(count+1)*j+1
             elementConn(count8+6)=base1+(count+1)*(j-1)+1
             elementConn(count8+7)=base1+(count+1)*(j-1)+k
             count8=count8+8
	   endif
        enddo
      else ! (ii==2)
        do j=1,lastlat-1
          do k=1,count
             elementConn(count8+4)=base1+(count+1)*(j-1)+k
             elementConn(count8+5)=base1+(count+1)*(j-1)+k+1
             elementConn(count8+6)=base1+(count+1)*j+k+1
             elementConn(count8+7)=base1+(count+1)*j+k
             count8=count8+8
           enddo
	   if (PetCnt==1) then
             elementConn(count8+4)=base1+(count+1)*(j-1)+k
             elementConn(count8+5)=base1+(count+1)*(j-1)+1
             elementConn(count8+6)=base1+(count+1)*j+1
             elementConn(count8+7)=base1+(count+1)*j+k
             count8=count8+8
	   endif
        enddo
      endif ! ii
      even=.TRUE.
      jj=lastlat
      do j=lastlat,totallats(i)-1
	if (even) then ! collapse the top into a line
          do k=1,count
            elementConn(count8+4)=base1+(count+1)*(jj-1)+k
            elementConn(count8+5)=base1+(count+1)*(jj-1)+k+1
            elementConn(count8+6)=base1+(count+1)*(jj-1)+k+1
            elementConn(count8+7)=base1+(count+1)*(jj-1)+k
            count8=count8+8
	  enddo
          if (PetCnt==1) then
            elementConn(count8+4)=base1+(count+1)*(jj-1)+k
            elementConn(count8+5)=base1+(count+1)*(jj-1)+1
            elementConn(count8+6)=base1+(count+1)*(jj-1)+1
            elementConn(count8+7)=base1+(count+1)*(jj-1)+k
            count8=count8+8
	  endif
	  even = .FALSE.            
        else ! odd, the top is a quad
	   if (ii==1) then
             do k=1,count
               elementConn(count8+4)=base1+(count+1)*jj+k
               elementConn(count8+5)=base1+(count+1)*jj+k+1
               elementConn(count8+6)=base1+(count+1)*(jj-1)+k+1
               elementConn(count8+7)=base1+(count+1)*(jj-1)+k
               count8=count8+8
             enddo
	     if (PetCnt==1) then
               elementConn(count8+4)=base1+(count+1)*jj+k
               elementConn(count8+5)=base1+(count+1)*jj+1
               elementConn(count8+6)=base1+(count+1)*(jj-1)+1
               elementConn(count8+7)=base1+(count+1)*(jj-1)+k
               count8=count8+8
             endif	       
	     jj=jj+1
           else ! (ii==2)
             do k=1,count
               elementConn(count8+4)=base1+(count+1)*(jj-1)+k
               elementConn(count8+5)=base1+(count+1)*(jj-1)+k+1
               elementConn(count8+6)=base1+(count+1)*jj+k+1
               elementConn(count8+7)=base1+(count+1)*jj+k
               count8=count8+8
             enddo
             if (PetCnt==1) then
               elementConn(count8+4)=base1+(count+1)*(jj-1)+k
               elementConn(count8+5)=base1+(count+1)*(jj-1)+1
               elementConn(count8+6)=base1+(count+1)*jj+1
               elementConn(count8+7)=base1+(count+1)*jj+k
               count8=count8+8
	     endif
             jj=jj+1
	   endif ! ii
	   even = .TRUE.
         endif ! even
      enddo ! j=lastlat,totallats(i)-1
      ! The extra row of cells
      if (ii==2) then
        do k=1,count
           elementConn(count8+4)=base1+(count+1)*(totallats(i+1)-1)+k
           elementConn(count8+5)=base1+(count+1)*(totallats(i+1)-1)+k+1
           elementConn(count8+6)=connectbase(northind(1))+k+1
           elementConn(count8+7)=connectbase(northind(1))+k
           count8=count8+8
        enddo      
        ! If PetCnt==1, connect the last column back to the first column
        if (PetCnt == 1) then
           elementConn(count8+4)=base1+(count+1)*(totallats(i+1)-1)+k
           elementConn(count8+5)=base1+(count+1)*(totallats(i+1)-1)+1
           elementConn(count8+6)=connectbase(northind(1))+1
           elementConn(count8+7)=connectbase(northind(1))+k
           count8=count8+8
        endif
      endif
      base1 = base1+totallats(i+1)*(count+1)
    endif ! (totallats(i)==totallats(i+1)) 

    ! Need to add the top level to base
    base = base + (totallats(northind(1)-1)-1)*count*8
  enddo ! ii=1,2
  
  do i=1,totalelements*8
    if (elementConn(i) > totalnodes .OR. elementConn(i)==0) then
       print *, PetNo, ' node index out of range', i/8, i, elementConn(i)
    endif
  enddo
  call ESMF_VMBarrier(vm) 
  call ESMF_MeshAddElements(mesh, elementIds, elementTypes, elementConn,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  is%wrap%ipedims = ipedims
  is%wrap%maxlevels => maxlevs
  is%wrap%highestlevel = northind(1)
  is%wrap%totallats => totallats
  is%wrap%localnodes = localnodes
  is%wrap%startlon = start
  is%wrap%countlon = count

  deallocate(southind, northind, gap)
  deallocate(nodeIds, nodeOwners, nodeCoords)
  deallocate(elementIds, elementTypes, elementConn)
  deallocate(ipelon, ipelat, ipehgt)
  deallocate(connectbase)
  return

#else
    call ESMF_LogSetError(ESMF_RC_LIB_NOT_PRESENT, & 
                 msg="- ESMF_NETCDF not defined when lib was compiled") 
    return
#endif
    end subroutine createIPEgrid

    subroutine realizeConnectedFields(state, mesh, rc)
      ! TODO: this method may move into the NUOPC_ utility layer
      type(ESMF_State)                :: state
      type(ESMF_Mesh)                 :: mesh
      integer, intent(out), optional  :: rc
      ! local variables
      character(len=80), allocatable  :: fieldNameList(:)
      integer                         :: i, itemCount, k
      type(ESMF_Field)                :: field

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
        ! create a Field
        field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
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
      deallocate(fieldNameList)

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
    type(ESMF_Mesh)                 :: mesh
    type(ESMF_Array)                :: array
    character(len=80), allocatable  :: fieldNameList(:)
    integer                         :: itemCount, count
    integer                         :: i, ii, j, k, index
    type(ESMF_Field)                :: field
    real(ESMF_KIND_R8), pointer     :: dataptr(:)
    real(ESMF_KIND_R8),allocatable  :: coords(:), varbuf(:,:,:)
    integer, save                   :: slice=1
    integer                         :: totalnodes, ownnodes, lbnd(1), ubnd(1)
    real(ESMF_KIND_R8)              :: x, y, z, errsum, var, deg2rad
    real(ESMF_KIND_R8)              :: sendbuf(2), recvbuf(2)
    type(ESMF_VM)                   :: vm
    integer                         :: PetNo, ipedims(3)
    integer                         :: countlon
    type(InternalState)             :: is
    real, parameter                 :: PI=3.1415927
    real, parameter                 :: earthradius = 6371.0
  
    deg2rad = PI/180.0
   
    rc = ESMF_SUCCESS

    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ipedims=is%wrap%ipedims
    countlon = is%wrap%countlon
    
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up local pet info
    call ESMF_VMGet(vm, localPet=PetNo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing DATAIPE from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, stopTime=stopTime, currTime=currTime, &
      timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
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
      
      ! Get the mesh from the first item
      call ESMF_StateGet(importState, field=field, itemName=fieldNameList(1), rc=rc)   
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_FieldGet(field, mesh=mesh, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

#if 0
      call ESMF_FieldGet(field, array=array, rc=rc)
      if (slice > 1) then 
        call ESMF_ArrayWrite(array, 'ipefield.nc', timeslice=slice-1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
#endif

      call ESMF_FieldGet(field, farrayPtr=dataptr, computationalLBound=lbnd, &
          computationalUBound=ubnd, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      totalnodes = size(dataptr,1)
      allocate(coords(totalnodes*3))
      
      call ESMF_MeshGet(mesh, ownedNodeCoords = coords, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
     errsum = 0.0
     count = 0
     do i=1,totalnodes
        ! if the destination point is mapped
        if (dataptr(i)>0) then
          x=coords((i-1)*3+1)*deg2rad
          y=coords((i-1)*3+2)*deg2rad
          z=(coords((i-1)*3+3)-1)*earthradius
          var = 2.0+cos(y)*cos(y)*cos(2*x)
          errsum = errsum + (dataptr(i)-var)/var
          count = count+1
        endif
     enddo
     if (slice > 1) then
       ! Sum up  all the errors and counts and calculate an average relative error
       sendbuf(1) = errsum
       sendbuf(2) = count
       call ESMF_VMReduce(vm, sendbuf, recvbuf, 2, ESMF_REDUCE_SUM, 0, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
       if (PetNo==0) then
      	 print *, 'Average relative error: ', slice, recvbuf(1)/recvbuf(2)
       endif
     endif

    ! Copy the array to a IPE 3D Array with dimensions (level, lat, lon)
    ! IPE dimension order:  nlevels, lats and lons (1115, 170, 80)
    allocate(varbuf(ipedims(1), ipedims(2), countlon))    
    ! Set to null value
    varbuf=-999.0
    index=1
    do i=1,is%wrap%highestlevel
       do j=1,is%wrap%totallats(i)
          do k=1,countlon
             varbuf(i,j,k)=dataptr(index)
	     index = index+1   
       enddo
       enddo
    enddo
    do ii=1,is%wrap%highestlevel
      do j=1,is%wrap%totallats(ii)
        !adjust the height index based on the maxlevel of the latitude
        i=is%wrap%maxlevels(j)-ii+1
        do k=1, countlon
          varbuf(i,j,k)=dataptr(index)
          index=index+1
        enddo
      enddo
    enddo
    if (index-1 /= is%wrap%localnodes) then
        print *, "Mismatching localnodes ", index-1, is%wrap%localnodes
        call ESMF_LogSetError(rcToCheck=ESMF_FAILURE, &
	      msg="Total number of localnodes do not match", &
              line=__LINE__, file=__FILE__, &
	      rcToReturn=rc)
        return
    endif
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! advance the time slice counter
    slice = slice + 1

  end subroutine

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

  !-----------------------------------------------------------------------------

end module
