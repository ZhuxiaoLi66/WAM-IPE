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
module ipeCap

  !-----------------------------------------------------------------------------
  ! IPE Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS    => SetServices, &
    model_label_DataInitialize  => label_DataInitialize, &
    model_label_Advance => label_Advance, &
    model_label_Finalize => label_Finalize
  use module_finalize_IPE,only   :  finalize_IPE      
  use module_sub_myIPE_Init, only:  myIPE_Init
  use module_update_IPE, only    :  update_IPE
  USE module_IPE_dimension,only  :  NMP,NLP,NPTS2D
  use module_FIELD_LINE_GRID_MKS, only: JMIN_ING, JMAX_ISG, plasma_grid_3d, &
           plasma_grid_Z, IGCOLAT, IGLON, MaxFluxTube, WamField, &
	   JMIN_IN,JMAX_IS

#ifdef ESMF_NETCDF
  use netcdf
#endif

  implicit none
  
  private
  integer, parameter :: MAXNAMELEN = 128

  ! private internal state to keep instance data
  type InternalStateStruct
    integer, pointer :: maxlevels(:)
    integer, pointer :: totalsouthlats(:), totalnorthlats(:)
    real(ESMF_KIND_R8), pointer :: ipelon(:,:,:), ipelat(:,:,:)
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
    use module_input_parameters,only:mype,swEsmfTime
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
      phaseLabelList=(/"IPDv02p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, specRoutine=InitializeData, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
!nm20161003 esmf timing library
    if(swEsmfTime) CALL ESMF_VMWtime(beg_time)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if(swEsmfTime) then
      CALL ESMF_VMWtime(end_time)    
      write(unit=9999,FMT=*)mype,"ModelAdvance endT=",(end_time-beg_time)
    end if !(swEsmfTime) then
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    ! attach specializing method(s)
    !---finalize IPE
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
     specRoutine=finalize_IPE, rc=rc)
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
      acceptStringList=(/"IPDv02p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! Enabeling the following macro, i.e. renaming it to WITHIMPORTFIELDS,
    ! will result in a model component that advertise import Field dependencies.
    ! In the single model case, where there isn't another model to satisfy these
    ! dependencies, it is expected to be caught by the compatability checking.

    ! import fields from WAM
    call NUOPC_Advertise(importState, StandardNames=(/ &
      "northward_wind_neutral         ", &
      "eastward_wind_neutral          ", &
      "upward_wind_neutral            ", &
      "temp_neutral                   ", &
      "O_Density                      ", &
      "O2_Density                     ", &
      "N2_Density                     " &
      /), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#if 0
    ! set Component name so it becomes identifiable in the output
    call ESMF_GridCompSet(gcomp, name="IPE", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif
      
#if 0    
    ! exportable field: air_pressure_at_sea_level
    call NUOPC_Advertise(exportState, &
      StandardName="air_pressure_at_sea_level", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! exportable field: surface_net_downward_shortwave_flux
    call NUOPC_Advertise(exportState, &
      StandardName="surface_net_downward_shortwave_flux", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field)        :: field
    type(ESMF_Grid)         :: gridIn
    type(ESMF_Grid)         :: gridOut
    type(ESMF_Mesh)         :: mesh

    rc = ESMF_SUCCESS
    
    ! initialize IPE 
    call myIPE_Init (gComp, importState, exportState, clock, rc)

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
      integer :: totalhalo, recvfm, sendto
      real(ESMF_KIND_R8), pointer    :: ipelon(:,:,:), ipelat(:,:,:), ipehgt(:,:)
      real(ESMF_KIND_R8), pointer    :: sendbuf1(:), recvbuf1(:)
      integer(ESMF_KIND_I4), pointer :: totalheight(:)
      integer(ESMF_KIND_I4), pointer :: elementIds(:), elementTypes(:), elementConn(:)
      integer(ESMF_KIND_I4), pointer :: nodeIds(:), nodeOwners(:)
      real(ESMF_KIND_R8), pointer    :: nodeCoords(:)
      integer(ESMF_KIND_I4), pointer :: southind(:), northind(:)
      integer(ESMF_KIND_I4), pointer :: totalnorthlats(:), totalsouthlats(:)
      integer(ESMF_KIND_I4), pointer :: totallats(:)
      integer(ESMF_KIND_I4), pointer :: maxlevs(:)
      integer(ESMF_KIND_I4), pointer :: elementCnt(:), nodeCnt(:), sendbuf(:), recvbuf(:)
      real(ESMF_KIND_R8), pointer :: conntbl(:), globalCoords(:,:), fptr2d(:,:), fptr1d(:)
      integer, pointer :: allCounts(:), connectbase(:)
      integer :: count, remind, start, neightbor, base, increment
      integer :: totalnodes, localnodes, totalelements, globalTotal, globalTotalelmt
      real(ESMF_KIND_R8) :: maxheight, lon, lat, hgt
      integer :: i, j, k, ii, jj, count1, count3, count8, localcount
      integer :: startid, nodestartid
      integer :: neighbor, elmtcount, countup, diff, lastlat, base1, save
      logical :: even
      real(ESMF_KIND_R8) :: rad2deg
      type(InternalState)  :: is
 
      rc = ESMF_SUCCESS
    
      rad2deg = 180.0/3.14159265
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

      maxheight = 800.0

#ifdef ESMF_NETCDF
#if 0
      ipefilename = 'data/ipe3dgrid.nc'

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
#endif

      ! MaxFluxTube: 1115, NLP: 170, NMP: 80
      ipedims(1)=MaxFluxTube
      ipedims(2)=NLP
      ipedims(3)=NMP

      ! Step 1: calculate number of longitude lines in the local processor


      if (ipedims(3) >= PetCnt) then
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
      ! print *, PetNo, ' start from ', start, ' with ', count, ' columns'

      ! IPE dimension order:  nlevels, lats and lons (1115, 170, 80)
      if (PetCnt > 1) then
          allocate(ipelon(ipedims(1), ipedims(2), count+1), &
     	      ipelat(ipedims(1), ipedims(2), count+1))
      else
          allocate(ipelon(ipedims(1), ipedims(2), count), &
     	      ipelat(ipedims(1), ipedims(2), count))
      endif
      allocate(ipehgt(ipedims(1), ipedims(2)))
      allocate(maxlevs(ipedims(2)))

      maxlevs = JMAX_ISG - JMIN_ING + 1
      
      do i=1,count
         do j=1,ipedims(2)
            do k=1,maxlevs(j)
               ipelon(k,j,i)=rad2deg*plasma_grid_3d(k,j,i+start-1,IGLON)
               ipelat(k,j,i)=90.0-rad2deg*plasma_grid_3d(k,j,i+start-1,IGCOLAT)
            enddo
         enddo
      enddo
      ipehgt = plasma_grid_Z/1000.0

      if (PetCnt > 1) then
        ! Get the coordinates for the halo column
        totalhalo = ipedims(1)*ipedims(2)*2
        allocate(sendbuf1(totalhalo), recvbuf1(totalhalo))
        ! copy the first column to sendbuf, first lon, then lat
        i=1
        do j=1,ipedims(2)
           do k=1,ipedims(1)
             sendbuf1(i)=ipelon(k,j,1)
             sendbuf1(i+1)=ipelat(k,j,1)
             i=i+2
           enddo
        enddo
        sendto = PetNo-1
        recvfm = PetNo+1
        if (sendto == -1) sendto = PetCnt-1
        if (recvfm == PetCnt) recvfm = 0
        call ESMF_VMSendRecv(vm, sendbuf1, totalhalo, sendto, &
      	     recvbuf1, totalhalo, recvfm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, &
             file=__FILE__)) &
             return  ! bail out
        i=1
        do j=1,ipedims(2)
           do k=1,ipedims(1)
             ipelon(k,j,count+1) = recvbuf1(i)
             ipelat(k,j,count+1) = recvbuf1(i+1)
             i=i+2
           enddo
        enddo
        deallocate(sendbuf1, recvbuf1)
      endif
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

      !Find the total node count and cell count in the local PE
      ! hardcode the maxheight to 800km, find how many height levels are within the range
      ! for each latitute tube
      ! Need to separate it into south hemisphere and north hemisphere 
      allocate(totalheight(ipedims(2)),northind(ipedims(2)), southind(ipedims(2)))
      do j=1,ipedims(2)
        do i=1,(maxlevs(j)-1)/2
          if (ipehgt(i,j)> maxheight) then
             northind(j)=i-1
	     exit
          endif
        enddo
        if (i==(maxlevs(j)/2+1)) then 
	   northind(j)=i-1 
           ! include the highest level in north hemisphere if it is below 800km 
           if (ipehgt(i,j) < maxheight .and. northind(j)<northind(1)) then
	      northind(j)=northind(j)+1
       	      if (PetNo==0) print *, 'Add highest level:', j,i,ipehgt(i,j)
           endif
        endif
        do i=maxlevs(j),(maxlevs(j)+1)/2+1,-1
          if (ipehgt(i,j)> maxheight) then
             southind(j)=i+1
	     exit
          endif
        enddo
        if (i==(maxlevs(j)+1)/2) southind(j)=i+1 
        if (southind(j)==northind(j)) southind(j)=southind(j)+1
        totalheight(j) = northind(j)+maxlevs(j)-southind(j)+1
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

      ! count number of latitudes and local elements at each height for both north and south hemispheres
      ! should allocate totallats for all the heights including south and north hemisphere
      ! and calculate the number of lats separately for north and south hemisphere
      totalelements=0
      allocate(totalnorthlats(northind(1)), totalsouthlats(northind(1)))
      allocate(totallats(northind(1)))
      ! North Hemisphere
      do j=1,northind(1)
         do i=1,ipedims(2)
            if (northind(i)<j) exit
         enddo
         totalnorthlats(j)=i-1
         ! do not count the toppest level
         if (j<northind(1)) then
            totalelements=totalelements + count*(totalnorthlats(j)-1)
         endif
      enddo
      ! South Hemisphere
      do j=1,northind(1)
         do i=1,ipedims(2)
            jj=maxlevs(i)-j+1
            if (southind(i)>jj) exit
         enddo
         totalsouthlats(j)=i-1
         ! do not count the toppest level
         if (j<northind(1)) then
            totalelements=totalelements + count*(totalsouthlats(j)-1)
         endif
      enddo

      print *, "Total nodes/elements ", PetNo, startid, totalnodes, localnodes, totalelements
#if 1
      ! for debugging only
      if (PetNo == 0) then      
      write(90,*) maxheight, totalnodes
      do j=1,ipedims(2)
         write(90,*) j, maxlevs(j), northind(j), maxlevs(j)-southind(j)+1, totalheight(j)
      enddo
      do j=1,northind(1)
         write(90,*) j, totalnorthlats(j), totalsouthlats(j)
      enddo
      endif
#endif

      ! totalelements include one extra row connecting the north and south hemisphere
      ! should not multiple totalelements by 2 assuming north and south has same number of elements
      totalelements = totalelements + count*(northind(1)-1)

      ! Do the halo longitude
      if (PetNo == PetCnt-1) then
         neighbor=0
         base=1
         increment=nodeCnt(2)
      else
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
         do j=1,totalnorthlats(i)
            do k=1,count
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
               lon=ipelon(i,j,count+1)
               lat=ipelat(i,j,count+1)
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
    do j=1,totalsouthlats(ii)
      !adjust the height index based on the maxlevel of the latitude
      i=maxlevs(j)-ii+1
      do k=1,count
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
       nodeIds(count1)=base
       base=base+increment
       lon=ipelon(i,j,count+1)
       lat=ipelat(i,j,count+1)
       hgt=ipehgt(i,j)
       call convert2Cart(lon,lat,hgt,nodeCoords(count3:count3+2))
       nodeOwners(count1)=neighbor
       count1=count1+1
       count3=count3+3
       endif 
     enddo
  enddo
  nodestartid = startid+1
  print *, PetNo, 'totalnodes, localnodes ', count1-1, localcount-1
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
    do j=1,totalnorthlats(i)-1
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
    base=base+totalnorthlats(i)*(count+1)
    connectbase(i)=base-(count+1)
  enddo
  ! add the top level to base
  base = base+totalnorthlats(i)*(count+1)
  connectbase(i)=base-(count+1)
  ! South Hemisphere
  do i=1,northind(1)-1
    do j=1,totalsouthlats(i)-1
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

    !!!!!!!!!!!!!!!  
    !!! TBD: Need to check if using totalsouthlats() to calculate the local node id is correct
    !!!!!!!!!!!!!!!  
    ! Add one extra row of cells to connect the north and south hemisphere
    ! Need to figure out the local node id of the last row in North Hemisphere
    do k=1,count
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*(totalsouthlats(i)-1)+k
         elementConn(count8+1)=base+(count+1)*(totalsouthlats(i)-1)+k+1
         elementConn(count8+2)=connectbase(i)+k+1
         elementConn(count8+3)=connectbase(i)+k
	 count1=count1+1
	 count8=count8+8
    enddo      
    ! If PetCnt==1, connect the last column back to the first column
    if (PetCnt == 1) then
         elementIds(count1)=count1+startid
         elementConn(count8)=base+(count+1)*(totalsouthlats(i)-1)+k
         elementConn(count8+1)=base+(count+1)*(totalsouthlats(i)-1)+1
         elementConn(count8+2)=connectbase(i)+1
         elementConn(count8+3)=connectbase(i)+k
	 count1=count1+1
	 count8=count8+8
    endif

    base=base+totalsouthlats(i)*(count+1)
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   Need to take account of the difference between totalnorthlats and totalsouthlats
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ii=1,2  ! ii=1 is North and 2 is South
    if (ii==1) then 
       totallats = totalnorthlats
    else
       totallats = totalsouthlats
    endif
    do i=1, northind(1)-2
      ! base is the starting local index to the elementConn table for the level i+1
      if (ii==1) then
         base = base+count*(totallats(i)-1)*8
      else
         ! For south hemisphere, we have one extra row to connect the south to north
         !!!!!!!!!!!  Should I use totalsouthlats here?????
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

  is%wrap%maxlevels => maxlevs
  is%wrap%ipelon => ipelon
  is%wrap%ipelat => ipelat
  is%wrap%highestlevel = northind(1)
  is%wrap%totalnorthlats => totalnorthlats
  is%wrap%totalsouthlats => totalsouthlats
  is%wrap%localnodes = localnodes
  is%wrap%startlon = start
  is%wrap%countlon = count
  is%wrap%ipedims = ipedims

  deallocate(southind, northind)
  deallocate(totallats)
  deallocate(nodeIds, nodeOwners, nodeCoords)
  deallocate(elementIds, elementTypes, elementConn)
  !deallocate(ipelon, ipelat, ipehgt)
  deallocate(ipehgt)
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
      real(ESMF_KIND_R8), pointer     :: fptr(:)

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
        if (NUOPC_IsConnected(state, fieldName=fieldNameList(i))) then
          ! create a Field
          field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
            name=fieldNameList(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          fptr = 0.d0
          ! realize the connected Field using the just created Field
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        else
          ! remove a not connected Field from State
          call ESMF_StateRemove(state, (/fieldNameList(i)/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif
      enddo
      deallocate(fieldNameList)

    end subroutine realizeConnectedFields

  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine InitializeData(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_State)       :: importState

    rc = ESMF_SUCCESS
    
    ! query the Component for its importState
    call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
        
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
  use module_input_parameters,only:mype,swEsmfTime
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
    integer                         :: itemCount, count, count1
    integer                         :: i, ii, iii, j, k, n, index
    type(ESMF_Field)                :: field
    real(ESMF_KIND_R8), pointer     :: dataptr(:)
    real(ESMF_KIND_R8),allocatable  :: varbuf(:,:,:)
    integer, save                   :: slice=1
    integer                         :: totalnodes, ownnodes, lbnd(1), ubnd(1)
    real(ESMF_KIND_R8)              :: x, y, z, errsum, var, deg2rad
    real(ESMF_KIND_R8)              :: sendbuf(2), recvbuf(2)
    type(ESMF_VM)                   :: vm
    integer                         :: PetNo, PetCnt, ipedims(3)
    integer                         :: countlon, startlon
    type(InternalState)             :: is
    real, parameter                 :: PI=3.14159265
    real, parameter                 :: earthradius = 6371.0
    
    character(len=256)              :: filename
    integer                         :: nc, varid, status
    integer                         :: start3(3), count3(3)  
    integer                         :: midway

    integer, parameter              :: maxItemCount=7
    real(ESMF_KIND_R8)              :: min_field(maxItemCount)
    real(ESMF_KIND_R8)              :: max_field(maxItemCount)
    real(ESMF_KIND_R8)              :: global_min_field(maxItemCount)
    real(ESMF_KIND_R8)              :: global_max_field(maxItemCount)

    real(ESMF_KIND_R8),parameter    :: BAD_VALUE = -999.0

    logical :: checkFields
    character (len=10) :: value

!---
!nm20161003: esmf timing lib
    real(ESMF_KIND_R8) :: beg_time, end_time
!---
    deg2rad = PI/180.0

    ! Init to success   
    rc = ESMF_SUCCESS

    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    startlon = is%wrap%startlon
    countlon = is%wrap%countlon
    ipedims = is%wrap%ipedims

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up local pet info
    call ESMF_VMGet(vm, localPet=PetNo, PetCount=PetCnt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Look for CheckFields attribute to see if we should be
    ! checking the fields for errors
    call ESMF_AttributeGet(gcomp, &
      name="CheckFields", value=value, defaultvalue="false", &
      convention="NUOPC", purpose="Instance", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! convert value to checkfields
    checkFields=.false.
    if (trim(value)=="true") then
       checkFields=.true.
    endif


    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing IPE from: ", rc=rc)
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

    ! Make sure itemCount is not bigger than maxItemCount
    if (itemCount > maxItemCount) then
       call ESMF_LogSetError(rcToCheck=ESMF_FAILURE, &
     msg="actual number of items is bigger than the  maximum allowed", &
             line=__LINE__, file=__FILE__, &
	     rcToReturn=rc)
        return
    endif

    ! Allocate array to hold field names
    allocate(fieldNameList(itemCount))

    ! arrange the fields in the same order as defined in wamfield
    if (itemCount == 7) then
        fieldNameList(1) = 'temp_neutral'
        fieldNameList(2) = 'eastward_wind_neutral'
        fieldNameList(3) = 'northward_wind_neutral'
        fieldNameList(4) = 'upward_wind_neutral'
        fieldNameList(5) = 'O_Density'
        fieldNameList(6) = 'O2_Density'
        fieldNameList(7) = 'N2_Density'
    endif

    ! Init field to bad value
    WamField = BAD_VALUE

    ! loop through items 
    do n=1,itemCount

      ! Get the mesh from the first item
      call ESMF_StateGet(importState, field=field, itemName=fieldNameList(n), rc=rc)   
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

    ! Get pointer from Field
    call ESMF_FieldGet(field, farrayPtr=dataptr, &
                       computationalLBound=lbnd, &
                       computationalUBound=ubnd, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! Copy the array to a IPE 3D Array with dimensions (level, lat, lon)
    ! IPE dimension order:  nlevels, lats and lons (1115, 170, 80)

    ! Init index into dataptr
    index=1

    ! Northern hemisphere
    do i=1,is%wrap%highestlevel
       do j=1,is%wrap%totalnorthlats(i)
          do k=startlon,startlon+countlon-1
             WamField(i,j,k,n)=dataptr(index)
 	     index = index+1     
          enddo
       enddo
    enddo

    ! Southern hemisphere
    do ii=1,is%wrap%highestlevel
      do j=1,is%wrap%totalsouthlats(ii)
        !adjust the height index based on the maxlevel of the latitude
        i=is%wrap%maxlevels(j)-ii+1
        do k=startlon, startlon+countlon-1
           wamField(i,j,k,n)=dataptr(index)
           index=index+1
        enddo
      enddo
    enddo

    ! Make sure that we're ending at the correct place
    if (index-1 /= is%wrap%localnodes) then
        print *, "Mismatching localnodes ", index-1, is%wrap%localnodes
        call ESMF_LogSetError(rcToCheck=ESMF_FAILURE, &
	      msg="Total number of localnodes do not match", &
              line=__LINE__, file=__FILE__, &
	      rcToReturn=rc)
        return
    endif

    ! Check for problem data
    if (checkFields) then
       ! TODO: if there's a problem point perhaps have an option to stop with an error? 

       ! Skip the first timestep, because the info coming from WAM isn't correct then
       if (slice > 1) then

          ! Init to values bigger and smaller than expected
          min_field(n)=HUGE(wamfield)
          max_field(n)=-HUGE(wamfield)

          ! Loop checking values
       	  do j=1,ipedims(2)
               do ii=JMIN_IN(j), JMAX_IS(j)
	         do k=1,countlon

	       	  ! Check points under 783km, since IPE starts its extrap at 782km
                  ! (The next level starts at 786.7km)
                  if (plasma_grid_Z(ii,j) <= 783.0*1000.0) then  

		        ! Make sure that there aren't any uninitialized points
                        if (wamfield(ii,j,startlon+k-1,n) .eq. BAD_VALUE) then
		        	 write(*,*) "ERROR ipeCap BAD_VALUE found (", wamfield(ii,j,startlon+k-1,n), &
                                 ") at ",plasma_grid_Z(ii,j),ii,j,startlon+k-1,trim(fieldNameList(n))
		        endif

		        ! Make sure that there aren't any approx uninitialized points
                        if (abs(wamfield(ii,j,startlon+k-1,n) - (BAD_VALUE)) < 0.1) then
		        	 write(*,*) "ERROR ipeCap approx BAD_VALUE found (", wamfield(ii,j,startlon+k-1,n), &
                                 ") at ",plasma_grid_Z(ii,j),ii,j,startlon+k-1,trim(fieldNameList(n))
		        endif

		        ! Make sure that there aren't any nans
                        if (.not. (wamfield(ii,j,startlon+k-1,n) .eq. wamfield(ii,j,startlon+k-1,n))) then
		        	 write(*,*) "ERROR ipeCap NaN found (", wamfield(ii,j,startlon+k-1,n), &
                                 ") at ",plasma_grid_Z(ii,j),ii,j,startlon+k-1,trim(fieldNameList(n))
		        endif

		        ! Make sure that there aren't any Infinities
                        if (wamfield(ii,j,startlon+k-1,n) .gt. HUGE(wamfield)) then
		        	 write(*,*) "ERROR ipeCap +Infinity found (", wamfield(ii,j,startlon+k-1,n), &
                                 ") at ",plasma_grid_Z(ii,j),ii,j,startlon+k-1,trim(fieldNameList(n))
		        endif

                        if (wamfield(ii,j,startlon+k-1,n) .lt. -HUGE(wamfield)) then
		        	 write(*,*) "ERROR ipeCap -Infinity found (", wamfield(ii,j,startlon+k-1,n), &
                                 ") at ",plasma_grid_Z(ii,j),ii,j,startlon+k-1,trim(fieldNameList(n))
		        endif


			! Check temps
			if (n .eq. 1) then
                           if (wamfield(ii,j,startlon+k-1,n) <= 0.0) then
		        	 write(*,*) "ERROR ipeCap temp <= 0.0 (", wamfield(ii,j,startlon+k-1,n),") found in ", &
                                 " wamfield"
        		   endif
                        endif

			! Check O density
			if (n .eq. 5) then
                           if (wamfield(ii,j,startlon+k-1,n) <= 0.0) then
		        	 write(*,*) "ERROR ipeCap O density <= 0.0 (", wamfield(ii,j,startlon+k-1,n),") found in ", &
                                 " wamfield"
        		   endif
                        endif

			! Check O2 density
			if (n .eq. 6) then
                           if (wamfield(ii,j,startlon+k-1,n) <= 0.0) then
		        	 write(*,*) "ERROR ipeCap O2 density <= 0.0 (", wamfield(ii,j,startlon+k-1,n),") found in ", &
                                 " wamfield"
        		   endif
                        endif

			! Check N2 density
			if (n .eq. 7) then
                           if (wamfield(ii,j,startlon+k-1,n) <= 0.0) then
		        	 write(*,*) "ERROR ipeCap N2 density <= 0.0 (", wamfield(ii,j,startlon+k-1,n),") found in ", &
                                 " wamfield"
        		   endif
                        endif

                        ! compute min/max
                        if (wamfield(ii,j,startlon+k-1,n) < min_field(n)) then
                            min_field(n)=wamfield(ii,j,startlon+k-1,n)
                        endif

                        if (wamfield(ii,j,startlon+k-1,n) > max_field(n)) then
                            max_field(n)=wamfield(ii,j,startlon+k-1,n)
                        endif
                   endif
               enddo
           enddo
         enddo
       endif

       ! For validation purpose, write the field out at the second slice
       if (slice==2) then
          filename = 'ipe3dgrid2.nc'
          ! write wamfield to ipe3dgrid.nc
          do i=0, PetCnt-1
             if (PetNo == i) then
                status = nf90_open(filename, NF90_WRITE, nc)
                call CheckNCError(status, filename)
                status = nf90_inq_varId(nc, trim(fieldNameList(n)), varid)
                call CheckNCError(status, 'wamfields')
                start3(1)=1
                start3(2)=1
                start3(3)=startlon
                count3(1)=ipedims(1)
                count3(2)=ipedims(2)
                count3(3)=countlon
	        if (n==1) allocate(varbuf(count3(1), count3(2), count3(3)))
                do ii=1,ipedims(1)
	        do j=1,ipedims(2)
	        do k=1,count3(3)
                  varbuf(ii,j,k) = wamfield(ii,j,startlon+k-1,n)
                enddo
                enddo
                enddo
	        status = nf90_put_var(nc, varid, varbuf, & 
	  	  start3, count3)
                call CheckNCError(status, 'wamfields')
                status = nf90_close(nc)
                call CheckNCError(status, filename)
             endif
             call ESMF_VMBarrier(vm)
          enddo  !i=0,PetCnt-1
          if (PetNo == 0) then
       	     print *, 'Write out ', trim(fieldNameList(n))
          endif
       endif !slice==2
     endif  ! if (checkFields)
    enddo ! itemCount

    ! If we're checking fields then output field min/max
    if (checkFields) then
       if (slice > 1) then
          ! Compute global min
          call ESMF_VMReduce(vm, min_field, global_min_field, itemCount, &
                            ESMF_REDUCE_MIN, 0, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out

          ! Compute global max
          call ESMF_VMReduce(vm, max_field, global_max_field, itemCount, &
                      ESMF_REDUCE_MAX, 0, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out

          ! Only write out if Pet 0
          if (PetNo==0) then
             do n=1,itemCount
                write(*,*) slice,trim(fieldNameList(n))," ipeCap min=",global_min_field(n)
                write(*,*) slice,trim(fieldNameList(n))," ipeCap max=",global_max_field(n)
             enddo
          endif
       endif
    endif	

#if 0  
    ! The error check if for the analytical values only
    ! Check IPE fields againt its analytical values
    errsum = 0.0
    count = 0
    do i=1,countlon
      do j=1,ipedims(2)
         do k=1,is%wrap%maxlevels(j)
            x=is%wrap%ipelon(k,j,i)*deg2rad
            y=is%wrap%ipelat(k,j,i)*deg2rad
            var = 1.0+plasma_grid_Z(k,j)*0.000004+cos(y)*cos(y)*cos(2*x)
            ! only count the values < 500km
	    if (wamfield(k,j,startlon+i-1,1) > 0 .and. &
	        plasma_grid_Z(k,j)<500000) then
              errsum = errsum + (wamfield(k,j,startlon+i-1,1)-var)/var
              count = count+1
	    endif  
         enddo
       enddo
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
#endif

!nm20161003 esmf timing library
    if(swEsmfTime) CALL ESMF_VMWtime(beg_time)
    call update_IPE ( clock , rc )  
    if(swEsmfTime) then
      CALL ESMF_VMWtime(end_time)    
      !if(mype==0.or.mype==8)
      write(unit=9999,FMT=*)mype,"update_IPE endT=",(end_time-beg_time)
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

#if 0
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE IPE ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing IPE from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call update_IPE ( clock , rc )  

    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
 
#endif
  !-----------------------------------------------------------------------------

    !
end module
