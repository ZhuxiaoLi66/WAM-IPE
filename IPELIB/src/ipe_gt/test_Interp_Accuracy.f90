!---------------------------------
!purpose: test the GT for IPE
!---------------------------------

PROGRAM  test_Interp_Accuracy

USE module_precision

USE modSizeFluxTube, ONLY : NPTS, NMP, NLP  ! sizes of flux tube grid


USE moduleInterfaceThermo2Iono, ONLY : INTERFACE__thermosphere_to_FIXED_GEO, &
                            INTERFACE__FIXED_GEO_to_IONOSPHERE

USE modSizeFixedGridIono, ONLY : nFixedGridIonoHeights, &  !  183 heights
                                 nFixedGridIonoLats, &
                                 nFixedGridIonoLons

USE modSizeFixedGridThermo, ONLY : nFixedGridThermoHeights, fixedThermoHeights_km

! sizes of thermosphere pressure grid
USE modSizeThermo, ONLY : GT_ht_dim, GT_lat_dim, GT_lon_dim


USE moduleInterfaceIono2Thermo, ONLY : INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO, &
                                       readIPEtoGeoGrid, INTERFACE__FIXED_GRID_to_THERMO

IMPLICIT NONE
      

! check status of fortran open
INTEGER :: FileOpenStatus


!------------------------------------------------------
! File unit number for checking startup IPE values :
!------------------------------------------------------
INTEGER, parameter :: uniterroripetofixed = 500
INTEGER, parameter :: uniterrorfixedtothermo = 501
INTEGER, parameter :: uniterroripetothermo = 502
INTEGER, parameter :: uniterrorthermotofixedlinear = 503
INTEGER, parameter :: uniterrorthermotofixedlog = 504
INTEGER, parameter :: uniterrorfixedtoipelinear = 505
INTEGER, parameter :: uniterrorfixedtoipelog = 506
INTEGER, parameter :: uniterrorthermotoipelinear = 507
INTEGER, parameter :: uniterrorthermotoipelog = 508

!------------------------------------------------------
! File unit number for checking test grid values :
!------------------------------------------------------
INTEGER, parameter :: unittestionofixeddata = 509
INTEGER, parameter :: unittestthermotofixedgrid = 510
INTEGER, parameter :: unittestipedata = 511
INTEGER, parameter :: unitfixedthermotoiono = 512


!------------------------------------------------------
! File unit number for checking thermosphere values :
!------------------------------------------------------
INTEGER, parameter :: unitCheckThermo =7700

!--------------------------------------------------------------------
! File unit number for checking Thermospheric interpolation values :
!--------------------------------------------------------------------
INTEGER, parameter :: unitCheckThermoInterp = 7800
!---------------------------------------------------
! Write out the Thermospheric interpolated values??
!---------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: debugThermoInterpFileName = 'interpOut.dat'
INTEGER, parameter :: unitCheckThermoInterpBefore = 7600  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckThermoInterpAfter = 7500  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckThermoInterpFluxT = 7400  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckPressureHeightThermo = 8800

!-----------------------------------------------------------------
! File unit number for checking Ionospheric interpolation values :
!-----------------------------------------------------------------
! THIS IS NOT BEING USED RIGHT NOW,  NEED TO WRITE TO THIS FILE
INTEGER, parameter :: unitCheckIonoInterpBefore = 9800
INTEGER, parameter :: unitCheckIonoInterpAfter = 9900

INTEGER, parameter :: numIonoVars = 10

CHARACTER(LEN=30) :: ionoVarName(numIonoVars)

!-----------------------------------------------------------------------------
! File unit number for checking Ionospheric values to the pressure grid
! and file unit number for getting the average height for each pressure level
!-----------------------------------------------------------------------------
INTEGER, parameter :: unitCheckPressureInterp = 8900
INTEGER, parameter :: unitCheckPressureHeight = 7810


!----------------------------------------
! File unit number for IPE startup files
!----------------------------------------
INTEGER, parameter :: unitStartUp = 211

!----------------------------------
! Output variables from GT model
!----------------------------------
REAL(kind=8) :: testthermogrid(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: fixedgridtothermo(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: errorfixedgridtothermo(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: ipetothermo(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: erroripetothermo(GT_ht_dim, GT_lat_dim, GT_lon_dim)

! transposed lon/lat gt grid 
REAL(kind=8) :: transthermogrid(GT_ht_dim, GT_lon_dim, GT_lat_dim)
REAL(kind=8) :: logtransthermogrid(GT_ht_dim, GT_lon_dim, GT_lat_dim) ! using exponential height

REAL(kind=8) :: Temperature_K_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: ht_FROM_GT_km(GT_ht_dim, GT_lat_dim, GT_lon_dim)
!REAL(kind=8) :: O_density_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


!-------------------------
! For GT_thermosphere  :
!-------------------------
REAL*8 :: magnetic_latitude_degrees(91,20)
REAL*8 :: magnetic_longitude_degrees(91,20)


INTEGER :: ii, jj, kk, ll, mm


!----------------------------------------------
! Ionospheric Variables in thermospheric grid
!----------------------------------------------
REAL*8 :: Oplus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL*8 :: Te_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL*8 :: junk_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim) = 0.0

!-------------------------------------------------------------------------------
! Ionospheric Variables from IPE model Ne_, NOplus_, O2plus_, Ti_Oplus
!-------------------------------------------------------------------------------
REAL(kind=8) :: Oplus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: junk_density_from_IPE(NPTS, NMP) =  0.0


!-------------------------------------------------------------
! Ionospheric Variables interpolated to fixed grid (mid_lat)
!-------------------------------------------------------------
REAL(kind=8) :: Oplus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: Te_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: junk_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons) = 0.0
REAL(kind=8) :: testionofixeddata(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: ipetofixeddata(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: erroripetofixeddata(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)


!--------------------------------------------------
! For interpolating from Thermosphere to IPE grid
!--------------------------------------------------
real*8 :: therm_model_geo_long_deg(GT_lon_dim), therm_model_geo_long_rad(GT_lon_dim)
real*8 :: therm_model_geo_lat_deg(GT_lat_dim),  therm_model_geo_lat_rad(GT_lat_dim) 
real*8 :: Altitude_m_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)


real*8 :: Tn_K_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: O_density_m3_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: junk_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim) = 10.0


!-------------------------------- 
! Fixed grid neutral atmosphere
!--------------------------------                                              
REAL(kind=8) :: logthermotofixedgrid(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: thermotofixedgrid(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: testthermotofixedgrid(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: logtestthermotofixedgrid(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: errorthermotofixedgrid(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: logerrorthermotofixedgrid(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: junk_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim) = 0.0

INTEGER :: fileNPTS, fileNMP, fileNLP  ! read from IPE grid file, but not used
REAL (KIND=8), DIMENSION(1:NMP+1) :: mlon_rad
REAL (KIND=8), DIMENSION(NPTS) :: ipeGL  ! read from IPE grid file, but not used

INTEGER :: inGIP1d(NLP), isGIP1d(NLP)   ! 1d in, read from the IPE grid file
! in - n means northpole, is - s means southpole
INTEGER :: inGIP(NMP, NLP), isGIP(NMP, NLP)  ! 2d inGIP, isGIP grids to match up w/ old GIP subroutine call

! ipe grid info
REAL(kind=8) :: pz_plasma_1d(NPTS)  ! 1d in, read from the IPE grid file
REAL(kind=8) :: glat_plasma_3d(npts,nmp), glond_plasma_3d(npts,nmp), pz_plasma_3d(npts,nmp)

! glats in co-lat radians
REAL(kind=8) :: glat_plasma_3d_crad(npts,nmp), glond_plasma_3d_crad(npts,nmp)

! ipe test grid variable
REAL(kind=8) :: testipedata(npts,nmp)
REAL(kind=8) :: logtestipedata(npts,nmp)


REAL(kind=8) :: fixedthermotoiono(npts, nmp)
REAL(kind=8) :: errorfixedthermotoiono(npts, nmp)
REAL(kind=8) :: errorlogfixedthermotoiono(npts, nmp)
REAL(kind=8) :: errorthermotoiono(npts, nmp)
REAL(kind=8) :: errorlogthermotoiono(npts, nmp)
REAL(kind=8) :: logfixedthermotoiono(npts, nmp)
REAL(kind=8) :: junk_plasma(npts, nmp) = 0.0


!----------------------------------------------------
! Output from fixed_geo 
! lat, lon, ht points to be used for 
! interpolation from thermo fixed grid to ipe grid
!-----------------------------------------------------
INTEGER :: ilon1_3d_fixed_ht(npts,nmp), ilon2_3d_fixed_ht(npts,nmp)
INTEGER :: ilat1_3d_fixed_ht(npts,nmp), ilat2_3d_fixed_ht(npts,nmp)
INTEGER :: ispecial_3d_fixed_ht(npts,nmp)
INTEGER :: ihl_3d_fixed_ht(npts,nmp), ihu_3d_fixed_ht(npts,nmp)


LOGICAL :: isFirstCallFixedHeight = .TRUE.

!---------------------
! Physical Parameters
!---------------------
REAL(kind=8), PARAMETER  :: PI = 3.141592654
REAL(kind=8), PARAMETER  :: DTR = PI/180.0 ! Degrees to Radians conversion
REAL(kind=8), PARAMETER  :: R0 = 6.370E06  ! in meters

!-------------------------------
! For converting apex grid data
!-------------------------------
!REAL(kind=8) :: gr(NPTS,nmp)
!REAL(kind=8) :: gcol(NPTS,nmp)
!REAL(kind=8) :: glon(NPTS,nmp)

!----------------------------------------------------------------------
! Debug results from INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE ?
!----------------------------------------------------------------------
INTEGER, parameter :: unitFixedGeo = 15
CHARACTER(LEN=*), PARAMETER :: debugFixedGeoFileName = 'CheckFixedGeo.dat'

!------------------------------------
! File with the IPE grid coordinates
!------------------------------------
CHARACTER(LEN=*), PARAMETER :: plasmaGridFileName = 'plasma_grid.new.ascii'
INTEGER, PARAMETER :: unitGridIPE = 16

!------------------------------------------------------
! File for checking IPE grid values :
!------------------------------------------------------
INTEGER, parameter :: unitCheckIPE = 17
CHARACTER(LEN=*), PARAMETER :: debugIPEFileName = 'checkIPEgrid.txt'
!-------------------------------
! Debug the IPE grid values???
!-------------------------------
LOGICAL, parameter :: debugGridIPE = .FALSE.

CHARACTER*10, parameter :: thermospheric_model_name = 'CMAT2'

logical ::  GIP_switches(GT_lon_dim)

!----------------------
! Files & Directories
!----------------------
CHARACTER(len=200) :: GT_input_dataset , GT_output_dataset
character(len=140) :: ionoStartDir
character(len=140) :: staticFileDir
character(len=140) :: debugDir, pressureDir
CHARACTER(len=140) :: giptogeoFileName



REAL(kind=8) :: high_res_long(nFixedGridIonoLons) 
REAL(kind=8) :: high_res_lat(nFixedGridIonoLats) 
REAL(kind=8) :: high_res_height(nFixedGridIonoHeights)

!-------------------------
! For timing the code :
!-------------------------
REAL :: startTime, endTime
REAL :: timeFluxGridtoFixedGrid = 0.
REAL :: timeFixedGridtoPressureGrid = 0.
REAL :: timePressureGridtoFixedGrid = 0.
REAL :: timeFixedGridtoFluxGrid = 0.


! BEGIN CODE ====================================================================================
staticFileDir = '/scratch1/portfolios/NCEPDEV/swpc/noscrub/Leslie.Mayer/DATA/IONO/static_files/'
debugDir = '/scratch1/portfolios/NCEPDEV/swpc/noscrub/Leslie.Mayer/DATA/IONO/gtgipIPE/DEBUG/interpAccur/TEST/'
! has the file w/ the pressure grid heights
pressureDir = '/scratch1/portfolios/NCEPDEV/swpc/noscrub/Leslie.Mayer/DATA/IONO/gtgipIPE/DEBUG/interpAccur/'

debugDir = TRIM(debugDir)

!-----------------------------------------------------------
! open files for writing out interpolation error arrays
!-----------------------------------------------------------
OPEN (uniterroripetofixed, FILE=TRIM(debugDir)//TRIM("erroripetofixed"), STATUS='REPLACE')
OPEN (uniterrorfixedtothermo, FILE=TRIM(debugDir)//TRIM("errorfixedtothermo"), STATUS='REPLACE')
OPEN (uniterroripetothermo, FILE=TRIM(debugDir)//TRIM("erroripetothermo"), STATUS='REPLACE')
OPEN (uniterrorthermotofixedlinear, FILE=TRIM(debugDir)//TRIM("errorthermotofixedlinear"), STATUS='REPLACE')
OPEN (uniterrorthermotofixedlog, FILE=TRIM(debugDir)//TRIM("errorthermotofixedlog"), STATUS='REPLACE')
OPEN (uniterrorfixedtoipelinear, FILE=TRIM(debugDir)//TRIM("errorfixedtoipelinear"), STATUS='REPLACE')
OPEN (uniterrorfixedtoipelog, FILE=TRIM(debugDir)//TRIM("errorfixedtoipelog"), STATUS='REPLACE')
OPEN (uniterrorthermotoipelinear, FILE=TRIM(debugDir)//TRIM("errorthermotoipelinear"), STATUS='REPLACE')
OPEN (uniterrorthermotoipelog, FILE=TRIM(debugDir)//TRIM("errorthermotoipelog"), STATUS='REPLACE')


!-----------------------------------------------------------
! open files for writing out grid arrays
!-----------------------------------------------------------
OPEN (unittestionofixeddata, FILE=TRIM(debugDir)//TRIM("gridtestionofixeddata"), STATUS='REPLACE')
OPEN (unittestthermotofixedgrid, FILE=TRIM(debugDir)//TRIM("gridtestthermotofixedgrid"), STATUS='REPLACE')
OPEN (unittestipedata, FILE=TRIM(debugDir)//TRIM("gridtestipedata"), STATUS='REPLACE')
OPEN (unitfixedthermotoiono, FILE=TRIM(debugDir)//TRIM("gridfixedthermotoiono"), STATUS='REPLACE')

			    
!--------------------------------------
! Define name of the geo to mag file
!--------------------------------------
giptogeoFileName = TRIM(staticFileDir)//'GIP_Fixed_GEO_grid_lowres'


    !------------------------
    ! Initialize arrays :
    !------------------------
    ht_FROM_GT = 0. ! height at each grid point in 3D                      ! OUTPUT



    !------------------
    ! Get the IPE grid 
    !------------------
    OPEN (unitGridIPE, FILE=TRIM(staticFileDir)//TRIM(plasmaGridFileName), STATUS='OLD')
    READ (UNIT=unitGridIPE,FMT="(i10)") fileNMP
    READ (UNIT=unitGridIPE,FMT="(i10)") fileNLP
    READ (UNIT=unitGridIPE,FMT="(i10)") fileNPTS
    READ (UNIT=unitGridIPE,FMT="(20i10)") inGIP1d(1:NLP)
    READ (UNIT=unitGridIPE,FMT="(20i10)") isGIP1d(1:NLP)
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") mlon_rad( 1: NMP+1 ) !rad
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") Pz_plasma_1d(1:NPTS) ! meters
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") ipeGL(1:NPTS) ! radians
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") GLAT_plasma_3d(1:NPTS, 1:NMP)  ! radians
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") GLONd_plasma_3d(1:NPTS, 1:NMP) ! radians
    CLOSE (unitGridIPE)

    print *,'test_Interp_Accuracy : Done reading the IPE grid.'
    
   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------

       OPEN (unitCheckIonoInterpBefore, FILE=TRIM(debugDir)//TRIM('IonoInterpBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter, FILE=TRIM(debugDir)//TRIM('IonoInterpAfter.txt'), STATUS='REPLACE')

   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------
       OPEN (unitCheckPressureInterp, FILE=TRIM(debugDir)//TRIM('InterpPressureGrid.txt'), STATUS='REPLACE')

   ! Open to read the pressure grid height file
       OPEN (unitCheckPressureHeight, FILE=TRIM(pressureDir)//TRIM('PressureHeightGrid.txt'))

   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------
       OPEN (unitCheckThermoInterpBefore, FILE=TRIM(debugDir)//TRIM('ThermoInterpBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter, FILE=TRIM(debugDir)//TRIM('ThermoInterpAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT, FILE=TRIM(debugDir)//TRIM('ThermoInterpFluxTube.txt'), STATUS='REPLACE')

    !------------------------------------------------------------------------------
    ! Need to change 1d inGIP, isGIP, to 2d to match up with the old GIP subroutine 
    !------------------------------------------------------------------------------
    DO ii = 1, NLP
       inGIP(:,ii) = inGIP1d(ii)
       isGIP(:,ii) = isGIP1d(ii)
    ENDDO ! ii


    !------------------------------------------------------------------------------
    ! Convert heights in meters to height in kilometers
    !------------------------------------------------------------------------------
    DO ii = 1, NPTS
       Pz_plasma_3d(ii,:) = Pz_plasma_1d(ii)/1000.
    ENDDO ! ii


    ! keep co-latitude in radians for grid function interpolation accuracy test ***
    GLAT_plasma_3d_crad = GLAT_plasma_3d
    GLONd_plasma_3d_crad = GLONd_plasma_3d

    !--------------------------------------------------------
    ! Convert co-latitude in radians to latitude in degrees
    !--------------------------------------------------------
    GLAT_plasma_3d = 90. - (GLAT_plasma_3d/DTR)

    !--------------------------------------------------------
    ! Convert longitude in radians to degrees
    !--------------------------------------------------------
    GLONd_plasma_3d = GLONd_plasma_3d/DTR

    !------------------------------------
    ! Set up grid function for ipe grid    
    !------------------------------------

     testipedata = 1.0 + pz_plasma_3d*0.01 + &
                   (cos(GLAT_plasma_3d_crad)**2) * cos(2.0*GLONd_plasma_3d_crad)

     ! For testing logarithmic interpolations
     logtestipedata = 1.0 + 10.**(pz_plasma_3d*0.0001) + &
                   (cos(GLAT_plasma_3d_crad)**2) * cos(2.0*GLONd_plasma_3d_crad)

     print *,' '
     print *,'minval(pz_plasma_3d) in km = ',minval(pz_plasma_3d)
     print *,'maxval(pz_plasma_3d) in km = ',maxval(pz_plasma_3d)
     print *,'min, max of GLAT_plasma_3d_crad = ',minval(GLAT_plasma_3d_crad), maxval(GLAT_plasma_3d_crad)
     print *,'min, max of GLONd_plasma_3d_crad = ',minval(GLONd_plasma_3d_crad), maxval(GLONd_plasma_3d_crad)     
     print *,' '

    !--------------------------------------
    ! Define name of the geo to mag file
    !--------------------------------------
    !geoToMagFileName = TRIM(staticFileDir)//'/Geographic_to_Magnetic_91_20.2000.0.format'


    ! read in ht_FROM_GT
    READ (unitCheckPressureHeight,*) ht_from_gt
    ! convert to km
    ht_from_gt_km = ht_from_gt/1000.

    !print *,'test_Interp_Accuracy : ht_from_gt_km(:,1,1) = ', ht_from_gt_km(:,1,1)
    !print *,'test_Interp_Accuracy : ht_from_gt_km(1,:,1) = ', ht_from_gt_km(1,:,1)
    !print *,'test_Interp_Accuracy : ht_from_gt_km(1,1,:) = ', ht_from_gt_km(1,1,:)

     print *,'min, max height of ht_from_gt_km = ',minval(ht_from_gt_km), maxval(ht_from_gt_km)


      !--------------------------------------------------
      ! Initialize all of the GIP switches to be .FALSE.
      !--------------------------------------------------
      GIP_switches(:) = .FALSE.
      GIP_switches(7) = .TRUE.   ! input_Auroral_production_is_single_overall_rate


      !-----------------------------------------------------------
      ! Set up lat, longs, used for interpolation to the IPE grid
      !-----------------------------------------------------------
      do ii = 1 , GT_lon_dim
         therm_model_geo_long_deg(ii) = (ii - 1) * 18.
      enddo

      do ii = 1 , GT_lat_dim
         therm_model_geo_lat_deg(ii) = (ii - 46) * 2.
      enddo

      ! convert to radians
      therm_model_geo_long_rad = therm_model_geo_long_deg*DTR

      ! convert to co-latitude, then to radians
      therm_model_geo_lat_rad = (90. - therm_model_geo_lat_deg)*DTR

      print *,' '
      print *,'min,max of therm_model_geo_lat_rad = ',minval(therm_model_geo_lat_rad), maxval(therm_model_geo_lat_rad)
      print *,'min,max of therm_model_geo_lat_rad = ',minval(therm_model_geo_long_rad), maxval(therm_model_geo_long_rad)
      print *,' '


     !----------------------------------------------
     ! Read in the fixed grid ionospheric grid for
     ! nearest point finding
     ! (only need to do this once)
     !----------------------------------------------
     CALL readIPEtoGeoGrid( giptogeoFileName)



!------------------------------------------------
! Open IPE species files
! One for each species, with all time steps in it
!------------------------------------------------

      ! print *,'test_Interp_Accuracy : inGIP1d(1:NLP) = ', inGIP1d(1:NLP)
      ! print *,'test_Interp_Accuracy : isGIP1d(1:NLP) = ', isGIP1d(1:NLP)


     ! For checking ...........
     !do mm = 1, NMP
!
     !   do ll = 1,6
     !       jj = inGIP1d(ll)  
     !       kk = isGIP1d(ll)  

     !       print *, jj, kk

     !       do ii = kk, kk+3
     !          print *, ii, ll, mm, Oplus_density_from_IPE(ii,mm)
     !       end do !ii
     !   end do !ll

     !end do !mm


! ionosphere to fixed grid :
do ii = 1 , nFixedGridIonoHeights 
   high_res_height(ii) = (5.* (ii-1)) + 90.   ! height in km
enddo
!high_res_height = zkm_fixed_ht

print *,' '
print *,'high_res_height (iono fixed grid) = ',high_res_height


!--------------------------------------
! Set up high resolution longitudes
!--------------------------------------
do ii = 1, nFixedGridIonoLons
   high_res_long(ii) = (float(ii-1)) * 4.
enddo

! convert to radians
high_res_long = DTR*high_res_long

!-----------------------------------
! Set up  high resolution latitudes
!-----------------------------------
do ii = 1, nFixedGridIonoLats
   high_res_lat(ii) = (float(ii-46)) * 2.
enddo

! convert to co-lat, then to radians
high_res_lat = (90. - high_res_lat)*DTR

print *,'iono fixed grid :'
print *,'min, max high_res_lat = ',minval(high_res_lat), maxval(high_res_lat)
print *,'min, max high_res_long = ',minval(high_res_long), maxval(high_res_long)
print *,' '

!arr3 = 1.0 + hgt2*0.01 + (cosd(lat2).^2).*cosd(2*lon2)

! Set up iono fixed grid to the grid function
do ii = 1, nFixedGridIonoHeights
   do jj = 1, nFixedGridIonoLats
      do kk = 1, nFixedGridIonoLons
                         ! height, lat, lon
         testionofixeddata(ii,jj,kk) = 1.0 + high_res_height(ii)*0.01 + &
                                        (cos(high_res_lat(jj))**2) * cos(2.0*high_res_long(kk))
         !print *,'ii,jj, kk, testionofixeddata(ii,jj,kk) = ',ii,jj, kk, testionofixeddata(ii,jj,kk)
         !print *,'high_res_height(ii) = ',high_res_height(ii)
         !print *,'high_res_lat(jj) = ',high_res_lat(jj)
         !print *,'high_res_long(kk) = ',high_res_long(kk)

      enddo ! ii
   enddo ! jj
enddo ! kk


    !print *,'size(testionofixeddata) = ',size(testionofixeddata)
    !print *,'nFixedGridIonoHeights*nFixedGridIonoLats*nFixedGridIonoLons = ',&
    !         nFixedGridIonoHeights*nFixedGridIonoLats*nFixedGridIonoLons

  !---------------------------------------------------------------
  ! Check the cpu time to get timing of the interface subroutine
  !---------------------------------------------------------------
  !CALL CPU_TIME(startTime)

  !-------------------------------------------------------
  ! 1st interpolation :  Interpolate from flux tube grid
  ! to the fixed height grid for all IPE species
  !-------------------------------------------------------
  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   testipedata, &      ! inputs
                                   ipetofixeddata)   ! output

  !CALL CPU_TIME(endTime)
  !------------------------------------------------------------------------
  ! Add up how much time we've been spent in this routine, then avg at end
  !------------------------------------------------------------------------
  !timeFluxGridtoFixedGrid = (endTime-startTime) + timeFluxGridtoFixedGrid
  !WRITE(*,*) "cpu_time for INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO, Oplus   : ",(endTime-startTime)

   !  Calculate the errors in the grid interpolation 
    erroripetofixeddata = abs(testionofixeddata - ipetofixeddata)/testionofixeddata

    print *,'min(erroripetofixeddata) = ', minval(erroripetofixeddata)
    print *,'max(erroripetofixeddata) = ', maxval(erroripetofixeddata)
    print *,'average(erroripetofixeddata) = ', sum(erroripetofixeddata)/size(erroripetofixeddata)
    print *,' '


   ! Write out grid after interpolation  (or only write out the error grid????)
   write(uniterroripetofixed,*) erroripetofixeddata


   ! SET UP GRID FUNCTION HERE
   do ii = 1, GT_ht_dim
      do jj = 1, GT_lat_dim
         do kk = 1, GT_lon_dim

            testthermogrid(ii,jj,kk) = 1.0 + ht_from_GT_km(ii,jj,kk)*0.01 + &
                                       (cos(therm_model_geo_lat_rad(jj))**2) * cos(2.0*therm_model_geo_long_rad(kk))
            !print *,'ht, lat, lon = ',ht_from_GT(ii,jj,kk), therm_model_geo_lat_rad(jj), therm_model_geo_long_rad(kk)

         enddo ! kk
      enddo ! jj
   enddo ! ii


      !------------------------------------------------------
      ! Time the fixed grid to pressure grid interpolation
      !------------------------------------------------------
      !CALL cpu_time(startTime)

      ! This is in the loop b/c of changing heights of the pressure grid
      CALL INTERFACE__FIXED_GRID_to_THERMO ( &
         thermospheric_model_name , GT_ht_dim , GT_lat_dim , GT_lon_dim , &  ! inputs
         junk_high_res_fixed, junk_high_res_fixed, junk_high_res_fixed, &    ! inputs
         junk_high_res_fixed, junk_high_res_fixed, &                     ! inputs
         junk_high_res_fixed, junk_high_res_fixed, &                      ! inputs
         testionofixeddata, junk_high_res_fixed, junk_high_res_fixed, &          ! inputs
         therm_model_geo_long_deg, therm_model_geo_lat_deg, &            ! inputs long, lat in degrees
         Ht_FROM_GT, &                                               ! inputs   ht in meters
         junk_FOR_GT, junk_FOR_GT, junk_FOR_GT, &    ! outputs
         junk_FOR_GT, junk_FOR_GT, &                     ! outputs
         junk_FOR_GT, junk_FOR_GT, &                      ! outputs
         fixedgridtothermo, junk_FOR_GT, junk_FOR_GT)                        ! outputs

      !CALL cpu_time(endTime)

      !------------------------------------------------------------------------
      ! Add up how much time we've been spent in this routine, then avg at end
      !------------------------------------------------------------------------
      !timeFixedGridtoPressureGrid = (endTime-startTime) + timeFixedGridtoPressureGrid
      !WRITE(*,*) "cpu_time for  INTERFACE__FIXED_GRID_to_THERMO  : ",(endTime-startTime)


      ! Calculate the errors in the grid interpolation 

     errorfixedgridtothermo = abs(testthermogrid - fixedgridtothermo)/testthermogrid

    print *,'min(errorfixedgridtothermo) = ', minval(errorfixedgridtothermo)
    print *,'max(errorfixedgridtothermo) = ', maxval(errorfixedgridtothermo)
    print *,'average(errorfixedgridtothermo) = ', sum(errorfixedgridtothermo)/size(errorfixedgridtothermo)

    WRITE(uniterrorfixedtothermo,*) errorfixedgridtothermo


    ! Check ipe -> gt grid error
      ! This is in the loop b/c of changing heights of the pressure grid
      CALL INTERFACE__FIXED_GRID_to_THERMO ( &
         thermospheric_model_name , GT_ht_dim , GT_lat_dim , GT_lon_dim , &  ! inputs
         junk_high_res_fixed, junk_high_res_fixed, junk_high_res_fixed, &    ! inputs
         junk_high_res_fixed, junk_high_res_fixed, &                     ! inputs
         junk_high_res_fixed, junk_high_res_fixed, &                      ! inputs
         ipetofixeddata, junk_high_res_fixed, junk_high_res_fixed, &          ! inputs
         therm_model_geo_long_deg, therm_model_geo_lat_deg, &                ! inputs
         Ht_FROM_GT, &                                               ! inputs   ht, lat, lon
         junk_FOR_GT, junk_FOR_GT, junk_FOR_GT, &    ! outputs
         junk_FOR_GT, junk_FOR_GT, &                     ! outputs
         junk_FOR_GT, junk_FOR_GT, &                      ! outputs
         ipetothermo, junk_FOR_GT, junk_FOR_GT)                        ! outputs


     erroripetothermo = abs(testthermogrid - ipetothermo)/testthermogrid

    print *,' '
    print *,'min(erroripetothermo) = ', minval(erroripetothermo)
    print *,'max(erroripetothermo) = ', maxval(erroripetothermo)
    print *,'average(erroripetothermo) = ', sum(erroripetothermo)/size(erroripetothermo)

    WRITE(uniterroripetothermo,*) erroripetothermo

      !-------------------------------------------------------------------------
      ! For interpolation to IPE grid, switch (ht, lat, lon) to (ht, lon, lat)
      ! for input into the subroutine
      !-------------------------------------------------------------------------
      do ii = 1 , GT_ht_dim
         do jj = 1 , GT_lon_dim
            do kk = 1 , GT_lat_dim

               Altitude_m_FOR_IPE(ii,jj,kk) = ht_FROM_GT(ii,kk,jj)
               Tn_K_FOR_IPE(ii,jj,kk) = Temperature_K_FROM_GT(ii,kk,jj)
               transthermogrid(ii,jj,kk) = testthermogrid(ii,kk,jj)
               logtransthermogrid(ii,jj,kk) = 1.0 + 10.**(ht_from_GT_km(ii,kk,jj)*0.0001) + &
                                       (cos(therm_model_geo_lat_rad(kk))**2) * cos(2.0*therm_model_geo_long_rad(jj))
            enddo
         enddo
      enddo

      ! test grid - what output should look like - in ht, lat, lon
      do ii = 1 , nfixedgridthermoheights
         do jj = 1 , GT_lat_dim
            do kk = 1 , GT_lon_dim
               testthermotofixedgrid(ii,jj,kk) = 1.0 + fixedThermoHeights_km(ii)*0.01 + &
                                              (cos(therm_model_geo_lat_rad(jj))**2) * cos(2.0*therm_model_geo_long_rad(kk))
               !print *,'testthermotofixedgrid(ii,jj,kk) = ',testthermotofixedgrid(ii,jj,kk)
               logtestthermotofixedgrid(ii,jj,kk) = 1.0 + 10.**(fixedThermoHeights_km(ii)*0.0001) + &
                                              (cos(therm_model_geo_lat_rad(jj))**2) * cos(2.0*therm_model_geo_long_rad(kk))

            enddo ! kk
         enddo ! jj
      enddo ! ii
            
print *,' '
print *,'fixedThermoHeights_km = ',fixedThermoHeights_km
print *,' '

      !-------------------------------------------------------------------------
      ! time how long it takes to do pressure grid to fixed grid interpolation
      !-------------------------------------------------------------------------
      !CALL cpu_time(startTime)

      !------------------------------------------------------
      ! interpolate from pressure grid to fixed height grid
      !------------------------------------------------------
        call INTERFACE__thermosphere_to_FIXED_GEO ( &
           GIP_switches, &
           thermospheric_model_name , &
           GT_ht_dim, &                    ! was : therm_model_ht_dim
           GT_lat_dim, &                   ! was : therm_model_lat_dim
           GT_lon_dim, &                   ! was : therm_model_lon_dim
           logtransthermogrid, & !  O_density_FROM_GT,  &             ! was : therm_model_o_density
           junk_FOR_IPE, & ! O2_density_FROM_GT, &             ! was : therm_model_o2_density
           junk_FOR_IPE, & ! N2_density_FROM_GT, &             ! was : therm_model_n2_density
           transthermogrid, & ! Temperature_K_FROM_GT, &        ! was : therm_model_Tn
           junk_FOR_IPE, & ! wind_eastwards_ms1_FROM_GT, &  ! was : therm_model_Vy
           junk_FOR_IPE, &  ! wind_southwards_ms1_FROM_GT, &   ! was : therm_model_Vx
           junk_FOR_IPE, &  ! wind_upwards_ms1_from_gt,  added lrm20121108
           junk_FOR_IPE, & ! was : therm_model_qion3d    ! input
           junk_FOR_IPE, & ! exns, &  ****                       ! was : therm_model_elx
           junk_FOR_IPE, & ! eyns, &  ****                      ! was : therm_model_ely
           therm_model_geo_long_deg, &  ! was : therm_model_geo_long
           therm_model_geo_lat_deg, &   ! was : therm_model_geo_lat 
           Altitude_m_FOR_IPE, &        ! was : therm_model_ht_m, &
           logthermotofixedgrid, &  ! O output -------------------------
           junk_fixed_ht, junk_fixed_ht, &     ! output
           junk_Fixed_ht, junk_Fixed_Ht, junk_Fixed_ht, &
           thermotofixedgrid, &             ! temperature output --------------------------
           junk_fixed_ht, junk_fixed_ht, junk_fixed_ht)                        ! output

           !CALL cpu_time(endTime)
           !timePressureGridtoFixedGrid = (endTime - startTime) + timePressureGridtoFixedGrid 
           !WRITE(*,*) "cpu_time for INTERFACE__thermosphere_to_FIXED_GEO   : ",(endTime-startTime)

           ! Calculate interpolation error

           errorthermotofixedgrid = abs(testthermotofixedgrid - thermotofixedgrid)/testthermotofixedgrid
    print *,' '
    print *,'min(errorthermotofixedgrid) = ', minval(errorthermotofixedgrid)
    print *,'max(errorthermotofixedgrid) = ', maxval(errorthermotofixedgrid)
    print *,'average(errorthermotofixedgrid) = ', sum(errorthermotofixedgrid)/size(errorthermotofixedgrid)

    WRITE(uniterrorthermotofixedlinear,*) errorthermotofixedgrid


           logerrorthermotofixedgrid = abs(logtestthermotofixedgrid - logthermotofixedgrid)/logtestthermotofixedgrid
    print *,' '
    print *,'min(logerrorthermotofixedgrid) = ', minval(logerrorthermotofixedgrid)
    print *,'max(logerrorthermotofixedgrid) = ', maxval(logerrorthermotofixedgrid)
    print *,'average(logerrorthermotofixedgrid) = ', sum(logerrorthermotofixedgrid)/size(logerrorthermotofixedgrid)

    !print *,'maxval(abs(testthermotofixedgrid - thermotofixedgrid)) = ',maxval(abs(testthermotofixedgrid - thermotofixedgrid))

    WRITE(uniterrorthermotofixedlog,*) logerrorthermotofixedgrid


  !-----------------------------------------------------------------
  ! Convert fixed height grid to flux tube grid with perfect inputs
  !-----------------------------------------------------------------

   !    testthermotofixedgrid = 1.0 !  JUST FOR TESTING !!!!!

  !---------------------------------------------------------------------
  ! Time how long it takes to interpolate fixed grid to flux tube grid
  !---------------------------------------------------------------------
  !CALL cpu_time(startTime)

  call INTERFACE__FIXED_GEO_to_IONOSPHERE( &
        therm_model_geo_long_deg, & ! was : geo_grid_longitudes_degrees, &
        therm_model_geo_lat_deg, &  ! was : geo_grid_latitudes_degrees, &
        logtestthermotofixedgrid, &  ! O ------------
        junk_fixed_ht, &  ! O2
        junk_fixed_ht, &  ! N2
        junk_Fixed_Ht, junk_Fixed_Ht, junk_Fixed_Ht, & ! v_east, v_south, v_up
        testthermotofixedgrid, &   ! temperature inputs
        inGIP, isGIP, &  ! was : IN, IS, &
        fixedthermotoiono, &  ! temperature output ---------
        logfixedthermotoiono, &  ! O ------------------
        junk_plasma, &  ! O2 out
        junk_plasma, &  ! N2 out
        GLAt_plasma_3d, &
        GLOnd_plasma_3d, &
        PZ_plasma_3d, &
        junk_plasma, junk_plasma, junk_plasma, &  ! wind east, south, up out 
        ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &  ! output
        ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, &  ! output
        ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht, ihu_3d_fixed_ht, &
        isFirstCallFixedHeight, &
        GIP_switches)

        !CALL cpu_time(endTime)
        !timeFixedGridtoFluxGrid = (endTime - startTime) + timeFixedGridtoFluxGrid 
        !WRITE(*,*) "cpu_time for INTERFACE__FIXED_GEO_to_IONOSPHERE  : ",(endTime-startTime)

        ! Write the grid functions out for plotting
        WRITE (unittestionofixeddata,*) testionofixeddata
        WRITE (unittestthermotofixedgrid,*) testthermotofixedgrid
        WRITE (unittestipedata,*) testipedata
        WRITE (unitfixedthermotoiono,*) fixedthermotoiono

        ! calculate errors and write them out
        errorfixedthermotoiono = abs(testipedata - fixedthermotoiono)/testipedata
    print *,' '
    print *,'min(errorfixedthermotoiono) = ', minval(errorfixedthermotoiono)
    print *,'max(errorfixedthermotoiono) = ', maxval(errorfixedthermotoiono)
    print *,'average(errorfixedthermotoiono) = ', sum(errorfixedthermotoiono)/size(errorfixedthermotoiono)

    WRITE (uniterrorfixedtoipelinear,*) errorfixedthermotoiono

        ! calculate errors and write them out (log data)
        errorlogfixedthermotoiono = abs(logtestipedata - logfixedthermotoiono)/logtestipedata
    print *,' '
    print *,'min(errorlogfixedthermotoiono) = ', minval(errorlogfixedthermotoiono)
    print *,'max(errorlogfixedthermotoiono) = ', maxval(errorlogfixedthermotoiono)
    print *,'average(errorlogfixedthermotoiono) = ', sum(errorlogfixedthermotoiono)/size(errorlogfixedthermotoiono)

    WRITE (uniterrorfixedtoipelog,*) errorlogfixedthermotoiono


! Test thermosphere grid -> ipe grid
  call INTERFACE__FIXED_GEO_to_IONOSPHERE( &
        therm_model_geo_long_deg, & ! was : geo_grid_longitudes_degrees, &
        therm_model_geo_lat_deg, &  ! was : geo_grid_latitudes_degrees, &
        logthermotofixedgrid, &  ! O  --------------------
        junk_fixed_ht, &
        junk_fixed_ht, &
        junk_Fixed_Ht, junk_Fixed_Ht, junk_Fixed_Ht, &
        thermotofixedgrid, &   ! temperature inputs  ----------------
        inGIP, isGIP, &  ! was : IN, IS, &
        fixedthermotoiono, &  ! temperature output --------
        logfixedthermotoiono, &  ! O output -------------------
        junk_plasma, &
        junk_plasma, &
        GLAt_plasma_3d, &
        GLOnd_plasma_3d, &
        PZ_plasma_3d, &
        junk_plasma, junk_plasma, junk_plasma, &  ! wind
        ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &  ! output
        ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, &  ! output
        ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht, ihu_3d_fixed_ht, &
        isFirstCallFixedHeight, &
        GIP_switches)

        ! calculate errors and write them out
        errorthermotoiono = abs(testipedata - fixedthermotoiono)/testipedata
    print *,' '
    print *,'min(errorthermotoiono) = ', minval(errorthermotoiono)
    print *,'max(errorthermotoiono) = ', maxval(errorthermotoiono)
    print *,'average(errorthermotoiono) = ', sum(errorthermotoiono)/size(errorthermotoiono)

    do ii = 1, nmp
       print *,'max(errorthermotoiono(:,ii) = ', maxval(errorthermotoiono(:,ii)), ii
    enddo

    WRITE (uniterrorthermotoipelinear,*) errorthermotoiono

        ! calculate errors and write them out (log data)
        errorlogthermotoiono = abs(logtestipedata - logfixedthermotoiono)/logtestipedata
    print *,' '
    print *,'min(errorlogthermotoiono) = ', minval(errorlogthermotoiono)
    print *,'max(errorlogthermotoiono) = ', maxval(errorlogthermotoiono)
    print *,'average(errorlogthermotoiono) = ', sum(errorlogthermotoiono)/size(errorlogthermotoiono)

    WRITE (uniterrorthermotoipelog,*) errorlogthermotoiono

print *,'test_Interp_Accuracy : END OF test_Interp_Accuracy '

!-----------------------------------------------
! Average the times spent in the subroutines
!-----------------------------------------------
!timeFluxGridtoFixedGrid = timeFluxGridtoFixedGrid/(bigLoop - 1)
!timeFixedGridtoPressureGrid = timeFixedGridtoPressureGrid/(littleLoop - 1)
!timePressureGridtoFixedGrid = timePressureGridtoFixedGrid/(bigLoop - 1)
!timeFixedGridtoFluxGrid = timeFixedGridtoFluxGrid/(bigLoop - 1)
!print *,'INTERPOLATION TIMES : -------------------------------------------------'
!print *,'timeFluxGridtoFixedGrid = ', timeFluxGridtoFixedGrid
!print *,'timeFixedGridtoPressureGrid = ', timeFixedGridtoPressureGrid
!print *,'timePressureGridtoFixedGrid = ', timePressureGridtoFixedGrid
!print *,'timeFixedGridtoFluxGrid = ', timeFixedGridtoFluxGrid
!print *,' '
!print *,'bigLoop, littleLoop = ', bigLoop, littleLoop

END PROGRAM  test_Interp_Accuracy
