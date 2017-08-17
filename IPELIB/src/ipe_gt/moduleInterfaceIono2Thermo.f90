MODULE moduleInterfaceIono2Thermo

USE modSizeFixedGridThermo
USE modSizeFixedGridIono


IMPLICIT NONE

PRIVATE  ! Set everything to private access, except those
         ! marked as public

PUBLIC :: INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO, readIPEtoGeoGrid, &
          INTERFACE__FIXED_GRID_to_THERMO


! lrm The names of nheights, mLats, lLons correspond to :  
! nFixedGridIonoHeights
! nFixedGridIonoLats
! nFixedGridIonoLons


  REAL(kind=8) :: facfac_interface(3,nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
  REAL(kind=8) ::     dd_interface(3,nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)

  INTEGER :: ii1_interface(3,nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
  INTEGER :: ii2_interface(3,nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
  INTEGER :: ii3_interface(3,nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
  INTEGER :: ii4_interface(3,nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)

  INTEGER ::  mlow_interface(nFixedGridIonoLons)
  INTEGER :: mhigh_interface(nFixedGridIonoLons)

  INTEGER, parameter :: lun2 = 101




CONTAINS



SUBROUTINE readIPEtoGeoGrid( giptogeoFileName)

  CHARACTER(*), intent(in) :: giptogeoFileName

!===============================================================================

  print *,'readIPEtoGeoGrid : giptogeoFileName = ', &
           giptogeoFileName


! WAITING FOR THE NEW FILE WHICH GOES ALL THE WAY TO THE POLES  2012/02/24
  open(UNIT=lun2, file=giptogeoFileName,  &
          form='formatted',status='unknown')

  read(lun2,*) facfac_interface
  read(lun2,*) dd_interface
  read(lun2,*) ii1_interface
  read(lun2,*) ii2_interface
  read(lun2,*) ii3_interface
  read(lun2,*) ii4_interface
  read(lun2,*) mlow_interface
  read(lun2,*) mhigh_interface
  close(lun2)


END SUBROUTINE readIPEtoGeoGrid



!---------------------------------------------------------------------------------

SUBROUTINE INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS2D, NMP, NLP, &     ! inputs
                                   nHeights, mLats, lLons, &         ! inputs
                                   ionoVar,  &      ! inputs
                                   ionoVar_high_res_fixed)   ! output

  IMPLICIT NONE

  INTEGER, intent(in) :: NPTS2D 
  INTEGER, intent(in) :: NMP !# of magnetic longitude sectors
  INTEGER, intent(in) :: NLP  !# of field line grids in one longitude sector


  INTEGER, intent(in) :: nHeights ! don't need this?
  integer, intent(in) :: mLats 
  integer, intent(in) :: lLons 


  REAL(kind=8) :: dtotinv, factor, d, fac
  INTEGER :: i , iheight , iup, ido, ih, ip
  INTEGER :: l , m ,  mp , lp , in1, in2
  INTEGER, parameter :: lun2 = 101

  INTEGER :: MHIgh(lLons) , MLOw(lLons)
  INTEGER,DIMENSION(3,nHeights,mLats,lLons) :: ii1, ii2, ii3, ii4

  !INTEGER,DIMENSION(3,nheights,mLats,lLons) :: ii1_interface, ii2_interface,&
  !                                             ii3_interface, ii4_interface

  !REAL(kind=8),Dimension(3,nheights,mLats,lLons) ::  facfac_interface, dd_interface
  !INTEGER,DIMENSION(lLons) :: mlow_interface, mhigh_interface  

  REAL(kind=8),Dimension(npts2d,nmp) :: ionoVar  ! Input

  REAL(kind=8),Dimension(nheights,mLats,lLons) :: ionoVar_high_res_fixed  ! Output

  REAL(kind=8),Dimension(3) :: ionoVar_interpolated

!=======================================================================================

!nm20120607dbg
!write(9998,*) ionoVar

  !print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : giptogeoFileName = ', &
  !         giptogeoFileName

!g ---------------------------------------------------
!g  loop over all heights and the lats/longs......
!g ---------------------------------------------------

! Already has coefficients from the height grid

!---------------------------------------------
!  Initialize our output parameters....
!  took out of the loop 20120508lrm
!---------------------------------------------
ionoVar_high_res_fixed = 0.0


  do iheight = 1, nheights
      do l = 1, lLons
         do m = mlow_interface(l), mhigh_interface(l)  ! latitudes
          !!do m = 1, nFixedGridIonoLats  In the future we want to do all latitudes

             !---------------------------------------------
             !  Initialize dtoinv
             !---------------------------------------------
                dtotinv = 0.0

             !g -------------------------------------------------
             !g  The number of contributing flux-tube points
             !g  is always 3....
             !g -------------------------------------------------
              do i = 1, 3

                 !g --------------------------------------
                 !g  The first interpolation uses facfac
                 !g --------------------------------------
                  factor = facfac_interface(i,iheight,m,l)
                  mp = ii1_interface(i,iheight,m,l)
                  lp = ii2_interface(i,iheight,m,l)
                  in2 = ii3_interface(i,iheight,m,l)
                  in1 = ii4_interface(i,iheight,m,l)


                  ionoVar_interpolated(i) = ((ionoVar(in2,mp)-ionoVar(in1,mp))*factor) &
                                            + ionoVar(in1,mp)

                 !g -----------------------------------------------------
                 !g Now we calculate the parameters at each point.....
                 !g -----------------------------------------------------

                  d = dd_interface(i,iheight,m,l)
                  dtotinv = (1./d) + dtotinv


                  ionoVar_high_res_fixed(iheight,m,l) = ionoVar_interpolated(i)/d &
                                            + ionoVar_high_res_fixed(iheight,m,l)

              ENDDO
          
              ionoVar_high_res_fixed(iheight,m,l) = ionoVar_high_res_fixed(iheight,m,l)/dtotinv

              !if (ionoVar_high_res_fixed(iheight,m,l) == 0 .and. iheight == 1) then
              !    print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : ionoVar_high_res_fixed(iheight,m,l) == 0', iheight, m, l
              !    print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : factor, mp, lp, in2, in1 = ', &
              !             factor, mp, lp, in2, in1
              !    print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : ionoVar(in1,mp) = ',ionoVar(in1,mp)
              !    !print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : ARTIFICIALLY SETTTING TO .1 ****'
              !    !ionoVar_high_res_fixed(iheight,m,l) = 0.1
              !endif

          ENDDO ! lats  m

          !----------------------------------------------------------------------------
          !nm20120517
          ! There are some zero values near the magnetic poles (north & south),
          ! we want to get rid of them by assigning the zero values to the lowest or
          ! highest available values.
          ! We need to solve this interpolation problem lrm20120717
          !----------------------------------------------------------------------------
          do m = 1, mlow_interface(l)-1
            ionoVar_high_res_fixed(iheight,m,l) = ionoVar_high_res_fixed(iheight, mlow_interface(l), l)
          end do

          do m = mhigh_interface(l)+1, nFixedGridIonoLats
             ionoVar_high_res_fixed(iheight,m,l) = ionoVar_high_res_fixed(iheight, mhigh_interface(l), l)
          end do


      ENDDO  ! lons   l
  ENDDO  ! iheight

  !print *,'ionoVar_high_res_fixed(1,8,:) = ', ionoVar_high_res_fixed(1,8,:)
  !print *,'ionoVar_high_res_fixed(2,8,:) = ', ionoVar_high_res_fixed(2,8,:)

  !-----------------------------------------------------
  ! Add fix for when values at iheight = 1 are 0
  ! MUST CHECK THIS *********************************
  !-----------------------------------------------------
  !print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : SETTING ionoVar_high_res_fixed(1,m,l) == 0.0 TO LEVEL ABOVE ***'
  !WHERE (ionoVar_high_res_fixed(1,:,:) == 0.0 ) 
  !       ionoVar_high_res_fixed(1,:,:) = ionoVar_high_res_fixed(2,:,:)
  !END WHERE


  !print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : MINVAL(ionoVar) = ',MINVAL(ionoVar)
  !print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : MINVAL(ionoVar_high_res_fixed) = ', MINVAL(ionoVar_high_res_fixed)

!  print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : mlow_interface = ',mlow_interface
!  print *,'INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO : mhigh_interface = ',mhigh_interface

 ! print *,'AFTER FIXING :  '
 ! print *,'ionoVar_high_res_fixed(1,8,:) = ', ionoVar_high_res_fixed(1,8,:)

!  STOP


  return


END SUBROUTINE INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO


!***********************************************************************************
!***********************************************************************************

SUBROUTINE INTERFACE__FIXED_GRID_to_THERMO ( &
         thermospheric_model_name , GT_ht_dim , lat_dim , lon_dim , &
         ne_high_res_fixed, oplus_high_res_fixed, hplus_high_res_fixed, &
         noplus_high_res_fixed, o2plus_high_res_fixed, &
         n2plus_high_res_fixed, nplus_high_res_fixed, &
         Te_high_res_fixed, Ti1_high_res_fixed, Ti2_high_res_fixed, &
         useIPEHeatingRates, heatingRate_high_res_fixed, &  
         therm_geo_long_input, therm_geo_lat_input, &
         therm_Z, &
         therm_Ne_density, therm_oplus_density, therm_hplus_density, &
         therm_noplus_density, therm_o2plus_density, &
         therm_n2plus_density, therm_nplus_density, &
         therm_Te, therm_Ti1, therm_Ti2, &
         neutralHeatingRates)

IMPLICIT NONE


! From modSizeFixedGridIono : -------------------------
! nFixedGridIonoHeights = 183
! nFixedGridIonoLats = 91
! nFixedGridIonoLons = 90
!----------------------------------
!----------------------------------------------------
! Input Variables
!----------------------------------------------------

character*10, intent(in) :: thermospheric_model_name
integer, intent(in) :: GT_ht_dim , lat_dim , lon_dim

REAL(kind=8), intent(in) :: ne_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: oplus_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: hplus_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: noplus_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: o2plus_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: n2plus_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: nplus_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: Te_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: Ti1_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)
REAL(kind=8), intent(in) :: Ti2_high_res_fixed(nFixedGridIonoHeights,nFixedGridIonoLats,nFixedGridIonoLons)

!------------------------------------------------------
! Neutral heating rates interpolated to fixed grid
!------------------------------------------------------
!INTEGER, intent(in) :: numHrate
LOGICAL, intent(in) :: useIPEheatingRates
REAL(kind=8), intent(in) :: heatingRate_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)


REAL(kind=8), intent(in) :: therm_geo_long_input(lon_dim)
REAL(kind=8), intent(in) :: therm_geo_lat_input(lat_dim)

! using ht_from_gt instead of altitude_m_for_ipe (only difference is ht, lat, lon order)
REAL(kind=8), intent(in) :: therm_Z(GT_ht_dim, lat_dim, lon_dim) ! lrm20120601  


!----------------------------------------------------
! Output Variables
!----------------------------------------------------

! Switched lon, lat to match input variables lrm20120601

REAL(kind=8), intent(out) :: therm_Ne_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_oplus_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_hplus_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_noplus_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_o2plus_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_n2plus_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_nplus_density(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_Te(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_Ti1(GT_ht_dim, lat_dim, lon_dim)
REAL(kind=8), intent(out) :: therm_Ti2(GT_ht_dim, lat_dim, lon_dim)

REAL(kind=8), intent(out) :: neutralHeatingRates(GT_ht_dim, lat_dim, lon_dim)



!-------------------------------------------
! Local Variables
!-------------------------------------------
REAL(kind=8) :: therm_geo_long(lon_dim)
REAL(kind=8) :: therm_Z_km


INTEGER :: ispecial, ii



REAL(kind=8) :: factor_lon_array(lon_dim)
REAL(kind=8) :: factor_lon
INTEGER :: ilon_west_array(lon_dim)
INTEGER :: ilon_east_array(lon_dim)
REAL(kind=8) :: factor_lat_array(lat_dim)
REAL(kind=8) :: factor_lat
INTEGER :: ilat_north_array(lat_dim)
INTEGER :: ilat_south_array(lat_dim)
REAL(kind=8) :: factor_ht

REAL(kind=8) :: ne_above_north_west
REAL(kind=8) :: ne_below_north_west 
REAL(kind=8) :: ne_above_north_east
REAL(kind=8) :: ne_below_north_east
REAL(kind=8) :: ne_above_south_west
REAL(kind=8) :: ne_below_south_west
REAL(kind=8) :: ne_above_south_east
REAL(kind=8) :: ne_below_south_east

REAL(kind=8) :: ne_north_west 
REAL(kind=8) :: ne_north_east
REAL(kind=8) :: ne_south_west 
REAL(kind=8) :: ne_south_east 

REAL(kind=8) :: oplus_above_north_west 
REAL(kind=8) :: oplus_below_north_west 
REAL(kind=8) :: oplus_above_north_east 
REAL(kind=8) :: oplus_below_north_east 
REAL(kind=8) :: oplus_above_south_west 
REAL(kind=8) :: oplus_below_south_west 
REAL(kind=8) :: oplus_above_south_east 
REAL(kind=8) :: oplus_below_south_east 

REAL(kind=8) :: oplus_north_west 
REAL(kind=8) :: oplus_north_east 
REAL(kind=8) :: oplus_south_west 
REAL(kind=8) :: oplus_south_east 

    REAL(kind=8) :: hplus_above_north_west 
    REAL(kind=8) :: hplus_below_north_west 
    REAL(kind=8) :: hplus_above_north_east 
    REAL(kind=8) :: hplus_below_north_east 
    REAL(kind=8) :: hplus_above_south_west 
    REAL(kind=8) :: hplus_below_south_west 
    REAL(kind=8) :: hplus_above_south_east 
    REAL(kind=8) :: hplus_below_south_east 

    REAL(kind=8) :: hplus_north_west 
    REAL(kind=8) :: hplus_north_east 
    REAL(kind=8) :: hplus_south_west 
    REAL(kind=8) :: hplus_south_east 

    REAL(kind=8) :: noplus_above_north_west 
    REAL(kind=8) :: noplus_below_north_west 
    REAL(kind=8) :: noplus_above_north_east 
    REAL(kind=8) :: noplus_below_north_east 
    REAL(kind=8) :: noplus_above_south_west 
    REAL(kind=8) :: noplus_below_south_west 
    REAL(kind=8) :: noplus_above_south_east 
    REAL(kind=8) :: noplus_below_south_east 

    REAL(kind=8) :: noplus_north_west 
    REAL(kind=8) :: noplus_north_east 
    REAL(kind=8) :: noplus_south_west 
    REAL(kind=8) :: noplus_south_east 

    REAL(kind=8) :: n2plus_above_north_west 
    REAL(kind=8) :: n2plus_below_north_west 
    REAL(kind=8) :: n2plus_above_north_east 
    REAL(kind=8) :: n2plus_below_north_east 
    REAL(kind=8) :: n2plus_above_south_west 
    REAL(kind=8) :: n2plus_below_south_west 
    REAL(kind=8) :: n2plus_above_south_east 
    REAL(kind=8) :: n2plus_below_south_east 

    REAL(kind=8) :: n2plus_north_west 
    REAL(kind=8) :: n2plus_north_east 
    REAL(kind=8) :: n2plus_south_west 
    REAL(kind=8) :: n2plus_south_east 

    REAL(kind=8) :: o2plus_above_north_west 
    REAL(kind=8) :: o2plus_below_north_west 
    REAL(kind=8) :: o2plus_above_north_east 
    REAL(kind=8) :: o2plus_below_north_east 
    REAL(kind=8) :: o2plus_above_south_west 
    REAL(kind=8) :: o2plus_below_south_west 
    REAL(kind=8) :: o2plus_above_south_east 
    REAL(kind=8) :: o2plus_below_south_east 

    REAL(kind=8) :: o2plus_north_west 
    REAL(kind=8) :: o2plus_north_east 
    REAL(kind=8) :: o2plus_south_west 
    REAL(kind=8) :: o2plus_south_east 

    REAL(kind=8) :: nplus_above_north_west 
    REAL(kind=8) :: nplus_below_north_west 
    REAL(kind=8) :: nplus_above_north_east 
    REAL(kind=8) :: nplus_below_north_east 
    REAL(kind=8) :: nplus_above_south_west 
    REAL(kind=8) :: nplus_below_south_west 
    REAL(kind=8) :: nplus_above_south_east 
    REAL(kind=8) :: nplus_below_south_east 

    REAL(kind=8) :: nplus_north_west 
    REAL(kind=8) :: nplus_north_east 
    REAL(kind=8) :: nplus_south_west 
    REAL(kind=8) :: nplus_south_east 

    REAL(kind=8) :: Te_above_north_west 
    REAL(kind=8) :: Te_below_north_west 
    REAL(kind=8) :: Te_above_north_east 
    REAL(kind=8) :: Te_below_north_east 
    REAL(kind=8) :: Te_above_south_west 
    REAL(kind=8) :: Te_below_south_west 
    REAL(kind=8) :: Te_above_south_east 
    REAL(kind=8) :: Te_below_south_east 

    REAL(kind=8) :: Te_north_west 
    REAL(kind=8) :: Te_north_east 
    REAL(kind=8) :: Te_south_west 
    REAL(kind=8) :: Te_south_east 

    REAL(kind=8) :: Ti1_above_north_west 
    REAL(kind=8) :: Ti1_below_north_west 
    REAL(kind=8) :: Ti1_above_north_east 
    REAL(kind=8) :: Ti1_below_north_east 
    REAL(kind=8) :: Ti1_above_south_west 
    REAL(kind=8) :: Ti1_below_south_west 
    REAL(kind=8) :: Ti1_above_south_east 
    REAL(kind=8) :: Ti1_below_south_east 

    REAL(kind=8) :: Ti1_north_west 
    REAL(kind=8) :: Ti1_north_east 
    REAL(kind=8) :: Ti1_south_west 
    REAL(kind=8) :: Ti1_south_east 

    REAL(kind=8) :: Ti2_above_north_west 
    REAL(kind=8) :: Ti2_below_north_west 
    REAL(kind=8) :: Ti2_above_north_east 
    REAL(kind=8) :: Ti2_below_north_east 
    REAL(kind=8) :: Ti2_above_south_west 
    REAL(kind=8) :: Ti2_below_south_west 
    REAL(kind=8) :: Ti2_above_south_east 
    REAL(kind=8) :: Ti2_below_south_east 

    REAL(kind=8) :: Ti2_north_west 
    REAL(kind=8) :: Ti2_north_east 
    REAL(kind=8) :: Ti2_south_west 
    REAL(kind=8) :: Ti2_south_east 


    REAL(kind=8) :: heatingRate_above_north_west
    REAL(kind=8) :: heatingRate_below_north_west
    REAL(kind=8) :: heatingRate_above_north_east
    REAL(kind=8) :: heatingRate_below_north_east
    REAL(kind=8) :: heatingRate_above_south_west
    REAL(kind=8) :: heatingRate_below_south_west 
    REAL(kind=8) :: heatingRate_above_south_east
    REAL(kind=8) :: heatingRate_below_south_east

    REAL(kind=8) :: heatingRate_north_west
    REAL(kind=8) :: heatingRate_north_east
    REAL(kind=8) :: heatingRate_south_west
    REAL(kind=8) :: heatingRate_south_east 



    REAL(kind=8) :: ne_east 
    REAL(kind=8) :: ne_west 
    REAL(kind=8) :: oplus_east 
    REAL(kind=8) :: oplus_west 
    REAL(kind=8) :: hplus_east
    REAL(kind=8) :: hplus_west 
    REAL(kind=8) :: noplus_east 
    REAL(kind=8) :: noplus_west 
    REAL(kind=8) :: o2plus_east 
    REAL(kind=8) :: o2plus_west 
    REAL(kind=8) :: n2plus_east 
    REAL(kind=8) :: n2plus_west 
    REAL(kind=8) :: nplus_east 
    REAL(kind=8) :: nplus_west 

    REAL(kind=8) :: Te_east 
    REAL(kind=8) :: Te_west 

    REAL(kind=8) :: Ti1_east 
    REAL(kind=8) :: Ti1_west 
    REAL(kind=8) :: Ti2_east 
    REAL(kind=8) :: Ti2_west 

    REAL(kind=8) :: heatingRate_east 
    REAL(kind=8) :: heatingRate_west 

    REAL(kind=8) :: high_res_long(nFixedGridIonoLons) 
    REAL(kind=8) :: high_res_lat(nFixedGridIonoLats) 
    REAL(kind=8) :: high_res_height(nFixedGridIonoHeights)

    INTEGER :: ilon_int
    INTEGER :: ilat_int
    INTEGER :: iht_int
    INTEGER :: ilon
    INTEGER :: ilat
    INTEGER :: iht

    INTEGER :: iht_above
    INTEGER :: iht_below
    INTEGER :: ilat_north
    INTEGER :: ilat_south
    INTEGER :: ilon_west
    INTEGER :: ilon_east


!------------------------------------
! High resolution fixed height grid
!------------------------------------
REAL(kind=8) :: zkm_fixed_ht(nFixedGridIonoHeights)


! NOT USING THIS ***************
!REAL(kind=8) :: fixed_heights_km(interface_hts)  ! was in module definition 20120430lrm
!DATA fixed_heights_km/90.,95.,100.,105.,110.,115.,120.,125., & ! was in module definition 20120430lrm
!       150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400., &
!       450.,500.,550.,600.,700.,800.,900.,1000., &
!       2000.,4000.,6370.,9000./

!===================================================================================
! Begin Code
!===================================================================================

!---------------------------------------------------------------
! Use high resolution fixed ionospheric grid and interpolate to 
! thermospheric grid - 20120508lrm
!---------------------------------------------------------------
!do n = 1 , interface_hts
do ii = 1 , nFixedGridIonoHeights
   zkm_fixed_ht(ii) = (5.* (ii-1)) + 90.  ! Naomi said use this 20120501lrm
   !zkm_fixed_ht(n) = fixed_heights_km(n)
enddo


therm_geo_long(:) = therm_geo_long_input(:)


!----------------------------------------
! Check for longitudes that are negative
!----------------------------------------
do ilon = 1, lon_dim
   if ( therm_geo_long(ilon) < 0.0 ) then
        therm_geo_long(ilon) = therm_geo_long(ilon) + 360.
   endif
enddo

!--------------------------------------
! Set up high resolution longitudes
!--------------------------------------
do ilon_int = 1, nFixedGridIonoLons
   high_res_long(ilon_int) = (float(ilon_int-1)) * 4.
enddo

!-----------------------------------
! Set up  high resolution latitudes
!-----------------------------------
do ilat_int = 1, nFixedGridIonoLats
   high_res_lat(ilat_int) = (float(ilat_int-46)) * 2.
enddo

!----------------------------------------------------------
! - from GIP_ionosphere_plasmasphere.fromresthermo.f90 
!----------------------------------------------------------
! These arrays are the same size (nFixedGridIonoHeights)
high_res_height = zkm_fixed_ht



!-------------------------------------------
! Loop over thermosphere longitudes....
!  ** (or ionospheric ?? they are the same)
!-------------------------------------------
do ilon = 1 , lon_dim 

   !----------------------------------------------------------
   ! longitude interpolation
   ! loop over high_res longs to find points east and west....
   !----------------------------------------------------------

   ispecial = 0

   do ilon_int = 1, nFixedGridIonoLons

      if (high_res_long(ilon_int) > therm_geo_long(ilon)) then
          ilon_east = ilon_int
          ilon_west = ilon_int - 1
          goto 1500  ! we found what we're looking for exit out of loop 
      endif

   enddo ! ilon_int

   ilon_east = 1
   ilon_west = nFixedGridIonoLons
   ispecial = 1
       
1500 continue

   if (ilon_east == 1) then
       ilon_east = 1
       ilon_west = nFixedGridIonoLons
       ispecial = 1
   endif

   if ( ispecial == 0 ) then
        factor_lon = (therm_geo_long(ilon) - high_res_long(ilon_west)) /  &
                     (high_res_long(ilon_east) - high_res_long(ilon_west))
   else

        factor_lon = (therm_geo_long(ilon) - high_res_long(ilon_west)) /  &
                     (high_res_long(ilon_east)+360. - high_res_long(ilon_west))
   endif

   if ( factor_lon > 1.0 .or. factor_lon < 0.0 ) then
        write(6,*) 'INTERFACE__FIXED_GRID_to_THERMO factor lon ', &
                   factor_lon, therm_geo_long(ilon), high_res_long(ilon_west), &
                   high_res_long(ilon_east)
   endif

   factor_lon_array(ilon) = factor_lon
   ilon_west_array(ilon) = ilon_west
   ilon_east_array(ilon) = ilon_east

   !--------------------------------------
   ! Loop over thermoshperic latitudes....
   !--------------------------------------
   do ilat = 1 , lat_dim

      !--------------------------------
      ! latitude interpolation....
      !--------------------------------
      do ilat_int = 1 , nFixedGridIonoLats

         if (high_res_lat(ilat_int) > therm_geo_lat_input(ilat)) then 
             ilat_north = ilat_int
             ilat_south = ilat_int - 1
             goto 2000  ! we found what we're looking for, exit out of loop
         endif
      enddo

2000 continue

      if (ilat_north == 1) then
          ilat_north = 2
          ilat_south = 1
      endif
      factor_lat = (therm_geo_lat_input(ilat) - high_res_lat(ilat_south)) /  &
                   (high_res_lat(ilat_north) - high_res_lat(ilat_south))

      factor_lat_array(ilat) = factor_lat
      ilat_north_array(ilat) = ilat_north
      ilat_south_array(ilat) = ilat_south


      !-------------------------------------
      ! Loop over thermosphere heights....
      !-------------------------------------
      do iht = 1 , GT_ht_dim           

         therm_Z_km = therm_Z(iht, ilat, ilon) / 1000.

         !-----------------------------------------------------------------------------
         ! height interpolation....
         ! Find the index of the two vertical ionospheric grid cells we are in between
         !-----------------------------------------------------------------------------
         do iht_int = 1, nFixedGridIonoHeights 

            if (high_res_height(iht_int) > therm_Z_km) then
               iht_above = iht_int
               iht_below = iht_int - 1
               goto 2500    ! we found what we're looking for, exit out of loop
            endif
         enddo

         ! should be an else statement instead ***********

         !print *,' HERE AT iht_above = nFixedGridIonoHeights .......'
         !print *,'iht_int, iht_above, iht_below = ',iht_int, iht_above, iht_below
         !print *,'therm_Z_km, high_res_height(iht_int) = ',therm_Z_km, high_res_height(iht_int)
         !STOP

         iht_above = nFixedGridIonoHeights         ! fixed lrm20130102
         iht_below = nFixedGridIonoHeights - 1
  
2500 continue


         if (iht_above == 1) then
             iht_above = 2
             iht_below = 1
         endif

         factor_ht = (therm_Z_km - high_res_height(iht_below)) /  &
                     (high_res_height(iht_above) - high_res_height(iht_below))

         !---------------------------------------------------------
         ! cg - make sure factor_ht doesn't get smaller than 0.0 
         !      which will happen for thermospheric
         ! cg - heights below 90km.....
         !---------------------------------------------------------

         if (factor_ht .lt. 0.0) factor_ht = 0.0


         ne_above_north_west = ne_high_res_fixed(iht_above,ilat_north,ilon_west)
         ne_below_north_west = ne_high_res_fixed(iht_below,ilat_north,ilon_west)
         ne_above_north_east = ne_high_res_fixed(iht_above,ilat_north,ilon_east)
         ne_below_north_east = ne_high_res_fixed(iht_below,ilat_north,ilon_east)
         ne_above_south_west = ne_high_res_fixed(iht_above,ilat_south,ilon_west)
         ne_below_south_west = ne_high_res_fixed(iht_below,ilat_south,ilon_west)
         ne_above_south_east = ne_high_res_fixed(iht_above,ilat_south,ilon_east)
         ne_below_south_east = ne_high_res_fixed(iht_below,ilat_south,ilon_east)

         ne_north_west = ((ne_above_north_west - ne_below_north_west)*factor_ht) + ne_below_north_west
         ne_north_east = ((ne_above_north_east - ne_below_north_east)*factor_ht) + ne_below_north_east
         ne_south_west = ((ne_above_south_west - ne_below_south_west)*factor_ht) + ne_below_south_west
         ne_south_east = ((ne_above_south_east - ne_below_south_east)*factor_ht) + ne_below_south_east

         oplus_above_north_west = oplus_high_res_fixed(iht_above,ilat_north,ilon_west)
         oplus_below_north_west = oplus_high_res_fixed(iht_below,ilat_north,ilon_west)
         oplus_above_north_east = oplus_high_res_fixed(iht_above,ilat_north,ilon_east)
         oplus_below_north_east = oplus_high_res_fixed(iht_below,ilat_north,ilon_east)
         oplus_above_south_west = oplus_high_res_fixed(iht_above,ilat_south,ilon_west)
         oplus_below_south_west = oplus_high_res_fixed(iht_below,ilat_south,ilon_west)
         oplus_above_south_east = oplus_high_res_fixed(iht_above,ilat_south,ilon_east)
         oplus_below_south_east = oplus_high_res_fixed(iht_below,ilat_south,ilon_east)

         oplus_north_west = ((oplus_above_north_west - oplus_below_north_west)*factor_ht) + oplus_below_north_west
         oplus_north_east = ((oplus_above_north_east - oplus_below_north_east)*factor_ht) + oplus_below_north_east
         oplus_south_west = ((oplus_above_south_west - oplus_below_south_west)*factor_ht) + oplus_below_south_west
         oplus_south_east = ((oplus_above_south_east - oplus_below_south_east)*factor_ht) + oplus_below_south_east

         hplus_above_north_west = hplus_high_res_fixed(iht_above,ilat_north,ilon_west)
         hplus_below_north_west = hplus_high_res_fixed(iht_below,ilat_north,ilon_west)
         hplus_above_north_east = hplus_high_res_fixed(iht_above,ilat_north,ilon_east)
         hplus_below_north_east = hplus_high_res_fixed(iht_below,ilat_north,ilon_east)
         hplus_above_south_west = hplus_high_res_fixed(iht_above,ilat_south,ilon_west)
         hplus_below_south_west = hplus_high_res_fixed(iht_below,ilat_south,ilon_west)
         hplus_above_south_east = hplus_high_res_fixed(iht_above,ilat_south,ilon_east)
         hplus_below_south_east = hplus_high_res_fixed(iht_below,ilat_south,ilon_east)

         hplus_north_west = ((hplus_above_north_west - hplus_below_north_west)*factor_ht) + hplus_below_north_west
         hplus_north_east = ((hplus_above_north_east - hplus_below_north_east)*factor_ht) + hplus_below_north_east
         hplus_south_west = ((hplus_above_south_west - hplus_below_south_west)*factor_ht) + hplus_below_south_west
         hplus_south_east = ((hplus_above_south_east - hplus_below_south_east)*factor_ht) + hplus_below_south_east

         noplus_above_north_west = noplus_high_res_fixed(iht_above,ilat_north,ilon_west)
         noplus_below_north_west = noplus_high_res_fixed(iht_below,ilat_north,ilon_west)
         noplus_above_north_east = noplus_high_res_fixed(iht_above,ilat_north,ilon_east)
         noplus_below_north_east = noplus_high_res_fixed(iht_below,ilat_north,ilon_east)
         noplus_above_south_west = noplus_high_res_fixed(iht_above,ilat_south,ilon_west)
         noplus_below_south_west = noplus_high_res_fixed(iht_below,ilat_south,ilon_west)
         noplus_above_south_east = noplus_high_res_fixed(iht_above,ilat_south,ilon_east)
         noplus_below_south_east = noplus_high_res_fixed(iht_below,ilat_south,ilon_east)

         noplus_north_west = ((noplus_above_north_west - noplus_below_north_west)*factor_ht) + noplus_below_north_west
         noplus_north_east = ((noplus_above_north_east - noplus_below_north_east)*factor_ht) + noplus_below_north_east
         noplus_south_west = ((noplus_above_south_west - noplus_below_south_west)*factor_ht) + noplus_below_south_west
         noplus_south_east = ((noplus_above_south_east - noplus_below_south_east)*factor_ht) + noplus_below_south_east

         n2plus_above_north_west = n2plus_high_res_fixed(iht_above,ilat_north,ilon_west)
         n2plus_below_north_west = n2plus_high_res_fixed(iht_below,ilat_north,ilon_west)
         n2plus_above_north_east = n2plus_high_res_fixed(iht_above,ilat_north,ilon_east)
         n2plus_below_north_east = n2plus_high_res_fixed(iht_below,ilat_north,ilon_east)
         n2plus_above_south_west = n2plus_high_res_fixed(iht_above,ilat_south,ilon_west)
         n2plus_below_south_west = n2plus_high_res_fixed(iht_below,ilat_south,ilon_west)
         n2plus_above_south_east = n2plus_high_res_fixed(iht_above,ilat_south,ilon_east)
         n2plus_below_south_east = n2plus_high_res_fixed(iht_below,ilat_south,ilon_east)

         n2plus_north_west = ((n2plus_above_north_west - n2plus_below_north_west)*factor_ht) + n2plus_below_north_west
         n2plus_north_east = ((n2plus_above_north_east - n2plus_below_north_east)*factor_ht) + n2plus_below_north_east
         n2plus_south_west = ((n2plus_above_south_west - n2plus_below_south_west)*factor_ht) + n2plus_below_south_west
         n2plus_south_east = ((n2plus_above_south_east - n2plus_below_south_east)*factor_ht) + n2plus_below_south_east

         o2plus_above_north_west = o2plus_high_res_fixed(iht_above,ilat_north,ilon_west)
         o2plus_below_north_west = o2plus_high_res_fixed(iht_below,ilat_north,ilon_west)
         o2plus_above_north_east = o2plus_high_res_fixed(iht_above,ilat_north,ilon_east)
         o2plus_below_north_east = o2plus_high_res_fixed(iht_below,ilat_north,ilon_east)
         o2plus_above_south_west = o2plus_high_res_fixed(iht_above,ilat_south,ilon_west)
         o2plus_below_south_west = o2plus_high_res_fixed(iht_below,ilat_south,ilon_west)
         o2plus_above_south_east = o2plus_high_res_fixed(iht_above,ilat_south,ilon_east)
         o2plus_below_south_east = o2plus_high_res_fixed(iht_below,ilat_south,ilon_east)

         o2plus_north_west = ((o2plus_above_north_west - o2plus_below_north_west)*factor_ht) + o2plus_below_north_west
         o2plus_north_east = ((o2plus_above_north_east - o2plus_below_north_east)*factor_ht) + o2plus_below_north_east
         o2plus_south_west = ((o2plus_above_south_west - o2plus_below_south_west)*factor_ht) + o2plus_below_south_west
         o2plus_south_east = ((o2plus_above_south_east - o2plus_below_south_east)*factor_ht) + o2plus_below_south_east

         nplus_above_north_west = nplus_high_res_fixed(iht_above,ilat_north,ilon_west)
         nplus_below_north_west = nplus_high_res_fixed(iht_below,ilat_north,ilon_west)
         nplus_above_north_east = nplus_high_res_fixed(iht_above,ilat_north,ilon_east)
         nplus_below_north_east = nplus_high_res_fixed(iht_below,ilat_north,ilon_east)
         nplus_above_south_west = nplus_high_res_fixed(iht_above,ilat_south,ilon_west)
         nplus_below_south_west = nplus_high_res_fixed(iht_below,ilat_south,ilon_west)
         nplus_above_south_east = nplus_high_res_fixed(iht_above,ilat_south,ilon_east)
         nplus_below_south_east = nplus_high_res_fixed(iht_below,ilat_south,ilon_east)

         nplus_north_west = ((nplus_above_north_west - nplus_below_north_west)*factor_ht) + nplus_below_north_west
         nplus_north_east = ((nplus_above_north_east - nplus_below_north_east)*factor_ht) + nplus_below_north_east
         nplus_south_west = ((nplus_above_south_west - nplus_below_south_west)*factor_ht) + nplus_below_south_west
         nplus_south_east = ((nplus_above_south_east - nplus_below_south_east)*factor_ht) + nplus_below_south_east

         Te_above_north_west = Te_high_res_fixed(iht_above,ilat_north,ilon_west)
         Te_below_north_west = Te_high_res_fixed(iht_below,ilat_north,ilon_west)
         Te_above_north_east = Te_high_res_fixed(iht_above,ilat_north,ilon_east)
         Te_below_north_east = Te_high_res_fixed(iht_below,ilat_north,ilon_east)
         Te_above_south_west = Te_high_res_fixed(iht_above,ilat_south,ilon_west)
         Te_below_south_west = Te_high_res_fixed(iht_below,ilat_south,ilon_west)
         Te_above_south_east = Te_high_res_fixed(iht_above,ilat_south,ilon_east)
         Te_below_south_east = Te_high_res_fixed(iht_below,ilat_south,ilon_east)

         Te_north_west = ((Te_above_north_west - Te_below_north_west)*factor_ht) + Te_below_north_west
         Te_north_east = ((Te_above_north_east - Te_below_north_east)*factor_ht) + Te_below_north_east
         Te_south_west = ((Te_above_south_west - Te_below_south_west)*factor_ht) + Te_below_south_west
         Te_south_east = ((Te_above_south_east - Te_below_south_east)*factor_ht) + Te_below_south_east

         Ti1_above_north_west = Ti1_high_res_fixed(iht_above,ilat_north,ilon_west)
         Ti1_below_north_west = Ti1_high_res_fixed(iht_below,ilat_north,ilon_west)
         Ti1_above_north_east = Ti1_high_res_fixed(iht_above,ilat_north,ilon_east)
         Ti1_below_north_east = Ti1_high_res_fixed(iht_below,ilat_north,ilon_east)
         Ti1_above_south_west = Ti1_high_res_fixed(iht_above,ilat_south,ilon_west)
         Ti1_below_south_west = Ti1_high_res_fixed(iht_below,ilat_south,ilon_west)
         Ti1_above_south_east = Ti1_high_res_fixed(iht_above,ilat_south,ilon_east)
         Ti1_below_south_east = Ti1_high_res_fixed(iht_below,ilat_south,ilon_east)

         Ti1_north_west = ((Ti1_above_north_west - Ti1_below_north_west)*factor_ht) + Ti1_below_north_west
         Ti1_north_east = ((Ti1_above_north_east - Ti1_below_north_east)*factor_ht) + Ti1_below_north_east
         Ti1_south_west = ((Ti1_above_south_west - Ti1_below_south_west)*factor_ht) + Ti1_below_south_west
         Ti1_south_east = ((Ti1_above_south_east - Ti1_below_south_east)*factor_ht) + Ti1_below_south_east

         Ti2_above_north_west = Ti2_high_res_fixed(iht_above,ilat_north,ilon_west)
         Ti2_below_north_west = Ti2_high_res_fixed(iht_below,ilat_north,ilon_west)
         Ti2_above_north_east = Ti2_high_res_fixed(iht_above,ilat_north,ilon_east)
         Ti2_below_north_east = Ti2_high_res_fixed(iht_below,ilat_north,ilon_east)
         Ti2_above_south_west = Ti2_high_res_fixed(iht_above,ilat_south,ilon_west)
         Ti2_below_south_west = Ti2_high_res_fixed(iht_below,ilat_south,ilon_west)
         Ti2_above_south_east = Ti2_high_res_fixed(iht_above,ilat_south,ilon_east)
         Ti2_below_south_east = Ti2_high_res_fixed(iht_below,ilat_south,ilon_east)

         Ti2_north_west = ((Ti2_above_north_west - Ti2_below_north_west)*factor_ht) + Ti2_below_north_west
         Ti2_north_east = ((Ti2_above_north_east - Ti2_below_north_east)*factor_ht) + Ti2_below_north_east
         Ti2_south_west = ((Ti2_above_south_west - Ti2_below_south_west)*factor_ht) + Ti2_below_south_west
         Ti2_south_east = ((Ti2_above_south_east - Ti2_below_south_east)*factor_ht) + Ti2_below_south_east

         !--------------------------------------------
         ! interpolate neutral heating rates
         !--------------------------------------------
      if (useIPEHeatingRates) then
         heatingRate_above_north_west = heatingRate_high_res_fixed(iht_above,ilat_north,ilon_west)
         heatingRate_below_north_west = heatingRate_high_res_fixed(iht_below,ilat_north,ilon_west)
         heatingRate_above_north_east = heatingRate_high_res_fixed(iht_above,ilat_north,ilon_east)
         heatingRate_below_north_east = heatingRate_high_res_fixed(iht_below,ilat_north,ilon_east)
         heatingRate_above_south_west = heatingRate_high_res_fixed(iht_above,ilat_south,ilon_west)
         heatingRate_below_south_west = heatingRate_high_res_fixed(iht_below,ilat_south,ilon_west)
         heatingRate_above_south_east = heatingRate_high_res_fixed(iht_above,ilat_south,ilon_east)
         heatingRate_below_south_east = heatingRate_high_res_fixed(iht_below,ilat_south,ilon_east)

         heatingRate_north_west = ((heatingRate_above_north_west - heatingRate_below_north_west)*factor_ht) + heatingRate_below_north_west
         heatingRate_north_east = ((heatingRate_above_north_east - heatingRate_below_north_east)*factor_ht) + heatingRate_below_north_east
         heatingRate_south_west = ((heatingRate_above_south_west - heatingRate_below_south_west)*factor_ht) + heatingRate_below_south_west
         heatingRate_south_east = ((heatingRate_above_south_east - heatingRate_below_south_east)*factor_ht) + heatingRate_below_south_east
     end if


         !------------------------------
         ! latitude interpolation....
         !------------------------------
         ne_east =    ((ne_north_east - ne_south_east) * factor_lat) + ne_south_east
         ne_west =     ((ne_north_west - ne_south_west) * factor_lat) + ne_south_west
         oplus_east =  ((oplus_north_east - oplus_south_east) * factor_lat) + oplus_south_east
         oplus_west =  ((oplus_north_west - oplus_south_west) * factor_lat) + oplus_south_west
         hplus_east =  ((hplus_north_east - hplus_south_east) * factor_lat) + hplus_south_east
         hplus_west =  ((hplus_north_west - hplus_south_west) * factor_lat) + hplus_south_west
         noplus_east = ((noplus_north_east - noplus_south_east) * factor_lat) + noplus_south_east
         noplus_west = ((noplus_north_west - noplus_south_west) * factor_lat) + noplus_south_west
         o2plus_east = ((o2plus_north_east - o2plus_south_east) * factor_lat) + o2plus_south_east
         o2plus_west = ((o2plus_north_west - o2plus_south_west) * factor_lat) + o2plus_south_west
         n2plus_east = ((n2plus_north_east - n2plus_south_east) * factor_lat) + n2plus_south_east
         n2plus_west = ((n2plus_north_west - n2plus_south_west) * factor_lat) + n2plus_south_west
         nplus_east =  ((nplus_north_east - nplus_south_east) * factor_lat) + nplus_south_east
         nplus_west =  ((nplus_north_west - nplus_south_west) * factor_lat) + nplus_south_west
         Te_east =     ((Te_north_east - Te_south_east) * factor_lat) + Te_south_east
         Te_west =     ((Te_north_west - Te_south_west) * factor_lat) + Te_south_west
         Ti1_east =    ((Ti1_north_east - Ti1_south_east) * factor_lat) + Ti1_south_east
         Ti1_west =    ((Ti1_north_west - Ti1_south_west) * factor_lat) + Ti1_south_west
         Ti2_east =    ((Ti2_north_east - Ti2_south_east) * factor_lat) + Ti2_south_east
         Ti2_west =    ((Ti2_north_west - Ti2_south_west) * factor_lat) + Ti2_south_west


       if (useIPEHeatingRates) then
         heatingRate_east =     ((heatingRate_north_east - heatingRate_south_east) * factor_lat) + heatingRate_south_east
         heatingRate_west =     ((heatingRate_north_west - heatingRate_south_west) * factor_lat) + heatingRate_south_west
       end if


         !--------------------------------
         ! longitude interpolation....
         !--------------------------------
         therm_Ne_density(iht, ilat, ilon) =     ((ne_east - ne_west) * factor_lon) + ne_west
         therm_oplus_density(iht, ilat, ilon) =  ((oplus_east - oplus_west) * factor_lon) + oplus_west
         therm_hplus_density(iht, ilat, ilon) =  ((hplus_east - hplus_west) * factor_lon) + hplus_west
         therm_noplus_density(iht, ilat, ilon) = ((noplus_east - noplus_west) * factor_lon) + noplus_west
         therm_o2plus_density(iht, ilat, ilon) = ((o2plus_east - o2plus_west) * factor_lon) + o2plus_west
         therm_n2plus_density(iht, ilat, ilon) = ((n2plus_east - n2plus_west) * factor_lon) + n2plus_west
         therm_nplus_density(iht, ilat, ilon) =  ((nplus_east - nplus_west) * factor_lon) + nplus_west
         therm_Te(iht, ilat, ilon) =             ((Te_east - Te_west) * factor_lon) + Te_west
         therm_Ti1(iht, ilat, ilon) =            ((Ti1_east - Ti1_west) * factor_lon) + Ti1_west
         therm_Ti2(iht, ilat, ilon) =            ((Ti2_east - Ti2_west) * factor_lon) + Ti2_west

       if (useIPEHeatingRates) then
         neutralHeatingRates(iht, ilat, ilon) =             ((heatingRate_east - heatingRate_west) * factor_lon) + heatingRate_west
       endif 

       enddo ! iht
   enddo ! ilat
enddo ! ilon


return




end SUBROUTINE INTERFACE__FIXED_GRID_to_THERMO

!********************************************************************************************************





END MODULE moduleInterfaceIono2Thermo
