MODULE IONOSPHERE_PLASMASPHERE

  IMPLICIT NONE

  include 'netcdf.inc'

  PUBLIC :: GIP_init
  PUBLIC :: GIP_calculation
  PRIVATE

  SAVE

  INTEGER NPTS
  INTEGER NMP
  INTEGER NLP
  PARAMETER (NPTS = 13813)
  PARAMETER (NMP  = 80)
  PARAMETER (NLP  = 67)

  INTEGER which_ion
  INTEGER high_lat_hts
  INTEGER high_lat_lats
  INTEGER high_lat_lons
  INTEGER interface_hts
  INTEGER fixed_height_lats
  INTEGER fixed_height_lons
  PARAMETER (which_ion = 2)
  PARAMETER (high_lat_hts  = 90)
  PARAMETER (high_lat_lats = 91)
  PARAMETER (high_lat_lons = 20)
  PARAMETER (interface_hts  = 31)
  PARAMETER (fixed_height_lats = 91)
  PARAMETER (fixed_height_lons = 90)

  INTEGER :: jion(6) 
  INTEGER :: ii1_interface(3,interface_hts,91,90)
  INTEGER :: ii2_interface(3,interface_hts,91,90)
  INTEGER :: ii3_interface(3,interface_hts,91,90)
  INTEGER :: ii4_interface(3,interface_hts,91,90)
  REAL(kind=8) facfac_interface(3,interface_hts,91,90)
  REAL(kind=8) dd_interface(3,interface_hts,91,90)

  REAL(kind=8) fixed_heights_km(interface_hts)

  INTEGER :: MLOw_interface(90)
  INTEGER :: MHIgh_interface(90)
  INTEGER :: in(nmp,nlp) , is(nmp,nlp)
  INTEGER :: L300 , L1000
  INTEGER :: nhgt , j

  REAL(kind=8) q_coordinate_plasma(NPTS,NMP)
  REAL(kind=8) gr_plasma(NPTS,NMP)
  REAL(kind=8) gcol_plasma(NPTS,NMP)
  REAL(kind=8) glon_plasma(NPTS,NMP)
  REAL(kind=8) bcol_plasma(NPTS,NMP)
  REAL(kind=8) blon_plasma(NMP,NLP)
  REAL(kind=8) Apex_D1(3,NPTS,NMP)
  REAL(kind=8) Apex_D2(3,NPTS,NMP)
  REAL(kind=8) Apex_E1(3,NPTS,NMP)
  REAL(kind=8) Apex_E2(3,NPTS,NMP)
  REAL(kind=8) Apex_BE3(NPTS,NMP)
  REAL(kind=8) Apex_grdlbm2(3,NPTS,NMP)
  REAL(kind=8) ETA_apex_3d(NPTS,NMP)
  REAL(kind=8) Apex_BHAT(3,npts,nmp)
  REAL(kind=8) Apex_D(npts,nmp)
  REAL(kind=8) Apex_d1d1(npts,nmp)
  REAL(kind=8) Apex_d1d2(npts,nmp)
  REAL(kind=8) Apex_d2d2(npts,nmp)
  REAL(kind=8) Apex_BMAG(npts,nmp)
  REAL(kind=8) re_apex_plasma(nmp,nlp)
  REAL(kind=8) m_plasma(6)
  REAL(kind=8) hprof(20,19)
  REAL(kind=8) Glat_plasma_3d(npts,nmp)
  REAL(kind=8) Glond_plasma_3d(npts,nmp)
  REAL(kind=8) Pz_plasma_3d(npts,nmp)
  REAL(kind=8) km_plasma(6)
  REAL(kind=8) r0
  REAL(kind=8) zkm(90)
  REAL(kind=8) hnew(90) , h2(90) , h2diff(90) , hprod2(90) , hprod3(90)
  REAL(kind=8) hz_metres(90)
  REAL(kind=8) tarea(90)
  REAL(kind=8) dth_radians
  REAL(kind=8) gravity(90)
  REAL(kind=8) btotal(91,20)
  REAL(kind=8) dip_angle_radians(91,20)
  REAL(kind=8) declination(91,20)

  INTEGER ::   ilon1_3d(npts,nmp),ilon2_3d(npts,nmp)
  INTEGER ::   ilat1_3d(npts,nmp),ilat2_3d(npts,nmp)
  INTEGER ::   ispecial_3d(npts,nmp)
  INTEGER ::   ihl11_3d(npts,nmp),ihu11_3d(npts,nmp)
  INTEGER ::   ihl12_3d(npts,nmp),ihu12_3d(npts,nmp)
  INTEGER ::   ihl21_3d(npts,nmp),ihu21_3d(npts,nmp)
  INTEGER ::   ihl22_3d(npts,nmp),ihu22_3d(npts,nmp)
  LOGICAL ::   sw_1st_call_int_fixed_ht
  INTEGER ::   ilon1_3d_fixed_ht(npts,nmp),ilon2_3d_fixed_ht(npts,nmp)
  INTEGER ::   ilat1_3d_fixed_ht(npts,nmp),ilat2_3d_fixed_ht(npts,nmp)
  INTEGER ::   ispecial_3d_fixed_ht(npts,nmp)
  INTEGER ::   ihl_3d_fixed_ht(npts,nmp),ihu_3d_fixed_ht(npts,nmp)

  REAL(kind=8) ti_plasma(NPTS,NMP,2)
  REAL(kind=8) te_plasma(NPTS,NMP)
  REAL(kind=8) ni_plasma(NPTS,NMP,2)
  REAL(kind=8) vi_plasma(NPTS,NMP,2)
  REAL(kind=8) ne_plasma(NPTS,NMP)
  REAL(kind=8) yni(NPTS,NMP,2)
  REAL(kind=8) yvi(NPTS,NMP,2)
  REAL(kind=8) yne(NPTS,NMP)
  REAL(kind=8) yqe(NPTS,NMP)
  REAL(kind=8) d13d(90,91,20,2)
  REAL(kind=8) d23d(90,91,20,2)
  REAL(kind=8) v13d(90,91,20,2)
  REAL(kind=8) v23d(90,91,20,2)
  REAL(kind=8) ti3d(90,91,20)
  REAL(kind=8) ti1_high_lat_ions(90,91,20)
  REAL(kind=8) ti2_high_lat_ions(90,91,20)
  REAL(kind=8) te_high_lat_ions(90,91,20)
  REAL(kind=8) no_plus_3d(NPTS,NMP)
  REAL(kind=8) o2_plus_3d(NPTS,NMP)
  REAL(kind=8) n2_plus_3d(NPTS,NMP)
  REAL(kind=8) n_plus_3d(NPTS,NMP)
  REAL(kind=8) vpeq(NMP,NLP)
  REAL(kind=8) vzon(NMP,NLP)

  INTEGER :: i_first_call_of_plasma
  REAL(kind=8) factor_ht_3d(90,91,20)
  INTEGER iht_above_3d(90,91,20)
  INTEGER iht_below_3d(90,91,20)
  REAL(kind=8) factor_ht_inverse_3d(interface_hts,91,20)
  INTEGER iht_above_inverse_3d(interface_hts,91,20)
  INTEGER iht_below_inverse_3d(interface_hts,91,20)
  REAL(kind=8) altitude_metres_3d(90,91,20)

  REAL(kind=8) UT_previous_call_plasma_secs
  REAL(kind=8) UT_previous_call_polar_secs
  INTEGER :: MASs
  INTEGER :: NDT
  INTEGER :: NIOns
  INTEGER :: midpoint(nlp)
  INTEGER :: ITDay
  REAL(kind=8) :: QB

  REAL(kind=8) :: angdif(91,20)

  REAL(kind=8) :: Hyd_grid_400km_m3(25,19)
  real(kind=8) :: Hyd_grid_lats(19)
  real(kind=8) :: Hyd_grid_local_times(25)

data Hyd_grid_lats/-90.,-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10.,0.,10.,20.,30.,40.,50.,60.,70.,80.,90./
data Hyd_grid_local_times/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20., &
                          21.,22.,23.,24./

  PARAMETER (R0=6.370E06)

  DATA fixed_heights_km/90.,95.,100.,105.,110.,115.,120.,125., &
         150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400., &
         450.,500.,550.,600.,700.,800.,900.,1000., &
         2000.,4000.,6370.,9000./

CONTAINS





!***********************************************************************************



SUBROUTINE INTERFACE__IPE_to_thermosphere ( &
         thermospheric_model_name , ht_dim , lat_dim , lon_dim , &
         ne_high_res_fixed,oplus_high_res_fixed,hplus_high_res_fixed, &
         noplus_high_res_fixed,o2plus_high_res_fixed, &
         n2plus_high_res_fixed,nplus_high_res_fixed, &
         Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed, &
         therm_geo_long_input,therm_geo_lat,therm_Z, &
         therm_Ne_density,therm_oplus_density,therm_hplus_density, &
         therm_noplus_density,therm_o2plus_density, &
         therm_n2plus_density,therm_nplus_density, &
         therm_Te,therm_Ti1,therm_Ti2)

IMPLICIT NONE

    character*10 thermospheric_model_name
    integer ht_dim , lat_dim , lon_dim
    INTEGER :: ispecial
    REAL(kind=8) :: ne_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: oplus_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: hplus_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: noplus_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: o2plus_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: n2plus_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: nplus_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: Te_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: Ti1_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: Ti2_high_res_fixed(interface_hts,91,90)
    REAL(kind=8) :: therm_geo_long_input(lon_dim)
    REAL(kind=8) :: therm_geo_long(lon_dim)
    REAL(kind=8) :: therm_geo_lat(lat_dim)
    REAL(kind=8) :: therm_Z(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_Z_km

    REAL(kind=8) :: therm_Ne_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_oplus_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_hplus_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_noplus_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_o2plus_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_n2plus_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_nplus_density(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_Te(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_Ti1(ht_dim,lon_dim,lat_dim)
    REAL(kind=8) :: therm_Ti2(ht_dim,lon_dim,lat_dim)

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

    REAL(kind=8) :: high_res_long(90) 
    REAL(kind=8) :: high_res_lat(91) 
    REAL(kind=8) :: high_res_height(interface_hts)

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

    therm_geo_long(:) = therm_geo_long_input(:)

    do ilon = 1 , lon_dim
      if ( therm_geo_long(ilon) < 0.0 ) then
         therm_geo_long(ilon) = therm_geo_long(ilon) + 360.
      endif
    enddo

    do ilon_int = 1 , 90
      high_res_long(ilon_int)=(float(ilon_int-1)) * 4.
    enddo
    do ilat_int = 1 , 91
      high_res_lat(ilat_int)=(float(ilat_int-46)) * 2.
    enddo
    do iht_int = 1 , interface_hts
      high_res_height(iht_int)=fixed_heights_km(iht_int)
    enddo

! Loop over therm longitudes....

    do ilon = 1 , lon_dim 

!longitude interpolation

! loop over high_res longs to find points east and west....

       ispecial = 0
    do ilon_int = 1 , 90

      if (high_res_long(ilon_int) > therm_geo_long(ilon)) then
        ilon_east = ilon_int
        ilon_west = ilon_int - 1
        goto 1500
      endif
    enddo
       ilon_east = 1
       ilon_west = 90
       ispecial = 1
       
1500 continue

    if(ilon_east == 1) then
       ilon_east = 1
       ilon_west = 90
       ispecial = 1
    endif
    if ( ispecial == 0 ) then
    factor_lon = (therm_geo_long(ilon) - high_res_long(ilon_west)) /  &
                 (high_res_long(ilon_east) - high_res_long(ilon_west))
    else

    factor_lon = (therm_geo_long(ilon) - high_res_long(ilon_west)) /  &
                 (high_res_long(ilon_east)+360. - high_res_long(ilon_west))
    endif

    if ( factor_lon > 1.0.or. factor_lon < 0.0 ) then
        write(6,*) ' factor lon ',factor_lon,therm_geo_long(ilon),high_res_long(ilon_west), &
        high_res_long(ilon_east)
    endif

    factor_lon_array(ilon) = factor_lon
    ilon_west_array(ilon) = ilon_west
    ilon_east_array(ilon) = ilon_east

! Loop over therm latitudes....

    do ilat = 1 , lat_dim

! latitude interpolation....

    do ilat_int = 1 , 91

      if(high_res_lat(ilat_int) > therm_geo_lat(ilat)) then 
        ilat_north = ilat_int
        ilat_south = ilat_int - 1
        goto 2000
      endif
    enddo
2000 continue

    if(ilat_north == 1) then
       ilat_north = 2
       ilat_south = 1
    endif
    factor_lat = (therm_geo_lat(ilat) - high_res_lat(ilat_south)) /  &
                 (high_res_lat(ilat_north) - high_res_lat(ilat_south))

    factor_lat_array(ilat) = factor_lat
    ilat_north_array(ilat) = ilat_north
    ilat_south_array(ilat) = ilat_south



! Loop over therm heights....

  do iht = 1 , ht_dim           

  therm_Z_km = therm_Z(iht,ilon,ilat) / 1000.

! height interpolation....

  do iht_int = 1 , interface_hts         
    if(high_res_height(iht_int) > therm_Z_km) then
      iht_above = iht_int
      iht_below = iht_int - 1
      goto 2500    
    endif
  enddo
  
  iht_above = interface_hts
  iht_below = interface_hts - 1
  
2500 continue


if (iht_above == 1) then
    iht_above = 2
    iht_below = 1
endif

factor_ht = (therm_Z_km - high_res_height(iht_below)) /  &
            (high_res_height(iht_above) - high_res_height(iht_below))

! cg - make sure factor_ht doesn't get smaller than 0.0 which will happen for thermospheric
! cg - heights below 90km.....

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


!--------------------------------
! longitude interpolation....
!--------------------------------
therm_Ne_density(iht,ilon,ilat) =     ((ne_east - ne_west) * factor_lon) + ne_west
therm_oplus_density(iht,ilon,ilat) =  ((oplus_east - oplus_west) * factor_lon) + oplus_west
therm_hplus_density(iht,ilon,ilat) =  ((hplus_east - hplus_west) * factor_lon) + hplus_west
therm_noplus_density(iht,ilon,ilat) = ((noplus_east - noplus_west) * factor_lon) + noplus_west
therm_o2plus_density(iht,ilon,ilat) = ((o2plus_east - o2plus_west) * factor_lon) + o2plus_west
therm_n2plus_density(iht,ilon,ilat) = ((n2plus_east - n2plus_west) * factor_lon) + n2plus_west
therm_nplus_density(iht,ilon,ilat) =  ((nplus_east - nplus_west) * factor_lon) + nplus_west
therm_Te(iht,ilon,ilat) =             ((Te_east - Te_west) * factor_lon) + Te_west
therm_Ti1(iht,ilon,ilat) =            ((Ti1_east - Ti1_west) * factor_lon) + Ti1_west
therm_Ti2(iht,ilon,ilat) =            ((Ti2_east - Ti2_west) * factor_lon) + Ti2_west

enddo
enddo
enddo


return




end SUBROUTINE INTERFACE__IPE_to_thermosphere

!********************************************************************************************************



