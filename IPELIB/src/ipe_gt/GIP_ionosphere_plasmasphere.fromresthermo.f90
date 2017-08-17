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





SUBROUTINE gip_init(GIP_switches, iday, uthr, f107, &
                        GIP_input_dataset, GIP_output_dataset, &
                        static_file_location, GIP_Apex_coords_static_file, &
                        thermospheric_model_name, &
                        therm_model_ht_dim, &
                        therm_model_lat_dim, &
                        therm_model_lon_dim, &
                        therm_model_geo_long, &
                        therm_model_geo_lat, &
                        therm_model_ht_m, &
                        therm_model_Ne_density, &
                        therm_model_oplus_density, &
                        therm_model_hplus_density, &
                        therm_model_noplus_density, &
                        therm_model_o2plus_density, &
                        therm_model_n2plus_density, &
                        therm_model_nplus_density, &
                        therm_model_Te, &
                        therm_model_Ti1, &
                        therm_model_Ti2, &
                        dynamo_sigma_phph_dsi, &
                        dynamo_sigma_lmlm_msi, &
                        dynamo_sigma_h, &
                        dynamo_sigma_c, &
                        dynamo_Kdmph_dsi, &
                        dynamo_Kdmlm)
!
! global variables from TGCM      
    implicit none

    integer, intent(in) :: iday
    real(kind=8), intent(in) :: uthr,f107
    character(100) :: GIP_input_dataset , GIP_output_dataset
    character(100) :: static_file_location
    character(100) :: GIP_Apex_coords_static_file
    logical :: GIP_switches(20)
    INTEGER :: idump_GIP
    INTEGER :: it

    REAL(kind=8) :: potential_field(81,97)
    REAL(kind=8) :: ed1(81,97)
    REAL(kind=8) :: ed2(81,97)

    REAL(kind=8) :: dynamo_sigma_phph_dsi(81,97), &
            dynamo_sigma_lmlm_msi(81,97), &
            dynamo_sigma_h(81,97),dynamo_sigma_c(81,97), &
            dynamo_Kdmph_dsi(81,97),dynamo_Kdmlm(81,97)

    CHARACTER*40 header
    LOGICAL :: sw_initialisation_call
    INTEGER :: i_call_polar
    INTEGER :: i_call_plasma
    INTEGER :: iday_number
    INTEGER :: istop
    REAL(kind=8) :: universal_time_seconds
    REAL(kind=8) :: geo_grid_longitudes_degrees(20)
    REAL(kind=8) :: geo_grid_latitudes_degrees(91)
    REAL(kind=8) :: dummy_1d_20(20)
    REAL(kind=8) :: dummy_1d_91(91)
    REAL(kind=8) :: ex2d(91,20)
    REAL(kind=8) :: ey2d(91,20)

    REAL(kind=8) ::  o_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  O2_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  N2_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  NO_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  N4S_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  N2D_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  Vx_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  Vy_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  Wvz_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  tts_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qion3d_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  elx_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  ely_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qo2p_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qop_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qn2p_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qnp_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qtef_aurora_fixed_ht(interface_hts,91,20)


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


    real(kind=8) :: NmF2(91,90)
    real(kind=8) :: HmF2_km(91,90)
    real(kind=8) :: TEC(91,90)
    real(kind=8) :: Te_400(91,90)
    real(kind=8) :: Ti_400(91,90)

    real(kind=8) :: this_Ne_profile(interface_hts)

      integer :: therm_model_ht_dim , therm_model_lat_dim , therm_model_lon_dim
      character*10 :: thermospheric_model_name
      real(kind=8) :: therm_model_geo_long(therm_model_lon_dim) , therm_model_geo_lat(therm_model_lat_dim)
      real(kind=8), dimension(therm_model_ht_dim,therm_model_lon_dim,therm_model_lat_dim) ::  &
                    therm_model_ht_m, &
                    therm_model_Ne_density,therm_model_oplus_density, &
                    therm_model_hplus_density, therm_model_noplus_density,therm_model_o2plus_density, &
                    therm_model_n2plus_density, therm_model_nplus_density, &
                    therm_model_Te, therm_model_Ti1, therm_model_Ti2

!

     iday_number            = iday
     universal_time_seconds = uthr*3600.

    OPEN (22,FILE=TRIM(static_file_location)//'hprof',STATUS='old')
    OPEN (174,FILE=TRIM(static_file_location)//'angdif_for_GIP',STATUS='old')
    OPEN (175,FILE=TRIM(static_file_location)//'btotal_dip_declination_old_dipole',STATUS='old')

    OPEN (51,FILE=GIP_Apex_coords_static_file,STATUS='OLD')

    call get_HYD400km_for_this_run(static_file_location,&
                                   iday_number,f107)

!                                  Hyd_grid_400km_m3,Hyd_grid_lats, &
!                                  Hyd_grid_local_times)

     write(6,*) '***************** HYDROGEN GRID AT 400KM ***********************************'
     write(6,1888) Hyd_grid_400km_m3 / 1.e11
     1888 format(25f4.1)
     write(6,*) '****************************************************************************'

    sw_initialisation_call = .TRUE.
    idump_GIP = 0
    i_call_polar = 1
    i_call_plasma = 1

    CALL GLOBAL_IONOSPHERE_PLASMASPHERE ( &
         GIP_switches, &
         i_call_polar,i_call_plasma, &
         idump_GIP, &
         GIP_input_dataset, GIP_output_dataset, &
         sw_initialisation_call, &
         iday_number, &
         universal_time_seconds, &
         geo_grid_longitudes_degrees, geo_grid_latitudes_degrees, &
         ex2d,ey2d, &
         potential_field,ed1,ed2, &
         f107, &
         o_density_fixed_ht,o2_density_fixed_ht,n2_density_fixed_ht, &
         NO_density_fixed_ht, &
         N4S_density_fixed_ht,N2D_density_fixed_ht, &
         vx_fixed_ht,vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
         qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
         qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
         oplus_high_res_fixed,hplus_high_res_fixed, &
         noplus_high_res_fixed,o2plus_high_res_fixed, &
         n2plus_high_res_fixed,nplus_high_res_fixed, &
         Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed, &
         ne_high_res_fixed, &
         NmF2,HmF2_km,TEC,Te_400,Ti_400,this_Ne_profile, &
         dynamo_sigma_phph_dsi,dynamo_sigma_lmlm_msi, &
         dynamo_sigma_h,dynamo_sigma_c, &
         dynamo_Kdmph_dsi,dynamo_Kdmlm)

    call INTERFACE__GIP_to_thermosphere ( &
         thermospheric_model_name , therm_model_ht_dim , therm_model_lat_dim , therm_model_lon_dim, &
         ne_high_res_fixed,oplus_high_res_fixed,hplus_high_res_fixed, &
         noplus_high_res_fixed,o2plus_high_res_fixed, &
         n2plus_high_res_fixed,nplus_high_res_fixed, &
         Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed, &
         therm_model_geo_long,therm_model_geo_lat,therm_model_ht_m, &
         therm_model_Ne_density,therm_model_oplus_density,therm_model_hplus_density, &
         therm_model_noplus_density,therm_model_o2plus_density, &
         therm_model_n2plus_density,therm_model_nplus_density, &
         therm_model_Te,therm_model_Ti1,therm_model_Ti2)

    return



end SUBROUTINE gip_init



SUBROUTINE GIP_CALCULATION ( &
                              GIP_switches, &
                              GIP_input_dataset, &
                              GIP_output_dataset, &
                              thermospheric_model_name, &
                              therm_model_ht_dim, &
                              therm_model_lat_dim, &
                              therm_model_lon_dim, &
                              therm_model_geo_long, &
                              therm_model_geo_lat, &
                              therm_model_ht_m, &
                              therm_model_O_density, &
                              therm_model_O2_density, &
                              therm_model_N2_density,  &
                              therm_model_NO_density,  &
                              therm_model_N4S_density,  &
                              therm_model_N2D_density,  &
                              therm_model_Tn, &
                              therm_model_Vx, &
                              therm_model_Vy, &
                              therm_model_Vz, &
                              therm_model_qion3d, &
                              therm_model_elx, &
                              therm_model_ely, &
                              therm_model_qo2p_aurora,  & ! ionization rates [1/m^3/s]
                              therm_model_qop_aurora, &
                              therm_model_qn2p_aurora, &
                              therm_model_qnp_aurora, &
                              therm_model_qtef_aurora, &
                              potential_field, &
                              ed1, &
                              ed2, &
                              iday_number, &
                              UT_hours, &
                              f107, &
                              idump_GIP, &
                              therm_model_Ne_density, &
                              therm_model_oplus_density, &
                              therm_model_hplus_density, &
                              therm_model_noplus_density, &
                              therm_model_o2plus_density, &
                              therm_model_n2plus_density, &
                              therm_model_nplus_density, &
                              therm_model_Te, &
                              therm_model_Ti1, &
                              therm_model_Ti2, &
                              dynamo_sigma_phph_dsi, &
                              dynamo_sigma_lmlm_msi, &
                              dynamo_sigma_h, &
                              dynamo_sigma_c, &
                              dynamo_Kdmph_dsi, &
                              dynamo_Kdmlm, &
                              ne_high_res_fixed)

    implicit none

    integer :: l,m,n
    integer :: therm_model_ht_dim , therm_model_lat_dim , therm_model_lon_dim
    character*10 :: thermospheric_model_name
    character*100 :: GIP_input_dataset
    character*100 :: GIP_output_dataset

    real(kind=8), intent(in) ::  UT_hours
!     integer, intent(in) :: iday,
!     real(kind=8), intent(in) :: f107
    integer :: UT_hours_part_integer
    integer :: UT_mins_part_integer
    character(2) :: hours_string
    character(2) :: mins_string

    logical :: polar,plasma
    logical :: GIP_switches(20)
    INTEGER :: idump_GIP

    real(kind=8), dimension(therm_model_ht_dim, therm_model_lon_dim, therm_model_lat_dim) ::  &
        therm_model_o_density, therm_model_o2_density, &
        therm_model_n2_density , therm_model_tn , therm_model_Vx , therm_model_Vy , therm_model_Vz , &
        therm_model_ht_m , therm_model_qion3d , therm_model_elx , therm_model_ely, &
        therm_model_qo2p_aurora, therm_model_qop_aurora, therm_model_qn2p_aurora, &
        therm_model_qnp_aurora, therm_model_qtef_aurora, &
        therm_model_NO_density, &
        therm_model_N4S_density, therm_model_N2D_density

    real(kind=8) :: therm_model_geo_long(therm_model_lon_dim) , therm_model_geo_lat(therm_model_lat_dim)
    character*40 header
    real(kind=8), dimension(therm_model_ht_dim,therm_model_lon_dim,therm_model_lat_dim) ::  &
                  therm_model_Ne_density,therm_model_oplus_density, &
                  therm_model_hplus_density, therm_model_noplus_density,therm_model_o2plus_density, &
                  therm_model_n2plus_density, therm_model_nplus_density, &
                  therm_model_Te, therm_model_Ti1, therm_model_Ti2

    REAL(kind=8) :: potential_field(81,97)
    REAL(kind=8) :: ed1(81,97)
    REAL(kind=8) :: ed2(81,97)

    REAL(kind=8) :: dynamo_sigma_phph_dsi(81,97), &
            dynamo_sigma_lmlm_msi(81,97), &
            dynamo_sigma_h(81,97),dynamo_sigma_c(81,97), &
            dynamo_Kdmph_dsi(81,97),dynamo_Kdmlm(81,97)

    real(kind=8) :: ne_big(interface_hts,91,90)
    real(kind=8) :: NmF2(91,90)
    real(kind=8) :: HmF2_km(91,90)
    real(kind=8) :: TEC(91,90)
    real(kind=8) :: this_Ne_profile(interface_hts)
    real(kind=8) :: Te_400(91,90)
    real(kind=8) :: Ti_400(91,90)
      REAL(kind=8) :: Upar_apex_at_F2(91,90)
      REAL(kind=8) :: Vx_at_F2(91,90)
      REAL(kind=8) :: Vy_at_F2(91,90)
      REAL(kind=8) :: O_density_at_F2(91,90)
      REAL(kind=8) :: N2_density_at_F2(91,90)
      REAL(kind=8) :: O_N2_ratio_F2(91,90)
    REAL(kind=8) :: Vx_300km(91,90)
    REAL(kind=8) :: Vy_300km(91,90)
    REAL(kind=8) :: O_density_300km(91,90)
    REAL(kind=8) :: N2_density_300km(91,90)
    REAL(kind=8) :: Upar_apex_300km(91,90)


    LOGICAL :: sw_initialisation_call
    INTEGER :: i_call_polar
    INTEGER :: i_call_plasma
    INTEGER :: iday_number
    INTEGER :: istop
    REAL(kind=8) :: universal_time_seconds
    REAL(kind=8) :: f107
    REAL(kind=8) :: geo_grid_longitudes_degrees(20)
    REAL(kind=8) :: geo_grid_latitudes_degrees(91)
    REAL(kind=8) :: dummy_1d_20(20)
    REAL(kind=8) :: dummy_1d_91(91)
    REAL(kind=8) :: ex2d(91,20)
    REAL(kind=8) :: ey2d(91,20)

    REAL(kind=8) ::  O_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  O2_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  N2_density_fixed_ht(interface_hts,91,20)

    REAL(kind=8) ::  NO_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  N4S_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  N2D_density_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  Vx_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  Vy_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  Wvz_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  tts_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qion3d_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qo2p_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qop_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qn2p_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qnp_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  qtef_aurora_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  elx_fixed_ht(interface_hts,91,20)
    REAL(kind=8) ::  ely_fixed_ht(interface_hts,91,20)

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


! lrm20111002
!------------------------------------------------------
! File unit number for checking thermosphere values :
!------------------------------------------------------
INTEGER, parameter :: unitNumber = 13
! Debug the Thermospheric values???
LOGICAL, parameter :: debugThermo = .FALSE.

! lrm added Oct 2, 2011 to compare these interpolated values
!------------------------------------------------------
! File unit number for checking interpolation values :
!------------------------------------------------------
INTEGER, parameter :: unitInterp = 14
! Write out the interpolated values??
LOGICAL, parameter :: writeInterp = .FALSE.
! check status of fortran open
INTEGER :: FileOpenStatus
INTEGER :: ii, jj

! BEGIN CODE ===========================================================================================

    polar = .TRUE.
    plasma = .TRUE.

    universal_time_seconds = UT_hours * 3600.

    UT_hours_part_integer = int(universal_time_seconds/3600.)
    UT_mins_part_integer = nint((((universal_time_seconds/3600.) - real(UT_hours_part_integer)) * 60.))

    if ( UT_hours_part_integer < 10 ) then
      write(hours_string,fmt='(i1)') UT_hours_part_integer
      hours_string = '0' // hours_string
    else
      write(hours_string,fmt='(i2)') UT_hours_part_integer
    endif

    if ( UT_mins_part_integer < 10 ) then
      write(mins_string,fmt='(i1)') UT_mins_part_integer
      mins_string = '0' // mins_string
    else
      write(mins_string,fmt='(i2)') UT_mins_part_integer
    endif

    write(6,*) '********* GIP CALLED AT ' // hours_string // ':' // mins_string // ' UT *********'

!      
    i_call_polar = 0    ! not necessary since at least the polar region is done
    i_call_plasma = 0
    if(polar) i_call_polar = 1
    if(plasma) i_call_plasma = 1
    sw_initialisation_call = .FALSE.

    call INTERFACE__thermosphere_to_GIP ( &
           GIP_switches, &
           thermospheric_model_name , &
           therm_model_ht_dim , &
           therm_model_lat_dim , &
           therm_model_lon_dim, &
           therm_model_o_density, &
           therm_model_o2_density, &
           therm_model_n2_density, &
           therm_model_NO_density, &
           therm_model_N4S_density, &
           therm_model_N2D_density, &
           therm_model_Tn, &
           therm_model_Vy, &
           therm_model_Vx, &
           therm_model_qion3d, &
           therm_model_elx, &
           therm_model_ely, &
           therm_model_qo2p_aurora, therm_model_qop_aurora, therm_model_qn2p_aurora, &
           therm_model_qnp_aurora, therm_model_qtef_aurora, &
           therm_model_geo_long, &
           therm_model_geo_lat, &
           therm_model_ht_m, &
           o_density_fixed_ht, o2_density_fixed_ht, n2_density_fixed_ht, &
           NO_density_fixed_ht, &
           N4S_density_fixed_ht, &
           N2D_density_fixed_ht, &
           Vx_fixed_ht, Vy_fixed_ht, wvz_fixed_ht, tts_fixed_ht, &
           qion3d_fixed_ht, elx_fixed_ht, ely_fixed_ht, &
           qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
           qnp_aurora_fixed_ht, qtef_aurora_fixed_ht)

!     write(6,*) '1/3 done'


      !-------------------------------------------------------------------
      ! If in Thermospheric debug mode, write out values to the a file
      ! & do some prints to the screen
      !-------------------------------------------------------------------
      if (debugThermo) then

           !========================================
           ! Open for checking values from GT 
           !========================================
           OPEN (UNIT=unitNumber, FILE='CheckGTGIP.dat', &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
           IF (FileOpenStatus > 0) STOP "Cannot open CheckGTGIP file ********"

          WRITE (unitNumber, '(A)' ) ' '
          WRITE(unitNumber,'(A)') 'GIP_CALCULATION : '
          write(unitNumber,'(A, 1F20.7)') 'GIP_CALCULATION : Universal_Time_seconds = ' , Universal_Time_seconds
          write(unitNumber,'(A)') 'GIP_CALCULATION : ----------------------------------------------------'

          WRITE (unitNumber, '(A)' ) ' '
          WRITE (unitNumber, '(A, 20F10.5)' ) 'therm_model_Vy(1,1,:) = ', &
                      therm_model_Vy(1,1,:)
          WRITE (unitNumber, '(A, 20F10.5)' ) 'therm_model_Vy(11,51,:) = ', &
                      therm_model_Vy(11,51,:)
          WRITE (unitNumber, '(A, 20F10.5)' ) 'therm_model_Vy(15,91,:) = ', &
                      therm_model_Vy(15,91,:)
          WRITE (unitNumber, '(A)' ) ' '
          WRITE (unitNumber, '(A, 20F10.5)' ) 'therm_model_Vx(1,1,:) = ', &
                      therm_model_Vx(1,1,:)
          WRITE (unitNumber, '(A, 20F10.5)' ) 'therm_model_Vx(11,51,:) = ', &
                      therm_model_Vx(11,51,:)
          WRITE (unitNumber, '(A, 20F10.5)' ) 'therm_model_Vx(15,91,:) = ', &
                      therm_model_Vx(15,91,:)

          !WRITE (unitNumber, '(A)' ) ' '
          !WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(1,1,:) = ', wvz_FROM_GT(1,1,:)
          !WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(11,51,:) = ', wvz_FROM_GT(11,51,:)
          !WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(15,91,:) = ', wvz_FROM_GT(15,91,:)

          !WRITE (unitNumber, '(A)' ) ' '
          !WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(1,1,:) = ',rmt_FROM_GT(1,1,:)
          !WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(11,51,:) = ',rmt_FROM_GT(11,51,:)
          !WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(15,91,:) = ',rmt_FROM_GT(15,91,:)

          WRITE (unitNumber, '(A)' ) ' '
          WRITE (unitNumber, '(A, 20F15.5)' ) 'therm_model_Tn(1,1,:) = ',therm_model_Tn(1,1,:)
          WRITE (unitNumber, '(A, 20F15.5)' ) 'therm_model_Tn(11,51,:) = ',therm_model_Tn(11,51,:)
          WRITE (unitNumber, '(A, 20F15.5)' ) 'therm_model_Tn(15,91,:) = ',therm_model_Tn(15,91,:)


          !WRITE (unitNumber, '(A)' ) ' '
          !WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(1,1,:) = ',ht_FROM_GT(1,1,:)
          !WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(11,51,:) = ',ht_FROM_GT(11,51,:)
          !WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(15,91,:) = ',ht_FROM_GT(15,91,:)


        endif ! debugThermo



           !-------------------------------------------------------------------------
           ! Write out results of the interpolation to an ascii file for examination
           !-------------------------------------------------------------------------
           if (writeInterp) then

              print *,'GIP_CALCULATION : Writing out the interpolated values to a file..........'

              !========================================
              ! Open for checking values from the interpolation
              !========================================
              OPEN (UNIT=unitInterp, FILE='interpOut.dat', &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
              IF (FileOpenStatus > 0) STOP "Cannot open  interpOut.dat file ********"


               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'O_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) O_density_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii


               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats

                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'O2_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) O2_density_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii


               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'N2_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) N2_density_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii

               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'Vx_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) Vx_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii


               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'Vy_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) Vy_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii

               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'wvz_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) wvz_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii

               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'tts_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) tts_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii


               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'qion3d_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) qion3d_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii

               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'elx_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) elx_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii

               do ii = 1, interface_hts
                  do jj = 1, high_lat_lats


                     WRITE (unitInterp, '(A,I2,A,I2,A)' ) 'ely_fixed_ht(',  ii, ',' ,jj, ',:) = '
                     WRITE (unitInterp, '(20ES20.10)' ) ely_fixed_ht(ii,jj,:)

                  end do  ! jj
               end do ! ii



               print *,'GIP_CALCULATION : STOPPING.............'
               STOP
           endif
















    ex2d(:,:) = elx_fixed_ht(10,:,:)
    ey2d(:,:) = ely_fixed_ht(10,:,:)

    CALL GLOBAL_IONOSPHERE_PLASMASPHERE ( &
         GIP_switches, &
         i_call_polar,i_call_plasma, &
         idump_GIP, &
         GIP_input_dataset, GIP_output_dataset, &
         sw_initialisation_call, &
         iday_number, &
         universal_time_seconds, &
         geo_grid_longitudes_degrees, geo_grid_latitudes_degrees, &
         ex2d,ey2d, &
         potential_field,ed1,ed2, &
         f107, &
         o_density_fixed_ht,o2_density_fixed_ht,n2_density_fixed_ht, &
         NO_density_fixed_ht, &
         N4S_density_fixed_ht,N2D_density_fixed_ht, &
         vx_fixed_ht,vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
         qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
         qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
         oplus_high_res_fixed,hplus_high_res_fixed, &
         noplus_high_res_fixed,o2plus_high_res_fixed, &
         n2plus_high_res_fixed,nplus_high_res_fixed, &
         Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed, &
         ne_high_res_fixed, &
         NmF2,HmF2_km,TEC,Te_400,Ti_400,this_Ne_profile, &
         dynamo_sigma_phph_dsi,dynamo_sigma_lmlm_msi, &
         dynamo_sigma_h,dynamo_sigma_c, &
         dynamo_Kdmph_dsi,dynamo_Kdmlm)

!     write(6,*) '2/3 done'

    call INTERFACE__GIP_to_thermosphere ( &
         thermospheric_model_name , therm_model_ht_dim , therm_model_lat_dim , therm_model_lon_dim, &
         ne_high_res_fixed,oplus_high_res_fixed,hplus_high_res_fixed, &
         noplus_high_res_fixed,o2plus_high_res_fixed, &
         n2plus_high_res_fixed,nplus_high_res_fixed, &
         Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed, &
         therm_model_geo_long,therm_model_geo_lat,therm_model_ht_m, &
         therm_model_Ne_density,therm_model_oplus_density,therm_model_hplus_density, &
         therm_model_noplus_density,therm_model_o2plus_density, &
         therm_model_n2plus_density,therm_model_nplus_density, &
         therm_model_Te,therm_model_Ti1,therm_model_Ti2)

!     write(6,*) 'All done'

!
    return



end SUBROUTINE gip_calculation      















!r========================================================
!r=            Global_Ionosphere_Plasmasphere            =
!r========================================================





SUBROUTINE GLOBAL_IONOSPHERE_PLASMASPHERE ( &
         GIP_switches, &
         i_call_polar,i_call_plasma, &
         idump, &
         GIP_input_dataset, GIP_output_dataset, &
         sw_initialisation_call, &
         iday_number, &
         universal_time_seconds, &
         geo_grid_longitudes_degrees, geo_grid_latitudes_degrees, &
         ex2d,ey2d, &
         potential_field,ed1,ed2, &
         f107, &
         o_density_fixed_ht,o2_density_fixed_ht,n2_density_fixed_ht, &
         NO_density_fixed_ht, &
         N4S_density_fixed_ht,N2D_density_fixed_ht, &
         vx_fixed_ht,vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
         qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
         qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
         Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
         NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
         N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
         Te_fixed_ht_com,Ti1_fixed_ht_com,Ti2_fixed_ht_com, &
         Ne_density_fixed_ht_com, &
         NmF2,HmF2_km,TEC,Te_400,Ti_400,this_Ne_profile, &
         dynamo_sigma_phph_dsi,dynamo_sigma_lmlm_msi, &
         dynamo_sigma_h,dynamo_sigma_c, &
         dynamo_Kdmph_dsi,dynamo_Kdmlm)




  IMPLICIT NONE

    INTEGER N_Pressure_Levels,N_Latitudes,N_longitudes
      PARAMETER(N_Pressure_Levels=15)
      PARAMETER(N_Latitudes=91)
      PARAMETER(N_Longitudes=20)

  LOGICAL :: sw_use_EUVAC_solar_spectrum
  LOGICAL :: sw_initialisation_call
  LOGICAL :: sw_reset_GIP_ion_densities_on_startup
  LOGICAL :: GIP_switches(20)
  character(100) :: GIP_input_dataset
  character(100) :: GIP_output_dataset
  INTEGER :: istop
  INTEGER :: iday_number 
  INTEGER :: i_total_no_days , i_graphics_out_start , i_no_day, &
             idump , ilat , ilon , idummy , i , iilat , iilon
  INTEGER :: ioutput_counter , iplasma_call_frequency_mins , &
             ioutput_high_res_counter , &
             ioutput_frequency_mins , &
             ioutput_high_res_frequency_mins
  INTEGER :: n , m , l
  INTEGER :: mp , lp
  INTEGER ::  file_res(27)
  INTEGER ::  i_call_polar,i_call_plasma
  INTEGER ::  iwrite_plasma_interface
  INTEGER ::  N_DYN_LAT
  INTEGER :: ihigh_lat_call_frequency_mins
  PARAMETER (N_DYN_LAT=97)
  REAL(kind=8) dynamo_sigma_phph_dsi(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_sigma_lmlm_msi(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_sigma_h(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_sigma_c(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_Kdmph_dsi(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_Kdmlm(NMP+1,N_DYN_LAT)
  REAL(kind=8) universal_time_seconds
  REAL(kind=8) :: potential_field(NMP+1,N_DYN_LAT)
  REAL(kind=8) :: ed1(NMP+1,N_DYN_LAT)
  REAL(kind=8) :: ed2(NMP+1,N_DYN_LAT)
  REAL(kind=8) Tn_plasma_input_3d(npts,nmp)
  REAL(kind=8) O_plasma_input_3d(npts,nmp)
  REAL(kind=8) O2_plasma_input_3d(npts,nmp)
  REAL(kind=8) N2_plasma_input_3d(npts,nmp)

  REAL(kind=8) NO_plasma_input_3d(npts,nmp)
  REAL(kind=8) N4S_plasma_input_3d(npts,nmp)
  REAL(kind=8) N2D_plasma_input_3d(npts,nmp)
  REAL(kind=8) Um_plasma_input_3d(npts,nmp)
  REAL(kind=8) Uz_plasma_input_3d(npts,nmp)
  REAL(kind=8) Uv_plasma_input_3d(npts,nmp)
  REAL(kind=8) Te_dum_plasma_input_3d(npts,nmp)
  REAL(kind=8) O_percent_failed , H_percent_failed
  REAL(kind=8) solar_declination_Angle_degrees
  REAL(kind=8) F107
  REAL(kind=8) ex2d(N_Latitudes,N_longitudes) , ey2d(N_Latitudes,N_longitudes)
  REAL(kind=8)  e1_not_used(2,nmp,nlp)
  REAL(kind=8)  e2_not_used(2,nmp,nlp)
  REAL(kind=8) geo_grid_longitudes_degrees(N_longitudes)
  REAL(kind=8) geo_grid_longitudes_radians(N_longitudes)
  REAL(kind=8) geo_grid_latitudes_degrees(N_Latitudes)
  REAL(kind=8) geo_grid_latitudes_radians(N_Latitudes)
  REAL(kind=8) geo_grid_colatitudes_degrees(N_Latitudes)
  REAL(kind=8) geo_grid_colatitudes_radians(N_Latitudes)
  REAL(kind=8) PI , DTR

  REAL(kind=8) :: O_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: O2_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N2_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: NO_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N4S_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N2D_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Vx_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Vy_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Wvz_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: tts_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qion3d_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qo2p_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qop_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qn2p_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qnp_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qtef_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: telec_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: tn_1d_fixed_ht(interface_hts) 
  REAL(kind=8) :: te_1d_fixed_ht(interface_hts)
  REAL(kind=8) :: zkm_fixed_ht(interface_hts)

  REAL(kind=8)  btotal_2000(91,20)
  REAL(kind=8)  dip_2000(91,20)
  REAL(kind=8)  decl_2000(91,20)
  REAL(kind=8)  btot_nanoTesla(91,20)
  REAL(kind=8)  btotal_Tesla(91,20)
  REAL(kind=8)  dip_angle_degrees(91,20)
  REAL(kind=8)  dipav1 , dipav91
  REAL(kind=8)  bav1 , bav91
  REAL(kind=8)  declination_degrees(91,20)
  REAL(kind=8)  ty
  REAL(kind=8)  geo_local_time_degrees(N_longitudes)
  REAL(kind=8)  fye
  REAL(kind=8)  solar_zenith_angle_radians , rlt , rlat

  REAL(kind=8) :: Oplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Hplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: NOplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: O2plus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: N2plus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Nplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Te_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Ti1_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Ti2_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Ne_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Oplus_density_at_300(90,91)

  REAL(kind=8) :: Oplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Hplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: NOplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: O2plus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: N2plus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Nplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Te_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Ti1_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Ti2_fixed_ht_plasma(interface_hts,91,90)

  REAL(kind=8) :: Oplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Hplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: NOplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: O2plus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: N2plus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Nplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Te_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Ti1_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Ti2_fixed_ht_polar(interface_hts,91,20)

  REAL(kind=8) :: NmF2(91,90)
  REAL(kind=8) :: HmF2_km(91,90)
  REAL(kind=8) :: TEC(91,90)
  REAL(kind=8) :: Te_400(91,90)
  REAL(kind=8) :: Ti_400(91,90)
  REAL(kind=8) :: this_Ne_profile(interface_hts)

  REAL(kind=8) :: high_lat_time_step_seconds
  REAL(kind=8) :: plasma_time_step_seconds

  PARAMETER (PI=3.14159,DTR=PI/180.0)

  ! print out results from INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE ?
  LOGICAL :: checkFixedGeo = .FALSE.
  INTEGER :: FileOpenStatus = 0
  INTEGER, parameter :: unitFixedGeo = 15

! BEGIN CODE =================================================


! set a few switches here to not bother the outside folks...
!
    sw_use_EUVAC_solar_spectrum = .FALSE.
    ioutput_counter = 0
    ioutput_high_res_counter = 0
    ioutput_frequency_mins = 60
    ioutput_high_res_frequency_mins = 15
    i_no_day = 0
    i_graphics_out_start = 1
    file_res(:) = 0
    iplasma_call_frequency_mins = 15
    ihigh_lat_call_frequency_mins = 5


  do l = 1 , N_longitudes
    geo_grid_longitudes_degrees(l) = (float(l-1))*18.
  enddo

  do m = 1 , N_latitudes
    geo_grid_latitudes_degrees(m) = (float(m - 46))*2.
  enddo

    do l = 1 , N_longitudes
      geo_grid_longitudes_radians(l) = geo_grid_longitudes_degrees(l) * DTR
    enddo
    do m = 1 , N_Latitudes
      geo_grid_latitudes_radians(m) = geo_grid_latitudes_degrees(m) * DTR
      geo_grid_colatitudes_degrees(m) = 90. - geo_grid_latitudes_degrees(m)
      geo_grid_colatitudes_radians(m) = geo_grid_colatitudes_degrees(m) * DTR
    enddo

    ty = (iday_number+15.5)*12./365.
    IF ( ty.GT.12.0 ) ty = ty - 12.0
    solar_declination_angle_degrees = (ATAN(0.434*SIN(PI/6.0*(ty-3.17)))) / DTR

  M_plasma(1) = 16.
  M_plasma(2) = 1.
  M_plasma(3) = 4.
  M_plasma(4) = 28.
  M_plasma(5) = 32.
  M_plasma(6) = 30.

! boltzmanns constant (kg).....

  DO j = 1 , 6
      KM_plasma(j) = 1.381E-23/(M_plasma(j)*1.673E-27)
  ENDDO

!g
!g At initialisation we need to read in the '.gip' file......
!g
  if( sw_initialisation_call ) then

      WRITE(6,*) '**************** GIP INITIALISATION CALL ******************* '

      i_first_call_of_plasma = 1
      write(6,*) '******** reading GIP startup file ********'

!     sw_reset_GIP_ion_densities_on_startup = .FALSE.
      sw_reset_GIP_ion_densities_on_startup = .TRUE.
      IF(sw_reset_GIP_ion_densities_on_startup) then

         d13d(:,:,:,:) = 1.e9
         d23d(:,:,:,:) = 1.e9
         v13d(:,:,:,:) = 0.0
         v23d(:,:,:,:) = 0.0
         vpeq(:,:) = 0.0
         vzon(:,:) = 0.0
         ni_plasma(:,:,:) = 1.e09
         vi_plasma(:,:,:) = 0.
         ti_plasma(:,:,:) = 1000.
         te_plasma(:,:) = 2500.
         ne_plasma(:,:) = 2.e09
         no_plus_3d(:,:) = 1.e3
         o2_plus_3d(:,:) = 1.e3
         Oplus_density_fixed_ht_com(:,:,:) = 1.e09
         Hplus_density_fixed_ht_com(:,:,:) = 1.e09
         NOplus_density_fixed_ht_com(:,:,:) = 1.e03
         O2plus_density_fixed_ht_com(:,:,:) = 1.e03
         N2plus_density_fixed_ht_com(:,:,:) = 1.e03
         Nplus_density_fixed_ht_com(:,:,:) = 1.e03
         Te_fixed_ht_com(:,:,:) = 2500.
         Ti1_fixed_ht_com(:,:,:) = 1000.
         Ti2_fixed_ht_com(:,:,:) = 1000.
         Ne_density_fixed_ht_com(:,:,:) = 2.e09

         write(6,*) '********** Using COLD STARTUP for GIP **********'

      else

         call IO__read_gip_netcdf_history(GIP_input_dataset, &
                                      ni_plasma, vi_plasma,ne_plasma, ti_plasma,te_plasma, &
                                      no_plus_3d, o2_plus_3d, d13d, d23d, v13d, v23d, &
                                      vpeq, vzon, &
                                      Oplus_density_fixed_ht_com, Hplus_density_fixed_ht_com, &
                                      NOplus_density_fixed_ht_com, O2plus_density_fixed_ht_com, &
                                      N2plus_density_fixed_ht_com, Nplus_density_fixed_ht_com, &
                                      Te_fixed_ht_com, Ti1_fixed_ht_com, Ti2_fixed_ht_com, &
                                      Ne_density_fixed_ht_com)

      endif

      write(6,*) '******** finished reading GIP startup file ********'

     !OPEN (22,FILE='../static_files/hprof',STATUS='old')
!
! Magnetic field stuff from Tucan....
!
      !OPEN (174,FILE='../static_files/angdif_for_GIP',STATUS='old')
           read(174,*) angdif
       close (174)

      !OPEN (175,FILE='../static_files/btotal_dip_declination_old_dipole',STATUS='old')
       do ilat = 1 , 91
       do ilon = 1 , 20
         read(175,*) idummy, idummy, btotal_2000(ilat,ilon), dip_2000(ilat,ilon), decl_2000(ilat,ilon)
       enddo
       enddo
       close (175)
!g
       DO ilon = 1 , 20
          DO ilat = 2 , 90
             btot_nanoTesla(ilat,ilon) = btotal_2000(ilat,ilon)
             dip_angle_degrees(ilat,ilon) = dip_2000(ilat,ilon)
             declination_degrees(ilat,ilon) = decl_2000(ilat,ilon)
             IF ( ABS(dip_angle_degrees(ilat,ilon)) .LT.5. ) THEN
             IF ( dip_angle_degrees(ilat,ilon).LT.0.0 ) dip_angle_degrees(ilat,ilon) = -5.
             IF ( dip_angle_degrees(ilat,ilon).GT.0.0 ) dip_angle_degrees(ilat,ilon) = 5.
             ENDIF
          ENDDO   
       ENDDO   

!g
!g  define btot and dip_angle_degrees over the poles....
!g
       bav1 = 0.0
       dipav1 = 0.0
       bav91 = 0.0
       dipav91 = 0.0
       DO ilon = 1 , 20
          bav1 = bav1 + btot_nanoTesla(2,ilon)/20.
          dipav1 = dipav1 + dip_angle_degrees(2,ilon)/20.
          bav91 = bav91 + btot_nanoTesla(90,ilon)/20.
          dipav91 = dipav91 + dip_angle_degrees(90,ilon)/20.
       ENDDO    
       DO ilon = 1 , 20
          dip_angle_degrees(1,ilon) = dipav1
          btot_nanoTesla(1,ilon) = bav1
          dip_angle_degrees(91,ilon) = dipav91
          btot_nanoTesla(91,ilon) = bav91
       ENDDO   
       DO ilon = 1 , 20
          DO ilat = 1 , 91
            btotal_Tesla(ilat,ilon) = btot_nanoTesla(ilat,ilon)*1.E-09
            dip_angle_radians(ilat,ilon) = dip_angle_degrees(ilat,ilon)*DTR
            btotal(ilat,ilon) = btotal_Tesla(ilat,ilon)
            declination(ilat,ilon) = declination_degrees(ilat,ilon)
          ENDDO
       ENDDO
!
! end of the mag field stuff from tucan
!



      sw_1st_call_int_fixed_ht = .TRUE.
      UT_previous_call_plasma_secs = universal_time_seconds
      UT_previous_call_polar_secs = universal_time_seconds
  endif
      if ( .NOT. sw_initialisation_call ) then

      !WRITE(6,*) '**************** GIP NORMAL CALL ******************* '

       DO l = 1 , 20
        geo_local_time_degrees(l) = geo_grid_longitudes_degrees(l) + (universal_time_seconds - 43200.)/240.0
        IF ( geo_local_time_degrees(l).GE.360.0 ) geo_local_time_degrees(l) = geo_local_time_degrees(l) - 360.0
        IF ( geo_local_time_degrees(l).lt.0.0 ) geo_local_time_degrees(l) = geo_local_time_degrees(l) + 360.0
        fye = geo_local_time_degrees(l)*dtr
        rlt = 180.0 + geo_local_time_degrees(l)
        IF ( rlt.GT.360.0 ) rlt = rlt - 360.0
        rlt = rlt*DTR
          DO m = 1 , 91
            rlat = geo_grid_latitudes_radians(m)
            solar_zenith_angle_radians =  &
                ACOS(-COS(rlat)*COS(solar_declination_angle_degrees*DTR)*COS(rlt)+SIN(rlat) &
                *SIN(solar_declination_angle_degrees*DTR))

            do n = 1 , interface_hts
               tn_1d_fixed_ht(n) = tts_fixed_ht(n,m,l)
!              zkm_fixed_ht(n) = (5.* (n-1)) + 90.
               zkm_fixed_ht(n) = fixed_heights_km(n)
            enddo

            do n = 1 , interface_hts
              telec_fixed_ht(n,m,l) = te_1d_fixed_ht(n)
            enddo
          ENDDO
       ENDDO

!  sort out the hygrogen density at 400km which gets fed into both Polar and Plasma codes....

!   write(6,*) '***************** HYDROGEN GRID AT 400KM ***********************************'
!   write(6,1888) Hyd_grid_400km_m3 / 1.e11
!   1888 format(25f4.1)
!   write(6,*) '****************************************************************************'


       endif  ! .NOT. sw_initialisation IF 


  if (i_call_polar == 1) then

      high_lat_time_step_seconds = universal_time_seconds - UT_previous_call_polar_secs
      if ( high_lat_time_step_seconds < 0.0 ) then
         high_lat_time_step_seconds = high_lat_time_step_seconds + 86400.
      endif


      CALL HL__POLAR_IONOSPHERE ( &
      GIP_switches, &
      geo_grid_longitudes_degrees,geo_grid_longitudes_radians, & 
      geo_grid_latitudes_degrees,geo_grid_latitudes_radians, &
      geo_grid_colatitudes_degrees,geo_grid_colatitudes_radians, &
      high_lat_time_step_seconds, &
      angdif, &
      btotal,dip_angle_radians,f107,hprof,dth_radians, &
      ex2d,ey2d,v13d,v23d,tarea,d13d,d23d, &
      gravity,nhgt,hz_metres,hnew,h2,h2diff,hprod2,hprod3,zkm,l300,l1000,r0,declination, &
      Solar_Declination_Angle_degrees,km_plasma,universal_time_seconds, &
      ioutput_counter,ioutput_high_res_counter, &
      ioutput_frequency_mins,ioutput_high_res_frequency_mins, &
      ihigh_lat_call_frequency_mins,file_res, &
      sw_initialisation_call, &
      O_density_fixed_ht,O2_density_fixed_ht,N2_density_fixed_ht,&
      NO_density_fixed_ht,&
      N4S_density_fixed_ht,N2D_density_fixed_ht,&
      Vx_fixed_ht,Vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
      qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
      qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
      factor_ht_3d,iht_above_3d,iht_below_3d, &
      factor_ht_inverse_3d,iht_above_inverse_3d,iht_below_inverse_3d, &
      altitude_metres_3d, &
      Oplus_density_fixed_ht_polar,Hplus_density_fixed_ht_polar,NOplus_density_fixed_ht_polar, &
      O2plus_density_fixed_ht_polar,N2plus_density_fixed_ht_polar,Nplus_density_fixed_ht_polar, &
      Te_fixed_ht_polar,Ti1_fixed_ht_polar,Ti2_fixed_ht_polar,iday_number)

      UT_previous_call_polar_secs = universal_time_seconds
  endif
!g
!g
  if(i_call_plasma == 1) then

      if ( .NOT. sw_initialisation_call ) then
!       write(6,*) '*********************** HERE 3 *******************',sw_initialisation_call
          istop = 0
          if (istop == 1) stop
          iwrite_plasma_interface = 0

  call INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE( &
        geo_grid_longitudes_degrees, &
        geo_grid_latitudes_degrees, &
        O_density_fixed_ht, O2_density_fixed_ht, N2_density_fixed_ht, &
        NO_density_fixed_ht, &
        N4S_density_fixed_ht, N2D_density_fixed_ht, &
        VX_fixed_ht, VY_fixed_ht, WVZ_fixed_ht, TTS_fixed_ht,telec_fixed_ht, &
        IN, IS, iwrite_plasma_interface, TN_plasma_input_3d, &
        O_plasma_input_3d, &
        O2_plasma_input_3d, N2_plasma_input_3d, &
        NO_plasma_input_3d, N4S_plasma_input_3d, N2D_plasma_input_3d, &
        GLAt_plasma_3d, &
        GLOnd_plasma_3d, &
        PZ_plasma_3d, &
        um_plasma_input_3d,uz_plasma_input_3d,uv_plasma_input_3d, &
        te_dum_plasma_input_3d, &
        ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &
        ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht, ihu_3d_fixed_ht, &
        sw_1st_call_int_fixed_ht, &
        GIP_switches)

checkFixedGeo = .TRUE.
!---------------------------------------
! print out results to check vs GT-IPE
!---------------------------------------
If (checkFixedGeo) then

    print *,'GLOBAL_IONOSPHERE_PLASMASPHERE : AFTER CALLING INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE.............'
    !========================================
    ! Open for checking values from GT 
    !========================================
    OPEN (UNIT=unitFixedGeo, FILE='./data/CheckFixedGeo.dat', &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
    IF (FileOpenStatus > 0) STOP "Cannot open  ./data/CheckFixedGeo.dat file ********"

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(1,:) = ', TN_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(700,:) = ', TN_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(NPTS,:) = ', TN_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'O_plasma_input_3d(1,:) = ', O_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'O_plasma_input_3d(700,:) = ', O_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'O_plasma_input_3d(NPTS,:) = ', O_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(1,:) = ', O2_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(700,:) = ', O2_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(NPTS,:) = ', O2_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(1,:) = ', N2_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(700,:) = ', N2_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(NPTS,:) = ', N2_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(1,:) = ', GLAt_plasma_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(700,:) = ', GLAt_plasma_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(NPTS,:) = ', GLAt_plasma_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(1,:) = ', GLOnd_plasma_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(700,:) = ', GLOnd_plasma_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(NPTS,:) = ', GLOnd_plasma_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'pz_plasma_3d(1,:) = ', pz_plasma_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'pz_plasma_3d(700,:) = ', pz_plasma_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'pz_plasma_3d(NPTS,:) = ', pz_plasma_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'um_plasma_input_3d(1,:) = ', um_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'um_plasma_input_3d(700,:) = ', um_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'um_plasma_input_3d(NPTS,:) = ', um_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'uz_plasma_input_3d(1,:) = ', uz_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'uz_plasma_input_3d(700,:) = ', uz_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'uz_plasma_input_3d(NPTS,:) = ', uz_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'uv_plasma_input_3d(1,:) = ', uv_plasma_input_3d(1,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'uv_plasma_input_3d(700,:) = ', uv_plasma_input_3d(700,:)
    WRITE (unitFixedGeo, '(A, 80ES20.10)' ) 'uv_plasma_input_3d(NPTS,:) = ', uv_plasma_input_3d(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(1,:) = ', ilon1_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(700,:) = ', ilon1_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(NPTS,:) = ', ilon1_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(1,:) = ', ilon2_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(700,:) = ', ilon2_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(NPTS,:) = ', ilon2_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(1,:) = ', ilat1_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(700,:) = ', ilat1_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(NPTS,:) = ', ilat1_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(1,:) = ', ilat2_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(700,:) = ', ilat2_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(NPTS,:) = ', ilat2_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(1,:) = ', ispecial_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(700,:) = ', ispecial_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(NPTS,:) = ', ispecial_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ihl_3d_fixed_ht(1,:) = ', ihl_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ihl_3d_fixed_ht(700,:) = ', ihl_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ihl_3d_fixed_ht(NPTS,:) = ', ihl_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ihu_3d_fixed_ht(1,:) = ', ihu_3d_fixed_ht(1,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ihu_3d_fixed_ht(700,:) = ', ihu_3d_fixed_ht(700,:)
    WRITE (unitFixedGeo, '(A, 80I7)' ) 'ihu_3d_fixed_ht(NPTS,:) = ', ihu_3d_fixed_ht(NPTS,:)

    WRITE (unitFixedGeo, '(A)' ) ' '
    WRITE (unitFixedGeo, '(A, I3)' ) 'sw_1st_call_int_fixed_ht = ', sw_1st_call_int_fixed_ht




    print *,'STOPPING AFTER CHECKING FIXED GEO ******************************************'
    STOP

endif ! checkFixedGeo


      plasma_time_step_seconds = universal_time_seconds - UT_previous_call_plasma_secs
      if ( plasma_time_step_seconds < 0.0 ) then
         plasma_time_step_seconds = plasma_time_step_seconds + 86400.
      endif

      endif

      CALL ML__MID_AND_LOW_LATITUDE_IONOSPHERE ( &
      GIP_switches, &
          geo_grid_longitudes_degrees,geo_grid_longitudes_radians, &
          geo_grid_latitudes_degrees,geo_grid_latitudes_radians, &
          geo_grid_colatitudes_degrees,geo_grid_colatitudes_radians, &
      iday_number, &
      i_first_call_of_plasma,universal_time_seconds,plasma_time_step_seconds, &
      idump,ioutput_counter,ioutput_high_res_counter,iplasma_call_frequency_mins, &
      ioutput_frequency_mins,ioutput_high_res_frequency_mins, &
      i_no_day,i_graphics_out_start, &
      file_res, &
      sw_initialisation_call, &
      F107,Solar_Declination_Angle_degrees, &
      ETA_Apex_3D, &
      Apex_D1,Apex_D2, &
      Apex_E1,Apex_E2, &
      Apex_BE3,Apex_bhat,Apex_Bmag, &
      Apex_D,Apex_d1d1,Apex_d1d2,Apex_d2d2,Apex_grdlbm2, &
      bcol_plasma,blon_plasma,q_coordinate_plasma,Re_apex_plasma, &
      gr_plasma,gcol_plasma,glon_plasma, &
      IN,IS, &
      glat_plasma_3d,glond_plasma_3d,pz_plasma_3d, &
      TN_plasma_input_3d,O_plasma_input_3d,O2_plasma_input_3d,N2_plasma_input_3d, &
      NO_plasma_input_3d,N4S_plasma_input_3d,N2D_plasma_input_3d, &
      te_dum_plasma_input_3d, &
      um_plasma_input_3d,uz_plasma_input_3d,uv_plasma_input_3d, &
      TI_plasma,TE_plasma,NI_plasma,VI_plasma,NE_plasma,YNI,YVI,YNE,YQE,no_plus_3d,o2_plus_3d, &
      n2_plus_3d,n_plus_3d, &
      JIOn,mass,ndt,nions,midpoint,itday,qb,KM_plasma,M_plasma, &
      dynamo_sigma_phph_dsi,dynamo_sigma_lmlm_msi, &
      dynamo_sigma_h,dynamo_sigma_c, &
      dynamo_Kdmph_dsi,dynamo_Kdmlm, &
      potential_field,ed1,ed2, &
      vpeq,vzon)


      i_first_call_of_plasma = 0
        !write(6,*) '*********************** HERE 5 *******************',sw_initialisation_call


      CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
      ni_plasma, &
      no_plus_3d,o2_plus_3d,n2_plus_3d,n_plus_3d, &
      TE_plasma,TI_plasma, &
      Oplus_density_fixed_ht_plasma,Hplus_density_fixed_ht_plasma, &
      Noplus_density_fixed_ht_plasma,O2plus_density_fixed_ht_plasma, &
      N2plus_density_fixed_ht_plasma,Nplus_density_fixed_ht_plasma, &
      Te_fixed_ht_plasma,Ti1_fixed_ht_plasma,Ti2_fixed_ht_plasma)

  !g

      UT_previous_call_plasma_secs = universal_time_seconds
  endif

      if ( .NOT. sw_initialisation_call ) then

          call INTERFACE__COMBINE_MID_LAT_AND_POLAR_PARAMS ( &
                     Oplus_density_fixed_ht_plasma,Hplus_density_fixed_ht_plasma, &
                     NOplus_density_fixed_ht_plasma,O2plus_density_fixed_ht_plasma, &
                     N2plus_density_fixed_ht_plasma,Nplus_density_fixed_ht_plasma, &
                     Te_fixed_ht_plasma,Ti1_fixed_ht_plasma,Ti2_fixed_ht_plasma, &
                     Oplus_density_fixed_ht_polar,Hplus_density_fixed_ht_polar, &
                     NOplus_density_fixed_ht_polar,O2plus_density_fixed_ht_polar, &
                     N2plus_density_fixed_ht_polar,Nplus_density_fixed_ht_polar, &
                     Te_fixed_ht_polar,Ti1_fixed_ht_polar,Ti2_fixed_ht_polar, &
                     Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
                     NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
                     N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
                     Te_fixed_ht_com,Ti1_fixed_ht_com,Ti2_fixed_ht_com)


          call INTERFACE__calculate_Ne_from_ion_densities ( &
                     Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
                     NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
                     N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
                     Ne_density_fixed_ht_com)

      endif


!g
!g If we are at the end of the run then we need to write out the '.gip' file......
!g
  if(idump == 1) then
      write(6,*) '******** writing GIP startup file ********'


        call IO__write_gip_netcdf_history(GIP_output_dataset, &
                                      ni_plasma,vi_plasma,ne_plasma,ti_plasma,te_plasma, &
                                      no_plus_3d,o2_plus_3d,d13d,d23d,v13d,v23d, &
                                      vpeq,vzon, &
                                      Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
                                      NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
                                      N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
                                      Te_fixed_ht_com,Ti1_fixed_ht_com,Ti2_fixed_ht_com, &
                                      Ne_density_fixed_ht_com)

      write(6,*) '******** finished writing GIP startup file ********'
  endif

  return




end SUBROUTINE GLOBAL_IONOSPHERE_PLASMASPHERE






SUBROUTINE INTERFACE__thermosphere_to_GIP ( &
           GIP_switches, &
           thermospheric_model_name,ht_dim,lat_dim,lon_dim, &
           therm_o_density,therm_o2_density,therm_n2_density, &
           therm_NO_density, &
           therm_N4S_density,therm_N2D_density, &
           therm_Tn,therm_Un,therm_Vn, &
           therm_qion3d, therm_elx, therm_ely, &
           therm_qo2p_aurora, therm_qop_aurora, therm_qn2p_aurora, &
           therm_qnp_aurora, therm_qtef_aurora, &
           therm_long,therm_lat,therm_Z, &
           o_density_fixed_ht,o2_density_fixed_ht,n2_density_fixed_ht, &
           NO_density_fixed_ht, &
           N4S_density_fixed_ht,N2D_density_fixed_ht, &
           Vx_fixed_ht,Vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht, &
           qion3d_fixed_ht, elx_fixed_ht, ely_fixed_ht, &
           qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
           qnp_aurora_fixed_ht, qtef_aurora_fixed_ht)

IMPLICIT NONE

character*10 thermospheric_model_name
integer ht_dim , lat_dim , lon_dim
INTEGER ilon  , ilat , iht
INTEGER ilon_therm , ilat_therm , iht_therm 
INTEGER ispecial 
LOGICAL :: GIP_switches(20)
LOGICAL :: sw_External_model_provides_NO_N4S_densities
LOGICAL :: sw_input_Auroral_production_is_single_overall_rate
REAL(kind=8) therm_long(lon_dim)
REAL(kind=8) therm_lat(lat_dim)
REAL(kind=8) therm_Z(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_height(ht_dim)
REAL(kind=8) therm_o_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_o2_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_n2_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_NO_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_N4S_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_N2D_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_Tn(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_Un(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_Vn(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_qion3d(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_elx(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_ely(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_qo2p_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_qop_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_qn2p_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_qnp_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) therm_qtef_aurora(ht_dim,lon_dim,lat_dim)

!                               ( hts, lats, lons )
REAL(kind=8) interface_o_density(interface_hts,91,20)
REAL(kind=8) interface_o2_density(interface_hts,91,20)
REAL(kind=8) interface_n2_density(interface_hts,91,20)
REAL(kind=8) interface_NO_density(interface_hts,91,20)
REAL(kind=8) interface_N4S_density(interface_hts,91,20)
REAL(kind=8) interface_N2D_density(interface_hts,91,20)
REAL(kind=8) interface_Tn(interface_hts,91,20)
REAL(kind=8) interface_Un(interface_hts,91,20)
REAL(kind=8) interface_Vn(interface_hts,91,20)
REAL(kind=8) interface_qion3d(interface_hts,91,20)
REAL(kind=8) interface_elx(interface_hts,91,20)
REAL(kind=8) interface_ely(interface_hts,91,20)
REAL(kind=8) interface_qo2p_aurora(interface_hts,91,20)
REAL(kind=8) interface_qop_aurora(interface_hts,91,20)
REAL(kind=8) interface_qn2p_aurora(interface_hts,91,20)
REAL(kind=8) interface_qnp_aurora(interface_hts,91,20)
REAL(kind=8) interface_qtef_aurora(interface_hts,91,20)


REAL(kind=8) o_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) o2_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) n2_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) NO_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) N4S_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) N2D_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) Vx_fixed_ht(interface_hts,91,20)
REAL(kind=8) Vy_fixed_ht(interface_hts,91,20)
REAL(kind=8) wvz_fixed_ht(interface_hts,91,20)
REAL(kind=8) tts_fixed_ht(interface_hts,91,20)
REAL(kind=8) qion3d_fixed_ht(interface_hts,91,20)
REAL(kind=8) elx_fixed_ht(interface_hts,91,20)
REAL(kind=8) ely_fixed_ht(interface_hts,91,20)
REAL(kind=8) qo2p_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qop_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qn2p_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qnp_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qtef_aurora_fixed_ht(interface_hts,91,20)


INTEGER ilon_west_array(20)
INTEGER ilon_east_array(20)
INTEGER ilat_north_array(91)
INTEGER ilat_south_array(91)
REAL(kind=8) factor_lat_array(91)
REAL(kind=8) factor_lon_array(20)
INTEGER ilon_west
INTEGER ilon_east
INTEGER ilat_north
INTEGER ilat_south
REAL(kind=8) factor_lat
REAL(kind=8) factor_lon

INTEGER iht_above_west_north(interface_hts,91,20)
INTEGER iht_below_west_north(interface_hts,91,20)
INTEGER iht_above_east_north(interface_hts,91,20)
INTEGER iht_below_east_north(interface_hts,91,20)
INTEGER iht_above_west_south(interface_hts,91,20)
INTEGER iht_below_west_south(interface_hts,91,20)
INTEGER iht_above_east_south(interface_hts,91,20)
INTEGER iht_below_east_south(interface_hts,91,20)
INTEGER iht_above
INTEGER iht_below

REAL(kind=8) factor_ht_west_north(interface_hts,91,20)
REAL(kind=8) factor_ht_east_north(interface_hts,91,20)
REAL(kind=8) factor_ht_west_south(interface_hts,91,20)
REAL(kind=8) factor_ht_east_south(interface_hts,91,20)

REAL(kind=8) factor_ht12
REAL(kind=8) factor_ht22
REAL(kind=8) factor_ht11
REAL(kind=8) factor_ht21

REAL(kind=8) O_den_p_u12_log
REAL(kind=8) O_den_p_l12_log
REAL(kind=8) O_den_p_12_log
REAL(kind=8) O2_den_p_u12_log
REAL(kind=8) O2_den_p_l12_log
REAL(kind=8) O2_den_p_12_log
REAL(kind=8) N2_den_p_u12_log
REAL(kind=8) N2_den_p_l12_log
REAL(kind=8) N2_den_p_12_log
REAL(kind=8) HYD_den_p_u12_log
REAL(kind=8) HYD_den_p_l12_log
REAL(kind=8) HYD_den_p_12_log
REAL(kind=8) HEL_den_p_u12_log
REAL(kind=8) HEL_den_p_l12_log
REAL(kind=8) HEL_den_p_12_log
REAL(kind=8) NO_den_p_u12_log
REAL(kind=8) NO_den_p_l12_log
REAL(kind=8) NO_den_p_12_log
REAL(kind=8) N4S_den_p_u12_log
REAL(kind=8) N4S_den_p_l12_log
REAL(kind=8) N4S_den_p_12_log
REAL(kind=8) N2D_den_p_u12_log
REAL(kind=8) N2D_den_p_l12_log
REAL(kind=8) N2D_den_p_12_log
REAL(kind=8) Tn_p_u12
REAL(kind=8) Tn_p_l12
REAL(kind=8) Tn_p_12
REAL(kind=8) Un_p_u12
REAL(kind=8) Un_p_l12
REAL(kind=8) Un_p_12
REAL(kind=8) Vn_p_u12
REAL(kind=8) Vn_p_l12
REAL(kind=8) Vn_p_12
REAL(kind=8) qion3d_p_u12
REAL(kind=8) qion3d_p_l12
REAL(kind=8) qion3d_p_12
REAL(kind=8) elx_p_u12
REAL(kind=8) elx_p_l12
REAL(kind=8) elx_p_12
REAL(kind=8) ely_p_u12
REAL(kind=8) ely_p_l12
REAL(kind=8) ely_p_12
REAL(kind=8) qo2p_aurora_p_u12
REAL(kind=8) qo2p_aurora_p_l12
REAL(kind=8) qo2p_aurora_p_12
REAL(kind=8) qop_aurora_p_u12
REAL(kind=8) qop_aurora_p_l12
REAL(kind=8) qop_aurora_p_12
REAL(kind=8) qn2p_aurora_p_u12
REAL(kind=8) qn2p_aurora_p_l12
REAL(kind=8) qn2p_aurora_p_12
REAL(kind=8) qnp_aurora_p_u12
REAL(kind=8) qnp_aurora_p_l12
REAL(kind=8) qnp_aurora_p_12
REAL(kind=8) qtef_aurora_p_u12
REAL(kind=8) qtef_aurora_p_l12
REAL(kind=8) qtef_aurora_p_12


REAL(kind=8) O_den_p_u22_log
REAL(kind=8) O_den_p_l22_log
REAL(kind=8) O_den_p_22_log
REAL(kind=8) O2_den_p_u22_log
REAL(kind=8) O2_den_p_l22_log
REAL(kind=8) O2_den_p_22_log
REAL(kind=8) N2_den_p_u22_log
REAL(kind=8) N2_den_p_l22_log
REAL(kind=8) N2_den_p_22_log
REAL(kind=8) HYD_den_p_u22_log
REAL(kind=8) HYD_den_p_l22_log
REAL(kind=8) HYD_den_p_22_log
REAL(kind=8) HEL_den_p_u22_log
REAL(kind=8) HEL_den_p_l22_log
REAL(kind=8) HEL_den_p_22_log
REAL(kind=8) NO_den_p_u22_log
REAL(kind=8) NO_den_p_l22_log
REAL(kind=8) NO_den_p_22_log
REAL(kind=8) N4S_den_p_u22_log
REAL(kind=8) N4S_den_p_l22_log
REAL(kind=8) N4S_den_p_22_log
REAL(kind=8) N2D_den_p_u22_log
REAL(kind=8) N2D_den_p_l22_log
REAL(kind=8) N2D_den_p_22_log
REAL(kind=8) Tn_p_u22
REAL(kind=8) Tn_p_l22
REAL(kind=8) Tn_p_22
REAL(kind=8) Un_p_u22
REAL(kind=8) Un_p_l22
REAL(kind=8) Un_p_22
REAL(kind=8) Vn_p_u22
REAL(kind=8) Vn_p_l22
REAL(kind=8) Vn_p_22
REAL(kind=8) qion3d_p_u22
REAL(kind=8) qion3d_p_l22
REAL(kind=8) qion3d_p_22
REAL(kind=8) elx_p_u22
REAL(kind=8) elx_p_l22
REAL(kind=8) elx_p_22
REAL(kind=8) ely_p_u22
REAL(kind=8) ely_p_l22
REAL(kind=8) ely_p_22
REAL(kind=8) qo2p_aurora_p_u22
REAL(kind=8) qo2p_aurora_p_l22
REAL(kind=8) qo2p_aurora_p_22
REAL(kind=8) qop_aurora_p_u22
REAL(kind=8) qop_aurora_p_l22
REAL(kind=8) qop_aurora_p_22
REAL(kind=8) qn2p_aurora_p_u22
REAL(kind=8) qn2p_aurora_p_l22
REAL(kind=8) qn2p_aurora_p_22
REAL(kind=8) qnp_aurora_p_u22
REAL(kind=8) qnp_aurora_p_l22
REAL(kind=8) qnp_aurora_p_22
REAL(kind=8) qtef_aurora_p_u22
REAL(kind=8) qtef_aurora_p_l22
REAL(kind=8) qtef_aurora_p_22


REAL(kind=8) O_den_p_u11_log
REAL(kind=8) O_den_p_l11_log
REAL(kind=8) O_den_p_11_log
REAL(kind=8) O2_den_p_u11_log
REAL(kind=8) O2_den_p_l11_log
REAL(kind=8) O2_den_p_11_log
REAL(kind=8) N2_den_p_u11_log
REAL(kind=8) N2_den_p_l11_log
REAL(kind=8) N2_den_p_11_log
REAL(kind=8) HYD_den_p_u11_log
REAL(kind=8) HYD_den_p_l11_log
REAL(kind=8) HYD_den_p_11_log
REAL(kind=8) HEL_den_p_u11_log
REAL(kind=8) HEL_den_p_l11_log
REAL(kind=8) HEL_den_p_11_log
REAL(kind=8) NO_den_p_u11_log
REAL(kind=8) NO_den_p_l11_log
REAL(kind=8) NO_den_p_11_log
REAL(kind=8) N4S_den_p_u11_log
REAL(kind=8) N4S_den_p_l11_log
REAL(kind=8) N4S_den_p_11_log
REAL(kind=8) N2D_den_p_u11_log
REAL(kind=8) N2D_den_p_l11_log
REAL(kind=8) N2D_den_p_11_log
REAL(kind=8) Tn_p_u11
REAL(kind=8) Tn_p_l11
REAL(kind=8) Tn_p_11
REAL(kind=8) Un_p_u11
REAL(kind=8) Un_p_l11
REAL(kind=8) Un_p_11
REAL(kind=8) Vn_p_u11
REAL(kind=8) Vn_p_l11
REAL(kind=8) Vn_p_11
REAL(kind=8) qion3d_p_u11
REAL(kind=8) qion3d_p_l11
REAL(kind=8) qion3d_p_11
REAL(kind=8) elx_p_u11
REAL(kind=8) elx_p_l11
REAL(kind=8) elx_p_11
REAL(kind=8) ely_p_u11
REAL(kind=8) ely_p_l11
REAL(kind=8) ely_p_11
REAL(kind=8) qo2p_aurora_p_u11
REAL(kind=8) qo2p_aurora_p_l11
REAL(kind=8) qo2p_aurora_p_11
REAL(kind=8) qop_aurora_p_u11
REAL(kind=8) qop_aurora_p_l11
REAL(kind=8) qop_aurora_p_11
REAL(kind=8) qn2p_aurora_p_u11
REAL(kind=8) qn2p_aurora_p_l11
REAL(kind=8) qn2p_aurora_p_11
REAL(kind=8) qnp_aurora_p_u11
REAL(kind=8) qnp_aurora_p_l11
REAL(kind=8) qnp_aurora_p_11
REAL(kind=8) qtef_aurora_p_u11
REAL(kind=8) qtef_aurora_p_l11
REAL(kind=8) qtef_aurora_p_11


REAL(kind=8) O_den_p_u21_log
REAL(kind=8) O_den_p_l21_log
REAL(kind=8) O_den_p_21_log
REAL(kind=8) O2_den_p_u21_log
REAL(kind=8) O2_den_p_l21_log
REAL(kind=8) O2_den_p_21_log
REAL(kind=8) N2_den_p_u21_log
REAL(kind=8) N2_den_p_l21_log
REAL(kind=8) N2_den_p_21_log
REAL(kind=8) HYD_den_p_u21_log
REAL(kind=8) HYD_den_p_l21_log
REAL(kind=8) HYD_den_p_21_log
REAL(kind=8) HEL_den_p_u21_log
REAL(kind=8) HEL_den_p_l21_log
REAL(kind=8) HEL_den_p_21_log
REAL(kind=8) NO_den_p_u21_log
REAL(kind=8) NO_den_p_l21_log
REAL(kind=8) NO_den_p_21_log
REAL(kind=8) N4S_den_p_u21_log
REAL(kind=8) N4S_den_p_l21_log
REAL(kind=8) N4S_den_p_21_log
REAL(kind=8) N2D_den_p_u21_log
REAL(kind=8) N2D_den_p_l21_log
REAL(kind=8) N2D_den_p_21_log
REAL(kind=8) Tn_p_u21
REAL(kind=8) Tn_p_l21
REAL(kind=8) Tn_p_21
REAL(kind=8) Un_p_u21
REAL(kind=8) Un_p_l21
REAL(kind=8) Un_p_21
REAL(kind=8) Vn_p_u21
REAL(kind=8) Vn_p_l21
REAL(kind=8) Vn_p_21
REAL(kind=8) qion3d_p_u21
REAL(kind=8) qion3d_p_l21
REAL(kind=8) qion3d_p_21
REAL(kind=8) elx_p_u21
REAL(kind=8) elx_p_l21
REAL(kind=8) elx_p_21
REAL(kind=8) ely_p_u21
REAL(kind=8) ely_p_l21
REAL(kind=8) ely_p_21
REAL(kind=8) qo2p_aurora_p_u21
REAL(kind=8) qo2p_aurora_p_l21
REAL(kind=8) qo2p_aurora_p_21
REAL(kind=8) qop_aurora_p_u21
REAL(kind=8) qop_aurora_p_l21
REAL(kind=8) qop_aurora_p_21
REAL(kind=8) qn2p_aurora_p_u21
REAL(kind=8) qn2p_aurora_p_l21
REAL(kind=8) qn2p_aurora_p_21
REAL(kind=8) qnp_aurora_p_u21
REAL(kind=8) qnp_aurora_p_l21
REAL(kind=8) qnp_aurora_p_21
REAL(kind=8) qtef_aurora_p_u21
REAL(kind=8) qtef_aurora_p_l21
REAL(kind=8) qtef_aurora_p_21


REAL(kind=8) O_den_p_2 
REAL(kind=8) O_den_p_1 
REAL(kind=8) O2_den_p_2 
REAL(kind=8) O2_den_p_1 
REAL(kind=8) N2_den_p_2
REAL(kind=8) N2_den_p_1
REAL(kind=8) HYD_den_p_2
REAL(kind=8) HYD_den_p_1
REAL(kind=8) HEL_den_p_2
REAL(kind=8) HEL_den_p_1
REAL(kind=8) NO_den_p_2
REAL(kind=8) NO_den_p_1
REAL(kind=8) N4S_den_p_2
REAL(kind=8) N4S_den_p_1
REAL(kind=8) N2D_den_p_2
REAL(kind=8) N2D_den_p_1
REAL(kind=8) Tn_p_2 
REAL(kind=8) Tn_p_1
REAL(kind=8) Un_p_2 
REAL(kind=8) Un_p_1
REAL(kind=8) Vn_p_2 
REAL(kind=8) Vn_p_1
REAL(kind=8) qion3d_p_2 
REAL(kind=8) qion3d_p_1
REAL(kind=8) elx_p_2
REAL(kind=8) elx_p_1
REAL(kind=8) ely_p_2
REAL(kind=8) ely_p_1
REAL(kind=8) qo2p_aurora_p_2
REAL(kind=8) qo2p_aurora_p_1
REAL(kind=8) qop_aurora_p_2
REAL(kind=8) qop_aurora_p_1
REAL(kind=8) qn2p_aurora_p_2
REAL(kind=8) qn2p_aurora_p_1
REAL(kind=8) qnp_aurora_p_2
REAL(kind=8) qnp_aurora_p_1
REAL(kind=8) qtef_aurora_p_2
REAL(kind=8) qtef_aurora_p_1


REAL(kind=8) interface_long(20)
REAL(kind=8) interface_lat(91)
REAL(kind=8) interface_height(interface_hts)


sw_External_model_provides_NO_N4S_densities = GIP_switches(5)
sw_input_Auroral_production_is_single_overall_rate = GIP_switches(7)

do ilon = 1 , 20 
  interface_long(ilon) = float(ilon-1)*18.
  if ( thermospheric_model_name == 'TGCM' ) then
     if ( interface_long(ilon) >  180. ) interface_long(ilon) = interface_long(ilon) - 360.
  endif
enddo

do ilat = 1 , 91
  interface_lat(ilat) = float(ilat-46) * 2.
enddo

do iht = 1 , interface_hts
!  interface_height(iht) = (float(iht-1)*5.) + 90.
  interface_height(iht) = fixed_heights_km(iht)
enddo



! Loop over interface longitudes....

do ilon = 1 , 20

!longitude interpolation

! loop over thermosphere longs to find points east and west....
   ispecial = 0
do ilon_therm = 1 , lon_dim
 if ( therm_long(ilon_therm) > interface_long(ilon)) then
   ilon_east = ilon_therm
   ilon_west = ilon_therm - 1
   goto 1500
 endif
enddo
   ilon_east = 1
   ilon_west = lon_dim
   ispecial = 1
1500 continue

if (ispecial == 0) then
factor_lon = (interface_long(ilon) - therm_long(ilon_west)) / (therm_long(ilon_east) - therm_long(ilon_west))
else
factor_lon = (interface_long(ilon) - therm_long(ilon_west)) / (therm_long(ilon_east)+360. - therm_long(ilon_west))
endif
factor_lon_array(ilon) = factor_lon
ilon_west_array(ilon) = ilon_west
ilon_east_array(ilon) = ilon_east

! Loop over interface latitudes....

do ilat = 1 , 91

! latitude interpolation....

do ilat_therm = 1 , lat_dim
 if ( therm_lat(ilat_therm) > interface_lat(ilat) ) then 
     ilat_north = ilat_therm
     ilat_south = ilat_therm - 1
   goto 2000
 endif
enddo 
2000 continue
if ( ilat_north == 1  ) then
ilat_north = 2
ilat_south = 1
endif

factor_lat = (interface_lat(ilat) - therm_lat(ilat_south)) / (therm_lat(ilat_north) - therm_lat(ilat_south))

if(factor_lat < 0.0) factor_lat = 0.0
if(factor_lat > 1.0) factor_lat = 1.0

factor_lat_array(ilat) = factor_lat
ilat_north_array(ilat) = ilat_north
ilat_south_array(ilat) = ilat_south







! Loop over interface heights....

do iht = 1 , interface_hts

! height interpolation....

do iht_therm = 1 , ht_dim
therm_height(iht_therm) = therm_Z(iht_therm,ilon_west,ilat_north) / 1000.
if ( therm_height(iht_therm) > interface_height(iht) ) then 
  iht_above = iht_therm
  iht_below = iht_therm - 1
  goto 2500
endif
enddo 
2500 continue
if ( iht_above == 1 ) then 
iht_above = 2
iht_below = 1
endif

factor_ht12 = (interface_height(iht) - therm_height(iht_below)) / (therm_height(iht_above) - therm_height(iht_below))

if( factor_ht12 < 0.0 ) factor_ht12 = 0.0

factor_ht_west_north(iht,ilat,ilon) = factor_ht12
iht_above_west_north(iht,ilat,ilon) = iht_above
iht_below_west_north(iht,ilat,ilon) = iht_below



do iht_therm = 1 , ht_dim
therm_height(iht_therm) = therm_Z(iht_therm,ilon_east,ilat_north) / 1000.
if ( therm_height(iht_therm) > interface_height(iht) ) then
  iht_above = iht_therm
  iht_below = iht_therm - 1
  goto 3000      
endif
enddo 
3000 continue
if ( iht_above == 1  ) then
iht_above = 2
iht_below = 1
endif

factor_ht22 = (interface_height(iht) - therm_height(iht_below)) / (therm_height(iht_above) - therm_height(iht_below))

if( factor_ht22 < 0.0 ) factor_ht22 = 0.0

factor_ht_east_north(iht,ilat,ilon) = factor_ht22
iht_above_east_north(iht,ilat,ilon) = iht_above
iht_below_east_north(iht,ilat,ilon) = iht_below



do iht_therm = 1 , ht_dim
therm_height(iht_therm) = therm_Z(iht_therm,ilon_west,ilat_south) / 1000.
if ( therm_height(iht_therm) > interface_height(iht) ) then
  iht_above = iht_therm
  iht_below = iht_therm - 1
  goto 3500     
endif
enddo 
3500 continue
if ( iht_above == 1 ) then
iht_above = 2
iht_below = 1
endif

factor_ht11 = (interface_height(iht) - therm_height(iht_below)) / (therm_height(iht_above) - therm_height(iht_below))

if( factor_ht11 < 0.0 ) factor_ht11 = 0.0

factor_ht_west_south(iht,ilat,ilon) = factor_ht11
iht_above_west_south(iht,ilat,ilon) = iht_above
iht_below_west_south(iht,ilat,ilon) = iht_below


do iht_therm = 1 , ht_dim
therm_height(iht_therm) = therm_Z(iht_therm,ilon_east,ilat_south) / 1000.
if ( therm_height(iht_therm) > interface_height(iht) ) then
  iht_above = iht_therm
  iht_below = iht_therm - 1
  goto 4000     
endif
enddo 
4000 continue
if ( iht_above == 1 ) then
iht_above = 2
iht_below = 1
endif

factor_ht21 = (interface_height(iht) - therm_height(iht_below)) / (therm_height(iht_above) - therm_height(iht_below))

if( factor_ht21 < 0.0 ) factor_ht21 = 0.0

factor_ht_east_south(iht,ilat,ilon) = factor_ht21
iht_above_east_south(iht,ilat,ilon) = iht_above
iht_below_east_south(iht,ilat,ilon) = iht_below



enddo 
enddo 
enddo 


! Loop over interface longitudes....

do ilon = 1 , 20
 !write(6,*) ilon
 factor_lon = factor_lon_array(ilon)
 ilon_west = ilon_west_array(ilon)
 ilon_east = ilon_east_array(ilon)

! Loop over interface latitudes....

do ilat = 1 , 91
  factor_lat = factor_lat_array(ilat)
  ilat_north = ilat_north_array(ilat)
  ilat_south = ilat_south_array(ilat)

! Loop over interface heights....

do iht = 1 , interface_hts


!    write(6,*) 'ted ',ilon,ilat,iht
! height interpolation....

   factor_ht12 = factor_ht_west_north(iht,ilat,ilon)
   iht_above = iht_above_west_north(iht,ilat,ilon)
   iht_below = iht_below_west_north(iht,ilat,ilon)

   O_den_p_u12_log=log10(therm_o_density(iht_above,ilon_west,ilat_north))
   O_den_p_l12_log=log10(therm_o_density(iht_below,ilon_west,ilat_north))
   O_den_p_12_log = ((O_den_p_u12_log - O_den_p_l12_log)*factor_ht12) + O_den_p_l12_log
   O_den_p_12_log = 10**(O_den_p_12_log)
   O2_den_p_u12_log=log10(therm_o2_density(iht_above,ilon_west,ilat_north))
   O2_den_p_l12_log=log10(therm_o2_density(iht_below,ilon_west,ilat_north))
   O2_den_p_12_log = ((O2_den_p_u12_log - O2_den_p_l12_log)*factor_ht12) + O2_den_p_l12_log
   O2_den_p_12_log = 10**(O2_den_p_12_log)
   N2_den_p_u12_log=log10(therm_n2_density(iht_above,ilon_west,ilat_north))
   N2_den_p_l12_log=log10(therm_n2_density(iht_below,ilon_west,ilat_north))
   N2_den_p_12_log = ((N2_den_p_u12_log - N2_den_p_l12_log)*factor_ht12) + N2_den_p_l12_log
   N2_den_p_12_log = 10**(N2_den_p_12_log)

if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_p_u12_log=log10(therm_NO_density(iht_above,ilon_west,ilat_north))
   NO_den_p_l12_log=log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_12_log = ((NO_den_p_u12_log - NO_den_p_l12_log)*factor_ht12) + NO_den_p_l12_log
   NO_den_p_12_log = 10**(NO_den_p_12_log)
   N4S_den_p_u12_log=log10(therm_N4S_density(iht_above,ilon_west,ilat_north))
   N4S_den_p_l12_log=log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_12_log = ((N4S_den_p_u12_log - N4S_den_p_l12_log)*factor_ht12) + N4S_den_p_l12_log
   N4S_den_p_12_log = 10**(N4S_den_p_12_log)
   N2D_den_p_u12_log=log10(therm_N2D_density(iht_above,ilon_west,ilat_north))
   N2D_den_p_l12_log=log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_12_log = ((N2D_den_p_u12_log - N2D_den_p_l12_log)*factor_ht12) + N2D_den_p_l12_log
   N2D_den_p_12_log = 10**(N2D_den_p_12_log)
endif

   Tn_p_u12=therm_Tn(iht_above,ilon_west,ilat_north)
   Tn_p_l12=therm_Tn(iht_below,ilon_west,ilat_north)
   Tn_p_12 = ((Tn_p_u12 - Tn_p_l12)*factor_ht12) + Tn_p_l12
   Un_p_u12=therm_Un(iht_above,ilon_west,ilat_north)
   Un_p_l12=therm_Un(iht_below,ilon_west,ilat_north)
   Un_p_12 = ((Un_p_u12 - Un_p_l12)*factor_ht12) + Un_p_l12
   Vn_p_u12=therm_Vn(iht_above,ilon_west,ilat_north)
   Vn_p_l12=therm_Vn(iht_below,ilon_west,ilat_north)
   Vn_p_12 = ((Vn_p_u12 - Vn_p_l12)*factor_ht12) + Vn_p_l12

if (sw_input_Auroral_production_is_single_overall_rate) then
   qion3d_p_u12=therm_qion3d(iht_above,ilon_west,ilat_north)
   qion3d_p_l12=therm_qion3d(iht_below,ilon_west,ilat_north)
   qion3d_p_12 = ((qion3d_p_u12 - qion3d_p_l12)*factor_ht12) + qion3d_p_l12
   if (qion3d_p_12 < 0.0) qion3d_p_12 = 0.0
else
   qo2p_aurora_p_u12=therm_qo2p_aurora(iht_above,ilon_west,ilat_north)
   qo2p_aurora_p_l12=therm_qo2p_aurora(iht_below,ilon_west,ilat_north)
   qo2p_aurora_p_12 = ((qo2p_aurora_p_u12 - qo2p_aurora_p_l12)*factor_ht12) + qo2p_aurora_p_l12
   if (qo2p_aurora_p_12 < 0.0) qo2p_aurora_p_12 = 0.0
   qop_aurora_p_u12=therm_qop_aurora(iht_above,ilon_west,ilat_north)
   qop_aurora_p_l12=therm_qop_aurora(iht_below,ilon_west,ilat_north)
   qop_aurora_p_12 = ((qop_aurora_p_u12 - qop_aurora_p_l12)*factor_ht12) + qop_aurora_p_l12
   if (qop_aurora_p_12 < 0.0) qop_aurora_p_12 = 0.0
   qn2p_aurora_p_u12=therm_qn2p_aurora(iht_above,ilon_west,ilat_north)
   qn2p_aurora_p_l12=therm_qn2p_aurora(iht_below,ilon_west,ilat_north)
   qn2p_aurora_p_12 = ((qn2p_aurora_p_u12 - qn2p_aurora_p_l12)*factor_ht12) + qn2p_aurora_p_l12
   if (qn2p_aurora_p_12 < 0.0) qn2p_aurora_p_12 = 0.0
   qnp_aurora_p_u12=therm_qnp_aurora(iht_above,ilon_west,ilat_north)
   qnp_aurora_p_l12=therm_qnp_aurora(iht_below,ilon_west,ilat_north)
   qnp_aurora_p_12 = ((qnp_aurora_p_u12 - qnp_aurora_p_l12)*factor_ht12) + qnp_aurora_p_l12
   if (qnp_aurora_p_12 < 0.0) qnp_aurora_p_12 = 0.0
   qtef_aurora_p_u12=therm_qtef_aurora(iht_above,ilon_west,ilat_north)
   qtef_aurora_p_l12=therm_qtef_aurora(iht_below,ilon_west,ilat_north)
   qtef_aurora_p_12 = ((qtef_aurora_p_u12 - qtef_aurora_p_l12)*factor_ht12) + qtef_aurora_p_l12
   if (qtef_aurora_p_12 < 0.0) qtef_aurora_p_12 = 0.0
endif

   elx_p_u12=therm_elx(iht_above,ilon_west,ilat_north)
   elx_p_l12=therm_elx(iht_below,ilon_west,ilat_north)
   elx_p_12 = ((elx_p_u12 - elx_p_l12)*factor_ht12) + elx_p_l12

   ely_p_u12=therm_ely(iht_above,ilon_west,ilat_north)
   ely_p_l12=therm_ely(iht_below,ilon_west,ilat_north)
   ely_p_12 = ((ely_p_u12 - ely_p_l12)*factor_ht12) + ely_p_l12


   factor_ht22 = factor_ht_east_north(iht,ilat,ilon) 
   iht_above = iht_above_east_north(iht,ilat,ilon) 
   iht_below = iht_below_east_north(iht,ilat,ilon)

   O_den_p_u22_log=log10(therm_o_density(iht_above,ilon_east,ilat_north))
   O_den_p_l22_log=log10(therm_o_density(iht_below,ilon_east,ilat_north))
   O_den_p_22_log = ((O_den_p_u22_log - O_den_p_l22_log)*factor_ht22) + O_den_p_l22_log
   O_den_p_22_log = 10**(O_den_p_22_log)
   O2_den_p_u22_log=log10(therm_o2_density(iht_above,ilon_east,ilat_north))
   O2_den_p_l22_log=log10(therm_o2_density(iht_below,ilon_east,ilat_north))
   O2_den_p_22_log = ((O2_den_p_u22_log - O2_den_p_l22_log)*factor_ht22) + O2_den_p_l22_log
   O2_den_p_22_log = 10**(O2_den_p_22_log)
   N2_den_p_u22_log=log10(therm_n2_density(iht_above,ilon_east,ilat_north))
   N2_den_p_l22_log=log10(therm_n2_density(iht_below,ilon_east,ilat_north))
   N2_den_p_22_log = ((N2_den_p_u22_log - N2_den_p_l22_log)*factor_ht22) + N2_den_p_l22_log
   N2_den_p_22_log = 10**(N2_den_p_22_log)

if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_p_u22_log=log10(therm_NO_density(iht_above,ilon_west,ilat_north))
   NO_den_p_l22_log=log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_22_log = ((NO_den_p_u22_log - NO_den_p_l22_log)*factor_ht22) + NO_den_p_l22_log
   NO_den_p_22_log = 10**(NO_den_p_22_log)
   N4S_den_p_u22_log=log10(therm_N4S_density(iht_above,ilon_west,ilat_north))
   N4S_den_p_l22_log=log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_22_log = ((N4S_den_p_u22_log - N4S_den_p_l22_log)*factor_ht22) + N4S_den_p_l22_log
   N4S_den_p_22_log = 10**(N4S_den_p_22_log)
   N2D_den_p_u22_log=log10(therm_N2D_density(iht_above,ilon_west,ilat_north))
   N2D_den_p_l22_log=log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_22_log = ((N2D_den_p_u22_log - N2D_den_p_l22_log)*factor_ht22) + N2D_den_p_l22_log
   N2D_den_p_22_log = 10**(N2D_den_p_22_log)
endif

   Tn_p_u22=therm_Tn(iht_above,ilon_east,ilat_north)
   Tn_p_l22=therm_Tn(iht_below,ilon_east,ilat_north)
   Tn_p_22 = ((Tn_p_u22 - Tn_p_l22)*factor_ht22) + Tn_p_l22
   Un_p_u22=therm_Un(iht_above,ilon_east,ilat_north)
   Un_p_l22=therm_Un(iht_below,ilon_east,ilat_north)
   Un_p_22 = ((Un_p_u22 - Un_p_l22)*factor_ht22) + Un_p_l22
   Vn_p_u22=therm_Vn(iht_above,ilon_east,ilat_north)
   Vn_p_l22=therm_Vn(iht_below,ilon_east,ilat_north)
   Vn_p_22 = ((Vn_p_u22 - Vn_p_l22)*factor_ht22) + Vn_p_l22

if (sw_input_Auroral_production_is_single_overall_rate) then
   qion3d_p_u22=therm_qion3d(iht_above,ilon_east,ilat_north)
   qion3d_p_l22=therm_qion3d(iht_below,ilon_east,ilat_north)
   qion3d_p_22 = ((qion3d_p_u22 - qion3d_p_l22)*factor_ht22) + qion3d_p_l22
   if (qion3d_p_22 < 0.0) qion3d_p_22 = 0.0
else
   qo2p_aurora_p_u22=therm_qo2p_aurora(iht_above,ilon_east,ilat_north)
   qo2p_aurora_p_l22=therm_qo2p_aurora(iht_below,ilon_east,ilat_north)
   qo2p_aurora_p_22 = ((qo2p_aurora_p_u22 - qo2p_aurora_p_l22)*factor_ht22) + qo2p_aurora_p_l22
   if (qo2p_aurora_p_22 < 0.0) qo2p_aurora_p_22 = 0.0
   qop_aurora_p_u22=therm_qop_aurora(iht_above,ilon_east,ilat_north)
   qop_aurora_p_l22=therm_qop_aurora(iht_below,ilon_east,ilat_north)
   qop_aurora_p_22 = ((qop_aurora_p_u22 - qop_aurora_p_l22)*factor_ht22) + qop_aurora_p_l22
   if (qop_aurora_p_22 < 0.0) qop_aurora_p_22 = 0.0
   qn2p_aurora_p_u22=therm_qn2p_aurora(iht_above,ilon_east,ilat_north)
   qn2p_aurora_p_l22=therm_qn2p_aurora(iht_below,ilon_east,ilat_north)
   qn2p_aurora_p_22 = ((qn2p_aurora_p_u22 - qn2p_aurora_p_l22)*factor_ht22) + qn2p_aurora_p_l22
   if (qn2p_aurora_p_22 < 0.0) qn2p_aurora_p_22 = 0.0
   qnp_aurora_p_u22=therm_qnp_aurora(iht_above,ilon_east,ilat_north)
   qnp_aurora_p_l22=therm_qnp_aurora(iht_below,ilon_east,ilat_north)
   qnp_aurora_p_22 = ((qnp_aurora_p_u22 - qnp_aurora_p_l22)*factor_ht22) + qnp_aurora_p_l22
   if (qnp_aurora_p_22 < 0.0) qnp_aurora_p_22 = 0.0
   qtef_aurora_p_u22=therm_qtef_aurora(iht_above,ilon_east,ilat_north)
   qtef_aurora_p_l22=therm_qtef_aurora(iht_below,ilon_east,ilat_north)
   qtef_aurora_p_22 = ((qtef_aurora_p_u22 - qtef_aurora_p_l22)*factor_ht22) + qtef_aurora_p_l22
   if (qtef_aurora_p_22 < 0.0) qtef_aurora_p_22 = 0.0
endif

   elx_p_u22=therm_elx(iht_above,ilon_east,ilat_north)
   elx_p_l22=therm_elx(iht_below,ilon_east,ilat_north)
   elx_p_22 = ((elx_p_u22 - elx_p_l22)*factor_ht22) + elx_p_l22
   ely_p_u22=therm_ely(iht_above,ilon_east,ilat_north)
   ely_p_l22=therm_ely(iht_below,ilon_east,ilat_north)
   ely_p_22 = ((ely_p_u22 - ely_p_l22)*factor_ht22) + ely_p_l22

   factor_ht11 = factor_ht_west_south(iht,ilat,ilon) 
   iht_above = iht_above_west_south(iht,ilat,ilon) 
   iht_below = iht_below_west_south(iht,ilat,ilon)

   O_den_p_u11_log=log10(therm_o_density(iht_above,ilon_west,ilat_south))
   O_den_p_l11_log=log10(therm_o_density(iht_below,ilon_west,ilat_south))
   O_den_p_11_log = ((O_den_p_u11_log - O_den_p_l11_log)*factor_ht11) + O_den_p_l11_log
   O_den_p_11_log = 10**(O_den_p_11_log)
   O2_den_p_u11_log=log10(therm_o2_density(iht_above,ilon_west,ilat_south))
   O2_den_p_l11_log=log10(therm_o2_density(iht_below,ilon_west,ilat_south))
   O2_den_p_11_log = ((O2_den_p_u11_log - O2_den_p_l11_log)*factor_ht11) + O2_den_p_l11_log
   O2_den_p_11_log = 10**(O2_den_p_11_log)
   N2_den_p_u11_log=log10(therm_n2_density(iht_above,ilon_west,ilat_south))
   N2_den_p_l11_log=log10(therm_n2_density(iht_below,ilon_west,ilat_south))
   N2_den_p_11_log = ((N2_den_p_u11_log - N2_den_p_l11_log)*factor_ht11) + N2_den_p_l11_log
   N2_den_p_11_log = 10**(N2_den_p_11_log)

if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_p_u11_log=log10(therm_NO_density(iht_above,ilon_west,ilat_north))
   NO_den_p_l11_log=log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_11_log = ((NO_den_p_u11_log - NO_den_p_l11_log)*factor_ht11) + NO_den_p_l11_log
   NO_den_p_11_log = 10**(NO_den_p_11_log)
   N4S_den_p_u11_log=log10(therm_N4S_density(iht_above,ilon_west,ilat_north))
   N4S_den_p_l11_log=log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_11_log = ((N4S_den_p_u11_log - N4S_den_p_l11_log)*factor_ht11) + N4S_den_p_l11_log
   N4S_den_p_11_log = 10**(N4S_den_p_11_log)
   N2D_den_p_u11_log=log10(therm_N2D_density(iht_above,ilon_west,ilat_north))
   N2D_den_p_l11_log=log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_11_log = ((N2D_den_p_u11_log - N2D_den_p_l11_log)*factor_ht11) + N2D_den_p_l11_log
   N2D_den_p_11_log = 10**(N2D_den_p_11_log)
endif

   Tn_p_u11=therm_Tn(iht_above,ilon_west,ilat_south)
   Tn_p_l11=therm_Tn(iht_below,ilon_west,ilat_south)
   Tn_p_11 = ((Tn_p_u11 - Tn_p_l11)*factor_ht11) + Tn_p_l11
   Un_p_u11=therm_Un(iht_above,ilon_west,ilat_south)
   Un_p_l11=therm_Un(iht_below,ilon_west,ilat_south)
   Un_p_11 = ((Un_p_u11 - Un_p_l11)*factor_ht11) + Un_p_l11
   Vn_p_u11=therm_Vn(iht_above,ilon_west,ilat_south)
   Vn_p_l11=therm_Vn(iht_below,ilon_west,ilat_south)
   Vn_p_11 = ((Vn_p_u11 - Vn_p_l11)*factor_ht11) + Vn_p_l11

if (sw_input_Auroral_production_is_single_overall_rate) then
   qion3d_p_u11=therm_qion3d(iht_above,ilon_west,ilat_south)
   qion3d_p_l11=therm_qion3d(iht_below,ilon_west,ilat_south)
   qion3d_p_11 = ((qion3d_p_u11 - qion3d_p_l11)*factor_ht11) + qion3d_p_l11
   if (qion3d_p_11 < 0.0) qion3d_p_11 = 0.0
else
   qo2p_aurora_p_u11=therm_qo2p_aurora(iht_above,ilon_west,ilat_south)
   qo2p_aurora_p_l11=therm_qo2p_aurora(iht_below,ilon_west,ilat_south)
   qo2p_aurora_p_11 = ((qo2p_aurora_p_u11 - qo2p_aurora_p_l11)*factor_ht11) + qo2p_aurora_p_l11
   if (qo2p_aurora_p_11 < 0.0) qo2p_aurora_p_11 = 0.0
   qop_aurora_p_u11=therm_qop_aurora(iht_above,ilon_west,ilat_south)
   qop_aurora_p_l11=therm_qop_aurora(iht_below,ilon_west,ilat_south)
   qop_aurora_p_11 = ((qop_aurora_p_u11 - qop_aurora_p_l11)*factor_ht11) + qop_aurora_p_l11
   if (qop_aurora_p_11 < 0.0) qop_aurora_p_11 = 0.0
   qn2p_aurora_p_u11=therm_qn2p_aurora(iht_above,ilon_west,ilat_south)
   qn2p_aurora_p_l11=therm_qn2p_aurora(iht_below,ilon_west,ilat_south)
   qn2p_aurora_p_11 = ((qn2p_aurora_p_u11 - qn2p_aurora_p_l11)*factor_ht11) + qn2p_aurora_p_l11
   if (qn2p_aurora_p_11 < 0.0) qn2p_aurora_p_11 = 0.0
   qnp_aurora_p_u11=therm_qnp_aurora(iht_above,ilon_west,ilat_south)
   qnp_aurora_p_l11=therm_qnp_aurora(iht_below,ilon_west,ilat_south)
   qnp_aurora_p_11 = ((qnp_aurora_p_u11 - qnp_aurora_p_l11)*factor_ht11) + qnp_aurora_p_l11
   if (qnp_aurora_p_11 < 0.0) qnp_aurora_p_11 = 0.0
   qtef_aurora_p_u11=therm_qtef_aurora(iht_above,ilon_west,ilat_south)
   qtef_aurora_p_l11=therm_qtef_aurora(iht_below,ilon_west,ilat_south)
   qtef_aurora_p_11 = ((qtef_aurora_p_u11 - qtef_aurora_p_l11)*factor_ht11) + qtef_aurora_p_l11
   if (qtef_aurora_p_11 < 0.0) qtef_aurora_p_11 = 0.0
endif

   elx_p_u11=therm_elx(iht_above,ilon_west,ilat_south)
   elx_p_l11=therm_elx(iht_below,ilon_west,ilat_south)
   elx_p_11 = ((elx_p_u11 - elx_p_l11)*factor_ht11) + elx_p_l11
   ely_p_u11=therm_ely(iht_above,ilon_west,ilat_south)
   ely_p_l11=therm_ely(iht_below,ilon_west,ilat_south)
   ely_p_11 = ((ely_p_u11 - ely_p_l11)*factor_ht11) + ely_p_l11


   factor_ht21 = factor_ht_east_south(iht,ilat,ilon) 
   iht_above = iht_above_east_south(iht,ilat,ilon)
   iht_below = iht_below_east_south(iht,ilat,ilon)

   O_den_p_u21_log=log10(therm_o_density(iht_above,ilon_east,ilat_south))
   O_den_p_l21_log=log10(therm_o_density(iht_below,ilon_east,ilat_south))
   O_den_p_21_log = ((O_den_p_u21_log - O_den_p_l21_log)*factor_ht21) + O_den_p_l21_log
   O_den_p_21_log = 10**(O_den_p_21_log)
   O2_den_p_u21_log=log10(therm_o2_density(iht_above,ilon_east,ilat_south))
   O2_den_p_l21_log=log10(therm_o2_density(iht_below,ilon_east,ilat_south))
   O2_den_p_21_log = ((O2_den_p_u21_log - O2_den_p_l21_log)*factor_ht21) + O2_den_p_l21_log
   O2_den_p_21_log = 10**(O2_den_p_21_log)
   N2_den_p_u21_log=log10(therm_n2_density(iht_above,ilon_east,ilat_south))
   N2_den_p_l21_log=log10(therm_n2_density(iht_below,ilon_east,ilat_south))
   N2_den_p_21_log = ((N2_den_p_u21_log - N2_den_p_l21_log)*factor_ht21) + N2_den_p_l21_log
   N2_den_p_21_log = 10**(N2_den_p_21_log)

if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_p_u21_log=log10(therm_NO_density(iht_above,ilon_west,ilat_north))
   NO_den_p_l21_log=log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_21_log = ((NO_den_p_u21_log - NO_den_p_l21_log)*factor_ht21) + NO_den_p_l21_log
   NO_den_p_21_log = 10**(NO_den_p_21_log)
   N4S_den_p_u21_log=log10(therm_N4S_density(iht_above,ilon_west,ilat_north))
   N4S_den_p_l21_log=log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_21_log = ((N4S_den_p_u21_log - N4S_den_p_l21_log)*factor_ht21) + N4S_den_p_l21_log
   N4S_den_p_21_log = 10**(N4S_den_p_21_log)
   N2D_den_p_u21_log=log10(therm_N2D_density(iht_above,ilon_west,ilat_north))
   N2D_den_p_l21_log=log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_21_log = ((N2D_den_p_u21_log - N2D_den_p_l21_log)*factor_ht21) + N2D_den_p_l21_log
   N2D_den_p_21_log = 10**(N2D_den_p_21_log)
endif

   Tn_p_u21=therm_Tn(iht_above,ilon_east,ilat_south)
   Tn_p_l21=therm_Tn(iht_below,ilon_east,ilat_south)
   Tn_p_21 = ((Tn_p_u21 - Tn_p_l21)*factor_ht21) + Tn_p_l21
   Un_p_u21=therm_Un(iht_above,ilon_east,ilat_south)
   Un_p_l21=therm_Un(iht_below,ilon_east,ilat_south)
   Un_p_21 = ((Un_p_u21 - Un_p_l21)*factor_ht21) + Un_p_l21
   Vn_p_u21=therm_Vn(iht_above,ilon_east,ilat_south)
   Vn_p_l21=therm_Vn(iht_below,ilon_east,ilat_south)
   Vn_p_21 = ((Vn_p_u21 - Vn_p_l21)*factor_ht21) + Vn_p_l21

if (sw_input_Auroral_production_is_single_overall_rate) then
   qion3d_p_u21=therm_qion3d(iht_above,ilon_east,ilat_south)
   qion3d_p_l21=therm_qion3d(iht_below,ilon_east,ilat_south)
   qion3d_p_21 = ((qion3d_p_u21 - qion3d_p_l21)*factor_ht21) + qion3d_p_l21
   if (qion3d_p_21 < 0.0) qion3d_p_21 = 0.0
else
   qo2p_aurora_p_u21=therm_qo2p_aurora(iht_above,ilon_east,ilat_south)
   qo2p_aurora_p_l21=therm_qo2p_aurora(iht_below,ilon_east,ilat_south)
   qo2p_aurora_p_21 = ((qo2p_aurora_p_u21 - qo2p_aurora_p_l21)*factor_ht21) + qo2p_aurora_p_l21
   if (qo2p_aurora_p_21 < 0.0) qo2p_aurora_p_21 = 0.0
   qop_aurora_p_u21=therm_qop_aurora(iht_above,ilon_east,ilat_south)
   qop_aurora_p_l21=therm_qop_aurora(iht_below,ilon_east,ilat_south)
   qop_aurora_p_21 = ((qop_aurora_p_u21 - qop_aurora_p_l21)*factor_ht21) + qop_aurora_p_l21
   if (qop_aurora_p_21 < 0.0) qop_aurora_p_21 = 0.0
   qn2p_aurora_p_u21=therm_qn2p_aurora(iht_above,ilon_east,ilat_south)
   qn2p_aurora_p_l21=therm_qn2p_aurora(iht_below,ilon_east,ilat_south)
   qn2p_aurora_p_21 = ((qn2p_aurora_p_u21 - qn2p_aurora_p_l21)*factor_ht21) + qn2p_aurora_p_l21
   if (qn2p_aurora_p_21 < 0.0) qn2p_aurora_p_21 = 0.0
   qnp_aurora_p_u21=therm_qnp_aurora(iht_above,ilon_east,ilat_south)
   qnp_aurora_p_l21=therm_qnp_aurora(iht_below,ilon_east,ilat_south)
   qnp_aurora_p_21 = ((qnp_aurora_p_u21 - qnp_aurora_p_l21)*factor_ht21) + qnp_aurora_p_l21
   if (qnp_aurora_p_21 < 0.0) qnp_aurora_p_21 = 0.0
   qtef_aurora_p_u21=therm_qtef_aurora(iht_above,ilon_east,ilat_south)
   qtef_aurora_p_l21=therm_qtef_aurora(iht_below,ilon_east,ilat_south)
   qtef_aurora_p_21 = ((qtef_aurora_p_u21 - qtef_aurora_p_l21)*factor_ht21) + qtef_aurora_p_l21
   if (qtef_aurora_p_21 < 0.0) qtef_aurora_p_21 = 0.0
endif
   elx_p_u21=therm_elx(iht_above,ilon_east,ilat_south)
   elx_p_l21=therm_elx(iht_below,ilon_east,ilat_south)
   elx_p_21 = ((elx_p_u21 - elx_p_l21)*factor_ht21) + elx_p_l21
   ely_p_u21=therm_ely(iht_above,ilon_east,ilat_south)
   ely_p_l21=therm_ely(iht_below,ilon_east,ilat_south)
   ely_p_21 = ((ely_p_u21 - ely_p_l21)*factor_ht21) + ely_p_l21


! latitude interpolation....

   O_den_p_2 = ((O_den_p_22_log - O_den_p_21_log) * factor_lat) + O_den_p_21_log
   O_den_p_1 = ((O_den_p_12_log - O_den_p_11_log) * factor_lat) + O_den_p_11_log
   O2_den_p_2 = ((O2_den_p_22_log - O2_den_p_21_log) * factor_lat) + O2_den_p_21_log
   O2_den_p_1 = ((O2_den_p_12_log - O2_den_p_11_log) * factor_lat) + O2_den_p_11_log
   N2_den_p_2 = ((N2_den_p_22_log - N2_den_p_21_log) * factor_lat) + N2_den_p_21_log
   N2_den_p_1 = ((N2_den_p_12_log - N2_den_p_11_log) * factor_lat) + N2_den_p_11_log

if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_p_2 = ((NO_den_p_22_log - NO_den_p_21_log) * factor_lat) + NO_den_p_21_log
   NO_den_p_1 = ((NO_den_p_12_log - NO_den_p_11_log) * factor_lat) + NO_den_p_11_log
   N4S_den_p_2 = ((N4S_den_p_22_log - N4S_den_p_21_log) * factor_lat) + N4S_den_p_21_log
   N4S_den_p_1 = ((N4S_den_p_12_log - N4S_den_p_11_log) * factor_lat) + N4S_den_p_11_log
   N2D_den_p_2 = ((N2D_den_p_22_log - N2D_den_p_21_log) * factor_lat) + N2D_den_p_21_log
   N2D_den_p_1 = ((N2D_den_p_12_log - N2D_den_p_11_log) * factor_lat) + N2D_den_p_11_log
endif
   Tn_p_2 = ((Tn_p_22 - Tn_p_21) * factor_lat) + Tn_p_21
   Tn_p_1 = ((Tn_p_12 - Tn_p_11) * factor_lat) + Tn_p_11
   Un_p_2 = ((Un_p_22 - Un_p_21) * factor_lat) + Un_p_21
   Un_p_1 = ((Un_p_12 - Un_p_11) * factor_lat) + Un_p_11
   Vn_p_2 = ((Vn_p_22 - Vn_p_21) * factor_lat) + Vn_p_21
   Vn_p_1 = ((Vn_p_12 - Vn_p_11) * factor_lat) + Vn_p_11

if (sw_input_Auroral_production_is_single_overall_rate) then
   qion3d_p_2 = ((qion3d_p_22 - qion3d_p_21) * factor_lat) + qion3d_p_21
   qion3d_p_1 = ((qion3d_p_12 - qion3d_p_11) * factor_lat) + qion3d_p_11
   if (qion3d_p_2 < 0.0) qion3d_p_2 = 0.0
   if (qion3d_p_1 < 0.0) qion3d_p_1 = 0.0
else
   qo2p_aurora_p_2 = ((qo2p_aurora_p_22 - qo2p_aurora_p_21) * factor_lat) + qo2p_aurora_p_21
   qo2p_aurora_p_1 = ((qo2p_aurora_p_12 - qo2p_aurora_p_11) * factor_lat) + qo2p_aurora_p_11
   if (qo2p_aurora_p_2 < 0.0) qo2p_aurora_p_2 = 0.0
   if (qo2p_aurora_p_1 < 0.0) qo2p_aurora_p_1 = 0.0
   qop_aurora_p_2 = ((qop_aurora_p_22 - qop_aurora_p_21) * factor_lat) + qop_aurora_p_21
   qop_aurora_p_1 = ((qop_aurora_p_12 - qop_aurora_p_11) * factor_lat) + qop_aurora_p_11
   if (qop_aurora_p_2 < 0.0) qop_aurora_p_2 = 0.0
   if (qop_aurora_p_1 < 0.0) qop_aurora_p_1 = 0.0
   qn2p_aurora_p_2 = ((qn2p_aurora_p_22 - qn2p_aurora_p_21) * factor_lat) + qn2p_aurora_p_21
   qn2p_aurora_p_1 = ((qn2p_aurora_p_12 - qn2p_aurora_p_11) * factor_lat) + qn2p_aurora_p_11
   if (qn2p_aurora_p_2 < 0.0) qn2p_aurora_p_2 = 0.0
   if (qn2p_aurora_p_1 < 0.0) qn2p_aurora_p_1 = 0.0
   qnp_aurora_p_2 = ((qnp_aurora_p_22 - qnp_aurora_p_21) * factor_lat) + qnp_aurora_p_21
   qnp_aurora_p_1 = ((qnp_aurora_p_12 - qnp_aurora_p_11) * factor_lat) + qnp_aurora_p_11
   if (qnp_aurora_p_2 < 0.0) qnp_aurora_p_2 = 0.0
   if (qnp_aurora_p_1 < 0.0) qnp_aurora_p_1 = 0.0
   qtef_aurora_p_2 = ((qtef_aurora_p_22 - qtef_aurora_p_21) * factor_lat) + qtef_aurora_p_21
   qtef_aurora_p_1 = ((qtef_aurora_p_12 - qtef_aurora_p_11) * factor_lat) + qtef_aurora_p_11
   if (qtef_aurora_p_2 < 0.0) qtef_aurora_p_2 = 0.0
   if (qtef_aurora_p_1 < 0.0) qtef_aurora_p_1 = 0.0
endif

   elx_p_2 = ((elx_p_22 - elx_p_21) * factor_lat) + elx_p_21
   elx_p_1 = ((elx_p_12 - elx_p_11) * factor_lat) + elx_p_11
   ely_p_2 = ((ely_p_22 - ely_p_21) * factor_lat) + ely_p_21
   ely_p_1 = ((ely_p_12 - ely_p_11) * factor_lat) + ely_p_11


! longitude interpolation....

   interface_o_density(iht,ilat,ilon) = ((O_den_p_2 - O_den_p_1) * factor_lon) + O_den_p_1
   interface_o2_density(iht,ilat,ilon) = ((O2_den_p_2 - O2_den_p_1) * factor_lon) + O2_den_p_1
   interface_n2_density(iht,ilat,ilon) = ((N2_den_p_2 - N2_den_p_1) * factor_lon) + N2_den_p_1
   if(interface_o_density(iht,ilat,ilon).lt.1.e-20)interface_o_density(iht,ilat,ilon)=1.e-20
   if(interface_o2_density(iht,ilat,ilon).lt.1.e-20)interface_o2_density(iht,ilat,ilon)=1.e-20
   if(interface_n2_density(iht,ilat,ilon).lt.1.e-20)interface_n2_density(iht,ilat,ilon)=1.e-20

if (sw_External_model_provides_NO_N4S_densities) then
   interface_NO_density(iht,ilat,ilon) = ((NO_den_p_2 - NO_den_p_1) * factor_lon) + NO_den_p_1
   interface_N4S_density(iht,ilat,ilon) = ((N4S_den_p_2 - N4S_den_p_1) * factor_lon) + N4S_den_p_1
   interface_N2D_density(iht,ilat,ilon) = ((N2D_den_p_2 - N2D_den_p_1) * factor_lon) + N2D_den_p_1
   if(interface_NO_density(iht,ilat,ilon).lt.1.e-20)interface_NO_density(iht,ilat,ilon)=1.e-20
   if(interface_N4S_density(iht,ilat,ilon).lt.1.e-20)interface_N4S_density(iht,ilat,ilon)=1.e-20
   if(interface_N2D_density(iht,ilat,ilon).lt.1.e-20)interface_N2D_density(iht,ilat,ilon)=1.e-20
endif
   interface_Tn(iht,ilat,ilon) = ((Tn_p_2 - Tn_p_1) * factor_lon) + Tn_p_1
   interface_Un(iht,ilat,ilon) = ((Un_p_2 - Un_p_1) * factor_lon) + Un_p_1
   interface_Vn(iht,ilat,ilon) = ((Vn_p_2 - Vn_p_1) * factor_lon) + Vn_p_1
if (sw_input_Auroral_production_is_single_overall_rate) then
   interface_qion3d(iht,ilat,ilon) = ((qion3d_p_2 - qion3d_p_1) * factor_lon) + qion3d_p_1
   if (interface_qion3d(iht,ilat,ilon) < 0.0) interface_qion3d(iht,ilat,ilon) = 0.0
else
   interface_qo2p_aurora(iht,ilat,ilon) = ((qo2p_aurora_p_2 - qo2p_aurora_p_1) * factor_lon) +qo2p_aurora_p_1
   if (interface_qo2p_aurora(iht,ilat,ilon) < 0.0) interface_qo2p_aurora(iht,ilat,ilon) = 0.0
   interface_qop_aurora(iht,ilat,ilon) = ((qop_aurora_p_2 - qop_aurora_p_1) * factor_lon) + qop_aurora_p_1
   if (interface_qop_aurora(iht,ilat,ilon) < 0.0) interface_qop_aurora(iht,ilat,ilon) = 0.0
   interface_qn2p_aurora(iht,ilat,ilon) = ((qn2p_aurora_p_2 - qn2p_aurora_p_1) * factor_lon) + qn2p_aurora_p_1
   if (interface_qn2p_aurora(iht,ilat,ilon) < 0.0) interface_qn2p_aurora(iht,ilat,ilon) = 0.0
   interface_qnp_aurora(iht,ilat,ilon) = ((qnp_aurora_p_2 - qnp_aurora_p_1) * factor_lon) + qnp_aurora_p_1
   if (interface_qnp_aurora(iht,ilat,ilon) < 0.0) interface_qnp_aurora(iht,ilat,ilon) = 0.0
   interface_qtef_aurora(iht,ilat,ilon) = ((qtef_aurora_p_2 - qtef_aurora_p_1) * factor_lon) + qtef_aurora_p_1
   if (interface_qtef_aurora(iht,ilat,ilon) < 0.0) interface_qtef_aurora(iht,ilat,ilon) = 0.0
endif
   interface_elx(iht,ilat,ilon) = ((elx_p_2 - elx_p_1) * factor_lon)+ elx_p_1
   interface_ely(iht,ilat,ilon) = ((ely_p_2 - ely_p_1) * factor_lon)+ ely_p_1

enddo
enddo
enddo

      o_density_fixed_ht(:,:,:) = interface_o_density(:,:,:) 
      o2_density_fixed_ht(:,:,:) = interface_o2_density(:,:,:)
      n2_density_fixed_ht(:,:,:) = interface_n2_density(:,:,:)

if (sw_External_model_provides_NO_N4S_densities) then
      NO_density_fixed_ht(:,:,:) = interface_NO_density(:,:,:)
      N4S_density_fixed_ht(:,:,:) = interface_N4S_density(:,:,:)
      N2D_density_fixed_ht(:,:,:) = interface_N2D_density(:,:,:)
endif
      Vx_fixed_ht(:,:,:) = interface_Vn(:,:,:)
      Vy_fixed_ht(:,:,:) = interface_Un(:,:,:)
      wvz_fixed_ht(:,:,:) = 0.0
      tts_fixed_ht(:,:,:) = interface_Tn(:,:,:)
if (sw_input_Auroral_production_is_single_overall_rate) then
      qion3d_fixed_ht(:,:,:) = interface_qion3d(:,:,:)
else
      qo2p_aurora_fixed_ht(:,:,:) = interface_qo2p_aurora(:,:,:)
      qop_aurora_fixed_ht(:,:,:)  = interface_qop_aurora(:,:,:)
      qn2p_aurora_fixed_ht(:,:,:) = interface_qn2p_aurora(:,:,:)
      qnp_aurora_fixed_ht(:,:,:)  = interface_qnp_aurora(:,:,:)
      qtef_aurora_fixed_ht(:,:,:) = interface_qtef_aurora(:,:,:)
endif
      elx_fixed_ht(:,:,:) = interface_elx(:,:,:)
      ely_fixed_ht(:,:,:) = interface_ely(:,:,:)

      do iht = 1 , interface_hts
      do ilat = 1 , 91
      do ilon = 1 , 20
      if(o2_density_fixed_ht(iht,ilat,ilon) < 0.0 ) then
         write(6,*) 'o2_density negative ',iht,ilat,ilon,o2_density_fixed_ht(iht,ilat,ilon)
         stop
      endif
      enddo
      enddo
      enddo

return




end SUBROUTINE INTERFACE__thermosphere_to_GIP









SUBROUTINE INTERFACE__GIP_to_thermosphere ( &
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
      !write(6,*) 'interface longitudes ',ilon_int , high_res_long(ilon_int)
    enddo
    do ilat_int = 1 , 91
      high_res_lat(ilat_int)=(float(ilat_int-46)) * 2.
      !write(6,*) 'interface latitudes ',ilat_int , high_res_lat(ilat_int)
    enddo
    do iht_int = 1 , interface_hts
!      high_res_height(iht_int)=((float(iht_int-1)) * 5.) + 90.
      high_res_height(iht_int)=fixed_heights_km(iht_int)
      !write(6,*) 'interface heights ',iht_int , high_res_height(iht_int)
    enddo

! Loop over therm longitudes....

    do ilon = 1 , lon_dim 
      !write(6,*) 'ilon ',ilon

!longitude interpolation

! loop over high_res longs to find points east and west....

       ispecial = 0
    do ilon_int = 1 , 90
      !write(6,*) '     ilon_int ',ilon_int
      !write(6,*) '  longs ',high_res_long(ilon_int),therm_geo_long(ilon)
      if (high_res_long(ilon_int) > therm_geo_long(ilon)) then
        ilon_east = ilon_int
        ilon_west = ilon_int - 1
      !write(6,*) '     ieast iwest ',ilon_east,ilon_west
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

!write(6,*) 'ilon ' , ilon,ilon_east,ilon_west,factor_lon
!write(6,*) 'lon confusion ',therm_geo_long(ilon),high_res_long(ilon_east),high_res_long(ilon_west)
! Loop over therm latitudes....

    do ilat = 1 , lat_dim
        !write(6,*) 'ilat ' , ilat

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


       !write(6,*) 'ilat ' , ilat,ilat_north,ilat_south,factor_lat





! Loop over therm heights....

  do iht = 1 , ht_dim           
       !write(6,*) '     iht ',iht

  therm_Z_km = therm_Z(iht,ilon,ilat) / 1000.

! height interpolation....

  do iht_int = 1 , interface_hts         
      !write(6,*) '     iht_int ',iht_int
      !write(6,*) '  heights ',high_res_height(iht_int),therm_Z(iht,ilon,ilat)
    if(high_res_height(iht_int) > therm_Z_km) then
      iht_above = iht_int
      iht_below = iht_int - 1
      goto 2500    
    endif
  enddo
  iht_above = interface_hts
  iht_below = interface_hts - 1
2500 continue

       !write(6,*) 'iht 1' , iht,iht_above,iht_below
if (iht_above == 1) then
iht_above = 2
iht_below = 1
endif
factor_ht = (therm_Z_km - high_res_height(iht_below)) /  &
          (high_res_height(iht_above) - high_res_height(iht_below))

! cg - make sure factor_ht doesn't get smaller than 0.0 which will happen for thermospheric
! cg - heights below 90km.....

if (factor_ht .lt. 0.0) factor_ht = 0.0

       !write(6,*) 'iht 2' , iht,iht_above,iht_below,factor_ht



    !write(6,*) iht_above,iht_below,ilat_north,ilat_south,ilon_east,ilon_west

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

! latitude interpolation....

ne_east = ((ne_north_east - ne_south_east) * factor_lat) + ne_south_east
ne_west = ((ne_north_west - ne_south_west) * factor_lat) + ne_south_west
oplus_east = ((oplus_north_east - oplus_south_east) * factor_lat) + oplus_south_east
oplus_west = ((oplus_north_west - oplus_south_west) * factor_lat) + oplus_south_west
hplus_east = ((hplus_north_east - hplus_south_east) * factor_lat) + hplus_south_east
hplus_west = ((hplus_north_west - hplus_south_west) * factor_lat) + hplus_south_west
noplus_east = ((noplus_north_east - noplus_south_east) * factor_lat) + noplus_south_east
noplus_west = ((noplus_north_west - noplus_south_west) * factor_lat) + noplus_south_west
o2plus_east = ((o2plus_north_east - o2plus_south_east) * factor_lat) + o2plus_south_east
o2plus_west = ((o2plus_north_west - o2plus_south_west) * factor_lat) + o2plus_south_west
n2plus_east = ((n2plus_north_east - n2plus_south_east) * factor_lat) + n2plus_south_east
n2plus_west = ((n2plus_north_west - n2plus_south_west) * factor_lat) + n2plus_south_west
nplus_east = ((nplus_north_east - nplus_south_east) * factor_lat) + nplus_south_east
nplus_west = ((nplus_north_west - nplus_south_west) * factor_lat) + nplus_south_west
Te_east = ((Te_north_east - Te_south_east) * factor_lat) + Te_south_east
Te_west = ((Te_north_west - Te_south_west) * factor_lat) + Te_south_west
Ti1_east = ((Ti1_north_east - Ti1_south_east) * factor_lat) + Ti1_south_east
Ti1_west = ((Ti1_north_west - Ti1_south_west) * factor_lat) + Ti1_south_west
Ti2_east = ((Ti2_north_east - Ti2_south_east) * factor_lat) + Ti2_south_east
Ti2_west = ((Ti2_north_west - Ti2_south_west) * factor_lat) + Ti2_south_west

! longitude interpolation....

therm_Ne_density(iht,ilon,ilat) = ((ne_east - ne_west) * factor_lon) + ne_west
therm_oplus_density(iht,ilon,ilat) = ((oplus_east - oplus_west) * factor_lon) + oplus_west
therm_hplus_density(iht,ilon,ilat) = ((hplus_east - hplus_west) * factor_lon) + hplus_west
therm_noplus_density(iht,ilon,ilat) = ((noplus_east - noplus_west) * factor_lon) + noplus_west
therm_o2plus_density(iht,ilon,ilat) = ((o2plus_east - o2plus_west) * factor_lon) + o2plus_west
therm_n2plus_density(iht,ilon,ilat) = ((n2plus_east - n2plus_west) * factor_lon) + n2plus_west
therm_nplus_density(iht,ilon,ilat) = ((nplus_east - nplus_west) * factor_lon) + nplus_west
therm_Te(iht,ilon,ilat) = ((Te_east - Te_west) * factor_lon) + Te_west
therm_Ti1(iht,ilon,ilat) = ((Ti1_east - Ti1_west) * factor_lon) + Ti1_west
therm_Ti2(iht,ilon,ilat) = ((Ti2_east - Ti2_west) * factor_lon) + Ti2_west

enddo
enddo
enddo


return




end SUBROUTINE INTERFACE__GIP_to_thermosphere




SUBROUTINE INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE( &
        geo_grid_longitudes_degrees, &
        geo_grid_latitudes_degrees, &
        O_density, O2_density, N2_density, &
        NO_density, &
        N4S_density, &
        N2D_density, &
        VX, VY, WVZ, TTS, telec, &
        IN, IS, IWRite2, &
        TN_plasma_input_3d, O_plasma_input_3d, &
        O2_plasma_input_3d, N2_plasma_input_3d, &
        NO_plasma_input_3d, N4S_plasma_input_3d, N2D_plasma_input_3d, &
        GLAt_3d, &
        GLOnd_3d, &
        PZ_3d, &
        um_plasma_input_3d, uz_plasma_input_3d, uv_plasma_input_3d, &
        te_plasma_input_3d, &
        ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &
        ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht, ihu_3d_fixed_ht, &
        sw_1st_call_int_fixed_ht, &
        GIP_switches)





! this calculates a neutral background for the
! plasmasphere code by interpolating values from the
! fixed height interface.....

  IMPLICIT NONE
  INTEGER :: N_heights,N_Latitudes,N_longitudes
  PARAMETER (N_heights=interface_hts,N_Latitudes=91,N_longitudes=20)

  LOGICAL :: sw_1st_call_int_fixed_ht
  LOGICAL :: GIP_switches(20)
  LOGICAL :: sw_External_model_provides_NO_N4S_densities

  REAL(kind=8) :: dNO1 , dNO11 , dNO12 , dNO2 , dNO21 , dNO22 , &
  dNOl11 , dNOl12 , dNOl21 , dNOl22 , dNOu11 , dNOu12, &
  dNOu21 , dNOu22
  REAL(kind=8) :: dN4S1 , dN4S11 , dN4S12 , dN4S2 , dN4S21 , dN4S22 , &
  dN4Sl11 , dN4Sl12 , dN4Sl21 , dN4Sl22 , dN4Su11 , dN4Su12, &
  dN4Su21 , dN4Su22
  REAL(kind=8) :: dN2D1 , dN2D11 , dN2D12 , dN2D2 , dN2D21 , dN2D22 , &
  dN2Dl11 , dN2Dl12 , dN2Dl21 , dN2Dl22 , dN2Du11 , dN2Du12, &
  dN2Du21 , dN2Du22

  REAL(kind=8) :: dnn1 , dnn11 , dnn12 , dnn2 , dnn21 , dnn22 , &
  dnnl11 , dnnl12 , dnnl21 , dnnl22 , dnnu11 , dnnu12
  REAL(kind=8) :: &
  dnnu21 , dnnu22 , fach,  &
  faclat , faclon , &
  glond2 

  REAL(kind=8) :: ol11 , ol12 , ol21 , ol22 ,  ou11 , ou12 , ou21 , ou22 , o11 , o12 , o21 , o22 , do1 , do2
  REAL(kind=8) :: ool11 , ool12 , ool21 , ool22 ,  oou11 , oou12 , oou21 , oou22 , oo11 , oo12 , oo21 , oo22  , doo1 , doo2
  REAL(kind=8) :: tnl11 , tnl12 , tnl21 , tnl22 ,  tnu11 , tnu12 , tnu21 , tnu22 , tn11 , tn12 , tn21 , tn22  , tn1 , tn2
  REAL(kind=8) :: tel11 , tel12 , tel21 , tel22 ,  teu11 , teu12 , teu21 , teu22 , te11 , te12 , te21 , te22 , te1 , te2
  REAL(kind=8) :: uml11 , uml12 , uml21 , uml22 ,  umu11 , umu12 , umu21 , umu22 , um11 , um12 , um21 , um22 , um1 , um2
  REAL(kind=8) :: uzl11 , uzl12 , uzl21 , uzl22 ,  uzu11 , uzu12 , uzu21 , uzu22 , uz11 , uz12 , uz21 , uz22 , uz1 , uz2
  REAL(kind=8) :: uvl11 , uvl12 , uvl21 , uvl22 ,  uvu11 , uvu12 , uvu21 , uvu22 , uv11 , uv12 , uv21 , uv22 , uv1 , uv2
  REAL(kind=8) :: topl11 , topl12 , topl21 , topl22 ,  topu11 , topu12 , topu21 , topu22 , top11 , top12 , top21 , top22  , top1 , top2
  REAL(kind=8) :: tnopl11 , tnopl12 , tnopl21 , tnopl22 ,  tnopu11 , tnopu12 , tnopu21 , tnopu22 , tnop11 , tnop12 , tnop21 , tnop22 , tnop1 , tnop2
  REAL(kind=8) :: to2pl11 , to2pl12 , to2pl21 , to2pl22 ,  to2pu11 , to2pu12 , to2pu21 , to2pu22 , to2p11 , to2p12 , to2p21 , to2p22 , to2p1 , to2p2

  INTEGER :: i , ih , ihl , ihu, ii , ilat , ilat1 , ilat2 , ilon , ilon1 , ilon2 , iprob
  INTEGER :: ispecial , IWRite , l , m , n ,  ifault , itube , iwrite1 , iwrite2 , istop

  REAL(kind=8) :: fixed_hts_in_km(N_heights)
  REAL(kind=8) :: O_density(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: O2_density(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: N2_density(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: NO_density(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: N4S_density(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: N2D_density(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: VX(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: VY(N_heights,N_Latitudes,N_longitudes) 
  REAL(kind=8) :: WVZ(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: TTS(N_heights,N_Latitudes,N_longitudes)
  REAL(kind=8) :: telec(N_heights,N_Latitudes,N_longitudes)

  REAL(kind=8) :: pzh(N_heights)

  integer :: in(nmp,nlp), is(nmp,nlp), mp, lp

  ! Field line grid , npts # of points along meriodional planes,
  ! nmp # of planes

  REAL(kind=8) tn_plasma_input_3d(npts,nmp), &
  um_plasma_input_3d(npts,nmp),uz_plasma_input_3d(npts,nmp),uv_plasma_input_3d(npts,nmp), &
  o_plasma_input_3d(npts,nmp), o2_plasma_input_3d(npts,nmp),n2_plasma_input_3d(npts,nmp), &
  te_plasma_input_3d(npts,nmp),te_dum(npts), &
  NO_plasma_input_3d(npts,nmp), &
  N4S_plasma_input_3d(npts,nmp), N2D_plasma_input_3d(npts,nmp)

  REAL(kind=8) ::  glat_3d(npts,nmp),glond_3d(npts,nmp),pz_3d(npts,nmp)
  REAL(kind=8) ::  NO(npts)
  REAL(kind=8) ::  N4S(npts)
  REAL(kind=8) ::  N2D(npts)

  REAL(kind=8) ::  uz(NPTS)
  REAL(kind=8) ::  um(NPTS)
  REAL(kind=8) ::  uv(NPTS)

  REAL(kind=8) ::  TN(NPTS) , O(NPTS) , O2(NPTS) , N2(NPTS) , GLAt(NPTS) , &
  PZ(NPTS) , GLOnd(NPTS)

  INTEGER ::   ilon1_3d_fixed_ht(npts,nmp),ilon2_3d_fixed_ht(npts,nmp)
  INTEGER ::   ilat1_3d_fixed_ht(npts,nmp),ilat2_3d_fixed_ht(npts,nmp)
  INTEGER ::   ispecial_3d_fixed_ht(npts,nmp)
  INTEGER ::   ihl_3d_fixed_ht(npts,nmp),ihu_3d_fixed_ht(npts,nmp)
  REAL(kind=8) pz_1000(npts)
  REAL(kind=8) geo_grid_longitudes_degrees(N_longitudes)
  REAL(kind=8) geo_grid_latitudes_degrees(N_Latitudes)
  REAL(kind=8) small_power,small_number
  small_power = -20.
  small_number = 1.d-20

  sw_External_model_provides_NO_N4S_densities = GIP_switches(5) 

!   do l = 1 , N_longitudes
!     geo_grid_longitudes_degrees(l) = (float(l-1))*18.
!   enddo

!   do m = 1 , N_latitudes
!     geo_grid_latitudes_degrees(m) = (float(m - 46))*2.
!   enddo

  !write(6,*) '************************************'
  !write(6,*) 'longs ',geo_grid_longitudes_degrees
  !write(6,*) '************************************'
  !write(6,*) 'lats ',geo_grid_latitudes_degrees
  !write(6,*) '************************************'

!g
  iwrite1 = 0
  iwrite=0
  istop = 0
  if (istop == 1) stop
  do i = 1, N_heights
!    fixed_hts_in_km(i) = (float(i-1)*5.) + 90.
    fixed_hts_in_km(i) = fixed_heights_km(i)
    pzh(i) = fixed_hts_in_km(i)
  enddo
!g
!g  Big loop over all flux tubes (nmp and nlp).....
!g
!       write(6,*) '***** calling interface to plasma ***** '
  do mp = 1 , nmp
      !write(6,*) 'interface to plasma fixed_ht',mp,sw_1st_call_int_fixed_ht
      do lp = 1 , nlp
      !write(6,*) '     lp ',lp
      !g
      !g  calculate the 1D geographic inputs......
      !g
          do i=in(mp,lp), is(mp,lp)
              glat(i) = glat_3d(i,mp)
              glond(i) = glond_3d(i,mp)
              pz(i) = pz_3d(i,mp)
              pz_1000(i) = pz(i)*1000.
          enddo
      !g

      ! loop over all points on the tube...

!print *,'INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE : IN = ',IN
!print *,'INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE : IS = ',IS
!STOP


          DO 900 i = IN(mp,lp) , IS(mp,lp)
              !write(6,*) 'points ',i
              glond2 = GLOnd(i)
              IF ( glond2 > 360. ) glond2 = glond2 - 360.
              IF ( glond2 < 0. ) glond2 = glond2 + 360.

              if(sw_1st_call_int_fixed_ht) then

                  ispecial = 0
                  DO 200 ilon = 1 , N_longitudes
                      IF ( geo_grid_longitudes_degrees(ilon) > glond2 ) THEN
                          ilon2 = ilon
                          ilon1 = ilon - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , glond2 , &
                          geo_grid_longitudes_degrees(ilon2) , geo_grid_longitudes_degrees(ilon1)
                          99001 FORMAT ('this 1 ',i3,3(2x,f10.1))
                          GOTO 250
                      ENDIF
                  200 ENDDO
                  ilon2 = 1
                  ilon1 = N_longitudes
                  ispecial = 1

                  250 DO 300 ilat = 1 , N_Latitudes
                      IF ( geo_grid_latitudes_degrees(ilat) > GLAt(i) ) THEN
                          ilat2 = ilat
                          ilat1 = ilat - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , GLAt(i) , &
                          geo_grid_latitudes_degrees(ilat2) , geo_grid_latitudes_degrees(ilat1)
                          GOTO 350
                      ENDIF
                  300 ENDDO
                  350 CONTINUE
              !g
              !g thus the required field line point lies
              !g within the square with latitudes ilat2,ilat1
              !g and longitudes ilon2,ilon1....at these four
              !g points there are two heights which are above
              !g and below the point.
              !g

                  DO 400 ih = 1 , N_heights
                      pzh(ih) = fixed_hts_in_km(ih)
                      IF ( pzh(ih) > PZ(i) ) THEN
                          ihu = ih
                          ihl = ih - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , PZ(i) , &
                          pzh(ihu) , pzh(ihl)
                          GOTO 450
                      ENDIF
                  400 ENDDO

                  ihu = N_heights
                  ihl = N_heights - 1
                  IF ( IWRite == 1 ) WRITE(6,99001) i , PZ(i) , pzh(ihu) , &
                  pzh(ihl)
                  450 IF ( ihu == 1 ) THEN
                      WRITE(6,*) 'ihu ' , ihu , pzh(ihu) , pzh(ihl)
                      DO 460 ii = IN(mp,lp) , IS(mp,lp)
                          WRITE(6,*) ii , PZ(ii)
                      460 ENDDO
                  ENDIF

              ! thus point lies in a box surrounded by the 8 points...

              ! ilon1 ilat1 ihl
              ! ilon1 ilat1 ihu
              ! ilon1 ilat2 ihl
              ! ilon1 ilat2 ihu
              ! ilon2 ilat1 ihl
              ! ilon2 ilat1 ihu
              ! ilon2 ilat2 ihl
              ! ilon2 ilat2 ihu

                  ilon1_3d_fixed_ht(i,mp) = ilon1
                  ilon2_3d_fixed_ht(i,mp) = ilon2
                  ilat1_3d_fixed_ht(i,mp) = ilat1
                  ilat2_3d_fixed_ht(i,mp) = ilat2
                  ispecial_3d_fixed_ht(i,mp) = ispecial
                  ihl_3d_fixed_ht(i,mp) = ihl
                  ihu_3d_fixed_ht(i,mp) = ihu

              else

                  ilon1 = ilon1_3d_fixed_ht(i,mp)
                  ilon2 = ilon2_3d_fixed_ht(i,mp)
                  ilat1 = ilat1_3d_fixed_ht(i,mp)
                  ilat2 = ilat2_3d_fixed_ht(i,mp)
                  ispecial = ispecial_3d_fixed_ht(i,mp)
                  ihl = ihl_3d_fixed_ht(i,mp)
                  ihu = ihu_3d_fixed_ht(i,mp)
!                   write(6,*) 'arggh ',ilon1,ilon2,ilat1,ilat2,ispecial,ihl,ihu

              endif

          ! neutral temperature on the eight surrounding points......

              tnu11 = TTS(ihu,ilat1,ilon1)
              tnl11 = TTS(ihl,ilat1,ilon1)
              tnu12 = TTS(ihu,ilat1,ilon2)
              tnl12 = TTS(ihl,ilat1,ilon2)
              tnu21 = TTS(ihu,ilat2,ilon1)
              tnl21 = TTS(ihl,ilat2,ilon1)
              tnu22 = TTS(ihu,ilat2,ilon2)
              tnl22 = TTS(ihl,ilat2,ilon2)

          ! electron temperature on the eight surrounding points......

!             teu11 = Telec(ihu,ilat1,ilon1)
!             tel11 = Telec(ihl,ilat1,ilon1)
!             teu12 = Telec(ihu,ilat1,ilon2)
!             tel12 = Telec(ihl,ilat1,ilon2)
!             teu21 = Telec(ihu,ilat2,ilon1)
!             tel21 = Telec(ihl,ilat2,ilon1)
!             teu22 = Telec(ihu,ilat2,ilon2)
!             tel22 = Telec(ihl,ilat2,ilon2)

          ! zonal wind on the eight surrounding points......

              uzu11 = VY(ihu,ilat1,ilon1)
              uzl11 = VY(ihl,ilat1,ilon1)
              uzu12 = VY(ihu,ilat1,ilon2)
              uzl12 = VY(ihl,ilat1,ilon2)
              uzu21 = VY(ihu,ilat2,ilon1)
              uzl21 = VY(ihl,ilat2,ilon1)
              uzu22 = VY(ihu,ilat2,ilon2)
              uzl22 = VY(ihl,ilat2,ilon2)

          ! meridional wind on the eight surrounding points......

              umu11 = VX(ihu,ilat1,ilon1)
              uml11 = VX(ihl,ilat1,ilon1)
              umu12 = VX(ihu,ilat1,ilon2)
              uml12 = VX(ihl,ilat1,ilon2)
              umu21 = VX(ihu,ilat2,ilon1)
              uml21 = VX(ihl,ilat2,ilon1)
              umu22 = VX(ihu,ilat2,ilon2)
              uml22 = VX(ihl,ilat2,ilon2)

          ! vertical wind on the eight surrounding points......

              uvu11 = WVZ(ihu,ilat1,ilon1)
              uvl11 = WVZ(ihl,ilat1,ilon1)
              uvu12 = WVZ(ihu,ilat1,ilon2)
              uvl12 = WVZ(ihl,ilat1,ilon2)
              uvu21 = WVZ(ihu,ilat2,ilon1)
              uvl21 = WVZ(ihl,ilat2,ilon1)
              uvu22 = WVZ(ihu,ilat2,ilon2)
              uvl22 = WVZ(ihl,ilat2,ilon2)

          ! atomic oxygen density on the eight surrounding points......

              ou11 = log10(O_density(ihu,ilat1,ilon1))
              ol11 = log10(O_density(ihl,ilat1,ilon1))
              ou12 = log10(O_density(ihu,ilat1,ilon2))
              ol12 = log10(O_density(ihl,ilat1,ilon2))
              ou21 = log10(O_density(ihu,ilat2,ilon1))
              ol21 = log10(O_density(ihl,ilat2,ilon1))
              ou22 = log10(O_density(ihu,ilat2,ilon2))
              ol22 = log10(O_density(ihl,ilat2,ilon2))

          ! molecular oxygen density on the eight surrounding points......

              oou11 = log10(O2_density(ihu,ilat1,ilon1))
              ool11 = log10(O2_density(ihl,ilat1,ilon1))
              oou12 = log10(O2_density(ihu,ilat1,ilon2))
              ool12 = log10(O2_density(ihl,ilat1,ilon2))
              oou21 = log10(O2_density(ihu,ilat2,ilon1))
              ool21 = log10(O2_density(ihl,ilat2,ilon1))
              oou22 = log10(O2_density(ihu,ilat2,ilon2))
              ool22 = log10(O2_density(ihl,ilat2,ilon2))

          ! molecular nitrogen density on the eight surrounding points......

              dnnu11 = log10(N2_density(ihu,ilat1,ilon1))
              dnnl11 = log10(N2_density(ihl,ilat1,ilon1))
              dnnu12 = log10(N2_density(ihu,ilat1,ilon2))
              dnnl12 = log10(N2_density(ihl,ilat1,ilon2))
              dnnu21 = log10(N2_density(ihu,ilat2,ilon1))
              dnnl21 = log10(N2_density(ihl,ilat2,ilon1))
              dnnu22 = log10(N2_density(ihu,ilat2,ilon2))
              dnnl22 = log10(N2_density(ihl,ilat2,ilon2))

          if(sw_External_model_provides_NO_N4S_densities) then
          ! NO density on the eight surrounding points......

              dNOu11 = log10(NO_density(ihu,ilat1,ilon1))
              dNOl11 = log10(NO_density(ihl,ilat1,ilon1))
              dNOu12 = log10(NO_density(ihu,ilat1,ilon2))
              dNOl12 = log10(NO_density(ihl,ilat1,ilon2))
              dNOu21 = log10(NO_density(ihu,ilat2,ilon1))
              dNOl21 = log10(NO_density(ihl,ilat2,ilon1))
              dNOu22 = log10(NO_density(ihu,ilat2,ilon2))
              dNOl22 = log10(NO_density(ihl,ilat2,ilon2))

          ! N4S density on the eight surrounding points......

              dN4Su11 = log10(N4S_density(ihu,ilat1,ilon1))
              dN4Sl11 = log10(N4S_density(ihl,ilat1,ilon1))
              dN4Su12 = log10(N4S_density(ihu,ilat1,ilon2))
              dN4Sl12 = log10(N4S_density(ihl,ilat1,ilon2))
              dN4Su21 = log10(N4S_density(ihu,ilat2,ilon1))
              dN4Sl21 = log10(N4S_density(ihl,ilat2,ilon1))
              dN4Su22 = log10(N4S_density(ihu,ilat2,ilon2))
              dN4Sl22 = log10(N4S_density(ihl,ilat2,ilon2))

          ! N2D density on the eight surrounding points......

              dN2Du11 = log10(N2D_density(ihu,ilat1,ilon1))
              dN2Dl11 = log10(N2D_density(ihl,ilat1,ilon1))
              dN2Du12 = log10(N2D_density(ihu,ilat1,ilon2))
              dN2Dl12 = log10(N2D_density(ihl,ilat1,ilon2))
              dN2Du21 = log10(N2D_density(ihu,ilat2,ilon1))
              dN2Dl21 = log10(N2D_density(ihl,ilat2,ilon1))
              dN2Du22 = log10(N2D_density(ihu,ilat2,ilon2))
              dN2Dl22 = log10(N2D_density(ihl,ilat2,ilon2))

          endif 

          ! now the 8 point interpolation.........

              fach = (PZ(i)-pzh(ihl))/(pzh(ihu)-pzh(ihl))
             if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'ZZZZZZZ ',ihl,ihu
!               write(6,*) 'XXXXXXX ',pz(i),pzh(ihl),pzh(ihu)
             endif

              IF ( fach > 1. ) THEN
                  tn11 = tnu11
                  uz11 = uzu11
                  um11 = umu11
                  uv11 = 0.
                  tn12 = tnu12
                  uz12 = uzu12
                  um12 = umu12
                  uv12 = 0.
                  tn21 = tnu21
                  uz21 = uzu21
                  um21 = umu21
                  uv21 = 0.
                  tn22 = tnu22
                  uz22 = uzu22
                  um22 = umu22
                  uv22 = 0.
!                 te11 = teu11+pz(i)-pzh(ihu)
!                 te12 = teu12+pz(i)-pzh(ihu)
!                 te21 = teu21+pz(i)-pzh(ihu)
!                 te22 = teu22+pz(i)-pzh(ihu)
              ELSE
                  tn11 = ((tnu11-tnl11)*fach) + tnl11
                  uz11 = ((uzu11-uzl11)*fach) + uzl11
                  um11 = ((umu11-uml11)*fach) + uml11
                  uv11 = ((uvu11-uvl11)*fach) + uvl11
                  tn12 = ((tnu12-tnl12)*fach) + tnl12
                  uz12 = ((uzu12-uzl12)*fach) + uzl12
                  um12 = ((umu12-uml12)*fach) + uml12
                  uv12 = ((uvu12-uvl12)*fach) + uvl12
                  tn21 = ((tnu21-tnl21)*fach) + tnl21
                  uz21 = ((uzu21-uzl21)*fach) + uzl21
                  um21 = ((umu21-uml21)*fach) + uml21
                  uv21 = ((uvu21-uvl21)*fach) + uvl21
                  tn22 = ((tnu22-tnl22)*fach) + tnl22
                  uz22 = ((uzu22-uzl22)*fach) + uzl22
                  um22 = ((umu22-uml22)*fach) + uml22
                  uv22 = ((uvu22-uvl22)*fach) + uvl22
!                 te11 = ((teu11-tel11)*fach) + tel11
!                 te12 = ((teu12-tel12)*fach) + tel12
!                 te21 = ((teu21-tel21)*fach) + tel21
!                 te22 = ((teu22-tel22)*fach) + tel22
              ENDIF

              o11 = (((ou11-ol11)*fach)+ol11)
              o12 = (((ou12-ol12)*fach)+ol12)
              o21 = (((ou21-ol21)*fach)+ol21)
              o22 = (((ou22-ol22)*fach)+ol22)
              oo11 =(((oou11-ool11)*fach)+ool11)
              oo12 = (((oou12-ool12)*fach)+ool12)
              oo21 = (((oou21-ool21)*fach)+ool21)
              oo22 = (((oou22-ool22)*fach)+ool22)
              dnn11 =(((dnnu11-dnnl11)*fach)+dnnl11)
              dnn12 = (((dnnu12-dnnl12)*fach)+dnnl12)
              dnn21 = (((dnnu21-dnnl21)*fach)+dnnl21)
              dnn22 = (((dnnu22-dnnl22)*fach)+dnnl22)

              top11 = (((topu11-topl11)*fach)+topl11)
              top12 = (((topu12-topl12)*fach)+topl12)
              top21 = (((topu21-topl21)*fach)+topl21)
              top22 = (((topu22-topl22)*fach)+topl22)
              tnop11 = (((tnopu11-tnopl11)*fach)+tnopl11)
              tnop12 = (((tnopu12-tnopl12)*fach)+tnopl12)
              tnop21 = (((tnopu21-tnopl21)*fach)+tnopl21)
              tnop22 = (((tnopu22-tnopl22)*fach)+tnopl22)
              to2p11 = (((to2pu11-to2pl11)*fach)+to2pl11)
              to2p12 = (((to2pu12-to2pl12)*fach)+to2pl12)
              to2p21 = (((to2pu21-to2pl21)*fach)+to2pl21)
              to2p22 = (((to2pu22-to2pl22)*fach)+to2pl22)

          if(sw_External_model_provides_NO_N4S_densities) then
              dNO11 =(((dNOu11-dNOl11)*fach)+dNOl11)
              dNO12 = (((dNOu12-dNOl12)*fach)+dNOl12)
              dNO21 = (((dNOu21-dNOl21)*fach)+dNOl21)
              dNO22 = (((dNOu22-dNOl22)*fach)+dNOl22)
              dN4S11 =(((dN4Su11-dN4Sl11)*fach)+dN4Sl11)
              dN4S12 = (((dN4Su12-dN4Sl12)*fach)+dN4Sl12)
              dN4S21 = (((dN4Su21-dN4Sl21)*fach)+dN4Sl21)
              dN4S22 = (((dN4Su22-dN4Sl22)*fach)+dN4Sl22)
              dN2D11 =(((dN2Du11-dN2Dl11)*fach)+dN2Dl11)
              dN2D12 = (((dN2Du12-dN2Dl12)*fach)+dN2Dl12)
              dN2D21 = (((dN2Du21-dN2Dl21)*fach)+dN2Dl21)
              dN2D22 = (((dN2Du22-dN2Dl22)*fach)+dN2Dl22)
          endif
          !g
!              if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'YYYYYY',ou11,ol11,fach
!               write(6,*) 'yabs 1',o11,fach
!              endif
              if(o11 > small_power) then
                  o11=10**o11
              else
                  o11=small_number
              endif
              if(o12 > small_power) then
                  o12=10**o12
              else
                  o12=small_number
              endif
              if(o21 > small_power) then
                  o21=10**o21
              else
                  o21=small_number
              endif
              if(o22 > small_power) then
                  o22=10**o22
              else
                  o22=small_number
              endif
!              if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'yabs ',o11
!              endif
          !g
              if(oo11 > small_power) then
                  oo11=10**oo11
              else
                  oo11=small_number
              endif
              if(oo12 > small_power) then
                  oo12=10**oo12
              else
                  oo12=small_number
              endif
              if(oo21 > small_power) then
                  oo21=10**oo21
              else
                  oo21=small_number
              endif
              if(oo22 > small_power) then
                  oo22=10**oo22
              else
                  oo22=small_number
              endif
          !g
              if(dnn11 > small_power) then
                  dnn11=10**dnn11
              else
                  dnn11=small_number
              endif
              if(dnn12 > small_power) then
                  dnn12=10**dnn12
              else
                  dnn12=small_number
              endif
              if(dnn21 > small_power) then
                  dnn21=10**dnn21
              else
                  dnn21=small_number
              endif
              if(dnn22 > small_power) then
                  dnn22=10**dnn22
              else
                  dnn22=small_number
              endif
          !g
          if(sw_External_model_provides_NO_N4S_densities) then
              if(dNO11 > small_power) then
                  dNO11=10**dNO11
              else
                  dNO11=small_number
              endif
              if(dNO12 > small_power) then
                  dNO12=10**dNO12
              else
                  dNO12=small_number
              endif
              if(dNO21 > small_power) then
                  dNO21=10**dNO21
              else
                  dNO21=small_number
              endif
              if(dNO22 > small_power) then
                  dNO22=10**dNO22
              else
                  dNO22=small_number
              endif
          !g
              if(dN4S11 > small_power) then
                  dN4S11=10**dN4S11
              else
                  dN4S11=small_number
              endif
              if(dN4S12 > small_power) then
                  dN4S12=10**dN4S12
              else
                  dN4S12=small_number
              endif
              if(dN4S21 > small_power) then
                  dN4S21=10**dN4S21
              else
                  dN4S21=small_number
              endif
              if(dN4S22 > small_power) then
                  dN4S22=10**dN4S22
              else
                  dN4S22=small_number
              endif
          !g
              if(dN2D11 > small_power) then
                  dN2D11=10**dN2D11
              else
                  dN2D11=small_number
              endif
              if(dN2D12 > small_power) then
                  dN2D12=10**dN2D12
              else
                  dN2D12=small_number
              endif
              if(dN2D21 > small_power) then
                  dN2D21=10**dN2D21
              else
                  dN2D21=small_number
              endif
              if(dN2D22 > small_power) then
                  dN2D22=10**dN2D22
              else
                  dN2D22=small_number
              endif

            endif
          !g

              IF ( ispecial == 0 ) THEN
                  faclon = (glond2-geo_grid_longitudes_degrees(ilon1)) / &
                           (geo_grid_longitudes_degrees(ilon2)-geo_grid_longitudes_degrees(ilon1))
              ELSEIF ( ispecial == 1 ) THEN
                  faclon = (glond2-geo_grid_longitudes_degrees(ilon1)) &
                  /(geo_grid_longitudes_degrees(ilon2)+360.-geo_grid_longitudes_degrees(ilon1))
              ENDIF
          ! write(6,*) 'faclon  ',faclon
              tn2 = ((tn22-tn21)*faclon) + tn21
              tn1 = ((tn12-tn11)*faclon) + tn11
!             te2 = ((te22-te21)*faclon) + te21
!             te1 = ((te12-te11)*faclon) + te11
              uz2 = ((uz22-uz21)*faclon) + uz21
              uz1 = ((uz12-uz11)*faclon) + uz11
              um2 = ((um22-um21)*faclon) + um21
              um1 = ((um12-um11)*faclon) + um11
              uv2 = ((uv22-uv21)*faclon) + uv21
              uv1 = ((uv12-uv11)*faclon) + uv11
              do2 = ((o22-o21)*faclon) + o21
              do1 = ((o12-o11)*faclon) + o11
              doo2 = ((oo22-oo21)*faclon) + oo21
              doo1 = ((oo12-oo11)*faclon) + oo11
              dnn2 = ((dnn22-dnn21)*faclon) + dnn21
              dnn1 = ((dnn12-dnn11)*faclon) + dnn11

              top2 = ((top22-top21)*faclon) + top21
              top1 = ((top12-top11)*faclon) + top11
              tnop2 = ((tnop22-tnop21)*faclon) + tnop21
              tnop1 = ((tnop12-tnop11)*faclon) + tnop11
              to2p2 = ((to2p22-to2p21)*faclon) + to2p21
              to2p1 = ((to2p12-to2p11)*faclon) + to2p11

          if(sw_External_model_provides_NO_N4S_densities) then
              dNO2 = ((dNO22-dNO21)*faclon) + dNO21
              dNO1 = ((dNO12-dNO11)*faclon) + dNO11
              dN4S2 = ((dN4S22-dN4S21)*faclon) + dN4S21
              dN4S1 = ((dN4S12-dN4S11)*faclon) + dN4S11
              dN2D2 = ((dN2D22-dN2D21)*faclon) + dN2D21
              dN2D1 = ((dN2D12-dN2D11)*faclon) + dN2D11
          endif

              faclat = (GLAt(i)-geo_grid_latitudes_degrees(ilat1)) / &
                       (geo_grid_latitudes_degrees(ilat2)-geo_grid_latitudes_degrees(ilat1))
          ! write(6,*) 'faclat  ',faclat
              TN(i) = ((tn2-tn1)*faclat) + tn1
!             Te_dum(i) = ((te2-te1)*faclat) + te1
              !if(te_dum(i) > 5000.) then
              ! write(6,*) 'TE_DUM ',i,mp,lp,te_dum(i)
              !   te_dum(i) = 5000.
              !endif
              uz(i) = ((uz2-uz1)*faclat) + uz1
              um(i) = ((um2-um1)*faclat) + um1
              uv(i) = ((uv2-uv1)*faclat) + uv1
              O(i) = ((do2-do1)*faclat) + do1
              if(o(i) < small_number) o(i)=small_number
              O2(i) = ((doo2-doo1)*faclat) + doo1
              if(o2(i) < small_number) o2(i)=small_number
              N2(i) = ((dnn2-dnn1)*faclat) + dnn1
              if(n2(i) < small_number) n2(i)=small_number
          if(sw_External_model_provides_NO_N4S_densities) then
              NO(i) = ((dNO2-dNO1)*faclat) + dNO1
              if(NO(i) < small_number) NO(i)=small_number
              N4S(i) = ((dN4S2-dN4S1)*faclat) + dN4S1
              if(N4S(i) < small_number) N4S(i)=small_number
              N2D(i) = ((dN2D2-dN2D1)*faclat) + dN2D1
              if(N2D(i) < small_number) N2D(i)=small_number
          endif

          !g THe above does not let any of the densities get lower
          !g than 'small_number'.  This can be a problem at low solar activity
          !g at the top of the larger flux-tubes.......
          !g
          ! write(199,*) i,mp,lp,faclon,faclat,fach
          900 ENDDO

          do i= in(mp,lp), is(mp,lp)
              TN_plasma_input_3d(i,mp) = tn(i)
              TN_plasma_input_3d(i,mp) = tn(i)
              if( TN_plasma_input_3d(i,mp) .lt. 0.) then
                  write(6,*) 'TN_plasma_input_3d lt 0)',i,mp,lp,tn(i)
                  TN_plasma_input_3d(i,mp) = 50.
              endif
              O_plasma_input_3d(i,mp) = o(i)
              O2_plasma_input_3d(i,mp) = o2(i)
              N2_plasma_input_3d(i,mp) = n2(i)
          if(sw_External_model_provides_NO_N4S_densities) then
              NO_plasma_input_3d(i,mp) = NO(i)
              N4S_plasma_input_3d(i,mp) = N4S(i)
              N2D_plasma_input_3d(i,mp) = N2D(i)
          endif
              um_plasma_input_3d(i,mp) = um(i)
              uz_plasma_input_3d(i,mp) = uz(i)
              uv_plasma_input_3d(i,mp) = uv(i)
!             te_plasma_input_3d(i,mp) = te_dum(i)

          enddo
      !g
      !g  end of the big flux tubes loop....
      !g
      enddo
  enddo
!       write(6,*) '***** done with interface to plasma ***** '

  sw_1st_call_int_fixed_ht = .FALSE.
  RETURN



end SUBROUTINE INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE











SUBROUTINE INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
  ni,no_plus_3d,o2_plus_3d, &
  n2_plus_3d,n_plus_3d, &
  Te,Ti, &
  oplus_high_res_fixed,hplus_high_res_fixed, &
  noplus_high_res_fixed,o2plus_high_res_fixed, &
  n2plus_high_res_fixed,nplus_high_res_fixed, &
  Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed)

  IMPLICIT NONE
  REAL(kind=8) :: dtotinv  , &
  factor , d
  REAL(kind=8) :: &
  fac
  INTEGER :: i , iheight , &
  iup , ido , ih , ip
  INTEGER ::  l , m ,  mp , lp ,  &
  in1 , in2 , &
  istop , NmF2_index , maxpos(1)
  INTEGER :: max_ic

  parameter (max_ic = 2*nmp*nlp)

  REAL(kind=8) :: ni(npts,nmp,2)
  REAL(kind=8) :: no_plus_3d(npts,nmp)
  REAL(kind=8) :: o2_plus_3d(npts,nmp)
  REAL(kind=8) :: n2_plus_3d(npts,nmp)
  REAL(kind=8) :: n_plus_3d(npts,nmp)
  REAL(kind=8) :: Te(npts,nmp)
  REAL(kind=8) :: Ti(npts,nmp,2)

  REAL(kind=8) :: oplus_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: hplus_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: noplus_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: o2plus_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: n2plus_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: nplus_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: Te_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: Ti1_high_res_fixed(interface_hts,91,90)
  REAL(kind=8) :: Ti2_high_res_fixed(interface_hts,91,90)

  REAL(kind=8) :: oplus_interpolated(3)
  REAL(kind=8) :: hplus_interpolated(3)
  REAL(kind=8) :: noplus_interpolated(3)
  REAL(kind=8) :: o2plus_interpolated(3)
  REAL(kind=8) :: n2plus_interpolated(3)
  REAL(kind=8) :: nplus_interpolated(3)
  REAL(kind=8) :: Te_interpolated(3)
  REAL(kind=8) :: Ti1_interpolated(3)
  REAL(kind=8) :: Ti2_interpolated(3)
!g
!g  loop over all heights and the lats/longs......
!g
  do 100 iheight=1,interface_hts
      do 200 l=1,90
      ! write(6,*) 'l       ',l
          do 300 m=mlow_interface(l),mhigh_interface(l)
          ! write(6,*) 'm       ',m
          !g
          !g  Initialise our output parameters....
          !g
              oplus_high_res_fixed(iheight,m,l) = 0.0
              hplus_high_res_fixed(iheight,m,l) = 0.0
              noplus_high_res_fixed(iheight,m,l) = 0.0
              o2plus_high_res_fixed(iheight,m,l) = 0.0
              n2plus_high_res_fixed(iheight,m,l) = 0.0
              nplus_high_res_fixed(iheight,m,l) = 0.0
              Te_high_res_fixed(iheight,m,l) = 0.0
              Ti1_high_res_fixed(iheight,m,l) = 0.0
              Ti2_high_res_fixed(iheight,m,l) = 0.0
              dtotinv = 0.0
          !g
          !g  The number of contributing flux-tube points
          !g  is always 3....
          !g
              do 400 i=1,3
              !g
              !g  The first interpolation uses facfac
              !g
                  factor=facfac_interface(i,iheight,m,l)
                  mp=ii1_interface(i,iheight,m,l)
                  lp=ii2_interface(i,iheight,m,l)
                  in2=ii3_interface(i,iheight,m,l)
                  in1=ii4_interface(i,iheight,m,l)

                  oplus_interpolated(i) = ((Ni(in2,mp,1)-Ni(in1,mp,1))*factor) &
                  + Ni(in1,mp,1)

                  hplus_interpolated(i) = ((Ni(in2,mp,2)-Ni(in1,mp,2))*factor) &
                  + Ni(in1,mp,2)

                  noplus_interpolated(i) = ((no_plus_3d(in2,mp)-no_plus_3d(in1,mp))*factor) &
                  + no_plus_3d(in1,mp)

                  o2plus_interpolated(i) = ((o2_plus_3d(in2,mp)-o2_plus_3d(in1,mp))*factor) &
                  + o2_plus_3d(in1,mp)

                  n2plus_interpolated(i) = ((n2_plus_3d(in2,mp)-n2_plus_3d(in1,mp))*factor) &
                  + n2_plus_3d(in1,mp)

                  nplus_interpolated(i) = ((n_plus_3d(in2,mp)-n_plus_3d(in1,mp))*factor) &
                  + n_plus_3d(in1,mp)

                  Te_interpolated(i) = ((Te(in2,mp)-Te(in1,mp))*factor) &
                  + Te(in1,mp)

                  Ti1_interpolated(i) = ((Ti(in2,mp,1)-Ti(in1,mp,1))*factor) &
                  + Ti(in1,mp,1)

                  Ti2_interpolated(i) = ((Ti(in2,mp,2)-Ti(in1,mp,2))*factor) &
                  + Ti(in1,mp,2)

              !g
              !g Now we calculate the parameters at each point.....
              !g

                  d = dd_interface(i,iheight,m,l)
                  dtotinv = (1./d) + dtotinv

                  oplus_high_res_fixed(iheight,m,l) = oplus_interpolated(i)/d &
                  + oplus_high_res_fixed(iheight,m,l)

                  hplus_high_res_fixed(iheight,m,l) = hplus_interpolated(i)/d &
                  + hplus_high_res_fixed(iheight,m,l)

                  noplus_high_res_fixed(iheight,m,l) = noplus_interpolated(i)/d &
                  + noplus_high_res_fixed(iheight,m,l)

                  o2plus_high_res_fixed(iheight,m,l) = o2plus_interpolated(i)/d &
                  + o2plus_high_res_fixed(iheight,m,l)

                  n2plus_high_res_fixed(iheight,m,l) = n2plus_interpolated(i)/d &
                  + n2plus_high_res_fixed(iheight,m,l)

                  nplus_high_res_fixed(iheight,m,l) = nplus_interpolated(i)/d &
                  + nplus_high_res_fixed(iheight,m,l)

                  Te_high_res_fixed(iheight,m,l) = Te_interpolated(i)/d &
                  + Te_high_res_fixed(iheight,m,l)

                  Ti1_high_res_fixed(iheight,m,l) = Ti1_interpolated(i)/d &
                  + Ti1_high_res_fixed(iheight,m,l)

                  Ti2_high_res_fixed(iheight,m,l) = Ti2_interpolated(i)/d &
                  + Ti2_high_res_fixed(iheight,m,l)

              400 ENDDO
          !g
              oplus_high_res_fixed(iheight,m,l) = oplus_high_res_fixed(iheight,m,l)/dtotinv

              hplus_high_res_fixed(iheight,m,l) = hplus_high_res_fixed(iheight,m,l)/dtotinv

              noplus_high_res_fixed(iheight,m,l) = noplus_high_res_fixed(iheight,m,l)/dtotinv

              o2plus_high_res_fixed(iheight,m,l) = o2plus_high_res_fixed(iheight,m,l)/dtotinv

              n2plus_high_res_fixed(iheight,m,l) = n2plus_high_res_fixed(iheight,m,l)/dtotinv

              nplus_high_res_fixed(iheight,m,l) = nplus_high_res_fixed(iheight,m,l)/dtotinv

              Te_high_res_fixed(iheight,m,l) = Te_high_res_fixed(iheight,m,l)/dtotinv

              Ti1_high_res_fixed(iheight,m,l) = Ti1_high_res_fixed(iheight,m,l)/dtotinv

              Ti2_high_res_fixed(iheight,m,l) = Ti2_high_res_fixed(iheight,m,l)/dtotinv

          300 ENDDO
      200 ENDDO
  100 ENDDO


  return




end SUBROUTINE INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO










SUBROUTINE interface__init_for_high_lat_ions (nhgt,altitude_metres,iht_above,iht_below,factor_ht, &
                                          iht_above_inverse,iht_below_inverse,factor_ht_inverse)

IMPLICIT NONE

INTEGER iht , iht_interface , nhgt
REAL(kind=8) interface_altitude_metres(interface_hts)
REAL(kind=8) altitude_metres(90)
REAL(kind=8) factor_ht(90)
INTEGER iht_above(90)
INTEGER iht_below(90)
REAL(kind=8) factor_ht_inverse(interface_hts)
INTEGER iht_above_inverse(interface_hts)
INTEGER iht_below_inverse(interface_hts)

do iht_interface = 1 , interface_hts
!  interface_altitude_metres(iht_interface) = ((float(iht_interface - 1)*5.) + 90.) * 1000.
  interface_altitude_metres(iht_interface) = fixed_heights_km(iht_interface) * 1000.
enddo

! Loop over high_lat_ions heights....

do iht = 1 , NHGT

!write(6,*) 'iht ',iht
! height interpolation....

do iht_interface = 1 , interface_hts
 !write(6,*) 'iht_interface ',iht_interface,interface_altitude_metres(iht_interface) , altitude_metres(iht)
if ( interface_altitude_metres(iht_interface) > altitude_metres(iht) ) then 
  iht_above(iht) = iht_interface
  iht_below(iht) = iht_interface - 1
  goto 2500
endif
enddo 
iht_above(iht) = interface_hts
iht_below(iht) = interface_hts - 1
2500 continue
if ( iht_above(iht) == 1 ) then 
iht_above(iht) = 2
iht_below(iht) = 1
endif

factor_ht(iht) = (altitude_metres(iht) - interface_altitude_metres(iht_below(iht))) /  &
               (interface_altitude_metres(iht_above(iht)) - interface_altitude_metres(iht_below(iht)))

enddo

! now need to calculate the inverse interpolation,
! ie, from the high_lat levels to the fixed height interface (183) levels....

do iht_interface = 1 , interface_hts

!write(6,*) 'iht ',iht
! height interpolation....

do iht = 1 , nhgt
 !write(6,*) 'iht_interface ',iht_interface,interface_altitude_metres(iht_interface) , altitude_metres(iht)
if ( altitude_metres(iht) > interface_altitude_metres(iht_interface) ) then 
  iht_above_inverse(iht_interface) = iht
  iht_below_inverse(iht_interface) = iht - 1
  goto 3000
endif
enddo 
iht_above_inverse(iht_interface) = nhgt
iht_below_inverse(iht_interface) = nhgt - 1
3000 continue
if ( iht_above_inverse(iht_interface) == 1 ) then 
iht_above_inverse(iht_interface) = 2
iht_below_inverse(iht_interface) = 1
endif

factor_ht_inverse(iht_interface) = (interface_altitude_metres(iht_interface) &
                               - altitude_metres(iht_below_inverse(iht_interface))) /  &
               (altitude_metres(iht_above_inverse(iht_interface)) - altitude_metres(iht_below_inverse(iht_interface)))

enddo
return




end SUBROUTINE interface__init_for_high_lat_ions







SUBROUTINE INTERFACE__FIXED_GEO_to_high_lat_ions ( &
           o_density_fixed_ht,o2_density_fixed_ht,n2_density_fixed_ht, &
           NO_density_fixed_ht, &
           N4S_density_fixed_ht,N2D_density_fixed_ht, &
           Vx_fixed_ht,Vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
           qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
           qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
           o_density_high_lat_ions,o2_density_high_lat_ions,n2_density_high_lat_ions, &
           NO_density_high_lat_ions, &
           N4S_density_high_lat_ions,N2D_density_high_lat_ions, &
           Vx_high_lat_ions,Vy_high_lat_ions,wvz_high_lat_ions,tn_high_lat_ions,qion3d_high_lat_ions, &
           qo2p_aurora_high_lat_ions, qop_aurora_high_lat_ions, qn2p_aurora_high_lat_ions, &
           qnp_aurora_high_lat_ions, qtef_aurora_high_lat_ions, &
           factor_ht_array,iht_above_array,iht_below_array,nhgt, &
           GIP_switches)

IMPLICIT NONE

INTEGER ilon  , ilat , iht , nhgt
LOGICAL :: sw_External_model_provides_NO_N4S_densities
LOGICAL :: sw_input_Auroral_production_is_single_overall_rate
LOGICAL :: GIP_switches(20)

REAL(kind=8) o_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) o2_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) n2_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) NO_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) N4S_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) N2D_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) Vx_fixed_ht(interface_hts,91,20)
REAL(kind=8) Vy_fixed_ht(interface_hts,91,20)
REAL(kind=8) wvz_fixed_ht(interface_hts,91,20)
REAL(kind=8) tts_fixed_ht(interface_hts,91,20)
REAL(kind=8) qion3d_fixed_ht(interface_hts,91,20)
REAL(kind=8) qo2p_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qop_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qn2p_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qnp_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) qtef_aurora_fixed_ht(interface_hts,91,20)
REAL(kind=8) interface_height(interface_hts)

REAL(kind=8) o_density_high_lat_ions(90,91,20)
REAL(kind=8) o2_density_high_lat_ions(90,91,20)
REAL(kind=8) n2_density_high_lat_ions(90,91,20)
REAL(kind=8) NO_density_high_lat_ions(90,91,20)
REAL(kind=8) N4S_density_high_lat_ions(90,91,20)
REAL(kind=8) N2D_density_high_lat_ions(90,91,20)
REAL(kind=8) Vx_high_lat_ions(90,91,20)
REAL(kind=8) Vy_high_lat_ions(90,91,20)
REAL(kind=8) Vz_high_lat_ions(90,91,20)
REAL(kind=8) wvz_high_lat_ions(90,91,20)
REAL(kind=8) Tn_high_lat_ions(90,91,20)
REAL(kind=8) qion3d_high_lat_ions(90,91,20)
REAL(kind=8) qo2p_aurora_high_lat_ions(90,91,20)
REAL(kind=8) qop_aurora_high_lat_ions(90,91,20)
REAL(kind=8) qn2p_aurora_high_lat_ions(90,91,20)
REAL(kind=8) qnp_aurora_high_lat_ions(90,91,20)
REAL(kind=8) qtef_aurora_high_lat_ions(90,91,20)
REAL(kind=8) zkm(90)
REAL(kind=8) O_den_u_log
REAL(kind=8) O_den_l_log
REAL(kind=8) O_den_log
REAL(kind=8) O_den 
REAL(kind=8) O2_den_u_log
REAL(kind=8) O2_den_l_log
REAL(kind=8) O2_den_log
REAL(kind=8) O2_den
REAL(kind=8) N2_den_u_log
REAL(kind=8) N2_den_l_log
REAL(kind=8) N2_den_log
REAL(kind=8) N2_den
REAL(kind=8) HYD_den_u_log
REAL(kind=8) HYD_den_l_log
REAL(kind=8) HYD_den_log
REAL(kind=8) HYD_den
REAL(kind=8) HEL_den_u_log
REAL(kind=8) HEL_den_l_log
REAL(kind=8) HEL_den_log
REAL(kind=8) HEL_den
REAL(kind=8) NO_den_u_log
REAL(kind=8) NO_den_l_log
REAL(kind=8) NO_den_log
REAL(kind=8) NO_den
REAL(kind=8) N4S_den_u_log
REAL(kind=8) N4S_den_l_log
REAL(kind=8) N4S_den_log
REAL(kind=8) N4S_den
REAL(kind=8) N2D_den_u_log
REAL(kind=8) N2D_den_l_log
REAL(kind=8) N2D_den_log
REAL(kind=8) N2D_den
REAL(kind=8) Tn_u
REAL(kind=8) Tn_l
REAL(kind=8) Tn
REAL(kind=8) Vx_u
REAL(kind=8) Vx_l
REAL(kind=8) Vx
REAL(kind=8) Vy_u
REAL(kind=8) Vy_l
REAL(kind=8) Vy
REAL(kind=8) Vz_u
REAL(kind=8) Vz_l
REAL(kind=8) Vz
REAL(kind=8) qion3d_u
REAL(kind=8) qion3d_l
REAL(kind=8) qion3d
REAL(kind=8) qo2p_aurora_u
REAL(kind=8) qo2p_aurora_l
REAL(kind=8) qo2p_aurora
REAL(kind=8) qop_aurora_u
REAL(kind=8) qop_aurora_l
REAL(kind=8) qop_aurora
REAL(kind=8) qn2p_aurora_u
REAL(kind=8) qn2p_aurora_l
REAL(kind=8) qn2p_aurora
REAL(kind=8) qnp_aurora_u
REAL(kind=8) qnp_aurora_l
REAL(kind=8) qnp_aurora
REAL(kind=8) qtef_aurora_u
REAL(kind=8) qtef_aurora_l
REAL(kind=8) qtef_aurora

REAL(kind=8) factor_ht_array(90,91,20)
INTEGER iht_above_array(90,91,20)
INTEGER iht_below_array(90,91,20)
REAL(kind=8) factor_ht
REAL(kind=8) factor_ht_for_tn
INTEGER iht_above
INTEGER iht_below

sw_External_model_provides_NO_N4S_densities = GIP_switches(5)
sw_input_Auroral_production_is_single_overall_rate = GIP_switches(7)

! Loop over high_lat_ions longitudes....

do ilon = 1 , 20

 !write(6,*) ' the ilon ******* ' ,ilon

! Loop over high_lat_ions latitudes....

do ilat = 1 , 91

 !write(6,*) ' the ilat ******* ' ,ilat

 ! only latitudes up to 40 or after 50...

if (ilat < 41 .or. ilat > 49) then

! Loop over high_lat_ions heights....

do iht = 1 , NHGT

 !write(6,*) ' the iht ******* ' ,iht

! height interpolation....

factor_ht = factor_ht_array(iht,ilat,ilon)
iht_above = iht_above_array(iht,ilat,ilon)
iht_below = iht_below_array(iht,ilat,ilon)

!write(6,*) '***** interpy  ',ilon,ilat,iht,iht_above,iht_below,factor_ht


   O_den_u_log=log10(O_density_fixed_ht(iht_above,ilat,ilon))
   O_den_l_log=log10(O_density_fixed_ht(iht_below,ilat,ilon))
   O_den_log = ((O_den_u_log - O_den_l_log)*factor_ht) + O_den_l_log
   O_den = 10**(O_den_log)

!    write(6,*) 'fgfg ',O2_density_fixed_ht(iht_above,ilat,ilon),iht_above,ilat,ilon
   O2_den_u_log=log10(O2_density_fixed_ht(iht_above,ilat,ilon))
   O2_den_l_log=log10(O2_density_fixed_ht(iht_below,ilat,ilon))
   O2_den_log = ((O2_den_u_log - O2_den_l_log)*factor_ht) + O2_den_l_log
   O2_den = 10**(O2_den_log)

   N2_den_u_log=log10(N2_density_fixed_ht(iht_above,ilat,ilon))
   N2_den_l_log=log10(N2_density_fixed_ht(iht_below,ilat,ilon))
   N2_den_log = ((N2_den_u_log - N2_den_l_log)*factor_ht) + N2_den_l_log
   N2_den = 10**(N2_den_log)

   if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_u_log=log10(NO_density_fixed_ht(iht_above,ilat,ilon))
   NO_den_l_log=log10(NO_density_fixed_ht(iht_below,ilat,ilon))
   NO_den_log = ((NO_den_u_log - NO_den_l_log)*factor_ht) + NO_den_l_log
   NO_den = 10**(NO_den_log)

   N4S_den_u_log=log10(N4S_density_fixed_ht(iht_above,ilat,ilon))
   N4S_den_l_log=log10(N4S_density_fixed_ht(iht_below,ilat,ilon))
   N4S_den_log = ((N4S_den_u_log - N4S_den_l_log)*factor_ht) + N4S_den_l_log
   N4S_den = 10**(N4S_den_log)

   N2D_den_u_log=log10(N2D_density_fixed_ht(iht_above,ilat,ilon))
   N2D_den_l_log=log10(N2D_density_fixed_ht(iht_below,ilat,ilon))
   N2D_den_log = ((N2D_den_u_log - N2D_den_l_log)*factor_ht) + N2D_den_l_log
   N2D_den = 10**(N2D_den_log)
   endif

   factor_ht_for_tn = factor_ht
   if ( factor_ht_for_tn .gt. 1.0) factor_ht_for_tn = 1.0

   Tn_u=tts_fixed_ht(iht_above,ilat,ilon)
   Tn_l=tts_fixed_ht(iht_below,ilat,ilon)
   Tn = ((Tn_u - Tn_l)*factor_ht_for_tn) + Tn_l

   Vx_u=Vx_fixed_ht(iht_above,ilat,ilon)
   Vx_l=Vx_fixed_ht(iht_below,ilat,ilon)
   Vx = ((Vx_u - Vx_l)*factor_ht) + Vx_l

   Vy_u=Vy_fixed_ht(iht_above,ilat,ilon)
   Vy_l=Vy_fixed_ht(iht_below,ilat,ilon)
   Vy = ((Vy_u - Vy_l)*factor_ht) + Vy_l

   Vz_u=wvz_fixed_ht(iht_above,ilat,ilon)
   Vz_l=wvz_fixed_ht(iht_below,ilat,ilon)
   Vz = ((Vz_u - Vz_l)*factor_ht) + Vz_l

if (sw_input_Auroral_production_is_single_overall_rate) then

   qion3d_u=qion3d_fixed_ht(iht_above,ilat,ilon)
   qion3d_l=qion3d_fixed_ht(iht_below,ilat,ilon)
   qion3d = ((qion3d_u - qion3d_l)*factor_ht) + qion3d_l

else

   qo2p_aurora_u=qo2p_aurora_fixed_ht(iht_above,ilat,ilon)
   qo2p_aurora_l=qo2p_aurora_fixed_ht(iht_below,ilat,ilon)
   qo2p_aurora = ((qo2p_aurora_u - qo2p_aurora_l)*factor_ht) + qo2p_aurora_l

   qop_aurora_u=qop_aurora_fixed_ht(iht_above,ilat,ilon)
   qop_aurora_l=qop_aurora_fixed_ht(iht_below,ilat,ilon)
   qop_aurora = ((qop_aurora_u - qop_aurora_l)*factor_ht) + qop_aurora_l

   qn2p_aurora_u=qn2p_aurora_fixed_ht(iht_above,ilat,ilon)
   qn2p_aurora_l=qn2p_aurora_fixed_ht(iht_below,ilat,ilon)
   qn2p_aurora = ((qn2p_aurora_u - qn2p_aurora_l)*factor_ht) + qn2p_aurora_l

   qnp_aurora_u=qnp_aurora_fixed_ht(iht_above,ilat,ilon)
   qnp_aurora_l=qnp_aurora_fixed_ht(iht_below,ilat,ilon)
   qnp_aurora = ((qnp_aurora_u - qnp_aurora_l)*factor_ht) + qnp_aurora_l

   qtef_aurora_u=qtef_aurora_fixed_ht(iht_above,ilat,ilon)
   qtef_aurora_l=qtef_aurora_fixed_ht(iht_below,ilat,ilon)
   qtef_aurora = ((qtef_aurora_u - qtef_aurora_l)*factor_ht) + qtef_aurora_l

endif

   o_density_high_lat_ions(iht,ilat,ilon) = O_den
   o2_density_high_lat_ions(iht,ilat,ilon) = O2_den
   n2_density_high_lat_ions(iht,ilat,ilon) = N2_den
   if (sw_External_model_provides_NO_N4S_densities) then 
   NO_density_high_lat_ions(iht,ilat,ilon) = NO_den
   N4S_density_high_lat_ions(iht,ilat,ilon) = N4S_den
   N2D_density_high_lat_ions(iht,ilat,ilon) = N2D_den
   endif
   Vx_high_lat_ions(iht,ilat,ilon) = Vx
   Vy_high_lat_ions(iht,ilat,ilon) = Vy
   wvz_high_lat_ions(iht,ilat,ilon) = Vz
   Tn_high_lat_ions(iht,ilat,ilon) = Tn

if (sw_input_Auroral_production_is_single_overall_rate) then
   qion3d_high_lat_ions(iht,ilat,ilon) = qion3d
else
   qo2p_aurora_high_lat_ions(iht,ilat,ilon) = qo2p_aurora
   qop_aurora_high_lat_ions(iht,ilat,ilon) = qop_aurora
   qn2p_aurora_high_lat_ions(iht,ilat,ilon) = qn2p_aurora
   qnp_aurora_high_lat_ions(iht,ilat,ilon) = qnp_aurora
   qtef_aurora_high_lat_ions(iht,ilat,ilon) = qtef_aurora
endif

enddo

endif

enddo
enddo

return




end SUBROUTINE INTERFACE__FIXED_GEO_to_high_lat_ions












SUBROUTINE INTERFACE__high_lat_ions_to_FIXED_GEO ( &
           Oplus_density_high_lat_ions,Hplus_density_high_lat_ions,NOplus_density_high_lat_ions, &
           O2plus_density_high_lat_ions,N2plus_density_high_lat_ions,Nplus_density_high_lat_ions, &
           Te_high_lat_ions,Ti1_high_lat_ions,Ti2_high_lat_ions, &
           Oplus_density_fixed_ht,Hplus_density_fixed_ht,NOplus_density_fixed_ht, &
           O2plus_density_fixed_ht,N2plus_density_fixed_ht,Nplus_density_fixed_ht, &
           Te_fixed_ht,Ti1_fixed_ht,Ti2_fixed_ht, &
           factor_ht_inverse_3d,iht_above_inverse_3d,iht_below_inverse_3d)

IMPLICIT NONE

INTEGER ilon  , ilat , iht_interface

REAL(kind=8) Oplus_density_high_lat_ions(90,91,20)
REAL(kind=8) Hplus_density_high_lat_ions(90,91,20)
REAL(kind=8) NOplus_density_high_lat_ions(90,91,20)
REAL(kind=8) O2plus_density_high_lat_ions(90,91,20)
REAL(kind=8) N2plus_density_high_lat_ions(90,91,20)
REAL(kind=8) Nplus_density_high_lat_ions(90,91,20)
REAL(kind=8) Te_high_lat_ions(90,91,20)
REAL(kind=8) Ti1_high_lat_ions(90,91,20)
REAL(kind=8) Ti2_high_lat_ions(90,91,20)

REAL(kind=8) Oplus_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) Hplus_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) NOplus_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) O2plus_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) N2plus_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) Nplus_density_fixed_ht(interface_hts,91,20)
REAL(kind=8) Te_fixed_ht(interface_hts,91,20)
REAL(kind=8) Ti1_fixed_ht(interface_hts,91,20)
REAL(kind=8) Ti2_fixed_ht(interface_hts,91,20)

REAL(kind=8) :: Oplus_den_u
REAL(kind=8) :: Oplus_den_l
REAL(kind=8) :: Oplus_den

REAL(kind=8) :: Hplus_den_u
REAL(kind=8) :: Hplus_den_l
REAL(kind=8) :: Hplus_den

REAL(kind=8) :: NOplus_den_u
REAL(kind=8) :: NOplus_den_l
REAL(kind=8) :: NOplus_den

REAL(kind=8) :: O2plus_den_u
REAL(kind=8) :: O2plus_den_l
REAL(kind=8) :: O2plus_den

REAL(kind=8) :: N2plus_den_u
REAL(kind=8) :: N2plus_den_l
REAL(kind=8) :: N2plus_den

REAL(kind=8) :: Nplus_den_u
REAL(kind=8) :: Nplus_den_l
REAL(kind=8) :: Nplus_den

REAL(kind=8) :: Te_u
REAL(kind=8) :: Te_l
REAL(kind=8) :: Te

REAL(kind=8) :: Ti1_u
REAL(kind=8) :: Ti1_l
REAL(kind=8) :: Ti1

REAL(kind=8) :: Ti2_u
REAL(kind=8) :: Ti2_l
REAL(kind=8) :: Ti2

REAL(kind=8) factor_ht_inverse_3d(interface_hts,91,20)
INTEGER iht_above_inverse_3d(interface_hts,91,20)
INTEGER iht_below_inverse_3d(interface_hts,91,20)
REAL(kind=8) factor_ht
INTEGER iht_above
INTEGER iht_below


 Oplus_density_fixed_ht(:,:,:) = 0.0
 Hplus_density_fixed_ht(:,:,:) = 0.0
 NOplus_density_fixed_ht(:,:,:) = 0.0
 O2plus_density_fixed_ht(:,:,:) = 0.0
 N2plus_density_fixed_ht(:,:,:) = 0.0
 Nplus_density_fixed_ht(:,:,:) = 0.0
 Te_fixed_ht(:,:,:) = 0.0
 Ti1_fixed_ht(:,:,:) = 0.0
 Ti2_fixed_ht(:,:,:) = 0.0

! Loop over interface longitudes....

do ilon = 1 , 20



! Loop over interface latitudes....

!do ilat = 1 , 91
do ilat = 3 , 89     ! ISSUE: need to check at some point
                     ! that this is correct - I assume so -
                     ! because the polar points (1,2,90 and 91) should 
                     ! all be filled in later (so we only need to loop from 3 to 89)

if ( ilat < 41 .or. ilat > 49 ) then

! Loop over interface heights....

do iht_interface = 1 , interface_hts

! height interpolation....

factor_ht = factor_ht_inverse_3d(iht_interface,ilat,ilon)

if(factor_ht < 0.0) then
! write(6,*) 'factor ht ', iht_interface , factor_ht
factor_ht = 0.0
endif

iht_above = iht_above_inverse_3d(iht_interface,ilat,ilon)
iht_below = iht_below_inverse_3d(iht_interface,ilat,ilon)

   Oplus_den_u=Oplus_density_high_lat_ions(iht_above,ilat,ilon)
   Oplus_den_l=Oplus_density_high_lat_ions(iht_below,ilat,ilon)
   Oplus_den = ((Oplus_den_u - Oplus_den_l)*factor_ht) + Oplus_den_l
   Oplus_density_fixed_ht(iht_interface,ilat,ilon) = Oplus_den

   Hplus_den_u=Hplus_density_high_lat_ions(iht_above,ilat,ilon)
   Hplus_den_l=Hplus_density_high_lat_ions(iht_below,ilat,ilon)
   Hplus_den = ((Hplus_den_u - Hplus_den_l)*factor_ht) + Hplus_den_l
   Hplus_density_fixed_ht(iht_interface,ilat,ilon) = Hplus_den

   NOplus_den_u=NOplus_density_high_lat_ions(iht_above,ilat,ilon)
   NOplus_den_l=NOplus_density_high_lat_ions(iht_below,ilat,ilon)
   NOplus_den = ((NOplus_den_u - NOplus_den_l)*factor_ht) + NOplus_den_l
   NOplus_density_fixed_ht(iht_interface,ilat,ilon) = NOplus_den

   O2plus_den_u=O2plus_density_high_lat_ions(iht_above,ilat,ilon)
   O2plus_den_l=O2plus_density_high_lat_ions(iht_below,ilat,ilon)
   O2plus_den = ((O2plus_den_u - O2plus_den_l)*factor_ht) + O2plus_den_l
   O2plus_density_fixed_ht(iht_interface,ilat,ilon) = O2plus_den

   N2plus_den_u=N2plus_density_high_lat_ions(iht_above,ilat,ilon)
   N2plus_den_l=N2plus_density_high_lat_ions(iht_below,ilat,ilon)
   N2plus_den = ((N2plus_den_u - N2plus_den_l)*factor_ht) + N2plus_den_l
   N2plus_density_fixed_ht(iht_interface,ilat,ilon) = N2plus_den

   Nplus_den_u=Nplus_density_high_lat_ions(iht_above,ilat,ilon)
   Nplus_den_l=Nplus_density_high_lat_ions(iht_below,ilat,ilon)
   Nplus_den = ((Nplus_den_u - Nplus_den_l)*factor_ht) + Nplus_den_l
   Nplus_density_fixed_ht(iht_interface,ilat,ilon) = Nplus_den

   Te_u=Te_high_lat_ions(iht_above,ilat,ilon)
   Te_l=Te_high_lat_ions(iht_below,ilat,ilon)
   Te = ((Te_u - Te_l)*factor_ht) + Te_l
   Te_fixed_ht(iht_interface,ilat,ilon) = Te

   Ti1_u=Ti1_high_lat_ions(iht_above,ilat,ilon)
   Ti1_l=Ti1_high_lat_ions(iht_below,ilat,ilon)
   Ti1 = ((Ti1_u - Ti1_l)*factor_ht) + Ti1_l
   Ti1_fixed_ht(iht_interface,ilat,ilon) = Ti1

   Ti2_u=Ti2_high_lat_ions(iht_above,ilat,ilon)
   Ti2_l=Ti2_high_lat_ions(iht_below,ilat,ilon)
   Ti2 = ((Ti2_u - Ti2_l)*factor_ht) + Ti2_l
   Ti2_fixed_ht(iht_interface,ilat,ilon) = Ti2

enddo

endif

enddo
enddo

return




end SUBROUTINE INTERFACE__high_lat_ions_to_FIXED_GEO











SUBROUTINE INTERFACE__COMBINE_MID_LAT_AND_POLAR_PARAMS ( &
                     Oplus_density_fixed_ht_plasma,Hplus_density_fixed_ht_plasma, &
                     NOplus_density_fixed_ht_plasma,O2plus_density_fixed_ht_plasma, &
                     N2plus_density_fixed_ht_plasma,Nplus_density_fixed_ht_plasma, &
                     Te_fixed_ht_plasma,Ti1_fixed_ht_plasma,Ti2_fixed_ht_plasma, &
                     Oplus_density_fixed_ht_polar,Hplus_density_fixed_ht_polar, &
                     NOplus_density_fixed_ht_polar,O2plus_density_fixed_ht_polar, &
                     N2plus_density_fixed_ht_polar,Nplus_density_fixed_ht_polar, &
                     Te_fixed_ht_polar,Ti1_fixed_ht_polar,Ti2_fixed_ht_polar, &
                     Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
                     NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
                     N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
                     Te_fixed_ht_com,Ti1_fixed_ht_com,Ti2_fixed_ht_com)


  IMPLICIT NONE
  INTEGER :: n , m , l
  INTEGER :: ilon , ilon20
  INTEGER :: ilon_east , ilon_west
  INTEGER :: ispecial
  REAL(kind=8) :: long_20(20)
  REAL(kind=8) :: long_90(90)
  REAL(kind=8) :: factor_lon
  REAL(kind=8) :: high_fac
  REAL(kind=8) :: plas_fac

  REAL(kind=8) :: Oplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Hplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: NOplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: O2plus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: N2plus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Nplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Te_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Ti1_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Ti2_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8) :: Oplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Hplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: NOplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: O2plus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: N2plus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Nplus_density_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Te_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Ti1_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Ti2_fixed_ht_plasma(interface_hts,91,90)
  REAL(kind=8) :: Oplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Hplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: NOplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: O2plus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: N2plus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Nplus_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Te_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Ti1_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Ti2_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Ne_density_fixed_ht_polar(interface_hts,91,20)
  REAL(kind=8) :: Oplus_den_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: Hplus_den_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: NOplus_den_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: O2plus_den_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: N2plus_den_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: Nplus_den_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: Te_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: Ti1_fixed_polar_upsamp(interface_hts,91,90)
  REAL(kind=8) :: Ti2_fixed_polar_upsamp(interface_hts,91,90)


  do ilon = 1 , 90
   long_90(ilon) = (float(ilon - 1)) * 4.
  enddo

  do ilon20 = 1 , 20
   long_20(ilon20) = (float(ilon20 - 1)) * 18.
  enddo



  ispecial = 0
  do ilon = 1 , 90

  do ilon20 = 1 , 20
    if ( long_20(ilon20) > long_90(ilon) ) then
        ilon_east = ilon20 
        ilon_west = ilon20  - 1
        goto 2500
    endif
  enddo
  ilon_east = 1
  ilon_west = 20
  ispecial = 1
2500 continue
  if ( ilon_east == 1 ) then
     ilon_east = 1 
     ilon_west = 20
     ispecial = 1
  endif

if ( ispecial == 0 ) then 
factor_lon = (long_90(ilon) - long_20(ilon_west)) / (long_20(ilon_east) - long_20(ilon_west))
else
factor_lon = (long_90(ilon) - long_20(ilon_west)) / (long_20(ilon_east) + 360. - long_20(ilon_west))
endif


if ( factor_lon < 0.0 .or. factor_lon > 1.0 ) then
 write(6,*) 'ARGHHHHH ', factor_lon, ilon_east , ilon_west
endif

Oplus_den_fixed_polar_upsamp(:,:,ilon) =  &
(( Oplus_density_fixed_ht_polar(:,:,ilon_east) - Oplus_density_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + Oplus_density_fixed_ht_polar(:,:,ilon_west)

Hplus_den_fixed_polar_upsamp(:,:,ilon) =  &
(( Hplus_density_fixed_ht_polar(:,:,ilon_east) - Hplus_density_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + Hplus_density_fixed_ht_polar(:,:,ilon_west)

NOplus_den_fixed_polar_upsamp(:,:,ilon) =  &
(( NOplus_density_fixed_ht_polar(:,:,ilon_east) - NOplus_density_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + NOplus_density_fixed_ht_polar(:,:,ilon_west)

O2plus_den_fixed_polar_upsamp(:,:,ilon) =  &
(( O2plus_density_fixed_ht_polar(:,:,ilon_east) - O2plus_density_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + O2plus_density_fixed_ht_polar(:,:,ilon_west)

N2plus_den_fixed_polar_upsamp(:,:,ilon) =  &
(( N2plus_density_fixed_ht_polar(:,:,ilon_east) - N2plus_density_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + N2plus_density_fixed_ht_polar(:,:,ilon_west)

Nplus_den_fixed_polar_upsamp(:,:,ilon) =  &
(( Nplus_density_fixed_ht_polar(:,:,ilon_east) - Nplus_density_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + Nplus_density_fixed_ht_polar(:,:,ilon_west)

Te_fixed_polar_upsamp(:,:,ilon) =  &
(( Te_fixed_ht_polar(:,:,ilon_east) - Te_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + Te_fixed_ht_polar(:,:,ilon_west)

Ti1_fixed_polar_upsamp(:,:,ilon) =  &
(( Ti1_fixed_ht_polar(:,:,ilon_east) - Ti1_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + Ti1_fixed_ht_polar(:,:,ilon_west)

Ti2_fixed_polar_upsamp(:,:,ilon) =  &
(( Ti2_fixed_ht_polar(:,:,ilon_east) - Ti2_fixed_ht_polar(:,:,ilon_west)) * &
                 factor_lon) + Ti2_fixed_ht_polar(:,:,ilon_west)

  enddo



  do l = 1 , 90

  do m = 1 , MLOw_interface(l) - 1
    Oplus_density_fixed_ht_com(:,m,l) = Oplus_den_fixed_polar_upsamp(:,m,l)
    Hplus_density_fixed_ht_com(:,m,l) = Hplus_den_fixed_polar_upsamp(:,m,l)
    NOplus_density_fixed_ht_com(:,m,l) = NOplus_den_fixed_polar_upsamp(:,m,l)
    O2plus_density_fixed_ht_com(:,m,l) = O2plus_den_fixed_polar_upsamp(:,m,l)
    N2plus_density_fixed_ht_com(:,m,l) = N2plus_den_fixed_polar_upsamp(:,m,l)
    Nplus_density_fixed_ht_com(:,m,l) = Nplus_den_fixed_polar_upsamp(:,m,l)
    Te_fixed_ht_com(:,m,l) = Te_fixed_polar_upsamp(:,m,l)
    Ti1_fixed_ht_com(:,m,l) = Ti1_fixed_polar_upsamp(:,m,l)
    Ti2_fixed_ht_com(:,m,l) = Ti2_fixed_polar_upsamp(:,m,l)
  enddo

  do m = Mlow_interface(l),  Mlow_interface(l) + 3
    plas_fac = float(m-Mlow_interface(l)+1) / 5.
    high_fac = float(Mlow_interface(l)-m+4) / 5.
    Oplus_density_fixed_ht_com(:,m,l) = (Oplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (Oplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    Hplus_density_fixed_ht_com(:,m,l) = (Hplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (Hplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    NOplus_density_fixed_ht_com(:,m,l) = (NOplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (NOplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    O2plus_density_fixed_ht_com(:,m,l) = (O2plus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (O2plus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    N2plus_density_fixed_ht_com(:,m,l) = (N2plus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (N2plus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    Nplus_density_fixed_ht_com(:,m,l) = (Nplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (Nplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    Te_fixed_ht_com(:,m,l) = (Te_fixed_polar_upsamp(:,m,l) * high_fac) + &
                             (Te_fixed_ht_plasma(:,m,l) * plas_fac)
    Ti1_fixed_ht_com(:,m,l) = (Ti1_fixed_polar_upsamp(:,m,l) * high_fac) + &
                             (Ti1_fixed_ht_plasma(:,m,l) * plas_fac)
    Ti2_fixed_ht_com(:,m,l) = (Ti2_fixed_polar_upsamp(:,m,l) * high_fac) + &
                             (Ti2_fixed_ht_plasma(:,m,l) * plas_fac)
  enddo

  do m = MLOw_interface(l) + 4 , MHIgh_interface(l) - 4
    Oplus_density_fixed_ht_com(:,m,l) = Oplus_density_fixed_ht_plasma(:,m,l)
    Hplus_density_fixed_ht_com(:,m,l) = Hplus_density_fixed_ht_plasma(:,m,l)
    NOplus_density_fixed_ht_com(:,m,l) = NOplus_density_fixed_ht_plasma(:,m,l)
    O2plus_density_fixed_ht_com(:,m,l) = O2plus_density_fixed_ht_plasma(:,m,l)
    N2plus_density_fixed_ht_com(:,m,l) = N2plus_density_fixed_ht_plasma(:,m,l)
    Nplus_density_fixed_ht_com(:,m,l) = Nplus_density_fixed_ht_plasma(:,m,l)
    Te_fixed_ht_com(:,m,l) = Te_fixed_ht_plasma(:,m,l)
    Ti1_fixed_ht_com(:,m,l) = Ti1_fixed_ht_plasma(:,m,l)
    Ti2_fixed_ht_com(:,m,l) = Ti2_fixed_ht_plasma(:,m,l)
  enddo

  do m = MHIgh_interface(l) - 3 , MHIgh_interface(l)
    plas_fac = float(MHIgh_interface(l)+1-m) / 5.
    high_fac = float(m-MHIgh_interface(l)+4) / 5.
    Oplus_density_fixed_ht_com(:,m,l) = (Oplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (Oplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    Hplus_density_fixed_ht_com(:,m,l) = (Hplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (Hplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    NOplus_density_fixed_ht_com(:,m,l) = (NOplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (NOplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    O2plus_density_fixed_ht_com(:,m,l) = (O2plus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (O2plus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    N2plus_density_fixed_ht_com(:,m,l) = (N2plus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (N2plus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    Nplus_density_fixed_ht_com(:,m,l) = (Nplus_den_fixed_polar_upsamp(:,m,l) * high_fac) + &
                                        (Nplus_density_fixed_ht_plasma(:,m,l) * plas_fac)
    Te_fixed_ht_com(:,m,l) = (Te_fixed_polar_upsamp(:,m,l) * high_fac) + &
                             (Te_fixed_ht_plasma(:,m,l) * plas_fac)
    Ti1_fixed_ht_com(:,m,l) = (Ti1_fixed_polar_upsamp(:,m,l) * high_fac) + &
                             (Ti1_fixed_ht_plasma(:,m,l) * plas_fac)
    Ti2_fixed_ht_com(:,m,l) = (Ti2_fixed_polar_upsamp(:,m,l) * high_fac) + &
                             (Ti2_fixed_ht_plasma(:,m,l) * plas_fac)
  enddo

  do m = MHIgh_interface(l) + 1 , 91
    Oplus_density_fixed_ht_com(:,m,l) = Oplus_den_fixed_polar_upsamp(:,m,l)
    Hplus_density_fixed_ht_com(:,m,l) = Hplus_den_fixed_polar_upsamp(:,m,l)
    NOplus_density_fixed_ht_com(:,m,l) = NOplus_den_fixed_polar_upsamp(:,m,l)
    O2plus_density_fixed_ht_com(:,m,l) = O2plus_den_fixed_polar_upsamp(:,m,l)
    N2plus_density_fixed_ht_com(:,m,l) = N2plus_den_fixed_polar_upsamp(:,m,l)
    Nplus_density_fixed_ht_com(:,m,l) = Nplus_den_fixed_polar_upsamp(:,m,l)
    Te_fixed_ht_com(:,m,l) = Te_fixed_polar_upsamp(:,m,l)
    Ti1_fixed_ht_com(:,m,l) = Ti1_fixed_polar_upsamp(:,m,l)
    Ti2_fixed_ht_com(:,m,l) = Ti2_fixed_polar_upsamp(:,m,l)
  enddo

  enddo



  return




end SUBROUTINE INTERFACE__COMBINE_MID_LAT_AND_POLAR_PARAMS












SUBROUTINE INTERFACE__calculate_Ne_from_ion_densities ( &
                     Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
                     NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
                     N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
                     Ne_density_fixed_ht_com)
  IMPLICIT NONE

  REAL(kind=8), intent(in) :: Oplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: Hplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: NOplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: O2plus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: N2plus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: Nplus_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(out) :: Ne_density_fixed_ht_com(interface_hts,91,90)


  Ne_density_fixed_ht_com(:,:,:) = Oplus_density_fixed_ht_com(:,:,:)  + &
                                   Hplus_density_fixed_ht_com(:,:,:)  + &
                                   NOplus_density_fixed_ht_com(:,:,:) + &
                                   O2plus_density_fixed_ht_com(:,:,:) + &
                                   N2plus_density_fixed_ht_com(:,:,:) + &
                                   Nplus_density_fixed_ht_com(:,:,:) 

  return


end SUBROUTINE INTERFACE__calculate_Ne_from_ion_densities





SUBROUTINE INTERFACE__calculate_NMF2_HMF2_TEC_etc ( &
                     Ne_density_fixed_ht_com, Te_fixed_ht_com, Ti1_fixed_ht_com, &
                     NmF2,HmF2_km,TEC,Te_400,Ti_400,this_Ne_profile)


  IMPLICIT NONE

  REAL(kind=8), intent(in) :: Ne_density_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: Te_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(in) :: Ti1_fixed_ht_com(interface_hts,91,90)
  REAL(kind=8), intent(out) :: NmF2(91,90)
  REAL(kind=8), intent(out) :: HmF2_km(91,90)
  REAL(kind=8), intent(out) :: TEC(91,90)
  REAL(kind=8), intent(out) :: Te_400(91,90)
  REAL(kind=8), intent(out) :: Ti_400(91,90)

  REAL(kind=8), intent(out) :: this_Ne_profile(interface_hts)

!   INTEGER, intent(out) :: NmF2_index_array(91,90)
  INTEGER :: NmF2_index_array(91,90)
  REAL(kind=8)  :: NE_1d(interface_hts)
  REAL(kind=8)  :: a , b , c 
  REAL(kind=8)  :: x1 , x2 , x3 , y1 , y2 , y3 , x 
  INTEGER :: m , l , iheight , max_position(1) , NmF2_index
  INTEGER :: iht, max_position2


  this_Ne_profile(:) = Ne_density_fixed_ht_com(:,66,31)

  DO 1250 l = 1 , 90
      DO 1220 m = 1 , 91

          NE_1d(:) = Ne_density_fixed_ht_com(:,m,l)


! loop downwards from interface_hts to 1
! find the first turning point - that's ya nmf2.....
!
          max_position2 = interface_hts - 3
!            do iht = 182,10,-1
! looping downwards from 103 to 23 (600km to 200km)
! ie, hmf2 must be 200 < hmf2 < 600
!
          do iht = interface_hts - 1,23,-1
!           if ((NE_1d(iht+2) .lt. NE_1d(iht+1)) .and. (NE_1d(iht) .lt. NE_1d(iht+1)) .and. (NE_1d(iht-1) .lt. NE_1d(iht))) then
           if ((NE_1d(iht) .lt. NE_1d(iht+1)) .and. (NE_1d(iht-1) .lt. NE_1d(iht))) then
            max_position2 = iht + 1
            goto 2500
           endif
          enddo
!           max_position = maxloc(NE_1d)

!           NmF2_index = max_position(1)

2500   continue

          NmF2_index = max_position2
          NmF2_index_array(m,l) = NmF2_index

!            if (NmF2_index < 183 .and. NmF2_index > 1) then
          if (NmF2_index < 103 .and. NmF2_index > 23) then

!         x2 = FLOAT(NmF2_index-1)*5. + 90.
          x2 = fixed_heights_km(NmF2_index)
          x1 = x2 - 5.
          x3 = x2 + 5.
          y2 = NE_1d(NmF2_index)
          y1 = NE_1d(NmF2_index-1)
          y3 = NE_1d(NmF2_index+1)

          a  = (y1*x3 - y2*x3 - y1*x2 - y3*x1 + y2*x1 + y3*x2 ) / &
               (x2*x2*x1 - x3*x3*x1 + x3*x3*x2 - x2*x2*x3 + x1*x1*x3 - x1*x1*x2)

          b = (y1 - y2 - a*(x1*x1 - x2*x2)) / (x1 - x2)

          c = y1 - a*x1*x1 - b*x1

          x = -b / (2 * a)

          hmF2_km(m,l) = x
          NmF2(m,l) = (a * x * x) + (b * x) + c

          else

!         hmF2_km(m,l) = FLOAT(NmF2_index-1)*5. + 90.
          hmF2_km(m,l) = fixed_heights_km(NmF2_index)
          NmF2(m,l) = NE_1d(NmF2_index)

          endif

          TEC(m,l) = 0.0
          DO iheight = 1,interface_hts
              TEC(m,l) = TEC(m,l) + NE_1d(iheight) * 5000.
          enddo

          Te_400(m,l) = Te_fixed_ht_com(19,m,l)
          Ti_400(m,l) = Ti1_fixed_ht_com(19,m,l)

      1220 ENDDO
  1250 ENDDO

  return




end SUBROUTINE INTERFACE__calculate_NMF2_HMF2_TEC_etc




SUBROUTINE IO__read_apex_mag_field_coords( &
  ETA_Apex_3D, &
  Apex_D1,Apex_D2, &
  Apex_E1,Apex_E2, &
  Apex_BE3,Apex_bhat,Apex_Bmag, &
  Apex_D,Apex_d1d1,Apex_d1d2,Apex_d2d2,Apex_grdlbm2, &
  bcol,blon,q_coordinate,Re_apex, &
  gr,gcol,glon, &
  IN, IS, &
  glat_plasma_3d,glond_plasma_3d,pz_plasma_3d, &
  iday_number)

  IMPLICIT NONE


  INTEGER :: NPTS_dum , NMP_dum , NLP_dum
  INTEGER :: istop

  INTEGER :: IN(NMP,NLP) , IS(NMP,NLP)
  INTEGER :: IN_dum(NMP,NLP) , IS_dum(NMP,NLP)
  INTEGER :: i , j
  INTEGER :: mp , lp
  INTEGER :: iday_number
  INTEGER :: mgtype_dum
  REAL(kind=8) R0
  REAL(kind=8) PI
  REAL(kind=8) DTR
  REAL(kind=8) RE_apex(NMP,NLP)
  REAL(kind=8) Apex_D1(3,NPTS,NMP)
  REAL(kind=8) Apex_D2(3,NPTS,NMP)
  REAL(kind=8) Apex_E1(3,NPTS,NMP)
  REAL(kind=8) Apex_E2(3,NPTS,NMP)
  REAL(kind=8) Apex_BE3(NPTS,NMP)
  REAL(kind=8) Apex_grdlbm2(3,NPTS,NMP)
  REAL(kind=8) ETA_APEX_3D(NPTS,NMP)
  REAL(kind=8) ETA_APEX_1D(NPTS)
  REAL(kind=8) magnitude_e2_at_300
  REAL(kind=8) magnitude_e2_at_apex
  REAL(kind=8) Apex_BHAT(3,npts,nmp)
  REAL(kind=8) Apex_D(npts,nmp)
  REAL(kind=8) Apex_d1d1(npts,nmp)
  REAL(kind=8) Apex_d1d2(npts,nmp)
  REAL(kind=8) Apex_d2d2(npts,nmp)
  REAL(kind=8) Apex_BMAG(npts,nmp)
  REAL(kind=8) gr(NPTS,nmp)
  REAL(kind=8) gcol(NPTS,nmp)
  REAL(kind=8) glon(NPTS,nmp)
  REAL(kind=8) bcol(NPTS,NMP)
  REAL(kind=8) BLOn(NMP,NLP)
  REAL(kind=8) q_coordinate(NPTS,NMP)
  real(kind=8) Glat_plasma_3d(npts,nmp)
  real(kind=8) Glond_plasma_3d(npts,nmp)
  real(kind=8) Pz_plasma_3d(npts,nmp)
  integer ijump , iwrite_binary_to_114

  PARAMETER (PI=3.141592654,DTR=PI/180.0,R0=6.370E06)

  write(6,*) '******* Reading in APEX mag field parameters for GIP *******'

    write(6,*) '5'
      read(51,*) IN , IS
      read(51,*) gr, gcol, glon
      read(51,*) RE_apex , q_coordinate , BLOn
      read(51,*) bcol
    write(6,*) '4'
      read(51,*) Apex_D1
      read(51,*) Apex_D2
    write(6,*) '3'
      read(51,*) Apex_E1
      read(51,*) Apex_E2
      read(51,*) Apex_BE3
      read(51,*) Apex_grdlbm2
    write(6,*) '2'
      read(51,*) Eta_Apex_3d
      read(51,*) apex_D
      read(51,*) apex_BHAT
      read(51,*) apex_d1d1
    write(6,*) '1'
      read(51,*) apex_d1d2
      read(51,*) apex_d2d2
      read(51,*) apex_BMAG
    write(6,*) '....and all the interpolation indices'
      read(51,*) facfac_interface
      read(51,*) dd_interface
      read(51,*) ii1_interface
      read(51,*) ii2_interface
      read(51,*) ii3_interface
      read(51,*) ii4_interface
      read(51,*) mlow_interface
      read(51,*) mhigh_interface
    CLOSE (51)
    write(6,*) 'Blast Off'

   write(6,*) '******* Finished reading in APEX mag field parameters for GIP *******'

      iwrite_binary_to_114 = 0
      if(iwrite_binary_to_114 == 1) then

      OPEN (114,FILE='../static_files/apex_coordinates.written_out.binary', &
            FORM='UNFORMATTED',STATUS='UNKNOWN')

      write(114) IN , IS
      write(114) gr, gcol, glon
      write(114) RE_apex , q_coordinate , BLOn
      write(114) bcol
      write(114) Apex_D1
      write(114) Apex_D2
      write(114) Apex_E1
      write(114) Apex_E2
      write(114) Apex_BE3
      write(114) Apex_grdlbm2
      write(114) Eta_Apex_3d
      write(114) apex_D
      write(114) apex_BHAT
      write(114) apex_d1d1
      write(114) apex_d1d2
      write(114) apex_d2d2
      write(114) apex_BMAG

      close(114)

      write(6,*) '***** Written out Apex coordinates as a binary file (unit 114) : Stopping **'
      stop

      endif
!g
!g calculate 3D parameters for lat, lon, pz...
!g we can get rid of this later.....
!g
  do mp=1,nmp
      do lp=1,nlp
          do i=in(mp,lp),is(mp,lp)
              Pz_plasma_3d(i,mp) = (gr(i,mp)-R0)/1000.
              Glat_plasma_3d(i,mp) = 90. - gcol(i,mp)/DTR
              Glond_plasma_3d(i,mp) = glon(i,mp)/DTR
          enddo
      enddo
  enddo

  RETURN




end SUBROUTINE IO__read_apex_mag_field_coords




SUBROUTINE IO__write_gip_netcdf_history(filename, &
             ni,vi,ne,ti,te,no_plus_3d,o2_plus_3d,d13d,d23d,v13d,v23d, &
             vpeq,vzon, &
             Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
             NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
             N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
             Te_fixed_ht_com,Ti1_fixed_ht_com,Ti2_fixed_ht_com, &
             Ne_density_fixed_ht_com)


!     include 'netcdf.inc'

      character(len=*),intent(in) :: filename
      integer :: ncid,istat,id
      integer :: id_npts, id_nmp, id_nlp, id_which_ion, id_highlat_hts, &
        id_highlat_lats, id_highlat_lons, id_fixht_hts, id_fixht_lats, &
        id_fixht_lons
      integer ::  &
        idv_ti     ,idv_te   ,idv_ni   ,idv_vi   ,idv_ne   ,idv_noplus, &
        idv_o2plus ,idv_d13d ,idv_d23d ,idv_v13d ,idv_v23d ,idv_vpeq, &
        idv_vzon   , &
        idv_fixht_oplus  ,idv_fixht_hplus  ,idv_fixht_noplus, & 
        idv_fixht_o2plus ,idv_fixht_n2plus ,idv_fixht_nplus , &
        idv_fixht_te     ,idv_fixht_ti1    ,idv_fixht_ti2   , &
        idv_fixht_ne
      integer :: ids1(1),ids2(2),ids3(3),ids4(4) ! vectors of dim id's
      character(len=120) :: long_name

      real(kind=8),dimension(npts,nmp) :: te,ne,no_plus_3d,o2_plus_3d
      real(kind=8),dimension(npts,nmp,which_ion) :: ti,ni,vi
      real(kind=8),dimension &
        (high_lat_hts, high_lat_lats, high_lat_lons, which_ion) :: &
        d13d, d23d, v13d, v23d
      real(kind=8),dimension(nmp,nlp) :: vpeq,vzon
      real(kind=8),dimension &
        (interface_hts, fixed_height_lats, fixed_height_lons) :: &
        oplus_density_fixed_ht_com , hplus_density_fixed_ht_com, & 
        noplus_density_fixed_ht_com, o2plus_density_fixed_ht_com, &
        n2plus_density_fixed_ht_com, nplus_density_fixed_ht_com, &
        te_fixed_ht_com, ti1_fixed_ht_com, ti2_fixed_ht_com, &
        ne_density_fixed_ht_com

!
! Create new data set (overwrite any pre-existing file):
!
      istat = nf_create(filename,NF_CLOBBER,ncid) 
      if (istat /= NF_NOERR)  &
        call IO__handle_ncerr(istat,'Error from nf_create',1)
      write(6,"(/,'write_hist: Created netcdf file ',a)") trim(filename) 
!
! Define dimensions:
!
      istat = nf_def_dim(ncid,"npts",npts,id_npts)
      istat = nf_def_dim(ncid,"nmp" ,nmp,id_nmp)
      istat = nf_def_dim(ncid,"nlp" ,nlp,id_nlp)
      istat = nf_def_dim(ncid,"which_ion" ,which_ion,id_which_ion)
      istat = nf_def_dim(ncid,"high_lat_hts" ,high_lat_hts , &
        id_highlat_hts)
      istat = nf_def_dim(ncid,"high_lat_lats",high_lat_lats, &
        id_highlat_lats)
      istat = nf_def_dim(ncid,"high_lat_lons",high_lat_lons, &
        id_highlat_lons)
      istat = nf_def_dim(ncid,"interface_hts" ,interface_hts, &
        id_fixht_hts)
      istat = nf_def_dim(ncid,"fixed_height_lats",fixed_height_lats, &
        id_fixht_lats)
      istat = nf_def_dim(ncid,"fixed_height_lons",fixed_height_lons, &
        id_fixht_lons)
      write(6,"('Defined dimensions on file ',a)") trim(filename)
!
! Define variables:
!
! 2d and 3d variables:
      ids2 = (/id_npts,id_nmp/)
      ids3 = (/id_npts,id_nmp,id_which_ion/)

      istat = nf_def_var(ncid,"TI",NF_DOUBLE,3,ids3,idv_ti)
      long_name = "Ion Temperature"
      istat = nf_put_att_text(ncid,idv_ti,"long_name", &
        len_trim(long_name),trim(long_name))
      istat = nf_put_att_text(ncid,idv_ti,"units",5,"deg K")

      istat = nf_def_var(ncid,"TE",NF_DOUBLE,2,ids2,idv_te)
      long_name = "Electron Temperature"
      istat = nf_put_att_text(ncid,idv_te,"long_name", &
        len_trim(long_name),trim(long_name))
      istat = nf_put_att_text(ncid,idv_te,"units",5,"deg K")

      istat = nf_def_var(ncid,"NI",NF_DOUBLE,3,ids3,idv_ni)
      istat = nf_def_var(ncid,"VI",NF_DOUBLE,3,ids3,idv_vi)
      istat = nf_def_var(ncid,"NE",NF_DOUBLE,2,ids2,idv_ne)
      istat = nf_def_var(ncid,"no_plus_3d",NF_DOUBLE,2,ids2,idv_noplus)
      istat = nf_def_var(ncid,"o2_plus_3d",NF_DOUBLE,2,ids2,idv_o2plus)
!
! High-lat variables (4d):
!
      ids4 = (/id_highlat_hts, id_highlat_lats, id_highlat_lons, &
               id_which_ion/)
      istat = nf_def_var(ncid,"d13d",NF_DOUBLE,4,ids4,idv_d13d)
      istat = nf_def_var(ncid,"d23d",NF_DOUBLE,4,ids4,idv_d23d)
      istat = nf_def_var(ncid,"v13d",NF_DOUBLE,4,ids4,idv_v13d)
      istat = nf_def_var(ncid,"v23d",NF_DOUBLE,4,ids4,idv_v23d)
!
! 2d (nmp,nlp):
      ids2 = (/id_nmp,id_nlp/)
      istat = nf_def_var(ncid,"vpeq",NF_DOUBLE,2,ids2,idv_vpeq)
      istat = nf_def_var(ncid,"vzon",NF_DOUBLE,2,ids2,idv_vzon)
!
! Fixed height variables (3d):
!
      ids3 = (/id_fixht_hts, id_fixht_lats, id_fixht_lons/)
      istat = nf_def_var(ncid,"oplus_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_oplus)
      istat = nf_def_var(ncid,"hplus_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_hplus)
      istat = nf_def_var(ncid,"noplus_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_noplus)
      istat = nf_def_var(ncid,"o2plus_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_o2plus)
      istat = nf_def_var(ncid,"n2plus_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_n2plus)
      istat = nf_def_var(ncid,"nplus_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_nplus)
      istat = nf_def_var(ncid,"te_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_te)
      istat = nf_def_var(ncid,"ti1_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_ti1)
      istat = nf_def_var(ncid,"ti2_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_ti2)
      istat = nf_def_var(ncid,"ne_density_fixed_ht_com", &
        NF_DOUBLE,3,ids3,idv_fixht_ne)

!     write(6,"('Defined variables on file ',a)") trim(filename)
!
! Take out of define mode:
      istat = nf_enddef(ncid)
!
! Write variables to the file:
!     write(6,"('Writing variables to file ',a,'..')") trim(filename)

      istat = nf_put_var_double(ncid,idv_ti,ti)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var TI',0)

      istat = nf_put_var_double(ncid,idv_te,te)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var TE',0)

      istat = nf_put_var_double(ncid,idv_ni,ni)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var NI',0)

      istat = nf_put_var_double(ncid,idv_vi,vi)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var VI',0)

      istat = nf_put_var_double(ncid,idv_ne,ne)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var NE',0)

      istat = nf_put_var_double(ncid,idv_noplus,no_plus_3d)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var no_plus_3d',0)

      istat = nf_put_var_double(ncid,idv_o2plus,o2_plus_3d)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var o2_plus_3d',0)

      istat = nf_put_var_double(ncid,idv_d13d,d13d)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var d13d',0)

      istat = nf_put_var_double(ncid,idv_d23d,d23d)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var d23d',0)

      istat = nf_put_var_double(ncid,idv_v13d,v13d)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var v13d',0)

      istat = nf_put_var_double(ncid,idv_v23d,v23d)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var v23d',0)

      istat = nf_put_var_double(ncid,idv_vpeq,vpeq)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var vpeq',0)

      istat = nf_put_var_double(ncid,idv_vzon,vzon)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var vzon',0)

      istat = nf_put_var_double(ncid,idv_fixht_oplus, &
        oplus_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var oplus_density_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_hplus, &
        hplus_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var hplus_density_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_noplus, &
        noplus_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var noplus_density_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_o2plus, &
        o2plus_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var o2plus_density_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_n2plus, &
        n2plus_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var n2plus_density_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_nplus, &
        nplus_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var nplus_density_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_te, &
        te_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var te_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_ti1, &
        ti1_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var ti1_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_ti2, &
        ti2_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var ti2_fixed_ht_com',0)

      istat = nf_put_var_double(ncid,idv_fixht_ne, &
        ne_density_fixed_ht_com)
      if (istat /= NF_NOERR) call IO__handle_ncerr(istat, &
        'Error writing var ne_density_fixed_ht_com',0)
!
! Close data set:
!
      istat = nf_close(ncid)
      if (istat /= NF_NOERR) &
        call IO__handle_ncerr(istat,'Error from nf_close',0)

end SUBROUTINE IO__write_gip_netcdf_history

!-----------------------------------------------------------------------



SUBROUTINE IO__read_gip_netcdf_history(filename, &
             ni,vi,ne,ti,te,no_plus_3d,o2_plus_3d,d13d,d23d,v13d,v23d, &
             vpeq,vzon, &
             Oplus_density_fixed_ht_com,Hplus_density_fixed_ht_com, &
             NOplus_density_fixed_ht_com,O2plus_density_fixed_ht_com, &
             N2plus_density_fixed_ht_com,Nplus_density_fixed_ht_com, &
             Te_fixed_ht_com,Ti1_fixed_ht_com,Ti2_fixed_ht_com, &
             Ne_density_fixed_ht_com)


!     include 'netcdf.inc'

      character(len=*),intent(in) :: filename
      integer :: ncid,istat,id
      integer :: id_npts, id_nmp, id_nlp, id_which_ion, id_highlat_hts, &
        id_highlat_lats, id_highlat_lons, id_fixht_hts, id_fixht_lats, &
        id_fixht_lons
      integer ::  &
        idv_ti     ,idv_te   ,idv_ni   ,idv_vi   ,idv_ne   ,idv_noplus, &
        idv_o2plus ,idv_d13d ,idv_d23d ,idv_v13d ,idv_v23d ,idv_vpeq, &
        idv_vzon   , &
        idv_fixht_oplus  ,idv_fixht_hplus  ,idv_fixht_noplus, & 
        idv_fixht_o2plus ,idv_fixht_n2plus ,idv_fixht_nplus , &
        idv_fixht_te     ,idv_fixht_ti1    ,idv_fixht_ti2   , &
        idv_fixht_ne
      integer :: ids1(1),ids2(2),ids3(3),ids4(4) ! vectors of dim id's
      character(len=120) :: long_name

      real(kind=8),dimension(npts,nmp) :: te,ne,no_plus_3d,o2_plus_3d
      real(kind=8),dimension(npts,nmp,which_ion) :: ti,ni,vi
      real(kind=8),dimension &
        (high_lat_hts, high_lat_lats, high_lat_lons, which_ion) :: &
        d13d, d23d, v13d, v23d
      real(kind=8),dimension(nmp,nlp) :: vpeq,vzon
      real(kind=8),dimension &
        (interface_hts, fixed_height_lats, fixed_height_lons) :: &
        oplus_density_fixed_ht_com , hplus_density_fixed_ht_com, & 
        noplus_density_fixed_ht_com, o2plus_density_fixed_ht_com, &
        n2plus_density_fixed_ht_com, nplus_density_fixed_ht_com, &
        te_fixed_ht_com, ti1_fixed_ht_com, ti2_fixed_ht_com, &
        ne_density_fixed_ht_com



!
! Open file for reading:
      istat = nf_open(filename,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) &
        call IO__handle_ncerr(istat,'Error from nf_open',1)
      write(6,"(/,'read_hist: Opened file ',a,' for reading.')") &
        trim(filename)

      istat = nf_inq_varid(ncid,"TE",id)
      istat = nf_get_var_double(ncid,id,te)
!     write(6,"('read_hist: TE min,max=',2e12.4)") minval(te),maxval(te)

      istat = nf_inq_varid(ncid,"TI",id)
      istat = nf_get_var_double(ncid,id,ti)
!     write(6,"('read_hist: TI min,max=',2e12.4)") minval(ti),maxval(ti)

      istat = nf_inq_varid(ncid,"NI",id)
      istat = nf_get_var_double(ncid,id,ni)
!     write(6,"('read_hist: NI min,max=',2e12.4)") minval(ni),maxval(ni)

      istat = nf_inq_varid(ncid,"VI",id)
      istat = nf_get_var_double(ncid,id,vi)
!     write(6,"('read_hist: VI min,max=',2e12.4)") minval(vi),maxval(vi)

      istat = nf_inq_varid(ncid,"NE",id)
      istat = nf_get_var_double(ncid,id,ne)
!     write(6,"('read_hist: NE min,max=',2e12.4)") minval(ne),maxval(ne)

      istat = nf_inq_varid(ncid,"no_plus_3d",id)
      istat = nf_get_var_double(ncid,id,no_plus_3d)
!     write(6,"('read_hist: no_plus_3d min,max=',2e12.4)") &
!       minval(no_plus_3d),maxval(no_plus_3d)

      istat = nf_inq_varid(ncid,"o2_plus_3d",id)
      istat = nf_get_var_double(ncid,id,o2_plus_3d)
!     write(6,"('read_hist: o2_plus_3d min,max=',2e12.4)") &
!       minval(o2_plus_3d),maxval(o2_plus_3d)

      istat = nf_inq_varid(ncid,"d13d",id)
      istat = nf_get_var_double(ncid,id,d13d)
!     write(6,"('read_hist: d13d min,max=',2e12.4)") &
!       minval(d13d),maxval(d13d)

      istat = nf_inq_varid(ncid,"d23d",id)
      istat = nf_get_var_double(ncid,id,d23d)
!     write(6,"('read_hist: d23d min,max=',2e12.4)") &
!       minval(d23d),maxval(d23d)

      istat = nf_inq_varid(ncid,"v13d",id)
      istat = nf_get_var_double(ncid,id,v13d)
!     write(6,"('read_hist: v13d min,max=',2e12.4)") &
!       minval(v13d),maxval(v13d)

      istat = nf_inq_varid(ncid,"v23d",id)
      istat = nf_get_var_double(ncid,id,v23d)
!     write(6,"('read_hist: v23d min,max=',2e12.4)") &
!       minval(v23d),maxval(v23d)

      istat = nf_inq_varid(ncid,"vpeq",id)
      istat = nf_get_var_double(ncid,id,vpeq)
!     write(6,"('read_hist: vpeq min,max=',2e12.4)") &
!       minval(vpeq),maxval(vpeq)

      istat = nf_inq_varid(ncid,"vzon",id)
      istat = nf_get_var_double(ncid,id,vzon)
!     write(6,"('read_hist: vzon min,max=',2e12.4)") &
!       minval(vzon),maxval(vzon)

      istat = nf_inq_varid(ncid,"oplus_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,oplus_density_fixed_ht_com)
!     write(6,"('read_hist: oplus_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(oplus_density_fixed_ht_com), &
!       maxval(oplus_density_fixed_ht_com)

      istat = nf_inq_varid(ncid,"hplus_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,hplus_density_fixed_ht_com)
!     write(6,"('read_hist: hplus_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(hplus_density_fixed_ht_com), &
!       maxval(hplus_density_fixed_ht_com)

      istat = nf_inq_varid(ncid,"noplus_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,noplus_density_fixed_ht_com)
!     write(6,"('read_hist: noplus_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(noplus_density_fixed_ht_com), &
!       maxval(noplus_density_fixed_ht_com)

      istat = nf_inq_varid(ncid,"o2plus_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,o2plus_density_fixed_ht_com)
!     write(6,"('read_hist: o2plus_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(o2plus_density_fixed_ht_com), &
!       maxval(o2plus_density_fixed_ht_com)

      istat = nf_inq_varid(ncid,"n2plus_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,n2plus_density_fixed_ht_com)
!     write(6,"('read_hist: n2plus_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(n2plus_density_fixed_ht_com), &
!       maxval(n2plus_density_fixed_ht_com)

      istat = nf_inq_varid(ncid,"nplus_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,nplus_density_fixed_ht_com)
!     write(6,"('read_hist: nplus_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(nplus_density_fixed_ht_com), &
!       maxval(nplus_density_fixed_ht_com)

      istat = nf_inq_varid(ncid,"te_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,te_fixed_ht_com)
!     write(6,"('read_hist: te_fixed_ht_com min,max=',2e12.4)") &
!       minval(te_fixed_ht_com),maxval(te_fixed_ht_com)

      istat = nf_inq_varid(ncid,"ti1_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,ti1_fixed_ht_com)
!     write(6,"('read_hist: ti1_fixed_ht_com min,max=',2e12.4)") &
!       minval(ti1_fixed_ht_com),maxval(ti1_fixed_ht_com)

      istat = nf_inq_varid(ncid,"ti2_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,ti2_fixed_ht_com)
!     write(6,"('read_hist: ti2_fixed_ht_com min,max=',2e12.4)") &
!       minval(ti2_fixed_ht_com),maxval(ti2_fixed_ht_com)

      istat = nf_inq_varid(ncid,"ne_density_fixed_ht_com",id)
      istat = nf_get_var_double(ncid,id,ne_density_fixed_ht_com)
!     write(6,"('read_hist: ne_density_fixed_ht_com min,max=',2e12.4)") &
!       minval(ne_density_fixed_ht_com), &
!       maxval(ne_density_fixed_ht_com)

      istat = nf_close(ncid)
      write(6,"('Completed read of file ',a)") trim(filename)



end SUBROUTINE IO__read_gip_netcdf_history


!-----------------------------------------------------------------------

SUBROUTINE IO__handle_ncerr(istat,msg,ifatal)

!     include 'netcdf.inc'
!
! Handle a netcdf lib error:
!
      integer,intent(in) :: istat,ifatal
      character(len=*),intent(in) :: msg
!
      write(6,"(/72('-'))")
      write(6,"('>>> Error from netcdf library:')")
      write(6,"(a)") trim(msg)
      write(6,"('istat=',i5)") istat
      write(6,"(a)") nf_strerror(istat)
      write(6,"(72('-')/)")
      if (ifatal > 0) stop 'IO__handle_ncerr'
      return
end SUBROUTINE IO__handle_ncerr









SUBROUTINE get_HYD400km_for_this_run(static_file_location,&
                                     iday_of_the_year_input,F107_input)
!                                     Hyd_grid_400km_m3,Hyd_grid_lats, &
!                                     Hyd_grid_local_times)

implicit none

character(100) :: static_file_location
integer f_low , f_high
integer idaylow , idayhigh
integer if10 , iday
integer :: iday_of_the_year_input
real(kind=8) :: F107_input
real(kind=8) :: myday
real(kind=8) :: myf107
real(kind=8) :: the_days(10)
real(kind=8) :: f10s(4)
real(kind=8) :: f_factor
real(kind=8) :: day_factor

data the_days/1.,40.,80.,126.,172.,219.,266.,310.,355.,366./
data f10s/80.,120.,160.,200./

real(kind=8) :: Hden(4,9,24,19)
real(kind=8) :: Hden2(4,10,24,19)
real(kind=8) :: Hden_daylow_flow(24,19)
real(kind=8) :: Hden_daylow_fhigh(24,19)
real(kind=8) :: Hden_dayhigh_flow(24,19)
real(kind=8) :: Hden_dayhigh_fhigh(24,19)
real(kind=8) :: Hden_daylow(24,19)
real(kind=8) :: Hden_dayhigh(24,19)
real(kind=8) :: Hden_final(24,19)
!real(kind=8) :: Hyd_grid_400km_m3(25,19)

!real(kind=8) :: Hyd_grid_lats(19)
!real(kind=8) :: Hyd_grid_local_times(25)




OPEN(140,FILE=TRIM(static_file_location)//'Hydrogen_file',STATUS='old')
READ(140,*) Hden
CLOSE(140)

Hden2(:,1:9,:,:) = Hden(:,1:9,:,:)
Hden2(:,10,:,:) = Hden2(:,1,:,:)

myday = dble(iday_of_the_year_input)
myf107 = F107_input

if (myf107.lt.80.) myf107 = 80.
if (myf107.gt.200.) myf107 = 200.


do if10 = 1 , 3 
if (myf107.ge.f10s(if10).and.myf107.lt.f10s(if10+1)) then
f_low = if10
f_high = if10+1
goto 1234
endif
enddo

f_low = 3
f_high = 4

1234 continue


f_factor = (myf107 - f10s(f_low)) / 40.

if (myday.ge.1..and.myday.le.366.) then

do iday = 1 , 10
if (myday.gt.the_days(iday).and.myday.le.the_days(iday+1)) then
idaylow = iday
idayhigh = iday+1
goto 2345
endif
enddo

2345 continue


day_factor = (myday - the_days(idaylow)) / (the_days(idayhigh) - the_days(idaylow))



Hden_daylow_flow(:,:) = Hden2(f_low,idaylow,:,:)
Hden_daylow_fhigh(:,:) = Hden2(f_high,idaylow,:,:)
Hden_dayhigh_flow(:,:) = Hden2(f_low,idayhigh,:,:)
Hden_dayhigh_fhigh(:,:) = Hden2(f_high,idayhigh,:,:)

Hden_daylow(:,:) = ((Hden_daylow_fhigh(:,:)-Hden_daylow_flow(:,:)) * f_factor) + Hden_daylow_flow(:,:)
Hden_dayhigh(:,:) = ((Hden_dayhigh_fhigh(:,:)-Hden_dayhigh_flow(:,:)) * f_factor) + Hden_dayhigh_flow(:,:)

Hden_final(:,:) = ((Hden_dayhigh(:,:) - Hden_daylow(:,:)) * day_factor) +  Hden_daylow(:,:)

Hyd_grid_400km_m3(1:24,:) = Hden_final(1:24,:)
Hyd_grid_400km_m3(25,:) = Hyd_grid_400km_m3(1,:)

Hyd_grid_400km_m3(:,:) = Hyd_grid_400km_m3(:,:) * 1.e11


else

write(6,*) 'day needs to be between 1 and 366 innit ?'

endif


return

END SUBROUTINE get_HYD400km_for_this_run










SUBROUTINE HL__POLAR_IONOSPHERE( &
      GIP_switches, &
      geo_grid_longitudes_degrees,geo_grid_longitudes_radians, &
      geo_grid_latitudes_degrees,geo_grid_latitudes_radians, &
      geo_grid_colatitudes_degrees,geo_grid_colatitudes_radians, &
  high_lat_time_step_seconds, &
  angdif, &
  btotal,dip_angle_radians,f107,hprof,dth_radians, &
  ex2d,ey2d,v13d,v23d,tarea,d13d,d23d, &
  gravity,nhgt,hz_metres,hnew,h2,h2diff,hprod2,hprod3,zkm,l300,l1000,r0,declination, &
  Solar_Declination_Angle_degrees,km_plasma,universal_time_seconds, &
  ioutput_counter,ioutput_high_res_counter, &
  ioutput_frequency_mins,ioutput_high_res_frequency_mins, &
  ihigh_lat_call_frequency_mins,file_res, &
  sw_initialisation_call, &
  O_density_fixed_ht,O2_density_fixed_ht,N2_density_fixed_ht, &
  NO_density_fixed_ht, &
  N4S_density_fixed_ht,N2D_density_fixed_ht, &
  Vx_fixed_ht,Vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
  qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
  qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
  factor_ht_3d,iht_above_3d,iht_below_3d, &
  factor_ht_inverse_3d,iht_above_inverse_3d,iht_below_inverse_3d, &
  altitude_metres_3d, &
  Oplus_density_fixed_ht,Hplus_density_fixed_ht,NOplus_density_fixed_ht, &
  O2plus_density_fixed_ht,N2plus_density_fixed_ht,Nplus_density_fixed_ht, &
  Te_fixed_ht,Ti1_fixed_ht,Ti2_fixed_ht,iday_number)
!g
  IMPLICIT NONE

    INTEGER N_Pressure_Levels,N_Latitudes,N_longitudes
      PARAMETER(N_Pressure_Levels=15)
      PARAMETER(N_Latitudes=91)
      PARAMETER(N_Longitudes=20)

  LOGICAL :: sw_initialisation_call
  LOGICAL :: GIP_switches(20)
  INTEGER :: l
  INTEGER :: m
  INTEGER :: n
  INTEGER :: iday_number
  INTEGER :: ihemi
  INTEGER :: nhgt
  INTEGER :: l300
  INTEGER :: l1000
  INTEGER :: istop
  INTEGER :: ioutput_counter
  INTEGER :: ioutput_high_res_counter
  INTEGER :: ioutput_frequency_mins
  INTEGER :: ioutput_high_res_frequency_mins
  INTEGER :: file_res(27)
  INTEGER :: mhigh(20)
  INTEGER :: mlow(20)
  INTEGER :: i , j
  INTEGER :: ik
  INTEGER :: ihigh_lat_call_frequency_mins
  INTEGER :: ilt
  INTEGER :: ilt_east
  INTEGER :: ilt_west
  INTEGER :: ilat
  INTEGER :: ilat_north
  INTEGER :: ilat_south
  REAL(kind=8) factor_lt
  REAL(kind=8) factor_lat
  REAL(kind=8) H_NE
  REAL(kind=8) H_NW
  REAL(kind=8) H_SE
  REAL(kind=8) H_SW
  REAL(kind=8) H_N
  REAL(kind=8) H_S
  REAL(kind=8) H_400
  REAL(kind=8) dtr
  REAL(kind=8) pi
  REAL(kind=8) qion_1d(90)
  REAL(kind=8) qo2p_aurora_1d(90)
  REAL(kind=8) qop_aurora_1d(90)
  REAL(kind=8) qn2p_aurora_1d(90)
  REAL(kind=8) qnp_aurora_1d(90)
  REAL(kind=8) qtef_aurora_1d(90)
  REAL(kind=8) high_lat_time_step_seconds
  REAL(kind=8) ssa
  REAL(kind=8) geo_local_time_degrees(N_longitudes)
  REAL(kind=8) geo_local_time_hours(N_longitudes)
  REAL(kind=8) angdif(N_Latitudes,N_longitudes)
  REAL(kind=8) btotal(N_Latitudes,N_longitudes)
  REAL(kind=8) dip_angle_radians(N_Latitudes,N_longitudes)
  REAL(kind=8) solar_zenith_angle_radians(N_Latitudes,N_longitudes)
  REAL(kind=8) coschi
  REAL(kind=8) f107
  REAL(kind=8) hprof(N_longitudes,19)
  REAL(kind=8) dth_radians
  REAL(kind=8) ti1_1D_output(N_Pressure_Levels)
  REAL(kind=8) ex2d(N_Latitudes,N_longitudes)
  REAL(kind=8) ey2d(N_Latitudes,N_longitudes)
  REAL(kind=8) d13d(90,N_Latitudes,N_longitudes,2)
  REAL(kind=8) d23d(90,N_Latitudes,N_longitudes,2)
  REAL(kind=8) v13d(90,N_Latitudes,N_longitudes,2)
  REAL(kind=8) v23d(90,N_Latitudes,N_longitudes,2)
  REAL(kind=8) tarea(90)
  REAL(kind=8) gravity(90)
  REAL(kind=8) hz_metres(90)
  REAL(kind=8) h2(90)
  REAL(kind=8) h2diff(90)
  REAL(kind=8) hnew(90)
  REAL(kind=8) hprod2(90)
  REAL(kind=8) hprod3(90)
  REAL(kind=8) zkm(90)
  REAL(kind=8) r0
  REAL(kind=8) r0_km
  REAL(kind=8) g0
  REAL(kind=8) declination(N_Latitudes,N_longitudes)
  REAL(kind=8) solar_declination_angle_degrees
  REAL(kind=8) km_plasma(6)
  REAL(kind=8) universal_time_seconds
  REAL(kind=8) topht
  REAL(kind=8) height_in
  REAL(kind=8) georat
  REAL(kind=8) geo_grid_longitudes_degrees(N_longitudes)
  REAL(kind=8) geo_grid_longitudes_radians(N_longitudes)
  REAL(kind=8) geo_grid_latitudes_degrees(N_Latitudes)
  REAL(kind=8) geo_grid_latitudes_radians(N_Latitudes)
  REAL(kind=8) geo_grid_colatitudes_degrees(N_Latitudes)
  REAL(kind=8) geo_grid_colatitudes_radians(N_Latitudes)
  REAL(kind=8) :: O_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: O2_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N2_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: NO_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N4S_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N2D_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Vx_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Vy_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: wvz_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: tts_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qion3d_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qo2p_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qop_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qn2p_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qnp_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: qtef_aurora_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: o_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: o2_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: n2_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: NO_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: N4S_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: N2D_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: Vx_high_lat_ions(90,91,20)
  REAL(kind=8) :: Vy_high_lat_ions(90,91,20)
  REAL(kind=8) :: wvz_high_lat_ions(90,91,20)
  REAL(kind=8) :: tn_high_lat_ions(90,91,20)
  REAL(kind=8) :: qion3d_high_lat_ions(90,91,20)
  REAL(kind=8) :: qo2p_aurora_high_lat_ions(90,91,20)
  REAL(kind=8) :: qop_aurora_high_lat_ions(90,91,20)
  REAL(kind=8) :: qn2p_aurora_high_lat_ions(90,91,20)
  REAL(kind=8) :: qnp_aurora_high_lat_ions(90,91,20)
  REAL(kind=8) :: qtef_aurora_high_lat_ions(90,91,20)
  REAL(kind=8) :: O_density_1D(90)
  REAL(kind=8) :: O2_density_1D(90)
  REAL(kind=8) :: N2_density_1D(90)
  REAL(kind=8) :: NO_density_1D(90)
  REAL(kind=8) :: N4S_density_1D(90)
  REAL(kind=8) :: N2D_density_1D(90)
  REAL(kind=8) :: Vx_1d(90)
  REAL(kind=8) :: Vy_1D(90)
  REAL(kind=8) :: wvz_1D(90)
  REAL(kind=8) :: tn_1D(90)
  REAL(kind=8) :: Oplus_density_1d(90)
  REAL(kind=8) :: Hplus_density_1d(90)
  REAL(kind=8) :: NOplus_density_1d(90)
  REAL(kind=8) :: O2plus_density_1d(90)
  REAL(kind=8) :: N2plus_density_1d(90)
  REAL(kind=8) :: Nplus_density_1d(90)
  REAL(kind=8) :: Te_1d(90)
  REAL(kind=8) :: Ti1_1d(90)
  REAL(kind=8) :: Ti2_1d(90)
  REAL(kind=8) :: Oplus_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: Hplus_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: NOplus_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: O2plus_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: N2plus_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: Nplus_density_high_lat_ions(90,91,20)
  REAL(kind=8) :: Te_high_lat_ions(90,91,20)
  REAL(kind=8) :: Ti1_high_lat_ions(90,91,20)
  REAL(kind=8) :: Ti2_high_lat_ions(90,91,20)
  REAL(kind=8) :: Oplus_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Hplus_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: NOplus_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: O2plus_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: N2plus_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Nplus_density_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Te_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Ti1_fixed_ht(interface_hts,91,20)
  REAL(kind=8) :: Ti2_fixed_ht(interface_hts,91,20)
  REAL(kind=8) ::  rlt , rlat
  REAL(kind=8) :: oplus_polar(90)
  REAL(kind=8) :: hplus_polar(90)
  REAL(kind=8) :: noplus_polar(90)
  REAL(kind=8) :: o2plus_polar(90)
  REAL(kind=8) :: n2plus_polar(90)
  REAL(kind=8) :: nplus_polar(90)
  REAL(kind=8) :: Te_polar(90)
  REAL(kind=8) :: Ti1_polar(90)
  REAL(kind=8) :: Ti2_polar(90)

  REAL(kind=8) factor_ht_3d(90,91,20)
  INTEGER iht_above_3d(90,91,20)
  INTEGER iht_below_3d(90,91,20)
  REAL(kind=8) factor_ht_1d(90)
  INTEGER iht_above_1d(90)
  INTEGER iht_below_1d(90)

  REAL(kind=8) factor_ht_inverse_3d(interface_hts,91,20)
  INTEGER iht_above_inverse_3d(interface_hts,91,20)
  INTEGER iht_below_inverse_3d(interface_hts,91,20)
  REAL(kind=8) factor_ht_inverse_1d(interface_hts)
  INTEGER iht_above_inverse_1d(interface_hts)
  INTEGER iht_below_inverse_1d(interface_hts)

  REAL(kind=8) altitude_metres_3d(90,91,20)
  REAL(kind=8) altitude_metres(90)

  PARAMETER (PI=3.14159,DTR=PI/180.0)
!g
  data mlow /14,16,18,20,22,24,25,25,23,21,20,17,15,13,12,12,12,12,12,13/
  data mhigh /75,75,75,76,76,76,76,76,76,76,76,75,74,73,71,69,70,71,73,74/

  if(sw_initialisation_call) then

  ! INITIALISATION PART........

      DATA r0_km , g0/6371. , 9.81/


      topht = 10000.
      CALL HL__calculate_points_along_the_tube(topht,NHGt,ZKM,TARea,HNEw,H2,H2Diff,HPRod2,HPRod3)

      DO 1001 i = 1 , NHGt
          georat = 1. + ZKM(i)/r0_km
          gravity(i) = g0/(georat*georat)
          hz_metres(i) = ZKM(i)*1.E+03
      1001 ENDDO

      do l = 1 , 20

      do m = 1 , 40

        call HL__calculate_altitudes (NHGT,hz_metres,dip_angle_radians(m,l),altitude_metres)
        altitude_metres_3d(1:nhgt,m,l) = altitude_metres(1:nhgt)

        call interface__init_for_high_lat_ions (nhgt,altitude_metres,iht_above_1d,iht_below_1d,factor_ht_1d, &
                                                     iht_above_inverse_1d,iht_below_inverse_1d,factor_ht_inverse_1d)

        iht_above_3d(1:nhgt,m,l) = iht_above_1d(1:nhgt)
        iht_below_3d(1:nhgt,m,l) = iht_below_1d(1:nhgt)
        factor_ht_3d(1:nhgt,m,l) = factor_ht_1d(1:nhgt)
        iht_above_inverse_3d(1:nhgt,m,l) = iht_above_inverse_1d(1:nhgt)
        iht_below_inverse_3d(1:nhgt,m,l) = iht_below_inverse_1d(1:nhgt)
        factor_ht_inverse_3d(1:nhgt,m,l) = factor_ht_inverse_1d(1:nhgt)

      enddo

      do m= 50,91

        call HL__calculate_altitudes (NHGT,hz_metres,dip_angle_radians(m,l),altitude_metres)
        altitude_metres_3d(1:nhgt,m,l) = altitude_metres(1:nhgt)

        call interface__init_for_high_lat_ions (nhgt,altitude_metres,iht_above_1d,iht_below_1d,factor_ht_1d, &
                                                     iht_above_inverse_1d,iht_below_inverse_1d,factor_ht_inverse_1d)

        iht_above_3d(1:nhgt,m,l) = iht_above_1d(1:nhgt)
        iht_below_3d(1:nhgt,m,l) = iht_below_1d(1:nhgt)
        factor_ht_3d(1:nhgt,m,l) = factor_ht_1d(1:nhgt)
        iht_above_inverse_3d(1:nhgt,m,l) = iht_above_inverse_1d(1:nhgt)
        iht_below_inverse_3d(1:nhgt,m,l) = iht_below_inverse_1d(1:nhgt)
        factor_ht_inverse_3d(1:nhgt,m,l) = factor_ht_inverse_1d(1:nhgt)

      enddo

      enddo


  ! we need the following height levels:
  ! 300 & 1000 km for HL__ETEMP

      ik = 1
      height_in=400.
      CALL HL__NEARHT(ZKM,height_in,ik,NHGt,L300)
      ik = L300 + 1
      height_in=1000.
      CALL HL__NEARHT(ZKM,height_in,ik,NHGt,L1000)


      READ (22,99002) ((hprof(i,j),i=1,20),j=1,19)
      99002 FORMAT (1x,8F8.0)
      CLOSE (22)

      dth_radians = 2.0 * dtr


  else  ! else we have a normal call (not initialisation) ....


  call INTERFACE__FIXED_GEO_to_high_lat_ions ( &
           o_density_fixed_ht,o2_density_fixed_ht,n2_density_fixed_ht, &
           NO_density_fixed_ht, &
           N4S_density_fixed_ht,N2D_density_fixed_ht, &
           Vx_fixed_ht,Vy_fixed_ht,wvz_fixed_ht,tts_fixed_ht,qion3d_fixed_ht, &
           qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &
           qnp_aurora_fixed_ht, qtef_aurora_fixed_ht, &
           o_density_high_lat_ions,o2_density_high_lat_ions,n2_density_high_lat_ions, &
           NO_density_high_lat_ions, &
           N4S_density_high_lat_ions,N2D_density_high_lat_ions, &
           Vx_high_lat_ions,Vy_high_lat_ions,wvz_high_lat_ions,tn_high_lat_ions,qion3d_high_lat_ions, &
           qo2p_aurora_high_lat_ions, qop_aurora_high_lat_ions, qn2p_aurora_high_lat_ions, &
           qnp_aurora_high_lat_ions, qtef_aurora_high_lat_ions, &
           factor_ht_3d,iht_above_3d,iht_below_3d,nhgt, &
           GIP_switches)


!   write(6,*) '***************** HYDROGEN GRID AT 400KM **INSIDE POLAR*********************'
!   write(6,1888) Hyd_grid_400km_m3 / 1.e11
!   1888 format(25f4.1)
!   write(6,*) '****************************************************************************'

      do 100 l=1,N_longitudes
        !write(6,*) 'longitude ',l

        geo_local_time_degrees(l) = geo_grid_longitudes_degrees(l) + (universal_time_seconds - 43200.)/240.0
        IF ( geo_local_time_degrees(l).GE.360.0 ) geo_local_time_degrees(l) = geo_local_time_degrees(l) - 360.0
        IF ( geo_local_time_degrees(l).lt.0.0 ) geo_local_time_degrees(l) = geo_local_time_degrees(l) + 360.0
        ssa = geo_local_time_degrees(l)
        rlt = 180.0 + geo_local_time_degrees(l)
        IF ( rlt.GT.360.0 ) rlt = rlt - 360.0
        rlt = rlt*DTR

        geo_local_time_hours(l) = (universal_time_seconds / 3600.) + (geo_grid_longitudes_degrees(l) / 15.)
        if ( geo_local_time_hours(l) .gt. 24.) geo_local_time_hours(l) = geo_local_time_hours(l) - 24.
        if ( geo_local_time_hours(l) .lt. 0.) geo_local_time_hours(l) = geo_local_time_hours(l) + 24.

        do ilt = 1 , 25
         if (hyd_grid_local_times(ilt) .gt. geo_local_time_hours(l)) then
           ilt_east = ilt
           ilt_west = ilt - 1
          goto 1822
         endif
        enddo

1822    factor_lt = geo_local_time_hours(l) - hyd_grid_local_times(ilt_west)

!       write(6,*) 'longitude ',l , geo_local_time_hours(l) , ilt_east , ilt_west , factor_lt

          do 200 m=3,N_Latitudes-2
             !write(6,*) 'latitude ',m
          !g
!               IF (m > mhigh(l)-2 .OR. m < mlow(l)+2) THEN
              IF (m > mhigh(l)-7 .OR. m < mlow(l)+7) THEN
              !g
                  if(m < 45) then
                      ihemi = 2
                  else
                      ihemi = 1
                  endif

                  do ilat = 1 , 19
                    if (hyd_grid_lats(ilat) .gt. geo_grid_latitudes_degrees(m)) then
                    ilat_north = ilat
                    ilat_south = ilat - 1
                    goto 1823
                    endif
                  enddo
1823              factor_lat = ( geo_grid_latitudes_degrees(m) - hyd_grid_lats(ilat_south) ) /10.
!                  write(6,*) 'factor lat ',m , ilat_north , ilat_south , factor_lat
                  H_NE = Hyd_grid_400km_m3(ilt_east,ilat_north)
                  H_NW = Hyd_grid_400km_m3(ilt_west,ilat_north)
                  H_SE = Hyd_grid_400km_m3(ilt_east,ilat_south)
                  H_SW = Hyd_grid_400km_m3(ilt_west,ilat_south)

                  H_N = ((H_NE - H_NW) * factor_lt) + H_NW
                  H_S = ((H_SE - H_SW) * factor_lt) + H_SW

                  H_400 = ((H_N - H_S) * factor_lat) + H_S
!                 write(168,*) geo_local_time_hours(l),geo_grid_latitudes_degrees(m),H_400

                  rlat = geo_grid_latitudes_radians(m)
                  solar_zenith_angle_radians(m,l) =  &
                            ACOS(-COS(rlat)*COS(solar_declination_angle_degrees*DTR)*COS(rlt)+SIN(rlat) &
                           *SIN(solar_declination_angle_degrees*DTR))
                  coschi = COS(solar_zenith_angle_radians(m,l))

                  O_density_1D(1:nhgt) = o_density_high_lat_ions(1:nhgt,m,l)
                  O2_density_1D(1:nhgt) = o2_density_high_lat_ions(1:nhgt,m,l)
                  N2_density_1D(1:nhgt) = n2_density_high_lat_ions(1:nhgt,m,l)
                  NO_density_1D(1:nhgt) = NO_density_high_lat_ions(1:nhgt,m,l)
                  N4S_density_1D(1:nhgt) = N4S_density_high_lat_ions(1:nhgt,m,l)
                  N2D_density_1D(1:nhgt) = N2D_density_high_lat_ions(1:nhgt,m,l)
                  Vx_1d(1:nhgt) = Vx_high_lat_ions(1:nhgt,m,l)
                  Vy_1D(1:nhgt) = Vy_high_lat_ions(1:nhgt,m,l)
                  wvz_1D(1:nhgt) = wvz_high_lat_ions(1:nhgt,m,l)
                  tn_1D(1:nhgt) = tn_high_lat_ions(1:nhgt,m,l)
                  qion_1d(1:nhgt) = qion3d_high_lat_ions(1:nhgt,m,l)
                  qo2p_aurora_1d(1:nhgt) = qo2p_aurora_high_lat_ions(1:nhgt,m,l)
                  qop_aurora_1d(1:nhgt) = qop_aurora_high_lat_ions(1:nhgt,m,l)
                  qn2p_aurora_1d(1:nhgt) = qn2p_aurora_high_lat_ions(1:nhgt,m,l)
                  qnp_aurora_1d(1:nhgt) = qnp_aurora_high_lat_ions(1:nhgt,m,l)
                  qtef_aurora_1d(1:nhgt) = qtef_aurora_high_lat_ions(1:nhgt,m,l)

                  CALL HL__HIGH_LAT_IONS( &
                  GIP_switches, &
                  qion_1d, &
                  qo2p_aurora_1d, qop_aurora_1d,&
                  qn2p_aurora_1d, qnp_aurora_1d,qtef_aurora_1d,&
                  high_lat_time_step_seconds, &
                  m,l,ssa,angdif(m,l), &
                  geo_grid_colatitudes_radians(m),geo_grid_longitudes_radians(l),btotal(m,l), &
                  dip_angle_radians(m,l),solar_zenith_angle_radians(m,l), &
                  ihemi,f107,hprof,dth_radians, &
                  ti1_1D_output, &
                  ex2d,ey2d, &
                  v13d,v23d,tarea, &
                  d13d,d23d,gravity,nhgt,hz_metres, &
                  coschi,hnew,h2,h2diff,hprod2, &
                  hprod3,zkm,l300,l1000,r0, &
                  geo_grid_longitudes_degrees, &
                  declination,Solar_Declination_Angle_degrees,km_plasma, &
                  universal_time_seconds,ioutput_counter,ioutput_high_res_counter, &
                  ioutput_frequency_mins,ioutput_high_res_frequency_mins, &
                  file_res, &
                  O_density_1D,O2_density_1D,N2_density_1D, &
                  NO_density_1D, &
                  N4S_density_1D,N2D_density_1D, &
                  Vx_1d,Vy_1D,wvz_1D,tn_1D, &
                  H_400, &
                  Oplus_density_1d,Hplus_density_1d,NOplus_density_1d,O2plus_density_1d, &
                  N2plus_density_1d,Nplus_density_1d, &
                  Te_1d,Ti1_1d,Ti2_1d,iday_number)
              !g
              !g
              !g

              ENDIF

           Oplus_density_high_lat_ions(1:nhgt,m,l) = Oplus_density_1d(1:nhgt)
           Hplus_density_high_lat_ions(1:nhgt,m,l) = Hplus_density_1d(1:nhgt)
           NOplus_density_high_lat_ions(1:nhgt,m,l) = NOplus_density_1d(1:nhgt)
           O2plus_density_high_lat_ions(1:nhgt,m,l) = O2plus_density_1d(1:nhgt)
           N2plus_density_high_lat_ions(1:nhgt,m,l) = N2plus_density_1d(1:nhgt)
           Nplus_density_high_lat_ions(1:nhgt,m,l) = Nplus_density_1d(1:nhgt)
           Te_high_lat_ions(1:nhgt,m,l) = Te_1d(1:nhgt)
           Ti1_high_lat_ions(1:nhgt,m,l) = Ti1_1d(1:nhgt)
           Ti2_high_lat_ions(1:nhgt,m,l) = Ti2_1d(1:nhgt)

          200 ENDDO
      100 ENDDO

! cg
! cg ensure that we have the polar values properly defined.....
! cg
! cg First the north pole.....

           oplus_polar(1:nhgt) = 0.0
           hplus_polar(1:nhgt) = 0.0
           noplus_polar(1:nhgt) = 0.0
           o2plus_polar(1:nhgt) = 0.0
           n2plus_polar(1:nhgt) = 0.0
           nplus_polar(1:nhgt) = 0.0
           Te_polar(1:nhgt) = 0.0
           Ti1_polar(1:nhgt) = 0.0
           Ti2_polar(1:nhgt) = 0.0
      do 150 l=1,N_longitudes
           oplus_polar(1:nhgt) = oplus_polar(1:nhgt) + Oplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           hplus_polar(1:nhgt) = hplus_polar(1:nhgt) + hplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           noplus_polar(1:nhgt) = noplus_polar(1:nhgt) + noplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           o2plus_polar(1:nhgt) = o2plus_polar(1:nhgt) + o2plus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           n2plus_polar(1:nhgt) = n2plus_polar(1:nhgt) + n2plus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           nplus_polar(1:nhgt) = nplus_polar(1:nhgt) + nplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           Te_polar(1:nhgt) = Te_polar(1:nhgt) + Te_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           Ti1_polar(1:nhgt) = Ti1_polar(1:nhgt) + Ti1_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
           Ti2_polar(1:nhgt) = Ti2_polar(1:nhgt) + Ti2_high_lat_ions(1:nhgt,N_Latitudes-2,l) / float(N_longitudes)
150    enddo
      do 160 l=1,N_longitudes
           Oplus_density_high_lat_ions(1:nhgt,N_Latitudes,l) = oplus_polar(1:nhgt)
           hplus_density_high_lat_ions(1:nhgt,N_Latitudes,l) = hplus_polar(1:nhgt)
           noplus_density_high_lat_ions(1:nhgt,N_Latitudes,l) = noplus_polar(1:nhgt)
           o2plus_density_high_lat_ions(1:nhgt,N_Latitudes,l) = o2plus_polar(1:nhgt)
           n2plus_density_high_lat_ions(1:nhgt,N_Latitudes,l) = n2plus_polar(1:nhgt)
           nplus_density_high_lat_ions(1:nhgt,N_Latitudes,l) = nplus_polar(1:nhgt)
           Te_high_lat_ions(1:nhgt,N_Latitudes,l) = Te_polar(1:nhgt)
           Ti1_high_lat_ions(1:nhgt,N_Latitudes,l) = Ti1_polar(1:nhgt)
           Ti2_high_lat_ions(1:nhgt,N_Latitudes,l) = Ti2_polar(1:nhgt)
160    enddo
      do 170 l=1,N_longitudes
           Oplus_density_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (Oplus_density_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + Oplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           hplus_density_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (hplus_density_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + hplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           noplus_density_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (noplus_density_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + noplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           o2plus_density_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (o2plus_density_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + o2plus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           n2plus_density_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (n2plus_density_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + n2plus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           nplus_density_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (nplus_density_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + nplus_density_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           Te_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (Te_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + Te_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           Ti1_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (Ti1_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + Ti1_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
           Ti2_high_lat_ions(1:nhgt,N_Latitudes-1,l) = (Ti2_high_lat_ions(1:nhgt,N_Latitudes,l)  &
                                                          + Ti2_high_lat_ions(1:nhgt,N_Latitudes-2,l)) / 2.0
170    enddo

! cg and the south north pole.....

           oplus_polar(1:nhgt) = 0.0
           hplus_polar(1:nhgt) = 0.0
           noplus_polar(1:nhgt) = 0.0
           o2plus_polar(1:nhgt) = 0.0
           n2plus_polar(1:nhgt) = 0.0
           nplus_polar(1:nhgt) = 0.0
           Te_polar(1:nhgt) = 0.0
           Ti1_polar(1:nhgt) = 0.0
           Ti2_polar(1:nhgt) = 0.0
      do 250 l=1,N_longitudes
           oplus_polar(1:nhgt) = oplus_polar(1:nhgt) + Oplus_density_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           hplus_polar(1:nhgt) = hplus_polar(1:nhgt) + hplus_density_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           noplus_polar(1:nhgt) = noplus_polar(1:nhgt) + noplus_density_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           o2plus_polar(1:nhgt) = o2plus_polar(1:nhgt) + o2plus_density_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           n2plus_polar(1:nhgt) = n2plus_polar(1:nhgt) + n2plus_density_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           nplus_polar(1:nhgt) = nplus_polar(1:nhgt) + nplus_density_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           Te_polar(1:nhgt) = Te_polar(1:nhgt) + Te_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           Ti1_polar(1:nhgt) = Ti1_polar(1:nhgt) + Ti1_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
           Ti2_polar(1:nhgt) = Ti2_polar(1:nhgt) + Ti2_high_lat_ions(1:nhgt,3,l) / float(N_longitudes)
250    enddo
      do 260 l=1,N_longitudes
           Oplus_density_high_lat_ions(1:nhgt,1,l) = oplus_polar(1:nhgt)
           hplus_density_high_lat_ions(1:nhgt,1,l) = hplus_polar(1:nhgt)
           noplus_density_high_lat_ions(1:nhgt,1,l) = noplus_polar(1:nhgt)
           o2plus_density_high_lat_ions(1:nhgt,1,l) = o2plus_polar(1:nhgt)
           n2plus_density_high_lat_ions(1:nhgt,1,l) = n2plus_polar(1:nhgt)
           nplus_density_high_lat_ions(1:nhgt,1,l) = nplus_polar(1:nhgt)
           Te_high_lat_ions(1:nhgt,1,l) = Te_polar(1:nhgt)
           Ti1_high_lat_ions(1:nhgt,1,l) = Ti1_polar(1:nhgt)
           Ti2_high_lat_ions(1:nhgt,1,l) = Ti2_polar(1:nhgt)
260    enddo
      do 270 l=1,N_longitudes
           Oplus_density_high_lat_ions(1:nhgt,2,l) = (Oplus_density_high_lat_ions(1:nhgt,1,l)  &
                                                          + Oplus_density_high_lat_ions(1:nhgt,3,l)) / 2.0
           hplus_density_high_lat_ions(1:nhgt,2,l) = (hplus_density_high_lat_ions(1:nhgt,1,l)  &
                                                          + hplus_density_high_lat_ions(1:nhgt,3,l)) / 2.0
           noplus_density_high_lat_ions(1:nhgt,2,l) = (noplus_density_high_lat_ions(1:nhgt,1,l)  &
                                                          + noplus_density_high_lat_ions(1:nhgt,3,l)) / 2.0
           o2plus_density_high_lat_ions(1:nhgt,2,l) = (o2plus_density_high_lat_ions(1:nhgt,1,l)  &
                                                          + o2plus_density_high_lat_ions(1:nhgt,3,l)) / 2.0
           n2plus_density_high_lat_ions(1:nhgt,2,l) = (n2plus_density_high_lat_ions(1:nhgt,1,l)  &
                                                          + n2plus_density_high_lat_ions(1:nhgt,3,l)) / 2.0
           nplus_density_high_lat_ions(1:nhgt,2,l) = (nplus_density_high_lat_ions(1:nhgt,1,l)  &
                                                          + nplus_density_high_lat_ions(1:nhgt,3,l)) / 2.0
           Te_high_lat_ions(1:nhgt,2,l) = (Te_high_lat_ions(1:nhgt,1,l)  &
                                                          + Te_high_lat_ions(1:nhgt,3,l)) / 2.0
           Ti1_high_lat_ions(1:nhgt,2,l) = (Ti1_high_lat_ions(1:nhgt,1,l)  &
                                                          + Ti1_high_lat_ions(1:nhgt,3,l)) / 2.0
           Ti2_high_lat_ions(1:nhgt,2,l) = (Ti2_high_lat_ions(1:nhgt,1,l)  &
                                                          + Ti2_high_lat_ions(1:nhgt,3,l)) / 2.0
270    enddo


      CALL HL__ITRPD1D2(v13d,v23d,d13d,d23d,nhgt)


call INTERFACE__high_lat_ions_to_FIXED_GEO ( &
           Oplus_density_high_lat_ions,Hplus_density_high_lat_ions,NOplus_density_high_lat_ions, &
           O2plus_density_high_lat_ions,N2plus_density_high_lat_ions,Nplus_density_high_lat_ions, &
           Te_high_lat_ions,Ti1_high_lat_ions,Ti2_high_lat_ions, &
           Oplus_density_fixed_ht,Hplus_density_fixed_ht,NOplus_density_fixed_ht, &
           O2plus_density_fixed_ht,N2plus_density_fixed_ht,Nplus_density_fixed_ht, &
           Te_fixed_ht,Ti1_fixed_ht,Ti2_fixed_ht, &
           factor_ht_inverse_3d,iht_above_inverse_3d,iht_below_inverse_3d)

  endif  ! initialisation call endif
  Return



end SUBROUTINE HL__POLAR_IONOSPHERE





SUBROUTINE HL__HIGH_LAT_IONS( &
  GIP_switches, &
  qion_1d, &
  qo2p_aurora_1d, qop_aurora_1d ,&
  qn2p_aurora_1d, qnp_aurora_1d,qtef_aurora_1d, &
  high_lat_time_step_seconds,MNEut,LNEut,SSA, &
  ANGdif,geo_colatitude_radians,geo_longitude_radians, &
  BTH,DIP,CHI, &
  IHEmi,F107,HPRof,dth_radians, &
  TI1_dum, &
  EX2d,EY2d,V13d,V23d,TARea, &
  D13d,D23d, &
  Gravity,NHGt,hz_metres,COSchi,HNEw,H2, &
  H2Diff,HPRod2,HPRod3,ZKM,L300,L1000, &
  R0, &
  geo_grid_longitudes_degrees, &
  declination,Solar_Declination_Angle_degrees,km,sec,ioutput, &
  iout_high,ipint,ipint_high,file_res, &
  O_density_1D,O2_density_1D,N2_density_1D, &
  NO_density_1D, &
  N4S_density_1D,N2D_density_1D, &
  Vx_1d,Vy_1D,wvz_1D,tn_1D, &
  H_400, &
  Oplus_density_1d,Hplus_density_1d,NOplus_density_1d,O2plus_density_1d, &
  N2plus_density_1d,Nplus_density_1d, &
  Te_1d,Ti1_1d,Ti2_1d,iday_number)

  IMPLICIT NONE

  INTEGER N_Pressure_Levels,N_Latitudes,N_longitudes
      PARAMETER(N_Pressure_Levels=15)
      PARAMETER(N_Latitudes=91)
      PARAMETER(N_Longitudes=20)

  integer :: iday_number
  integer :: ioutprod,i_ion , j
  integer :: i_molecular_ions_failed
  integer :: ifail_HL__ION_TEMP_HIGH_LAT
  logical :: GIP_switches(20)
  logical :: sw_input_Auroral_production_is_single_overall_rate

  REAL(kind=8) :: &
  AKN2 , AKO , AKO2 , ANGdif , aurqn2 , aurqo , aurqo2 , &
  beta2(90) , BTH , cf21(90) , cf1n(90) , CHI , geo_colatitude_degrees , COSchi , &
  cosi ,  D13d , geo_latitude_degrees , geo_longitude_degrees
  REAL(kind=8) :: &
  D23d , DIP , dmol1d , dte(90) , high_lat_time_step_seconds , dth_radians , &
  dti(90) , DTR , eld 
  REAL(kind=8) :: EX2d , EY2d , F107 , Gravity , gsin(90) , atomic_hydrogen_density(90),&
    atomic_hydrogen_density2(90), &
    ha , hb , hbeta , &
  hbetasub(90) , height , hp , HPLus , hprod , &
  hprodsub(90) , HPRof(20,19)
  REAL(kind=8) :: atomic_hydrogen_density_new(90)
  REAL(kind=8) :: hz_metres(90) , altitude_metres(90) , altitude_km(90) , &
  o , o2 , &
  o2prod , oa , ob , &
  obeta , op , oprod , phi , &
  geo_longitude_radians , &
  PHI_stepped_back_radians , PI , &
  qion_1d(90) , qout(90) , ramin , rincf,&
  qo2p_aurora_1d(90), qop_aurora_1d(90) ,&
  qn2p_aurora_1d(90), qnp_aurora_1d(90),qtef_aurora_1d(90)

  REAL(kind=8) :: &
  rk1 , rk2 , RTD , rteff , rtti(90) , rttin(90) , &
  sini , SSA , TARea(90) , &
  geo_colatitude_radians , thcorr(90)
  REAL(kind=8) :: &
  th_stepped_back_radians , ti , tin , tinlog(90) , topvh , &
  topvo , u2dif , ucos , V13d , V23d , VHTOP
  REAL(kind=8) :: vi , vix(90) , viy(90) , viz(90) , vrsini ,  &
  ZKM(90)
  REAL(kind=8) :: HNEw(90) , H2(90) , H2Diff(90) , HPRod2(90) , HPRod3(90)
  REAL(kind=8) :: AGR , &
  AGT , AGP , R0 , &
  declination , sin_decl , cos_decl ,sdecl , cdecl , &
  solar_Declination_Angle_degrees, &
  time(90), &
  peuvn(90,6),gcol(90), &
  km(6),gr(90),sec , chid(90) , &
  xtime(90) , rx , ry
  REAL(kind=8) H_400
  INTEGER :: i , IHEmi , ihghtout , il1 , iouty , ix1 , ix2 , &
  iy1 , iy2 , k , l , L1000 , L300 , len , lmaxm ,ioutput , &
  ioutputnow , iout_high , ipint_high , file_res(27)
  INTEGER :: LNEut , men , MNEut , NHGt , istop , i_in , ipint

! *********************************************
! *                                           *
! *         sheffield high-latitude           *
! *            ionosphere model.              *
! *                                           *
! *********************************************

  PARAMETER (PI=3.14159,AKO=1.8551E-3,AKO2=3.7101E-3,AKN2=3.2464E-3, &
  DTR=PI/180.0,RTD=180.0/PI)

  LOGICAL :: sw_use_EUVAC_solar_spectrum
  REAL(kind=8) :: nit , n2prod 
  REAL(kind=8) :: &
  ti1_dum(15) 
  DIMENSION Gravity(90) , nit(90) , o(90) , o2(90) , &
  rk1(90) , rk2(90) , oa(90) , ob(90) &
  , obeta(90) , oprod(90) , ha(90) , hb(90) , hbeta(90) , &
  hprod(90) , rincf(90) , rteff(90) , &
  ramin(90) , vi(3) , &
  dmol1d(90) , o2prod(90) , n2prod(90) , &
  D13d(90,91,20,2) , D23d(90,91,20,2) , &
  tin(15) 
  DIMENSION &
  eld(90) , hp(90) , &
  op(90) 

  DIMENSION &
  EX2d(91,20) , EY2d(91,20) , V13d(90,91,20,2) , &
  V23d(90,91,20,2)
  DIMENSION aurqo(90) , aurqo2(90) , aurqn2(90)
  DIMENSION declination(91,20)
  DIMENSION u2dif(90) , ucos(90)
  REAL(kind=8) :: O_density_1D(90)
  REAL(kind=8) :: O2_density_1D(90)
  REAL(kind=8) :: N2_density_1D(90)
  REAL(kind=8) :: NO_density_1D(90)
  REAL(kind=8) :: N4S_density_1D(90)
  REAL(kind=8) :: N2D_density_1D(90)
  REAL(kind=8) :: Vx_1d(90)
  REAL(kind=8) :: Vy_1D(90)
  REAL(kind=8) :: wvz_1D(90)
  REAL(kind=8) :: tn_1D(90)
  REAL(kind=8) :: Oplus_density_1d(90)
  REAL(kind=8) :: Hplus_density_1d(90)
  REAL(kind=8) :: Oplus_density_1d_back(90)
  REAL(kind=8) :: Hplus_density_1d_back(90)
  REAL(kind=8) :: Oplus_velocity_1d(90)
  REAL(kind=8) :: Hplus_velocity_1d(90)
  REAL(kind=8) :: Oplus_velocity_1d_back(90)
  REAL(kind=8) :: Hplus_velocity_1d_back(90)
  REAL(kind=8) :: NOplus_density_1d(90)
  REAL(kind=8) :: O2plus_density_1d(90)
  REAL(kind=8) :: N2plus_density_1d(90)
  REAL(kind=8) :: Nplus_density_1d(90)
  REAL(kind=8) :: Te_1d(90)
  REAL(kind=8) :: Ti1_1d(90)
  REAL(kind=8) :: Ti2_1d(90)
  REAL(kind=8) :: Electron_density_1d(90)
  REAL(kind=8) :: cos_sza_1d(90)
  REAL(kind=8) :: eccentric

  REAL(kind=8) :: auroral_production_rate_Oplus(90)
  REAL(kind=8) :: auroral_production_rate_N2plus(90)
  REAL(kind=8) :: auroral_production_rate_O2plus(90)
  REAL(kind=8) :: euv_production_rate_Oplus(90)
  REAL(kind=8) :: euv_production_rate_N2plus(90)
  REAL(kind=8) :: euv_production_rate_O2plus(90)
  REAL(kind=8) :: total_production_rate_Oplus(90)
  REAL(kind=8) :: total_production_rate_N2plus(90)
  REAL(kind=8) :: total_production_rate_O2plus(90)
  REAL(kind=8) :: geo_grid_longitudes_degrees(20)

  INTEGER :: i_diagnose
  i_diagnose = 0 
  if(i_diagnose == 1) write(6,*) ' mneut lneut ', mneut,lneut

  sw_input_Auroral_production_is_single_overall_rate = GIP_switches(7)

           sw_use_EUVAC_solar_spectrum = .FALSE.  ! ISSUE : this isn't being passed through correctly
                                                  ! at present


  altitude_metres(:) = 0.0   ! just an initialisation

! call the path routine step which steps back one time step
! along the ion drift velocity in the geographic frame.
  call HL__calculate_altitudes (NHGT,hz_metres,dip,altitude_metres)

  altitude_km(:) = altitude_metres(:) / 1000.

  CALL HL__STEP_BACKWARDS(geo_colatitude_radians,geo_longitude_radians,IHEmi, &
            th_stepped_back_radians,PHI_stepped_back_radians, &
            EX2d(MNEut,LNEut), EY2d(MNEut,LNEut),BTH,DIP,high_lat_time_step_seconds,vi,ANGdif)

if (i_diagnose == 1) write(6,*) 'here 2'

  CALL HL__PREVIOUS_POSITION_CALC_INDEXES(th_stepped_back_radians,PHI_stepped_back_radians,dth_radians,ry,rx,iy1,iy2,ix1,ix2)

  DO 100 k = 1 , NHGt

      Oplus_density_1d_back(k) = D13d(k,ix1,iy1,1)*(1.-rx)*(1.-ry) + D13d(k,ix2,iy1,1) &
      *rx*(1.-ry) + D13d(k,ix1,iy2,1)*(1.-rx) &
      *ry + D13d(k,ix2,iy2,1)*rx*ry

      Hplus_density_1d_back(k) = D23d(k,ix1,iy1,1)*(1.-rx)*(1.-ry) + D23d(k,ix2,iy1,1) &
      *rx*(1.-ry) + D23d(k,ix1,iy2,1)*(1.-rx) &
      *ry + D23d(k,ix2,iy2,1)*rx*ry

  100 ENDDO

  DO 200 k = 1 , NHGt

      Oplus_velocity_1d_back(k) = V13d(k,ix1,iy1,1)*(1.-rx)*(1.-ry) + V13d(k,ix2,iy1,1) &
      *rx*(1.-ry) + V13d(k,ix1,iy2,1)*(1.-rx) &
      *ry + V13d(k,ix2,iy2,1)*rx*ry

      Hplus_velocity_1d_back(k) = V23d(k,ix1,iy1,1)*(1.-rx)*(1.-ry) + V23d(k,ix2,iy1,1) &
      *rx*(1.-ry) + V23d(k,ix1,iy2,1)*(1.-rx) &
      *ry + V23d(k,ix2,iy2,1)*rx*ry

  200 ENDDO


if (i_diagnose == 1) write(6,*) 'here 3'


  geo_colatitude_degrees = geo_colatitude_radians/dtr
  geo_latitude_degrees = 90. - geo_colatitude_degrees
  geo_longitude_degrees = geo_longitude_radians/dtr

 do i = 1 , nhgt
 CALL HL_ML__HYDEQ_BASIC(H_400,altitude_km(i),tn_1d(i),gravity(i), &
                atomic_hydrogen_density_new(i))
       
 atomic_hydrogen_density(i) = atomic_hydrogen_density_new(i)

3421 format(3i4,f9.0,2e12.4)
 enddo

!
! these were interpolated from the pressure levels....
! need to sort this out...
! set to zero for now...
!
  do k = 1 , nhgt
    qout(k) = 0.0
  enddo

if (i_diagnose == 1) write(6,*) 'here 4'
!c

  DO 400 k = 1 , NHGt
!g
!g  ...need to fix this - need the Ti from the previous time step....
!g  ... in the meantime just use tn....
!g      rteff(k) = (tn_1D(k)+TI3d(k,MNEut,LNEut))/2.
      rteff(k) = tn_1D(k)
      IF ( rteff(k) < tn_1D(k) ) rteff(k) = tn_1D(k)
  400 ENDDO
 
! ISSUE : the ion_neutral collisions routine needs the NO+ and O2+
!         from the previous timestep - but we don't have this at present.
!         So, I am just initialising the arrays for now.....

  noplus_density_1d(:) = 0.0
  o2plus_density_1d(:) = 0.0

  CALL HL__ION_NEUTRAL_COLLISIONS(o_density_1D,o2_density_1D,n2_density_1d, &
                              Oplus_density_1d_back,noplus_density_1d,o2plus_density_1d,rteff,rincf,ramin,NHGt,1)

  CALL HL__IDRIFT(vx_1d,vy_1d,rincf,EX2d(MNEut,LNEut),EY2d(MNEut,LNEut),ramin,vix,viy,viz,DIP,BTH,NHGt)

if (i_diagnose == 1) write(6,*) 'here 5'


! wind differences

  sini = ABS(SIN(DIP))
  cosi = COS(DIP)
  sin_decl = sin(0.0-declination(mneut,lneut)*dtr)
  cos_decl = cos(0.0-declination(mneut,lneut)*dtr)
  IF ( sini < 0.18 ) THEN
      sini = 0.18
      cosi = 0.983666
  ENDIF

  CALL HL__WINDIF(vx_1d,vy_1d,wvz_1d,DIP,vix,viy,viz,NHGt,cosi,u2dif, &
  ucos,AGR,AGT,AGP,sin_decl,cos_decl,sini)

!c  **
!c  calculate solar ionization rates from reference spectrum.
!c  **

  iouty = 0
  do i=1,nhgt
      gcol(i) = geo_colatitude_radians
      gr(i) = (1.0E5+(hz_metres(i)-1.0E5)*sini)+r0
      time(i) = SEC + geo_longitude_radians*240./dtr
      IF ( time(i) > 86400. ) THEN
          time(i) = time(i) - 86400.
      ELSEIF ( time(i) < 0. ) THEN
          time(i) = time(i) + 86400.
      ENDIF
      xtime(i) = time(i)/3600.
  enddo
if (i_diagnose == 1) write(6,*) 'here 8'

  if(sw_use_EUVAC_solar_spectrum) then
! write(6,*) ' high_lat switch 1 ',sw_use_EUVAC_solar_spectrum
!*****************************
! ALD March 09: Put HL_ML__EUV_ION_PRODUCTION2 call in here:
! This routine will use TIEGCM EUVAC solar spectrum to calculate major
! species photoionisation rates.

  eccentric=(1.0+0.0167*COS(2.0*pi*(iday_number-3)/366.0))**2

  CALL HL_ML__EUV_ION_PRODUCTION_2(90,1,nhgt,o_density_1d,o2_density_1d, &
                 n2_density_1d, Solar_Declination_Angle_degrees, &
                 time, f107, gcol, R0, tn_1d, gravity, km, gr,   &
                 eccentric, chid, peuvn)

  else
! write(6,*) ' high_lat switch 2 ',sw_use_EUVAC_solar_spectrum

  CALL HL_ML__EUV_ION_PRODUCTION(90,1,nhgt,o_density_1d,o2_density_1d, &
                     n2_density_1d,                             &
                     Solar_Declination_Angle_degrees,time,f107, &
                     peuvn,gcol,r0,tn_1d,gravity, &
                     km,gr,chid,0)

! *************************
  endif

if (i_diagnose == 1) write(6,*) 'here 9'

  do i=1,nhgt
      euv_production_rate_Oplus(i)=peuvn(i,1)
      euv_production_rate_N2plus(i)=peuvn(i,4)
      euv_production_rate_O2plus(i)=peuvn(i,5)
  enddo

! calculate auroral ionisation rates for O+, N2+ and O2+....

if (sw_input_Auroral_production_is_single_overall_rate) then
  CALL HL__PARTITION_AURORAL_IONISATION(qion_1d,o_density_1d,o2_density_1d,n2_density_1d,NHGt, &
  auroral_production_rate_Oplus,auroral_production_rate_N2plus,auroral_production_rate_O2plus)
else
  auroral_production_rate_Oplus  = qop_aurora_1d*1.e6
  auroral_production_rate_N2plus = qn2p_aurora_1d*1.e6
  auroral_production_rate_O2plus = qo2p_aurora_1d*1.e6
endif

  do i=1,nhgt
      total_production_rate_Oplus(i) = euv_production_rate_Oplus(i) + auroral_production_rate_Oplus(i)
      total_production_rate_N2plus(i) = euv_production_rate_N2plus(i) + auroral_production_rate_N2plus(i)
      total_production_rate_O2plus(i) = euv_production_rate_O2plus(i) + auroral_production_rate_O2plus(i)
  enddo

if (i_diagnose == 1) write(6,*) 'here 10'
  phi = SSA*DTR
  IF ( SSA > 180. ) phi = (SSA-360.)*DTR

  CALL HL__ETEMP(phi,CHI,tn_1D,te_1d,NHGt,ZKM,L300,L1000)

  do i=1,nhgt
      if(Oplus_density_1d_back(i) < 1.e-5) Oplus_density_1d_back(i) = 1.e-5
      if(Hplus_density_1d_back(i) < 1.e-5) Hplus_density_1d_back(i) = 1.e-5
  enddo

  CALL HL__ION_TEMP_HIGH_LAT(tn_1D,te_1d,Oplus_density_1d_back,Hplus_density_1d_back,ti1_1d,cf1n,rttin,tinlog, &
                         n2_density_1d,o_density_1d,o2_density_1d,NHGt,u2dif,ifail_HL__ION_TEMP_HIGH_LAT)

  do i=1,nhgt
      ti2_1d(i) =  ti1_1d(i)
  enddo

  CALL HL__RATEK(tn_1D,ti1_1d,u2dif,rk1,rk2,NHGt)

 if (i_diagnose == 1) write(6,*) 'here 11'

!c  **
!c  check we have the right component here
!c  **
  vrsini = vi(1)*sini
!g
!g vrsini allows for non-vertical flux-tubes (during a time step the
!g plasma moves up or down a bit due to EXB - this needs to be accounted
!g for since we are solving on a fixed height grid)
!g
 

  CALL HL__OCOEFF(Oplus_density_1d_back,Hplus_density_1d_back,Hplus_velocity_1d_back,tn_1D,ti1_1d,te_1d,rk1,rk2,vrsini,oa,ob,obeta, &
  total_production_rate_Oplus,hbeta,dte,dti,hprodsub,hbetasub,gsin,cf21,beta2, &
  thcorr,rtti,cf1n,n2_density_1d,o_density_1d,o2_density_1d,atomic_hydrogen_density, &
  Gravity,NHGt,HNEw, &
  H2,H2Diff, &
  HPRod3,sini,ucos)

if (i_diagnose == 1) write(6,*) 'here 12'
!g
  CALL HL__HCOEFF(Hplus_density_1d_back,Oplus_density_1d_back,Oplus_velocity_1d_back,ti2_1d,te_1d,vrsini,ha,hb,hprod,dte,dti,hprodsub, &
  gsin,cf21,beta2,thcorr,rtti,rttin,tinlog,n2_density_1d,o_density_1d, &
  o2_density_1d,atomic_hydrogen_density, &
  NHGt,HNEw,H2,H2Diff,HPRod3,sini,ucos)

if (i_diagnose == 1) write(6,*) 'here 13'
  topvo = vrsini
!c  **
!g  VHTOP sets the field-aligned velocity of the H+ ions at the top of
!g  the flux-tube as a function of magnetic latitude and longitude.
!g  in VHTOP though the final velcotiy has been set to zero...
!g  so we don't need to call it do we....
!g
! topvh = VHTOP(th,PHM)
  topvh=0.0
!c         topvo=0.0
!g
!g Solve the triadiagonal matrix for O+
!g
  i_ion=0

  Oplus_density_1d(:) = Oplus_density_1d_back(:)
  Oplus_velocity_1d(:) = Oplus_velocity_1d_back(:)

  CALL HL__IONCHNG(oa,ob,obeta,total_production_rate_Oplus,Oplus_density_1d,Oplus_velocity_1d, &
               topvo,high_lat_time_step_seconds,TARea,NHGt,HNEw,H2, &
               HPRod2,sini,i_ion)

if (i_diagnose == 1) write(6,*) 'here 14'
  lmaxm = 0
  DO 600 l = 1 , NHGt
      IF ( Oplus_density_1d(l) <= 0.0 ) lmaxm = l
  600 ENDDO
  IF ( lmaxm /= 0 ) THEN
      if(lmaxm > nhgt) then
          lmaxm=nhgt
          write(6,*) 'lmaxm ',mneut,lneut,lmaxm
      endif
      DO 650 l = 1 , lmaxm
          Oplus_density_1d(l) = total_production_rate_Oplus(l)/obeta(l)
      650 ENDDO
  !g
  !g A latest attempt to get the point above the last negative point to
  !g behave - here goes .......
  !g
      Oplus_density_1d(lmaxm+1)=(Oplus_density_1d(lmaxm)+Oplus_density_1d(lmaxm+2))/2.0
  !g
  !g If we have had to bodge up the O+ density then we should really
  !g do the same with the field-aligned velocity to make it consistent
  !g
  ! do l=2,lmaxm
  ! V1(l) = oA(l) + oB(l)*((D1(l-1)/D1(l))-1.)*SINi/Hnew(l-1)
  ! enddo
  ENDIF

!g
!g Solve the triadiagonal matrix for H+
!g
  i_ion=0

  Hplus_density_1d(:) = Hplus_density_1d_back(:)
  Hplus_velocity_1d(:) = Hplus_velocity_1d_back(:)

  CALL HL__IONCHNG(ha,hb,hbeta,hprod,Hplus_density_1d,Hplus_velocity_1d,topvh,high_lat_time_step_seconds,TARea,NHGt,HNEw,H2, &
  HPRod2,sini,i_ion)

  IF ( lmaxm /= 0 ) THEN
  ! write(6,*) 'lmaxm problem ',lmaxm
      DO 700 l = 1 , lmaxm
      ! write(6,*) 'list ',i,hprod(i),hbeta(i)
          Hplus_density_1d(l) = hprod(l)/hbeta(l)
      700 ENDDO
  ENDIF
  lmaxm = 0
  DO 800 l = 1 , NHGt
      IF ( Hplus_density_1d(l) <= 0.0 ) lmaxm = l
  800 ENDDO
  IF ( lmaxm /= 0 ) THEN
      DO 850 l = 1 , lmaxm
          Hplus_density_1d(l) = hprod(l)/hbeta(l)
          Oplus_density_1d(l) = total_production_rate_Oplus(l)/obeta(l)
      850 ENDDO
  ENDIF

if (i_diagnose == 1) write(6,*) 'here 15'

!g
      D13d(:,MNEut,LNEut,2) = Oplus_density_1d(:)
      D23d(:,MNEut,LNEut,2) = Hplus_density_1d(:)
      V13d(:,MNEut,LNEut,2) = Oplus_velocity_1d(:)
      V23d(:,MNEut,LNEut,2) = Hplus_velocity_1d(:)

if (i_diagnose == 1) write(6,*) 'here 16'

  cos_sza_1d(:) = coschi

  call HL_ML__MOLECULAR_IONS_ON_TUBES( &
  GIP_switches, &
  90,1,nhgt,f107,O_density_1D,O2_density_1D,N2_density_1D, &
  no_density_1d,n4s_density_1d,n2d_density_1d, &
  tn_1d,te_1d,ti1_1d,Oplus_density_1d, &
  altitude_km,cos_sza_1d, &
  euv_production_rate_Oplus,euv_production_rate_N2plus,euv_production_rate_O2plus, &
  auroral_production_rate_Oplus,auroral_production_rate_N2plus,auroral_production_rate_O2plus, &
  N2plus_density_1d,NOplus_density_1d, &
  O2plus_density_1d,Nplus_density_1d, &
  i_molecular_ions_failed)

if (i_diagnose == 1) write(6,*) 'here 18'

  RETURN




end SUBROUTINE HL__HIGH_LAT_IONS




SUBROUTINE HL__TRIDIAG(SUB,DIAg,SUP,RHS,M,ANS)
  IMPLICIT NONE
  REAL(KIND=8) :: ANS , DIAg , RHS , SUB , SUP , x
  INTEGER :: i , i1 , j , M

! we solve a TRIDIAGONAL system whose subscripts run from 2 to m.
! in practice m=n-1.the values of ans-number density-at 1 and n are
! supplied as boundary conditions in the main program.

  DIMENSION SUB(90) , DIAg(90) , SUP(90) , RHS(90) , ANS(90)
  DO 100 i = 3 , M
      i1 = i - 1
  ! write(6,*) 'yaargh ',i1,diag(i1)
      x = SUB(i)/DIAg(i1)
  ! write(6,*) 'diag 1 ',i,diag(i),x,sup(i1)
      DIAg(i) = DIAg(i) - x*SUP(i1)
      if(diag(i) == 0.0) then
          diag(i) = diag(i-1)
          write(6,*) 'using this fix..... '
      endif
  ! write(6,*) 'diag 2 ',i,diag(i)
      RHS(i) = RHS(i) - x*RHS(i1)
  100 ENDDO
! write(6,*) 'yaargh2 ',m,diag(m)
  ANS(M) = RHS(M)/DIAg(M)
  DO 200 i = 1 , M - 2
      j = M - i
  ! write(6,*) 'yaargh3 ',j,diag(j)
      ANS(j) = (RHS(j)-SUP(j)*ANS(j+1))/DIAg(j)
  200 ENDDO
  RETURN



end SUBROUTINE HL__TRIDIAG



SUBROUTINE HL__DIFFARR(ARR,DARr,N,H,H2,H2Diff,HPRod3)
  IMPLICIT NONE
  REAL(KIND=8) :: ARR , DARr , H , H2 , h21 , H2Diff , h2n , HPRod3 , hsum , &
  hsum2
  INTEGER :: i , N

  DIMENSION H(90) , H2(90) , H2Diff(90) , HPRod3(90)
  DIMENSION ARR(90) , DARr(90)
  hsum = H(1) + H(2)
  hsum2 = hsum*hsum
  h21 = H2(1)
  DARr(1) = (-ARR(1)*(hsum2-h21)+ARR(2)*hsum2-ARR(3)*h21)/HPRod3(1)
  h2n = H2(N-1)
  hsum = H(N-1) + H(N-2)
  hsum2 = hsum*hsum
  DARr(N) = (ARR(N-2)*h2n-ARR(N-1)*hsum2+ARR(N)*(hsum2-h2n)) &
  /HPRod3(N-2)
  DO 100 i = 2 , N - 1
      DARr(i) = (ARR(i+1)*H2(i-1)+ARR(i)*H2Diff(i-1)-ARR(i-1)*H2(i)) &
      /HPRod3(i-1)
  100 ENDDO
  RETURN

end SUBROUTINE HL__DIFFARR



SUBROUTINE HL__CFIONS(D1,D2,TI,CF12,CF21,BETa2,RTTi,N)
  IMPLICIT NONE
  REAL(KIND=8) :: &
  BETa2 , CF12 , CF21 , D1 , D2 , ft , rat1 , rat2 , RTTi(90) , &
  TI
  INTEGER :: i , N

! this routine has as inputs the 2 ion number densities and
! ti;we are assuming t1=t2.
! the collision frequencies are returned,and also the
! thermal coefficient beta2.

  DIMENSION D1(90) , D2(90) , TI(90) , CF12(90) , CF21(90) , BETa2(90)

  RTTi(1) = SQRT(TI(1))

  DO 100 i = 2 , N
      RTTi(i) = SQRT(TI(i))
      ft = TI(i)*RTTi(i)
      rat1 = D1(i)/ft
      rat2 = D2(i)/ft
      CF12(i) = 8.4E-8*rat2
      CF21(i) = 1.35E-6*rat1
  !c  **
  !c  modified by tjfr jan 86
  !c  **
      BETa2(i) = 2.6475/(D2(i)/D1(i)+2.34)
  100 ENDDO
  RETURN



end SUBROUTINE HL__CFIONS


SUBROUTINE HL__RATEK(TN,TI,U2Dif,RK1,RK2,Nhgt)
  IMPLICIT NONE
  REAL(kind=8) :: &
  ako , RK1 , RK2 , teff1 , teff2 , TI , TN , tni , tr , tsum , &
  U2Dif , vterm
  INTEGER :: i , Nhgt

! the rate coefficients are as st-maurice and torr-jgr,1978

  DIMENSION TN(90) , TI(90) , RK1(90) , RK2(90) , U2Dif(90)
  DATA ako/1.924E-3/

  DO 100 i = 1 , Nhgt
      vterm = ako*U2Dif(i)/3.
      tni = TN(i)
      tsum = vterm + TI(i) - tni
      teff1 = 0.6364*tsum + tni
      teff2 = 0.6667*tsum + tni
      tr = teff2/300.
      RK2(i) = 2.82E-17 + &
      tr*(-7.74E-18+tr*(1.073E-18+tr*(-5.17E-20+tr*9.65E-22) &
      ))
      tr = teff1/300.
      IF ( teff1 > 1700. ) THEN
          RK1(i) = 2.73E-18 + tr*(-1.155E-18+tr*1.483E-19)
      ELSE
          RK1(i) = 1.533E-18 + tr*(-5.92E-19+tr*8.6E-20)
      ENDIF
  100 ENDDO
  RETURN



end SUBROUTINE HL__RATEK



SUBROUTINE HL__OCOEFF(D1,D2,V2,TN,TI,TE,RK1,RK2,VRSini,OA,OB,OBEta, &
  OPRod,HBEta,DTE,DTI,HPRodsub,HBEtasub,GSIn,CF21, &
  BETa2,THCorr,RTTi,CF1n,NIT,O,O2,atomic_hydrogen,G,Nhgt,HNEw, &
  H2,H2Diff,HPRod3,SINi,UCOs)
  IMPLICIT NONE
  REAL(kind=8) :: ako , akosin , beta12 , BETa2(90) , cf12 , CF1n(90) , &
  cf1tot , CF21(90) , D1 , D2 , dd2 , delta , DTE , DTI , &
  dtot , G(90) , GSIn(90) , atomic_hydrogen(90) , HBEta , HBEtasub(90) , &
  HPRodsub(90) , O(90) , O2(90) , OA , OB , OBEta , OPRod , &
  RK1 , RK2 , RTTi(90) , SINi , TE , THCorr(90) , TI , TN , &
  UCOs(90) , V2
  REAL(kind=8) :: VRSini , wind , y
  INTEGER :: i , Nhgt

  REAL(kind=8) :: NIT(90) , HNEw(90) , H2(90) , H2Diff(90) , HPRod3(90)
  DIMENSION D1(90) , D2(90) , V2(90) , TN(90) , TI(90) , TE(90) , &
  DTI(90) , DTE(90) , RK1(90) , RK2(90) , cf12(90) , &
  OA(90) , OB(90) , OBEta(90) , OPRod(90) , HBEta(90) , &
  dd2(90)
  DATA ako/1.8551E-3/

! the photoionization is given as the input 0q;
! on output oprod includes the effect of o-h charge exchange.
! hbeta is returned as part of the calculation.
! cf1n is passed from itemp,sqrt(ti) is passed from HL__CFIONS.
  CALL HL__DIFFARR(TI,DTI,Nhgt,HNEw,H2,H2Diff,HPRod3)
  CALL HL__DIFFARR(TE,DTE,Nhgt,HNEw,H2,H2Diff,HPRod3)
  DO 100 i = 1 , Nhgt
      GSIn(i) = G(i)*SINi
  100 ENDDO
  akosin = SINi/ako
!g
!g  make sure our ions are not set to zero .....
!g
  do i=1,nhgt
      if(d1(i) < 1.e-5) d1(i) = 1.e-5
      if(d2(i) < 1.e-5) d2(i) = 1.e-5
      if(TI(i) < 10.) ti(i) = 1000.
  enddo

  CALL HL__CFIONS(D1,D2,TI,cf12,CF21,BETa2,RTTi,Nhgt)
  CALL HL__DIFFARR(D2,dd2,Nhgt,HNEw,H2,H2Diff,HPRod3)
  DO 200 i = 2 , Nhgt
      delta = 0.565*BETa2(i)
      beta12 = D2(i)*BETa2(i)/D1(i)
      THCorr(i) = 1. - delta
      cf1tot = CF1n(i) + cf12(i)*THCorr(i)
      dtot = D1(i) + D2(i)
      y = (TE(i)*dd2(i)/dtot+DTE(i)+DTI(i)*(1.-beta12))*akosin
      wind = CF1n(i)*UCOs(i) + VRSini*cf1tot
      OA(i) = (cf12(i)*THCorr(i)*V2(i)+wind-GSIn(i)-y)/cf1tot
      OB(i) = (TI(i)+TE(i)*D1(i)/dtot)/(ako*cf1tot)
  200 ENDDO
  DO 300 i = 1 , Nhgt
      HPRodsub(i) = 2.5E-17*atomic_hydrogen(i)*SQRT(TN(i))
      HBEtasub(i) = 2.3E-17*O(i)
      OBEta(i) = RK1(i)*NIT(i) + RK2(i)*O2(i) + HPRodsub(i)
      HBEta(i) = HBEtasub(i)*RTTi(i)
      OPRod(i) = OPRod(i) + HBEta(i)*D2(i)
  300 ENDDO
  RETURN



end SUBROUTINE HL__OCOEFF



SUBROUTINE HL__HCOEFF(D2,D1,V1,TI,TE,VRSini,HA,HB,HPRod,DTE,DTI, &
  HPRodsub,GSIn,CF21,BETa2,THCorr,RTTi,RTTin, &
  TINlog,NIT,O,O2,atomic_hydrogen,Nhgt,HNEw,H2,H2Diff,HPRod3, &
  SINi,UCOs)
  IMPLICIT NONE

  REAL(kind=8) :: akh , akhsin , BETa2(90) , CF21(90) , cf2n , cf2nk1 , &
  cf2nk2 , cf2tot , D1 , D2 , dd1 , DTE , DTI , dtot , GSIn(90) &
  , atomic_hydrogen(90) , HA , HB , HPRod , HPRodsub(90) , HNEw(90) , H2(90) &
  , H2Diff(90) , HPRod3(90)
  REAL(kind=8) :: O(90) , O2(90) , RTTi(90) , RTTin(90) , SINi , TE , &
  THCorr(90) , TI , TINlog(90) , UCOs(90) , V1 , VRSini , &
  wind , x1 , x2 , y
  INTEGER :: i , Nhgt

  REAL(kind=8) :: NIT(90)
  DIMENSION D2(90) , D1(90) , V1(90) , dd1(90) , TI(90) , DTI(90) , &
  TE(90) , DTE(90) , HA(90) , HB(90) , HPRod(90)
  DATA akh/1.203E-4/

! the necessary functions of t=ti+tn are all calculated
! in itemp and passed across.sqrt(ti) is passed from HL__CFIONS.

  akhsin = SINi/akh
  CALL HL__DIFFARR(D1,dd1,Nhgt,HNEw,H2,H2Diff,HPRod3)
  DO 100 i = 2 , Nhgt
      cf2nk1 = 6.61E-17*O(i)
      cf2nk2 = 2.0E-16*atomic_hydrogen(i)
      x1 = cf2nk1*RTTi(i)*(1.0-.047*log10(TI(i)))**2
      x2 = cf2nk2*RTTin(i)*(1.0-0.082*TINlog(i))**2
      cf2n = x1 + x2 + 3.2E-15*O2(i) + 3.36E-15*NIT(i)
      cf2tot = cf2n + CF21(i)*THCorr(i)
      dtot = D1(i) + D2(i)
      y = (TE(i)*dd1(i)/dtot+DTE(i)+DTI(i)*(1.+BETa2(i)))*akhsin
      wind = cf2n*UCOs(i) + VRSini*cf2tot
      HA(i) = (CF21(i)*THCorr(i)*V1(i)+wind-GSIn(i)-y)/cf2tot
      HB(i) = (TI(i)+TE(i)*D2(i)/dtot)/(akh*cf2tot)
      HPRod(i) = HPRodsub(i)*D1(i)
  100 ENDDO
  HPRod(1) = HPRodsub(1)*D1(1)
  RETURN



end SUBROUTINE HL__HCOEFF



SUBROUTINE HL__IONCHNG(A,B,BETa,Q,DN,V,VTOp,TCHange,TARea,Nhgt,H,H2, &
  HPRod2,SINi,i_ion)
  IMPLICIT NONE
  REAL(kind=8) :: A , aa , B , bb , BETa , diag , DN , H , H2 , HPRod2 , Q , &
  rhs , SINi , sini2 , sub , sup
  REAL(kind=8) :: TARea(90) , TCHange , topratio , V , VTOp , x , y
  INTEGER :: i , i1 , Nhgt , n1 , i_ion

  DIMENSION H(90) , H2(90) , HPRod2(90)
  DIMENSION DN(90) , sub(90) , diag(90) , sup(90) , rhs(90) , &
  A(90) , B(90) , BETa(90) , Q(90) , aa(90) , bb(90) , &
  V(90)
  sini2 = SINi*SINi
  DO 100 i = 2 , Nhgt
      aa(i) = A(i)*SINi
      bb(i) = B(i)*sini2
  100 ENDDO
  DO 200 i = 2 , Nhgt - 1
      x = bb(i)/HPRod2(i-1)
      y = bb(i+1)/H2(i)
      sub(i) = x
      diag(i) = aa(i)/H(i) - x - y*TARea(i) - BETa(i) - 1./TCHange
      sup(i) = (y-aa(i+1)/H(i))*TARea(i)
      rhs(i) = -Q(i) - DN(i)/TCHange
      if(i_ion == 1) then
          write(60,1234) i,aa(i),h(i),x,y,tarea(i),beta(i),sini
      endif
      1234 format(i4,7e12.4)
  200 ENDDO
  DN(1) = Q(1)/BETa(1)
  rhs(2) = rhs(2) - sub(2)*DN(1)
  n1 = Nhgt - 1
  topratio = bb(Nhgt)/(bb(Nhgt)+H(n1)*(VTOp*SINi-aa(Nhgt)))
  diag(n1) = diag(n1) + sup(n1)*topratio
!g
  if (i_ion == 1) then
      write(60,*) ' before tridiag'
      do i=1,90
          if(diag(i) > -0.01) diag(i) = -0.01
          write(60,6822) i,sub(i),diag(i),sup(i),rhs(i)
          6822 format(i4,e12.4,f20.10,2e12.4)
      enddo
  endif
  CALL HL__TRIDIAG(sub,diag,sup,rhs,n1,DN)
  DN(Nhgt) = DN(n1)*topratio

  if (i_ion == 1) then
      write(60,*) ' after tridiag'
      do i=1,90
          write(60,6822) i,dn(i)
      enddo
  endif

  DO 300 i = 2 , Nhgt
      i1 = i - 1
      V(i) = A(i) + B(i)*((DN(i1)/DN(i))-1.)*SINi/H(i1)
  300 ENDDO
  RETURN



end SUBROUTINE HL__IONCHNG



SUBROUTINE HL__WINDIF(VX,VY,VZ,DIP,VIX,VIY,VIZ,Nhgt,COSi,U2Dif,UCOs, &
  AGR,AGT,AGP,sin_decl,cos_decl,Sini)
  IMPLICIT NONE
  REAL(kind=8) :: &
  COSi , DIP , U2Dif(90) , UCOs(90) , udifx , udify , udifz , &
  VIX(90) , VIY(90) , VIZ(90) , VX(90) , VY(90) , VZ(90)
  REAL(kind=8) :: AGR , AGT , AGP , sin_decl , cos_decl , sini
  INTEGER :: j , Nhgt

! calculate wind difference

  DO 100 j = 1 , Nhgt
      udifx = VZ(j) - VIZ(j)
      udify = VX(j) - VIX(j)
      udifz = VY(j) - VIY(j)
      U2Dif(j) = udifx*udifx + udify*udify + udifz*udifz
  100 ENDDO
!c  **
  IF ( DIP > 0.0 ) THEN
      DO 150 j = 1 , Nhgt
          UCOs(j) = (VX(j)*COSi*cos_decl)+(VY(j)*COSi*sin_decl) &
          + (VZ(j)*SINi)
      150 ENDDO
  ELSE
      DO 200 j = 1 , Nhgt
          UCOs(j) = -(VX(j)*COSi*cos_decl)-(VY(j)*COSi*sin_decl) &
          + (VZ(j)*SINi)
      200 ENDDO
  ENDIF
!c  **
  RETURN

end SUBROUTINE HL__WINDIF



SUBROUTINE HL__STEP_BACKWARDS(TH1,PHI1,IHEmi,TH2,PHI2,EX,EY,BTH,DIP,high_lat_time_step_seconds,VI, &
  ANGdif)
  IMPLICIT NONE
  REAL(kind=8) :: ANGdif , BTH , DIP , dph , high_lat_time_step_seconds , dth , EX , EY , PHI1 , &
  PHI2 , PI , R0 , TH1 , TH2 , VI
  INTEGER :: IHEmi
  PARAMETER (PI=3.14159,R0=6.370E06)
  DIMENSION VI(3)
!c  **
!c  th1,phi1 are the geographic lat and long of start point.
!c  th2,phi2 are the geographic lat and long of end point,
!c  after stepping back along ion velocity track for one
!c  time step.
!c  **
!c  ex,ey are southward and eastward electric fields
!c  **
!c  evaluate ion velocities
!c  1=upward
!c  2=southward
!c  3=eastward
!c  **
  VI(1) = (EY*COS(ANGdif)-EX*SIN(ANGdif))*COS(DIP)/BTH
  VI(2) = -EY*SIN(DIP)/BTH
  VI(3) = EX/(BTH*SIN(DIP))

  IF ( IHEmi == 2 ) THEN
      TH1 = PI - TH1
      VI(2) = -VI(2)
  ENDIF

  dth = VI(2)*high_lat_time_step_seconds/R0
  dph = VI(3)*high_lat_time_step_seconds/(R0*SIN(TH1))
  TH2 = TH1 - dth
  PHI2 = PHI1 - dph
  IF ( PHI2 < 0.0 ) PHI2 = PHI2 + 2.*PI
  IF ( PHI2 >= 2.*PI ) PHI2 = PHI2 - 2.*PI
  IF ( IHEmi == 2 ) THEN
      TH2 = PI - TH2
      TH1 = PI - TH1
      VI(2) = -VI(2)
  ENDIF
  RETURN



end SUBROUTINE HL__STEP_BACKWARDS



SUBROUTINE HL__PREVIOUS_POSITION_CALC_INDEXES(TH_stepped_back_radians,PHI_stepped_back_radians,dth_radians,ry,rx,iy1,iy2,ix1,ix2)

  IMPLICIT NONE
  REAL(kind=8) :: TH_stepped_back_radians , PHI_stepped_back_radians , dth_radians , rtd , pi
  REAL(kind=8) :: adif , ry , rx
  INTEGER :: iy , iy1 , iy2 , ix , ix1 , ix2
  PARAMETER (PI=3.14159,RTD=180./PI)

!c  **
!c  find the index for interpolating back into the d and v arrays for
!c  initializing d1 and d2, and v1 and v2.
!c  **
!c  **

   adif = MOD(NINT(PHI_stepped_back_radians*RTD),NINT(360.))
   IF ( adif < 0.0 ) adif = adif + 360.
   ry = adif/18. + 1.E-06
   iy = ry
   ry = ry - iy
   iy1 = iy + 1
   IF ( iy1 == 21 ) iy1 = 1
   iy2 = iy1 + 1
   IF ( iy2 == 21 ) iy2 = 1

  IF ( iy1 < 1 .OR. iy1 > 20 ) WRITE(6,99002) iy1
  99002 FORMAT ('0',10x, &
  'array index out of bounds in high-lat ionosphere,  iy1=', &
  i5)

   rx = 90. - TH_stepped_back_radians/dth_radians
   ix = rx + 1.E-06
   rx = rx - ix
   ix = ix + 1
   ix1 = ix
   ix2 = ix1 + 1

  RETURN



end SUBROUTINE HL__PREVIOUS_POSITION_CALC_INDEXES



SUBROUTINE HL__PARTCLQ(OPRod,O2Prod,N2Prod,Nhgt,AURqo,AURqo2,AURqn2)

  IMPLICIT NONE

  REAL(kind=8) :: AURqn2 , AURqo , AURqo2 , O2Prod , OPRod
  INTEGER :: i , Nhgt

! this routine returns the original ionization profile
! +any additional profile due to particle precipitation.

  REAL(kind=8) :: N2Prod
  DIMENSION AURqo(90) , AURqo2(90) , AURqn2(90)
  DIMENSION OPRod(90) , O2Prod(90) , N2Prod(90)
!*
! ADF  Originally there were Cusp contributions here too - Cuspqo(i),
! Cuspqo2 and Cuspqn2, but these were not set to anything
! ADF
  DO 100 i = 1 , Nhgt
      OPRod(i) = OPRod(i) + AURqo(i)
      O2Prod(i) = O2Prod(i) + AURqo2(i)
      N2Prod(i) = N2Prod(i) + AURqn2(i)
  100 ENDDO
  RETURN



end SUBROUTINE HL__PARTCLQ







SUBROUTINE HL__ETEMP(FYE,SZA,TN,TE,Nhgt,ZKM,L300,L1000)
  IMPLICIT NONE
  REAL(kind=8) :: dte , dted , dten , fh , fl , fract , FYE , grad , rad100 , &
  rad80 , rad90 , SZA , szaday , szanight , TE , tefhd , tefhn
  REAL(kind=8) :: tefld , tefln , teh , tel , TN , zh , ZKM(90) , zl
  INTEGER :: i , L1000 , L300 , Nhgt

! we calculate te as a function of tn and solar zenith angle.

  DIMENSION TN(90) , TE(90)
  DATA rad80 , rad90 , rad100/1.3963 , 1.5708 , 1.7453/ , tefhd , &
  tefld/4.2 , 2.4/ , tefhn , tefln/3.0 , 1.35/ , dted , &
  dten/1.0 , 0.75/
  IF ( SZA <= rad80 ) THEN

  ! daytime conditions

      fh = tefhd
      fl = tefld
      dte = dted
  ELSE
      IF ( SZA <= rad100 ) THEN
          IF ( SZA <= rad90 .OR. FYE >= 0. ) THEN
              IF ( FYE < 0. ) THEN

              ! sunrise conditions for 80<sza<90.

                  szaday = rad80
                  szanight = rad90
              ELSE

              ! sunset conditions for 80<sza<100.

                  szaday = rad80
                  szanight = rad100
              ENDIF

          ! put r=te/tn.r is fixed at night and during the day.for the transitio
          ! periods we make r a linear function of solar zenith angle-
          ! this is slightly different from gb who uses hour-angle.

              fract = (SZA-szaday)/(szanight-szaday)
              fh = tefhd + fract*(tefhn-tefhd)
              fl = tefld + fract*(tefln-tefld)
              dte = dted + fract*(dten-dted)
              GOTO 100
          ENDIF
      ENDIF

  ! nighttime conditions

      fh = tefhn
      fl = tefln
      dte = dten
  ENDIF

! the above has fixed the values of te/tn at 1000 &300 km,and also
! the temperature gradient above 1000 km.
! using the known values of tn we make te a piecewise linear functio
! of z,with the extra condition that te>tn(relevant at low altitudes

  100 zh = ZKM(L1000)
  zl = ZKM(L300)
  teh = fh*TN(L1000)
  tel = fl*TN(L300)
  grad = (teh-tel)/(zh-zl)
  DO 200 i = 1 , L1000
      TE(i) = tel + grad*(ZKM(i)-zl)
      IF ( TE(i) < TN(i) ) TE(i) = TN(i)
  200 ENDDO
  DO 300 i = L1000 + 1 , Nhgt
      TE(i) = teh + dte*(ZKM(i)-zh)
  300 ENDDO
  RETURN



end SUBROUTINE HL__ETEMP






SUBROUTINE HL__ION_TEMP_HIGH_LAT(TN,TE,N1,N2,TI,CF1n,RTTin,TINlog, &
  NIT,O,O2,Nhgt, &
  U2Dif,ifail)
  IMPLICIT NONE
  REAL(kind=8) :: &
  cf1a , CF1n(90) , cf1nk1 , cf1nk2 , cf1nk3 , factor , O(90) , &
  O2(90) , RTTin(90) , rttp , t1 , t2 , TE , TI , TINlog(90) , &
  TN , tp , tplog
  REAL(kind=8) :: u2 , U2Dif(90) , w , w1 , w2 , w3 , w4 , w5 , w6 , w8 , x
  INTEGER :: l , Nhgt , icount , ifail

! for the purposes of this routine it is necessary to work in the cg
! system.we alter the coefficients to express the number densities a
! velocities in the right units.

  REAL(kind=8) :: N1 , N2 , NIT(90)
  DIMENSION TN(90) , TE(90) , N1(90) , N2(90) , TI(90)
!c  **
!c  **
  factor = 1.0
  ifail = 0
!c  **
  DO 100 l = 1 , Nhgt
      u2 = U2Dif(l)*1.0E4

      if(o(l) < 1.d-50) o(l)=1.d-50
      if(o2(l) < 1.d-50) o2(l)=1.d-50
      if(nit(l) < 1.d-50) nit(l)=1.d-50
      cf1nk1 = 3.42E-17*O(l)*factor
      cf1nk2 = 6.66E-16*O2(l)
      cf1nk3 = 6.82E-16*NIT(l)
      w1 = 4.8E-13*(N1(l)+N2(l))
      w2 = SQRT(TE(l))
      w3 = w1/w2
      w4 = w1/(w2*w2*w2)
      w5 = 2.1E-21*O(l)
      w6 = (6.6E-20*NIT(l)+5.8E-20*O2(l))

  ! we now solve the heat balance equation for ti,using newton-raphson

      t1 = TE(l)
      icount = 0
      50 w = w5*SQRT(t1+TN(l)) + w6
      icount = icount + 1

      if (icount == 10 ) then
         ifail=1
         return
      endif

  ! frictional heating terms.

      tp = (t1+TN(l))/2.
      rttp = SQRT(tp)
      tplog = log10(tp)
      x = 1.04 - 0.067*tplog
      cf1a = cf1nk1*rttp*x*x
      w8 = (8.37E-12*cf1a+1.12E-11*cf1nk2+1.06E-11*cf1nk3)*u2
      t2 = (w*TN(l)+w3+w8)/(w+w4)
      IF ( ABS(t2-t1) < 0.5 ) THEN
      !c  **
      !c  modified by tjfr jan 86
      !c  **
          IF ( t2 < TN(l) ) t2 = TN(l)
          TI(l) = t2
          tp = (t2+TN(l))/2.
          rttp = SQRT(tp)
          tplog = log10(tp)
          RTTin(l) = rttp
          TINlog(l) = tplog
          CF1n(l) = cf1a + cf1nk2 + cf1nk3
      ELSE
          t1 = t2
          GOTO 50
      ENDIF
  100 ENDDO
  RETURN



end SUBROUTINE HL__ION_TEMP_HIGH_LAT


SUBROUTINE HL__PARTITION_AURORAL_IONISATION(QT,O,O2,N2,NHGt,OPRod,O2Prod,N2Prod)

  IMPLICIT NONE
  INTEGER :: n , NHGt
  real (kind=8) :: q , QT(90) , O(90) , O2(90) , N2(90) , OPRod(90) , O2Prod(90) , N2Prod(90)

! routine to break down ionization rate given in qt
! into different species

! tjfr may 85


  DO 100 n = 1 , NHGt
      q = QT(n)/(0.92*N2(n)+1.5*O2(n)+0.56*O(n))
      OPRod(n) = (0.5*O2(n)+0.56*O(n))*q
      O2Prod(n) = O2(n)*q
      N2Prod(n) = 0.92*N2(n)*q
  100 ENDDO

  RETURN



end SUBROUTINE HL__PARTITION_AURORAL_IONISATION

SUBROUTINE HL__ION_NEUTRAL_COLLISIONS(P1,P2,P3,PI1,PI2,PI3,T,VIN,AMIn,NMAx,iout)
  IMPLICIT NONE
  REAL(kind=8) :: a , AMIn , amu , b , factor , P1 , P2 , P3 , PI1 , PI2 , &
  sum , summol , T , v1 , v2 , VIN , PI3
  INTEGER :: n , NM , NMAx , iout
  PARAMETER (NM=90)
  DIMENSION P1(NM) , P2(NM) , P3(NM) , T(NM) , VIN(NM) , AMIn(NM) , &
  a(3) , b(3) , PI1(NM) , PI2(NM) , PI3(NM)
  REAL(kind=8) :: mi1 , mi2 , mi3
  DATA mi1 , mi2 , mi3/16. , 30. , 32./
  DATA a/3.42E-11 , 6.66E-10 , 6.82E-10/
  DATA b/2.44E-10 , 4.28E-10 , 4.34E-10/
  amu = 1.66E-27
!c  **
!c  **
  factor=1.0
!c  **
!c  **
  DO 100 n = 1 , NMAx
      summol = PI2(n) + PI3(n)
      sum = PI1(n) + PI2(n) + PI3(n)
      v2 = b(1)*P1(n) + b(2)*P2(n) + b(3)*P3(n)
      v1 = a(3)*P3(n) + a(2)*P2(n) + a(1)*P1(n)*factor*SQRT(T(n)) &
      *(1.08-0.139*log10(T(n))+4.51E-03*log10(T(n))**2)
      if(summol < 1.d-90) summol=0.0
      if(v1 < 1.d-90) v1=0.0
      if(v2 < 1.d-90) v2=0.0
  ! if(pi1(n).lt.1.d-90) pi1(n)=0.0
  ! if(iout.eq.1) write(6,*) 'here 5',n
  ! if(iout.eq.1) write(6,*) 'here 5.5',PI1(n) , PI2(n) , PI3(n)
  ! if(iout.eq.1) stop
      VIN(n) = (v1*PI1(n)+v2*summol)*1.E-06/sum
  ! if(iout.eq.1) write(6,*) 'here 6',n
      AMIn(n) = (PI1(n)*mi1+PI2(n)*mi2+PI3(n)*mi3)*amu/sum
  100 ENDDO
  RETURN



end SUBROUTINE HL__ION_NEUTRAL_COLLISIONS















SUBROUTINE HL__IDRIFT(VX,VY,RINcf,EX,EY,RAMin,VIX,VIY,VIZ,DIP,BTH, &
  NHGt)
  IMPLICIT NONE
  REAL(kind=8) :: &
  BTH , bx , DIP , ELCH , EX , EY , r , r2 , RAMin , RINcf , &
  sdip , VIX , VIY , VIZ , vperp , VX
  REAL(kind=8) :: VY , wi
  INTEGER :: n , NHGt
  PARAMETER (ELCH=1.602E-19)
  DIMENSION VX(90) , VY(90) , RAMin(90) , RINcf(90) , VIX(90) , &
  VIY(90) , VIZ(90) , vperp(90)
  sdip = SIN(DIP)
  bx = BTH*sdip
  DO 100 n = 1 , NHGt
      wi = ELCH*BTH/RAMin(n)
      r = RINcf(n)/wi
      r2 = r*r
      vperp(n) = (EY/BTH-r*EX/bx+r*VY(n)-r2*VX(n)*sdip)/(1.+r2)
      VIY(n) = (EX/bx+r*EY/BTH+r*VX(n)*sdip+r2*VY(n))/(1.+r2)
  100 ENDDO
  DO 200 n = 1 , NHGt
      VIX(n) = -vperp(n)*sdip
      VIZ(n) = 0.0
  200 ENDDO
  RETURN



end SUBROUTINE HL__IDRIFT




SUBROUTINE HL__calculate_altitudes (NHGT,hz_metres,dip_radians,altitude_metres)

  IMPLICIT NONE
  INTEGER :: N , NHGT
  REAL(kind=8) :: hz_metres(90)
  REAL(kind=8) :: dip_radians , sini
  REAL(kind=8) :: altitude_metres(90)

  sini = ABS(SIN(DIP_radians))
  IF ( sini < 0.18 ) sini = 0.18
  DO n = 1 , NHGt
      altitude_metres(n) = 1.0E5 + (hz_metres(n)-1.0E5) * sini
  ENDDO

  RETURN

end SUBROUTINE HL__calculate_altitudes




SUBROUTINE HL__calculate_points_along_the_tube(TOPht,NSTeps,ZKM,TARea,H,H2,H2Diff,HPRod2, &
  HPRod3)
  IMPLICIT NONE
  REAL(kind=8) :: &
  diff , H , H2 , H2Diff , hi , hi1 , HPRod2 , HPRod3 , prod , &
  re , TARea(90) , TOPht , x , z , ZKM
  INTEGER :: i , NSTeps
  DIMENSION H(90) , H2(90) , H2Diff(90) , HPRod2(90) , HPRod3(90)

! this routine sets up the heightsteps along a plasma tube(presumed
! straight).as r=1,2,....,nsteps,we have the following:
! zkm(r)=f(r)
! z(r)=(zkm(r)+re)*1000.
! h(r)=z(r+1)-z(r)
! h2(r)=h(r)*h(r)
! h2diff(r)=h2(r+1)-h2(r)
! hprod2(r)=h(r)*h(r+1)
! hprod3(r)=h(r)*h(r+1)*(h(r)+h(r+1))

! all functions of h are passed in common blocks.

  DATA re/6371./
  DIMENSION ZKM(90) , z(90)

! set lower height

  ZKM(1) = 100.
  z(1) = (ZKM(1)+re)*1.E3

  DO 100 i = 1 , 57
  !c       x=zkm(i) + 5.*exp(0.0036*(zkm(i)-100.))
      x = ZKM(i) + 5.*EXP(0.0036*(ZKM(i)-100.))
      IF ( x > TOPht ) GOTO 200
      ZKM(i+1) = x
      z(i+1) = (x+re)*1.E3
  100 ENDDO
  200 DO 300 i = 59 , 90
      x = ZKM(i-1) + 300.
      IF ( x > TOPht ) GOTO 400
      ZKM(i) = x
      z(i) = (x+re)*1.E3
  300 ENDDO
  400 NSTeps = i - 1

  DO 500 i = 1 , NSTeps - 1
      TARea(i) = (z(i+1)/z(i))**3
  500 ENDDO
  DO 600 i = 1 , NSTeps - 1
      diff = z(i+1) - z(i)
      H(i) = diff
      H2(i) = diff*diff
  600 ENDDO
  DO 700 i = 1 , NSTeps - 2
      H2Diff(i) = H2(i+1) - H2(i)
      hi = H(i)
      hi1 = H(i+1)
      prod = hi*hi1
      HPRod2(i) = prod
      HPRod3(i) = prod*(hi+hi1)
  700 ENDDO
! do i = 1 , nsteps
! write(6,*) i, zkm(i)
! enddo
  RETURN



end SUBROUTINE HL__calculate_points_along_the_tube







SUBROUTINE HL__NEARHT(ARR,VALue,M,N,NEArest)
  IMPLICIT NONE
  REAL(kind=8) :: ARR , diff1 , diff2 , VALue
  INTEGER :: i , M , N , NEArest

! this routine finds which element of arr,assumed monotonic increasi
! is nearest to value.


  DIMENSION ARR(90)
  DO 100 i = M , N
      IF ( ARR(i) > VALue ) GOTO 200
  100 ENDDO
  200 diff1 = ARR(i) - VALue
  diff2 = VALue - ARR(i-1)
  NEArest = i
  IF ( diff2 < diff1 ) NEArest = i - 1
  RETURN



end SUBROUTINE HL__NEARHT














SUBROUTINE HL__ITRPD1D2(V13d,V23d,D13d,D23d,NHGt)
  IMPLICIT NONE
  REAL(kind=8) :: D13d , d1n , d1s , D23d , d2n , d2s , V13d , v1n , v1s , &
  V23d , v2n , v2s
  INTEGER :: l , m , n , NHGt

!- define parameters

  DIMENSION D13d(90,91,20,2) , D23d(90,91,20,2)

  DIMENSION V13d(90,91,20,2) , V23d(90,91,20,2)

  DO 200 n = 1 , NHGt
      d1n = 0.0
      d1s = 0.0
      d2n = 0.0
      d2s = 0.0
      v1n = 0.0
      v1s = 0.0
      v2n = 0.0
      v2s = 0.0
      DO 50 l = 1 , 20
          DO 20 m = 2 , 90
              V13d(n,m,l,1) = V13d(n,m,l,2)
              V23d(n,m,l,1) = V23d(n,m,l,2)
              D13d(n,m,l,1) = D13d(n,m,l,2)
              D23d(n,m,l,1) = D23d(n,m,l,2)
          20 ENDDO
          d1n = d1n + D13d(n,89,l,1)
          d2n = d2n + D23d(n,89,l,1)
          d1s = d1s + D13d(n,3,l,1)
          d2s = d2s + D23d(n,3,l,1)
          v1n = v1n + V13d(n,89,l,1)
          v1s = v1s + V13d(n,3,l,1)
          v2n = v2n + V23d(n,89,l,1)
          v2s = v2s + V23d(n,3,l,1)

      !c  set boundary conditions for cut off at 40 lat
      !c       d13d(n,26,l,1)=d13d(n,25,l,1)
      !c       d23d(n,26,l,1)=d23d(n,25,l,1)
      !c       d13d(n,66,l,1)=d13d(n,67,l,1)
      !c       d23d(n,66,l,1)=d23d(n,67,l,1)

      !c       v13d(n,26,l,1)=v13d(n,25,l,1)
      !c       v13d(n,66,l,1)=v13d(n,67,l,1)
      !c       v23d(n,26,l,1)=v23d(n,25,l,1)
      !c       v23d(n,66,l,1)=v23d(n,67,l,1)
      !c
      !c  set boundary conditions for cut off at 20 lat
      !c       d13d(n,36,l,1)=d13d(n,35,l,1)
      !c       d23d(n,36,l,1)=d23d(n,35,l,1)
      !c       d13d(n,56,l,1)=d13d(n,57,l,1)
      !c       d23d(n,56,l,1)=d23d(n,57,l,1)

      !c       v13d(n,36,l,1)=v13d(n,35,l,1)
      !c       v13d(n,56,l,1)=v13d(n,57,l,1)
      !c       v23d(n,36,l,1)=v23d(n,35,l,1)
      !c       v23d(n,56,l,1)=v23d(n,57,l,1)
      50 ENDDO
      d1n = d1n/20.
      d2n = d2n/20.
      d1s = d1s/20.
      d2s = d2s/20.
      v1n = v1n/20.
      v2n = v2n/20.
      v1s = v1s/20.
      v2s = v2s/20.

      DO 100 l = 1 , 20
          D13d(n,91,l,1) = d1n
          D23d(n,91,l,1) = d2n
          D13d(n,1,l,1) = d1s
          D23d(n,1,l,1) = d2s
          V13d(n,91,l,1) = v1n
          V13d(n,1,l,1) = v1s
          V23d(n,91,l,1) = v1n
          V23d(n,1,l,1) = v1s
          D13d(n,90,l,1) = (d1n+D13d(n,89,l,1))/2.0
          D23d(n,90,l,1) = (d2n+D23d(n,89,l,1))/2.0
          D13d(n,2,l,1) = (d1s+D13d(n,3,l,1))/2.0
          D23d(n,2,l,1) = (d2s+D23d(n,3,l,1))/2.0
          V13d(n,90,l,1) = (v1n+V13d(n,89,l,1))/2.
          V13d(n,2,l,1) = (v1s+V13d(n,3,l,1))/2.
          V23d(n,90,l,1) = (v2n+V23d(n,89,l,1))/2.
          V23d(n,2,l,1) = (v2s+V23d(n,3,l,1))/2.
      100 ENDDO
  200 ENDDO

  RETURN



end SUBROUTINE HL__ITRPD1D2






SUBROUTINE HL_ML__HYDEQ_BASIC(H_400,PZ,TN,G,hydrogen_out)

  IMPLICIT NONE

  REAL(kind=8), intent(in) :: PZ , TN , G , H_400
  REAL(kind=8), intent(out) :: Hydrogen_out

  REAL(kind=8) :: SHT
  REAL(kind=8) :: GASCON
  REAL(kind=8) :: dh

  PARAMETER (GASCON=8.3141E+03)

  dh = 400.E+03 - PZ*1.0E03
  sht = GASCON*TN/G
! Hydrogen_out = 5.e10 * EXP(dh/sht)
  Hydrogen_out = H_400 * EXP(dh/sht)

  RETURN


end SUBROUTINE HL_ML__HYDEQ_BASIC






SUBROUTINE HL_ML__EUV_ION_PRODUCTION(n_array_size,i1,i2,o,o2,n2,Solar_Declination_Angle_degrees,time,f107, &
  peuvi,gcol,r0,tn,gravity,km,gr,CHID,iout)

!***********************************************************************
! ROUTINE TO CALCULATE THE PHOTOIONIZATION PRODUCTION RATE

! SOLAR FLUXES FROM TOBISKA MODEL
! (JATP XXVIII COSPAR EDITION 1991)

! CROSS SECTIONS FROM TORR AND TORR
! (TORR AND TORR, 1982, REV. GEOPHYS. SPACE RES., 20, 91-144)

! NIGHTTIME PRODUCTION DUE TO RING CURRENT
! (PROLSS, IAGA93)
!***********************************************************************

  implicit none
  integer :: mm,k,j,i,i1,i2,iout,n_array_size
  REAL(kind=8) :: CSIO,CSIO2,CSIN2,CSIHE,PO,PO2,PN2,PHE, &
  sf1,sf2,r0,TN,gravity,km(6),gr
  REAL(kind=8) :: pye,dtr,Solar_Declination_Angle_degrees,f107,cschi,chid,csfye,emtau, &
  emtau_tot
  REAL(kind=8) :: O,O2,N2,peuvi,time,gcol,sf,chi,w,f107a1,f107a2
  REAL(kind=8) :: pngt , pngto2 , pngtn2
  PARAMETER(PYE=3.14159)
  PARAMETER(DTR=PYE/180.)

  PARAMETER (PNGT=1.0E-17, &
  PNGTO2=2.0E-17,PNGTN2=1.8E-17)

! PARAMETER (PNGT=1.0E-23, &
!   PNGTO2=2.0E-23,PNGTN2=1.8E-23)

! PARAMETER (PNGT=1.0E-26, &
!   PNGTO2=2.0E-26,PNGTN2=1.8E-26)

! PARAMETER (PNGT=0.0, &
!   PNGTO2=0.0,PNGTN2=0.0)

  SAVE CSIO,CSIO2,CSIN2,CSIHE,MM,PO,PO2,PN2,PHE,SF

  DIMENSION CSIO(37),CSIO2(37),CSIN2(37),CSIHE(37) &
  ,PO(37),PO2(37),PN2(37),PHE(37) &
  ,EMTAU(37),CHID(n_array_size),CSCHI(n_array_size)
  dimension O(n_array_size),O2(n_array_size),N2(n_array_size),TN(n_array_size),gravity(n_array_size) &
  ,PEUVI(n_array_size,6),gr(n_array_size) &
  ,TIME(n_array_size),GCOL(n_array_size),SF(37),CHI(n_array_size),sf1(37),sf2(37)

!--solar flux spectrum sc 21refw for solar minimum
!--in 10**13 photons m-2 sec-1
  data f107a1/68./,sf1/ &
  0.3834, 0.1346, 1.8418, 0.9235, 0.2713, 0.1000, 0.8405, 0.2350, &
  6.0000, 0.8661, 0.7394, 0.2121, 0.3926, 0.1800, 0.3063, 0.5085, &
  0.7992, 1.5800, 0.4843, 0.4500, 1.5000, 0.1746, 0.2223, 0.3915, &
  0.1667, 0.1997, 0.2425, 0.7931, 0.8728, 1.9311, 4.4325, 4.2170, &
  5.9570, 1.7850, 4.3750, 3.1840, 3.6401/

!--solar flux spectrum f79050n for solar maximum
!--in 10**13 photons m-2 sec-1
  data f107a2/243./,sf2/ &
  1.1487, 0.3433, 4.8498, 3.7013, 0.5947, 3.1675, 4.1358, 2.4995, &
  11.2800, 5.6326, 1.3949, 2.1965, 0.9932, 0.3621, 1.6716, 1.5468, &
  1.5904, 4.8664, 1.0213, 1.4621, 3.0180, 0.4820, 0.4554, 0.7165, &
  0.4256, 0.4318, 0.6709, 1.5869, 2.1809, 5.0135,13.2975,12.0345, &
  13.1770, 4.4204,13.1251, 9.0426, 8.6669/

!***********************************************************************
! ALL CROSS SECTIONS ARE GIVEN IN UNITS 10**-22 M2
!***********************************************************************

!--PHOTOIONIZATION CROSS SECTIONS - O
  DATA CSIO/ &
  1.06, 3.53, 5.96, 7.55, 8.43, 9.26, 8.78, 9.70, &
  9.72,10.03,10.84,10.70,11.21,11.25,11.64,11.91, &
  12.13,12.17,11.90,12.23,12.22,12.21,10.04,11.35, &
  8.00, 4.18, 4.18, 4.28, 4.23, 4.38, 4.18, 2.12, &
  0.00, 0.00, 0.00, 0.00, 0.00/

!--PHOTOIONIZATION CROSS SECTIONS - O2
  DATA CSIO2/ &
  1.18, 3.61, 7.27,10.50,12.80,14.80,13.65,15.98, &
  16.00,17.19,18.40,18.17,19.39,20.40,21.59,24.06, &
  25.59,22.00,25.04,26.10,25.80,25.94,22.05,23.00, &
  23.81, 8.59, 9.69,11.05, 9.39, 6.12, 4.69, 9.34, &
  2.50,12.22, 1.00, 0.00, 0.27/

!--PHOTOIONIZATION CROSS SECTIONS - N2
  DATA CSIN2/ &
  0.60, 2.32, 5.40, 8.15, 9.65,10.60,10.08,11.58, &
  11.60,14.60,18.00,17.51,21.07,21.80,21.85,24.53, &
  24.69,23.20,22.38,23.10,23.20,23.22,25.06,23.00, &
  23.20,23.77,18.39,10.18,16.75, 0.00, 0.00, 0.00, &
  0.00, 0.00, 0.00, 0.00, 0.00/

!--PHOTOIONIZATION CROSS SECTIONS - HE
  DATA CSIHE/ &
  0.21, 0.53, 1.02, 1.71, 2.16, 2.67, 2.38, 3.05, &
  3.05, 3.65, 4.35, 4.25, 5.51, 6.53, 7.09, 0.72, &
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
  0.00, 0.00, 0.00, 0.00, 0.00/

  DATA MM/0/

!--FIRST TIME ROUTINE IS ENTERED CALCULATE THE SOLAR EUV FLUX
!--FOR EACH WAVELENGTH BAND AND SET THE IONIZATION & ABSORPTION CROSS
!--SECTIONS TO THE CORRECT MAGNITUDE

  IF(MM == 0) THEN
      w=(f107-f107a1)/(f107a2-f107a1)
      DO 1 K=1,37
          sf(k)=(sf1(k)+w*(sf2(k)-sf1(k)))*1.e13
          PO(K)=SF(K)*CSIO(K)*1.E-22
          PO2(K)=SF(K)*CSIO2(K)*1.E-22
          PN2(K)=SF(K)*CSIN2(K)*1.E-22
          PHE(K)=SF(K)*CSIHE(K)*1.E-22
      1 ENDDO
      MM=1
  ENDIF

!--SET THE PRODUCTION RATE, PEUVI, EQUAL TO ZERO AND CALCULATE THE SOLAR
!--ZENITH ANGLE CHID AND THE SCALE HEIGHT FOR HYDROGEN, H.

  DO J=1,6
     DO I = i1 , i2
        PEUVI(I,J)=0.
     ENDDO
  ENDDO


  DO 101 I = i1 , i2

  ! FOR DAYTIME PRODUCTION

      if(o(i) < 1.d-50) o(i)=1.d-50
      if(o2(i) < 1.d-50) o2(i)=1.d-50
      if(n2(i) < 1.d-50) n2(i)=1.d-50
      CSFYE=COS(PYE*(TIME(I)/43200.-1.))
      CSCHI(I)=Sin(Solar_declination_Angle_degrees * DTR)*COS(GCOL(I)) &
      +Cos(Solar_declination_Angle_degrees * DTR)*SIN(GCOL(I))*CSFYE
      CHI(I)=ACOS(CSCHI(I))
      CHID(I)=CHI(I)/dtr
      if(iout == 1) then
      ! write(6,*) i,chi(i),gravity(i),tn(i),gr(i),o(i),o2(i),n2(i)
      ! write(6,*) i,chi(i),gravity(i),tn(i),gr(i),o(i),o2(i),n2(i)
      endif

      CALL HL_ML__OPTICAL_DEPTH(CHI(i),r0,tn(i),gravity(i),km,gr(i),o(i),o2(i),n2(i),EMTAU)

      if(iout == 1) then
      ! write(6,*) 'finished '
      endif
      emtau_tot = 0.
      do k=1,37
      ! emtau(k)=0.0
          emtau_tot = emtau_tot + emtau(k)
      enddo

  !--O+ PHOTOIONIZATION RATE ( WAVELENGTHS .LT. 910 A )

      DO 10 K=1,32
          PEUVI(I,1)=PEUVI(I,1)+PO(K)*EMTAU(K)
      10 ENDDO
      PEUVI(I,1)=(PEUVI(I,1)+pngt)*O(I)

  !--HE+ PHOTOIONIZATION RATE ( WAVELENGTHS .LT. 504 A )

  ! IF(JION(3).NE.0) THEN
  ! DO 30 K=1,16
  ! PEUVI(I,3)=PEUVI(I,3)+PHE(K)*EMTAU(K)
  ! PEUVN(I,3)=PEUVN(I,3)+PHE(K)*EMTAUN(K)
  ! 30    CONTINUE
  ! PEUVI(I,3)=PEUVI(I,3)*HEL(I)
  ! PEUVN(I,3)=PEUVN(I,3)*HEL(I)/PFAC3
  ! IF(PEUVI(I,3).GT.PEUVN(I,3)) THEN
  ! PEUVN(I,3)=PEUVI(I,3)
  ! ENDIF
  ! ENDIF

  !--N2+ PHOTOIONIZATION RATE ( WAVELENGTHS .LT. 800 A )

  ! IF(JION(4).NE.0) THEN
      DO 40 K=1,29
          PEUVI(I,4)=PEUVI(I,4)+PN2(K)*EMTAU(K)
      40 ENDDO
      PEUVI(I,4)=(PEUVI(I,4)+pngtn2)*N2(I)
  ! ENDIF

  !--O2+ PHOTOIONIZATION RATE ( WAVELENGTHS .LT. 1050 A )

  ! IF(JION(5).NE.0) THEN
      DO 50 K=1,37
          PEUVI(I,5)=PEUVI(I,5)+PO2(K)*EMTAU(K)
      50 ENDDO
      PEUVI(I,5)=(PEUVI(I,5)+pngto2)*O2(I)
  ! ENDIF
  ! if(iout.eq.1) then
  ! write(6,2424) i,peuvi(i,1),o(i)
  ! endif
  ! 424       format(' ion_prod',i4,3e12.4)

  101 ENDDO

  RETURN




end SUBROUTINE HL_ML__EUV_ION_PRODUCTION




SUBROUTINE HL_ML__MOLECULAR_IONS_ON_TUBES( &
   GIP_switches, &
   n_array_size,i1,i2,f107,O_number_density_m3,O2_number_density_m3,N2_number_density_m3, &
   NO_number_density_m3,N4S_number_density_m3,N2D_number_density_m3, &
   tn,te,ti,OPLUS_number_density_m3,ht,cos_SZA, &
   euv_production_rate_Oplus,euv_production_rate_N2plus,euv_production_rate_O2plus, &
   auroral_production_rate_Oplus,auroral_production_rate_N2plus,auroral_production_rate_O2plus, &
   N2PLUS_number_density_m3,NOPLUS_number_density_m3,O2PLUS_number_density_m3,NPLUS_number_density_m3, &
   ifailed)


  IMPLICIT NONE

! Inputs and Outputs....................................................

  integer, intent(in) :: n_array_size
  logical :: GIP_switches(20)
  logical :: sw_External_model_provides_NO_N4S_densities

  real (kind=8), intent(in) :: TE(n_array_size)
  real (kind=8), intent(in) :: TI(n_array_size)
  real (kind=8), intent(in) :: TN(n_array_size)
  real (kind=8), intent(in) :: HT(n_array_size)
  real (kind=8), intent(in) :: O_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: O2_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: N2_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: NO_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: N4S_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: N2D_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: OPLUS_number_density_m3(n_array_size)
  real (kind=8), intent(in) :: cos_SZA(n_array_size)
  real (kind=8), intent(in) :: F107
  real (kind=8), intent(in) :: euv_production_rate_Oplus(n_array_size)
  real (kind=8), intent(in) :: euv_production_rate_N2plus(n_array_size)
  real (kind=8), intent(in) :: euv_production_rate_O2plus(n_array_size)
  real (kind=8), intent(in) :: auroral_production_rate_Oplus(n_array_size)
  real (kind=8), intent(in) :: auroral_production_rate_N2plus(n_array_size)
  real (kind=8), intent(in) :: auroral_production_rate_O2plus(n_array_size)
  integer, intent(in) :: i1 , i2

  real (kind=8), intent(out) :: NOPLUS_number_density_m3(n_array_size)
  real (kind=8), intent(out) :: O2PLUS_number_density_m3(n_array_size)
  real (kind=8), intent(out) :: N2PLUS_number_density_m3(n_array_size)
  real (kind=8), intent(out) :: NPLUS_number_density_m3(n_array_size)
  integer, intent(out) :: ifailed

! ......................................................................

  real (kind=8) :: O_density_cm3(n_array_size)
  real (kind=8) :: O2_density_cm3(n_array_size)
  real (kind=8) :: N2_density_cm3(n_array_size)
  real (kind=8) :: OPLUS_density_cm3(n_array_size)
  real (kind=8) :: ELECTRON_density_cm3(n_array_size)
  real (kind=8) :: NOPLUS_density_cm3(n_array_size)
  real (kind=8) :: O2PLUS_density_cm3(n_array_size)
  real (kind=8) :: N2PLUS_density_cm3(n_array_size)
  real (kind=8) :: NPLUS_density_cm3(n_array_size)
  real (kind=8) :: sht(n_array_size)
  real (kind=8) :: term
  real (kind=8) :: uvo(n_array_size)
  real (kind=8) :: uvo2(n_array_size)
  real (kind=8) :: uvn2(n_array_size)
  real (kind=8) :: aurqo(n_array_size)
  real (kind=8) :: aurqo2(n_array_size)
  real (kind=8) :: aurqn2(n_array_size)

  real (kind=8) :: rno1(n_array_size)
  real (kind=8) :: rn4sl(n_array_size)
  real (kind=8) :: rn4sh(n_array_size)
  real (kind=8) :: rno1l(n_array_size)
  real (kind=8) :: rno1h(n_array_size)
  real (kind=8) :: ac(5)
  real (kind=8) :: z(2,4)
  real (kind=8) :: w(20)
  real (kind=8) :: RN4SL_15(15)
  real (kind=8) :: RNO1L_15(15)
  real (kind=8) :: RN4SH_15(15)
  real (kind=8) :: RNO1H_15(15)
  real (kind=8) :: EUVOP_15(15)
  real (kind=8) :: EUVNP_15(15)
  real (kind=8) :: PEOP_15(15)
  real (kind=8) :: PEO2_15(15)
  real (kind=8) :: PEN2_15(15)
  real (kind=8) :: rn4s(n_array_size)
  real (kind=8) :: euvop(n_array_size)
  real (kind=8) :: euvnp(n_array_size)
  real (kind=8) :: peop(n_array_size)
  real (kind=8) :: pen2(n_array_size)
  real (kind=8) :: peo2(n_array_size)

  real (kind=8) :: a , a0 , a1 , a2 , a3 , a4 , af1 , af2 , af3 , b , c , cg , ch
  real (kind=8) :: d , e , ee , f , g , p 
  real (kind=8) :: k5 , k6 , k7 , k10 , k15 , ly , lyl , lyh
  real (kind=8) :: pq , q , ri , rj , rl10 , rl11 , rl2 , rl4 , rl5 , rl7 , rl8 , rl9 , rlam
  real (kind=8) :: RMt 
  real (kind=8) :: root , rp10 , rp11 , rp14 , rp17 , rp18 , rp19 , rp20 , rp5 , rp8 , rp9
  INTEGER :: ifail , n , nc , l
  LOGICAL :: scale
  COMPLEX*16 :: cp


  sw_External_model_provides_NO_N4S_densities = GIP_switches(5)
  ifailed = 0
  scale = .TRUE.

  ee = 0.0E-40

! Chemical reaction rates in units of cm-3 s-1...

!C        LY=3.E-7
  lyl = 2.6E11
  lyh = 4.0E11
  ly = lyl + (F107-67.)/(243.-67.)*(lyh-lyl)
  k5 = 2.0E-10
  k6 = 1.2E-10
  k7 = 4.4E-10
  k10 = 4.0E-10
  k15 = 1.0E-12

!g
!g the following data is on 15 pressure levels.
!g Need to interpolate this onto our flux tube (from 1 to 90)...
!g

  if (.NOT. sw_External_model_provides_NO_N4S_densities) then 
!
! if we are not getting the NO and N4S from the external thermospheric model
! then we need to calculate these as simple 1D profiles.....
!
  DATA rn4sl_15/10.13 , 10.56 , 11.07 , 11.59 , 12.04 , 12.44 , 12.71 , &
                12.84 , 13.05 , 13.14 , 13.02 , 12.77 , 12.45 , 12.10 , 11.74/

  DATA rno1l_15/11.82 , 11.81 , 11.90 , 12.02 , 12.17 , 12.31 , 12.37 , &
                12.41 , 12.36 , 12.22 , 11.86 , 11.23 , 10.51 , 9.78 , 9.04/

  DATA rn4sh_15/10.56 , 10.72 , 10.83 , 11.12 , 11.55 , 11.94 , 12.24 , &
                12.51 , 12.94 , 13.21 , 13.17 , 12.91 , 12.54 , 12.18 , 11.81/

  DATA rno1h_15/12.64 , 12.78 , 12.98 , 13.16 , 13.23 , 13.26 , 13.14 , &
                13.01 , 12.78 , 12.45 , 11.89 , 11.21 , 10.51 , 9.76 , 9.02/

  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,rn4sl_15,rn4sl)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,rno1l_15,rno1l)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,rn4sh_15,rn4sh)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,rno1h_15,rno1h)

  do l = i1 , i2
      rn4s(l) = 10**(rn4sl(l)+(F107-67.)/(243.-67.)*(rn4sh(l)-rn4sl(l))-6.0)
      rno1(l) = 10**(rno1l(l)+(F107-67.)/(243.-67.)*(rno1h(l)-rno1l(l))-6.0)
  enddo

  else
!
! ...else we pass these numbers in from outside. NO and N4S are passed in as m-3 
! but are needed here in cm-3, hence the 1.0e-06 factor......
!
  do l = i1 , i2
      rn4s(l) = N4S_number_density_m3(l)*1.0E-6
      rno1(l) = NO_number_density_m3(l)*1.0E-6
  enddo

  endif


  DATA euvop_15/4*0.0 , 0.002 , 0.0055 , 0.0185 , 0.08 , 0.1475 , &
  0.1625 , 0.155 , 4*0.15/
  DATA euvnp_15/3*0.4 , 0.405 , 0.405 , 0.395 , 0.385 , 0.315 , 0.205 , &
  0.27 , 0.1 , 4*0.09/
  DATA peop_15/4*5.67 , 4.86 , 3.405 , 1.99 , 0.835 , 0.35 , 0.21 , &
  0.1625 , 0.1525 , 0.145 , 2*0.14/
  DATA peo2_15/4*.01 , 0.015 , 0.04 , 0.105 , 0.18 , 0.185 , 0.14 , &
  0.115 , 0.1085 , 0.1035 , 2*0.1/
  DATA pen2_15/81. , 52. , 35. , 20. , 11. , 6.9 , 4.2 , 1.7 , .47 , &
  .23 , .16 , .15 , .14 , .13 , .12/


  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,euvop_15,euvop)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,euvnp_15,euvnp)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,peop_15,peop)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,peo2_15,peo2)
  CALL HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,ht,pen2_15,pen2)

  do l = i1 , i2

      O_density_cm3(l) = O_number_density_m3(l)*1.0E-6
      O2_density_cm3(l) = O2_number_density_m3(l)*1.0E-6
      N2_density_cm3(l) = N2_number_density_m3(l)*1.0E-6
      OPLUS_density_cm3(l) = OPLUS_number_density_m3(l)*1.0E-6

      rmt = ((O_density_cm3(l)*16.)+(N2_density_cm3(l)*28.)+(O2_density_cm3(l)*32.))/ &
            (O_density_cm3(l)+N2_density_cm3(l)+O2_density_cm3(l))
      sht(l) = 8314.*TN(l)*100./RMt/9.5

      uvo(l)=euv_production_rate_Oplus(l)*1.0E-6
      uvo2(l)=euv_production_rate_O2plus(l)*1.0E-6
      uvn2(l)=euv_production_rate_N2plus(l)*1.0E-6
      aurqo(l)=auroral_production_rate_Oplus(l)*1.0E-6
      aurqo2(l)=auroral_production_rate_O2plus(l)*1.0E-6
      aurqn2(l)=auroral_production_rate_N2plus(l)*1.0E-6

  enddo



  DO 500 l = i1 , i2
  !g
  !g only do this calculation for altitudes of less than 1000 km
  !g (above this moleculars get stupidly small and sometimes NaN)..
  !g
      if (ht(l) <= 1000.) then

          rl7 = O2_density_cm3(l)*k10
          rl8 = O2_density_cm3(l)*k5
          rl9 = O_density_cm3(l)*k15
          rp17 = aurqn2(l)*0.24
          NPLUS_density_cm3(l) = (rp17+uvn2(l)*(euvnp(l)+0.24*pen2(l)))/(rl7+rl8+rl9)
          g = NPLUS_density_cm3(l)

          rl10 = O2_density_cm3(l)*K8(TI(l),TN(l))
          rl11 = N2_density_cm3(l)*K3(TI(l),TN(l))
          rp14 = rl9*NPLUS_density_cm3(l)
          rp18 = aurqo2(l)*0.33 + aurqo(l)

          IF ( OPLUS_density_cm3(l) <= 0. ) then
             OPLUS_density_cm3(l) = (uvo2(l)*euvop(l)+uvo(l)*(1.+peop(l))+0.33*peo2(l)*uvo2(l)+rp14+rp18)/(rl10+rl11)
          ENDIF

          IF ( OPLUS_density_cm3(l) <= 0.0 ) OPLUS_density_cm3(l) = 1.E-25

          f = OPLUS_density_cm3(l)
          rl2 = O_density_cm3(l)*K4(TI(l),TN(l))
          e = rl2
          rp19 = aurqn2(l)*0.76
          d = rp19 + uvn2(l)*((1.-euvnp(l))+0.76*pen2(l))
          rp5 = rno1(l)*ly*1.E-3*2.E-18*EXP(-1.E-20*O2_density_cm3(l)*sht(l))
          IF ( cos_SZA(l) > 0.001 ) rp5 = rno1(l) &
          *ly*2.E-18*EXP(-1.E-20*O2_density_cm3(l)*sht(l) &
          /cos_SZA(l))
          rp9 = OPLUS_density_cm3(l)*N2_density_cm3(l)*K3(TI(l),TN(l))
          rp10 = NPLUS_density_cm3(l)*rl8
          a = rp5 + rp9 + rp10
          rl4 = rn4s(l)*k6
          rl5 = rno1(l)*k7
          c = rl4 + rl5
          rp8 = OPLUS_density_cm3(l)*O2_density_cm3(l)*K8(TI(l),TN(l))
          rp11 = NPLUS_density_cm3(l)*rl7
          rp20 = aurqo2(l)*0.67
          b = uvo2(l)*((1.-euvop(l))+0.67*peo2(l)) + rp8 + rp11 + rp20
          af1 = 4.2E-7*(300.0/TE(l))**.85
          IF ( TE(l) >= 1200.0 ) af2 = 1.6E-7*(300.0/TE(l))**.55
          IF ( TE(l) < 1200.0 ) af2 = 2.7E-7*(300.0/TE(l))**.7
          af3 = 1.8E-7/((TE(l)/300.0)**0.39)
          a4 = -e*c*(a+b+d)
          a3 = (-af1*(e*c*(f+g)+d*c+b*e)-af2*e*(a+d)-af3*c*(a+b))/4.
          a2 = (af1*e*c-af1*(af2*e+af3*c)*(f+g)-af1*af2*d-af1*af3*b- &
          af2*a*af3)/6.
          a1 = (af1*(af2*e+af3*c)-af1*af2*af3*(f+g))/4.
          a0 = af1*af2*af3
          ac(1) = a0
          ac(2) = a1*4.
          ac(3) = a2*6.
          ac(4) = a3*4.
          ac(5) = a4
          ri = a0*a4 - 4.*a1*a3 + 3.*a2**2
          rj = a0*(a2*a4-a3**2) - a1*(a1*a4-a3*a2) + a2*(a1*a3-a2**2)
          ch = -ri/12.
          cg = rj/4.
          cp = cg**2 + 4.*ch**3 + ee
          cp = (.5*(cg+cp**.5)+ee)**(1./3.)
          rlam = -2.*REAL(cp)
          p = a0*rlam + a1**2 - a0*a2 + ee
          IF ( p <= 0.0 ) THEN
              GOTO 250
          ENDIF
          p = SQRT(p)
          q = (2.*rlam+a2)**2 - a0*a4 + ee
          IF ( q <= 0.0 ) THEN
          ! WRITE(6,99003)
          ! WRITE(6,99002) q , rlam , a2 , a0 , a4 , ee
          ! WRITE(6,99005) l , O_density_cm3(l) , O2_density_cm3(l) , N2_density_cm3(l) , OPLUS_density_cm3(l) ,
          ! &                      NPLUS_density_cm3(l)
          ! WRITE(6,99005) l , a , b , c , d , e
          ! WRITE(6,99004) rp19 , uvn2(l) , aurqn2(l) , 
              GOTO 250
          ENDIF
          pq = 2.*a1*rlam + a1*a2 - a0*a3 + ee
          q = SQRT(q)
          p = SIGN(p,q*pq)
          root = p - a1
          term = root**2 - a0*(a2+2.*rlam-q)
      ! IF ( term.LE.0.0 ) THEN
      ! WRITE(6,*) term , p , a1
      ! WRITE(6,99002) q , rlam , a2 , a0 , a4 , ee
      ! WRITE(6,99005) l , O_density_cm3(l) , O2_density_cm3(l) , N2_density_cm3(l) , OPLUS_density_cm3(l) ,
      ! &                      NPLUS_density_cm3(l)
      ! WRITE(6,99005) l , a , b , c , d , e
      ! WRITE(6,99004) rp19 , uvn2(l) , aurqn2(l) ,
      ! ENDIF
          ELECTRON_density_cm3(l) = (root+SQRT(term))/a0
      !c      ELECTRON_density_cm3(l)=(ROOT+SQRT(ROOT**2-A0*(A2+2.*RLAM-Q)))/A0
      ! IF ( ELECTRON_density_cm3(l).LE.0.0 ) goto 250
          IF ( ELECTRON_density_cm3(l) <= 0.0 ) then
              goto 250
          ENDIF
          GOTO 400
          250 nc = 4
          ifail = 0
          write(6,*) ' Here in MOLECULAR IONS ON TUBES - need CO2AGF'
!
!     if we get to here Molecular ions is going to fail.
!     just set a fail flag and carry on....
!
          ifailed = 1
          if(ifailed == 1) then
          goto 1500
          endif
!
!           CALL C02AGF(ac,nc,scale,z,w,ifail)
      !c        call c02age(ac,nc,scale,z,w,ifail)
          DO 300 n = 1 , 4
              IF ( ABS(z(2,n)) <= 0.0 ) THEN
                  IF ( z(1,n) > 0.0 ) THEN
                      ELECTRON_density_cm3(l) = z(1,n)
                      GOTO 350
                  ENDIF
              ENDIF
          300 ENDDO
          350 IF ( n == 5 ) WRITE(6,99006)
          IF ( ifail /= 0 ) WRITE(6,99007) ifail

400        CONTINUE

          N2PLUS_density_cm3(l) = d/(e+af3*ELECTRON_density_cm3(l))
          O2PLUS_density_cm3(l) = b/(c+af2*ELECTRON_density_cm3(l))
          NOPLUS_density_cm3(l) = (a+e*N2PLUS_density_cm3(l)+c*O2PLUS_density_cm3(l))/(af1*ELECTRON_density_cm3(l))

      else

      !g  .. if above 1000km just set our molecular densities to zero...
      !g  .... actually set them to a small number - but not zero
      !g  (can lead to crashes in fastinterp).....

          N2PLUS_density_cm3(l) = 1.e-6
          O2PLUS_density_cm3(l) = 1.e-6
          NOPLUS_density_cm3(l) = 1.e-6
          NPLUS_density_cm3(l) = 1.e-6

      !g end of the if below 1000km if....

      endif

  500 ENDDO

      N2PLUS_number_density_m3(i1:i2) = N2PLUS_density_cm3(i1:i2)*1.0E6
      O2PLUS_number_density_m3(i1:i2) = O2PLUS_density_cm3(i1:i2)*1.0E6
      NOPLUS_number_density_m3(i1:i2) = NOPLUS_density_cm3(i1:i2)*1.0E6
      NPLUS_number_density_m3(i1:i2) = NPLUS_density_cm3(i1:i2)*1.0E6

1500 continue

  RETURN
  99001 FORMAT (1x,'p is negative')
  99002 FORMAT (1x,1P6e12.4)
  99003 FORMAT (1x,'q is negative')
  99004 FORMAT (1x,1P5e12.4)
  99005 FORMAT (1x,i3,1P5e12.4)
  99006 FORMAT (1x,'no root found in molion')
  99007 FORMAT (1x,'ifail = ',i2)
  99008 FORMAT (1x,1P15e8.1)
  99009 FORMAT (1x,1P15e8.1,'ELECTRON_density_cm3 neg')




end SUBROUTINE HL_ML__MOLECULAR_IONS_ON_TUBES





  ! ------------------------------------------------------------------------
  !                          HL_ML__EUV_ION_PRODUCTION_2
  ! ------------------------------------------------------------------------
  ! Alison Dobbin: March 2009
  ! ali.dobbin@gmail.com
  ! Routine to calculate photo-ionisation rates 
  ! Uses solar fluxes from EUVAC, (TIEGCM). 
  ! Ref spectrum: Solomon and Qian, JGR, V110, 2005. 
  ! Fluxes and cross sections supplied by Qian.
  !
  ! Passed in:
  ! npts                 : vertical array dimension
  ! low_level            : lowest pressure or height level counter
  ! top_level            : top pressure or height level counter (15 or 90)
  ! O_ndensity_1d(npts)  : Number density O (m-3)
  ! O2_ndensity_1d(npts) : Number density O2 (m-3)
  ! N2_ndensity_1d(npts) : Number density N2 (m-3)
  ! solar_declination_angle_degrees : solar declination angle
  ! time (npts)          : current time (seconds)
  ! f107                 : F107a
  ! colatitude_geo(npts) : Geographic colatitude radians
  ! R0(npts)             : radius of earth (m)
  ! temperature_1d(npts) : Neutral temp (kelvin) 
  ! gravity(npts)        : acc due to gravity (m/s^2)
  ! km_plasma(6)         : =Boltzmann's const/(atomic_weight_X*1.673E-27)
  ! R0_eff(npts)         : radius of earth + altitude (m)
  ! eccentric            : eccentricity of earth's orbit
  !
  ! Passed out:
  ! chi_degrees(npts)    : solar zenith angle (degrees)
  ! ionisation_rates(npt,6) : Ionisation rates of major
  !                     species (O = 1, O2 = 5, N2 = 4) due to solar EUV 
  !                     (m-3s-1)
  !

    
SUBROUTINE HL_ML__EUV_ION_PRODUCTION_2(npts, low_level, top_level, O_ndensity_1d, &
                                O2_ndensity_1d, N2_ndensity_1d,                &
                                solar_declination_angle_degrees, time,         &
                                f107, colatitude_geo, R0, temperature_1d,      &
                                gravity, km_plasma, R0_eff,                    &
                                eccentric, chi_degrees, ionisation_rates)
    
    implicit none
    
    INTEGER, INTENT(IN)       :: npts ! height array dimension
    INTEGER, INTENT(IN)       :: low_level ! lowest pressure or height level counter
    INTEGER, INTENT(IN)       :: top_level ! top pressure or height counter 
    REAL(kind=8), INTENT(IN)  :: O_ndensity_1d(npts)  !m-3
    REAL(kind=8), INTENT(IN)  :: O2_ndensity_1d(npts) !m-3
    REAL(kind=8), INTENT(IN)  :: N2_ndensity_1d(npts) !m-3
    REAL(kind=8), INTENT(IN)  :: solar_declination_angle_degrees
    REAL(kind=8), INTENT(IN)  :: time(npts) ! seconds
    REAL(kind=8), INTENT(IN)  :: f107
    REAL(kind=8), INTENT(IN)  :: colatitude_geo(npts) ! radians
    REAL(kind=8), INTENT(IN)  :: R0 ! Radius of earth (m)

    REAL(kind=8), INTENT(IN)  :: temperature_1d(npts)
    REAL(kind=8), INTENT(IN)  :: gravity(npts)
    REAL(kind=8), INTENT(IN)  :: km_plasma(6) !=Boltzmann's const/(atomic_weight_X*1.673E-27)
    REAL(kind=8), INTENT(IN)  :: R0_eff(npts) ! radius of earth + altitude (m)
    REAL(kind=8), INTENT(IN)  :: eccentric

    REAL(kind=8), INTENT(OUT) :: ionisation_rates(npts,6)
    REAL(kind=8), INTENT(OUT) :: chi_degrees(npts)  ! solar zenith angle (degrees)


    INTEGER  :: n
    REAL(kind=8), PARAMETER :: PI    = 3.14159
    REAL(kind=8), PARAMETER :: DTR   = PI/180
    REAL(kind=8), PARAMETER :: PI2   = PI/2
    REAL(kind=8) :: cos_sda ! cos(solar declination angle in radians)
    REAL(kind=8) :: sin_sda ! sin(solar declination angle in radians)
    REAL(kind=8) :: cos_Hourangle    ! radians
    REAL(kind=8) :: csza(npts)  ! cos(solar zenith angle (radians))

    ! Night time ionisation factors
    REAL(kind=8) :: pngt , pngto2 , pngtn2
    PARAMETER (PNGT=1.0E-17, PNGTO2=2.0E-17,PNGTN2=1.8E-17)

    ! NWAVES is the total number of wavelength bands.
    ! The arrays are populated as
    ! Xrays in 1 -> NWAVES_XRAYS,
    ! EUV in NWAVES_XRAYS+1->NWAVES_EUV+NWAVES_XRAYS,
    ! SRC in NWAVES_EUV+NWAVES_XRAYS + 1->
    ! NWAVES_EUV+NWAVES_XRAYS+NWAVES_SRC

    ! NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY are returned from
    ! flux routines

    INTEGER :: NWAVES, NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY
    PARAMETER(NWAVES=200)

    ! wavelengths
    REAL(kind=8) :: WAVELS(NWAVES)

    ! ionisation and absorption cross sections in same bins as fluxes
    REAL(kind=8) :: CSAO(NWAVES),CSIO(NWAVES), CSAN2(NWAVES),           &
                    CSIN2(NWAVES), CSAO2(NWAVES), CSIO2(NWAVES)

    ! O2 EUV photon dissociation Xsection used if using TIEGCM EUVAC model
    REAL(kind=8) :: CSDO2(NWAVES)
    ! O2 EUV photoelectron dissociation Xsection used if using TIEGCM EUVAC model
    REAL(kind=8) :: CSDeO2(NWAVES)

    ! Generic flux array used for calculation
    REAL(kind=8) :: FLUX(NWAVES)
    REAL(kind=8) :: local_flux

    ! Lyman-a flux
    REAL(kind=8) :: lyman_a_flux
    INTEGER :: lyman_a_num

    ! Raw s2k fluxes and cross-sections in input file
    REAL(kind=8) :: S2K_ALLFLUXES(NWAVES,12), XSECTIONS(NWAVES,7)

    ! Plank*c. See below.
    REAL(kind=8) :: PCC
    ! column densities of level
    REAL(kind=8) :: WO, WO2, WN2

    ! Integrated column densities
    REAL(kind=8) :: TAUO, TAUO2, TAUN2, TAU, attenuation

    ! Some countv variables
    INTEGER :: wav, countv

    ! Consituent scale heights
    REAL(kind=8) :: HO(npts),HO2(npts),HN2(npts)

    ! Total height R0+ht1d for Chapmann function in cm
    REAL(kind=8) :: Z(npts)

    ! ionisation rate coefficients
    REAL(kind=8) :: PEUVO(npts), &
                    PEUVO2(npts),&
                    PEUVN2(npts)

    ! Constiuent densities cm-3
    REAL(kind=8) :: OALL(npts),O2ALL(npts),N2ALL(npts)

    ! seczenang for o2,o, n2
    REAL(kind=8) :: seco2,seco, secn2

    ! Ratio of nighttime to daytime ionization
    REAL(kind=8) :: rnight_o,rnight_o2, rnight_n2, nightfac

    ! Height between levels (cm)
    REAL(kind=8) :: zht

    ! local ionisation rates
    REAL(kind=8) :: PEUVO2_loc, PEUVO_loc, PEUVN2_loc

    REAL(kind=8) :: last_f107 = 0.

    ! Solar variability factor for SRB dissociation calculation
    REAL(kind=8) :: fmxfmn

    ! Save these variables between calls. Need this if not compiled with
    ! static flag.
    SAVE CSAO,CSIO, CSAN2,CSIN2, CSAO2, CSDO2, CSDeO2,                  &
         CSIO2, FLUX, WAVELS, lyman_a_flux,                             &
         lyman_a_num, NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY


    fmxfmn=(F107-71.)/(220-71.)

    ! Plank constant*c, with a units conversion factor (Jm * 1e6)
    PCC=1.985E-25*1.e6   !*1.e9

    ! Ratio of nightime ionisation to sec=1. This is to represent ionisation by
    ! Galactic cosmic rays
    !nightfac=1.e-6
    ! ald: taking out night facor and using CTIPe night time ionisation rates
    nightfac=0.

    ! Which fluxes/cross sections to use 1=CMAT1, 2=TIEGCM,
    ! >2 = SOLAR2000
    !DATA FLUX_SOURCE/1/


    ! --------------------- Initialise begin -----------------------------

    ! Get appropriate flux/cross-section data

    IF(f107 /= last_f107) THEN

      write(6,*) "F107 changed, from ",last_f107," to ",f107,             &
                 " :recalculating thermospheric fluxes."

      write(6,*) " Solar Fluxes from TIEGCM "
      CALL HL_ML__TIEGCM_FLUXES(f107, CSAO, CSIO, CSAO2, CSIO2,  &
         CSAN2, CSIN2, FLUX, WAVELS, Lyman_a_flux,        &
         NWAVES_XRAY, NWAVES_SRC, NWAVES_EUV, lyman_a_num,&
         CSDO2, CSDeO2)

    ENDIF

    !-----------------------------------------------------------------------


    ! Reset columns
    WO=0.
    WO2=0.
    WN2=0.
    ionisation_rates(:,:) = 0.


    ! Main height loop from top of atmosphere
    HEIGHTLOOP: DO n=low_level, top_level

       ! number densities cm-3
       Oall(n)=O_ndensity_1d(n)*1.e-6
       O2all(n)=O2_ndensity_1d(n)*1.e-6
       N2all(n)=N2_ndensity_1d(n)*1.e-6

       ! Scale heights, in cm

       HO = (temperature_1d(n)*km_plasma(1)/gravity(n))*1.E2
       HO2 = (temperature_1d(n)*km_plasma(5)/gravity(n))*1.E2
       HN2 = (temperature_1d(n)*km_plasma(4)/gravity(n))*1.E2

       ! Set total height in cm
       !R0_eff(n) = R0 + height_3d(n,m,l)
       z(n)=R0_eff(n)*1.E2

       PEUVO(n)  = 0.
       PEUVO2(n) = 0.
       PEUVN2(n) = 0.

       ! calculate sec zenith angle, incorporating Chapmann grazing incidence
       ! function

       sin_sda = sin(solar_declination_angle_degrees*DTR)
       cos_sda = cos(solar_declination_angle_degrees*DTR)

       cos_Hourangle = COS(PI*(TIME(n)/43200.-1.))
       csza(n) = sin_sda*cos(colatitude_geo(n))   &
                 + cos_sda*sin(colatitude_geo(n)) &
                 * cos_Hourangle
       chi_degrees(n) = (acos(csza(n)))/DTR

       SECO =CHAPMANN(csza(n), HO(n), Z(n))
       SECO2=CHAPMANN(csza(n), HO2(n), Z(n))
       SECN2=CHAPMANN(csza(n), HN2(n), Z(n))

       rnight_o=1.
       rnight_o2=1.
       rnight_n2=1.

      ! The stars are out ...
       IF(SECO < 0) THEN
          rnight_o=nightfac
          SECO=1.
       ENDIF
       IF(SECO2 < 0) THEN
          rnight_o2=nightfac
          SECO2=1.
       ENDIF
       IF(SECN2 < 0) THEN
          rnight_n2=nightfac
          SECN2=1.
       ENDIF

       ! Integrate Column densities

       WO = Oall(n)*HO(n)*SECO
       WO2= O2all(n)*HO2(n)*SECO2
       WN2= N2all(n)*HN2(n)*SECN2

       ! Loop over wavelength

       PEUVO2_loc  = 0.
       PEUVO_loc   = 0.
       PEUVN2_loc  = 0.

       TAU = 0.

       WAVELENGTH1: DO wav=1, NWAVES_EUV+NWAVES_XRAY+NWAVES_SRC


          TAUO = CSAO(wav)*WO
          TAUO2= CSAO2(wav)*WO2
          TAUN2= CSAN2(wav)*WN2
          TAU  = TAUO+TAUO2+TAUN2

          attenuation=0.
          IF (tau <= 200. .AND. tau /= 0.) attenuation = exp(-tau)
          IF (tau < 1.d-99 ) attenuation = 1.

          ! A tau of zero means all light gets through (none is absorbed)
          ! so attenuation = 1.

          local_flux=attenuation*FLUX(wav)

          !if(local_flux.le.0.)write(6,*)'local flux < 0 ' &
          !,n, wav, (acos(csza(n)))/DTR, attenuation, tau

          IF(local_flux < 1.e-20 ) local_flux=0.

          ! Calculate EUV ionization rate coefficients
          PEUVO_loc  = PEUVO_loc+local_flux*CSIO(wav)
          PEUVO2_loc = PEUVO2_loc+local_flux*CSIO2(wav)
          PEUVN2_loc = PEUVN2_loc+local_flux*CSIN2(wav)


       ENDDO WAVELENGTH1


       !Ionisation rates(m-3s-1)
       PEUVO2(n) = ((PEUVO2_loc*rnight_o2)+PNGTO2)*O2_ndensity_1d(n)
       PEUVO(n)  = ((PEUVO_loc*rnight_o)+PNGT)*O_ndensity_1d(n)
       PEUVN2(n) = ((PEUVN2_loc*rnight_n2)+PNGTN2)*N2_ndensity_1d(n)

       ionisation_rates(n,1) = PEUVO(n)
       ionisation_rates(n,4) = PEUVN2(n)
       ionisation_rates(n,5) = PEUVO2(n)

       if(PEUVO2(n).le.0.)write(6,*)'peuvo2 error ', &
                                n, (acos(csza(n)))/DTR, PEUVO2(n)
       if(PEUVO(n).le.0.)write(6,*)'peuvo error ', &
                                n, (acos(csza(n)))/DTR, PEUVO(n)
       if(PEUVN2(n).le.0.)write(6,*)'peuvn2 error ', &
                                n, (acos(csza(n)))/DTR, PEUVN2(n)


    ENDDO HEIGHTLOOP   ! End height loop

    last_f107 = f107

    RETURN

end SUBROUTINE HL_ML__EUV_ION_PRODUCTION_2


    !r========================================================
    !r=                   Chapmann                           =       
    !r========================================================
    !r   
    !r    Chapmann     Chapmann grazing incidence function from
    !r                 Smith and Smith 1972: JGR, V77, No. 19.
    !r                 Numerical Evauation of Chapmann's Grazing 
    !r                 Incidence Integral CH(X,SDA)
    !r                 (ALD 2005 replaced Risbeth and Garriot  with
    !r                 this as it's smoother)
    !r   
    !r    Passed in:
    !r    coschi       cos zenith angle
    !r    scale_ht     scht of constituent in question (cm)
    !r    ht           total altitude of point (cm)  
    !r   
    !r    Output :
    !r    Chapmann grazing incidence function
    !r   
    !r    Internal :
    !r    Y_ERR        Y in Smith and smith paper. Parameter for function
    !r    RAD_TO_Z     distance from centre of planet, normalised to scale height    
    !r    CHI          solar zenith angle in radians
    !r    CHID         solar zenith angle in degrees
    !r    ERFC         approxmation of complementary error function 
    !r                 (eq 12 and 13 in paper)

  REAL(kind=8) FUNCTION CHAPMANN(COSCHI,SCALE_HT,HT)

  REAL(kind=8), INTENT(IN) :: SCALE_HT, HT, COSCHI
  REAL(kind=8), PARAMETER :: PI    = 3.14159
  REAL(kind=8), PARAMETER :: PI2   = PI/2

  REAL*8 Y_ERR, RAD_TO_Z, CHI, CHID, ERFC

  CHI = ACOS(COSCHI)
  CHID = CHI*(180/PI)

  ! calculate chapman bit for angles over 75 degrees

  IF(chid.GT.75. .and. chid .LT. 96.) THEN
  !IF(chid.GT.75.) THEN

      RAD_TO_Z = HT/SCALE_HT

  ! Step1: calculate error function

      Y_ERR = SQRT(0.5*RAD_TO_Z) * ABS(COSCHI)
      if (Y_ERR.gt.100.0.and.HT.lt.(2*R0*1.e2))            &
        write(6,*) 'WARNING Problem in chapman function',  &
      Y_ERR,HT

      IF (Y_ERR.LE.8.0) THEN
        ERFC = (1.0606963 + 0.55643831* Y_ERR) /           &
               (1.0619898 + 1.7245609* Y_ERR +             &
               Y_ERR *Y_ERR)
      ELSE
        ERFC = 0.56498823 / (0.06651874 + Y_ERR)
      ENDIF

  ! step2. Calculate chapmann
  ! For solar zenith angles <= 90

      IF(CHID.LE.90.)THEN
         CHAPMANN = SQRT(PI2 * RAD_TO_Z) * ERFC
      ELSE

  ! For solar zenith angles > 90 (equation 15)
      CHAPMANN = SQRT(2*PI * RAD_TO_Z)*                     &
                 (SQRT(SIN(CHI))*EXP(RAD_TO_Z*(1-sin(CHI))) &
                 -0.5*ERFC)
  !         write(6,*)'chapman over 90', CHID,CHAPMANN
      ENDIF

  ELSE

  ! for angles below 80 degrees, use approximation that chapmann
  ! function is approx equal to sec(zenith angle)

      CHAPMANN = 1/COSCHI

  ENDIF

  END FUNCTION Chapmann


! **************************************************************
!   HL_ML__TIEGCM_FLUXES ROUTINE
! *************************************************************
! ALD Jan 2009
! Fluxes and cross sections taken from QRJ.F 
! (supplied by lqian@ucar.edu)
!
SUBROUTINE HL_ML__TIEGCM_FLUXES(F107, CSAO, CSIO, CSAO2, CSIO2,      &
                CSAN2, CSIN2, FLUX, WAVELS, Lyman_a_flux,       &
                NWAVES_XRAY_OUT, NWAVES_SRC_OUT, NWAVES_EUV_OUT,&
                lyman_a_num, CSDO2, CSDeO2)

    IMPLICIT NONE

    REAL(kind=8),INTENT(IN) :: f107  ! Daily value
    REAL(kind=8) F107a               ! 81 day average.
    REAL(kind=8) :: unit_conv = 1.e-18  ! Xsections into cm^2

    INTEGER   :: j, jinv
    INTEGER   :: NWAVES, NWAVES_XRAY, NWAVES_EUV, NWAVES_SRC
    PARAMETER (NWAVES = 37, NWAVES_XRAY = 0,     &
               NWAVES_SRC = 15, NWAVES_EUV = 22)

    ! FLUX, which is returned
    REAL(kind=8), INTENT(OUT) :: FLUX(NWAVES)

    ! cross sections (cm^2)
    REAL(kind=8), INTENT(OUT) :: CSAO(NWAVES),CSAO2(NWAVES), CSAN2(NWAVES)
    REAL(kind=8), INTENT(OUT) :: CSIO(NWAVES),CSIO2(NWAVES), CSIN2(NWAVES)
    REAL(kind=8), INTENT(OUT) :: CSDO2(NWAVES),CSDeO2(NWAVES)
    ! Wavelengths(m) (rlmeuv*1e-2)
    REAL(kind=8), INTENT(OUT) :: WAVELS(NWAVES)

    ! Lyman-a flux and record number
    REAL(kind=8), INTENT(OUT) :: lyman_a_flux
    INTEGER,      INTENT(OUT) :: lyman_a_num
    ! Output versions, to be returned to heating code
    INTEGER,      INTENT(OUT) :: NWAVES_EUV_OUT, NWAVES_SRC_OUT, &
                                 NWAVES_XRAY_OUT
    ! TIEGCM variables
    REAL(kind=8) sfmin(nwaves), afac(nwaves)
    REAL(kind=8) pind
    REAL(kind=8) rlmeuv(nwaves) ! wavelengths (cm)
    ! Total Absoption X sections (x1e18cm^2)
    REAL(kind=8) sigeuv_O(nwaves), sigeuv_O2(nwaves), sigeuv_N2(nwaves)
    ! ionization branching ratios (off absorption)
    REAL(kind=8) BphotonI_O(nwaves), BphotonI_O2(nwaves), BphotonI_N2(nwaves)
    ! O2 Photon dissociation branching ratios
    REAL(kind=8) bro2DPh(nwaves)
    ! O2 Photoelectron dissociation branching ratios
    REAL(kind=8) bro2Del(nwaves)

    DATA  rlmeuv/1.725e-05, 1.675e-05, 1.625e-05, 1.575e-05,  &
                 1.525e-05, 1.475e-05, 1.425e-05, 1.375e-05,  &
                 1.325e-05, 1.275e-05, 1.225e-05, 1.216e-05,  &
                 1.175e-05, 1.125e-05, 1.075e-05, 1.038e-05,  &
                 1.007e-05, 9.810e-06, 9.440e-06, 9.440e-06,  &
                 9.440e-06, 8.555e-06, 8.555e-06, 8.555e-06,  &
                 7.240e-06, 7.240e-06, 5.950e-06, 4.300e-06,  &
                 3.050e-06, 2.570e-06, 1.895e-06, 1.125e-06,  &
                 5.100e-07, 2.500e-07, 1.300e-07, 6.000e-08,  &
                 2.250e-08/

! O2 absorption coefficient:
    DATA  sigeuv_O2 /                                             &
                5.00e-01, 1.50e+00, 3.40e+00, 6.00e+00, 1.00e+01, &
                1.30e+01, 1.50e+01, 1.20e+01, 2.20e+00, 4.00e-01, &
                1.30e+01, 1.00e-02, 1.40e+00, 4.00e-01, 1.00e+00, &
                1.15e+00, 1.63e+00, 1.87e+01, 3.25e+01, 1.44e+01, &
                1.34e+01, 1.33e+01, 1.09e+01, 1.05e+01, 2.49e+01, &
                2.36e+01, 2.70e+01, 2.03e+01, 1.68e+01, 1.32e+01, &
                7.63e+00, 2.63e+00, 6.46e-01, 2.10e-01, 2.25e-01, &
                3.40e-02, 4.54e-03/

! O absorption coefficient:

    DATA  sigeuv_O /                                              &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 3.79e+00, 4.10e+00, 3.00e+00, 4.79e+00, &
                8.52e+00, 1.31e+01, 1.07e+01, 7.72e+00, 6.02e+00, &
                3.78e+00, 1.32e+00, 3.25e-01, 1.05e-01, 1.13e-01, &
                1.70e-02, 2.27e-03/

! N2 absorption coefficient:

    DATA  sigeuv_N2 /                                             &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 2.55e+00, 1.15e+02, 1.44e+01, &
                2.18e+00, 7.17e+01, 1.31e+01, 2.14e+00, 5.45e+01, &
                2.30e+01, 2.31e+01, 1.97e+01, 1.17e+01, 9.94e+00, &
                5.09e+00, 1.53e+00, 3.46e-01, 1.14e+00, 1.41e-01, &
                2.01e-02, 2.53e-03/

! The three major species' ionization branching ratio (off absorption):


! O2
    DATA  BPhotonI_O2 /                                           &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 6.13e-01, 8.30e-01, 6.20e-01, 7.86e-01, &
                7.56e-01, 5.34e-01, 5.74e-01, 5.49e-01, 4.76e-01, &
                6.73e-01, 9.83e-01, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00/
! O
    DATA  BPhotonI_O /                                            &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00/
! N2
    DATA  BPhotonI_N2 /                                           &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 4.29e-01, &
                6.80e-01, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00/

! O2 photon dissociation branching ratio
    DATA  bro2DPh  /                                              &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, &
                1.00e+00, 3.87e-01, 1.70e-01, 3.80e-01, 2.14e-01, &
                2.44e-01, 4.66e-01, 4.26e-01, 4.51e-01, 5.24e-01, &
                3.27e-01, 1.74e-02, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00/


! O2 photoelectron dissociation branching ratio
    DATA  bro2DEl  /                                              &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, &
                0.00e+00, 1.10e-02, 6.53e-01, 7.62e-01, 9.96e-01, &
                1.27e+00, 2.04e+00, 4.11e+00, 5.70e+01, 1.78e+01, &
                2.03e+01, 8.79e+01/


! Solar spectrum based on EUVAC and glow for wave length less than 1050 A
! and Woods for wavelength greater than 1050 A

! solar minimum flux (when P_index=80, unit:photon cm^-2 S^-1)

    DATA  sfmin /3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10,  &
                 5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10,  &
                 2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11,  &
                 8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09,  &
                 4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09,  &
                 3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09,  &
                 1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09,  &
                 8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09,  &
                 5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04,  &
                 5.010e+01/

! scaling factor A as defined in EUVAC model

    DATA  afac /5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03,   &
                1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03,   &
                2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03,   &
                2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03,   &
                5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03,   &
                4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03,   &
                3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02,   &
                7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03,   &
                1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01,   &
                6.240e-01/

    write(6,*)'Using Solar fluxes etc from TIEGCM: ', f107
    write(6,*)''

     F107a = F107 ! set 81 day average to be same as daily value

! solar model is the same as EUVAC, i.e., solar flux at
! (f107,f107a) is scale as: sflux=sfmin(1+afac(P-80.))
! where P=0.5(f107+f107a)

      pind=0.5*(f107+f107a)

! bins

      do J = 1, NWAVES

        jinv = nwaves-J+1
        WAVELS(J) = rlmeuv(jinv)*1.e-2 ! (convert to m-1)
        CSAO(J)   = sigeuv_O(jinv)*unit_conv
        CSAO2(J)  = sigeuv_O2(jinv)*unit_conv
        CSAN2(J)  = sigeuv_N2(jinv)*unit_conv
        CSIO(J)   = CSAO(J)*BphotonI_O(jinv)
        CSIO2(J)  = CSAO2(J)*BphotonI_O2(jinv)
        CSIN2(J)  = CSAN2(J)*BphotonI_N2(jinv)
        CSDO2(J)  = CSAO2(J)*bro2Dph(jinv)
        CSDeO2(J) = CSIO2(J)*bro2Del(jinv)

        flux(j)=sfmin(jinv)*(1+afac(jinv)*(pind-80.))

        IF (FLUX(J) .LT. 0.0) FLUX(J) = 0.0
        if(wavels(j) == 1.216e-07) then
                        lyman_a_num  = j
                        lyman_a_flux = flux(j)
        endif

        WRITE(6,"(4E14.6,A)") WAVELS(j), CSIO(j), &
                 CSIO2(j),CSIN2(j), "  tiegcm"
        WRITE(6,"(4E14.6,A)") WAVELS(j), CSAO(j), &
                 CSAO2(j),CSAN2(j), "  tiegcm"
        WRITE(6,"(I4, 2E14.6, A)") j, WAVELS(j), &
                 FLUX(j), "  tiegcm"

      enddo

  NWAVES_XRAY_out = nwaves_xray
  NWAVES_SRC_out  = nwaves_src
  NWAVES_EUV_out  = nwaves_euv



end SUBROUTINE HL_ML__TIEGCM_FLUXES





SUBROUTINE HL_ML__INTERP_PRESSURE_TO_TUBE(n_array_size,i1,i2,pz,param_pressure,param_tube)
!g
!g Interpolates a parameter in height from pressure coordinates to tube coordinates....
!g
  INTEGER :: n_array_size , i1 , i2
  INTEGER :: nhgt ,  i, i_pres, iupper, ilower
  REAL(kind=8) ::  ht(15), param_pressure(15), pz(n_array_size), param_tube(n_array_size)

ht(1)= 80.12440678759815
ht(2)= 86.34782667063907
ht(3)= 92.30040653307229
ht(4)= 99.1644693936756
ht(5)= 106.92273205414112
ht(6)= 115.87573705778922
ht(7)= 126.54786377322566
ht(8)= 142.36483990939157
ht(9)= 164.12432836837356
ht(10)= 193.81772599760085
ht(11)= 226.0844005019562
ht(12)= 263.9841685189012
ht(13)= 305.46160600421796
ht(14)= 349.51432004724006
ht(15)= 395.2760145837757
!g
!g loop over the points on the tube....
!g
  do 100 i = i1 , i2
  !g
  !g   find which points in pressure are above and below our chosen tube point....
  !g
      iupper=16
      do 200 i_pres = 1, 15

          ! Changed . gt . .gt. mjh
          if (ht(i_pres) .gt. pz(i)) then
              iupper = i_pres
              ilower = i_pres - 1
              goto  300
          endif

      200 ENDDO
  !g
  !g if we get to here - then pz is higher than ht(15) (the highest pressure level)
  !g Setting both iupper and ilower to 15 will make all values of param_tube above ht(15)
  !g equal to param_pres(15)....
  !g

      300 continue
  !g
  !g  we jump to here when we have found the upper and lower points.
  !g  Now do the interpolation.....
      if (iupper <= 15) then
          param_tube(i) = (((pz(i) - ht(ilower))/(ht(iupper) - ht(ilower))) * &
          (param_pressure(iupper) - param_pressure(ilower))) + param_pressure(ilower)
      else
          param_tube(i) = param_pressure(15)
      endif

  ! write(6,*) 'param ',i,param_tube(i)
  100 ENDDO

end SUBROUTINE HL_ML__INTERP_PRESSURE_TO_TUBE






SUBROUTINE HL_ML__OPTICAL_DEPTH(CHI,r0,tn,gravity,km,gr,o,o2,n2,EMTAU)

!***********************************************************************
! CALLED BY ROUTINE HL_ML__EUV_ION_PRODUCTION, CALCULATES THE OPTICAL DEPTH
! CROSS SECTIONS FROM TORR AND TORR
! (TORR AND TORR, 1982, REV. GEOPHYS. SPACE RES., 20, 91-144)

! CHAPMAN FUNCTION APPROXIMATION FROM RISHBETH AND GARRIOTT
! (INTRODUCTION TO IONOSPHERIC PHYSICS)
!***********************************************************************


  implicit none
  integer :: mm,k
  REAL(kind=8) :: pye,r0,dtr,tn,o,o2,n2,gr,gravity,km(6), &
  chap,chi,colum(3),snchi,cs2chi,w1
  REAL(kind=8) :: csao,csao2,csan2,pyeh,chid,w,ho,ho2,hn2,seco,seco2,secn2
  REAL(kind=8) :: ferf,wwo,w2,w3,w4,wwn2,tau,emtau,wwo2,cschi
  REAL(kind=8) :: height
  !Commented out mjh

  parameter(pye=3.14159)
  parameter(dtr=pye/180.)
  SAVE CSAO,CSAO2,CSAN2,PYEH,MM
  REAL(kind=8) :: M,NO
  DIMENSION CSAO(37),CSAO2(37),CSAN2(37),EMTAU(37) 

!***********************************************************************
! ALL CROSS SECTIONS ARE GIVEN IN UNITS 10**-22 M2
!***********************************************************************

!--PHOTOABSORPTION CROSS SECTIONS - O
  DATA CSAO/ &
  1.06, 3.53, 5.96, 7.55, 8.43, 9.26, 8.78, 9.70, &
  9.72,10.03,10.84,10.70,11.21,11.25,11.64,11.91, &
  12.13,12.17,11.90,12.23,12.22,12.21,10.04,11.35, &
  8.00, 4.18, 4.18, 4.28, 4.23, 4.38, 4.18, 2.12, &
  0.00, 0.00, 0.00, 0.00, 0.00/

!--PHOTOABSORPTION CROSS SECTIONS - O2
  DATA CSAO2/ &
  1.18, 3.61, 7.27,10.50,12.80,14.80,13.65,15.98, &
  16.00,17.19,18.40,18.17,19.39,20.40,21.59,24.06, &
  25.59,22.00,25.04,26.10,25.80,26.02,26.27,25.00, &
  29.05,21.98,25.18,26.66,27.09,20.87, 9.85,15.54, &
  4.00,16.53, 1.60, 1.00, 1.10/

!--PHOTOABSORPTION CROSS SECTIONS - N2
  DATA CSAN2/ &
  0.60, 2.32, 5.40, 8.15, 9.65,10.60,10.08,11.58, &
  11.60,14.60,18.00,17.51,21.07,21.80,21.85,24.53, &
  24.69,23.20,22.38,23.10,23.20,23.22,29.75,26.30, &
  30.94,35.46,26.88,19.26,30.71,15.05,46.63,16.99, &
  0.70,36.16, 0.00, 0.00, 0.00/

  DATA MM/0/

!--FIRST TIME SET THE IONIZATION & ABSORPTION CROSS
!--SECTIONS TO THE CORRECT MAGNITUDE

  IF(MM == 0) THEN
      DO 1 K=1,37
          CSAO(K)=CSAO(K)*1.E-22
          CSAO2(K)=CSAO2(K)*1.E-22
          CSAN2(K)=CSAN2(K)*1.E-22
      1 ENDDO
      PYEH=PYE*.5
      MM=1
  ENDIF

!--SET HO, HO2, HN2 EQUAL TO THE SCALE HEIGHTS FOR O, O2 AND N2 RESP.

  CHID=CHI/dtr
  CSCHI = COS(CHI)
  IF(CHID <= 115.0) THEN
      W=TN/Gravity
  ! write(6,*) 'yowzer ',i,tn,gravity,km(1)
      HO=W*KM(1)
      HO2=W*KM(5)
      HN2=W*KM(4)
      height = (gr-r0)/1000.

      IF((CHID <= 70.) .OR. (chid <= 89 .AND. height <= 140.)) THEN
          SECO=1./CSCHI
          SECO2=SECO
          SECN2=SECO

      !g
      !g  there is a bug here somewhere.   a point with chid = 113
      !g  and a height of 103km  has nowhere to go.. changed the next elseif
      !g  line to a simple else....
      !g
      ! ELSEIF(height.gt.140.) then
      ELSE

      !--CALCULATE THE CHAPMAN FUNCTION WHICH REPLACES SEC(CHI) IN THE
      !--PRODUCTION FUNCTION WHEN CHID IS GREATER THAN 70 DEGREES (PUT THESE
      !--VALUES INTO THE VARIABLES SECO, SECO2 AND SECN2)

          SNCHI=SIN(CHI)
          CS2CHI=COS(CHI)*COS(CHI)
          W1=gr
          IF(CHID > 90.) THEN
              FERF=1.
          ELSE
              FERF=-1.
          ENDIF

          WWO=W1/HO
          W2=SQRT(PYEH*WWO*SNCHI)
          W3=WWO*CS2CHI*.5
          W4=SQRT(W3)
      ! write(6,*) 'optical 3.91',i,wwo,ho
          SECO=W2*EXP(W3)*(1.+FERF*ERFCHE(W4))
      ! write(6,*) 'here optical 3.92'
      ! if(ERFCHE(W4).eq.1.0) seco = 1.e20

          WWO2=W1/HO2
          W2=SQRT(PYEH*WWO2*SNCHI)
          W3=WWO2*0.5*CS2CHI
          W4=SQRT(W3)
          SECO2=W2*EXP(W3)*(1.+FERF*ERFCHE(W4))
      ! if(ERFCHE(W4).eq.1.0) seco2 = 1.e20

          WWN2=W1/HN2
          W2=SQRT(PYEH*WWN2*SNCHI)
          W3=WWN2*0.5*CS2CHI
          W4=SQRT(W3)
          SECN2=W2*EXP(W3)*(1.+FERF*ERFCHE(W4))
      ! if(ERFCHE(W4).eq.1.0) secn2 = 1.e20

      ENDIF

  !***********************************************************************
  ! HAVING CALCULATED SECO, SECO2, SECN2 CONTINUE FOR ALL SOLAR ZENITH
  ! ANGLES LESS THEN 100 DEGREES
  !***********************************************************************

      COLUM(1)=O*HO*SECO
      CHAP=SECO
      COLUM(2)=O2*HO2*SECO2
      COLUM(3)=N2*HN2*SECN2

      DO 20 K=1,37
          TAU=CSAO(K)*COLUM(1)+CSAO2(K)*COLUM(2) &
          +CSAN2(K)*COLUM(3)
      !g
      !g I have set this up so that we only compute emtau if tau is less than
      !g 200.  e-200 is 10-87 so this is ok.  This came from the HP at SEL
      !g crashing because it couldn't do e-710 (my calculator gives the right
      !g answer (ie zero).   - George 12th June 1995
      !g
      !g same with the 1.e-99 bit - what is it with this computer - can't it do
      !g maths ????
      !g
          IF(TAU > 200.) THEN
              EMTAU(K)=0.
          ELSEIF(TAU < 1.d-99) THEN
              EMTAU(K)=1.
          ELSE
              EMTAU(K)=EXP(-TAU)
          ENDIF
      20 ENDDO

  !-- COLUM CONVERTED TO CGS UNITS FOR USE IN PHOTOELECTRON HEATING
  !--   CALCULATIONS
      COLUM(1)=COLUM(1)*1.E-4
      COLUM(2)=COLUM(2)*1.E-4
      COLUM(3)=COLUM(3)*1.E-4
  ELSE
  !..... a negative value of COLUM(1) means nighttime conditions in QE calcs.
      COLUM(1)=-1.
      CHAP=1.0E+10
      DO 25 K=1,37
          EMTAU(K)=0.
      25 ENDDO
  ENDIF

  RETURN




end SUBROUTINE HL_ML__OPTICAL_DEPTH






SUBROUTINE ML__MID_AND_LOW_LATITUDE_IONOSPHERE( &
  GIP_switches, &
  geo_grid_longitudes_degrees,geo_grid_longitudes_radians, &
  geo_grid_latitudes_degrees,geo_grid_latitudes_radians, &
  geo_grid_colatitudes_degrees,geo_grid_colatitudes_radians, &
  iday_number, &
  i_first_call_of_plasma,UT_in_seconds,plasma_time_step_seconds, &
  idump,iout,iout_high_res,ipcall,ipint,ipint_high, &
  i_no_day,i_graphics_out_start, &
  file_res, &
  sw_initialisation_call, &
  F107,Solar_Declination_Angle_degrees, &
  ETA_Apex_3D, &
  Apex_D1,Apex_D2, &
  Apex_E1,Apex_E2, &
  Apex_BE3,Apex_bhat,Apex_Bmag, &
  Apex_D,Apex_d1d1,Apex_d1d2,Apex_d2d2,Apex_grdlbm2, &
  bcol,blon,q_coordinate,Re_apex, &
  gr,gcol,glon, &
  IN,IS, &
  glat_plasma_3d,glond_plasma_3d,pz_plasma_3d, &
  TN_plasma_input_3d,O_plasma_input_3d,O2_plasma_input_3d,N2_plasma_input_3d, &
  NO_plasma_input_3d,N4S_plasma_input_3d, N2D_plasma_input_3d, &
  te_dum_plasma_input_3d, &
  um_plasma_input_3d,uz_plasma_input_3d,uv_plasma_input_3d, &
  Temp_ion,Temp_electron,NI,VI,NE,YNI,YVI,YNE,YQE,no_plus_3d,o2_plus_3d, &
  n2_plus_3d,n_plus_3d, &
  JIOn,mass,ndt,nions,midpoint,itday,qb,KM_plasma,M_plasma, &
  dynamo_sigma_phph_dsi,dynamo_sigma_lmlm_msi, &
  dynamo_sigma_h,dynamo_sigma_c, &
  dynamo_Kdmph_dsi,dynamo_Kdmlm, &
  potential_field,ed1_from_dynamo,ed2_from_dynamo, &
  vpeq,vzon)

!*********************************************
! *
! sheffield mid and low latitude 'global'  *
! plasmasphere/ionosphere program.      *
! *
!*********************************************

  IMPLICIT NONE

    INTEGER N_Pressure_Levels,N_Latitudes,N_longitudes
      PARAMETER(N_Pressure_Levels=15)
      PARAMETER(N_Latitudes=91)
      PARAMETER(N_Longitudes=20)

  LOGICAL :: sw_initialisation_call
  LOGICAL :: GIP_switches(20)
  LOGICAL :: sw_External_model_provides_low_lat_E_fields 
  INTEGER :: IN(NMP,NLP) , IS(NMP,NLP)
  INTEGER :: IN_dum(NMP,NLP) , IS_dum(NMP,NLP)
  INTEGER :: i , j
  INTEGER :: iday_number
  INTEGER :: istop
  INTEGER :: IDUmp
  INTEGER :: ifailed
  INTEGER :: iout
  INTEGER :: iout_high_res
  INTEGER :: iseasav , iutav
  INTEGER :: ite_failed
  INTEGER :: iti_failed
  INTEGER :: ipint
  INTEGER :: ipint_high
  INTEGER :: file_res(27)
  INTEGER :: i_no_day
  INTEGER :: i_graphics_out_start
  INTEGER :: iwrite
  INTEGER :: JIOn(6)
  INTEGER :: mp , lp
  INTEGER :: MASs
  INTEGER :: NDT
  INTEGER :: NIOns
  INTEGER :: midpoint(nlp)
  INTEGER :: i_first_call_of_plasma
  INTEGER :: ITDay
  INTEGER :: in1
  INTEGER :: in2
  INTEGER :: N_DYN_LAT
  integer :: ilon_dynamo
  PARAMETER (N_DYN_LAT=97)     ! latitude grid used by the new apex dynamo solver
  INTEGER :: ilat_dynamo
  INTEGER :: i_call_the_electrodynamics
  INTEGER :: IPCall
  INTEGER :: i300(nmp,nlp)
  INTEGER :: Ofailed(nmp,nlp)
  INTEGER :: Hfailed(nmp,nlp)
  INTEGER :: i_create_interpolation_indexes
  INTEGER :: i_create_high_res_indexes
  INTEGER :: ilon,iheight,ipoint
  INTEGER :: i_failed_molecular_ions
  INTEGER :: iwrite98
  INTEGER :: idiagnose
  INTEGER :: ilt
  INTEGER :: ilt_east
  INTEGER :: ilt_west
  INTEGER :: ilat
  INTEGER :: ilat_north
  INTEGER :: ilat_south
  REAL(kind=8) factor_lt
  REAL(kind=8) factor_lat
  REAL(kind=8) H_NE
  REAL(kind=8) H_NW
  REAL(kind=8) H_SE
  REAL(kind=8) H_SW
  REAL(kind=8) H_N
  REAL(kind=8) H_S
  REAL(kind=8) H_400
  REAL(kind=8) dayno
  REAL(kind=8) plasma_time_step_seconds
  REAL(kind=8) DTR
  REAL(kind=8) F107
  REAL(kind=8) gr300
  REAL(kind=8) plasma_energy_time_step_seconds
  REAL(kind=8) e2_factor_from_300km_to_apex
  REAL(kind=8) apex_height_km
  REAL(kind=8) geolon_apex_degrees
  REAL(kind=8) QB
  REAL(kind=8) R0
  REAL(kind=8) UT_in_seconds
  REAL(kind=8) UT_in_hours
  REAL(kind=8) w1
  REAL(kind=8) w2
  REAL(kind=8) wn2
  REAL(kind=8) wo
  REAL(kind=8) wo2
  REAL(kind=8) xmlat_300
  REAL(kind=8) xmlon_300
  REAL(kind=8) vi_limited
  REAL(kind=8) PI
  REAL(kind=8) Solar_Declination_Angle_degrees
  REAL(kind=8) atomic_mass_unit
  REAL(kind=8) RE_apex(NMP,NLP)
  REAL(kind=8) Burnside_factor
  REAL(kind=8) ne_save(npts)
  REAL(kind=8) ni_save(npts,nmp,2)
  REAL(kind=8) vi_save(npts,nmp,2)
  REAL(kind=8) ti_save(npts,nmp,2)
  REAL(kind=8) te_save(npts,nmp)
  REAL(kind=8) U_merid(npts)
  REAL(kind=8) U_zonal(npts)
  REAL(kind=8) U_vertical(npts)
  REAL(kind=8) mag_lat_tubeN(nmp,nlp)
  REAL(kind=8) mag_lat_tubeS(nmp,nlp)
  REAL(kind=8) xtimes(npts)
  REAL(kind=8) Apex_D1(3,NPTS,NMP)
  REAL(kind=8) Apex_D2(3,NPTS,NMP)
  REAL(kind=8) Apex_E1(3,NPTS,NMP)
  REAL(kind=8) Apex_E2(3,NPTS,NMP)
  REAL(kind=8) Apex_BE3(NPTS,NMP)
  REAL(kind=8) Apex_grdlbm2(3,NPTS,NMP)
  REAL(kind=8) ETA_APEX_3D(NPTS,NMP)
  REAL(kind=8) sigma_phph_dsi(2,NMP,NLP)
  REAL(kind=8) sigma_lmlm_msi(2,NMP,NLP)
  REAL(kind=8) sigma_h(2,NMP,NLP)
  REAL(kind=8) sigma_c(2,NMP,NLP)
  REAL(kind=8) Kdmph_dsi(2,NMP,NLP)
  REAL(kind=8) Kdmlm(2,NMP,NLP)
  REAL(kind=8) sigma_phph_dsi_1d(2)
  REAL(kind=8) sigma_lmlm_msi_1d(2)
  REAL(kind=8) sigma_h_1d(2)
  REAL(kind=8) sigma_c_1d(2)
  REAL(kind=8) Kdmph_dsi_1d(2)
  REAL(kind=8) Kdmlm_1d(2)
  REAL(kind=8) int_sigma_ped(nmp,nlp)
  REAL(kind=8) int_sigma_ped2
  REAL(kind=8) dynamo_sigma_phph_dsi(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_sigma_lmlm_msi(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_sigma_h(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_sigma_c(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_Kdmph_dsi(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_Kdmlm(NMP+1,N_DYN_LAT)
  REAL(kind=8) potential_field(NMP+1,N_DYN_LAT)
  REAL(kind=8) ed1_from_dynamo(NMP+1,N_DYN_LAT)
  REAL(kind=8) ed2_from_dynamo(NMP+1,N_DYN_LAT)
  REAL(kind=8) dynamo_latitude(N_DYN_LAT)
  REAL(kind=8) potential_on_tubes(nmp,nlp)
  REAL(kind=8) ed1_on_tube(nmp,nlp)
  REAL(kind=8) ed2_on_tube(nmp,nlp)
  REAL(kind=8) ve1
  REAL(kind=8) ve2
  REAL(kind=8) V_upwards_at_apex
  REAL(kind=8) V_upwards_at_apex_array(nmp,nlp)
  REAL(kind=8) V_upwards_at_apex_empirical
  REAL(kind=8) local_time_apex(nmp,nlp)
  REAL(kind=8) magnitude_e2_at_300
  REAL(kind=8) magnitude_e2_at_apex
  REAL(kind=8) Apex_BHAT(3,npts,nmp)
  REAL(kind=8) Apex_D(npts,nmp)
  REAL(kind=8) Apex_d1d1(npts,nmp)
  REAL(kind=8) Apex_d1d2(npts,nmp)
  REAL(kind=8) Apex_d2d2(npts,nmp)
  REAL(kind=8) Apex_BMAG(npts,nmp)
  REAL(kind=8) f_1d_parameter(NPTS)
  REAL(kind=8) df_1d_parameter(NPTS)
  REAL(kind=8) o(NPTS)
  REAL(kind=8) o2(NPTS)
  REAL(kind=8) n2(NPTS)
  REAL(kind=8) hyd(NPTS)
  REAL(kind=8) neutral_hydrogen_1d(NPTS)
  REAL(kind=8) hel(NPTS)
  REAL(kind=8) chi_degrees(NPTS)
  REAL(kind=8) cos_SZA(NPTS)
  REAL(kind=8) no_plus(NPTS)
  REAL(kind=8) o2_plus(NPTS)
  REAL(kind=8) n_plus(NPTS)
  REAL(kind=8) n2_plus(NPTS)
  REAL(kind=8) no_plus_saved(NPTS)
  REAL(kind=8) o2_plus_saved(NPTS)
  REAL(kind=8) n_plus_saved(NPTS)
  REAL(kind=8) n2_plus_saved(NPTS)
  REAL(kind=8) o_plus_1d_saved(NPTS)
  REAL(kind=8) no_plus_3d(npts,nmp)
  REAL(kind=8) o2_plus_3d(npts,nmp)
  REAL(kind=8) n2_plus_3d(npts,nmp)
  REAL(kind=8) n_plus_3d(npts,nmp)
  REAL(kind=8) nuin(NPTS)
  REAL(kind=8) fh(NPTS,2)
  REAL(kind=8) chemp(NPTS,2)
  REAL(kind=8) beta(NPTS,2)
  REAL(kind=8) peuvi(NPTS,2)
  REAL(kind=8) peuvn(NPTS,6)
  REAL(kind=8) gravity(NPTS)
  REAL(kind=8) o_plus_1d(NPTS)
  REAL(kind=8) ti_1d(npts)
  REAL(kind=8) qin(npts)
  REAL(kind=8) RNE(npts)
  REAL(kind=8) gr(NPTS,nmp)
  REAL(kind=8) gcol(NPTS,nmp)
  REAL(kind=8) glon(NPTS,nmp)
  REAL(kind=8) bcol(NPTS,NMP)
  REAL(kind=8) BLOn(NMP,NLP)
  REAL(kind=8) glat(NPTS)
  REAL(kind=8) altitude_pz_km(NPTS)
  REAL(kind=8) glond(NPTS)
  REAL(kind=8) time(NPTS)
  REAL(kind=8) angle(NPTS)
  REAL(kind=8) g_parallel(NPTS)
  REAL(kind=8) dte(NPTS)
  REAL(kind=8) dti(NPTS,2)
  REAL(kind=8) q_coordinate(NPTS,NMP)
  REAL(kind=8) upar_zonal_contribution(NPTS)
  REAL(kind=8) upar_meridional_contribution(NPTS)
  REAL(kind=8) upar_vertical_contribution(NPTS)
  REAL(kind=8) Temp_ion(NPTS,NMP,2)
  REAL(kind=8) Temp_electron(NPTS,NMP)
  REAL(kind=8) NI(NPTS,NMP,2)
  REAL(kind=8) dni(NPTS,2)
  REAL(kind=8) VI(NPTS,NMP,2)
  REAL(kind=8) NE(NPTS,NMP)
  REAL(kind=8) bb(NPTS)
  REAL(kind=8) div_vperp_3d(NPTS,NMP)
  REAL(kind=8) div_vperp_dipole(3,NPTS)
  REAL(kind=8) vp(npts)
  REAL(kind=8) M_plasma(6)
  REAL(kind=8) KM_plasma(6)
  REAL(kind=8) YNI(NPTS,NMP,2)
  REAL(kind=8) YVI(NPTS,NMP,2)
  REAL(kind=8) YNE(NPTS,NMP)
  REAL(kind=8) YQE(NPTS,NMP)
  REAL(kind=8) N1(NPTS)
  REAL(kind=8) nn2(NPTS)
  REAL(kind=8) U2DIF(NPTS)
  REAL(kind=8) TT(NPTS)
  REAL(kind=8) te_dum(npts)
  REAL(kind=8) dni_save(npts)
  REAL(kind=8) upar_apex(npts)
  REAL(kind=8) uperp_apex(npts)
  REAL(kind=8) Um_apex(npts)
  REAL(kind=8) vpeq(nmp,nlp)
  REAL(kind=8) vzon(nmp,nlp)
  REAL(kind=8) nuin_tot(npts)
  REAL(kind=8) nu_mol_n(npts)
  REAL(kind=8) ion_mass(npts)
  REAL(kind=8) e1(2,nmp,nlp)
  REAL(kind=8) e2(2,nmp,nlp)
  REAL(kind=8) vzon_coupled(nmp,nlp)
  real(kind=8) Tn_plasma_input_3d(npts,nmp)
  real(kind=8) O_plasma_input_3d(npts,nmp)
  real(kind=8) O2_plasma_input_3d(npts,nmp)
  real(kind=8) N2_plasma_input_3d(npts,nmp)

  real(kind=8) NO_plasma_input_3d(npts,nmp)
  real(kind=8) N4S_plasma_input_3d(npts,nmp)
  real(kind=8) N2D_plasma_input_3d(npts,nmp)
  real(kind=8) Glat_plasma_3d(npts,nmp)
  real(kind=8) Glond_plasma_3d(npts,nmp)
  real(kind=8) Pz_plasma_3d(npts,nmp)
  real(kind=8) Um_plasma_input_3d(npts,nmp)
  real(kind=8) Uz_plasma_input_3d(npts,nmp)
  real(kind=8) Uv_plasma_input_3d(npts,nmp)
  real(kind=8) Te_dum_plasma_input_3d(npts,nmp)
  REAL(kind=8) geo_grid_longitudes_degrees(N_longitudes)
  REAL(kind=8) geo_grid_longitudes_radians(N_longitudes)
  REAL(kind=8) geo_grid_latitudes_degrees(N_Latitudes)
  REAL(kind=8) geo_grid_latitudes_radians(N_Latitudes)
  REAL(kind=8) geo_grid_colatitudes_degrees(N_Latitudes)
  REAL(kind=8) geo_grid_colatitudes_radians(N_Latitudes)
  REAL(kind=8) euv_production_rate_Oplus(npts)
  REAL(kind=8) euv_production_rate_N2plus(npts)
  REAL(kind=8) euv_production_rate_O2plus(npts)
  REAL(kind=8) auroral_production_rate_Oplus(npts)
  REAL(kind=8) auroral_production_rate_N2plus(npts)
  REAL(kind=8) auroral_production_rate_O2plus(npts)

  REAL(kind=8) Apex_D1_1d(3,npts)
  REAL(kind=8) Apex_D2_1d(3,npts)
  REAL(kind=8) Apex_BE3_1d(npts)
  REAL(kind=8) Apex_D_1d(npts)
  REAL(kind=8) Apex_d1d1_1d(npts)
  REAL(kind=8) Apex_d1d2_1d(npts)
  REAL(kind=8) Apex_d2d2_1d(npts)
  REAL(kind=8) Apex_BMAG_1d(npts)
  REAL(kind=8) apex_bhat_1d(3,npts)
  REAL(kind=8) gr_1d(npts)
  REAL(kind=8) gcol_1d(npts)
  REAL(kind=8) glon_1d(npts)
  REAL(kind=8) tn_1d(npts)
  REAL(kind=8) o_1d(npts)
  REAL(kind=8) o2_1d(npts)
  REAL(kind=8) n2_1d(npts)

  REAL(kind=8) HYD_1d(npts)
  REAL(kind=8) HEL_1d(npts)
  REAL(kind=8) NO_1d(npts)
  REAL(kind=8) N4S_1d(npts)
  REAL(kind=8) N2D_1d(npts)
  REAL(kind=8) U_merid_1d(npts)
  REAL(kind=8) U_zonal_1d(npts)
  REAL(kind=8) U_vertical_1d(npts)
  REAL(kind=8) Te_dum_1d(npts)
  REAL(kind=8) eta_apex_1d(npts)
  REAL(kind=8) no_plus_1d(npts)
  REAL(kind=8) o2_plus_1d(npts)
  REAL(kind=8) n2_plus_1d(npts) 
  REAL(kind=8) n_plus_1d(npts)
  REAL(kind=8) div_vperp_1d(npts)
  REAL(kind=8) q_coordinate_1d(NPTS)

  REAL(kind=8) ni_oplus_1d(npts)
  REAL(kind=8) ni_hplus_1d(npts)
  REAL(kind=8) vi_oplus_1d(npts)
  REAL(kind=8) vi_hplus_1d(npts)
  REAL(kind=8) ne_1d(npts)
  REAL(kind=8) dni_oplus_1d(npts)
  REAL(kind=8) dni_hplus_1d(npts)
  REAL(kind=8) dne_1d(npts)
  REAL(kind=8) ti_oplus_1d(npts)
  REAL(kind=8) ti_hplus_1d(npts)
  REAL(kind=8) te_1d(npts)
  integer i_write_out_tube
  REAL(kind=8) te_tn_factor
  REAL(kind=8) glat_rad
  REAL(kind=8) geo_local_time_degrees
  REAL(kind=8) rlt
  REAL(kind=8) solar_zenith_angle_degrees
  LOGICAL :: morning

  real(kind=8) V_exb_east(npts,nmp)
  real(kind=8) V_exb_north(npts,nmp)
  real(kind=8) V_exb_up(npts,nmp)

  INTEGER i_attempt
  INTEGER i_attempt_array(nmp,nlp)

!g
  PARAMETER (PI=3.141592654,DTR=PI/180.0,R0=6.370E06)
  parameter (atomic_mass_unit=1.66e-27)

  idiagnose = 0

  if(sw_initialisation_call) then

  ! INITIALISATION PART........


      NIOns = 2
      JIOn(1) = 2
      JIOn(2) = 2
      JIOn(3) = 0
!       midpoint = (NPTS+1)/2
      ITDay = 1
      QB = 3.E-18
      Mass = 48
      NDT = 10



      CALL IO__read_apex_mag_field_coords ( &
        ETA_Apex_3D, &
        Apex_D1,Apex_D2, &
        Apex_E1,Apex_E2, &
        Apex_BE3,Apex_bhat,Apex_Bmag, &
        Apex_D,Apex_d1d1,Apex_d1d2,Apex_d2d2,Apex_grdlbm2, &
        bcol,blon,q_coordinate,Re_apex, &
        gr, gcol, glon, &
        IN, IS, &
        glat_plasma_3d, glond_plasma_3d, pz_plasma_3d, &
        iday_number)

  !g
  !g  Check here for any negative densities coming in (this shouldn't
  !g  happen but has caused me trouble before)......
  !g

      CALL ML__check_for_negative_densities( in, is, ni, vi, ne)

      iwrite98 = 0
      if (iwrite98 == 1) then

      DO 202 mp = 1 , NMP
          DO 150 lp = 1 , NLP
          !g
              if(lp == 1) then
                  write(6,3562) mp,glon(in(mp,lp),mp)/dtr, &
                  90.-(gcol(in(mp,lp),mp)/dtr)
                  3562 format('Outer flux-tube base(N): ',i4,2f7.2)
              endif
          !g
              if(mp == 1 .OR. mp == ((nmp/4)+1) .OR. mp == ((nmp/2)+1) .OR. &
              mp == ((nmp*3/4)+1)) then
                  write(98,*) mp,lp,in(mp,lp),is(mp,lp)
                  do i=in(mp,lp),is(mp,lp)
                      write(98,7877) (gr(i,mp)-r0)/1000.,90.-(gcol(i,mp)/dtr),glon(i,mp)/dtr, &
                      ni(i,mp,1)
                      7877 format(4e12.4)
                  enddo
              endif
          !g


          !g
          !g ... end of the Mp,Lp loop
          !g
          150 ENDDO
      202 ENDDO

      endif

  ! end of initial set up here....

!*******************************************************************************

  ELSE          ! .......else not initialisation, ie, normal call

!*******************************************************************************

  sw_External_model_provides_low_lat_E_fields = GIP_switches(6)


  div_vperp_3d(:,:) = 0.0

          UT_in_hours = UT_in_seconds/3600.
          do mp=1,nmp
          do lp=1,nlp
              midpoint(lp) = (in(mp,lp) + is(mp,lp)) / 2
              local_time_apex(mp,lp) = glon(midpoint(lp),mp)/dtr/15. + UT_in_hours
              if(local_time_apex(mp,lp) > 24.) local_time_apex(mp,lp) = local_time_apex(mp,lp) - 24.
          enddo
          enddo
  !g
      if(i_first_call_of_plasma == 0) then
          do mp=1,nmp

!
! ......But we need to remember that the dynamo grid (from the electrodynamics solver) goes from
!       -180 to 180 whereas GIP goes from 0 to 360... (this is the rotation bug spotted by Houjun)...
!
!       ....so we need the following 2 lines.....
!
          ilon_dynamo = mp + (nmp / 2)
          if(ilon_dynamo > nmp+1) ilon_dynamo = ilon_dynamo - nmp


          !g
          !g Loop over every second tube. These are the ones that correspond directly to the
          !g dynamo grid positions....
          !g
              do lp=1,nlp,2

                  ilat_dynamo = ((N_DYN_LAT - 1) / 2) - ((NLP-1)/2) + ((LP-1)/2)

                  potential_on_tubes(mp,lp) = potential_field(ilon_dynamo,ilat_dynamo)
                  ed1_on_tube(mp,lp) = ed1_from_dynamo(ilon_dynamo,ilat_dynamo)
                  ed2_on_tube(mp,lp) = ed2_from_dynamo(ilon_dynamo,ilat_dynamo)

              enddo
          !g
          !g For the other remaining tubes, we get our values by interpolation (averaging)...
          !g
              do lp=2,nlp-1,2

                  potential_on_tubes(mp,lp) = (potential_on_tubes(mp,lp+1)+potential_on_tubes(mp,lp-1))/2.0
                  ed1_on_tube(mp,lp) = (ed1_on_tube(mp,lp+1) + ed1_on_tube(mp,lp-1))/2.0
                  ed2_on_tube(mp,lp) = (ed2_on_tube(mp,lp+1) + ed2_on_tube(mp,lp-1))/2.0

              enddo
          enddo

      endif
  !g
  !g
!       midpoint = (NPTS+1)/2

      UT_in_hours = UT_in_seconds/3600.

      if(iwrite98 == 1) then

      WRITE (98,*) UT_in_seconds

      endif

      if (i_no_day >= i_graphics_out_start) then
          write(97,*) UT_in_hours
      endif

! Make sure we have defined our 'save' parameters here before we
! start with the ExB ionterpolation stuff......

      DO mp = 1 , NMP
      DO lp = 1 , NLP
         do i=in(mp,lp),is(mp,lp)
              ni_save(i,mp,1)=ni(i,mp,1)
              ni_save(i,mp,2)=ni(i,mp,2)
              vi_save(i,mp,1)=vi(i,mp,1)
              vi_save(i,mp,2)=vi(i,mp,2)
              ti_save(i,mp,1)=Temp_ion(i,mp,1)
              ti_save(i,mp,2)=Temp_ion(i,mp,2)
              te_save(i,mp)=Temp_electron(i,mp)
         enddo
      ENDDO
      ENDDO

!
! make sure the following params are defined for lp = 1 and nlp
!
      lp = 1
      DO mp = 1 , NMP
              do i=in(mp,lp),is(mp,lp)
                  gr_1d(i) = gr(i,mp)
                  altitude_PZ_km(i) = (GR_1D(i)-R0)/1000.
                  gcol_1d(i) = gcol(i,mp)
                  GLAt(i) = 90. - GCOl_1D(i)/DTR
                  glon_1d(i) = glon(i,mp)
                  GLOnd(i) = GLOn_1D(i)/DTR
                  div_vperp_3d(i,mp)=0.0
                  BB(i) = 0.0
                  Gravity(i) = 9.81*R0*R0/(GR_1D(i)*GR_1D(i))
              enddo
       ENDDO
      lp = nlp
      DO mp = 1 , NMP
              do i=in(mp,lp),is(mp,lp)
                  gr_1d(i) = gr(i,mp)
                  altitude_PZ_km(i) = (GR_1D(i)-R0)/1000.
                  gcol_1d(i) = gcol(i,mp)
                  GLAt(i) = 90. - GCOl_1D(i)/DTR
                  glon_1d(i) = glon(i,mp)
                  GLOnd(i) = GLOn_1D(i)/DTR
                  div_vperp_3d(i,mp)=0.0
                  BB(i) = 0.0
                  Gravity(i) = 9.81*R0*R0/(GR_1D(i)*GR_1D(i))
              enddo
       ENDDO
!
!
!       write(6,*) '***** calculating mid-latitude ExB drifts *****'



      DO 800 mp = 1 , NMP
           !write(6,*) 'mp = ',mp
      !g
      !g Only calculate an electric field and EXB drift for tubes 2 to 69.
      !g The outer and inner tubes just undergo corotation.......
      !g
          DO 750 lp = 2 , NLP-1
            !write(6,*) 'lp = ',lp
          !g
              do i=in(mp,lp),is(mp,lp)
                  gr_1d(i) = gr(i,mp)
                  altitude_PZ_km(i) = (GR_1D(i)-R0)/1000.
                  gcol_1d(i) = gcol(i,mp)
                  GLAt(i) = 90. - GCOl_1D(i)/DTR
                  glon_1d(i) = glon(i,mp)
                  GLOnd(i) = GLOn_1D(i)/DTR
                  div_vperp_3d(i,mp)=0.0
                  BB(i) = 0.0
                  Gravity(i) = 9.81*R0*R0/(GR_1D(i)*GR_1D(i))
              enddo
              if(i_first_call_of_plasma == 0) then
                  call ML__Calculate_div_vperp(Apex_e1,Apex_e2,Apex_Be3,Apex_grdlbm2, &
                  ed1_on_tube,ed2_on_tube,div_vperp_3d,in,is,mp,lp)
              endif

!cg
!cg need to know the midpoint of each tube..... 
!cg
              midpoint(lp) = (in(mp,lp) + is(mp,lp)) / 2

!             write(6,*) 'midpoint ',lp,in(mp,lp),is(mp,lp),midpoint(lp)

          !g
          !g  get an empirical value for the vertical ExB drift (for when we are not running the
          !g  full coupled electrodynamics - or at least not allowing the electrodynamic Efields to feedback)
          !g  .....
          !g
              IF ( .NOT. sw_External_model_provides_low_lat_E_fields ) THEN

              apex_height_km = altitude_PZ_km(midpoint(lp))
              geolon_apex_degrees = glond(midpoint(lp))
              local_time_apex(mp,lp) = glon(midpoint(lp),mp)/dtr/15. + UT_in_hours
              if(local_time_apex(mp,lp) > 24.) local_time_apex(mp,lp) = local_time_apex(mp,lp) - 24.
              if(apex_height_km > 1000.) then
              !g
              !g for apex_heights of greater than 1000km we use the Richmond model
              !g and therefore need to know a few things like which point on the tube is nearest
              !g 300km and such like...
              !g
                  gr300= r0 + 300.e3
                  call ML__FASTHL__NEARHTPLA(gr_1d,gr300,in(mp,lp),midpoint(lp),in1,in2,i300(mp,lp), &
                  ifailed)
                  xmlon_300 = blon(mp,lp)/dtr
                  xmlat_300 = 90. - (bcol(i300(mp,lp),mp)/dtr)
              !g
                  magnitude_e2_at_300 = sqrt((apex_e2(1,i300(mp,lp),mp)*apex_e2(1,i300(mp,lp),mp)) + &
                  (apex_e2(2,i300(mp,lp),mp)*apex_e2(2,i300(mp,lp),mp)) + &
                  (apex_e2(3,i300(mp,lp),mp)*apex_e2(3,i300(mp,lp),mp)))
              !g
                  magnitude_e2_at_apex = sqrt((apex_e2(1,midpoint(lp),mp)*apex_e2(1,midpoint(lp),mp)) + &
                  (apex_e2(2,midpoint(lp),mp)*apex_e2(2,midpoint(lp),mp)) + &
                  (apex_e2(3,midpoint(lp),mp)*apex_e2(3,midpoint(lp),mp)))

                  e2_factor_from_300km_to_apex = magnitude_e2_at_apex / magnitude_e2_at_300
              else
              !g
              !g ... for apex heights below 1000km we do not use the Richmond model and therefore these
              !g  parameters are not needed (thus set to zero).....
              !g
                  xmlon_300 = 0.0
                  xmlat_300 = 0.0
                  e2_factor_from_300km_to_apex = 0.0
              endif

              call ML__get_empirical_vertical_exb(apex_height_km,xmlat_300,xmlon_300,local_time_apex(mp,lp), &
              geolon_apex_degrees, &
              ut_in_hours,iday_number,f107,e2_factor_from_300km_to_apex, &
              v_upwards_at_apex_empirical)

              ENDIF

          !g
          !g The outer most tube and inner tubes with apex heights of less than 100km
          !g do not have any ExB convection......
          !g
                if(lp > 1 .AND. ((re_apex(mp,lp)-r0)/1000.) > 100.) then
!             if(lp > 10 .AND. ((re_apex(mp,lp)-r0)/1000.) > 100.) then
              !g
              ! vpeq(mp,lp)=vpeq_coupled(mp,lp)
              !g
              !g
              !g  The potential field is on the dynamo grid - related to our flux-tube grid
              !g  but it extends further north and south to each pole.
              !g  Convert between our 'lp' and the dynamo index 'ilat_dynamo' and define
              !g  a new parameter 'potential_on_tubes' ...
              !g
                  if(i_first_call_of_plasma == 0) then
                      Ve1 = ed2_on_tube(mp,lp) / apex_be3(in(mp,lp),mp)
                      Ve2 = 0.0 - ( ed1_on_tube(mp,lp) / apex_be3(in(mp,lp),mp) )
                  else
                      Ve1 = 0.0
                      Ve2 = 0.0
                  endif
              !g
              !g  Ve1 and Ve2 are constant along the magnetic field line.
              !g  We can now calculate the vertical component of velocity at the apex...
              !g

		 V_upwards_at_apex = (Ve1 * apex_e1(3,midpoint(lp),mp)) + (Ve2 * apex_e2(3,midpoint(lp),mp))
		 V_upwards_at_apex_array(mp,lp) = V_upwards_at_apex

                 do i = in(mp,lp) , is(mp,lp)
                   V_exb_east(i,mp) = (Ve1 * apex_e1(1,i,mp)) + (Ve2 * apex_e2(1,i,mp))
                   V_exb_north(i,mp) = (Ve1 * apex_e1(2,i,mp)) + (Ve2 * apex_e2(2,i,mp))
                   V_exb_up(i,mp) = (Ve1 * apex_e1(3,i,mp)) + (Ve2 * apex_e2(3,i,mp))
                 enddo

              !g
              !g  ...which then becomes our vpeq
              !g
              !g
                  IF ( sw_External_model_provides_low_lat_E_fields ) THEN

                      if(i_first_call_of_plasma == 1) then
                          vpeq(mp,lp)=vpeq(mp,lp)
                      else
                          vpeq(mp,lp)=V_upwards_at_apex
                      endif
!                     write(185,*) UT_in_seconds,mp,lp,V_upwards_at_apex

                  ELSE

                      vpeq(mp,lp) = v_upwards_at_apex_empirical

                  ENDIF
              !g


! cgm - start of the empirical electric fields insert

              apex_height_km = altitude_PZ_km(midpoint(lp))
              geolon_apex_degrees = glond(midpoint(lp))
              local_time_apex(mp,lp) = glon(midpoint(lp),mp)/dtr/15. + UT_in_hours
              if(local_time_apex(mp,lp) > 24.) local_time_apex(mp,lp) = local_time_apex(mp,lp) - 24.
              if(apex_height_km > 1000.) then
              !g
              !g for apex_heights of greater than 1000km we use the Richmond model
              !g and therefore need to know a few things like which point on the tube is nearest
              !g 300km and such like...
              !g
                  gr300= r0 + 300.e3
                  call ML__FASTHL__NEARHTPLA(gr_1d,gr300,in(mp,lp),midpoint(lp),in1,in2,i300(mp,lp), &
                  ifailed)
                  xmlon_300 = blon(mp,lp)/dtr
                  xmlat_300 = 90. - (bcol(i300(mp,lp),mp)/dtr)
              !g
                  magnitude_e2_at_300 = sqrt((apex_e2(1,i300(mp,lp),mp)*apex_e2(1,i300(mp,lp),mp)) + &
                  (apex_e2(2,i300(mp,lp),mp)*apex_e2(2,i300(mp,lp),mp)) + &
                  (apex_e2(3,i300(mp,lp),mp)*apex_e2(3,i300(mp,lp),mp)))
              !g
                  magnitude_e2_at_apex = sqrt((apex_e2(1,midpoint(lp),mp)*apex_e2(1,midpoint(lp),mp)) + &
                  (apex_e2(2,midpoint(lp),mp)*apex_e2(2,midpoint(lp),mp)) + &
                  (apex_e2(3,midpoint(lp),mp)*apex_e2(3,midpoint(lp),mp)))

                  e2_factor_from_300km_to_apex = magnitude_e2_at_apex / magnitude_e2_at_300
              else
              !g
              !g ... for apex heights below 1000km we do not use the Richmond model and therefore these
              !g  parameters are not needed (thus set to zero).....
              !g
                  xmlon_300 = 0.0
                  xmlat_300 = 0.0
                  e2_factor_from_300km_to_apex = 0.0
              endif

              call ML__get_empirical_vertical_exb(apex_height_km,xmlat_300,xmlon_300,local_time_apex(mp,lp), &
              geolon_apex_degrees, &
              ut_in_hours,iday_number,f107,e2_factor_from_300km_to_apex, &
              v_upwards_at_apex_empirical)




! cgm - end of the empirical electric fields insert




                  local_time_apex(mp,lp) = glon(midpoint(lp),mp)/dtr/15. + UT_in_hours
                  if(local_time_apex(mp,lp) > 24.) local_time_apex(mp,lp) = local_time_apex(mp,lp) - 24.
              !g
!                write(121,9611) ut_in_hours,mp,lp,altitude_pz_km(midpoint(lp)),local_time_apex(mp,lp),Vpeq(mp,lp), &
!                                ed1_on_tube(mp,lp),ed2_on_tube(mp,lp),apex_be3(in(mp,lp),mp),apex_e1(3,midpoint(lp),mp), &
!                                apex_e2(3,midpoint(lp),mp)

!      write(121,9611) ut_in_hours,mp,lp,altitude_pz_km(midpoint(lp)),local_time_apex(mp,lp),V_exb_east(midpoint(lp),mp), &
!            V_exb_north(midpoint(lp),mp),V_exb_up(midpoint(lp),mp),v_upwards_at_apex_empirical

!                9611 format(f7.2,2i4,f12.2,f7.2,f9.2,5e12.4)
                 9611 format(f7.2,2i4,f12.2,f7.2,4f9.2)

!                if (mp.eq.1.or.mp.eq.21.or.mp.eq.41.or.mp.eq.61) then
!                if (lp.eq.7.or.lp.eq.18.or.lp.eq.28.or.lp.eq.38) then
!                write(122,*) ut_in_hours,mp,lp,in(mp,lp),is(mp,lp),local_time_apex(mp,lp)
!                do i = in(mp,lp) , is(mp,lp)
!                   write(122,9612) i,altitude_pz_km(i),V_exb_east(i,mp),V_exb_north(i,mp),V_exb_up(i,mp)
!9612 format(i4,f12.2,3f9.2)
!                enddo
!                endif
!                endif
              !g

                  CALL ML__TUBES_LINEAR_INTERPOLATE(re_apex,vpeq,plasma_time_step_seconds,ni, &
                  q_coordinate,vi,Temp_ion,Temp_electron,ne,ni_save,vi_save,ti_save,te_save &
                  ,mp,lp,in,is)

              !g
              !g else we have a tube which hasn't drifted....
              !g
              else
              !g
                  do i=in(mp,lp),is(mp,lp)
                      ni_save(i,mp,1)=ni(i,mp,1)
                      ni_save(i,mp,2)=ni(i,mp,2)
                      vi_save(i,mp,1)=vi(i,mp,1)
                      vi_save(i,mp,2)=vi(i,mp,2)
                      ti_save(i,mp,1)=Temp_ion(i,mp,1)
                      ti_save(i,mp,2)=Temp_ion(i,mp,2)
                      te_save(i,mp)=Temp_electron(i,mp)
                  enddo
                  vpeq(mp,lp)=0.0
              !g
              endif
          !g
          750 ENDDO
      800 ENDDO
!       write(6,*) '***** done mid-latitude ExB drifts *****'
     do mp=1,NMP
! BUG this should be looped over all the tubes - fixed on the next line
!        do lp=2,NLP-1
! BUG
          do lp=1,NLP
              do i=in(mp,lp),is(mp,lp)
                  ni(i,mp,1)=ni_save(i,mp,1)
                  ni(i,mp,2)=ni_save(i,mp,2)
                  ne(i,mp)=ni(i,mp,1)+ni(i,mp,2)
                  vi(i,mp,1)=vi_save(i,mp,1)
                  vi(i,mp,2)=vi_save(i,mp,2)
                  Temp_ion(i,mp,1)=ti_save(i,mp,1)
                  Temp_ion(i,mp,2)=ti_save(i,mp,2)
                  Temp_electron(i,mp)=te_save(i,mp)
              enddo
          enddo
      enddo
  !g
  !g Now we know the previous position of each tube and have calculated all
  !g of the densities and temperatures and stuff at this previous position
  !g we can get on with calculating all these parameters for the present
  !g time step.  We loop over all the tubes again....
  !g
!           write(6,*) '******* main plasma calculation ********'


!   write(6,*) '***************** HYDROGEN GRID AT 400KM **INSIDE PLASMA*********************'
!   write(6,1888) Hyd_grid_400km_m3 / 1.e11
!   1888 format(25f4.1)
!   write(6,*) '****************************************************************************'

      DO 500 mp = 1 , NMP

          vpeq(mp,nlp)=0.0
          vpeq(mp,1)=0.0
!           write(6,*) 'mp = ',mp

          DO 450 lp = 1 , NLP
!             write(6,*) 'lp ',lp

              do i=in(mp,lp),is(mp,lp)
                  Apex_D1_1d(1,i) = Apex_D1(1,i,mp)
                  Apex_D1_1d(2,i) = Apex_D1(2,i,mp)
                  Apex_D1_1d(3,i) = Apex_D1(3,i,mp)
                  Apex_D2_1d(1,i) = Apex_D2(1,i,mp)
                  Apex_D2_1d(2,i) = Apex_D2(2,i,mp)
                  Apex_D2_1d(3,i) = Apex_D2(3,i,mp)
                  Apex_BE3_1d(i) = Apex_BE3(i,mp)
                  Apex_D_1d(i) = Apex_D(i,mp)
                  Apex_d1d1_1d(i) = Apex_d1d1(i,mp)
                  Apex_d1d2_1d(i) = Apex_d1d2(i,mp)
                  Apex_d2d2_1d(i) = Apex_d2d2(i,mp)
                  Apex_BMAG_1d(i) = Apex_BMAG(i,mp)
                  apex_bhat_1d(1,i) =  Apex_bhat(1,i,mp)
                  apex_bhat_1d(2,i) =  Apex_bhat(2,i,mp)
                  apex_bhat_1d(3,i) =  Apex_bhat(3,i,mp)
                  gr_1d(i) = gr(i,mp)
                  altitude_PZ_km(i) = (GR_1D(i)-R0)/1000.
                  gcol_1d(i) = gcol(i,mp)
                  GLAt(i) = 90. - GCOl_1D(i)/DTR
                  glon_1d(i) = glon(i,mp)
                  GLOnd(i) = GLOn_1D(i)/DTR
                  tn_1d(i) = Tn_plasma_input_3d(i,mp)
                  o_1d(i) = O_plasma_input_3d(i,mp)
                  o2_1d(i) = O2_plasma_input_3d(i,mp)
                  n2_1d(i) = N2_plasma_input_3d(i,mp)
!                 NO_1d(i) = NO_plasma_input_3d(i,mp)
!                 N4S_1d(i) = N4S_plasma_input_3d(i,mp)
!                 N2D_1d(i) = N2D_plasma_input_3d(i,mp)
                  U_merid_1d(i) = Um_plasma_input_3d(i,mp)
                  U_zonal_1d(i) = Uz_plasma_input_3d(i,mp)
                  U_vertical_1d(i) = Uv_plasma_input_3d(i,mp)
!                 Te_1d(i) = Te_dum_plasma_input_3d(i,mp)
                  eta_apex_1d(i) = eta_apex_3d(i,mp)
                  no_plus_1d(i) = no_plus_3d(i,mp)
                  o2_plus_1d(i) = o2_plus_3d(i,mp)
                  n2_plus_1d(i) = n2_plus_3d(i,mp) 
                  n_plus_1d(i) = n_plus_3d(i,mp)
                  div_vperp_1d(i) = div_vperp_3d(i,mp)
                  Gravity(i) = 9.81*R0*R0/(GR_1D(i)*GR_1D(i))
                  q_coordinate_1d(i) = q_coordinate(i,mp)
                  ni_oplus_1d(i) = ni(i,mp,1)
                  ni_hplus_1d(i) = ni(i,mp,2)
                  vi_oplus_1d(i) = vi(i,mp,1)
                  vi_hplus_1d(i) = vi(i,mp,2)
                  ne_1d(i) = ne(i,mp)

!    ISSUE:      xtimes isn't defined.....

!                 if(xtimes(i) .lt. 12.0 ) then
!                    morning = .true.
!                 else
!                    morning = .false.
!                 endif

                  glat_rad = glat(i)*dtr
                  geo_local_time_degrees = glond(i) + (UT_in_seconds - 43200.)/240.0
                  rlt = 180.0 + geo_local_time_degrees
                  IF ( rlt.GT.360.0 ) rlt = rlt - 360.0
                  rlt = rlt*DTR
                  solar_zenith_angle_degrees =  &
                            (ACOS(-COS(glat_rad)*COS(solar_declination_angle_degrees*DTR)*COS(rlt)+SIN(glat_rad) &
                           *SIN(solar_declination_angle_degrees*DTR))/dtr)

!                 call ML__ETEMP_day_night(solar_zenith_angle_degrees,morning,altitude_PZ_km(i),te_tn_factor)
!                 te_1d(i) = te_tn_factor * tn_1d(i)
                  te_1d(i) = 2.0 * tn_1d(i)

              enddo





        do ilt = 1 , 25
         if (hyd_grid_local_times(ilt) .gt. local_time_apex(mp,lp)) then
           ilt_east = ilt
           ilt_west = ilt - 1
          goto 1822
         endif
        enddo

1822    factor_lt = local_time_apex(mp,lp) - hyd_grid_local_times(ilt_west)

!        write(6,*) 'factor_lt ',mp,lp,local_time_apex(mp,lp) , hyd_grid_local_times(ilt_west),  factor_lt

                do i=in(mp,lp),is(mp,lp)

                  do ilat = 1 , 19
                    if (hyd_grid_lats(ilat) .gt. glat_plasma_3d(i,mp)) then
                    ilat_north = ilat
                    ilat_south = ilat - 1
                    goto 1823
                    endif
                  enddo
1823              factor_lat = ( glat_plasma_3d(i,mp) - hyd_grid_lats(ilat_south) ) /10.
!                  write(6,*) 'factor lat ',m , ilat_north , ilat_south , factor_lat
                  H_NE = Hyd_grid_400km_m3(ilt_east,ilat_north)
                  H_NW = Hyd_grid_400km_m3(ilt_west,ilat_north)
                  H_SE = Hyd_grid_400km_m3(ilt_east,ilat_south)
                  H_SW = Hyd_grid_400km_m3(ilt_west,ilat_south)

                  H_N = ((H_NE - H_NW) * factor_lt) + H_NW
                  H_S = ((H_SE - H_SW) * factor_lt) + H_SW

                  H_400 = ((H_N - H_S) * factor_lat) + H_S
!              if (lp.eq.1) then
!              write(169,2341) i,mp,lp,local_time_apex(mp,lp),glat_plasma_3d(i,mp),H_400
!               write(169,2341) i,local_time_apex(mp,lp),glon(midpoint(lp),mp),UT_in_hours
!2341 format(3i6,2f10.1,x3e12.4)
!              endif

                    CALL HL_ML__HYDEQ_BASIC(H_400,altitude_pz_km(i),tn_1d(i),gravity(i), &
                                   neutral_hydrogen_1d(i))
                enddo


              midpoint(lp) = (in(mp,lp) + is(mp,lp)) / 2

          i_write_out_tube = 0
!         if(lp.eq.1.and.mp.eq.56) then
!         i_write_out_tube = 1
!         write(6,*) 'calling single_flux_tube ',mp,lp,midpoint(lp)
!         endif
          call ML__single_flux_tube_1D_calculation( &
                                               GIP_switches, &
                                               mp,lp,in(mp,lp),is(mp,lp), &
                                               midpoint(lp), &
                                               glon_1d, &
                                               Apex_D1_1d, &
                                               Apex_D2_1d, &
                                               Apex_BE3_1d, &
                                               Apex_D_1d, &
                                               Apex_d1d1_1d, &
                                               Apex_d1d2_1d, &
                                               Apex_d2d2_1d, &
                                               Apex_BMAG_1d, &
                                               apex_bhat_1d, &
                                               gr_1d, &
                                               gcol_1d, &
                                               tn_1d, &
                                               o_1d, &
                                               o2_1d, &
                                               n2_1d, &
                                               neutral_hydrogen_1d, &
                                               NO_1d, &
                                               N4S_1d, &
                                               N2D_1d, &
                                               U_merid_1d, &
                                               U_zonal_1d, &
                                               U_vertical_1d, &
                                               gravity, &
                                               eta_apex_1d, &
                                               div_vperp_1d, &
                                               q_coordinate_1d, &
                                               no_plus_1d, &
                                               o2_plus_1d, &
                                               n2_plus_1d, & 
                                               n_plus_1d, &
                                               ni_oplus_1d, &
                                               ni_hplus_1d, &
                                               vi_oplus_1d, &
                                               vi_hplus_1d, &
                                               ne_1d, &
                                               ti_oplus_1d, &
                                               ti_hplus_1d, &
                                               te_1d, &
                                               sigma_phph_dsi_1d, &
                                               sigma_lmlm_msi_1d, &
                                               sigma_h_1d, &
                                               sigma_c_1d, &
                                               Kdmph_dsi_1d, &
                                               Kdmlm_1d, &
                                               plasma_time_step_seconds, &
                                               f107, &
                                               solar_declination_angle_degrees, &
                                               ut_in_seconds, &
                                               i_write_out_tube, &
                                               i_attempt, &
                                               iday_number)

          i_attempt_array(mp,lp) = i_attempt


          do i=in(mp,lp),is(mp,lp)
                no_plus_3d(i,mp) = no_plus_1d(i)
                o2_plus_3d(i,mp) = o2_plus_1d(i)
                n2_plus_3d(i,mp) = n2_plus_1d(i)
                n_plus_3d(i,mp) = n_plus_1d(i)
                NI(i,MP,1) =  ni_oplus_1d(i)
                NI(i,MP,2) =  ni_hplus_1d(i)
                Ne(i,MP) = ne_1d(i)
                VI(i,MP,1) = vi_oplus_1d(i)
                VI(i,MP,2) = vi_hplus_1d(i)
                Temp_ion(i,MP,1) = ti_oplus_1d(i)
                Temp_ion(i,MP,2) = ti_hplus_1d(i)
                Temp_electron(i,MP) = te_1d(i)
          enddo
          do i=1,2
                sigma_phph_dsi(i,mp,lp) = sigma_phph_dsi_1d(i)
                sigma_lmlm_msi(i,mp,lp) = sigma_lmlm_msi_1d(i)
                sigma_h(i,mp,lp) = sigma_h_1d(i)
                sigma_c(i,mp,lp) = sigma_c_1d(i)
                Kdmph_dsi(i,mp,lp) = Kdmph_dsi_1d(i)
                Kdmlm(i,mp,lp) = Kdmlm_1d(i)
          enddo

          450 ENDDO
      500 ENDDO

!        write(196,*) UT_in_hours
!        write(196,*) i_attempt_array 
!        write(196,*) local_time_apex 



!           write(6,*) '******* finished plasma calculation ********'
  !g
      i_call_the_electrodynamics = 1
      if(i_call_the_electrodynamics == 1) then
      !g
      !g
      !g  with our new dynamo code we now convert our ML__FIELD_LINE_INTEGRALS onto the (related) grid
      !g  used by the dynamo solver code.....
      !g
          !g
          !g Need to calculate the magnetic latitude of the tube end point in northern
          !g and southern hemispheres......
          !g
          do mp = 1 , NMP
          do lp = 1 , NLP
              mag_lat_tubeN(mp,lp)=90.-(bcol(in(mp,lp),mp)/dtr)
              mag_lat_tubeS(mp,lp)=90.-(bcol(is(mp,lp),mp)/dtr)
           enddo
           enddo

          call ML__convert_integral_to_dynamo_grid(sigma_phph_dsi , sigma_lmlm_msi , sigma_h , &
          sigma_c , Kdmph_dsi , Kdmlm , mag_lat_tubeN, mag_lat_tubeS, &
          dynamo_sigma_phph_dsi , dynamo_sigma_lmlm_msi , dynamo_sigma_h , &
          dynamo_sigma_c , dynamo_Kdmph_dsi , dynamo_Kdmlm , &
          dynamo_latitude)
      !g
      endif
  !g
  !g  Write the electrodynamic velocities to unit 96...
  !g
!     if (i_no_day >= i_graphics_out_start) then
!         write(96,*) UT_in_hours
!     !g
!     !g  writing out tube 36 ExB velocities here (an altitude of 310 km) ....
!     !g
!     ! write(96,3955) (vpeq(mp,20),vzon_coupled(mp,20),mp=1,nmp)
!         write(96,3955) (v_upwards_at_apex_array(mp,36),mp=1,nmp)
!         write(96,3955) (local_time_apex(mp,36),mp=1,nmp)
!     endif
      3955 format(10f8.2)


  endif !  nnloop eq zero endif

  RETURN




end SUBROUTINE ML__MID_AND_LOW_LATITUDE_IONOSPHERE















SUBROUTINE ML__ETEMP_day_night(solar_zenith_angle_degrees,morning,altitude_km,factor)

    implicit none
    real(kind=8) :: &
    rad80, rad90, rad100, tefhd, tefld, tefhn, tefln, &
    dted, dten, sza, szaday, szanight, fh, fl, &
    dte, fract, altitude_km, factor, &
    solar_zenith_angle_degrees

    logical :: morning

    rad80 = 80.
    rad90 = 90.
    rad100 = 100.
    tefhd = 4.2
    tefld = 2.4
    tefhn = 3.0
    tefln = 1.35
    dted = 1.0
    dten = 0.75

    sza = solar_zenith_angle_degrees

    IF ( SZA.LE.rad80 ) THEN
!
!     daytime conditions
!
       fh = tefhd
       fl = tefld
       dte = dted

    ENDIF

!     IF ( (SZA.gt.rad90.and.morning) .or. (SZA.gt.rad100.and.(.NOT.morning)) ) THEN
    IF ( (SZA.gt.rad100.and.morning) .or. (SZA.gt.rad100.and.(.NOT.morning)) ) THEN
!
!   nighttime conditions
!
       fh = tefhn
       fl = tefln
       dte = dten

    ENDIF

!     IF (SZA.le.rad90.and.SZA.gt.rad80.and.morning) THEN
    IF (SZA.le.rad100.and.SZA.gt.rad80.and.morning) THEN
!
!   sunrise conditions for 80<sza<90.
!
       szaday = rad80
       szanight = rad90
       fract = (SZA-szaday)/(szanight-szaday)
       fh = tefhd + fract*(tefhn-tefhd)
       fl = tefld + fract*(tefln-tefld)
       dte = dted + fract*(dten-dted)

    ENDIF

    IF (SZA.le.rad100.and.SZA.gt.rad80.and.(.NOT.morning)) THEN
!
!   sunset conditions for 80<sza<100.
!
       szaday = rad80
       szanight = rad100
       fract = (SZA-szaday)/(szanight-szaday)
       fh = tefhd + fract*(tefhn-tefhd)
       fl = tefld + fract*(tefln-tefld)
       dte = dted + fract*(dten-dted)

    ENDIF

!
!  now we have fh and fl - lets get f.....
!
    factor = (((altitude_km - 300.)/(1000. - 300.))*(fh - fl)) + fl

    if ( factor > fh ) factor = fh
    if ( factor < 1.0 ) factor = 1.0



    return

end SUBROUTINE ML__ETEMP_day_night









SUBROUTINE ML__single_flux_tube_1D_calculation( &
                                               GIP_switches, &
                                               mp,lp,in,is, &
                                               midpoint, &
                                               glon_1d, &
                                               Apex_D1_1d, &
                                               Apex_D2_1d, &
                                               Apex_BE3_1d, &
                                               Apex_D_1d, &
                                               Apex_d1d1_1d, &
                                               Apex_d1d2_1d, &
                                               Apex_d2d2_1d, &
                                               Apex_BMAG_1d, &
                                               apex_bhat_1d, &
                                               gr_1d, &
                                               gcol_1d, &
                                               tn, &
                                               o, &
                                               o2, &
                                               n2, &
                                               hyd, &
                                               NO_1d, &
                                               N4S_1d, &
                                               N2D_1d, &
                                               U_merid, &
                                               U_zonal, &
                                               U_vertical, &
                                               gravity, &
                                               eta_apex_1d, &
                                               div_vperp_1d, &
                                               q_coordinate_1d, &
                                               no_plus_1d, &
                                               o2_plus_1d, &
                                               n2_plus_1d, & 
                                               n_plus_1d, &
                                               ni_oplus_1d, &
                                               ni_hplus_1d, &
                                               vi_oplus_1d, &
                                               vi_hplus_1d, &
                                               ne_1d, &
                                               ti_oplus_1d, &
                                               ti_hplus_1d, &
                                               te_1d, &
                                               sigma_phph_dsi, &
                                               sigma_lmlm_msi, &
                                               sigma_h, &
                                               sigma_c, &
                                               Kdmph_dsi, &
                                               Kdmlm, &
                                               plasma_time_step_seconds, &
                                               f107, &
                                               solar_declination_angle_degrees, &
                                               ut_in_seconds, &
                                               i_write_out_tube, &
                                               i_attempt, &
                                               iday_number)

           IMPLICIT NONE
           INTEGER :: NPTS
           PARAMETER (NPTS = 13813)
           REAL(kind=8) :: atomic_mass_unit
           PARAMETER (atomic_mass_unit=1.66e-27)
           REAL(kind=8) :: PI , DTR
           PARAMETER (PI=3.141592654,DTR=PI/180.0)

           INTEGER :: iday_number
           INTEGER :: mp , lp
           INTEGER :: IN , IS , i , in1 , in2
           INTEGER :: midpoint
           INTEGER :: Ofailed , Hfailed
           INTEGER :: ite_failed , iti_failed
           INTEGER :: i_failed_molecular_ions
           INTEGER :: idiagnose
           INTEGER :: istop
           LOGICAL :: sw_use_EUVAC_solar_spectrum
           LOGICAL :: GIP_switches(20)

           REAL(kind=8) :: plasma_time_step_seconds
           REAL(kind=8) :: plasma_energy_time_step_seconds
           REAL(kind=8) :: ut_in_seconds
           REAL(kind=8) :: solar_declination_angle_degrees
           REAL(kind=8) :: f107
           REAL(kind=8) :: gr_1d(NPTS)
           REAL(kind=8) :: gcol_1d(NPTS)
           REAL(kind=8) :: glon_1d(NPTS)
           REAL(kind=8) :: tn(NPTS)
           REAL(kind=8) :: o(NPTS)
           REAL(kind=8) :: o2(NPTS)
           REAL(kind=8) :: n2(NPTS)
           REAL(kind=8) :: hyd(NPTS)
           REAL(kind=8) :: hel(NPTS)
           REAL(kind=8) :: NO_1d(NPTS)
           REAL(kind=8) :: N4S_1d(NPTS)
           REAL(kind=8) :: N2D_1d(NPTS)
           REAL(kind=8) :: U_merid(NPTS)
           REAL(kind=8) :: U_zonal(NPTS)
           REAL(kind=8) :: U_vertical(NPTS)
           REAL(kind=8) :: U_merid_for_integrals(NPTS)
           REAL(kind=8) :: U_zonal_for_integrals(NPTS)
           REAL(kind=8) :: U_vertical_for_integrals(NPTS)
           REAL(kind=8) :: upar_zonal_contribution(NPTS)
           REAL(kind=8) :: upar_meridional_contribution(NPTS)
           REAL(kind=8) :: upar_vertical_contribution(NPTS)
           REAL(kind=8) :: g_parallel(NPTS)
           REAL(kind=8) :: eta_apex_1d(NPTS)
           REAL(kind=8) :: apex_bhat_1d(3,NPTS)
           REAL(kind=8) :: no_plus_1d(NPTS)
           REAL(kind=8) :: o2_plus_1d(NPTS)
           REAL(kind=8) :: n2_plus_1d(NPTS)
           REAL(kind=8) :: n_plus_1d(NPTS)
           REAL(kind=8) :: DQ_1d(NPTS)
           REAL(kind=8) :: distance_ds(NPTS)
           REAL(kind=8) :: altitude_PZ_km(NPTS)
           REAL(kind=8) :: GLAt(NPTS)
           REAL(kind=8) :: GLOnd(NPTS)
           REAL(kind=8) :: div_vperp_1d(NPTS)
           REAL(kind=8) :: BB(NPTS)
           REAL(kind=8) :: Gravity(NPTS)
           REAL(kind=8) :: time(NPTS)
           REAL(kind=8) :: xtimes(NPTS)
           REAL(kind=8) :: angle(NPTS)
           REAL(kind=8) :: w1 , w2 , wo , wo2 , wn2
           REAL(kind=8) :: Burnside_factor
           REAL(kind=8) :: tt(NPTS)
           REAL(kind=8) :: nuin(NPTS)
           REAL(kind=8) :: vp(NPTS)
           REAL(kind=8) :: u2dif(NPTS)
           REAL(kind=8) :: uperp_apex(NPTS)
           REAL(kind=8) :: upar_apex(NPTS)
           REAL(kind=8) :: vi_limited
           REAL(kind=8) :: nu_mol_n(NPTS)
           REAL(kind=8) :: frictional_heating_1(NPTS)
           REAL(kind=8) :: frictional_heating_2(NPTS)
           REAL(kind=8) :: total_moleculars_1d(NPTS)
           REAL(kind=8) :: nuin_tot(NPTS)
           REAL(kind=8) :: beta_1(NPTS)
           REAL(kind=8) :: beta_2(NPTS)
           REAL(kind=8) :: chemp_1(NPTS)
           REAL(kind=8) :: chemp_2(NPTS)
           REAL(kind=8) :: ion_mass(NPTS)
           REAL(kind=8) :: ni_Oplus_1d(NPTS)
           REAL(kind=8) :: ni_Hplus_1d(NPTS)
           REAL(kind=8) :: dni_Oplus_1d(NPTS)
           REAL(kind=8) :: dni_Hplus_1d(NPTS)
           REAL(kind=8) :: vi_Oplus_1d(NPTS)
           REAL(kind=8) :: vi_Hplus_1d(NPTS)
           REAL(kind=8) :: ti_Oplus_1d(NPTS)
           REAL(kind=8) :: ti_Hplus_1d(NPTS)
           REAL(kind=8) :: dti_Oplus_1d(NPTS)
           REAL(kind=8) :: dti_Hplus_1d(NPTS)
           REAL(kind=8) :: ne_1d(NPTS)
           REAL(kind=8) :: te_1d(NPTS)
           REAL(kind=8) :: dte_1d(NPTS)
           REAL(kind=8) :: ni_Oplus_1d_saved(NPTS)
           REAL(kind=8) :: ni_Hplus_1d_saved(NPTS)
           REAL(kind=8) :: dni_Oplus_1d_saved(NPTS)
           REAL(kind=8) :: dni_Hplus_1d_saved(NPTS)
           REAL(kind=8) :: vi_Oplus_1d_saved(NPTS)
           REAL(kind=8) :: vi_Hplus_1d_saved(NPTS)
           REAL(kind=8) :: ne_1d_saved(NPTS)
           REAL(kind=8) :: no_plus_1d_saved(NPTS)
           REAL(kind=8) :: o2_plus_1d_saved(NPTS)
           REAL(kind=8) :: n2_plus_1d_saved(NPTS)
           REAL(kind=8) :: n_plus_1d_saved(NPTS)
           REAL(kind=8) :: ti_Oplus_1d_saved(NPTS)
           REAL(kind=8) :: ti_Hplus_1d_saved(NPTS)
           REAL(kind=8) :: te_1d_saved(NPTS)
           REAL(kind=8) :: euv_production_rate_Oplus(NPTS)
           REAL(kind=8) :: euv_production_rate_N2plus(NPTS)
           REAL(kind=8) :: euv_production_rate_O2plus(NPTS)
           REAL(kind=8) :: auroral_production_rate_Oplus(NPTS)
           REAL(kind=8) :: auroral_production_rate_N2plus(NPTS)
           REAL(kind=8) :: auroral_production_rate_O2plus(NPTS)
           REAL(kind=8) :: qin(NPTS)
           REAL(kind=8) :: cos_SZA(NPTS)
           REAL(kind=8) :: chi_degrees(NPTS)
           REAL(kind=8) :: sigma_phph_dsi(2)
           REAL(kind=8) :: sigma_lmlm_msi(2)
           REAL(kind=8) :: sigma_h(2)
           REAL(kind=8) :: sigma_c(2)
           REAL(kind=8) :: Kdmph_dsi(2)
           REAL(kind=8) :: Kdmlm(2)
           REAL(kind=8) :: peuvn(NPTS,6)
           REAL(kind=8) :: peuvi(NPTS,6)
           REAL(kind=8) :: Apex_D1_1d(3,NPTS)
           REAL(kind=8) :: Apex_D2_1d(3,NPTS)
           REAL(kind=8) :: Apex_BE3_1d(NPTS)
           REAL(kind=8) :: Apex_D_1d(npts)
           REAL(kind=8) :: Apex_d1d1_1d(npts)
           REAL(kind=8) :: Apex_d1d2_1d(npts)
           REAL(kind=8) :: Apex_d2d2_1d(npts)
           REAL(kind=8) :: Apex_BMAG_1d(npts)
           REAL(kind=8) :: q_coordinate_1d(npts)
           REAL(kind=8) :: height_factor
           REAL(kind=8) :: sigma_ped_1d(npts)
           REAL(kind=8) :: sigma_hall_1d(npts)
           REAL(kind=8) :: tiegcm_sigma_ped_1d(npts)
           REAL(kind=8) :: tiegcm_sigma_hall_1d(npts)
           REAL(kind=8) :: O_plus_production_fudge_factor

           integer i_write_out_tube
           integer i_use_tiegcm_ions_for_dynamo_calculation
           integer i_attempt

           integer IN_double
           integer IS_double
           REAL(kind=8) :: nuin_double(1000)
           REAL(kind=8) :: chemp_1_double(1000)
           REAL(kind=8) :: chemp_2_double(1000)   
           REAL(kind=8) :: beta_1_double(1000)
           REAL(kind=8) :: beta_2_double(1000)
           REAL(kind=8) :: peuvi_double(1000)
           REAL(kind=8) :: g_parallel_double(1000)
           REAL(kind=8) :: dte_1d_double(1000)
           REAL(kind=8) :: dti_Oplus_1d_double(1000)
           REAL(kind=8) :: dti_Hplus_1d_double(1000)
           REAL(kind=8) :: upar_apex_double(1000)
           REAL(kind=8) :: eta_apex_1d_double(1000)
           REAL(kind=8) :: dq_1d_double(1000)
           REAL(kind=8) :: TI_Oplus_1d_double(1000)
           REAL(kind=8) :: TI_Hplus_1d_double(1000)
           REAL(kind=8) :: TE_1d_double(1000)
           REAL(kind=8) :: ni_OPlus_1d_double(1000)
           REAL(kind=8) :: ni_Hplus_1d_double(1000)
           REAL(kind=8) :: dni_oplus_1d_double(1000)
           REAL(kind=8) :: dni_hplus_1d_double(1000)
           REAL(kind=8) :: Vi_oplus_1d_double(1000)
           REAL(kind=8) :: vi_hplus_1d_double(1000)
           REAL(kind=8) :: NE_1d_double(1000)
           REAL(kind=8) :: div_vperp_1d_double(1000)
           REAL(kind=8) :: O_double(1000)
           REAL(kind=8) :: eccentric ! eccentricity of earth's orbit


           sw_use_EUVAC_solar_spectrum = .FALSE.  ! ISSUE : this isn't begin passed through correctly
                                                  ! at present

           idiagnose = 0

              do i=in,is
                  altitude_PZ_km(i) = (GR_1D(i)-R0)/1000.
                  GLAt(i) = 90. - GCOl_1D(i)/DTR
                  GLOnd(i) = GLOn_1D(i)/DTR
                  BB(i) = 0.0
                  g_parallel(i) = apex_bhat_1d(3,i) * 9.81*r0*r0/(gr_1d(i)*gr_1d(i))
              end do
              DO i = IN , IS - 1
                  DQ_1d(i) = q_coordinate_1d(i+1) - q_coordinate_1d(i)
                  distance_ds(i) = dq_1d(i)/eta_apex_1d(i)
              ENDDO
                  distance_ds(IS) = distance_ds(IS-1)    ! Make sure we have a definition for the Southern end point IS
          !g
          !g  calculates LT for all points on tube (in seconds for some reason).....
          !g
              DO 20 i = IN , IS
                  time(i) = UT_in_seconds + glond(i)*240.
                  IF ( time(i) > 86400. ) THEN
                      time(i) = time(i) - 86400.
                  ELSEIF ( time(i) < 0. ) THEN
                      time(i) = time(i) + 86400.
                  ENDIF
                  xtimes(i) = time(i) / 3600.
              ! m
              ! m  ANGLE is the hour angle (degrees) corresponding to the LT above..
              ! m
                  angle(i) = UT_in_seconds/240. + glon_1d(i)/DTR
                  IF ( angle(i) >= 360. ) angle(i) = angle(i) - 360.
              20 ENDDO

          if (idiagnose.eq.1) write(6,*) 'here 1'

          !g calculate the neutral background by interpolating
          !g from the neutral code to give Tn, winds and composition
          !g
          !g Calculate upar_apex, uperp_apex, using the apex vectors.....
          !g
              CALL ML__get_U_apex(in,is,Apex_bhat_1d,U_zonal,U_merid,U_vertical, &
              upar_apex,uperp_apex, &
              upar_zonal_contribution,upar_meridional_contribution,upar_vertical_contribution)

          if (idiagnose.eq.1) write(6,*) 'here 2'
          !g
                  do i = in , is
                      vi_limited = vi_Oplus_1d(i)
                      if(vi_limited > 1000.) vi_limited = 1000.
                      if(vi_limited < -1000.) vi_limited = -1000.
                  !g
!g
!g this vp needs investigating.... it's a bug captain....
!g
vp(i) = 0.0


                      u2dif(i) = (vp(i)-uperp_apex(i))*(vp(i)-uperp_apex(i)) &
                      + ((vi_limited - upar_apex(i)) &
                      * (vi_limited - upar_apex(i)))
                  enddo
              !g

          if (idiagnose.eq.1) write(6,*) 'here 3'
 
                  CALL ML__ION_TEMP_PLASMA(TN,TE_1d,ni_Oplus_1d,ni_Hplus_1d,TT, &
                                       N2,O,O2,U2Dif,in,is)
              !g
                      do i = in, is
                          ti_Oplus_1d(i) = tt(i)
                          ti_Hplus_1d(i) = tt(i)
                      enddo

          if (idiagnose.eq.1) write(6,*) 'here 4'


          ! Calculate df/ds for n(h+),n(o+),t(h+),t(o+) and te

              CALL ML__DF_BY_DS(te_1d,dte_1d,IN,IS,eta_apex_1d,dq_1d)

              CALL ML__DF_BY_DS(TI_Oplus_1d,dti_Oplus_1d,IN,IS,eta_apex_1d,dq_1d)

              CALL ML__DF_BY_DS(TI_Hplus_1d,dti_Hplus_1d,IN,IS,eta_apex_1d,dq_1d)

              CALL ML__DF_BY_DS(NI_Oplus_1d,dni_Oplus_1d,IN,IS,eta_apex_1d,dq_1d)

              CALL ML__DF_BY_DS(NI_Hplus_1d,dni_Hplus_1d,IN,IS,eta_apex_1d,dq_1d)

          if (idiagnose.eq.1) write(6,*) 'here 5'
          !g
          !g  calculate parameters required for
          !g  high latitude production function......
          !g
          !g
          !g ....and calculate the solar production rate....
          !g

 if(sw_use_EUVAC_solar_spectrum) then
! write(6,*) ' plasma switch 1 ',sw_use_EUVAC_solar_spectrum
!*****************************
! ALD March 09: Put HL_ML__EUV_ION_PRODUCTION2 call in here:
! This routine will use TIEGCM EUVAC solar spectrum to calculate major
! species photoionisation rates.

      eccentric=(1.0+0.0167*COS(2.0*pi*(iday_number-3)/366.0))**2

      CALL HL_ML__EUV_ION_PRODUCTION_2(npts, in, is, o, o2, n2,     &
                  Solar_Declination_Angle_degrees, time,     &
                  f107, gcol_1d, R0, tn, gravity, km_plasma, &
                  gr_1d, eccentric, chi_degrees, peuvn)

!     if (mp.eq.1.and.lp.eq.5) then
!       write(196,*) 'new euv_ion ', UT_in_seconds/3600., midpoint, &
!                     xtimes(midpoint)
!       do i=in,is
!         write(196,2129) i,(gr_1d(i)-r0)/1000.,peuvn(i,1)
!       enddo
!     endif

  else
! write(6,*) ' plasma switch 2 ',sw_use_EUVAC_solar_spectrum
      CALL HL_ML__EUV_ION_PRODUCTION(npts,in,is,o,o2,n2 &
           ,Solar_Declination_Angle_degrees,time,f107, &
           peuvn,gcol_1d,r0,tn,gravity,km_plasma,gr_1d,chi_degrees,1)

!     if (mp.eq.1.and.lp.eq.5) then
!       write(190,*) 'old euv_ion ', UT_in_seconds/3600.,midpoint, &
!                    xtimes(midpoint)
!       do i=in,is
!         write(190,2129) i,(gr_1d(i)-r0)/1000.,peuvn(i,1)
!       enddo
!     endif
! *************************
2129   format(i6,2e12.4)
  endif


              DO 320 i = IN , IS
                  peuvi(i,1) = peuvn(i,1)
                  peuvi(i,2) = 0.0
              320 ENDDO
          !g
          ! calculate the ion- neutral collision frequency, nuin, the frictional
          ! heating term fh (in ev s-1) to be used in routine temp,
          ! and then the ion concentration and flux

          ! O+

          if (idiagnose.eq.1) write(6,*) 'here 6'

              CALL ML__CHEMISTRY_O_plus(chemp_1,beta_1,IN,IS,tn,o,o2,n2, &
              hyd,TI_Oplus_1d,Ti_Hplus_1D,ni_oplus_1d,ni_hplus_1d)

          !c  **
          !c  !include Burnside factor and adjust constant for O-O+
          !c  collision frequency to be compatible with high lat code and
          !c  neutrals
              Burnside_factor=1.0
          !c
              DO 340 i = IN , IS
                  w1 = (TI_Oplus_1d(i)+tn(i))*0.5
                  w2 = 1.04 - 0.067*log10(w1)
              !c  change O-O+ factor to be same as neutral and high lat code
                  wo = 3.42e-17*o(i)*sqrt(w1)*w2*w2*Burnside_factor
              ! wo = 4.45E-17*o(i)*SQRT(w1)*w2*w2
                  wo2 = 6.64E-16*o2(i)
                  wn2 = 6.82E-16*n2(i)
                  nuin(i) = wo + wo2 + wn2
                  frictional_heating_1(i) = (16./32.*wo+32./48.*wo2+28./44.*wn2) &
                  *16.*1.673E-27/1.602E-19
                  nu_mol_n(i)=4.34e-16*n2(i)+4.28e-16*o2(i)+ &
                  2.44e-16*o(i)

                  total_moleculars_1d(i) = no_plus_1d(i)+o2_plus_1d(i)+n2_plus_1d(i)+n_plus_1d(i)

                  nuin_tot(i)=(nuin(i) * ni_Oplus_1d(i) + nu_mol_n(i) * total_moleculars_1d(i)) &
                              / (ni_Oplus_1d(i) + total_moleculars_1d(i))

                  ion_mass(i) = ( 16. * ni_Oplus_1d(i) + 30. * no_plus_1d(i) + &
                                  32. * o2_plus_1d(i) + 28. * n2_plus_1d(i) + 14. * n_plus_1d(i)) & 
                                  * atomic_mass_unit / (ni_Oplus_1d(i)+total_moleculars_1d(i))
              340 ENDDO
          if (idiagnose.eq.1) write(6,*) 'here 7'

          !g ...and solve the diffusion equation for O+....
          !g
              do i=in , is
                  ni_Oplus_1d_saved(i)=ni_Oplus_1d(i)
                  vi_Oplus_1d_saved(i)=vi_Oplus_1d(i)
                  ne_1d_saved(i)=ne_1d(i)
                  dni_Oplus_1d_saved(i)=dni_Oplus_1d(i)
              enddo


              do 2316 i_attempt = 1 , 6

              if (i_attempt == 1) O_plus_production_fudge_factor = 0.0
              if (i_attempt == 2) O_plus_production_fudge_factor = 5.0e-12
              if (i_attempt == 3) O_plus_production_fudge_factor = 5.0e-11
              if (i_attempt == 4) O_plus_production_fudge_factor = 1.0e-10
              if (i_attempt == 5) O_plus_production_fudge_factor = 5.0e-10
              if (i_attempt == 6) O_plus_production_fudge_factor = 1.0e-9

!             if(i_write_out_tube.eq.1) write(6,*) 'i_attempt ',i_attempt

              CALL ML__DIFFUSION_EQUATION_O_PLUS(1,IN,IS,nuin,chemp_1, &
              beta_1,peuvi,g_parallel,dte_1d,dti_Oplus_1d,dti_Hplus_1d,upar_apex, &
              eta_apex_1d,dq_1d,TI_Oplus_1d,TI_Hplus_1d,TE_1d,ni_OPlus_1d,ni_Hplus_1d, &
              dni_oplus_1d,dni_hplus_1d,Vi_oplus_1d,vi_hplus_1d,NE_1d,div_vperp_1d, &
              O,M_plasma(1),M_plasma(2),KM_plasma,plasma_time_step_seconds,Ofailed, &
              i_write_out_tube,altitude_PZ_km,O_plus_production_fudge_factor)
 
              if (Ofailed == 1) then
                  if (i_attempt == 5) write(6,*) 'Ofailed  (5th attempt) ' , mp , lp
                  do i = in , is
                      ni_Oplus_1d(i) = ni_Oplus_1d_saved(i)
                      vi_Oplus_1d(i) = vi_Oplus_1d_saved(i)
                      ne_1d(i) = ne_1d_saved(i)
                      dni_Oplus_1d(i) = dni_Oplus_1d_saved(i)
                  enddo
              else
                  goto 2317
              endif

 2316         continue

 2317         continue




              if (idiagnose.eq.1) write(6,*) 'here 8'


          !g
          !g same as above but for H+......
          !g
              CALL ML__CHEMISTRY_H_plus(IN,IS,tn,o,hyd,Ti_Hplus_1d,ni_oplus_1d,chemp_2,beta_2)


              DO 360 i = IN , IS
                  w1 = 1. - 0.047*log10(TI_Hplus_1d(i))
                  wo = 6.61E-17*o(i)*SQRT(TI_Hplus_1d(i))*w1*w1
                  wo2 = 3.2E-15*o2(i)
                  wn2 = 3.36E-15*n2(i)
                  nuin(i) = wo + wo2 + wn2
                  frictional_heating_2(i) = (16./17.*wo+32./33.*wo2+28./29.*wn2) &
                  *1.673E-27/1.602E-19
              360 ENDDO
          if (idiagnose.eq.1) write(6,*) 'here 9'
          !g
          !g ...and then the diffusion equation for H+....
          !g
              do i= in , is
                  ni_Hplus_1d_saved(i)=ni_Hplus_1d(i)
                  vi_Hplus_1d_saved(i)=vi_Hplus_1d(i)
                  ne_1d_saved(i)=ne_1d(i)
                  dni_Hplus_1d_saved(i)=dni_Hplus_1d(i)
              enddo

              CALL ML__DIFFUSION_EQUATION_H_PLUS(2,IN,IS,nuin,chemp_1,chemp_2, &
              beta_1,beta_2,peuvi,g_parallel,dte_1d,dti_Oplus_1d,dti_Hplus_1d,upar_apex, &
              eta_apex_1d,dq_1d,TI_Oplus_1d,Ti_Hplus_1d,TE_1d,ni_oplus_1d,ni_hplus_1d, &
              dni_oplus_1d,dni_hplus_1d,VI_oplus_1d,vi_hplus_1d,NE_1d,div_vperp_1d, &
              O,M_plasma(1),M_plasma(2),KM_plasma,plasma_time_step_seconds,Hfailed, &
              i_write_out_tube,altitude_PZ_km,O_plus_production_fudge_factor)

          if (idiagnose.eq.1) write(6,*) 'here 10'
          !g
          !g  if tube has a negative density - set the ion density
          !g  to the saved values...
          !g
              if (Hfailed == 1) then
                  write(6,*) 'Hfailed ' , mp , lp
                  do i = in , is 
                      ni_Oplus_1d(i) = ni_Oplus_1d_saved(i)
                      ni_Hplus_1d(i) = ni_Hplus_1d_saved(i)
                      vi_Hplus_1d(i) = vi_Hplus_1d_saved(i)
                      ne_1d(i) = ne_1d_saved(i)
                      dni_Hplus_1d(i) = dni_Hplus_1d_saved(i)
                  enddo
              endif

              ite_failed = 0
              iti_failed = 0
          !g
          !g  ENERGY_EQUATIONS now calculates Ti (for both O+ and H+) and Te....
          !g
              plasma_energy_time_step_seconds=plasma_time_step_seconds
          !g
          !g
          !g After having calculated the O+ and H+ densities, and the ion and electron temperatures
          !g calculate the molecular ion densities.....
          !g
          !g Need oplus, te, and ti as 1D parameters...
          !g
              do i = in , is
              !g
              !g originally there was no Oplus calculated below 130km (in plasma)
              !g and thus it was set to zero.  It was then given a value in molecular ions.
              !g The following simulates this to check if this explains the problem...
              !g

!                   if (altitude_pz_km(i) > 130) then
!                       ni_Oplus_1d(i) = ni_Oplus_1d(i)
!                   else
!                       ni_Oplus_1d(i) = 0.0
!                   endif
              !g
                  cos_SZA(i) = cos(chi_degrees(i) * dtr)
              !g
              !g qin is the TIROS ionisation rate. set to zero here because all of these
              !g tubes are sub-auroral - therefore no particle precipitation....
              !g
                  qin(i) = 0.0
              enddo
          !g
                do i=in, is
                  ni_Oplus_1d_saved(i) = ni_Oplus_1d(i)
                  no_plus_1d_saved(i) = no_plus_1d(i)
                  o2_plus_1d_saved(i) = o2_plus_1d(i)
                  n2_plus_1d_saved(i) = n2_plus_1d(i)
                  n_plus_1d_saved(i) = n_plus_1d(i)
                enddo

              do i = IN , IS
                  euv_production_rate_Oplus(i) = peuvn(i,1)
                  euv_production_rate_N2plus(i) = peuvn(i,4)
                  euv_production_rate_O2plus(i) = peuvn(i,5)
                  auroral_production_rate_Oplus(i) = 0.0
                  auroral_production_rate_N2plus(i) = 0.0
                  auroral_production_rate_O2plus(i) = 0.0
              enddo

          if (idiagnose.eq.1) write(6,*) 'here 12'
     call HL_ML__MOLECULAR_IONS_ON_TUBES( &
                                GIP_switches, &
                                npts,IN,IS,f107,O,O2,N2,no_1d,n4s_1d,n2d_1d, &
                                tn,te_1d,ti_Oplus_1d,ni_Oplus_1d, &
                                altitude_pz_km,cos_SZA, &
                                euv_production_rate_Oplus,euv_production_rate_N2plus,euv_production_rate_O2plus, &
                                auroral_production_rate_Oplus,auroral_production_rate_N2plus,auroral_production_rate_O2plus, &
                                N2_PLUS_1d,NO_PLUS_1d,O2_PLUS_1d,N_PLUS_1d, &
                                i_failed_molecular_ions)
          if (idiagnose.eq.1) write(6,*) 'here 13'


              if(i_failed_molecular_ions == 1 ) then
!                 write(6,*) 'molecular_ions failed '
                do i=in, is
                  ni_Oplus_1d(i) = ni_Oplus_1d_saved(i)
                  no_plus_1d(i) = no_plus_1d_saved(i)
                  o2_plus_1d(i) = o2_plus_1d_saved(i)
                  n_plus_1d(i) = n_plus_1d_saved(i)
                  n2_plus_1d(i) = n2_plus_1d_saved(i)
                enddo
              endif
          !g
          !g  Here is the new Apex ML__FIELD_LINE_INTEGRALS call....
          !g

              do i = in , is
           !    if (altitude_pz_km(i).lt.130..and.lp.ge.1.and.lp.le.7) then     ! lats 50 to 60
           !    if (altitude_pz_km(i).lt.130..and.lp.ge.8.and.lp.le.14) then    ! lats 40 to 50
           !    if (altitude_pz_km(i).lt.130..and.lp.ge.15.and.lp.le.20) then   ! lats 30 to 40
           !    if (altitude_pz_km(i).lt.130..and.lp.ge.21.and.lp.le.28) then   ! lats 20 to 30
           !    if (altitude_pz_km(i).lt.130..and.lp.ge.29.and.lp.le.39) then   ! lats 10 to 20
           !    if (altitude_pz_km(i).lt.130..and.lp.ge.40.and.lp.le.67) then   ! lats 0 to 10
                  U_zonal_for_integrals(i) = U_zonal(i)
                  U_merid_for_integrals(i) = U_merid(i)
                  U_vertical_for_integrals(i) = U_vertical(i)
           !    else
           !      U_zonal_for_integrals(i) = 0.0
           !      U_merid_for_integrals(i) = 0.0
           !      U_vertical_for_integrals(i) = 0.0
           !    endif
              enddo

              CALL ML__FIELD_LINE_INTEGRALS(in,is,distance_ds,midpoint,Apex_d1_1d, &
                                               Apex_d2_1d,Apex_BE3_1d, &
                                               apex_D_1d,apex_d1d1_1d,apex_d1d2_1d,apex_d2d2_1d,apex_BMAG_1d, &
                                               ni_oplus_1d,ni_hplus_1d,ti_oplus_1d, &
                                               no_plus_1d,o2_plus_1d,n2_plus_1d,n_plus_1d,tn,o,o2,n2, &
                                               U_zonal_for_integrals,U_merid_for_integrals,U_vertical_for_integrals, &
                                               sigma_ped_1d,sigma_hall_1d, &
                                               sigma_phph_dsi,sigma_lmlm_msi,sigma_h,sigma_c,Kdmph_dsi,Kdmlm)

          if (idiagnose.eq.1) write(6,*) 'here 14'
          !g
          !g
          !g
          !g  Thats it for this flux-tube....
          !g

  return

end SUBROUTINE ML__single_flux_tube_1D_calculation
















SUBROUTINE ML__check_for_negative_densities(in,is,ni,vi,ne)

  IMPLICIT NONE

  INTEGER ::  mp , lp , i

  INTEGER :: IN(NMP,NLP) , IS(NMP,NLP)
  REAL(kind=8) NI(NPTS,NMP,2)
  REAL(kind=8) VI(NPTS,NMP,2)
  REAL(kind=8) NE(NPTS,NMP)

  do mp=1,nmp
      do lp=1,nlp
          do i=in(mp,lp),is(mp,lp)
              if(ni(i,mp,1) < 0.0) then
                  write(6,*) 'Negative density O+',i,mp,lp
                  goto 6790
              elseif(ni(i,mp,2) < 0.0) then
                  write(6,*) 'Negative density H+',i,mp,lp
                  goto 6790
              endif
          enddo
          goto 6799
          6790 do i= in(mp,lp) , is(mp,lp)
              ni(i,mp,1) = 1.e11
              ni(i,mp,2) = 1.e5
              vi(i,mp,1) = 0.
              vi(i,mp,2) = 0.
              ne(i,mp) = 1.e11
          enddo
          6799 continue
      enddo
  enddo

  return

end SUBROUTINE ML__check_for_negative_densities








SUBROUTINE ML__FASTHL__NEARHTPLA(ARR,VALue,M_IN,N_IN,M,N,nearst,IFAiled)

! this routine finds which two elements of ARR
! surround Value.  Arr is a monotonically increasing array
! The start and end points for the search are M_IN and N_IN.
! The exit points surrounding Value are N and M
  IMPLICIT NONE
  REAL(kind=8) :: ARR , VALue , dist1 , dist2
  INTEGER :: i , IFAiled , M , N ,  &
  M_IN , N_IN , n1 , nearst
  DIMENSION ARR(NPTS)

  N=N_IN
  M=M_IN
  IFAiled = 0
  if((value <= arr(n) .AND. value >= arr(m)) .OR. &
  (value >= arr(n) .AND. value <= arr(m))) then
      234 n1 =(n - m)/2 + m
      IF ( ARR(n1) > VALue ) then
          n=n1
      else
          m=n1
      ENDIF
      if(abs(n-m) > 1) goto 234
      dist1=abs(arr(n)-value)
      dist2=abs(arr(m)-value)
      if(dist1 < dist2) then
          nearst=n
      else
          nearst=m
      endif
  else
      ifailed=1
  endif

  RETURN



end SUBROUTINE ML__FASTHL__NEARHTPLA





SUBROUTINE ML__get_empirical_vertical_exb(apex_height_km,xmlat_300,xmlon_300,lt,geolon_apex, &
  ut,iday,f107,e2_factor_300km_to_apex, &
  v_upwards_at_apex)

  implicit none
  REAL(kind=8), intent(in)  :: apex_height_km
  REAL(kind=8), intent(in)  :: xmlat_300
  REAL(kind=8), intent(in)  :: xmlon_300
  REAL(kind=8), intent(in)  :: lt
  REAL(kind=8), intent(in)  :: geolon_apex
  REAL(kind=8), intent(in)  :: ut
  REAL(kind=8), intent(in)  :: f107
  REAL(kind=8), intent(in)  :: e2_factor_300km_to_apex
  REAL(kind=8) v_upwards_at_apex_fejer
  REAL(kind=8) v_upwards_at_apex_richmond
  REAL(kind=8) ve
  REAL(kind=8) pot
  REAL(kind=8) dayno
  REAL(kind=8) factor_less_than_300
  REAL(kind=8) apex_factor
  REAL(kind=8) v_outwards_at_300km
  INTEGER :: iday
  INTEGER :: iseasav
  INTEGER :: iutav
  REAL(kind=8), intent(out)  :: v_upwards_at_apex
!g



!g  This combines the low-latitude electric field models
!g  of Richmond and Fejer.
!g  The Fejer model is used for flux-tubes with apex heights up to 1000km.
!g  The Richmond model is used for flux-tubes with apex heights greater than
!g  2000km.  Between 1000km and 2000km a linear interpolation between the two
!g  is used.
!g
  if(apex_height_km < 2000.) then
  !g
  !g  use the Fejer model....
  !g
      call ML__fejer_exb_model(LT,Geolon_apex,IDAY,F107,v_upwards_at_apex_fejer)
  !g
  !g  if the apex height is less than 300km then tail the upwards drift off
  !g  in height between a full value at 300km and zero at 100km....
  !g
      if(apex_height_km < 300.) then
          factor_less_than_300 = (apex_height_km - 100. ) / 200.
          v_upwards_at_apex_fejer = factor_less_than_300 * v_upwards_at_apex_fejer
          if(apex_height_km < 100.) v_upwards_at_apex_fejer = 0.0
      endif
  !g
      v_upwards_at_apex = v_upwards_at_apex_fejer
  !g
  endif

  if(apex_height_km > 1000.) then
  !g
  !g  use the Richmond model....
  !g
      Iseasav = 3
      Iutav = 0
      dayno = float(iday)
      call ML__Richmond_lowlat_Efield_model(XMLat_300,XMLon_300,dayno,UT,ISEasav,IUTav, &
      POT,v_outwards_at_300km,VE)
  !g
      v_upwards_at_apex_richmond = v_outwards_at_300km * e2_factor_300km_to_apex
  !g
      v_upwards_at_apex = v_upwards_at_apex_richmond
  !g
  endif

  if(apex_height_km > 1000. .AND. apex_height_km < 2000.) then
  !g
  !g  interpolate between the 2 models according to the apex height....
  !g
      apex_factor = (apex_height_km - 1000.) / 1000.
      v_upwards_at_apex = apex_factor * (v_upwards_at_apex_richmond - v_upwards_at_apex_fejer) &
      + v_upwards_at_apex_fejer
  endif
!g
  return



end SUBROUTINE ML__get_empirical_vertical_exb













SUBROUTINE ML__TUBES_LINEAR_INTERPOLATE(re_apex,vpeq,dt,ni &
  ,q,vi,ti,te,ne,ni_save,vi_save,ti_save,te_save &
  ,mp,lp,in,is)

!*********************************************
! *
! Interpolation routine for incorporating  *
! EXB drift into the 'fixed tubes' version *
! *
!*********************************************

  IMPLICIT NONE
  REAL(kind=8) :: BLON , DT , Q , TE , TI , VI , vpeq , vzeq , &
  vpeq_tim , vt300
  ! removedNMP, NLP, NPTS declaration mjh
  INTEGER :: i , lp , mp , istop

  REAL(kind=8) :: NI , RE_apex(NMP,NLP) , NE , &
  VU(120,NMP,NLP) , &
  REDummy(NMP,NLP)
  real(kind=8) :: reback, factor,pout,pin,factor2,ni1_in(NPTS),ni2_in(NPTS), &
  vi1_in(npts),vi2_in(npts)
  real(kind=8) :: ni1_out(NPTS),ni2_out(NPTS), &
  vi1_out(npts),vi2_out(npts), &
  ni1old(npts)
  real(kind=8) :: ti1_in(npts),ti2_in(npts),ti1_out(npts), &
  ti2_out(npts),te_in(npts),te_out(npts), &
  ni_save(npts,nmp,2),vi_save(npts,nmp,2), &
  ti_save(npts,nmp,2),te_save(npts,nmp) , &
  xtime_save(nmp,nlp), &
  lambda_m_inner,lambda_m_outer,lambda_m_stepped_backwards,r0
  PARAMETER (R0=6.370E06)
  integer :: lp_out,lp_in,lpdum,ip,isouth,inorth, &
  iww,ispecial
  INTEGER :: IACtive(NMP,NLP)
  INTEGER ::   IN(NMP,NLP) , IS(NMP,NLP)
  DIMENSION  BLOn(NMP,NLP)
  DIMENSION  Q(NPTS,NMP) , &
  TI(NPTS,NMP,2) , TE(NPTS,NMP) , &
  NI(NPTS,NMP,2) , VI(NPTS,NMP,2) , &
  NE(NPTS,NMP)
  DIMENSION &
  vpeq_tim(nmp,nlp) , vpeq(nmp,nlp) , &
  vt300(nmp,nlp) , vzeq(nmp,nlp)

  REAL(kind=8) :: sqrt_part

! step backwards to imagined previous flux-tube position....

! stop=1
! f(istop.eq.1) stop
  iww=0
! f(mp.eq.12.and.lp.eq.10) iww=1
!g
  reback=re_apex(mp,lp)-vpeq(mp,lp)*DT
!g
! if(iww.eq.1) write(88,*) mp,lp,nlp,
! &  	re(mp,lp),reback,vpeq(mp,lp),dt,ne(125,mp,lp),
! &  in(mp,lp),is(mp,lp)
  istop=1
! f(istop.eq.1) stop

! we now need parameters on this imagined tube (by interpolation).....

  do lpdum=1,NLP
      if(reback > re_apex(mp,lpdum)) then
          lp_out=lpdum-1
          lp_in=lpdum
          if(iww == 1) write(88,*) lp_out,lp_in
          goto 3456
      endif
  enddo
  lp_out=NLP-1
  lp_in=NLP
  3456 continue
  if(lp_out == 0) then
      lp_in = 2
      lp_out = 1
  endif
!g
! factor=(reback-re(mp,lp_in))/(re(mp,lp_out)-re(mp,lp_in))
!g
!g  For the apex code we use the modified apex latitude as our parameter for
!g  interpolating in the inwards/outwards direction....
!g

  lambda_m_inner = acos(sqrt((r0+90000.)/re_apex(mp,lp_in)))
  lambda_m_outer = acos(sqrt((r0+90000.)/re_apex(mp,lp_out)))

!cg   For the very small flux tubes it is possible that stepping backwards
!cg   leads to an apex height which is lower than the base height of 90km.
!cg   This leads to a NAN because the ACOS won't compute.
!cg   So we do the 'sqrt_part' bit below and check that this is less than 1.0

  sqrt_part = sqrt((r0+90000.)/reback)

  if(sqrt_part.le.0.999999) then
       lambda_m_stepped_backwards = acos(sqrt_part)
  else
       lambda_m_stepped_backwards = 0.0
  endif

  factor=(lambda_m_stepped_backwards - lambda_m_inner) / &
  (lambda_m_outer - lambda_m_inner)
!g
!g  Need to flag when factor gets too big or small.....
!g  Factor of less than zero means that the position of
!g  the imagined tubes is inside of the inner most real tube.
!g  The interpolation then becomes extrapolation (which isn't
!g  necessarily a problem in itself but can produce nasty
!g  negative densities if the gradients are large enough).
!g  In this case set the factor to zero which means that you
!g  end up just using the inner flux_tube values...
!g  (Should be safe)
!g
  if(factor < 0.0) then
  ! write(6,*) 'FACTOR ',factor
  ! write(6,*) mp,lp_in,lp_out
  ! write(6,*) reback,re(mp,lp_in),re(mp,lp_out)
  ! write(6,*) '     '
      factor = 0.0
  endif
!g
!g  The same could be true if factor is greater than 1.0
!g  The extrapolation isn't in itself a problem but large
!g  negative gradients in density (temp, whatever) could
!g  lead to negative values.  For the time being lets just
!g  flag this and see if it is ever a problem....
!g
  if(factor > 1.0) then
      write(6,*) 'FACTOR ',factor
      write(6,*) mp,lp_in,lp_out
      write(6,*) reback,re_apex(mp,lp_in),re_apex(mp,lp_out)
      write(6,*) '     '
  endif
!g
!g first interpolate in q for the inner tube.....
!g
  if(iww == 1) write(88,*) '++++++++++ INNER TUBE ++++++++++++'
  do 7000 ip=in(mp,lp),is(mp,lp)
      if(iww == 1) write(88,*) 'point ',ip
      ni1old(ip)=ni(ip,mp,1)
      ispecial=0
  !g
  !g loop over the inner tube ......
  !g
      do i=in(mp,lp_in),is(mp,lp_in)
          if(q(i,mp) < q(ip,mp)) then
              isouth=i
              inorth=i-1
              goto 3488
          endif
      enddo
  !g
      ispecial=1
      isouth=is(mp,lp_in)
      inorth=is(mp,lp_in)-1
      3488 continue
      if(isouth == in(mp,lp_in)) then
          isouth=in(mp,lp_in)+1
          inorth=in(mp,lp_in)
          ispecial=2
      endif
      24 format(i4,2x,3f9.5)
  !g
      if(iww == 1) write(88,*) 'ispecial ',ispecial
      if(ispecial == 0) then
          factor2=(q(ip,mp)-q(isouth,mp))/ &
          (q(inorth,mp)-q(isouth,mp))
          if(iww == 1) write(88,*) 'factor2 inner',factor2
          ni1_in(ip)=(factor2*(ni(inorth,mp,1) - &
          ni(isouth,mp,1))) + ni(isouth,mp,1)
          ni2_in(ip)=(factor2*(ni(inorth,mp,2) - &
          ni(isouth,mp,2))) + ni(isouth,mp,2)
          vi1_in(ip)=(factor2*(vi(inorth,mp,1) - &
          vi(isouth,mp,1))) + vi(isouth,mp,1)
          vi2_in(ip)=(factor2*(vi(inorth,mp,2) - &
          vi(isouth,mp,2))) + vi(isouth,mp,2)
          ti1_in(ip)=(factor2*(ti(inorth,mp,1) - &
          ti(isouth,mp,1))) + ti(isouth,mp,1)
          ti2_in(ip)=(factor2*(ti(inorth,mp,2) - &
          ti(isouth,mp,2))) + ti(isouth,mp,2)
          te_in(ip)=(factor2*(te(inorth,mp) - &
          te(isouth,mp))) + te(isouth,mp)
      elseif(ispecial == 1) then
          ni1_in(ip)=ni(isouth,mp,1)
          ni2_in(ip)=ni(isouth,mp,2)
          vi1_in(ip)=vi(isouth,mp,1)
          vi2_in(ip)=vi(isouth,mp,2)
          ti1_in(ip)=ti(isouth,mp,1)
          ti2_in(ip)=ti(isouth,mp,2)
          te_in(ip)=te(isouth,mp)
      elseif(ispecial == 2) then
          ni1_in(ip)=ni(inorth,mp,1)
          ni2_in(ip)=ni(inorth,mp,2)
          vi1_in(ip)=vi(inorth,mp,1)
          vi2_in(ip)=vi(inorth,mp,2)
          ti1_in(ip)=ti(inorth,mp,1)
          ti2_in(ip)=ti(inorth,mp,2)
          te_in(ip)=te(inorth,mp)
      endif
      if(iww == 1) write(88,*) factor2,isouth,inorth,ni1_in(ip)
  7000 ENDDO
!g
!g then interpolate in q for the outer tube.....
!g
  if(iww == 1) write(88,*) '++++++++++ OUTER TUBE ++++++++++++'
  do 8000 ip=in(mp,lp),is(mp,lp)
      if(iww == 1) write(88,*) 'point ',ip
      ispecial=0
      do i=in(mp,lp_out),is(mp,lp_out)
          if(q(i,mp) < q(ip,mp)) then
              isouth=i
              inorth=i-1
              goto 3489
          endif
      enddo
  !g
      ispecial=1
      isouth=is(mp,lp_out)
      inorth=is(mp,lp_out)-1
      3489 continue
      if(isouth == in(mp,lp_out)) then
          isouth=in(mp,lp_out)+1
          inorth=in(mp,lp_out)
          ispecial=2
      endif
  !g
      if(iww == 1) write(88,*) 'ispecial ',ispecial
      if(ispecial == 0) then
          factor2=(q(ip,mp)-q(isouth,mp))/ &
          (q(inorth,mp)-q(isouth,mp))
          if(iww == 1) write(88,*) 'factor2 outer',factor2
          ni1_out(ip)=(factor2*(ni(inorth,mp,1) - &
          ni(isouth,mp,1))) + ni(isouth,mp,1)
          ni2_out(ip)=(factor2*(ni(inorth,mp,2) - &
          ni(isouth,mp,2))) + ni(isouth,mp,2)
          vi1_out(ip)=(factor2*(vi(inorth,mp,1) - &
          vi(isouth,mp,1))) + vi(isouth,mp,1)
          vi2_out(ip)=(factor2*(vi(inorth,mp,2) - &
          vi(isouth,mp,2))) + vi(isouth,mp,2)
          ti1_out(ip)=(factor2*(ti(inorth,mp,1) - &
          ti(isouth,mp,1))) + ti(isouth,mp,1)
          ti2_out(ip)=(factor2*(ti(inorth,mp,2) - &
          ti(isouth,mp,2))) + ti(isouth,mp,2)
          te_out(ip)=(factor2*(te(inorth,mp) - &
          te(isouth,mp))) + te(isouth,mp)
      elseif(ispecial == 1) then
          ni1_out(ip)=ni(isouth,mp,1)
          ni2_out(ip)=ni(isouth,mp,2)
          vi1_out(ip)=vi(isouth,mp,1)
          vi2_out(ip)=vi(isouth,mp,2)
          ti1_out(ip)=ti(isouth,mp,1)
          ti2_out(ip)=ti(isouth,mp,2)
          te_out(ip)=te(isouth,mp)
      elseif(ispecial == 2) then
          ni1_out(ip)=ni(inorth,mp,1)
          ni2_out(ip)=ni(inorth,mp,2)
          vi1_out(ip)=vi(inorth,mp,1)
          vi2_out(ip)=vi(inorth,mp,2)
          ti1_out(ip)=ti(inorth,mp,1)
          ti2_out(ip)=ti(inorth,mp,2)
          te_out(ip)=te(inorth,mp)
      endif
      if(iww == 1) write(88,*) factor2,isouth,inorth,ni1_out(ip)
      if(iww == 1) write(88,*) ni(inorth,mp,1), &
      ni(isouth,mp,1),inorth,isouth,mp,lp_out
  8000 ENDDO
! if(istop.eq.1) stop
!g
!g linear interpolation ......
!g
  do 9000 i=in(mp,lp),is(mp,lp)
  ! pout=log10(ni1_out(i))
  ! pin=log10(ni1_in(i))
  ! ni_save(i,mp,lp,1)=10.**(((pout-pin)*factor)+pin)
  ! pout=log10(ni2_out(i))
  ! pin=log10(ni2_in(i))
  ! ni_save(i,mp,lp,2)=10.**(((pout-pin)*factor)+pin)
      pout=ni1_out(i)
      pin=ni1_in(i)
      ni_save(i,mp,1)=((pout-pin)*factor)+pin
      pout=ni2_out(i)
      pin=ni2_in(i)
      ni_save(i,mp,2)=((pout-pin)*factor)+pin
      pout=vi1_out(i)
      pin=vi1_in(i)
      vi_save(i,mp,1)=((pout-pin)*factor)+pin
      pout=vi2_out(i)
      pin=vi2_in(i)
      vi_save(i,mp,2)=((pout-pin)*factor)+pin
      pout=ti1_out(i)
      pin=ti1_in(i)
      ti_save(i,mp,1)=((pout-pin)*factor)+pin
      pout=ti2_out(i)
      pin=ti2_in(i)
      ti_save(i,mp,2)=((pout-pin)*factor)+pin
      pout=te_out(i)
      pin=te_in(i)
      te_save(i,mp)=((pout-pin)*factor)+pin
  !g
  !g  If any of these parameters have gone -ve we got some
  !g  trub.........
  !g
  ! if(ni_save(i,mp,lp,1).lt.0.0) then
  ! write(6,*) 'Lin Interp O+ ',i,mp,lp
  ! write(6,*) '   ',pout,pin,factor
  ! elseif(ni_save(i,mp,lp,2).lt.0.0) then
  ! write(6,*) 'Lin Interp H+ ',i,mp,lp
  ! write(6,*) '   ',pout,pin,factor
  ! endif
  9000 ENDDO
  return



end SUBROUTINE ML__TUBES_LINEAR_INTERPOLATE



SUBROUTINE ML__get_U_apex(in,is,Apex_bhat_1d,U_zonal,U_merid,U_vert,upar,uperp, &
  upar_zonal_contribution,upar_meridional_contribution,upar_vertical_contribution)

!-------------------------------------------------------------------------
  implicit none

  integer(kind=4), intent(in)  ::  in    ! North footpoint of a flux tube
  integer(kind=4), intent(in)  ::  is    ! South footpoint of a flux tube

! input Arguments
  real(kind=8), intent(in)  ::  Apex_bhat_1d(3,npts)  ! Calculated in initialise_plasma

  real(kind=8), intent(in)  ::  U_zonal(npts)    !neutral wind +ggeast  [m/s]
  real(kind=8), intent(in)  ::  U_merid(npts)    !neutral wind +ggsouth
  real(kind=8), intent(in)  ::  U_vert(npts)     !neutral wind +ggup

! Output Arguments
  real(kind=8),  intent(out) ::  Upar(npts)              !Un parallel [m/s]
  real(kind=8),  intent(out) ::  upar_zonal_contribution(npts)      ! [m/s]
  real(kind=8),  intent(out) ::  upar_meridional_contribution(npts) ! [m/s]
  real(kind=8),  intent(out) ::  upar_vertical_contribution(npts)   ! [m/s]
  real(kind=8),  intent(out) ::  Uperp(npts)             !Un perp

!----------------------------Local variables-----------------------------

  integer(kind=4) ::  ipts              ! npts increment
  real(kind=8) ::  w1                ! | Un |^2
  real(kind=8) ::  upar_squared      ! upar^2
!------------------------------------------------------------------------

  ipts_loop: do ipts=in , is

! Ue3=d3*u
!g  Here is Naomis original definition...
! upar(ipts)=  apex_bhat(1,ipts,mp,lp)*U_zonal(ipts)
! &            -apex_bhat(2,ipts,mp,lp)*U_merid(ipts)
! &            +apex_bhat(3,ipts,mp,lp)*U_vert(ipts)
!g  Looking at the output I reckon that the meridional contribution needs turning around...
!g  ..in the plasma code U-merid is southwards and Upar is also southwards....
!g  (not sure about the others yet)....
!g  I'm sure now - looked at the numbers and the definition above - is exactly the wrong way
!g  round.  - ie, upar = -b1 * Uzonal + b2 * Umerid - b3 * Uvert  .....
!g
  upar(ipts)=  0.0 - apex_bhat_1d(1,ipts)*U_zonal(ipts) &
  +apex_bhat_1d(2,ipts)*U_merid(ipts) &
  - apex_bhat_1d(3,ipts)*U_vert(ipts)
  upar_zonal_contribution(ipts) = - apex_bhat_1d(1,ipts)*U_zonal(ipts)
  upar_meridional_contribution(ipts) = apex_bhat_1d(2,ipts)*U_merid(ipts)
  upar_vertical_contribution(ipts) = - apex_bhat_1d(3,ipts)*U_vert(ipts)
!g
!g I think the above line is correct - but doing experiments for the moment with the
!g field-aligned wind set to zero....
!g
! upar(ipts) = 0.0
!g
!g  This takes upar from just the merdional wind contribution (for testing)....
!g  ....again Naomis original minus sign has been reversed....
!g
! upar(ipts)= apex_bhat(2,ipts,mp,lp)*U_merid(ipts)


! Uperp: exctracted from thermosphere_input2.f
  w1 = U_merid(ipts)*U_merid(ipts) &
  +U_zonal(ipts)*U_zonal(ipts) &
  + U_vert(ipts)*U_vert(ipts)
  upar_squared = upar(ipts) * upar(ipts)
  if(w1>upar_squared) then
      UPErp(ipts) = SQRT(w1 - upar_squared)
  else
      UPErp(ipts) = 0.0
  ! write(6,*)'uperp set to zero',ipts
  ! &        ,U_zonal(ipts),apex_bhat(1,ipts)
  ! &        ,U_merid(ipts),apex_bhat(2,ipts)
  ! &        ,U_vert(ipts) ,apex_bhat(3,ipts)
  endif


enddo  ipts_loop




end SUBROUTINE ML__get_U_apex






SUBROUTINE ML__ION_TEMP_PLASMA(TN,TE,N1,N2,TI, &
                             NIT,O,O2,U2Dif,in,is)
  IMPLICIT NONE
  INTEGER :: in , is 

  REAL(kind=8) :: &
  cf1a , CF1n(NPTS) , cf1nk1 , cf1nk2 , cf1nk3 , factor , O(NPTS) , &
  O2(NPTS) , RTTin(NPTS) , rttp , t1 , t2 , TE , TI , &
  TINLOG(NPTS) , TN , tp , tplog
  REAL(kind=8) :: u2 , U2Dif(NPTS) , w , w1 , w2 , w3 , w4 , w5 , w6 , w8 , x
  INTEGER :: l

! for the purposes of this routine it is necessary to work in the cg
! system.we alter the coefficients to express the number densities a
! velocities in the right units.

  REAL(kind=8) :: N1 , N2 , NIT(NPTS)
  DIMENSION TN(NPTS) , TE(NPTS) , N1(NPTS) , N2(NPTS) , TI(NPTS)
!c  **
!c  **
  factor = 1.0
!c  **
  DO 100 l = in , is
      u2 = U2Dif(l)*1.0E4

      if(o(l) < 1.d-50) o(l)=1.d-50
      if(o2(l) < 1.d-50) o2(l)=1.d-50
      if(nit(l) < 1.d-50) nit(l)=1.d-50
      cf1nk1 = 3.42E-17*O(l)*factor
      cf1nk2 = 6.66E-16*O2(l)
      cf1nk3 = 6.82E-16*NIT(l)
      w1 = 4.8E-13*(N1(l)+N2(l))
      w2 = SQRT(TE(l))
      w3 = w1/w2
      w4 = w1/(w2*w2*w2)
      w5 = 2.1E-21*O(l)
      w6 = (6.6E-20*NIT(l)+5.8E-20*O2(l))

  ! we now solve the heat balance equation for ti,using newton-raphson

      t1 = TE(l)
      50 w = w5*SQRT(t1+TN(l)) + w6

  ! frictional heating terms.

      tp = (t1+TN(l))/2.
      rttp = SQRT(tp)
      tplog = log10(tp)
      x = 1.04 - 0.067*tplog
      cf1a = cf1nk1*rttp*x*x
      w8 = (8.37E-12*cf1a+1.12E-11*cf1nk2+1.06E-11*cf1nk3)*u2
      t2 = (w*TN(l)+w3+w8)/(w+w4)
      IF ( ABS(t2-t1) < 0.5 ) THEN
      !c  **
      !c  modified by tjfr jan 86
      !c  **
          IF ( t2 < TN(l) ) t2 = TN(l)
          TI(l) = t2
          tp = (t2+TN(l))/2.
          rttp = SQRT(tp)
          tplog = log10(tp)
          RTTin(l) = rttp
          TINlog(l) = tplog
          CF1n(l) = cf1a + cf1nk2 + cf1nk3
      ELSE
          t1 = t2
          GOTO 50
      ENDIF
  100 ENDDO
  RETURN



end SUBROUTINE ML__ION_TEMP_PLASMA
































SUBROUTINE ML__DF_BY_DS(F,DF,IN,IS,ETA,DQ)

!***********************************************************************



! differentiates f with respect to arc length s and puts
! result in df
! df/ds = df/dq * dq/ds    and dq/ds=eta
! the data is not equally spaced (wrt q) and so the array dq contains
! the spaces between succesive q with dq(i) = q(i+1)-q(i)
!***********************************************************************

  IMPLICIT NONE
  REAL(kind=8) :: DF , DQ , ETA , F , w , w0 , w1 , w2 , xrs1
  INTEGER :: i , IN , in1 , in2 , IS , is1 , is2 

  DIMENSION F(NPTS) , DF(NPTS) , ETA(NPTS) , DQ(NPTS)

  in1 = IN + 1
  in2 = IN + 2
  is1 = IS - 1
  is2 = IS - 2

  w0 = -F(IN)*(2*DQ(IN)+DQ(in1))/(DQ(IN)*(DQ(IN)+DQ(in1)))
  w1 = F(in1)*(DQ(in1)+DQ(IN))/(DQ(IN)*DQ(in1))
  w2 = -F(in2)*DQ(IN)/(DQ(in1)*(DQ(in1)+DQ(IN)))
  DF(IN) = (w0+w1+w2)*ETA(IN)

  DO 100 i = in1 , is1
      xrs1 = DQ(i-1)*DQ(i)
      w0 = -F(i-1)*DQ(i)/(DQ(i-1)*DQ(i-1)+xrs1)
      w1 = F(i)*(DQ(i)-DQ(i-1))/xrs1
      w2 = F(i+1)*DQ(i-1)/(DQ(i)*DQ(i)+xrs1)
      DF(i) = (w0+w1+w2)*ETA(i)
  100 ENDDO

  w = DQ(is1)*DQ(is2)
  w0 = F(is2)*DQ(is1)/(DQ(is2)*DQ(is2)+w)
  w1 = -F(is1)*(DQ(is1)+DQ(is2))/w
  w2 = F(IS)*(DQ(is2)+2*DQ(is1))/(DQ(is1)*DQ(is1)+w)
  DF(IS) = (w0+w1+w2)*ETA(IS)

  RETURN



end SUBROUTINE ML__DF_BY_DS


















SUBROUTINE ML__CHEMISTRY_O_plus(chemp_1,beta_1,IN,IS,TN,O,O2,N2,HYD, &
                       ti_oplus_1d,ti_hplus_1d,ni_oplus_1d,ni_hplus_1d)

!***********************************************************************
! routine to calculate chemical production rates and loss coefficients
! these are stored in chemp and beta respectively
!***********************************************************************

  IMPLICIT NONE
  INTEGER :: i , IN , IS

  real(kind=8) :: w , wk1(NPTS) , wk2(NPTS) , wk3(NPTS) , wk4(NPTS)
  real(kind=8) :: tn(NPTS) , O(NPTS) , O2(NPTS) , N2(NPTS) , HYD(NPTS)
  real(kind=8) :: ti_oplus_1d(NPTS)
  real(kind=8) :: ti_hplus_1d(NPTS)
  real(kind=8) :: ni_oplus_1d(NPTS)
  real(kind=8) :: ni_hplus_1d(NPTS)
  real(kind=8) :: chemp_1(NPTS)
  real(kind=8) :: beta_1(NPTS)


  DO 100 i = IN , IS

  !--when this routine is called for o+ calculate the rate coefficients
  !--for the reactions  in m3/sec (since they depend only upon the o+ and
  !--neutral temperatures
  !--o+ on n2 and o+ on o2 reactions from
  !--st. maurice and torr j. geophys. res. 83,969,1978.

      w = ti_oplus_1d(i)/300.
      IF ( TI_oplus_1d(i) < 1700. ) THEN
          wk1(i) = 1.533E-18 - 5.920E-19*w + 8.600E-20*w*w
      ELSE
          wk1(i) = 2.730E-18 - 1.155E-18*w + 1.483E-19*w*w
      ENDIF

      wk2(i) = 2.82E-17 - 7.74E-18*w + 1.073E-18*w*w - &
      5.17E-20*w*w*w + 9.65E-22*w*w*w*w
      wk3(i) = 2.5E-17*SQRT(TN(i))
      wk4(i) = 2.5E-17/1.125*SQRT(ti_hplus_1d(i))


          CHEmp_1(i) = wk4(i)*O(i)*ni_hplus_1d(i)
          BETa_1(i) = wk1(i)*N2(i) + wk2(i)*O2(i) + wk3(i)*HYD(i)

  100 ENDDO

      RETURN


end SUBROUTINE ML__CHEMISTRY_O_plus







 
SUBROUTINE ML__CHEMISTRY_H_plus(IN,IS,TN,O,HYD,ti_hplus_1d,ni_oplus_1d,Chemp_2,beta_2)

!***********************************************************************
! routine to calculate chemical production rates and loss coefficients
! these are stored in chemp and beta respectively
!***********************************************************************

  IMPLICIT NONE
  INTEGER :: i , IN , IS

  real(kind=8) :: wk3(NPTS) , wk4(NPTS)
  real(kind=8) :: tn(NPTS)
  real(kind=8) :: O(NPTS)
  real(kind=8) :: HYD(NPTS)
  real(kind=8) :: ti_hplus_1d(NPTS)
  real(kind=8) :: ni_oplus_1d(NPTS)
  real(kind=8) :: chemp_2(NPTS)
  real(kind=8) :: beta_2(NPTS)


  DO 100 i = IN , IS

  !--when this routine is called for o+ calculate the rate coefficients
  !--for the reactions  in m3/sec (since they depend only upon the o+ and
  !--neutral temperatures
  !--o+ on n2 and o+ on o2 reactions from
  !--st. maurice and torr j. geophys. res. 83,969,1978.

      wk3(i) = 2.5E-17*SQRT(TN(i))
      wk4(i) = 2.5E-17/1.125*SQRT(ti_hplus_1d(i))


  !--h+
          CHEmp_2(i) = wk3(i)*HYD(i)*ni_oplus_1d(i)
          BETa_2(i) = wk4(i)*O(i)

  100 ENDDO

      RETURN




end SUBROUTINE ML__CHEMISTRY_H_plus
















SUBROUTINE ML__DIFFUSION_EQUATION_O_PLUS(J,IN,IS,NUIn,CHEmp_1,BETa_1,PEUvi, &
                                GPAr,DTE_1d,DTI_oplus_1d,dti_hplus_1d,UPAr,ETA,DQ, &
                                TI_oplus_1d,ti_hplus_1d,TE_1d,ni_oplus_1d,ni_hplus_1d, &
                                DNI_oplus_1d,dni_hplus_1d,VI_oplus_1d,vi_hplus_1d,NE_1d, &
                                DVP,O,Mass_Oplus,Mass_Hplus,km,dt,ifailed,i_write_out_tube, &
                                altitude_PZ_km,O_plus_production_fudge_factor)

!***********************************************************************
!          routine to evaluate O+ concentrations and fluxes 
!                 by solving the diffusion equation 
!***********************************************************************

  IMPLICIT NONE
  INTEGER :: i , IN , in1 , IS , is1 , J
  INTEGER :: ifailed
  INTEGER :: i_write_out_tube

  real(kind=8) :: w1 , w2 , w3 , w4 , w5
  real(kind=8) :: w51, w52
  real(kind=8) :: fn , fs
  real(kind=8) :: a(NPTS)
  real(kind=8) :: b(NPTS)
  real(kind=8) :: c(NPTS)
  real(kind=8) :: d(NPTS)
  real(kind=8) :: f(NPTS)
  real(kind=8) :: df(NPTS)
  real(kind=8) :: yp(NPTS)
  real(kind=8) :: yi(NPTS)
  real(kind=8) :: xi(NPTS)
  real(kind=8) :: beta_i(NPTS)
  real(kind=8) :: beta_ij(NPTS)
  real(kind=8) :: nuij(NPTS)
  real(kind=8) :: nuin(NPTS)
  real(kind=8) :: ww2(NPTS)
  real(kind=8) :: chemp_1(NPTS)
  real(kind=8) :: beta_1(NPTS)
  real(kind=8) :: O(NPTS)
  real(kind=8) :: upar(NPTS)
  real(kind=8) :: eta(NPTS)
  real(kind=8) :: dq(NPTS)
  real(kind=8) :: gpar(NPTS)
  real(kind=8) :: dvp(NPTS)
  real(kind=8) :: ni_Oplus_1d(NPTS)
  real(kind=8) :: ni_Hplus_1d(NPTS)
  real(kind=8) :: dni_Oplus_1d(NPTS)
  real(kind=8) :: dni_Hplus_1d(NPTS)
  real(kind=8) :: ne_1d(NPTS)
  real(kind=8) :: ne_1d_saved(NPTS)
  real(kind=8) :: vi_Oplus_1d(NPTS)
  real(kind=8) :: vi_Hplus_1d(NPTS)
  real(kind=8) :: te_1d(NPTS)
  real(kind=8) :: ti_Oplus_1d(NPTS)
  real(kind=8) :: ti_Hplus_1d(NPTS)
  real(kind=8) :: dte_1d(NPTS)
  real(kind=8) :: dti_Oplus_1d(NPTS)
  real(kind=8) :: dti_Hplus_1d(NPTS)
  real(kind=8) :: peuvi(npts,2)
  real(kind=8) :: Mass_Oplus , Mass_Hplus
  real(kind=8) :: km(6)
  real(kind=8) :: dt
  real(kind=8) :: altitude_PZ_km(NPTS)
  real(kind=8) :: O_plus_production_fudge_factor

  in1 = IN + 1
  is1 = IS - 1

    CALL ML__THDIFF(IN,IS,Mass_oplus,Mass_hplus,ni_oplus_1d,ni_hplus_1d, &
                ti_oplus_1d,ti_hplus_1d,beta_i,beta_ij,nuij)

! evaluate the coefficients for the diffusion equation

!--first find hij,hik and hin then xi and yi

  DO 100 i = IN , IS
  !g
  !g  Add a bit to the O+ production rate....the dreaded fudge bit....
  !g
            peuvi(i,j)=peuvi(i,j)+O_plus_production_fudge_factor*o(i)

      ww2(i) = 0.
  100 ENDDO

          DO i = IN , IS
              ww2(i) = ww2(i) + DNI_hplus_1d(i)
          ENDDO

  DO i = IN , IS
      w1 = 1./(nuij(i)+NUIn(i)+1.E-3)
      w2 = 1./NE_1d(i)
      w3 = -GPAr(i) + KM(J)*(TE_1d(i)*ww2(i)*w2+DTE_1d(i)+DTI_oplus_1d(i))
      w4 = KM(J)*(beta_i(i)*DTI_oplus_1d(i)-beta_ij(i)*DTI_hplus_1d(i))
      w5 = nuij(i)*VI_hplus_1d(i) + NUIn(i)*UPAr(i)
      xi(i) = (w3+w4-w5)*w1
      yi(i) = KM(J)*(TI_oplus_1d(i)+TE_1d(i)*NI_oplus_1d(i)*w2)*w1
      NE_1d_saved(i) = NE_1d(i)
      NE_1d(i) = NE_1d(i) - NI_oplus_1d(i)
  ENDDO



  DO i = IN , is1
      yp(i) = (yi(i)+yi(i+1))/DQ(i)
  ENDDO

!--set up the boundary conditions

    NI_oplus_1d(IN) = (PEUvi(IN,J)+CHEmp_1(IN)+NI_oplus_1d(IN)/DT) &
                      /(BETa_1(IN)+1./DT)
    NI_oplus_1d(IS) = (PEUvi(IS,J)+CHEmp_1(IS)+NI_oplus_1d(IS)/DT) &
                     /(BETa_1(IS)+1./DT)

! solve the diffusion equation
!--set up the coefficients of the ML__TRIDIAGONAL system

  DO i = in1 , is1
      d(i) = -PEUvi(i,J) - CHEmp_1(i) - NI_oplus_1d(i)/DT
      w1 = ETA(i)*ETA(i)/(DQ(i-1)+DQ(i))
      a(i) = w1*(yp(i-1)-xi(i-1)/ETA(i-1))
      b(i) = -BETa_1(i) - DVP(i) - w1*(yp(i)+yp(i-1)) - (1./DT)
      c(i) = w1*(yp(i)+xi(i+1)/ETA(i+1))
  ENDDO

    fn = NI_oplus_1d(IN)
    fs = NI_oplus_1d(IS)
 
  CALL ML__TRIDIAGONAL(a,b,c,d,f,fn,fs,IN,IS)

  if (i_write_out_tube.eq.1) then
    do i = in , is
!     write(6,5466) i, altitude_PZ_km(i),a(i),b(i),c(i),d(i),f(i)
      write(6,5466) i, altitude_PZ_km(i),BETa_1(i),DVP(i),yp(i),f(i)
    enddo
  5466 format(i5,f10.0,5e12.4)
  endif


  ifailed=0

! do i = in1 , is1
!     IF ( f(i) <= 0.0 ) THEN
!      write(6,*) 'ifailpoint ',i,altitude_PZ_km(i)
!     ENDIF
! enddo

  DO 600 i = in1 , is1
      IF ( f(i) <= 0.0 ) THEN
          ifailed=1
          !stop
          RETURN
      ENDIF
  600 ENDDO

    DO i = in1 , is1
        NI_oplus_1d(i) = f(i)
    ENDDO

  CALL ML__DF_BY_DS(ni_oplus_1d,dni_oplus_1d,IN,IS,ETA,DQ)

    DO i = IN , IS
        NE_1d(i) = NE_1d(i) + NI_oplus_1d(i)
    ENDDO

! find the velocity from the momentum equation

    DO i = IN , is1
        VI_oplus_1d(i) = -xi(i) - ETA(i)*yi(i) &
        *(NI_oplus_1d(i+1)/NI_oplus_1d(i)-1.)/DQ(i)
    ENDDO
    VI_oplus_1d(IS) = VI_oplus_1d(IS-1)



  RETURN



end SUBROUTINE ML__DIFFUSION_EQUATION_O_PLUS

























SUBROUTINE ML__DIFFUSION_EQUATION_H_PLUS(J,IN,IS,NUIn,CHEmp_1,Chemp_2,BETa_1,beta_2,PEUvi, &
                                GPAr,DTE_1d,DTI_oplus_1d,dti_hplus_1d,UPAr,ETA,DQ, &
                                TI_oplus_1d,ti_hplus_1d,TE_1d,ni_oplus_1d,ni_hplus_1d, &
                                DNI_oplus_1d,dni_hplus_1d,VI_oplus_1d,vi_hplus_1d,NE_1d, &
                                DVP,O,Mass_Oplus,Mass_Hplus,km,dt,ifailed,i_write_out_tube, &
                                altitude_PZ_km,O_plus_production_fudge_factor)

!***********************************************************************
!          routine to evaluate H+ concentrations and fluxes 
!                 by solving the diffusion equation 
!***********************************************************************

  IMPLICIT NONE
  INTEGER :: i , IN , in1 , IS , is1 , J
  INTEGER :: ifailed
  INTEGER :: i_write_out_tube

  real(kind=8) :: w1 , w2 , w3 , w4 , w5
  real(kind=8) :: w51, w52
  real(kind=8) :: fn , fs
  real(kind=8) :: a(NPTS)
  real(kind=8) :: b(NPTS)
  real(kind=8) :: c(NPTS)
  real(kind=8) :: d(NPTS)
  real(kind=8) :: f(NPTS)
  real(kind=8) :: df(NPTS)
  real(kind=8) :: yp(NPTS)
  real(kind=8) :: yi(NPTS)
  real(kind=8) :: xi(NPTS)
  real(kind=8) :: beta_i(NPTS)
  real(kind=8) :: beta_ij(NPTS)
  real(kind=8) :: nuij(NPTS)
  real(kind=8) :: nuin(NPTS)
  real(kind=8) :: ww2(NPTS)
  real(kind=8) :: chemp_1(NPTS)
  real(kind=8) :: chemp_2(NPTS)
  real(kind=8) :: beta_1(NPTS)
  real(kind=8) :: beta_2(NPTS)
  real(kind=8) :: O(NPTS)
  real(kind=8) :: upar(NPTS)
  real(kind=8) :: eta(NPTS)
  real(kind=8) :: dq(NPTS)
  real(kind=8) :: gpar(NPTS)
  real(kind=8) :: dvp(NPTS)
  real(kind=8) :: ni_Oplus_1d(NPTS)
  real(kind=8) :: ni_Hplus_1d(NPTS)
  real(kind=8) :: dni_Oplus_1d(NPTS)
  real(kind=8) :: dni_Hplus_1d(NPTS)
  real(kind=8) :: ne_1d(NPTS)
  real(kind=8) :: vi_Oplus_1d(NPTS)
  real(kind=8) :: vi_Hplus_1d(NPTS)
  real(kind=8) :: te_1d(NPTS)
  real(kind=8) :: ti_Oplus_1d(NPTS)
  real(kind=8) :: ti_Hplus_1d(NPTS)
  real(kind=8) :: dte_1d(NPTS)
  real(kind=8) :: dti_Oplus_1d(NPTS)
  real(kind=8) :: dti_Hplus_1d(NPTS)
  real(kind=8) :: peuvi(npts,2)
  real(kind=8) :: Mass_Oplus , Mass_Hplus
  real(kind=8) :: km(6)
  real(kind=8) :: dt
  real(kind=8) :: altitude_PZ_km(NPTS)
  real(kind=8) :: O_plus_production_fudge_factor

  in1 = IN + 1
  is1 = IS - 1

!g
    CALL ML__THDIFF(IN,IS,Mass_hplus,Mass_oplus,ni_hplus_1d,ni_oplus_1d, &
                ti_hplus_1d,ti_oplus_1d,beta_i,beta_ij,nuij)

! evaluate the coefficients for the diffusion equation

!--first find hij,hik and hin then xi and yi

  DO 100 i = IN , IS
      ww2(i) = 0.
  100 ENDDO

          DO i = IN , IS
              ww2(i) = ww2(i) + DNI_oplus_1d(i)
          ENDDO

  DO i = IN , IS
      w1 = 1./(nuij(i)+NUIn(i)+1.E-3)
      w2 = 1./NE_1d(i)
      w3 = -GPAr(i) + KM(J)*(TE_1d(i)*ww2(i)*w2+DTE_1d(i)+DTI_hplus_1d(i))
      w4 = KM(J)*(beta_i(i)*DTI_hplus_1d(i)-beta_ij(i)*DTI_oplus_1d(i))
      w5 = nuij(i)*VI_oplus_1d(i) + NUIn(i)*UPAr(i)
      w51 = nuij(i)*VI_oplus_1d(i)
      w52 = NUIn(i)*UPAr(i)
      xi(i) = (w3+w4-w5)*w1
      yi(i) = KM(J)*(TI_hplus_1d(i)+TE_1d(i)*NI_hplus_1d(i)*w2)*w1
      NE_1d(i) = NE_1d(i) - NI_hplus_1d(i)
  ENDDO



  DO i = IN , is1
      yp(i) = (yi(i)+yi(i+1))/DQ(i)
  ENDDO

!--set up the boundary conditions

    NI_hplus_1d(IN) = (PEUvi(IN,J)+CHEmp_2(IN)+NI_hplus_1d(IN)/DT) &
                      /(BETa_2(IN)+1./DT)
    NI_hplus_1d(IS) = (PEUvi(IS,J)+CHEmp_2(IS)+NI_hplus_1d(IS)/DT) &
                     /(BETa_2(IS)+1./DT)

! solve the diffusion equation
!--set up the coefficients of the ML__TRIDIAGONAL system

  DO i = in1 , is1
      d(i) = -PEUvi(i,J) - CHEmp_2(i) - NI_hplus_1d(i)/DT
      w1 = ETA(i)*ETA(i)/(DQ(i-1)+DQ(i))
      a(i) = w1*(yp(i-1)-xi(i-1)/ETA(i-1))
      b(i) = -BETa_2(i) - DVP(i) - w1*(yp(i)+yp(i-1)) - (1./DT)
      c(i) = w1*(yp(i)+xi(i+1)/ETA(i+1))
  ENDDO

    fn = NI_hplus_1d(IN)
    fs = NI_hplus_1d(IS)


  CALL ML__TRIDIAGONAL(a,b,c,d,f,fn,fs,IN,IS)

 !if (i_write_out_tube.eq.1) then
 !  do i = in , is
 !    write(6,5466) i, altitude_PZ_km(i),a(i),b(i),c(i),d(i),f(i)
 !  enddo
 !466 format(i5,f10.0,5e12.4)
 !endif

!cg  don't let the H+ density get lower than about 10.
!cg  - no need and it risks going negative....
 
! do i = in , is
!   if (f(i).lt.1.e-6) f(i) = 1.e-6
! enddo

!cg
!cg


  ifailed=0
  DO 600 i = in1 , is1
      IF ( f(i) <= 0.0 ) THEN
          ifailed=1
  !       stop
          RETURN
      ENDIF
  600 ENDDO

    DO i = in1 , is1
        NI_hplus_1d(i) = f(i)
    ENDDO

  CALL ML__DF_BY_DS(ni_hplus_1d,dni_hplus_1d,IN,IS,ETA,DQ)

    DO i = IN , IS
        NE_1d(i) = NE_1d(i) + NI_hplus_1d(i)
    ENDDO

! find the velocity from the momentum equation

    DO i = IN , is1
        VI_hplus_1d(i) = -xi(i) - ETA(i)*yi(i) &
        *(NI_hplus_1d(i+1)/NI_hplus_1d(i)-1.)/DQ(i)
    ENDDO
    VI_hplus_1d(IS) = VI_hplus_1d(IS-1)



  RETURN



end SUBROUTINE ML__DIFFUSION_EQUATION_H_PLUS



























SUBROUTINE ML__FIELD_LINE_INTEGRALS(in,is,ds,midpoint, &
                                         Apex_d1,Apex_d2,Apex_BE3, &
                                         apex_D,apex_d1d1,apex_d1d2,apex_d2d2,apex_BMAG, &
                                         ni_oplus_1d,ni_hplus_1d,ti_oplus_1d,no_plus,o2_plus, &
                                         n2_plus,n_plus, &
                                         tn,o,o2,n2, &
                                         U_zonal,U_merid,U_vert, &
                                         sigma_ped,sigma_hall, &
                                         sigma_phph_dsi,sigma_lmlm_msi, &
                                         sigma_h,sigma_c, &
                                         Kdmph_dsi,Kdmlm)

  implicit none

! Input Arguments:

  integer(kind=4), intent(in)  ::  in    ! North footpoint of a flux tube
  integer(kind=4), intent(in)  ::  is    ! South footpoint of a flux tube
  integer(kind=4), intent(in)  ::  midpoint
  real(kind=8), intent(in)  ::  ds(npts)            !ds 1D [???

  real(kind=8), intent(in)  ::  Apex_d1(3,npts) !d1 of reference above (3.8)
  real(kind=8), intent(in)  ::  Apex_d2(3,npts) !d2 of reference above (3.9)
  real(kind=8), intent(in)  ::  Apex_BE3( npts) !B_e3 of reference above (= Bmag/D),[nT] (4.13)
  real(kind=8), intent(in)  ::  Apex_D(npts)
  real(kind=8), intent(in)  ::  Apex_d1d1(npts)
  real(kind=8), intent(in)  ::  Apex_d1d2(npts)
  real(kind=8), intent(in)  ::  Apex_d2d2(npts)
  real(kind=8), intent(in)  ::  Apex_BMAG(npts)

  real(kind=8), intent(in)  ::  ni_oplus_1d(npts)  !O+ densities [m-6]
  real(kind=8), intent(in)  ::  ni_hplus_1d(npts)  !H+ densities [m-6]
  real(kind=8), intent(in)  ::  ti_oplus_1d(npts)  !O+ temperatures [K]
  real(kind=8), intent(in)  ::  no_plus(npts)  !NO+ density [m-6]
  real(kind=8), intent(in)  ::  o2_plus(npts)  !O2+ density [m-6]
  real(kind=8), intent(in)  ::  n2_plus(npts)  !N2+ density [m-6]
  real(kind=8), intent(in)  ::  n_plus(npts)  !N+ density [m-6]
  real(kind=8), intent(in)  ::  tn(npts)         !Tn [K]
  real(kind=8), intent(in)  ::  o(npts)          !atomic oxygen density []
  real(kind=8), intent(in)  ::  o2(npts)         !molecular oxygen density []
  real(kind=8), intent(in)  ::  n2(npts)         !molecular nitrogen density []

  real(kind=8), intent(in)  ::  U_zonal(npts)    !neutral wind +ggeast  [m/s]
  real(kind=8), intent(in)  ::  U_merid(npts)    !neutral wind +ggsouth
  real(kind=8), intent(in)  ::  U_vert(npts)     !neutral wind +ggup

! Output Arguments: 

  real(kind=8),    intent(out) ::    sigma_phph_dsi(2)      !(5.13) divided by |sin I_m |
  real(kind=8),    intent(out) ::    sigma_lmlm_msi(2)      !(5.14) multiplied by | sin I_m |

  real(kind=8),    intent(out) ::	    sigma_h(2)      !(5.17)
  real(kind=8),    intent(out) ::	    sigma_c(2)      !(5.18)
  real(kind=8),    intent(out) ::         Kdmph_dsi(2)      !(5.19) divided by |sin I_m |
  real(kind=8),    intent(out) ::	    Kdmlm(2)	  !(5.20) plus or minus ????


!----------------------------Local variables-----------------------------

  integer(kind=4) ::  ipts     !
  real(kind=8)    ::  effective_temp(npts)   !???
  real(kind=8)    ::  ion_neut_cf(npts)  ! collision frequency []
  real(kind=8)    ::  ion_mass_amu(npts) ! ion mass in [AMU]
  integer(kind=4) ::  iout   !???
  real(kind=8)    ::  integral513
  real(kind=8)    ::  integral514
  real(kind=8)    ::  integral517
  real(kind=8)    ::  integral518
  real(kind=8)    ::  integral519
  real(kind=8)    ::  integral520
  integer(kind=4) ::  ihem  !
  integer(kind=4) ::  istart
  integer(kind=4) ::  istop
  integer(kind=4) ::  istep

  real(kind=8)    ::  abs_ds     ! |ds(ipts)|
  real(kind=8)    ::  electron_density  ![m-6]
  real(kind=8)    ::  r_factor
  real(kind=8)    ::  sigma_ped(npts)  !pedersen conductivity [mho/m]
  real(kind=8)    ::  sigma_hall(npts) !hall conductivity
  real(kind=8)    ::  Ue1(npts)        !Un parallel to e1 [m/s]  (5.6)
  real(kind=8)    ::  Ue2(npts)        !Un parallel to e2 [m/s]

! na100504:  sinIm not needed
! real(kind=8)    ::  sinlm    ! sin(lam_m)
! real(kind=8)    ::  clm2     ! cos^2(lam_m)
! real(kind=8)    ::  sinIm    ! sin I_m of reference above (3.7)
! real(kind=8)    ::  abs_sinIm  ! | sin I_m |

  real(kind=8)    ::  electron_charge_Coulombs
  electron_charge_Coulombs=1.6022E-19   ! electronic charge [C]

!---------------------------------------------------------

  do ipts = in,is
    effective_temp(ipts) = (tn(ipts)+ti_oplus_1d(ipts))/2.
    IF ( effective_temp(ipts)<tn(ipts) ) effective_temp(ipts) = tn(ipts)
  enddo 

  iout = 0
! get ion_mass_amu & ion_neut_cf

  CALL ML__IONNEUT_PLAS(o,o2,n2, ni_oplus_1d , no_plus, o2_plus, &
                    effective_temp,ion_neut_cf,ion_mass_amu, &
                    in,is,iout)



!   midpoint = (NPTS+1)/2
! separate South(1) & Northern(2) hemisphere
  ihem_loop: do ihem=1,2

! southern hemisphere
  if (ihem==1) then
      istart=is
      istop=midpoint 
      istep=-1
  ! Northern hemisphere
  else if (ihem==2) then
      istart=in
      istop=midpoint - 1
      istep=+1
  endif


  integral513=0.0
  integral514=0.0
  integral517=0.0
  integral518=0.0
  integral519=0.0
  integral520=0.0


! get integral
  ipts_loop1: do  ipts=istart, istop, istep

! get pedersen & hall conductivities

!  electron_density = ni_oplus_1d(ipts)+ni_hplus_1d(ipts)+ no_plus(ipts)+o2_plus(ipts) + n2_plus(ipts)+n_plus(ipts)

!  ....just use the O+, NO+ and O2+ for the Electron density here (original equation above)........
  electron_density = ni_oplus_1d(ipts) + no_plus(ipts) + o2_plus(ipts)


  if (electron_density.gt.1.e-10) then


  r_factor         = ion_mass_amu(ipts)*ion_neut_cf(ipts)/electron_charge_Coulombs/apex_BMAG(ipts)
  sigma_ped(ipts)  = electron_density*electron_charge_Coulombs/apex_BMAG(ipts)*r_factor/(1+(r_factor*r_factor))
  sigma_hall(ipts) = sigma_ped(ipts)*r_factor

! get neutral wind vectors
! Ue1=d1*u    :e1:+east (5.6)
  Ue1(ipts)= apex_d1(1,ipts)*U_zonal(ipts) &
  -apex_d1(2,ipts)*U_merid(ipts) &
  +apex_d1(3,ipts)*U_vert(ipts)
!    Ue1(ipts)= -apex_d1(2,ipts)*50.    ! am 071608 test integrals

! Ue2=d2*u    :e2:+down/equatorward
  Ue2(ipts)= apex_d2(1,ipts)*U_zonal(ipts) &
  -apex_d2(2,ipts)*U_merid(ipts) &
  +apex_d2(3,ipts)*U_vert(ipts)
!    Ue2(ipts)= 0.  ! am 071608 test integrals


  abs_ds=ABS(ds(ipts))

! get integrals

!g
!g  The following integrals all come from page 203 and 204 of the paper.  They are numbered
!g  to match the equations in the paper...
!g  The integral parts are calculated here and then some additional factors are applied below.
!g  Note, however, we ignore the |sin I_m| factor wherever it appears since this is applied
!g  later on within the Dynamo solver...
!g

  integral513= integral513 + sigma_ped(ipts)*apex_d1d1(ipts)*abs_ds/apex_D(ipts)

  integral514= integral514 + sigma_ped(ipts)*apex_d2d2(ipts)*abs_ds/apex_D(ipts)

  integral517= integral517 + sigma_hall(ipts)*abs_ds

  integral518= integral518 + sigma_ped(ipts)*apex_d1d2(ipts)*abs_ds/apex_D(ipts)

  integral519= integral519 + (sigma_ped(ipts)*apex_d1d1(ipts)*Ue2(ipts)/apex_D(ipts) &
  + (sigma_hall(ipts)-sigma_ped(ipts)*apex_d1d2(ipts) &
  /apex_D(ipts))*Ue1(ipts))*abs_ds

  integral520=integral520 + ( (sigma_hall(ipts)+sigma_ped(ipts)*apex_d1d2(ipts) &
  /apex_D(ipts) )*Ue2(ipts) &
  - sigma_ped(ipts)*apex_d2d2(ipts)*Ue1(ipts)/apex_D(ipts) )*abs_ds

!    integral520=integral520 + ( Ue1(ipts)/apex_D(ipts) )*abs_ds! am 071608 test integrals

endif    ! electron density gt 1.e-10 if block (stops electron density = 0.0 producing NaNs)

enddo  ipts_loop1 !: do  ipts=istart, istop



! inputs to the dynamo solver

!g
!g  Integrals 5.13 and 5.14 do not !include the |sin I_m| or 1/|sin I_m| factors respectively
!g  as these are dealt with in the dynamo module....
!g
  sigma_phph_dsi(ihem)=integral513   !(5.13) divided by |sin I_m |
  sigma_lmlm_msi(ihem)=integral514   !(5.14) multiplied by | sin I_m |
!g
!g  Integrals 5.17 and 5.18.....
!g
  sigma_h(ihem)= integral517       !(5.17)
  sigma_c(ihem)= integral518       !(5.18)
!g
!g  integral 5.19 is multiplied by BE3.  However we have not multiplied by |sin I_m| because
!g  this is done within the dynamo module itself (I've said this enough yeh ?)
!g
  Kdmph_dsi(ihem)=  Apex_BE3(istart)*integral519  !(5.19) divided by |sin I_m |
!g
!g  The following is equation 5.20.  The integral is multiplied by BE3.  There is also a minus
!g  sign for the northern hemisphere part (and a plus sign for the southern hemisphere part).
!g  The +- statement is that the upper (lower) sign applies to the northern (southern) magnetic
!g  hemisphere.  This statement is written on page 200 of the paper - just below equation 3.22
!g
  if(ihem == 1) then
      Kdmlm(ihem)= +Apex_BE3(istart)*integral520  !(5.20) plus for southern hemi
  endif
  if(ihem == 2) then
      Kdmlm(ihem)= +Apex_BE3(istart)*integral520  !(5.20) minus for northern hemi
  endif

enddo ihem_loop  !: do ihem=1,2




end SUBROUTINE ML__FIELD_LINE_INTEGRALS




















SUBROUTINE ML__convert_integral_to_dynamo_grid(plasma_integral_1,plasma_integral_2,plasma_integral_3, &
  plasma_integral_4,plasma_integral_5,plasma_integral_6, &
  mag_lat_tubeN, mag_lat_tubeS, &
  dynamo_integral_1,dynamo_integral_2, &
  dynamo_integral_3,dynamo_integral_4, &
  dynamo_integral_5,dynamo_integral_6, &
  dynamo_latitude)
!g
  IMPLICIT NONE
  INTEGER :: N_dynamo_lats,mp,lp, &
  ilat_dynamo,ilon_dynamo,lp_from_every_second_tube

  parameter (N_dynamo_lats = 97)
  REAL(kind=8) :: plasma_integral_1(2,NMP,NLP),plasma_integral_2(2,NMP,NLP), &
  plasma_integral_3(2,NMP,NLP),plasma_integral_4(2,NMP,NLP), &
  plasma_integral_5(2,NMP,NLP),plasma_integral_6(2,NMP,NLP), &
  mag_lat_tubeN(nmp,nlp) , mag_lat_tubeS(nmp,nlp)
!g
  REAL(kind=8) :: dynamo_integral_1(NMP + 1,N_dynamo_lats),dynamo_integral_2(NMP + 1,N_dynamo_lats), &
  dynamo_integral_3(NMP + 1,N_dynamo_lats),dynamo_integral_4(NMP + 1,N_dynamo_lats), &
  dynamo_integral_5(NMP + 1,N_dynamo_lats),dynamo_integral_6(NMP + 1,N_dynamo_lats), &
  dynamo_latitude(N_dynamo_lats)

!g
!g  Loop over our flux tubes....
!g
  do 1000 mp = 1 , NMP

      ilon_dynamo = mp + (nmp / 2)
      if(ilon_dynamo > nmp+1) ilon_dynamo = ilon_dynamo - nmp


  ! write(6,*) 'dynamo convert ',mp,ilon_dynamo

  !g
  !g the dynamo solver needs values from every second flux_tube - because we are now
  !g solving the plasmasphere with twice the number of tubes as used for the dynamo.
  !g The conversion is that we loop over every second tube and then calculate the
  !g corresponding dynamo latitude - using the variable lp_from_every_second_tube
  !g as calulated below.....
  !g

      do 2000 lp = 1 , NLP, 2

          lp_from_every_second_tube = (lp + 1) / 2
      !g
      !g  Southern hemisphere integrals....
      !g
          ilat_dynamo = ((N_dynamo_lats - 1) / 2) - ((NLP + 1) / 2) + lp_from_every_second_tube

          dynamo_integral_1(ilon_dynamo,ilat_dynamo) = plasma_integral_1(1,MP,lp)
          dynamo_integral_2(ilon_dynamo,ilat_dynamo) = plasma_integral_2(1,MP,lp)
          dynamo_integral_3(ilon_dynamo,ilat_dynamo) = plasma_integral_3(1,MP,lp)
          dynamo_integral_4(ilon_dynamo,ilat_dynamo) = plasma_integral_4(1,MP,lp)
          dynamo_integral_5(ilon_dynamo,ilat_dynamo) = plasma_integral_5(1,MP,lp)
          dynamo_integral_6(ilon_dynamo,ilat_dynamo) = plasma_integral_6(1,MP,lp)
      !g
          if(mp == 1) then
              dynamo_latitude(ilat_dynamo) = mag_lat_tubeS(mp,lp)
          endif
      !g
      !g  Northern hemisphere integrals....
      !g
          ilat_dynamo = ((N_dynamo_lats + 1) / 2) + ((NLP + 1) / 2) - lp_from_every_second_tube + 1

          dynamo_integral_1(ilon_dynamo,ilat_dynamo) = plasma_integral_1(2,MP,LP)
          dynamo_integral_2(ilon_dynamo,ilat_dynamo) = plasma_integral_2(2,MP,LP)
          dynamo_integral_3(ilon_dynamo,ilat_dynamo) = plasma_integral_3(2,MP,LP)
          dynamo_integral_4(ilon_dynamo,ilat_dynamo) = plasma_integral_4(2,MP,LP)
          dynamo_integral_5(ilon_dynamo,ilat_dynamo) = plasma_integral_5(2,MP,LP)
          dynamo_integral_6(ilon_dynamo,ilat_dynamo) = plasma_integral_6(2,MP,LP)
      !g
          if(mp == 1) then
              dynamo_latitude(ilat_dynamo) = mag_lat_tubeN(mp,lp)
          endif

      2000 ENDDO
  1000 ENDDO
!g
!g  this covers all dynamo longitudes except ilon_dynamo = 1 which has the same values as 81....
!g
  do ilat_dynamo = 1 , N_dynamo_lats
      dynamo_integral_1(1,ilat_dynamo) = dynamo_integral_1(nmp+1,ilat_dynamo)
      dynamo_integral_2(1,ilat_dynamo) = dynamo_integral_2(nmp+1,ilat_dynamo)
      dynamo_integral_3(1,ilat_dynamo) = dynamo_integral_3(nmp+1,ilat_dynamo)
      dynamo_integral_4(1,ilat_dynamo) = dynamo_integral_4(nmp+1,ilat_dynamo)
      dynamo_integral_5(1,ilat_dynamo) = dynamo_integral_5(nmp+1,ilat_dynamo)
      dynamo_integral_6(1,ilat_dynamo) = dynamo_integral_6(nmp+1,ilat_dynamo)
  enddo



  RETURN



end SUBROUTINE ML__convert_integral_to_dynamo_grid
















SUBROUTINE ML__IONNEUT_PLAS(P1,P2,P3,PI1,PI2,PI3,T,VIN,AMIn, &
  IN,IS,iout)
  IMPLICIT NONE
  REAL(kind=8) :: a , AMIn , amu , b , factor , P1 , P2 , P3 , PI1 , PI2 , &
  sum , summol , T , v1 , v2 , VIN , PI3
  INTEGER :: n , NMAx , iout, in, is

  DIMENSION P1(npts) , P2(npts) , P3(npts) , T(npts) , &
  VIN(npts) , AMIn(npts) , &
  a(3) , b(3) , PI1(npts) , PI2(npts) , PI3(npts)
  REAL(kind=8) :: mi1 , mi2 , mi3
  DATA mi1 , mi2 , mi3/16. , 30. , 32./
  DATA a/3.42E-11 , 6.66E-10 , 6.82E-10/
  DATA b/2.44E-10 , 4.28E-10 , 4.34E-10/
  amu = 1.66E-27
!c  **
!c  **
  factor=1.0
!c  **
!c  **
  DO 100 n = IN , IS
      summol = PI2(n) + PI3(n)
      sum = PI1(n) + PI2(n) + PI3(n)
      v2 = b(1)*P1(n) + b(2)*P2(n) + b(3)*P3(n)
      v1 = a(3)*P3(n) + a(2)*P2(n) + a(1)*P1(n)*factor*SQRT(T(n)) &
      *(1.08-0.139*log10(T(n))+4.51E-03*log10(T(n))**2)
      if(summol < 1.d-90) summol=0.0
      if(v1 < 1.d-90) v1=0.0
      if(v2 < 1.d-90) v2=0.0
  ! if(pi1(n).lt.1.d-90) pi1(n)=0.0
  ! if(iout.eq.1) write(6,*) 'here 5',n
      VIN(n) = (v1*PI1(n)+v2*summol)*1.E-06/sum
      AMIn(n) = (PI1(n)*mi1+PI2(n)*mi2+PI3(n)*mi3)*amu/sum
  100 ENDDO
  RETURN



end SUBROUTINE ML__IONNEUT_PLAS




































SUBROUTINE ML__TRIDIAGONAL(A,B,C,D,F,FN,FS,IN,IS)

!***********************************************************************
! routine to solve a ML__TRIDIAGONAL system of equations with coefficients
! in arrays a,b,c,d, boundary condtions in fn and fs and the solution
! in array f
!***********************************************************************

  IMPLICIT NONE
  REAL(kind=8) :: &
  A , ac , arbd , B , bd , C , crbd , D , dd , ddrbd , F , FN , &
  FS , rbd
  INTEGER :: i , IN , in1 , in2 , IS , is1 , is2 

  DIMENSION A(NPTS) , B(NPTS) , C(NPTS) , D(NPTS) , F(NPTS) , &
  bd(NPTS) , dd(NPTS)
  DIMENSION ac(NPTS) , rbd(NPTS) , arbd(NPTS) , ddrbd(NPTS) , &
  crbd(NPTS)

  in1 = IN + 1
  in2 = IN + 2
  is1 = IS - 1
  is2 = IS - 2
  DO 100 i = in2 , is1
      ac(i) = A(i)*C(i-1)
  100 ENDDO

  bd(in1) = B(in1)
  DO 200 i = in2 , is1
      rbd(i-1) = 1./bd(i-1)
      bd(i) = B(i) - ac(i)*rbd(i-1)
  200 ENDDO
  rbd(is1) = 1./bd(is1)

  DO 300 i = in2 , is1
      arbd(i-1) = A(i)*rbd(i-1)
  300 ENDDO

  D(in1) = D(in1) - A(in1)*FN
  dd(in1) = D(in1)
  D(is1) = D(is1) - C(is1)*FS
  DO 400 i = in2 , is1
      dd(i) = D(i) - dd(i-1)*arbd(i-1)
  400 ENDDO

  F(is1) = dd(is1)/bd(is1)
  DO 500 i = is2 , in1 , -1
      ddrbd(i) = dd(i)*rbd(i)
      crbd(i) = C(i)*rbd(i)
  500 ENDDO

  DO 600 i = is2 , in1 , -1
      F(i) = ddrbd(i) - crbd(i)*F(i+1)
  600 ENDDO

  RETURN



end SUBROUTINE ML__TRIDIAGONAL















SUBROUTINE ML__THDIFF(IN,IS,AI,AJ,NI,NJ,TI,TJ,BETai,BETaij,NUIj)

! routine to evaluate thermal diffusion coefficients (three major ions)
! and ion-ion collision frequencies (quegan et al.,1981)
! a vectorized version.
! a evaluates coefficients and frequencies for altitudes from in to is.

  IMPLICIT NONE
  REAL(kind=8) :: AI , AJ , dij , rnuid , rnujd , w , yijk , zijk
  INTEGER :: i , IN , IS 

  REAL(kind=8) :: nuid , nujd
  REAL(kind=8) :: NI(NPTS) , NJ(NPTS) , TI(NPTS) , TJ(NPTS)
  REAL(kind=8) :: BETai(NPTS) , BETaij(NPTS)
  REAL(kind=8) :: NUIj(NPTS)
  REAL(kind=8) :: nuii(NPTS) , nujj(NPTS)
  REAL(kind=8) :: nuji(NPTS)
  REAL(kind=8) :: nuijd(NPTS)
  REAL(kind=8) :: nujid(NPTS)
  REAL(kind=8) :: nuijdd(NPTS)
  REAL(kind=8) :: nujidd(NPTS)
  REAL(kind=8) :: cij(NPTS) , cji(NPTS)
  REAL(kind=8) :: phiij(NPTS) , phiji(NPTS)

  CALL ML__THD(IN,IS,AI,AJ,NI,NJ,TI,TJ,cij,cji,NUIj,nuji,nuii,nujj, &
  nuijdd,nujidd,nuijd,nujid,phiij,phiji)

  DO 100 i = IN , IS
      nuid = nuii(i) + nuijdd(i)
      nujd = nujj(i) + nujidd(i)
      rnuid = 1./nuid
      rnujd = 1./nujd
      yijk = phiij(i)*(1.0-nujid(i)*rnujd)
      zijk = phiij(i)*(1.0-nuijd(i)*rnuid)
      w = 1.547744E4/(1.-nuijd(i)*rnuid*nujid(i)*rnujd)
      BETai(i) = w*TI(i)*yijk*rnuid/AI
      BETaij(i) = w*AI*TJ(i)*zijk/(AJ*AJ*nujd)
      dij = 0.4*(cij(i)*BETai(i)+cji(i)*BETaij(i)*NI(i)/NJ(i))
      NUIj(i) = NUIj(i)*(1.0-dij)
  100 ENDDO

  RETURN



end SUBROUTINE ML__THDIFF






SUBROUTINE ML__fejer_exb_model(SLT,GL,IDAY,F107,VPE)

  IMPLICIT NONE

    REAL(kind=8), INTENT(IN)  :: SLT, GL, F107
    INTEGER,      INTENT(IN)  :: IDAY

    REAL(kind=8), INTENT(OUT) :: VPE

!*******************************************************************************

! ROUTINE TO DETERMINE VERTICAL EXB DRIFT
! (SCHERLIESS, L., AND B.G. FEJER,
! J. GEOPHYS. RES., 104, 6829-6842, 1999)

!*******************************************************************************

! SLT:   SOLAR LOCAL TIME
! GL:    GEOGRAPHIC LONGITUDE (+ EAST)
! IDAY:  DAY OF YEAR
! F107:  F10.7
! VPE:   EQUATORIAL VERTICAL DRIFT

!*******************************************************************************

  INTEGER :: i, kk, ind, k , j
  REAL(kind=8) :: BSPL4, funct, coeff
  REAL(kind=8) :: BSPL4T , BSPL4L

  DIMENSION COEFF(624),FUNCT(6)



  DATA (COEFF(I),I=1,60)/ &
  -10.80592, -9.63722,-11.52666, -0.05716, -0.06288,  0.03564, &
  -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138, &
  2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171, &
  -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656, &
  -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419, &
  -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984, &
  -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062, &
  -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297, &
  1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023, &
  5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166/
  DATA (COEFF(I),I=61,120)/ &
  3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456, &
  7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008, &
  -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507, &
  -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477, &
  -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528, &
  -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075, &
  14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723, &
  12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484, &
  18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249, &
  4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394/
  DATA (COEFF(I),I=121,180)/ &
  14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739, &
  7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626, &
  7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039, &
  10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177, &
  21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715, &
  19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434, &
  26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631, &
  21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903, &
  28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622, &
  22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206/
  DATA (COEFF(I),I=181,240)/ &
  31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322, &
  46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982, &
  13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970, &
  21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589, &
  16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570, &
  18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991, &
  10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775, &
  12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851, &
  -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433, &
  -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385/
  DATA (COEFF(I),I=241,300)/ &
  2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932, &
  3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084, &
  17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516, &
  0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331, &
  15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934, &
  4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308, &
  9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124, &
  13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465, &
  5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222, &
  9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418/
  DATA (COEFF(I),I=301,360)/ &
  9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919, &
  13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801, &
  10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208, &
  10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751, &
  6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705, &
  5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151, &
  29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325, &
  17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398, &
  8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374, &
  -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187/
  DATA (COEFF(I),I=361,420)/ &
  8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156, &
  14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061, &
  14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129, &
  6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993, &
  7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631, &
  -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188, &
  -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214, &
  -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655, &
  2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650, &
  5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376/
  DATA (COEFF(I),I=421,480)/ &
  13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223, &
  -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585, &
  -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463, &
  -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679, &
  -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355, &
  -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593, &
  -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839, &
  -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226, &
  -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004, &
  -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699/
  DATA (COEFF(I),I=481,540)/ &
  -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796, &
  -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959, &
  -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792, &
  -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535, &
  -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143, &
  -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688, &
  -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422, &
  -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978, &
  -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144, &
  -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133/
  DATA (COEFF(I),I=541,600)/ &
  -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949, &
  -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583, &
  -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213, &
  -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053, &
  -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085, &
  -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726, &
  -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529, &
  -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544, &
  -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807, &
  -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214/
  DATA (COEFF(I),I=601,624)/ &
  -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652, &
  -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235, &
  -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214, &
  -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/

  CALL ML__fejer_G(IDAY,F107,FUNCT,GL)


  VPE=0.
  DO I=1,13
      DO J=1,8,1
          KK=8*(I-1)+J
          DO K=1,6
              IND=6*(KK-1)+K

              call ML__FUNCTION_BSPL4T(I,SLT,BSPL4T)
              call ML__FUNCTION_BSPL4L(J,GL,BSPL4L)

              BSPL4=BSPL4T*BSPL4L

              VPE=VPE+BSPL4*FUNCT(K)*COEFF(IND)
          ENDDO
      ENDDO
  ENDDO

  RETURN




end SUBROUTINE ML__fejer_exb_model



SUBROUTINE ML__fejer_G(IDAY,F107,FUNCT,T)

  IMPLICIT NONE

  REAL(kind=8) :: f107, T
  INTEGER      :: iday
  REAL(kind=8) :: FUNCT, flux
  REAL(kind=8) :: t1, sigma, gauss, cflux, a
  INTEGER      :: i, kk

  DIMENSION FUNCT(6)

  IF(F107 < 75.) THEN
      FLUX=75.
  ELSEIF(F107 > 230.) THEN
      FLUX=230.
  ELSE
      FLUX=F107
  ENDIF
  CFLUX=FLUX

  A=0.
  IF(IDAY >= 120 .AND. IDAY <= 240) THEN
      A=170.
      SIGMA=60.
  ENDIF
  IF(IDAY >= 60 .OR. IDAY >= 300.) THEN
      A=170.
      SIGMA=40.
  ENDIF
  IF(FLUX < 95. .AND. A > 0.) THEN
      GAUSS=EXP(-0.5*((T-A)**2)/SIGMA**2)
      CFLUX=GAUSS*95.+(1.-GAUSS)*FLUX
  ENDIF

  DO 11 I=1,6
      FUNCT(I)=0.
  11 ENDDO

  IF(IDAY >= 135 .AND. IDAY <= 230) THEN
      FUNCT(1)=1.
  ENDIF
  IF(IDAY <= 45 .OR. IDAY >= 320) THEN
      FUNCT(2)=1.
  ENDIF
  IF(IDAY >= 75 .AND. IDAY <= 105) THEN
      FUNCT(3)=1.
  ENDIF
  IF(IDAY >= 260 .AND. IDAY <= 290) THEN
      FUNCT(3)=1.
  ENDIF

  IF(IDAY >= 45 .AND. IDAY <= 75) THEN    ! W-E
      FUNCT(2)=1.-(IDAY-45)/30.
      FUNCT(3)=1.-FUNCT(2)
  ENDIF
  IF(IDAY >= 105 .AND. IDAY <= 135) THEN  ! E-S
      FUNCT(3)=1.-(IDAY-105)/30.
      FUNCT(1)=1.-FUNCT(3)
  ENDIF
  IF(IDAY >= 230 .AND. IDAY <= 260) THEN  ! S-E
      FUNCT(1)=1.-(IDAY-230)/30.
      FUNCT(3)=1.-FUNCT(1)
  ENDIF
  IF(IDAY >= 290 .AND. IDAY <= 320) THEN  ! E-W
      FUNCT(3)=1.-(IDAY-290)/30.
      FUNCT(2)=1.-FUNCT(3)
  ENDIF

  FUNCT(4)=(CFLUX-140.)*FUNCT(1)
  FUNCT(5)=(CFLUX-140.)*FUNCT(2)
  FUNCT(6)=(CFLUX-140.)*FUNCT(3)

  RETURN




end SUBROUTINE ML__fejer_G
SUBROUTINE ML__FUNCTION_BSPL4T(I,T1,BSPL4T)

  IMPLICIT NONE

  integer, intent(in) :: i
  REAL(kind=8), intent(in) :: T1
  REAL(kind=8), intent(out) :: bspl4t

  INTEGER :: j,k
  REAL(kind=8) :: T
  REAL(kind=8) :: tt,b

  DIMENSION TT(0:39),B(20,20)

  DATA TT/ 0.00, 2.75, 4.75, 5.50, 6.25, &
  7.25,10.00,14.00,17.25,18.00, &
  18.75,19.75,21.00,24.00,26.75, &
  28.75,29.50,30.25,31.25,34.00, &
  38.00,41.25,42.00,42.75,43.75, &
  45.00,48.00,50.75,52.75,53.50, &
  54.25,55.25,58.00,62.00,65.25, &
  66.00,66.75,67.75,69.00,72.00/


  T=T1
  IF(I >= 0 .AND. T < TT(I)) THEN
      T=T+24.
  ENDIF
  DO J=I,I+4-1
      IF(T >= TT(J) .AND. T < TT(J+1)) THEN
          B(J,1)=1.
      ELSE
          B(J,1)=0.
      ENDIF
  ENDDO
  DO J=2,4
      DO K=I,I+4-J
          B(K,J)=(T-TT(K))/(TT(K+J-1)-TT(K))*B(K,J-1)
          B(K,J)=B(K,J)+(TT(K+J)-T)/(TT(K+J)-TT(K+1))*B(K+1,J-1)
      ENDDO
  ENDDO

  BSPL4T=B(I,4)

  RETURN




end SUBROUTINE ML__FUNCTION_BSPL4T



SUBROUTINE ML__FUNCTION_BSPL4L(I,T1,bspl4l)

  IMPLICIT NONE

  integer, intent(in) :: i
  REAL(kind=8), intent(in) :: T1
  REAL(kind=8), intent(out) :: bspl4l

  INTEGER :: j,k
  REAL(kind=8) :: TL, T, B

  DIMENSION TL(0:24),B(20,20)

  DATA TL/  0, 10,100,190,200,250,280,310, &
  360,370,460,550,560,610,640,670, &
  720,730,820,910,920,970,1000,1030,1080/

  T=T1
  IF(I >= 0 .AND. T < TL(I)) THEN
      T=T+360.
  ENDIF
  DO J=I,I+4-1
      IF(T >= TL(J) .AND. T < TL(J+1)) THEN
          B(J,1)=1.
      ELSE
          B(J,1)=0.
      ENDIF
  ENDDO

  DO J=2,4
      DO K=I,I+4-J
          B(K,J)=(T-TL(K))/(TL(K+J-1)-TL(K))*B(K,J-1)
          B(K,J)=B(K,J)+(TL(K+J)-T)/(TL(K+J)-TL(K+1))*B(K+1,J-1)
      ENDDO
  ENDDO

  BSPL4L=B(I,4)

  RETURN




end SUBROUTINE ML__FUNCTION_BSPL4L



SUBROUTINE ML__Richmond_lowlat_Efield_model(XMLat,XMLon,DAYno,UT,ISEasav,IUTav,POT,VU,VE)
  IMPLICIT NONE

  REAL(kind=8) :: &
  a , ang , cl , cml , ct , cts , dang , DAYno , daynop , fs , &
  ft , fut , hrang , p , pa , pb , POT , q , rad , rb
  REAL(kind=8) :: &
  rbt , rr , rs , rsm , rsrs , sl , sml , sq2 , st , tl , tla , &
  tlp , tua , tup , UT , va , vb , VE , vp , VU
  REAL(kind=8) :: x , xm , XMLat , xmlatp , XMLon , xnms , xns , xz , z , zp
  INTEGER :: i , icpt , imax , ISEasav , IUTav , j , jf , k , kf , kp , &
  l , lend , lf , lp , lst , m , mf , mm
  INTEGER :: mmf , mp , mpf , mpp , n , nf , np
!g
!g  Equatorial Electric Field model.......
!g
! gives quiet-day ionospheric electrostatic pseudo-potential and e x b
! drifts at 300 km for solar minimum conditions.  see richmond et al. (jgr,
! 1980, p. 4658) for definitions of magnetic coordinates, pseudo-potential,
! and drift components.
!***********************************************************************
! 8/2/2 routine has been modified by deleting shortcuts that depend
! on the invalid assumption that variables are saved between
! successive calls to it.  to reactivate the shortcuts it will be
! necessary to use the save statement to retain all needed variables.
! input parameters -
! xmlat, xmlon are magnetic latitude and east longitude in degrees.
! dayno is day number of the year from 1. to 365.24, with 1. being jan. 1.
! ut is universal time in hours.
! iseasav is 0 if no seasonal averaging is desired.
! is 1 for average over nov. - feb.
! is 2 for average over may -aug.
! is 3 for average over mar., apr., sept., oct.
! is 4 for average over entire year.
! if iseasav.ne.0, dayno is ignored.
! iutav is 0 if no ut averaging is desired.
! is 1 for average over all ut at the fixed local time given by
! ut + (xmlon - 69.)/15.

! output parameters -
! pot is the electrostatic pseudo-potential in volts.
! vu is the drift velocity component perpendicular to the geomagnetic field
! in the upward/poleward direction in the magnetic meridian plane, in m/s.
! ve is the drift velocity component in the magnetic eastward direction, in
! m/s.
! the output values are geophysically meaningful only for latitudes between
! about -65 and +65 degrees.  if iseasav or iutav is out of range pot, vu,
! and ve are set to -1/0.
  DIMENSION kf(128) , lf(128) , mf(128) , nf(128) , jf(128) , q(5) , &
  rs(16) , fut(5) , rr(16) , rsrs(16) , p(16) , vp(16) , &
  pa(5,9) , va(5,9) , fs(3) , ft(3,3) , sml(4) , cml(4) , &
  a(128) , pb(9) , vb(9)
  DATA a/ - 70. , -183. , 31. , -112. , 19. , -39. , -2. , 2. , &
  -33. , 2. , 2. , -111. , 46. , -4. , -5. , 7. , 9. , -17. , &
  2. , 9. , -10. , 2. , -9. , 22. , 145. , -57. , -42. , -6. , &
  6. , -5. , -2. , 20. , 16. , 16. , -77. , -18. , 13. , -8. , &
  16. , -52. , -10. , 7. , 2. , 11. , -28. , 2. , -85. , -82. , &
  3. , -281. , -71. , -25. , -57. , -50. , 21. , -10. , 10. , &
  -81. , 24. , 7. , 5. , 30. , 32. , 5. , -5. , 11. , -31. , &
  8. , 10. , 20. , -15. , -42. , 32. , 7. , -19. , 7. , 34. , &
  -11. , -15. , 26. , 21. , 1. , 22. , 12. , -2. , 275. , &
  777. , -318. , -320. , -208. , 47. , 429. , -523. , 8. , &
  -35. , -224. , -450. , -66. , -7. , -8. , -231. , 55. , 6. , &
  -28. , -51. , -81. , 48. , 9. , 2. , -10. , 54. , 16. , &
  112. , 69. , -33. , 120. , -47. , 5. , -19. , -17. , -23. , &
  -40. , -22. , -21. , -7. , -30. , 15. , 3./
!***********************************************************************
! 8/2/2 deactivate shortcut, since it won't necessarily work right if
! variables are not retained with a save statement.
!***********************************************************************
! (through statement 510) set up constant parameters in first call to
! routine.
! set up values of xmlatp, tup, tlp, and daynop which are not equal to xmlat,
! ut, magnetic local time, and dayno, respectively.
  xmlatp = -361.
  IF ( xmlatp == XMLat ) xmlatp = 0.
  tup = -25.
  IF ( tup == UT ) tup = 0.
  tl = UT + (XMLon-69.)/15.
  tlp = -25.
  IF ( tlp == tl ) tlp = 0.
  daynop = -366.
  IF ( daynop == DAYno ) daynop = 0.
! (through statement 100) select only those terms in series for which
! coefficients a are defined to be non-zero.
  i = 0
  DO 100 kp = 1 , 3
      lst = 4 - kp
      lend = 2 + kp
      DO 50 lp = lst , lend
          l = lp - 3
          IF ( MOD(kp+lp,2) == 0 ) THEN
              DO 10 np = 2 , 8
                  n = np - 1
                  IF ( kp == 3 .OR. n <= 6 ) THEN
                      IF ( IABS(l) /= 2 .OR. n <= 5 ) THEN
                          DO 2 mp = 1 , 9
                              m = mp - 5
                              IF ( IABS(m) <= n ) THEN
                                  IF ( IABS(m) <= 3 .OR. IABS(l) /= 2 ) THEN
                                      IF ( MOD(n-m,2) == 0 ) THEN
                                          i = i + 1
                                          kf(i) = kp
                                          lf(i) = lp
                                          nf(i) = np
                                          mf(i) = mp
                                      ENDIF
                                  ENDIF
                              ENDIF
                          2 ENDDO
                      ENDIF
                  ENDIF
              10 ENDDO
          ENDIF
      50 ENDDO
  100 ENDDO
  imax = i
  ft(1,1) = .75*SQRT(6.E0)/3.1415926535898
  ft(1,2) = 2.E0*ft(1,1)
  ft(1,3) = 1.E0
  ft(2,1) = ft(1,1)
  ft(2,2) = -ft(1,2)
  ft(2,3) = 1.E0
  ft(3,1) = ft(2,2)
  ft(3,2) = 0.
  ft(3,3) = 1.E0
  hrang = 3.1415926535898/12.
  dang = 3.1415926535898/182.62
  sq2 = SQRT(2.E0)
  rad = 180./3.1415926535898
  rb = -6.671E6*5.2E-5
! rb is -(earth radius + 3.e5 m) times dipole magnetic field at pole at 300 km.
  DO 200 i = 1 , imax
      mm = IABS(mf(i)-5)
  ! jf gives appropriate index of legendre polynomials as ordered between
  ! statements 530 and 595.
      jf(i) = (2*(nf(i)+7*mm)-(mm-1)**2+4)/4
  200 ENDDO
! (through statement 500) compute rr (defined as r(n,m)*r(n-1,m)) and rsrs
! (defined as r(n-1,m)**2 + r(n-2,m)**2) needed for legendre polynomial
! generating recursion relations, where r(n,m) is defined as
! sqrt(n**2 - m**2)/sqrt(4*n**2 - 1).  ordering is same as for p(n,m).
  j = 0
  DO 300 mp = 1 , 5
      m = mp - 1
      xm = m
      IF ( m /= 0 ) q(mp) = SQRT((2.*xm+1.)/(2.*xm))
      DO 250 np = mp , 8 , 2
          n = np - 1
          xns = n*n
          xnms = (n-1)**2
          j = j + 1
          rs(j) = (xns-xm*xm)/(4.*xns-1.)
          rsm = (xnms-xm*xm)/(4.*xnms-1.)
          rr(j) = SQRT(rs(j)*rsm)
          IF ( np /= mp ) rsrs(j) = rsm + rs(j-1)
      250 ENDDO
  300 ENDDO
  IF ( IUTav >= 0 .AND. IUTav <= 1 ) THEN
      IF ( IABS(ISEasav-2) <= 2 ) THEN
          icpt = 1
      ! (through statement 530) fs(1), fs(2), fs(3) are factors for amplitude of
      ! semiannual, annual, yearly average components, respectively.
      ! if iseasav = 0, compute fs for given day of the year.
          IF ( ISEasav /= 0 ) THEN
          ! if iseasav = 4, use only yearly average component.
              IF ( ISEasav == 4 ) THEN
                  fs(1) = 0.
                  fs(2) = 0.
                  fs(3) = 1.
              ELSE
              ! if iseasav = 1 - 3, compute fs for appropriate seasonal average.
                  DO 305 k = 1 , 3
                      fs(k) = ft(ISEasav,k)
                  305 ENDDO
              ENDIF
              GOTO 320
          ENDIF
      !***********************************************************************
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (dayno.eq.daynop) go to 530
      !***********************************************************************
          daynop = DAYno
          icpt = 1
          ang = (DAYno+9.)*dang
          fs(1) = sq2*COS(2.*ang)
          fs(2) = sq2*COS(ang)
          fs(3) = 1.
      ! if magnetic latitude is same as in previous call to this routine, skip to
      ! statement 596.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (xmlat.eq.xmlatp) go to 596
      !***********************************************************************
          320 icpt = 1
          xmlatp = XMLat
          ct = SIN(XMLat/rad)
          cts = ct*ct
          st = SQRT(1.-cts)
          rbt = rb*SQRT(.25+.75*cts)
      ! (through statement 595) calculate legendre polynomials p(n,m) as well as vp,
      ! defined as (dp(n,m)/d(colatitude))/(rb*ct).
          j = 0
          DO 340 mp = 1 , 5
          ! for mp=1, p=p(n,m).  for mp.gt.1, p=p(n,m)/st.  this difference is so that
          ! program never divides by st.
              j = j + 1
              xm = mp - 1
              mpp = mp + 2
              IF ( mp > 1 ) THEN
                  p(j) = q(mp)*x
                  x = p(j)*st
                  vp(j) = xm*p(j)/rb
                  xz = 2.*st*st/rb
              ELSE
                  x = 1.
                  p(1) = 1.E0
                  xz = 2.*st/rb
                  vp(1) = 0.
              ENDIF
              DO 330 np = mpp , 8 , 2
                  j = j + 1
                  z = 0.
                  zp = 0.
                  IF ( np /= mpp ) THEN
                      z = rr(j-1)*p(j-2)
                      zp = rr(j-1)*vp(j-2)
                  ENDIF
                  p(j) = ((cts-rsrs(j))*p(j-1)-z)/rr(j)
                  vp(j) = ((cts-rsrs(j))*vp(j-1)-zp-xz*p(j-1))/rr(j)
              330 ENDDO
          340 ENDDO
      ! (through statement 600) calculate arrays of fourier coefficients, with lp
      ! indicating harmonic of ut and mpf indicating harmonic of magnetic local
      ! time.  if arrays are same as in previous call to this routine, skip to
      ! statement 601.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (icpt.eq.0) go to 601
      !***********************************************************************
          DO 360 mpf = 1 , 9
              DO 350 lp = 1 , 5
                  pa(lp,mpf) = 0.
                  va(lp,mpf) = 0.
              350 ENDDO
          360 ENDDO
          DO 380 i = 1 , imax
              lp = lf(i)
              mpf = mf(i)
              j = jf(i)
              x = a(i)*fs(kf(i))
              pa(lp,mpf) = pa(lp,mpf) + x*p(j)
              va(lp,mpf) = va(lp,mpf) + x*vp(j)
          380 ENDDO
      ! (through statement 607) calculate fourier coefficients pb and vb at given ut
      ! for harmonics of tl.
      ! if ut is same as in previous call to this routine, skip to statement 603.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (ut.eq.tup) go to 603
      !***********************************************************************
          tup = UT
          icpt = 1
          tua = UT*hrang
          sl = SIN(tua)
          cl = COS(tua)
          fut(3) = 1.
          fut(2) = sq2*cl
          fut(4) = sq2*sl
          fut(1) = cl*fut(2) - sl*fut(4)
          fut(5) = cl*fut(4) + sl*fut(2)
          IF ( icpt /= 0 ) THEN
              DO 390 mpf = 1 , 9
                  pb(mpf) = 0.
                  vb(mpf) = 0.
                  DO 385 lp = 1 , 5
                      IF ( IUTav == 0 .OR. lp == 3 ) THEN
                          pb(mpf) = pb(mpf) + fut(lp)*pa(lp,mpf)
                          vb(mpf) = vb(mpf) + fut(lp)*va(lp,mpf)
                      ENDIF
                  385 ENDDO
              390 ENDDO
          ENDIF
      ! tl is magnetic local time.
          tl = UT + (XMLon-69.)/15.
      ! if tl is same as in previous call to this routine, skip to statement 630.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (tl.eq.tlp) go to 630
      !***********************************************************************
          tlp = tl
      ! (through statement 610) calculate sines and cosines, times sq2, of harmonics
      ! of magnetic local time.
          tla = tl*hrang
          sl = SIN(tla)
          cl = COS(tla)
          sml(1) = sq2*sl
          cml(1) = sq2*cl
          DO 400 m = 2 , 4
              sml(m) = cl*sml(m-1) + sl*cml(m-1)
              cml(m) = cl*cml(m-1) - sl*sml(m-1)
          400 ENDDO
      ! calculate pot, ve, and vu by summing fourier coefficients multiplied by
      ! appropriate sines and cosines.
          POT = 0.
          VE = 0.
          VU = 0.
          DO 420 m = 1 , 4
              mpf = m + 5
              mmf = 5 - m
              xm = m
              POT = POT + pb(mpf)*sml(m) + pb(mmf)*cml(m)
              VE = VE + vb(mpf)*sml(m) + vb(mmf)*cml(m)
              VU = VU + xm*(pb(mpf)*cml(m)-pb(mmf)*sml(m))
          420 ENDDO
          POT = POT*st + pb(5)
          VE = VE + vb(5)
          VU = VU/rbt
          GOTO 500
      ENDIF
  ENDIF
  POT = 999.
  VU = 0.0
  VE = 0.0
  500 RETURN



end SUBROUTINE ML__Richmond_lowlat_Efield_model












SUBROUTINE ML__THD(IN,IS,AS,AT,NS,NT,TS,TT,CST,CTS,NUSt,NUTs,NUSs, &
  NUTt,NUStdd,NUTsdd,NUStd,NUTsd,PHIst,PHIts)

! routine to evaluate various terms for ML__THDIFF

  IMPLICIT NONE
  REAL(KIND=8) :: AS , ast , AT , b1st , b1ts , b2st , b2ts , b3st , b3ts , &
  d1st , d1ts , d4st , d4ts , ras , rat , rtst , rw1 , rw25 , &
  tst , w
  REAL(KIND=8) :: w1 , w10 , w11 , w12 , w13 , w14 , w15 , w2 , w20 , w21 , &
  w22 , w23 , w24 , w25 , w3 , w4 , w5 , w6 , w7 , w8
  REAL(KIND=8) :: w9 , yst , yts
  INTEGER :: i , IN , IS 

  REAL(KIND=8) :: NS(NPTS) , NT(NPTS) , TS(NPTS) , TT(NPTS) , CST(NPTS) , &
  CTS(NPTS) , NUSt(NPTS) , NUTs(NPTS) , NUSs(NPTS) , NUTt(NPTS) &
  , NUStdd(NPTS) , NUTsdd(NPTS) , NUStd(NPTS) , NUTsd(NPTS) , &
  PHIst(NPTS) , PHIts(NPTS)

  ras = 1./AS
  rat = 1./AT
  w1 = AS + AT
  rw1 = 1./w1
  ast = AS*AT*rw1
  w9 = AS*rat
  w12 = 1.27E-6*SQRT(ast)*ras
  w13 = w12*w9
  w2 = ast*ras
  w3 = ast*rat
  w4 = ast*ast*rw1
  w5 = w4*ras
  w6 = w4*rat
  w7 = w2*w2
  w8 = w3*w3
  w10 = AT*ras
  w11 = ast*rw1
  w14 = 1.27E-6*SQRT(AS*0.5)*ras
  w15 = 1.27E-6*SQRT(AT*0.5)*rat

  DO 100 i = IN , IS

      tst = (AS*TT(i)+AT*TS(i))*rw1
      rtst = 1./tst
      w25 = tst*SQRT(tst)
      rw25 = 1./w25
      NUSt(i) = w12*NT(i)*rw25
      NUTs(i) = w13*NS(i)*rw25

      w20 = (TT(i)-TS(i))*rtst
      w21 = TS(i)*rtst
      w22 = TT(i)*rtst
      w23 = tst/TS(i)
      w24 = tst/TT(i)
      b1st = -0.2*w5*w20*(4.+3.*w9*w22)
      b1ts = 0.2*w6*w20*(4.+3.*w10*w21)
      b2st = 0.2*w5*w20*(4.-3.*w21)
      b2ts = -0.2*w6*w20*(4.-3.*w22)
      b3st = w7*(1.+3.*w8*w20*w20)
      b3ts = w8*(1.+3.*w7*w20*w20)
      yst = w2*(4.-5.*w21-3.*w3*w20)
      yts = w3*(4.-5.*w22+3.*w2*w20)
      CST(i) = 2.5 - w3*(2.5*w22+w10*w23-yst*w20*w23)
      CTS(i) = 2.5 - w2*(2.5*w21+w9*w24+yts*w20*w24)
      d1st = 3.*w8*w22*w22 - b1st - 0.2*b3st + w11*w22*(1.6-1.5*w21)
      d1ts = 3.*w7*w21*w21 - b1ts - 0.2*b3ts + w11*w21*(1.6-1.5*w22)
      d4st = 3.*w7*w21*w21 + b2st - 0.2*b3st - &
      w11*w21*(1.6*w10+1.5*w22)
      d4ts = 3.*w8*w22*w22 + b2ts - 0.2*b3ts - &
      w11*w22*(1.6*w9+1.5*w21)
      NUSs(i) = w14*NS(i)/(TS(i)*SQRT(TS(i)))
      NUTt(i) = w15*NT(i)/(TT(i)*SQRT(TT(i)))
      NUStdd(i) = 1.25*NUSt(i)*(d1st+1.5*w2*w21)
      NUTsdd(i) = 1.25*NUTs(i)*(d1ts+1.5*w3*w22)
      NUStd(i) = 1.25*NUSt(i)*(d4st+1.5*w2*w21)
      NUTsd(i) = 1.25*NUTs(i)*(d4ts+1.5*w3*w22)
      w = 1.211441E-4*rtst*ast
      PHIst(i) = NUSt(i)*w
      PHIts(i) = NUTs(i)*w
  100 ENDDO

  RETURN




end SUBROUTINE ML__THD









SUBROUTINE ML__Calculate_div_vperp(Apex_e1,Apex_e2,Apex_Be3,Apex_grdlbm2, &
  ed1_on_tube,ed2_on_tube,div_vperp_3d,in,is,mp,lp)

  implicit none

  integer, intent(in) :: in(nmp,nlp),is(nmp,nlp), &
  mp,lp
  real(kind=8), intent(in) :: &
  apex_e1(3,npts,nmp), &
  apex_e2(3,npts,nmp), &
  apex_Be3(npts,nmp), &
  Apex_grdlbm2(3,npts,nmp), &
  ed1_on_tube(nmp,nlp), &
  ed2_on_tube(nmp,nlp)
  real(kind=8), intent(out) :: div_vperp_3d(npts,nmp)

  integer  :: i_component, i
  real(kind=8)  :: divvem, ve1, ve2


  do i = in(mp,lp) , is(mp,lp)

  !C BE3 is in nT; need to convert to T
  !g    ve1 =  ed2_on_tube(mp,lp)/(Apex_BE3(i,mp,lp)*1.e-9)
  !g    ve2 = -ed1_on_tube(mp,lp)/(Apex_BE3(i,mp,lp)*1.e-9)
  !g



  !g Apex_BE3s already in Tesla - not nT - it has already been changed in sub.. plasma.
  !g Therefore the above 2 lines become the below 2 lines...
  !g
      ve1 =  ed2_on_tube(mp,lp)/(Apex_BE3(i,mp))
      ve2 = -ed1_on_tube(mp,lp)/(Apex_BE3(i,mp))
  !C divvem is dot product of ExB/B^2 velocity [ve1*E1(i) + ve2*E2(i)]
  !C   with gradient of ALOG(B0**(-2))
      divvem = 0.
      do i_component=1,3
          divvem = divvem + (ve1*Apex_E1(i_component,i,mp) + ve2*Apex_E2(i_component,i,mp)) &
          *Apex_grdlbm2(i_component,i,mp)
      enddo

      div_vperp_3d(i,mp) = divvem

  enddo

end SUBROUTINE ML__Calculate_div_vperp















!-----------------------------------------------------------------------









  FUNCTION K8(TI,TN)
  IMPLICIT NONE
  REAL(KIND=8) :: K8 , TI , TN , temp
  temp = .667*TI + .333*TN
  K8 = 2.82E-11 - 7.74E-12*(temp/300.0) + 1.073E-12*(temp/300.0) &
  **2 - 5.17E-14*(temp/300.0)**3 + 9.65E-16*(temp/300.0)**4
  end FUNCTION K8

  FUNCTION K3(TI,TN)
  IMPLICIT NONE
  REAL(KIND=8) :: K3 , TI , TN , temp
  temp = .6363*TI + .3637*TN
  IF ( temp < 1700.0 ) K3 = 1.533E-12 - 5.92E-13*(temp/300.0) &
  + 8.6E-14*(temp/300.0)**2
  IF ( temp >= 1700.0 ) K3 = 2.73E-12 - 1.155E-12*(temp/300.0) &
  + 1.483E-13*(temp/300.0)**2
  end FUNCTION K3

  FUNCTION K4(TI,TN)
  IMPLICIT NONE
  REAL(KIND=8) :: K4 , TI , TN , temp
  temp = (TI+TN)/2
  IF ( temp < 1500.0 ) K4 = 1.4E-10*(300.0/temp)**.44
  IF ( temp >= 1500.0 ) K4 = 5.2E-11*(temp/300.0)**.2
  end FUNCTION K4


  FUNCTION ERFCHE(X)
  IMPLICIT NONE

  REAL(KIND=8) :: ERFCHE , t , X , z

  z = ABS(X)
  t = 1./(1.+0.5*z)
  ERFCHE = 1. - t*EXP(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*( &
  .09678418+ &
  t*(-.18628806+t*(.27886807+t*(-1.13520398+t*(1.48851587+ &
  t*(-.82215223+t*.17087277)))))))))
  IF ( X < 0. ) ERFCHE = -ERFCHE

  end FUNCTION ERFCHE





















END MODULE IONOSPHERE_PLASMASPHERE
