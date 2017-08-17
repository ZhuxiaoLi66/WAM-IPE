MODULE moduleInterfaceThermo2Iono

USE modSizeFluxTube

IMPLICIT NONE

PRIVATE  ! Set everything to private access, except those
         ! marked as public

PUBLIC :: INTERFACE__thermosphere_to_FIXED_GEO, &
          INTERFACE__FIXED_GEO_to_IONOSPHERE


SAVE

!--------------------------
! Flux tube array sizes
!--------------------------
! GIP grid
!INTEGER, PARAMETER ::  NPTS = 13813
!INTEGER, PARAMETER ::  NMP  = 80
!INTEGER, PARAMETER ::  NLP  = 67

!----------------
! IPE grid      -  in modSizeFluxTube  lrm 20120229
!----------------
!INTEGER, PARAMETER ::  NPTS = 44438
!INTEGER, PARAMETER ::  NMP  = 80
!INTEGER, PARAMETER ::  NLP  = 170



INTEGER, PARAMETER ::  interface_hts = 31
REAL(kind=8) :: fixed_heights_km(interface_hts)


DATA fixed_heights_km /90.,95.,100.,105.,110.,115.,120.,125., &
         150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400., &
         450.,500.,550.,600.,700.,800.,900.,1000., &
         2000.,4000.,6370.,9000./

CONTAINS

!--------------------------------------------------------------------------------------------
! This subroutine interpolates variables from the GT pressure grid to the fixed height grid

SUBROUTINE INTERFACE__thermosphere_to_FIXED_GEO ( &
           GIP_switches, & ! input
           thermospheric_model_name, &  ! input
           ht_dim,  &   ! input
           lat_dim, &   ! input
           lon_dim, &   ! input
           therm_o_density, &   ! input
           therm_o2_density, &   ! input
           therm_n2_density, &    ! input
           !! therm_NO_density, &   ! not needed now
           !! therm_N4S_density, &  ! not needed now
           !! therm_N2D_density, &  ! not needed now
           therm_Tn, &   ! input : temperature from GT
           Vn_east, &   ! input : eastward wind from GT
           Vn_South, &   ! input : southward wind from GT
           Vn_Upward, &   ! input : upward wind from GT
           therm_qion3d, &   ! input
           therm_elx, &   ! input
           therm_ely, &   ! input

           !! AURORA VARIABLES NOT NEEDED FOR IPE YET lrm20110929
           !! therm_qo2p_aurora, therm_qop_aurora, therm_qn2p_aurora, &
           !! therm_qnp_aurora, therm_qtef_aurora, &

           therm_long, &  ! input
           therm_lat, &   ! input
           therm_Z, &     ! input

           !o_density_fixed_ht, o2_density_fixed_ht, n2_density_fixed_ht, &  ! output
           interface_o_density, interface_o2_density, interface_n2_density, &  ! output

           !! NO_density_fixed_ht, &                                           ! output
           !! N4S_density_fixed_ht, &
           !! N2D_density_fixed_ht, &                    ! output

           !Vx_fixed_ht, Vy_fixed_ht, wvz_fixed_ht, & ! wind  ! output
           interface_East, interface_South, interface_Upward, & ! wind  ! output changed names lrm20121126

           !tts_fixed_ht, &    ! output
           interface_Tn, &    ! output  changed names lrm20121126

           qion3d_fixed_ht, elx_fixed_ht, ely_fixed_ht)         ! output

           ! AURORA VARIABLES NOT NEEDED FOR IPE YET lrm20110929
           !qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, & ! output
           !qnp_aurora_fixed_ht, qtef_aurora_fixed_ht)                       ! output

IMPLICIT NONE

!-----------
! Inputs
!-----------
LOGICAL, INTENT(IN) :: GIP_switches(20)
character*10, INTENT(IN)  :: thermospheric_model_name
integer, INTENT(IN)  :: ht_dim , lat_dim , lon_dim

REAL(kind=8), INTENT(IN) :: therm_o_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8), INTENT(IN) :: therm_o2_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8), INTENT(IN) :: therm_n2_density(ht_dim,lon_dim,lat_dim)

! These are used but not set *****
REAL(kind=8) :: therm_NO_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) :: therm_N4S_density(ht_dim,lon_dim,lat_dim)
REAL(kind=8) :: therm_N2D_density(ht_dim,lon_dim,lat_dim)

REAL(kind=8), INTENT(IN) :: therm_Tn(ht_dim,lon_dim,lat_dim)

! wind
REAL(kind=8), INTENT(IN) :: Vn_east(ht_dim,lon_dim,lat_dim)
REAL(kind=8), INTENT(IN) :: Vn_South(ht_dim,lon_dim,lat_dim)
REAL(kind=8), INTENT(IN) :: Vn_Upward(ht_dim,lon_dim,lat_dim)

REAL(kind=8), INTENT(IN) :: therm_qion3d(ht_dim,lon_dim,lat_dim)
REAL(kind=8), INTENT(IN) :: therm_elx(ht_dim,lon_dim,lat_dim)
REAL(kind=8), INTENT(IN) :: therm_ely(ht_dim,lon_dim,lat_dim)

REAL(kind=8), INTENT(IN) :: therm_long(lon_dim)
REAL(kind=8), INTENT(IN) :: therm_lat(lat_dim)
REAL(kind=8), INTENT(IN) :: therm_Z(ht_dim,lon_dim,lat_dim)

! to be passed in the future.........
REAL(kind=8) :: therm_qo2p_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) :: therm_qop_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) :: therm_qn2p_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) :: therm_qnp_aurora(ht_dim,lon_dim,lat_dim)
REAL(kind=8) :: therm_qtef_aurora(ht_dim,lon_dim,lat_dim)

! Outputs ----------------------------
REAL(kind=8), INTENT(OUT) :: interface_o_density(interface_hts,lat_dim,lon_dim)
REAL(kind=8), INTENT(OUT) ::  interface_o2_density(interface_hts,lat_dim,lon_dim)
REAL(kind=8), INTENT(OUT) ::  interface_n2_density(interface_hts,lat_dim,lon_dim)

!REAL(kind=8) :: NO_density_fixed_ht(interface_hts,lat_dim,lon_dim) ! was set but never used lrm20121108
REAL(kind=8) :: N4S_density_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: N2D_density_fixed_ht(interface_hts,lat_dim,lon_dim)


! Wind
REAL(kind=8), INTENT(OUT) :: interface_East(interface_hts,lat_dim,lon_dim)
REAL(kind=8), INTENT(OUT) :: interface_South(interface_hts,lat_dim,lon_dim)
REAL(kind=8), INTENT(OUT) :: interface_Upward(interface_hts,lat_dim,lon_dim)


! temperature
REAL(kind=8), INTENT(OUT) ::  interface_Tn(interface_hts,lat_dim,lon_dim)

REAL(kind=8) :: qion3d_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: elx_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: ely_fixed_ht(interface_hts,lat_dim,lon_dim)

REAL(kind=8) :: qo2p_aurora_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: qop_aurora_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: qn2p_aurora_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: qnp_aurora_fixed_ht(interface_hts,lat_dim,lon_dim)
REAL(kind=8) :: qtef_aurora_fixed_ht(interface_hts,lat_dim,lon_dim)


! Local variables -----------------------------------------



INTEGER :: ilon, ilat, iht
INTEGER :: ilon_therm , ilat_therm , iht_therm 
INTEGER :: ispecial 

LOGICAL :: sw_External_model_provides_NO_N4S_densities
LOGICAL :: sw_input_Auroral_production_is_single_overall_rate

REAL(kind=8) therm_height(ht_dim)


REAL(kind=8) interface_NO_density(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_N4S_density(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_N2D_density(interface_hts,lat_dim,lon_dim)


REAL(kind=8) interface_qion3d(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_elx(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_ely(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_qo2p_aurora(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_qop_aurora(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_qn2p_aurora(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_qnp_aurora(interface_hts,lat_dim,lon_dim)
REAL(kind=8) interface_qtef_aurora(interface_hts,lat_dim,lon_dim)


INTEGER ilon_west_array(lon_dim)
INTEGER ilon_east_array(lon_dim)
INTEGER ilat_north_array(lat_dim)
INTEGER ilat_south_array(lat_dim)
REAL(kind=8) factor_lat_array(lat_dim)
REAL(kind=8) factor_lon_array(lon_dim)
INTEGER ilon_west
INTEGER ilon_east
INTEGER ilat_north
INTEGER ilat_south
REAL(kind=8) factor_lat
REAL(kind=8) factor_lon

INTEGER iht_above_west_north(interface_hts,lat_dim,lon_dim)
INTEGER iht_below_west_north(interface_hts,lat_dim,lon_dim)
INTEGER iht_above_east_north(interface_hts,lat_dim,lon_dim)
INTEGER iht_below_east_north(interface_hts,lat_dim,lon_dim)
INTEGER iht_above_west_south(interface_hts,lat_dim,lon_dim)
INTEGER iht_below_west_south(interface_hts,lat_dim,lon_dim)
INTEGER iht_above_east_south(interface_hts,lat_dim,lon_dim)
INTEGER iht_below_east_south(interface_hts,lat_dim,lon_dim)
INTEGER iht_above
INTEGER iht_below

REAL(kind=8) factor_ht_west_north(interface_hts,lat_dim,lon_dim)
REAL(kind=8) factor_ht_east_north(interface_hts,lat_dim,lon_dim)
REAL(kind=8) factor_ht_west_south(interface_hts,lat_dim,lon_dim)
REAL(kind=8) factor_ht_east_south(interface_hts,lat_dim,lon_dim)

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
!REAL(kind=8) HYD_den_p_u12_log  these are not used lrm20121108
!REAL(kind=8) HYD_den_p_l12_log
!REAL(kind=8) HYD_den_p_12_log
!REAL(kind=8) HEL_den_p_u12_log
!REAL(kind=8) HEL_den_p_l12_log
!REAL(kind=8) HEL_den_p_12_log
REAL(kind=8) NO_den_p_u12_log
REAL(kind=8) NO_den_p_l12_log
REAL(kind=8) NO_den_p_12_log
REAL(kind=8) N4S_den_p_u12_log
REAL(kind=8) N4S_den_p_l12_log
REAL(kind=8) N4S_den_p_12_log
REAL(kind=8) N2D_den_p_u12_log
REAL(kind=8) N2D_den_p_l12_log
REAL(kind=8) N2D_den_p_12_log

! temperature
REAL(kind=8) Tn_p_u12
REAL(kind=8) Tn_p_l12
REAL(kind=8) Tn_p_12

! wind
REAL(kind=8) Vn_E_u12
REAL(kind=8) Vn_E_l12
REAL(kind=8) Vn_E_12

REAL(kind=8) Vn_S_u12
REAL(kind=8) Vn_S_l12
REAL(kind=8) Vn_S_12

REAL(kind=8) Vn_Up_u12
REAL(kind=8) Vn_Up_l12
REAL(kind=8) Vn_Up_12


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
!REAL(kind=8) HYD_den_p_u22_log  these are not used lrm20121108
!REAL(kind=8) HYD_den_p_l22_log
!REAL(kind=8) HYD_den_p_22_log
!REAL(kind=8) HEL_den_p_u22_log
!REAL(kind=8) HEL_den_p_l22_log
!REAL(kind=8) HEL_den_p_22_log
REAL(kind=8) NO_den_p_u22_log
REAL(kind=8) NO_den_p_l22_log
REAL(kind=8) NO_den_p_22_log
REAL(kind=8) N4S_den_p_u22_log
REAL(kind=8) N4S_den_p_l22_log
REAL(kind=8) N4S_den_p_22_log
REAL(kind=8) N2D_den_p_u22_log
REAL(kind=8) N2D_den_p_l22_log
REAL(kind=8) N2D_den_p_22_log

! temperature
REAL(kind=8) Tn_p_u22
REAL(kind=8) Tn_p_l22
REAL(kind=8) Tn_p_22

! wind
REAL(kind=8) Vn_E_u22
REAL(kind=8) Vn_E_l22
REAL(kind=8) Vn_E_22

REAL(kind=8) Vn_S_u22
REAL(kind=8) Vn_S_l22
REAL(kind=8) Vn_S_22

REAL(kind=8) Vn_Up_u22
REAL(kind=8) Vn_Up_l22
REAL(kind=8) Vn_Up_22

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
!REAL(kind=8) HYD_den_p_u11_log  these are not used lrm20121108
!REAL(kind=8) HYD_den_p_l11_log
!REAL(kind=8) HYD_den_p_11_log
!REAL(kind=8) HEL_den_p_u11_log
!REAL(kind=8) HEL_den_p_l11_log
!REAL(kind=8) HEL_den_p_11_log
REAL(kind=8) NO_den_p_u11_log
REAL(kind=8) NO_den_p_l11_log
REAL(kind=8) NO_den_p_11_log
REAL(kind=8) N4S_den_p_u11_log
REAL(kind=8) N4S_den_p_l11_log
REAL(kind=8) N4S_den_p_11_log
REAL(kind=8) N2D_den_p_u11_log
REAL(kind=8) N2D_den_p_l11_log
REAL(kind=8) N2D_den_p_11_log

! temperature
REAL(kind=8) Tn_p_u11
REAL(kind=8) Tn_p_l11
REAL(kind=8) Tn_p_11

! wind
REAL(kind=8) Vn_E_u11
REAL(kind=8) Vn_E_l11
REAL(kind=8) Vn_E_11

REAL(kind=8) Vn_S_u11
REAL(kind=8) Vn_S_l11
REAL(kind=8) Vn_S_11

REAL(kind=8) Vn_Up_u11
REAL(kind=8) Vn_Up_l11
REAL(kind=8) Vn_Up_11

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
!REAL(kind=8) HYD_den_p_u21_log  these are not used lrm20121108
!REAL(kind=8) HYD_den_p_l21_log
!REAL(kind=8) HYD_den_p_21_log
!REAL(kind=8) HEL_den_p_u21_log
!REAL(kind=8) HEL_den_p_l21_log
!REAL(kind=8) HEL_den_p_21_log
REAL(kind=8) NO_den_p_u21_log
REAL(kind=8) NO_den_p_l21_log
REAL(kind=8) NO_den_p_21_log
REAL(kind=8) N4S_den_p_u21_log
REAL(kind=8) N4S_den_p_l21_log
REAL(kind=8) N4S_den_p_21_log
REAL(kind=8) N2D_den_p_u21_log
REAL(kind=8) N2D_den_p_l21_log
REAL(kind=8) N2D_den_p_21_log

! temperature
REAL(kind=8) Tn_p_u21
REAL(kind=8) Tn_p_l21
REAL(kind=8) Tn_p_21

! wind
REAL(kind=8) Vn_E_u21
REAL(kind=8) Vn_E_l21
REAL(kind=8) Vn_E_21

REAL(kind=8) Vn_S_u21
REAL(kind=8) Vn_S_l21
REAL(kind=8) Vn_S_21

REAL(kind=8) Vn_Up_u21
REAL(kind=8) Vn_Up_l21
REAL(kind=8) Vn_Up_21



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
!REAL(kind=8) HYD_den_p_2  these are not used lrm20121108
!REAL(kind=8) HYD_den_p_1
!REAL(kind=8) HEL_den_p_2
!REAL(kind=8) HEL_den_p_1
REAL(kind=8) NO_den_p_2
REAL(kind=8) NO_den_p_1
REAL(kind=8) N4S_den_p_2
REAL(kind=8) N4S_den_p_1
REAL(kind=8) N2D_den_p_2
REAL(kind=8) N2D_den_p_1

! temperature
REAL(kind=8) Tn_p_2 
REAL(kind=8) Tn_p_1

! wind
REAL(kind=8) Vn_E_2 
REAL(kind=8) Vn_E_1
REAL(kind=8) Vn_S_2 
REAL(kind=8) Vn_S_1
REAL(kind=8) Vn_Up_2 
REAL(kind=8) Vn_Up_1

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


REAL(kind=8) interface_long(lon_dim)
REAL(kind=8) interface_lat(lat_dim)
REAL(kind=8) interface_height(interface_hts)

! BEGIN CODE ===========================================================================


sw_External_model_provides_NO_N4S_densities = GIP_switches(5)
sw_input_Auroral_production_is_single_overall_rate = GIP_switches(7)

do ilon = 1 , lon_dim 
  interface_long(ilon) = float(ilon-1)*18.
  if ( thermospheric_model_name == 'TGCM' ) then
     if ( interface_long(ilon) >  180. ) interface_long(ilon) = interface_long(ilon) - 360.
  endif
enddo

do ilat = 1 , lat_dim
  interface_lat(ilat) = float(ilat-46) * 2.
enddo

do iht = 1 , interface_hts
!  interface_height(iht) = (float(iht-1)*5.) + 90.
  interface_height(iht) = fixed_heights_km(iht)
enddo

!print *,'therm_long = ',therm_long
!print *,'therm_lat = ',therm_lat
!print *,'therm_Z(:,10,20) = ',therm_Z(:,10,20)
!print *,'therm_Z(:,11,20) = ',therm_Z(:,11,20)
!print *,'therm_Z(:,9,20) = ',therm_Z(:,9,20)


!print *,'interface_long = ',interface_long
!print *,'interface_lat = ',interface_lat
!print *,'interface_height = ',interface_height

!STOP



!------------------------------------------
! Loop over interface longitudes....
!------------------------------------------
do ilon = 1 , lon_dim

    !longitude interpolation

    ! loop over thermosphere longs to find points east and west....
   ! MAKE MORE EFFICIENT SEARCH ?????
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

   !print *,'factor_lon = ',factor_lon

   factor_lon_array(ilon) = factor_lon
   ilon_west_array(ilon) = ilon_west
   ilon_east_array(ilon) = ilon_east

   !---------------------------------------
   ! Loop over interface latitudes....
   !------------------------------------
   do ilat = 1 , lat_dim

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

     !if(factor_lat < 0.0) factor_lat = 0.0
     !if(factor_lat > 1.0) factor_lat = 1.0  lrm20130829
     if (factor_lat < 0.0) then
        factor_lat = 0.0
     else if (factor_lat > 1.0) then
        factor_lat = 1.0
     endif


     !print *,'factor_lat = ', factor_lat

     factor_lat_array(ilat) = factor_lat
     ilat_north_array(ilat) = ilat_north
     ilat_south_array(ilat) = ilat_south


     !------------------------------------
     ! Loop over interface heights....
     !------------------------------------
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

2500   continue

       if ( iht_above == 1 ) then 
          iht_above = 2
          iht_below = 1
       endif

       factor_ht12 = (interface_height(iht) - therm_height(iht_below)) / &
              (therm_height(iht_above) - therm_height(iht_below))

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

   factor_ht22 = (interface_height(iht) - therm_height(iht_below)) / &
              (therm_height(iht_above) - therm_height(iht_below))

   if( factor_ht22 < 0.0 ) factor_ht22 = 0.0

   !print *,'factor_ht22 = ',factor_ht22

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

   factor_ht11 = (interface_height(iht) - therm_height(iht_below)) / &
              (therm_height(iht_above) - therm_height(iht_below))

   if( factor_ht11 < 0.0 ) factor_ht11 = 0.0

   !print *,'factor_ht11 = ',factor_ht11

   factor_ht_west_south(iht,ilat,ilon) = factor_ht11
   iht_above_west_south(iht,ilat,ilon) = iht_above
   iht_below_west_south(iht,ilat,ilon) = iht_below

   ! REPLACE W/ MORE EFFICIENT SEARCH ????
   do iht_therm = 1 , ht_dim
      therm_height(iht_therm) = therm_Z(iht_therm, ilon_east, ilat_south) / 1000.
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

   !if( factor_ht21 < 0.0 ) factor_ht21 = 0.0  ! moved below to a where statement for the array

   factor_ht_east_south(iht,ilat,ilon) = factor_ht21
   iht_above_east_south(iht,ilat,ilon) = iht_above
   iht_below_east_south(iht,ilat,ilon) = iht_below


enddo  ! heights
enddo  ! lats
enddo  ! lons

WHERE (factor_ht_east_south < 0.) 
       factor_ht_east_south = 0.0
END WHERE


! Loop over interface longitudes....

do ilon = 1 , lon_dim

   factor_lon = factor_lon_array(ilon)
   ilon_west = ilon_west_array(ilon)
   ilon_east = ilon_east_array(ilon)

   ! Loop over interface latitudes....

   do ilat = 1 , lat_dim
     factor_lat = factor_lat_array(ilat)
     ilat_north = ilat_north_array(ilat)
     ilat_south = ilat_south_array(ilat)


     !--------------------------------
     ! Loop over interface heights....
     !--------------------------------
     do iht = 1 , interface_hts

        !-----------------------------------
        ! height interpolation....
        !-----------------------------------

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
   NO_den_p_u12_log = log10(therm_NO_density(iht_above,ilon_west,ilat_north))   ! therm_NO_density is used but never set ***
   NO_den_p_l12_log = log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_12_log = ((NO_den_p_u12_log - NO_den_p_l12_log)*factor_ht12) + NO_den_p_l12_log
   NO_den_p_12_log = 10**(NO_den_p_12_log)
   N4S_den_p_u12_log=log10(therm_N4S_density(iht_above,ilon_west,ilat_north)) ! therm_N4S_density is used but never set ***
   N4S_den_p_l12_log=log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_12_log = ((N4S_den_p_u12_log - N4S_den_p_l12_log)*factor_ht12) + N4S_den_p_l12_log
   N4S_den_p_12_log = 10**(N4S_den_p_12_log)
   N2D_den_p_u12_log=log10(therm_N2D_density(iht_above,ilon_west,ilat_north)) ! therm_N2D_density is used but never set ***
   N2D_den_p_l12_log=log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_12_log = ((N2D_den_p_u12_log - N2D_den_p_l12_log)*factor_ht12) + N2D_den_p_l12_log
   N2D_den_p_12_log = 10**(N2D_den_p_12_log)
endif

   ! Temperature
   Tn_p_u12 = therm_Tn(iht_above,ilon_west,ilat_north)
   Tn_p_l12 = therm_Tn(iht_below,ilon_west,ilat_north)
   Tn_p_12  = ((Tn_p_u12 - Tn_p_l12)*factor_ht12) + Tn_p_l12

   if (Tn_p_12 < 0) then
       print *,'Tn_p_u12, Tn_p_l12, factor_ht12, iht_above,ilon_west,ilat_north = ', &
                Tn_p_u12, Tn_p_l12, factor_ht12, iht_above,ilon_west,ilat_north
   endif 

   ! wind
   Vn_E_u12 = Vn_east(iht_above,ilon_west,ilat_north)
   Vn_E_l12 = Vn_east(iht_below,ilon_west,ilat_north)
   Vn_E_12  = ((Vn_E_u12 - Vn_E_l12)*factor_ht12) + Vn_E_l12

   Vn_S_u12 = Vn_South(iht_above,ilon_west,ilat_north)
   Vn_S_l12 = Vn_South(iht_below,ilon_west,ilat_north)
   Vn_S_12  = ((Vn_S_u12 - Vn_S_l12)*factor_ht12) + Vn_S_l12

   Vn_Up_u12 = Vn_Upward(iht_above,ilon_west,ilat_north)
   Vn_Up_l12 = Vn_Upward(iht_below,ilon_west,ilat_north)
   Vn_Up_12  = ((Vn_Up_u12 - Vn_Up_l12)*factor_ht12) + Vn_Up_l12


if (sw_input_Auroral_production_is_single_overall_rate) then

   qion3d_p_u12 = therm_qion3d(iht_above,ilon_west,ilat_north)
   qion3d_p_l12 = therm_qion3d(iht_below,ilon_west,ilat_north)
   qion3d_p_12  = ((qion3d_p_u12 - qion3d_p_l12)*factor_ht12) + qion3d_p_l12
   if (qion3d_p_12 < 0.0) qion3d_p_12 = 0.0

else

   qo2p_aurora_p_u12 = therm_qo2p_aurora(iht_above,ilon_west,ilat_north) ! therm_qo2p_aurora is used but never set
   qo2p_aurora_p_l12 = therm_qo2p_aurora(iht_below,ilon_west,ilat_north)
   qo2p_aurora_p_12  = ((qo2p_aurora_p_u12 - qo2p_aurora_p_l12)*factor_ht12) + qo2p_aurora_p_l12
   if (qo2p_aurora_p_12 < 0.0) qo2p_aurora_p_12 = 0.0

   qop_aurora_p_u12 = therm_qop_aurora(iht_above,ilon_west,ilat_north) ! therm_qop_aurora is used but never set
   qop_aurora_p_l12 = therm_qop_aurora(iht_below,ilon_west,ilat_north)
   qop_aurora_p_12  = ((qop_aurora_p_u12 - qop_aurora_p_l12)*factor_ht12) + qop_aurora_p_l12
   if (qop_aurora_p_12 < 0.0) qop_aurora_p_12 = 0.0

   qn2p_aurora_p_u12 = therm_qn2p_aurora(iht_above,ilon_west,ilat_north) ! therm_qn2p_aurora is used but never set
   qn2p_aurora_p_l12 = therm_qn2p_aurora(iht_below,ilon_west,ilat_north)
   qn2p_aurora_p_12  = ((qn2p_aurora_p_u12 - qn2p_aurora_p_l12)*factor_ht12) + qn2p_aurora_p_l12
   if (qn2p_aurora_p_12 < 0.0) qn2p_aurora_p_12 = 0.0

   qnp_aurora_p_u12 = therm_qnp_aurora(iht_above,ilon_west,ilat_north) ! therm_qnp_aurora is used but never set
   qnp_aurora_p_l12 = therm_qnp_aurora(iht_below,ilon_west,ilat_north)
   qnp_aurora_p_12  = ((qnp_aurora_p_u12 - qnp_aurora_p_l12)*factor_ht12) + qnp_aurora_p_l12
   if (qnp_aurora_p_12 < 0.0) qnp_aurora_p_12 = 0.0

   qtef_aurora_p_u12 = therm_qtef_aurora(iht_above,ilon_west,ilat_north) ! therm_qtef_aurora is used but never set
   qtef_aurora_p_l12 = therm_qtef_aurora(iht_below,ilon_west,ilat_north)
   qtef_aurora_p_12  = ((qtef_aurora_p_u12 - qtef_aurora_p_l12)*factor_ht12) + qtef_aurora_p_l12
   if (qtef_aurora_p_12 < 0.0) qtef_aurora_p_12 = 0.0

endif ! sw_input_Auroral_production_is_single_overall_rate

   elx_p_u12 = therm_elx(iht_above,ilon_west,ilat_north)
   elx_p_l12 = therm_elx(iht_below,ilon_west,ilat_north)
   elx_p_12  = ((elx_p_u12 - elx_p_l12)*factor_ht12) + elx_p_l12

   ely_p_u12 = therm_ely(iht_above,ilon_west,ilat_north)
   ely_p_l12 = therm_ely(iht_below,ilon_west,ilat_north)
   ely_p_12  = ((ely_p_u12 - ely_p_l12)*factor_ht12) + ely_p_l12


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

   ! temperature
   Tn_p_u22 = therm_Tn(iht_above,ilon_east,ilat_north)
   Tn_p_l22 = therm_Tn(iht_below,ilon_east,ilat_north)
   Tn_p_22  = ((Tn_p_u22 - Tn_p_l22)*factor_ht22) + Tn_p_l22

   if (Tn_p_22 < 0) then
       print *,'Tn_p_u22, Tn_p_l22, factor_ht22, iht_above,ilon_east,ilat_north = ', &
                Tn_p_u22, Tn_p_l22, factor_ht22, iht_above,ilon_east,ilat_north
   endif

   ! wind
   Vn_E_u22 = Vn_east(iht_above,ilon_east,ilat_north)
   Vn_E_l22 = Vn_east(iht_below,ilon_east,ilat_north)
   Vn_E_22  = ((Vn_E_u22 - Vn_E_l22)*factor_ht22) + Vn_E_l22

   Vn_S_u22 = Vn_South(iht_above,ilon_east,ilat_north)
   Vn_S_l22 = Vn_South(iht_below,ilon_east,ilat_north)
   Vn_S_22  = ((Vn_S_u22 - Vn_S_l22)*factor_ht22) + Vn_S_l22

   Vn_Up_u22 = Vn_Upward(iht_above,ilon_east,ilat_north)
   Vn_Up_l22 = Vn_Upward(iht_below,ilon_east,ilat_north)
   Vn_Up_22  = ((Vn_Up_u22 - Vn_Up_l22)*factor_ht22) + Vn_Up_l22


if (sw_input_Auroral_production_is_single_overall_rate) then

   qion3d_p_u22 = therm_qion3d(iht_above,ilon_east,ilat_north)
   qion3d_p_l22 = therm_qion3d(iht_below,ilon_east,ilat_north)
   qion3d_p_22 = ((qion3d_p_u22 - qion3d_p_l22)*factor_ht22) + qion3d_p_l22
   if (qion3d_p_22 < 0.0) qion3d_p_22 = 0.0

else

   qo2p_aurora_p_u22 = therm_qo2p_aurora(iht_above,ilon_east,ilat_north)
   qo2p_aurora_p_l22 = therm_qo2p_aurora(iht_below,ilon_east,ilat_north)
   qo2p_aurora_p_22  = ((qo2p_aurora_p_u22 - qo2p_aurora_p_l22)*factor_ht22) + qo2p_aurora_p_l22
   if (qo2p_aurora_p_22 < 0.0) qo2p_aurora_p_22 = 0.0

   qop_aurora_p_u22 = therm_qop_aurora(iht_above,ilon_east,ilat_north)
   qop_aurora_p_l22 = therm_qop_aurora(iht_below,ilon_east,ilat_north)
   qop_aurora_p_22  = ((qop_aurora_p_u22 - qop_aurora_p_l22)*factor_ht22) + qop_aurora_p_l22
   if (qop_aurora_p_22 < 0.0) qop_aurora_p_22 = 0.0

   qn2p_aurora_p_u22 = therm_qn2p_aurora(iht_above,ilon_east,ilat_north)
   qn2p_aurora_p_l22 = therm_qn2p_aurora(iht_below,ilon_east,ilat_north)
   qn2p_aurora_p_22  = ((qn2p_aurora_p_u22 - qn2p_aurora_p_l22)*factor_ht22) + qn2p_aurora_p_l22
   if (qn2p_aurora_p_22 < 0.0) qn2p_aurora_p_22 = 0.0

   qnp_aurora_p_u22 = therm_qnp_aurora(iht_above,ilon_east,ilat_north)
   qnp_aurora_p_l22 = therm_qnp_aurora(iht_below,ilon_east,ilat_north)
   qnp_aurora_p_22  = ((qnp_aurora_p_u22 - qnp_aurora_p_l22)*factor_ht22) + qnp_aurora_p_l22
   if (qnp_aurora_p_22 < 0.0) qnp_aurora_p_22 = 0.0

   qtef_aurora_p_u22 = therm_qtef_aurora(iht_above,ilon_east,ilat_north)
   qtef_aurora_p_l22 = therm_qtef_aurora(iht_below,ilon_east,ilat_north)
   qtef_aurora_p_22  = ((qtef_aurora_p_u22 - qtef_aurora_p_l22)*factor_ht22) + qtef_aurora_p_l22
   if (qtef_aurora_p_22 < 0.0) qtef_aurora_p_22 = 0.0

endif ! sw_input_Auroral_production_is_single_overall_rate

   elx_p_u22 = therm_elx(iht_above,ilon_east,ilat_north)
   elx_p_l22 = therm_elx(iht_below,ilon_east,ilat_north)
   elx_p_22  = ((elx_p_u22 - elx_p_l22)*factor_ht22) + elx_p_l22

   ely_p_u22 = therm_ely(iht_above,ilon_east,ilat_north)
   ely_p_l22 = therm_ely(iht_below,ilon_east,ilat_north)
   ely_p_22  = ((ely_p_u22 - ely_p_l22)*factor_ht22) + ely_p_l22

   factor_ht11 = factor_ht_west_south(iht,ilat,ilon) 
   iht_above = iht_above_west_south(iht,ilat,ilon) 
   iht_below = iht_below_west_south(iht,ilat,ilon)

   O_den_p_u11_log = log10(therm_o_density(iht_above,ilon_west,ilat_south))
   O_den_p_l11_log = log10(therm_o_density(iht_below,ilon_west,ilat_south))
   O_den_p_11_log = ((O_den_p_u11_log - O_den_p_l11_log)*factor_ht11) + O_den_p_l11_log
   O_den_p_11_log = 10**(O_den_p_11_log)

   O2_den_p_u11_log = log10(therm_o2_density(iht_above,ilon_west,ilat_south))
   O2_den_p_l11_log = log10(therm_o2_density(iht_below,ilon_west,ilat_south))
   O2_den_p_11_log = ((O2_den_p_u11_log - O2_den_p_l11_log)*factor_ht11) + O2_den_p_l11_log
   O2_den_p_11_log = 10**(O2_den_p_11_log)

   N2_den_p_u11_log = log10(therm_n2_density(iht_above,ilon_west,ilat_south))
   N2_den_p_l11_log = log10(therm_n2_density(iht_below,ilon_west,ilat_south))
   N2_den_p_11_log = ((N2_den_p_u11_log - N2_den_p_l11_log)*factor_ht11) + N2_den_p_l11_log
   N2_den_p_11_log = 10**(N2_den_p_11_log)

if (sw_External_model_provides_NO_N4S_densities) then
   NO_den_p_u11_log = log10(therm_NO_density(iht_above,ilon_west,ilat_north))
   NO_den_p_l11_log = log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_11_log = ((NO_den_p_u11_log - NO_den_p_l11_log)*factor_ht11) + NO_den_p_l11_log
   NO_den_p_11_log = 10**(NO_den_p_11_log)

   N4S_den_p_u11_log = log10(therm_N4S_density(iht_above,ilon_west,ilat_north))
   N4S_den_p_l11_log = log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_11_log = ((N4S_den_p_u11_log - N4S_den_p_l11_log)*factor_ht11) + N4S_den_p_l11_log
   N4S_den_p_11_log = 10**(N4S_den_p_11_log)

   N2D_den_p_u11_log = log10(therm_N2D_density(iht_above,ilon_west,ilat_north))
   N2D_den_p_l11_log = log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_11_log = ((N2D_den_p_u11_log - N2D_den_p_l11_log)*factor_ht11) + N2D_den_p_l11_log
   N2D_den_p_11_log = 10**(N2D_den_p_11_log)
endif

! temperature
Tn_p_u11 = therm_Tn(iht_above,ilon_west,ilat_south)
Tn_p_l11 = therm_Tn(iht_below,ilon_west,ilat_south)
Tn_p_11  = ((Tn_p_u11 - Tn_p_l11)*factor_ht11) + Tn_p_l11

! wind
Vn_E_u11 = Vn_east(iht_above,ilon_west,ilat_south)
Vn_E_l11 = Vn_east(iht_below,ilon_west,ilat_south)
Vn_E_11  = ((Vn_E_u11 - Vn_E_l11)*factor_ht11) + Vn_E_l11

Vn_S_u11 = Vn_South(iht_above,ilon_west,ilat_south)
Vn_S_l11 = Vn_South(iht_below,ilon_west,ilat_south)
Vn_S_11  = ((Vn_S_u11 - Vn_S_l11)*factor_ht11) + Vn_S_l11

Vn_Up_u11 = Vn_Upward(iht_above,ilon_west,ilat_south)
Vn_Up_l11 = Vn_Upward(iht_below,ilon_west,ilat_south)
Vn_Up_11  = ((Vn_Up_u11 - Vn_Up_l11)*factor_ht11) + Vn_Up_l11


if (sw_input_Auroral_production_is_single_overall_rate) then

   qion3d_p_u11 = therm_qion3d(iht_above,ilon_west,ilat_south)
   qion3d_p_l11 = therm_qion3d(iht_below,ilon_west,ilat_south)
   qion3d_p_11  = ((qion3d_p_u11 - qion3d_p_l11)*factor_ht11) + qion3d_p_l11
   if (qion3d_p_11 < 0.0) qion3d_p_11 = 0.0

else

   qo2p_aurora_p_u11 = therm_qo2p_aurora(iht_above,ilon_west,ilat_south)
   qo2p_aurora_p_l11 = therm_qo2p_aurora(iht_below,ilon_west,ilat_south)
   qo2p_aurora_p_11 = ((qo2p_aurora_p_u11 - qo2p_aurora_p_l11)*factor_ht11) + qo2p_aurora_p_l11
   if (qo2p_aurora_p_11 < 0.0) qo2p_aurora_p_11 = 0.0

   qop_aurora_p_u11 = therm_qop_aurora(iht_above,ilon_west,ilat_south)
   qop_aurora_p_l11 = therm_qop_aurora(iht_below,ilon_west,ilat_south)
   qop_aurora_p_11 = ((qop_aurora_p_u11 - qop_aurora_p_l11)*factor_ht11) + qop_aurora_p_l11
   if (qop_aurora_p_11 < 0.0) qop_aurora_p_11 = 0.0

   qn2p_aurora_p_u11 = therm_qn2p_aurora(iht_above,ilon_west,ilat_south)
   qn2p_aurora_p_l11 = therm_qn2p_aurora(iht_below,ilon_west,ilat_south)
   qn2p_aurora_p_11 = ((qn2p_aurora_p_u11 - qn2p_aurora_p_l11)*factor_ht11) + qn2p_aurora_p_l11
   if (qn2p_aurora_p_11 < 0.0) qn2p_aurora_p_11 = 0.0

   qnp_aurora_p_u11 = therm_qnp_aurora(iht_above,ilon_west,ilat_south)
   qnp_aurora_p_l11 = therm_qnp_aurora(iht_below,ilon_west,ilat_south)
   qnp_aurora_p_11 = ((qnp_aurora_p_u11 - qnp_aurora_p_l11)*factor_ht11) + qnp_aurora_p_l11
   if (qnp_aurora_p_11 < 0.0) qnp_aurora_p_11 = 0.0

   qtef_aurora_p_u11 = therm_qtef_aurora(iht_above,ilon_west,ilat_south)
   qtef_aurora_p_l11 = therm_qtef_aurora(iht_below,ilon_west,ilat_south)
   qtef_aurora_p_11 = ((qtef_aurora_p_u11 - qtef_aurora_p_l11)*factor_ht11) + qtef_aurora_p_l11
   if (qtef_aurora_p_11 < 0.0) qtef_aurora_p_11 = 0.0

endif ! sw_input_Auroral_production_is_single_overall_rate

elx_p_u11 = therm_elx(iht_above,ilon_west,ilat_south)
elx_p_l11 = therm_elx(iht_below,ilon_west,ilat_south)
elx_p_11  = ((elx_p_u11 - elx_p_l11)*factor_ht11) + elx_p_l11
ely_p_u11 = therm_ely(iht_above,ilon_west,ilat_south)
ely_p_l11 = therm_ely(iht_below,ilon_west,ilat_south)
ely_p_11  = ((ely_p_u11 - ely_p_l11)*factor_ht11) + ely_p_l11


factor_ht21 = factor_ht_east_south(iht,ilat,ilon) 
iht_above  = iht_above_east_south(iht,ilat,ilon)
iht_below  = iht_below_east_south(iht,ilat,ilon)

O_den_p_u21_log = log10(therm_o_density(iht_above,ilon_east,ilat_south))
O_den_p_l21_log = log10(therm_o_density(iht_below,ilon_east,ilat_south))
O_den_p_21_log = ((O_den_p_u21_log - O_den_p_l21_log)*factor_ht21) + O_den_p_l21_log
O_den_p_21_log = 10**(O_den_p_21_log)

O2_den_p_u21_log = log10(therm_o2_density(iht_above,ilon_east,ilat_south))
O2_den_p_l21_log = log10(therm_o2_density(iht_below,ilon_east,ilat_south))
O2_den_p_21_log = ((O2_den_p_u21_log - O2_den_p_l21_log)*factor_ht21) + O2_den_p_l21_log
O2_den_p_21_log = 10**(O2_den_p_21_log)

N2_den_p_u21_log = log10(therm_n2_density(iht_above,ilon_east,ilat_south))
N2_den_p_l21_log = log10(therm_n2_density(iht_below,ilon_east,ilat_south))
N2_den_p_21_log = ((N2_den_p_u21_log - N2_den_p_l21_log)*factor_ht21) + N2_den_p_l21_log
N2_den_p_21_log = 10**(N2_den_p_21_log)


if (sw_External_model_provides_NO_N4S_densities) then

   NO_den_p_u21_log = log10(therm_NO_density(iht_above,ilon_west,ilat_north))
   NO_den_p_l21_log = log10(therm_NO_density(iht_below,ilon_west,ilat_north))
   NO_den_p_21_log = ((NO_den_p_u21_log - NO_den_p_l21_log)*factor_ht21) + NO_den_p_l21_log
   NO_den_p_21_log = 10**(NO_den_p_21_log)

   N4S_den_p_u21_log = log10(therm_N4S_density(iht_above,ilon_west,ilat_north))
   N4S_den_p_l21_log = log10(therm_N4S_density(iht_below,ilon_west,ilat_north))
   N4S_den_p_21_log = ((N4S_den_p_u21_log - N4S_den_p_l21_log)*factor_ht21) + N4S_den_p_l21_log
   N4S_den_p_21_log = 10**(N4S_den_p_21_log)

   N2D_den_p_u21_log = log10(therm_N2D_density(iht_above,ilon_west,ilat_north))
   N2D_den_p_l21_log = log10(therm_N2D_density(iht_below,ilon_west,ilat_north))
   N2D_den_p_21_log = ((N2D_den_p_u21_log - N2D_den_p_l21_log)*factor_ht21) + N2D_den_p_l21_log
   N2D_den_p_21_log = 10**(N2D_den_p_21_log)

endif

! temperature   
Tn_p_u21 = therm_Tn(iht_above,ilon_east,ilat_south) 
Tn_p_l21 = therm_Tn(iht_below,ilon_east,ilat_south)
Tn_p_21  = ((Tn_p_u21 - Tn_p_l21)*factor_ht21) + Tn_p_l21

if (Tn_p_21 < 0) then
   print *,'Tn_p_u21, Tn_p_l21, factor_ht21, iht_above,ilon_east,ilat_south = ', &
            Tn_p_u21, Tn_p_l21, factor_ht21, iht_above,ilon_east,ilat_south
endif

! wind
Vn_E_u21 = Vn_east(iht_above,ilon_east,ilat_south)
Vn_E_l21 = Vn_east(iht_below,ilon_east,ilat_south)
Vn_E_21  = ((Vn_E_u21 - Vn_E_l21)*factor_ht21) + Vn_E_l21

Vn_S_u21 = Vn_South(iht_above,ilon_east,ilat_south)
Vn_S_l21 = Vn_South(iht_below,ilon_east,ilat_south)
Vn_S_21  = ((Vn_S_u21 - Vn_S_l21)*factor_ht21) + Vn_S_l21

Vn_Up_u21 = Vn_Upward(iht_above,ilon_east,ilat_south)
Vn_Up_l21 = Vn_Upward(iht_below,ilon_east,ilat_south)
Vn_Up_21  = ((Vn_Up_u21 - Vn_Up_l21)*factor_ht21) + Vn_Up_l21



if (sw_input_Auroral_production_is_single_overall_rate) then

   qion3d_p_u21 = therm_qion3d(iht_above,ilon_east,ilat_south)
   qion3d_p_l21 = therm_qion3d(iht_below,ilon_east,ilat_south)
   qion3d_p_21 = ((qion3d_p_u21 - qion3d_p_l21)*factor_ht21) + qion3d_p_l21
   if (qion3d_p_21 < 0.0) qion3d_p_21 = 0.0

else

   qo2p_aurora_p_u21 = therm_qo2p_aurora(iht_above,ilon_east,ilat_south)
   qo2p_aurora_p_l21 = therm_qo2p_aurora(iht_below,ilon_east,ilat_south)
   qo2p_aurora_p_21 = ((qo2p_aurora_p_u21 - qo2p_aurora_p_l21)*factor_ht21) + qo2p_aurora_p_l21
   if (qo2p_aurora_p_21 < 0.0) qo2p_aurora_p_21 = 0.0

   qop_aurora_p_u21 = therm_qop_aurora(iht_above,ilon_east,ilat_south)
   qop_aurora_p_l21 = therm_qop_aurora(iht_below,ilon_east,ilat_south)
   qop_aurora_p_21 = ((qop_aurora_p_u21 - qop_aurora_p_l21)*factor_ht21) + qop_aurora_p_l21
   if (qop_aurora_p_21 < 0.0) qop_aurora_p_21 = 0.0

   qn2p_aurora_p_u21 = therm_qn2p_aurora(iht_above,ilon_east,ilat_south)
   qn2p_aurora_p_l21 = therm_qn2p_aurora(iht_below,ilon_east,ilat_south)
   qn2p_aurora_p_21 = ((qn2p_aurora_p_u21 - qn2p_aurora_p_l21)*factor_ht21) + qn2p_aurora_p_l21
   if (qn2p_aurora_p_21 < 0.0) qn2p_aurora_p_21 = 0.0

   qnp_aurora_p_u21 = therm_qnp_aurora(iht_above,ilon_east,ilat_south)
   qnp_aurora_p_l21 = therm_qnp_aurora(iht_below,ilon_east,ilat_south)
   qnp_aurora_p_21 = ((qnp_aurora_p_u21 - qnp_aurora_p_l21)*factor_ht21) + qnp_aurora_p_l21
   if (qnp_aurora_p_21 < 0.0) qnp_aurora_p_21 = 0.0

   qtef_aurora_p_u21 = therm_qtef_aurora(iht_above,ilon_east,ilat_south)
   qtef_aurora_p_l21 = therm_qtef_aurora(iht_below,ilon_east,ilat_south)
   qtef_aurora_p_21 = ((qtef_aurora_p_u21 - qtef_aurora_p_l21)*factor_ht21) + qtef_aurora_p_l21
   if (qtef_aurora_p_21 < 0.0) qtef_aurora_p_21 = 0.0

endif ! sw_input_Auroral_production_is_single_overall_rate


elx_p_u21 = therm_elx(iht_above,ilon_east,ilat_south)
elx_p_l21 = therm_elx(iht_below,ilon_east,ilat_south)
elx_p_21 = ((elx_p_u21 - elx_p_l21)*factor_ht21) + elx_p_l21

ely_p_u21 = therm_ely(iht_above,ilon_east,ilat_south)
ely_p_l21 = therm_ely(iht_below,ilon_east,ilat_south)
ely_p_21 = ((ely_p_u21 - ely_p_l21)*factor_ht21) + ely_p_l21


!--------------------------------
! latitude interpolation....
!--------------------------------
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


! temperature
Tn_p_2 = ((Tn_p_22 - Tn_p_21) * factor_lat) + Tn_p_21
Tn_p_1 = ((Tn_p_12 - Tn_p_11) * factor_lat) + Tn_p_11

!-----------------------------------------------------
! For debugging negative temperature lrm20130520
!------------------------------------------------------
if ( (Tn_p_1 < 0) .or. (Tn_p_2  < 0.)) then
        write (6,*) 'INTERFACE__thermosphere_to_FIXED_GEO : Tn_p_11, Tn_p_12, Tn_p_21, Tn_p_22, factor_lat = ',&
                     Tn_p_11, Tn_p_12, Tn_p_21, Tn_p_22, factor_lat
        write (6,*) 'INTERFACE__thermosphere_to_FIXED_GEO : iht_below, iht_above, ilon_east, ilat_south = ', &
                     iht_below, iht_above, ilon_east, ilat_south
        write (6,*) 'INTERFACE__thermosphere_to_FIXED_GEO : MINVAL(therm_Tn), MINLOC(therm_Tn) = ', &
                     MINVAL(therm_Tn), MINLOC(therm_Tn)
        STOP
endif




! wind
Vn_E_2 = ((Vn_E_22 - Vn_E_21) * factor_lat) + Vn_E_21
Vn_E_1 = ((Vn_E_12 - Vn_E_11) * factor_lat) + Vn_E_11

Vn_S_2 = ((Vn_S_22 - Vn_S_21) * factor_lat) + Vn_S_21
Vn_S_1 = ((Vn_S_12 - Vn_S_11) * factor_lat) + Vn_S_11

Vn_Up_2 = ((Vn_Up_22 - Vn_Up_21) * factor_lat) + Vn_Up_21
Vn_Up_1 = ((Vn_Up_12 - Vn_Up_11) * factor_lat) + Vn_Up_11


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

endif ! sw_input_Auroral_production_is_single_overall_rate


elx_p_2 = ((elx_p_22 - elx_p_21) * factor_lat) + elx_p_21
elx_p_1 = ((elx_p_12 - elx_p_11) * factor_lat) + elx_p_11
ely_p_2 = ((ely_p_22 - ely_p_21) * factor_lat) + ely_p_21
ely_p_1 = ((ely_p_12 - ely_p_11) * factor_lat) + ely_p_11

!-------------------------------
! longitude interpolation....
!-------------------------------
   interface_o_density(iht,ilat,ilon) = ((O_den_p_2 - O_den_p_1) * factor_lon) + O_den_p_1
   interface_o2_density(iht,ilat,ilon) = ((O2_den_p_2 - O2_den_p_1) * factor_lon) + O2_den_p_1
   interface_n2_density(iht,ilat,ilon) = ((N2_den_p_2 - N2_den_p_1) * factor_lon) + N2_den_p_1

   if (interface_o_density(iht,ilat,ilon) .lt. 1.e-20) interface_o_density(iht,ilat,ilon) = 1.e-20
   if (interface_o2_density(iht,ilat,ilon) .lt. 1.e-20) interface_o2_density(iht,ilat,ilon) = 1.e-20
   if (interface_n2_density(iht,ilat,ilon) .lt. 1.e-20) interface_n2_density(iht,ilat,ilon) = 1.e-20


if (sw_External_model_provides_NO_N4S_densities) then

   interface_NO_density(iht,ilat,ilon) = ((NO_den_p_2 - NO_den_p_1) * factor_lon) + NO_den_p_1
   interface_N4S_density(iht,ilat,ilon) = ((N4S_den_p_2 - N4S_den_p_1) * factor_lon) + N4S_den_p_1
   interface_N2D_density(iht,ilat,ilon) = ((N2D_den_p_2 - N2D_den_p_1) * factor_lon) + N2D_den_p_1

   if(interface_NO_density(iht,ilat,ilon) .lt. 1.e-20)  interface_NO_density(iht,ilat,ilon)  = 1.e-20
   if(interface_N4S_density(iht,ilat,ilon) .lt. 1.e-20) interface_N4S_density(iht,ilat,ilon) = 1.e-20
   if(interface_N2D_density(iht,ilat,ilon) .lt. 1.e-20) interface_N2D_density(iht,ilat,ilon) = 1.e-20

endif

   interface_Tn(iht,ilat,ilon) = ((Tn_p_2 - Tn_p_1) * factor_lon) + Tn_p_1

   !-----------------------------------------------------
   ! For debugging negative temperature lrm20130520
   !------------------------------------------------------
   if ( interface_Tn(iht,ilat,ilon) < 0. ) then
        write (6,*) 'INTERFACE__thermosphere_to_FIXED_GEO : iht, ilat, ilon, Tn_p_2, Tn_p_1, factor_lon = ',&
                     iht, ilat, ilon, Tn_p_2, Tn_p_1, factor_lon
        write (6,*) 'INTERFACE__thermosphere_to_FIXED_GEO : MINVAL(therm_Tn), MINLOC(therm_Tn) = ', &
                     MINVAL(therm_Tn), MINLOC(therm_Tn)
        STOP
   endif
   !--------------------------------------------------------

   ! wind
   interface_East(iht,ilat,ilon) = ((Vn_E_2 - Vn_E_1) * factor_lon) + Vn_E_1
   interface_South(iht,ilat,ilon) = ((Vn_S_2 - Vn_S_1) * factor_lon) + Vn_S_1
   interface_Upward(iht,ilat,ilon) = ((Vn_Up_2 - Vn_Up_1) * factor_lon) + Vn_Up_1


if (sw_input_Auroral_production_is_single_overall_rate) then

   interface_qion3d(iht,ilat,ilon) = ((qion3d_p_2 - qion3d_p_1) * factor_lon) + qion3d_p_1
   if (interface_qion3d(iht,ilat,ilon) < 0.0) interface_qion3d(iht,ilat,ilon) = 0.0

else

   interface_qo2p_aurora(iht,ilat,ilon) = ((qo2p_aurora_p_2 - qo2p_aurora_p_1) * factor_lon) + qo2p_aurora_p_1
   if (interface_qo2p_aurora(iht,ilat,ilon) < 0.0) interface_qo2p_aurora(iht,ilat,ilon) = 0.0

   interface_qop_aurora(iht,ilat,ilon) = ((qop_aurora_p_2 - qop_aurora_p_1) * factor_lon) + qop_aurora_p_1
   if (interface_qop_aurora(iht,ilat,ilon) < 0.0) interface_qop_aurora(iht,ilat,ilon) = 0.0

   interface_qn2p_aurora(iht,ilat,ilon) = ((qn2p_aurora_p_2 - qn2p_aurora_p_1) * factor_lon) + qn2p_aurora_p_1
   if (interface_qn2p_aurora(iht,ilat,ilon) < 0.0) interface_qn2p_aurora(iht,ilat,ilon) = 0.0

   interface_qnp_aurora(iht,ilat,ilon) = ((qnp_aurora_p_2 - qnp_aurora_p_1) * factor_lon) + qnp_aurora_p_1
   if (interface_qnp_aurora(iht,ilat,ilon) < 0.0) interface_qnp_aurora(iht,ilat,ilon) = 0.0

   interface_qtef_aurora(iht,ilat,ilon) = ((qtef_aurora_p_2 - qtef_aurora_p_1) * factor_lon) + qtef_aurora_p_1
   if (interface_qtef_aurora(iht,ilat,ilon) < 0.0) interface_qtef_aurora(iht,ilat,ilon) = 0.0

endif ! sw_input_Auroral_production_is_single_overall_rate

   interface_elx(iht,ilat,ilon) = ((elx_p_2 - elx_p_1) * factor_lon)+ elx_p_1
   interface_ely(iht,ilat,ilon) = ((ely_p_2 - ely_p_1) * factor_lon)+ ely_p_1

enddo
enddo
enddo

! maybe able to remove this lrm20121126  - took these out lrm20121126
!o_density_fixed_ht(:,:,:) = interface_o_density(:,:,:) 
!o2_density_fixed_ht(:,:,:) = interface_o2_density(:,:,:)
!n2_density_fixed_ht(:,:,:) = interface_n2_density(:,:,:)

!if (sw_External_model_provides_NO_N4S_densities) then
    ! NO_density_fixed_ht(:,:,:) = interface_NO_density(:,:,:) not used lrm20121108
    !N4S_density_fixed_ht(:,:,:) = interface_N4S_density(:,:,:) not used lrm20121115
    !N2D_density_fixed_ht(:,:,:) = interface_N2D_density(:,:,:) not used lrm20121115
!endif

! wind
!  Why is this done?   does not seem needed **** lrm20121121
!Vx_fixed_ht(:,:,:) = interface_South(:,:,:) commented out, not needed lrm20121126
!Vy_fixed_ht(:,:,:) = interface_East(:,:,:)  commented out, not needed lrm20121126


!print *,'INTERFACE__thermosphere_to_FIXED_GEO : **************************************'
!print *,'INTERFACE__thermosphere_to_FIXED_GEO : wvz_fixed_ht IS SET TO 0 ***********'
!print *,'INTERFACE__thermosphere_to_FIXED_GEO : **************************************'
!wvz_fixed_ht(:,:,:) = 0.0 !   NEED TO FIX THIS ************************** 
!wvz_fixed_ht = interface_Upward  ! lrm20121109  commented out, not needed lrm20121126


! temperature
!tts_fixed_ht(:,:,:) = interface_Tn(:,:,:)  don't need to do this lrm20121126

if (sw_input_Auroral_production_is_single_overall_rate) then
    qion3d_fixed_ht(:,:,:) = interface_qion3d(:,:,:)
endif

!else
    !qo2p_aurora_fixed_ht(:,:,:) = interface_qo2p_aurora(:,:,:)  never used lrm20121115
    !qop_aurora_fixed_ht(:,:,:)  = interface_qop_aurora(:,:,:)  never used lrm20121115
    !qn2p_aurora_fixed_ht(:,:,:) = interface_qn2p_aurora(:,:,:)  never used lrm20121115
    !qnp_aurora_fixed_ht(:,:,:)  = interface_qnp_aurora(:,:,:)  never used lrm20121115
    !qtef_aurora_fixed_ht(:,:,:) = interface_qtef_aurora(:,:,:)  never used lrm20121115
!endif


elx_fixed_ht(:,:,:) = interface_elx(:,:,:)
ely_fixed_ht(:,:,:) = interface_ely(:,:,:)

!------------------------------------------
! for checking if O2 density is negative
!------------------------------------------
!do iht = 1 , interface_hts
!   do ilat = 1 , 91
!      do ilon = 1 , 20
!         if (interface_o2_density(iht,ilat,ilon) < 0.0 ) then
!           write(6,*) 'o2_density negative ', iht, ilat, ilon, &
!                       interface_o2_density(iht,ilat,ilon)
!           stop
!         endif
!      enddo
!   enddo
!enddo

if (MINVAL(interface_o2_density) < 0) then
    print *,'interface__thermosphere_to_FIXED_GEO : NEGATIVE interface_o2_density  ****'
    print *,'interface__thermosphere_to_FIXED_GEO : MINVAL(interface_o2_density), MINLOC(interface_o2_density) = ', &
             MINVAL(interface_o2_density), MINLOC(interface_o2_density)
endif


return

end SUBROUTINE INTERFACE__thermosphere_to_FIXED_GEO


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


! the original name is INTERFACE_TO_PLASMA_FIXED_HT  

SUBROUTINE INTERFACE__FIXED_GEO_to_IONOSPHERE( &
        geo_grid_longitudes_degrees, &
        geo_grid_latitudes_degrees, &
        O_density, &
        O2_density, &
        N2_density, &

        !! NO_density, &  not needed now
        !! N4S_density, &  not needed now
        !! N2D_density, &  not needed now

        V_East , V_South, V_Upward, TTS, &  ! wind, temperature

        !! telec, & not needed now

        IN, IS, &
        ! IWRite2, &  not used lrm20121115
        TN_plasma_input_3d, &
        O_plasma_input_3d, &
        O2_plasma_input_3d, &
        N2_plasma_input_3d, &

        !! NO_plasma_input_3d, N4S_plasma_input_3d, N2D_plasma_input_3d, &  not needed now

        GLAt_3d, &
        GLOnd_3d, &
        PZ_3d, &
        east_plasma_input_3d, south_plasma_input_3d, upward_plasma_input_3d, &

        !! te_plasma_input_3d, & not needed now


        ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &
        ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, &
        ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht, ihu_3d_fixed_ht, &
        sw_1st_call_int_fixed_ht, &
        GIP_switches)




!--------------------------------------------------------
! this calculates a neutral background for the
! plasmasphere code by interpolating values from the
! fixed height interface.....
!--------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: N_heights, N_Latitudes, N_longitudes
  PARAMETER (N_heights=interface_hts, N_Latitudes=91, N_longitudes=20)

  REAL(kind=8) :: geo_grid_longitudes_degrees(N_longitudes)
  REAL(kind=8) :: geo_grid_latitudes_degrees(N_Latitudes)
  REAL(kind=8) :: O_density(N_heights, N_Latitudes, N_longitudes) 
  REAL(kind=8) :: O2_density(N_heights, N_Latitudes, N_longitudes)
  REAL(kind=8) :: N2_density(N_heights, N_Latitudes, N_longitudes) 

  REAL(kind=8) :: V_East (N_heights, N_Latitudes, N_longitudes) ! meridional wind
  REAL(kind=8) :: V_South(N_heights, N_Latitudes, N_longitudes)  ! zonal wind
  REAL(kind=8) :: V_Upward(N_heights, N_Latitudes, N_longitudes)  ! vertical wind
  REAL(kind=8) :: TTS(N_heights, N_Latitudes, N_longitudes)  ! temperature

  ! Not being passed in right now :
  !REAL(kind=8) :: telec(N_heights, N_Latitudes, N_longitudes)  ! not needed yet lrm20121115

  integer :: in(nmp, nlp), is(nmp, nlp)

  !INTEGER :: iwrite2 not used lrm20121115

  REAL(kind=8), INTENT(OUT) :: TN_plasma_input_3d(npts, nmp), O_plasma_input_3d(npts, nmp), &
                  O2_plasma_input_3d(npts, nmp), N2_plasma_input_3d(npts, nmp)

  REAL(kind=8) ::  glat_3d(npts,nmp), glond_3d(npts,nmp)
  REAL(kind=8), INTENT(IN) ::  pz_3d(npts,nmp)

  REAL(kind=8), INTENT(OUT) :: east_plasma_input_3d(npts, nmp), &  ! meridional wind
                  south_plasma_input_3d(npts, nmp), &               ! zonal wind
                  upward_plasma_input_3d(npts, nmp)                  ! vertical wind


  ! Not being passed in right now :
  !REAL(kind=8) :: te_plasma_input_3d(npts, nmp) not needed yet lrm20121115

  INTEGER, INTENT(INOUT) :: ilon1_3d_fixed_ht(npts,nmp), ilon2_3d_fixed_ht(npts,nmp)
  INTEGER, INTENT(INOUT) :: ilat1_3d_fixed_ht(npts,nmp), ilat2_3d_fixed_ht(npts,nmp)

  INTEGER, INTENT(OUT) :: ispecial_3d_fixed_ht(npts,nmp)

  INTEGER, INTENT(OUT) :: ihl_3d_fixed_ht(npts,nmp), ihu_3d_fixed_ht(npts,nmp)

  LOGICAL, INTENT(INOUT)  :: sw_1st_call_int_fixed_ht
  LOGICAL, INTENT(IN) :: GIP_switches(20)


  ! Not being passed in right now :
  REAL(kind=8) :: NO_density(N_heights, N_Latitudes, N_longitudes) 
  REAL(kind=8) :: N4S_density(N_heights, N_Latitudes, N_longitudes) 
  REAL(kind=8) :: N2D_density(N_heights, N_Latitudes, N_longitudes) 

  !REAL(kind=8) :: NO_plasma_input_3d(npts, nmp), &  not used lrm20121115
  !                N4S_plasma_input_3d(npts, nmp), & not used lrm20121115
  !                N2D_plasma_input_3d(npts, nmp)    not used lrm20121115

! Local variables -----------------------------------------------------


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

  REAL(kind=8) :: dnnu21 , dnnu22 , fach,  &
                  faclat , faclon , glond2 

  REAL(kind=8) :: ol11 , ol12 , ol21 , ol22 ,  ou11 , ou12 , ou21 , ou22 , o11 , o12 , o21 , o22 , do1 , do2
  REAL(kind=8) :: ool11 , ool12 , ool21 , ool22 ,  oou11 , oou12 , oou21 , oou22 , oo11 , oo12 , oo21 , oo22  , doo1 , doo2

  REAL(kind=8) :: tnl11 , tnl12 , tnl21 , tnl22 ,  tnu11 , tnu12 , tnu21 , tnu22 , tn11 , tn12 , tn21 , tn22  , tn1 , tn2

  ! meridional wind
  REAL(kind=8) :: eastl11 , eastl12 , eastl21 , eastl22 ,  eastu11 , eastu12 , eastu21 , eastu22 , &
                  east11 , east12 , east21 , east22 , east1 , east2
  ! zonal wind
  REAL(kind=8) :: southl11 , southl12 , southl21 , southl22 ,  southu11 , southu12 , southu21 , southu22 , &
                  south11 , south12 , south21 , south22 , south1 , south2
  ! vertical wind
  REAL(kind=8) :: upwardl11 , upwardl12 , upwardl21 , upwardl22 ,  upwardu11 , upwardu12 , upwardu21 , upwardu22 , &
                  upward11 , upward12 , upward21 , upward22 , upward1 , upward2

 !---------------------------------------------------------------------------------------
 ! This is not used, but Naomi will double check - was used for tiegcm oplus density
 !---------------------------------------------------------------------------------------
 ! REAL(kind=8) :: topl11 , topl12 , topl21 , topl22 ,  topu11 , topu12 , topu21 , &
 !                 topu22 , top11 , top12 , top21 , top22  , top1 , top2
!
 ! REAL(kind=8) :: tnopl11 , tnopl12 , tnopl21 , tnopl22 ,  tnopu11 , tnopu12 , &
 !                 tnopu21 , tnopu22 , tnop11 , tnop12 , tnop21 , tnop22 , tnop1 , tnop2

 ! REAL(kind=8) :: to2pl11 , to2pl12 , to2pl21 , to2pl22 ,  to2pu11 , to2pu12 , &
  !                to2pu21 , to2pu22 , to2p11 , to2p12 , to2p21 , to2p22 , to2p1 , to2p2
  !----------------------------------------------------------------------------------------

  INTEGER :: i , ih , ihl , ihu, ii , ilat , ilat1 , ilat2 , ilon , ilon1 , ilon2 

  !INTEGER :: iprob not used lrm20121115

  INTEGER :: ispecial , IWRite, istop

  !INTEGER ::  l , m , n , itube, ifault, iwrite1  not used lrm20121115

  REAL(kind=8) :: fixed_hts_in_km(N_heights)

  REAL(kind=8) :: pzh(N_heights)

  integer :: mp, lp

  !REAL(kind=8) :: te_dum(npts)  not used lrm20121115

  REAL(kind=8) ::  NO(npts)
  REAL(kind=8) ::  N4S(npts)
  REAL(kind=8) ::  N2D(npts)

  REAL(kind=8) ::  south(NPTS)  ! zonal wind
  REAL(kind=8) ::  east(NPTS)  ! meridional wind
  REAL(kind=8) ::  upward(NPTS)  ! vertical wind

  REAL(kind=8) ::  TN(NPTS) , O(NPTS) , O2(NPTS) , N2(NPTS) , GLAt(NPTS) , &
                   PZ(NPTS) , GLOnd(NPTS)

  REAL(kind=8) :: small_power, small_number


! BEGIN CODE =======================================================================================


small_power = -20.
small_number = 1.d-20

sw_External_model_provides_NO_N4S_densities = GIP_switches(5) 


iwrite = 0
istop = 0
if (istop == 1) stop

do i = 1, N_heights
   fixed_hts_in_km(i) = fixed_heights_km(i)
   pzh(i) = fixed_hts_in_km(i)
enddo


!g -----------------------------------------------------
!g  Big loop over all flux tubes (nmp and nlp).....
!g -----------------------------------------------------

do mp = 1 , nmp   ! longitude sectors
   do lp = 1 , nlp   !  number of tubes for each longitude sector

      !g ---------------------------------------------
      !g  calculate the 1D geographic inputs......
      !l (lat, lon, height of each point we want to 
      !l  interpolate to
      !g ---------------------------------------------
       do i = in(mp,lp), is(mp,lp)  ! points along each tube
          glat(i) = glat_3d(i,mp)
          glond(i) = glond_3d(i,mp)
          pz(i) = pz_3d(i,mp)
       enddo ! i
  

      ! loop over all points on the tube...

      DO 900 i = IN(mp,lp) , IS(mp,lp)

              glond2 = GLOnd(i)
              IF ( glond2 > 360. ) glond2 = glond2 - 360.
              IF ( glond2 < 0. ) glond2 = glond2 + 360.

              if (sw_1st_call_int_fixed_ht) then

                  ispecial = 0

                  DO ilon = 1 , N_longitudes
                      IF ( geo_grid_longitudes_degrees(ilon) > glond2 ) THEN
                          ilon2 = ilon
                          ilon1 = ilon - 1
                          IF ( IWRite == 1 ) WRITE(6,99001) i , glond2 , &
                              geo_grid_longitudes_degrees(ilon2) , geo_grid_longitudes_degrees(ilon1)
                         99001 FORMAT ('this 1 ',i3,3(2x,f10.1))
                             GOTO 250
                      ENDIF
                  ENDDO

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


              !g ------------------------------------------------
              !g thus the required field line point lies
              !g within the square with latitudes ilat2,ilat1
              !g and longitudes ilon2,ilon1....at these four
              !g points there are two heights which are above
              !g and below the point.
              !g ------------------------------------------------

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

                  IF ( IWRite == 1 ) WRITE(6,99001) i , PZ(i) , pzh(ihu) , pzh(ihl)

                  450 IF ( ihu == 1 ) THEN
                          WRITE(6,*) 'Interface__fixed_ : ihu ' , ihu , pzh(ihu) , pzh(ihl)
                          DO ii = IN(mp,lp) , IS(mp,lp)
                              WRITE(6,*) ii , PZ(ii)
                          ENDDO
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

          ! zonal wind on the eight surrounding points......

              southu11 = V_South(ihu,ilat1,ilon1)
              southl11 = V_South(ihl,ilat1,ilon1)
              southu12 = V_South(ihu,ilat1,ilon2)
              southl12 = V_South(ihl,ilat1,ilon2)
              southu21 = V_South(ihu,ilat2,ilon1)
              southl21 = V_South(ihl,ilat2,ilon1)
              southu22 = V_South(ihu,ilat2,ilon2)
              southl22 = V_South(ihl,ilat2,ilon2)

          ! meridional wind on the eight surrounding points......

              eastu11 = V_East (ihu,ilat1,ilon1)
              eastl11 = V_East (ihl,ilat1,ilon1)
              eastu12 = V_East (ihu,ilat1,ilon2)
              eastl12 = V_East (ihl,ilat1,ilon2)
              eastu21 = V_East (ihu,ilat2,ilon1)
              eastl21 = V_East (ihl,ilat2,ilon1)
              eastu22 = V_East (ihu,ilat2,ilon2)
              eastl22 = V_East (ihl,ilat2,ilon2)

          ! vertical wind on the eight surrounding points......

              upwardu11 = V_Upward(ihu,ilat1,ilon1)
              upwardl11 = V_Upward(ihl,ilat1,ilon1)
              upwardu12 = V_Upward(ihu,ilat1,ilon2)
              upwardl12 = V_Upward(ihl,ilat1,ilon2)
              upwardu21 = V_Upward(ihu,ilat2,ilon1)
              upwardl21 = V_Upward(ihl,ilat2,ilon1)
              upwardu22 = V_Upward(ihu,ilat2,ilon2)
              upwardl22 = V_Upward(ihl,ilat2,ilon2)

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

          if (sw_External_model_provides_NO_N4S_densities) then

              print *,'INTERFACE__FIXED_GEO_to_IONOSPHERE : ***********************************************************'
              print *,'INTERFACE__FIXED_GEO_to_IONOSPHERE : ***********************************************************'
              print *,'INTERFACE__FIXED_GEO_to_IONOSPHERE : NO_density, N4S_density, N2D_density USED BUT NOT SET !!!!'
              print *,'INTERFACE__FIXED_GEO_to_IONOSPHERE : ***********************************************************'
              print *,'INTERFACE__FIXED_GEO_to_IONOSPHERE : ***********************************************************'

              ! NO density on the eight surrounding points......

              dNOu11 = log10(NO_density(ihu,ilat1,ilon1)) ! NO_density used by never set ! *******  lrm20121115
              dNOl11 = log10(NO_density(ihl,ilat1,ilon1))
              dNOu12 = log10(NO_density(ihu,ilat1,ilon2))
              dNOl12 = log10(NO_density(ihl,ilat1,ilon2))
              dNOu21 = log10(NO_density(ihu,ilat2,ilon1))
              dNOl21 = log10(NO_density(ihl,ilat2,ilon1))
              dNOu22 = log10(NO_density(ihu,ilat2,ilon2))
              dNOl22 = log10(NO_density(ihl,ilat2,ilon2))

              ! N4S density on the eight surrounding points......

              dN4Su11 = log10(N4S_density(ihu,ilat1,ilon1)) ! N4S_density used by never set ! *******  lrm20121115
              dN4Sl11 = log10(N4S_density(ihl,ilat1,ilon1))
              dN4Su12 = log10(N4S_density(ihu,ilat1,ilon2))
              dN4Sl12 = log10(N4S_density(ihl,ilat1,ilon2))
              dN4Su21 = log10(N4S_density(ihu,ilat2,ilon1))
              dN4Sl21 = log10(N4S_density(ihl,ilat2,ilon1))
              dN4Su22 = log10(N4S_density(ihu,ilat2,ilon2))
              dN4Sl22 = log10(N4S_density(ihl,ilat2,ilon2))

              ! N2D density on the eight surrounding points......

              dN2Du11 = log10(N2D_density(ihu,ilat1,ilon1)) ! N2D_density used by never set ! *******  lrm20121115
              dN2Dl11 = log10(N2D_density(ihl,ilat1,ilon1))
              dN2Du12 = log10(N2D_density(ihu,ilat1,ilon2))
              dN2Dl12 = log10(N2D_density(ihl,ilat1,ilon2))
              dN2Du21 = log10(N2D_density(ihu,ilat2,ilon1))
              dN2Dl21 = log10(N2D_density(ihl,ilat2,ilon1))
              dN2Du22 = log10(N2D_density(ihu,ilat2,ilon2))
              dN2Dl22 = log10(N2D_density(ihl,ilat2,ilon2))

          endif  ! sw_External_model_provides_NO_N4S_densities


          ! now the 8 point interpolation.........

          fach = (PZ(i)-pzh(ihl))/(pzh(ihu)-pzh(ihl))

!          if(.not. sw_1st_call_int_fixed_ht) then
!               write(6,*) 'ZZZZZZZ ',ihl,ihu
!               write(6,*) 'XXXXXXX ',pz(i),pzh(ihl),pzh(ihu)
!          endif

          IF ( fach > 1. ) THEN
                  tn11 = tnu11
                  south11 = southu11
                  east11 = eastu11
                  upward11 = 0.
                  tn12 = tnu12
                  south12 = southu12
                  east12 = eastu12
                  upward12 = 0.
                  tn21 = tnu21
                  south21 = southu21
                  east21 = eastu21
                  upward21 = 0.
                  tn22 = tnu22
                  south22 = southu22
                  east22 = eastu22
                  upward22 = 0.
!                 te11 = teu11+pz(i)-pzh(ihu)
!                 te12 = teu12+pz(i)-pzh(ihu)
!                 te21 = teu21+pz(i)-pzh(ihu)
!                 te22 = teu22+pz(i)-pzh(ihu)
          ELSE
                  tn11 = ((tnu11-tnl11)*fach) + tnl11
                  south11 = ((southu11-southl11)*fach) + southl11
                  east11 = ((eastu11-eastl11)*fach) + eastl11
                  upward11 = ((upwardu11-upwardl11)*fach) + upwardl11
                  tn12 = ((tnu12-tnl12)*fach) + tnl12
                  south12 = ((southu12-southl12)*fach) + southl12
                  east12 = ((eastu12-eastl12)*fach) + eastl12
                  upward12 = ((upwardu12-upwardl12)*fach) + upwardl12
                  tn21 = ((tnu21-tnl21)*fach) + tnl21
                  south21 = ((southu21-southl21)*fach) + southl21
                  east21 = ((eastu21-eastl21)*fach) + eastl21
                  upward21 = ((upwardu21-upwardl21)*fach) + upwardl21
                  tn22 = ((tnu22-tnl22)*fach) + tnl22
                  south22 = ((southu22-southl22)*fach) + southl22
                  east22 = ((eastu22-eastl22)*fach) + eastl22
                  upward22 = ((upwardu22-upwardl22)*fach) + upwardl22
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

          if (sw_External_model_provides_NO_N4S_densities) then

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
          if (sw_External_model_provides_NO_N4S_densities) then
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

              tn2 = ((tn22-tn21)*faclon) + tn21
              tn1 = ((tn12-tn11)*faclon) + tn11
!             te2 = ((te22-te21)*faclon) + te21
!             te1 = ((te12-te11)*faclon) + te11

              south2 = ((south22-south21)*faclon) + south21
              south1 = ((south12-south11)*faclon) + south11
              east2 = ((east22-east21)*faclon) + east21
              east1 = ((east12-east11)*faclon) + east11
              upward2 = ((upward22-upward21)*faclon) + upward21
              upward1 = ((upward12-upward11)*faclon) + upward11
              do2 = ((o22-o21)*faclon) + o21
              do1 = ((o12-o11)*faclon) + o11
              doo2 = ((oo22-oo21)*faclon) + oo21
              doo1 = ((oo12-oo11)*faclon) + oo11
              dnn2 = ((dnn22-dnn21)*faclon) + dnn21
              dnn1 = ((dnn12-dnn11)*faclon) + dnn11

              !top2 = ((top22-top21)*faclon) + top21
              !top1 = ((top12-top11)*faclon) + top11
              !tnop2 = ((tnop22-tnop21)*faclon) + tnop21
              !tnop1 = ((tnop12-tnop11)*faclon) + tnop11
              !to2p2 = ((to2p22-to2p21)*faclon) + to2p21
              !to2p1 = ((to2p12-to2p11)*faclon) + to2p11

          if (sw_External_model_provides_NO_N4S_densities) then
              dNO2 = ((dNO22-dNO21)*faclon) + dNO21
              dNO1 = ((dNO12-dNO11)*faclon) + dNO11
              dN4S2 = ((dN4S22-dN4S21)*faclon) + dN4S21
              dN4S1 = ((dN4S12-dN4S11)*faclon) + dN4S11
              dN2D2 = ((dN2D22-dN2D21)*faclon) + dN2D21
              dN2D1 = ((dN2D12-dN2D11)*faclon) + dN2D11
          endif

              faclat = (GLAt(i)-geo_grid_latitudes_degrees(ilat1)) / &
                       (geo_grid_latitudes_degrees(ilat2)-geo_grid_latitudes_degrees(ilat1))

              TN(i) = ((tn2-tn1)*faclat) + tn1
              ! check for negative temperature 
              if (TN(i) < 0.) then 
                  write(*,*) 'tn1, tn2, faclat = ',tn1, tn2, faclat
              endif

!             Te_dum(i) = ((te2-te1)*faclat) + te1
              !if(te_dum(i) > 5000.) then
              ! write(6,*) 'TE_DUM ',i,mp,lp,te_dum(i)
              !   te_dum(i) = 5000.
              !endif
              south(i) = ((south2-south1)*faclat) + south1
              east(i) = ((east2-east1)*faclat) + east1
              upward(i) = ((upward2-upward1)*faclat) + upward1

              O(i) = ((do2-do1)*faclat) + do1
              if(o(i) < small_number) o(i) = small_number
              O2(i) = ((doo2-doo1)*faclat) + doo1
              if(o2(i) < small_number) o2(i) = small_number
              N2(i) = ((dnn2-dnn1)*faclat) + dnn1
              if (n2(i) < small_number) n2(i) = small_number

          if (sw_External_model_provides_NO_N4S_densities) then
              NO(i) = ((dNO2-dNO1)*faclat) + dNO1
              if(NO(i) < small_number) NO(i)=small_number
              N4S(i) = ((dN4S2-dN4S1)*faclat) + dN4S1
              if(N4S(i) < small_number) N4S(i)=small_number
              N2D(i) = ((dN2D2-dN2D1)*faclat) + dN2D1
              if(N2D(i) < small_number) N2D(i)=small_number
          endif

          !----------------------------------------------------------------------
          !g THe above does not let any of the densities get lower
          !g than 'small_number'.  This can be a problem at low solar activity
          !g at the top of the larger flux-tubes.......
          !g----------------------------------------------------------------------

          900 ENDDO


          do i = in(mp,lp), is(mp,lp)

              TN_plasma_input_3d(i,mp) = tn(i)
              TN_plasma_input_3d(i,mp) = tn(i)
              if (TN_plasma_input_3d(i,mp) .lt. 0.) then
                  write(6,*) 'TN_plasma_input_3d lt 0)',i,mp,lp,tn(i)
                  TN_plasma_input_3d(i,mp) = 50.
              endif

              O_plasma_input_3d(i,mp) = o(i)
              O2_plasma_input_3d(i,mp) = o2(i)
              N2_plasma_input_3d(i,mp) = n2(i)

              if (sw_External_model_provides_NO_N4S_densities) then
                 !NO_plasma_input_3d(i,mp) = NO(i)  ! not used lrm20121115
                 !N4S_plasma_input_3d(i,mp) = N4S(i) ! not used lrm20121115
                 !N2D_plasma_input_3d(i,mp) = N2D(i) ! not used lrm20121115
              endif

              east_plasma_input_3d(i,mp) = east(i)
              south_plasma_input_3d(i,mp) = south(i)
              upward_plasma_input_3d(i,mp) = upward(i)
!             te_plasma_input_3d(i,mp) = te_dum(i)

          enddo
      !g
      !g  end of the big flux tubes loop....
      !g
      enddo
  enddo

  sw_1st_call_int_fixed_ht = .FALSE.

  RETURN



end SUBROUTINE INTERFACE__FIXED_GEO_to_IONOSPHERE












END MODULE moduleInterfaceThermo2Iono
