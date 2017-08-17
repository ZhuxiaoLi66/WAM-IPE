MODULE moduleTHERMOSPHERE

  USE moduleAmplitude
  USE moduleTidalPhase
  ! sizes of thermosphere pressure gridb
  USE modSizeThermo, ONLY : GT_ht_dim, GT_lat_dim, GT_lon_dim  ! = 15, 91, 20

  IMPLICIT NONE

  PUBLIC :: GT_thermosphere_INIT
  PUBLIC :: GT_thermosphere
  PUBLIC :: Foster
  PUBLIC :: calculate_magnetic_parameters_using_apex
  PUBLIC :: low_lat_efield
  PUBLIC :: setThermoLatLons
  PUBLIC :: therm_model_geo_lat_deg, therm_model_geo_long_deg

  PRIVATE

  ! Parameters
  integer, parameter :: n_levels = GT_ht_dim
  integer, parameter :: n_lats = GT_lat_dim
  integer, parameter :: n_lons = GT_lon_dim 

  real(kind=8), parameter :: secondsInDay = 86400.


  REAL(kind=8) :: wind_southwards_ms1(n_levels, n_lats, n_lons)
  REAL(kind=8) :: wind_eastwards_ms1(n_levels, n_lats, n_lons)
  REAL(kind=8) :: wvz(n_levels, n_lats, n_lons)
  REAL(kind=8) :: eps(n_levels, n_lats, n_lons)
  REAL(kind=8) :: rmt(n_levels, n_lats, n_lons)
  REAL(kind=8) :: Temperature_K(n_levels, n_lats, n_lons)
  REAL(kind=8) :: ht(n_levels, n_lats, n_lons)
  REAL(kind=8) :: psao(n_levels, n_lats, n_lons)
  REAL(kind=8) :: psmo(n_levels, n_lats, n_lons)
  REAL(kind=8) :: psmn(n_levels, n_lats, n_lons)

  REAL(kind=8) :: km(n_levels)
  REAL(kind=8) :: kt(n_levels)
  REAL(kind=8) :: tin(n_levels)
  REAL(kind=8) :: umut(n_levels)


  !REAL(kind=8) :: vy0(91)
  REAL(kind=8) :: vy0(n_lats) ! lrm20130829

  REAL(kind=8) :: temp0(91)
  REAL(kind=8) :: temp0av
  REAL(kind=8) :: dh0(91)
  REAL(kind=8) :: htold(n_lats,n_lons)
  REAL(kind=8) :: pressure(n_levels)

  REAL(kind=8) :: examp
  REAL(kind=8) :: fden

  INTEGER :: m_North_Pole
  INTEGER :: m_equator
  INTEGER :: ixx
  INTEGER :: mjk
  INTEGER :: mj1
  INTEGER :: mj2

  integer :: idv_vsouth2 , idv_veast2 , idv_wvz2
  integer :: idv_rmt2 , idv_tn2, idv_ht2
  integer :: idv_psao2 , idv_psmo2
  integer :: idv_Ne2
  integer :: idv_UT2

  integer :: irec_number  ! for writing to the graphics file
  integer :: iout_high  !  counter for # of times gt_thermosphere has been called

 CHARACTER(len=200) :: graphics_file

 REAL*8, parameter :: PI = 3.14159
 REAL*8, parameter :: DTR = PI/180.0 ! lrm20130822
 REAL*8, parameter :: RTD = 180.0/PI

!-------------------------------
! thermosphere lats & lons
!-------------------------------
real*8 :: therm_model_geo_long_deg(GT_lon_dim)
real*8 :: therm_model_geo_lat_deg(GT_lat_dim)


CONTAINS



SUBROUTINE GT_thermosphere( &
                 GT_output_dataset, &  ! input
                 debugDir, & ! input
                 timestep_in_seconds, &           ! input
                 idump_gt, &
                 solar_declination_angle_radians, &
                 nn, nnloop, & ! nn is minute loop number, nnloop is loop # * minute #
                 Universal_Time_seconds, &
                 nn_smoothing_counter, &
                 i_smoothing_frequency, & ! GT smoothing frequency
                 nn_composition_counter, &
                 i_neutral_composition_calling_frequency, &
                 hough11 , hough22 , hough23 , hough24 , hough25, & ! these are constant
                 amplitude, &  ! these are in the user inputs (constant)
                 tidalPhase, &  ! these are last in the namelist file (also constant)
                 Electron_density_m3, &
                 O_plus_density_m3, &
                 NO_plus_density_m3, &
                 O2_plus_density_m3, &
                 Ti_Oplus_K, &
                 useIPEHeatingRates, &  ! input
                 neutralHeatingRates, & ! input
                 exns, eyns, ezns, &   ! electric field, important for storm, will change dimensions eventually
                 Dip_angle_degrees, & ! doesn't change during the run (angle b/t horizon & magnetic field)
                 B_magnitude_nT, & ! magnitude of magnetic field (2D, no height variation for now)
                 Magnetic_latitude_degrees, Magnetic_longitude_degrees, & ! doesn't change in time
                 newl, & ! important for running the storm (index for the level of auroral precipition)
                 gw, & ! hemispheric power input due to auroral precip. (changed every time step for storm run)
                       ! in gigawatts
                 f107, & ! a daily value
                 emaps, cmaps, profile, & ! read at beginning, they don't change during run
                 wind_southwards_ms1_copy, &  ! outputs from GT
                 wind_eastwards_ms1_copy, & ! outputs from GT
                 wvz_copy, & ! outputs from GT
                 rmt_copy, & ! outputs from GT
                 Temperature_K_copy, & ! outputs from GT
                 ht_copy, & ! outputs from GT
                 O_density_copy, & ! outputs from GT
                 O2_density_copy, & ! outputs from GT
                 N2_density_copy, & ! outputs from GT
                 qion3d) ! outputs from GT

!M: move to module parameters 
!N: not needed
!NE: not essential
!NU: not used
!Y: yes needed


      IMPLICIT NONE

     
      CHARACTER*200, INTENT(IN) :: GT_output_dataset  ! Startup filename for the next run
      character(len=*), INTENT(IN)  :: debugDir
      !INTEGER, INTENT(IN) :: i_no_day                 !M:I: nday (input from unit5)
      INTEGER, INTENT(IN) :: timestep_in_seconds    !M:I:  (input from unit5)
      !INTEGER, INTENT(IN)  :: i_graphics_out_start   !NE:I:  (input from unit5)
      INTEGER, INTENT(IN)  :: idump_gt               !NE:I:flag to call write_gt_nc_history
      !INTEGER, INTENT(INOUT)  :: iout_high              !NE:IO:counter to call write_gt_netcdf_graphics
      !INTEGER, INTENT(IN)  :: ipint_high             !NE:I:output frequency for graphics(input from unit5)
      REAL*8, INTENT(IN) :: solar_declination_angle_radians !Y:I:
      INTEGER, INTENT(IN) :: nn , nnloop            !Y:I:
      REAL(kind=8), INTENT(IN) :: Universal_Time_seconds  !Y:I:
      INTEGER, INTENT(INOUT) :: nn_smoothing_counter   !Y:IO:counter for smoothing
      INTEGER, INTENT(IN) :: i_smoothing_frequency  !M:I:smoothing freq (input from unit5)
      INTEGER, INTENT(INOUT) :: nn_composition_counter !Y:IO:counter to call composition
      INTEGER, INTENT(IN) :: i_neutral_composition_calling_frequency !M:I: freq to call compo(input from unit5)
      REAL*8, INTENT(IN) :: hough11(181), hough22(181) !M:I:read in from static_files/hough(unit36/33/40)
      REAL*8, INTENT(IN) :: hough23(181), hough24(181), hough25(181)

      TYPE(amplitudeType), INTENT(IN) :: amplitude      !M:I:tidal amplitudes (input from namelist)
      TYPE(tidalPhaseType), INTENT(IN) :: tidalPhase  ! M:I:tidal phases (input from namelist.  these are constant)

      REAL*8, INTENT(INOUT) :: Electron_density_m3(15,91,20) ! Y:IO:electron density[m3] in pres coords
      REAL*8, INTENT(IN) :: O_plus_density_m3(15,91,20)   !Y:I:O+ density[m3] in pres coords
      REAL*8, INTENT(IN) :: NO_plus_density_m3(15,91,20)  !Y:I:NO+ density[m3] in pres coords
      REAL*8, INTENT(IN) :: O2_plus_density_m3(15,91,20)  !Y:I:O2+ density[m3] in pres coords
      REAL*8, INTENT(IN) :: Ti_Oplus_K(15,91,20)          !Y:I: O+ ion temperature[Kelvin]

      LOGICAL, INTENT(IN) :: useIPEHeatingRates  ! Use the ipe neutral heating rate values ???
      REAL(kind=8) :: neutralHeatingRates(GT_ht_dim, GT_lat_dim, GT_lon_dim)  ! neutral heating rates from ipe

      !Y:I: electric field X/Y/Z components[V/m](needs improvement)
      REAL*8, INTENT(IN) :: exns(2,45,20), eyns(2,45,20), ezns(2,45,20)  

      REAL*8, INTENT(IN) :: Dip_angle_degrees(91,20)      !M:I:(unit175) read in from gip
      REAL*8, INTENT(IN) :: B_magnitude_nT(91,20)
      REAL*8, INTENT(IN) :: Magnetic_latitude_degrees(91,20)  !M:I
      REAL*8, INTENT(IN) :: Magnetic_longitude_degrees(91,20) !M:I
      INTEGER, INTENT(IN) :: newl !Y:I
      REAL*8, INTENT(IN) :: gw    !Y:I
      REAL*8, INTENT(IN) :: f107  !Y:I
      REAL*8, INTENT(IN) :: emaps(21,20,7) , cmaps(21,20,7) , profile(15,21) !M:I

      ! OUTPUTS ----------------------------------------
      REAL(kind=8), INTENT(OUT) :: wind_southwards_ms1_copy(n_levels, n_lats, n_lons)  ! lrm20130731
      REAL(kind=8), INTENT(OUT) :: wind_eastwards_ms1_copy(n_levels, n_lats, n_lons)

      REAL(kind=8), INTENT(OUT) :: wvz_copy(15,91,20)
      REAL(kind=8), INTENT(OUT) :: rmt_copy(15,91,20)
      REAL(kind=8), INTENT(OUT) :: Temperature_K_copy(15,91,20)
      REAL(kind=8), INTENT(OUT) :: ht_copy(15,91,20)
      REAL(kind=8), INTENT(OUT) :: O_density_copy(15,91,20)
      REAL(kind=8), INTENT(OUT) :: O2_density_copy(15,91,20)
      REAL(kind=8), INTENT(OUT) :: N2_density_copy(15,91,20)

      REAL*8, INTENT(OUT) :: qion3d(n_levels,n_lats,n_lons)

      ! Local Variables -------------------------------------------------------
      !----------------
      ! Parameters
      !----------------
      !REAL*8, parameter :: PI = 3.14159
      REAL*8, parameter :: R0 = 6.370E06
      REAL*8, parameter :: R0SQ = R0*R0

      !REAL*8, parameter :: DTR = PI/180.0  moved to module head lrm20130822
      !REAL*8, parameter :: RTD = 180.0/PI  moved to module head lrm20130822
      REAL*8, parameter :: Electron_charge_Coulombs = 1.602E-19

      REAL*8, parameter :: GRAV = 9.5
      REAL*8, parameter :: GSCON = 8.3141E+03
      REAL*8, parameter :: PI2 = PI/2.0


      LOGICAL :: sw_smoothing_this_time_step ! initialized to false
      INTEGER :: ndd , j , mmin , mmax           
      REAL*8 :: s


      REAL(kind=8) :: Universal_Time_hours
      REAL(kind=8) :: realTimeStep
      REAL(kind=8) :: UT_minus_12UT_in_seconds

      REAL*8 :: ped(15,91,20), hall(15,91,20)
      REAL*8 :: dmsp(21,20,5,5) 

      REAL*8 :: htt
      REAL*8 :: a1 , a11 , a12 , a2 , a3 , a4 , a5 , a6 , aa 

      REAL*8 :: angdif 
      REAL*8 :: b1 , b11 , b12 , b2 , b3 , b4 , b5 , b6

      REAL*8, parameter :: ANROT = 7.29E-05
      REAL*8, parameter :: BZ = 1.38E-23

      REAL*8 :: bb , brad , &
                c1 , c2 , c3 , c4 , c5 , &
                c51 , c52 , btheta , bphi

      REAL*8 :: Magnetic_latitude_radians

      REAL*8 :: B_magnitude_Tesla

      REAL*8 :: Dip_angle_radians
      REAL*8 :: Sine_Dip_angle
      REAL*8 :: Cosine_Dip_angle
      REAL*8 :: c6 , c8 , magnetic_colatitude_degrees


      REAL*8 :: cosdif , csza , cth , &
                ddddom , dddom , ddom
      REAL*8 :: ddt , delphi , deltha , &
                dkmkt , &
                dom , &
                source1 , source2 , k3 , k8 , tr

      REAL*8 :: dt, dvydy, &
                elecx , elecz

      REAL*8 :: edep(n_levels)

      REAL*8 :: elecy , &
                magnetic_local_time_in_degrees 
      REAL*8 :: ex2d , exd , &
                ey2d , eyd , &
                fak , fakj , ezd , ez2d
      REAL*8 :: vionx(15,91,20), viony(15,91,20)

      REAL*8 :: fp2 

      REAL*8 :: pdn , pdnd , pdndd , pdnu
      REAL*8 :: pdnuu , longitude_degrees , &
                magnetic_longitude_radians , longitude_radians
                

      REAL*8 :: r_variable
      REAL*8 :: rj1 , latitude_radians , rlt , rm2 , &
                rnumden

      REAL*8 :: sighal , sigped , sindif , &
                smin , splus , ssa , sigpar
      REAL*8 :: sth , sth2 , &
                sx1 , sx3 , sx2 , sx4 , sx5 , sx6 , sx7 , sx8 , sza

      REAL*8 :: longitude_subsolar_point_degrees
      REAL*8 :: latitude_subsolar_point_degrees
      REAL*8 :: magnetic_longitude_subsolar_point_degrees
      REAL*8 :: magnetic_latitude_subsolar_point_degrees
      REAL*8 :: t14 , & ! t15 , &
                tepse , tepsn , tepso , tepss , tepsw , &
                colatitude_radians


      REAL*8 :: colatitude_degrees


      REAL*8 :: volume
      REAL*8 :: vxmd , vxmu
      REAL*8 :: x10 , x11 , x12 , &
                x13 , x8 , x9 , y10 , y11 , y12 , y13
      REAL*8 :: y8 , y9 , z10 , z8 , z9 , bq(15)


      INTEGER :: i1 , i2 , ib , &
                 ihemi , &
                 mgtype, &
                 j1 , j2 , jb , kb

      INTEGER :: l , lb , le , level , lw , m


      INTEGER :: mn , ms , n , &
                 nb , nd ,  nmax , &
                 nmin


      INTEGER :: nu , nxx , nyy

      INTEGER :: yy, ggg

      REAL*8 :: om1(15,91,20) , &
                aeuv(15,91,20) , &
                elx(91,20), ely(91,20), &
                sped(15), teuv(15), &
                elz(91,20)

 !-----------------------------------------------------------
 ! Arrays for checking values of source1, source2
 ! 20130610lrm
 !-----------------------------------------------------------
 REAL*8 :: source1Array(n_levels,n_lats,n_lons), source2Array(n_levels,n_lats,n_lons)
! REAL*8 :: c7Array(n_levels,n_lats,n_lons), k3Array(n_levels,n_lats,n_lons)
 REAL*8 :: a5Array(n_levels,n_lats,n_lons), b5Array(n_levels,n_lats,n_lons)

 REAL*8 :: cp(15) , &
           temp(15), &
           ! neutral_density_3d(15,91,20), & lrm20130731
           scht(15) , om(15) , &
           c7(15) , stfac(15) , fo2(15) , qiont(15) , &
           div(15) , c7_3d(15,91,20), &
           c77_3d(15,91,20)

 REAL*8 :: effqia(n_levels), neutral_density_1d(n_levels)
 REAL*8 :: neutral_density_3d(n_levels, n_lats, n_lons)


 REAL*8 :: factor !  for multipling source1, source2

 REAL*8 :: vxe(15) , vye(15) , epse(15) , hte(15) , vxw(15) , &
                vyw(15) , epsw(15) , htw(15) , vxs(15) , vys(15) , &
                epss(15) , hts(15) , a8(15) , b8(15) , c10(15) , &
                qeuv(15) , qir(15) 

 REAL*8 :: vx2(n_levels,n_lats) , vy2(n_levels,n_lats) , eps2(n_levels,n_lats) 
 REAL*8 :: vxl(n_levels,n_lats) , vyl(n_levels,n_lats) , epsl(n_levels,n_lats)

      REAL*8 :: Vx_1d_copy(15)
      REAL*8 :: Vy_1d_copy(15)
      REAL*8 :: Eps_1d_copy(15)

      REAL*8 :: vpx(15,2) , vpy(15,2) , &
                epsp(15,2), fp1 , fp, &
                phir , sphi , cphi
      INTEGER :: mm , mp , zz


      REAL*8 :: jth(15) , jphi(15) , dvx(15) , dvy(15) , ddvx(15) , &
                ddvy(15) , dumut(15) , &
                c77(15) , &
                jrad(15)

 !REAL*8 :: sum1(15) , sum2(15) , sum3(15)
 REAL*8 :: sum1(n_levels) , sum2(n_levels) , sum3(n_levels)


      REAL*8 :: p1(15) , p2(15) , p33(15) , &
                O_plus_1d(15), &
                NO_plus_1d(15), &
                O2_plus_1d(15), &
                teff(15) , rvin(15) , ramin(15)


      REAL*8 :: ti1(15,91,20)


      DIMENSION ex2d(91,20) , &
                ey2d(91,20) , &
                ez2d(91,20)


       REAL(kind=8) :: prsall , rnt

       REAL*8 :: jrad_3d(15,91,20), jth_3d(15,91,20), jphi_3d(15,91,20)
!
! more blocks for combined models
!
      REAL*8 :: elxr, elyr
      real*8 :: rjth(15), rjphi(15)
      INTEGER :: efield_var, efieldfreq

!-------------------------------------------------
! file units for checking aeuv, source1, source2 lrm20130610
!-------------------------------------------------
INTEGER, parameter :: aeuvUnit = 5000, source1Unit = 5010, source2Unit = 5020
INTEGER, parameter :: psmnUnit = 5030, rmtUnit = 5040, neutral_density_3dUnit = 5050
INTEGER, parameter :: OplusUnit = 5060
INTEGER, parameter :: a5Unit = 5070, b5Unit = 5080, c7Unit = 5090, k3Unit = 6000

LOGICAL, parameter :: checkHeatingSources = .FALSE.



! Small scale variation in E-Field
      DATA efield_var/0/  ! 0 old, 1 just X2 Joule heating, 2 random
                          ! number variation on E-field, 3 rms added
                          ! to efield (not right)
      DATA efieldfreq/1/  ! frequency randomness called (has to be 1)

!c  **
!    polar cap limits in magnetic latitude
      DATA effqia/8*.55 , .54 , .42 , .28 , .15 , .07 , .05 , .03/

      DATA fo2/0.2096 , 0.2095 , 0.2090 , 0.1997 , 0.1758 , 0.1567 , &
           0.1326 , 0.1085 , 0.08394 , 0.05991 , 0.03922 , 0.02250 , &
           0.01100 , 0.004809 , 0.001963/

       data bq/0., 0., 0., 0., 0., 0., 0., 0., 2.e+1, 2.3e+4, 8.5e5, &
           14.8e6, 8.7e6, 8.3e6, 3.9e6/

!=======================================================================
! BEGIN CODE
!=======================================================================  
 !------------------------------------------------
 ! Initialize source1Array, source2Array to 0 lrm20130610
 ! open output files
 !------------------------------------------------
 if (checkHeatingSources) then 
     source1Array = 0.
     source2Array = 0.
     !c7Array = 0.
     !k3Array = 0.
     !a5Array = 0.
     !b5Array = 0.
     !OPEN (aeuvUnit, FILE=TRIM(debugDir)//'aeuv.txt')
     !print *,'debugdir = ',debugdir
     OPEN (source1Unit, FILE=TRIM(debugDir)//'source1.txt')
     OPEN (source2Unit, FILE=TRIM(debugDir)//'source2.txt')
     !OPEN (psmnUnit, FILE=TRIM(debugDir)//'psmn.txt')
     !OPEN (rmtUnit, FILE=TRIM(debugDir)//'rmt.txt')
     !OPEN (neutral_density_3dUnit, FILE=TRIM(debugDir)//'neutral_density_3d.txt')
     OPEN (OplusUnit, FILE=TRIM(debugDir)//'OplusPressure.txt')
     !OPEN (a5Unit, FILE=TRIM(debugDir)//'a5.txt')
     !OPEN (b5Unit, FILE=TRIM(debugDir)//'b5.txt')
     !OPEN (c7Unit, FILE=TRIM(debugDir)//'c7.txt')
     !OPEN (k3Unit, FILE=TRIM(debugDir)//'k3.txt')


 endif


ndd = 0
j = 0
mmin = 0
mmax = 0         
s = 0
Universal_Time_hours = 0
realTimeStep = 0
UT_minus_12UT_in_seconds = 0
ped  = 0
hall = 0
dmsp = 0
htt = 0
a1 = 0  
a11 = 0  
a12 = 0  
a2 = 0  
a3 = 0  
a4 = 0  
a5 = 0  
a6  = 0 
aa  = 0
angdif  = 0
b1 = 0 
b11 = 0 
b12 = 0 
b2 = 0 
b3 = 0 
b4 = 0 
b5 = 0 
b6 = 0
bb = 0 
brad = 0                
c1 = 0 
c2 = 0 
c3 = 0 
c4 = 0 
c5 = 0                 
c51 = 0 
c52 = 0 
btheta = 0 
bphi = 0
Magnetic_latitude_radians = 0
B_magnitude_Tesla = 0
Dip_angle_radians = 0
Sine_Dip_angle = 0
Cosine_Dip_angle = 0
c6 = 0 
c8 = 0 
magnetic_colatitude_degrees = 0
cosdif = 0 
csza = 0 
cth = 0                 
ddddom = 0 
dddom = 0 
ddom = 0      
ddt = 0 
delphi = 0 
deltha = 0            
dkmkt  = 0                  
dom   = 0              
source1  = 0   
source2  = 0   
k3   = 0  
k8   = 0  
tr = 0   
dt  = 0  
dvydy   = 0              
edep   = 0   
elecx   = 0  
elecz = 0 
elecy      = 0              
magnetic_local_time_in_degrees  = 0     
ex2d  = 0   
exd   = 0                 
ey2d  = 0   
eyd  = 0                   
fak  = 0   
fakj  = 0   
ezd   = 0  
ez2d = 0      
vionx  = 0  
viony = 0      
fp2  = 0      
pdn   = 0  
pdnd  = 0   
pdndd  = 0   
pdnu = 0     
pdnuu   = 0  
longitude_degrees  = 0                   
magnetic_longitude_radians  = 0   
longitude_radians = 0                      
r_variable = 0      
rj1   = 0  
latitude_radians  = 0   
rlt   = 0  
rm2   = 0                  
rnumden = 0     
sighal   = 0  
sigped  = 0   
sindif  = 0                   
smin  = 0   
splus = 0    
ssa   = 0  
sigpar = 0      
sth  = 0   
sth2  = 0                   
sx1  = 0   
sx3  = 0   
sx2  = 0   
sx4 = 0    
sx5 = 0    
sx6  = 0   
sx7  = 0   
sx8  = 0   
sza = 0
longitude_subsolar_point_degrees = 0
latitude_subsolar_point_degrees = 0
magnetic_longitude_subsolar_point_degrees = 0
magnetic_latitude_subsolar_point_degrees = 0      
t14  = 0   
!t15  = 0               
tepse  = 0   
tepsn  = 0   
tepso  = 0   
tepss  = 0   
tepsw  = 0                  
colatitude_radians = 0
colatitude_degrees = 0    
volume = 0     
vxmd   = 0  
vxmu = 0    
x10 = 0    
x11 = 0    
x12 = 0                    
x13  = 0   
x8  = 0   
x9  = 0   
y10  = 0   
y11  = 0   
y12  = 0   
y13 = 0      
y8 = 0 
y9 = 0
z10 = 0 
z8 = 0
z9 = 0      
i1  = 0   
i2  = 0   
ib  = 0                    
ihemi  = 0                    
mgtype  = 0                   
j1  = 0   
j2  = 0   
jb  = 0   
kb = 0     
l  = 0   
lb  = 0   
le   = 0  
level  = 0   
lw   = 0  
m = 0   
mn   = 0  
ms  = 0   
n  = 0                  
nb  = 0   
nd  = 0    
nmax  = 0               
nmin = 0
nu  = 0 
nxx  = 0 
nyy = 0 
yy = 0  
om1   = 0              
aeuv = 0         
elx  = 0  
ely  = 0  
sped   = 0  
teuv      = 0           
elz = 0   
cp   = 0               
temp  = 0    
neutral_density_1d  = 0             
neutral_density_3d  = 0              
scht    = 0  
om  = 0                
c7   = 0   
stfac   = 0      
qiont   = 0                 
div    = 0  
c7_3d   = 0            
c77_3d = 0      
vxe  = 0    
vye  = 0    
epse  = 0    
hte  = 0    
vxw  = 0                    
vyw  = 0    
epsw = 0     
htw  = 0    
vxs  = 0    
vys  = 0                    
epss  = 0    
hts  = 0    
a8  = 0    
b8  = 0    
c10  = 0                    
qeuv  = 0    
qir   = 0                   
vx2  = 0   
vy2  = 0   
eps2 = 0                 
vxl  = 0   
vyl   = 0  
epsl = 0     
Vx_1d_copy  = 0
Vy_1d_copy  = 0
Eps_1d_copy  = 0      
vpx   = 0  
vpy   = 0                  
epsp  = 0  
fp1  = 0   
fp   = 0                 
phir  = 0   
sphi   = 0  
cphi = 0      
mm  = 0   
mp  = 0   
zz = 0   
jth  = 0    
jphi  = 0    
dvx  = 0    
dvy  = 0    
ddvx  = 0                    
ddvy   = 0   
dumut  = 0    
sum1  = 0    
sum2  = 0    
sum3  = 0                 
c77   = 0               
jrad  = 0      
p1  = 0    
p2 = 0     
p33  = 0                   
O_plus_1d  = 0                   
NO_plus_1d = 0                    
O2_plus_1d  = 0                   
teff  = 0    
rvin  = 0    
ramin  = 0      
ti1 = 0     
prsall  = 0  
rnt = 0     
jrad_3d = 0  
jth_3d  = 0  
jphi_3d = 0      
elxr  = 0  
elyr = 0      
rjth   = 0  
rjphi  = 0
      
!----------------------------------
! Initialize output variables to 0
!----------------------------------
wind_southwards_ms1_copy = 0
wind_eastwards_ms1_copy = 0
wvz_copy = 0
rmt_copy = 0
Temperature_K_copy = 0
ht_copy = 0
O_density_copy = 0
O2_density_copy = 0
N2_density_copy = 0
qion3d = 0

! -------------------------------------------

 !print *,'GT_thermoshpere : timestep_in_seconds = ',  timestep_in_seconds  
 realTimeStep = REAL( timestep_in_seconds )

 ! 20110927 lrm - changed 86400. to secondsInDay

 UT_minus_12UT_in_seconds = universal_time_seconds - 43200.

 if ( UT_minus_12UT_in_seconds < 0.0 ) then
      UT_minus_12UT_in_seconds = UT_minus_12UT_in_seconds + secondsInDay
 endif

 UT_minus_12UT_in_seconds = DMOD(UT_minus_12UT_in_seconds, secondsInDay)

 !print *,'moduleThermosphere : UT_minus_12UT_in_seconds = ',UT_minus_12UT_in_seconds

 mgtype = 2
 aa = 180.0
 yy = 2
 bb = 4.0
 deltha = bb*DTR
 delphi = 36.0*DTR
 sx1 = deltha*R0
 sx3 = (sx1*sx1)/4.0
 sx5 = sx3/R0SQ
 nyy = yy
 nxx = ixx

 IF ( bb == 4. ) nyy = yy + 1
 IF ( bb == 4. ) nxx = ixx - 1

 !************************
 !    longitude loop l   *
 !************************
 lonloop1740 : DO l = 1 , n_lons

    le = l + 1
    lw = l - 1

    !g
    !g Ingo's version of the code had these next two lines in it ??????
    !g
    IF (lw  ==  0)  lw = 20
    IF (le  ==  21) le = 1

    !longitude_degrees = (FLOAT(l)-1.0)*18.0
    longitude_degrees = therm_model_geo_long_deg(l)  ! lrm20130731

    longitude_radians = longitude_degrees*DTR
    ssa = longitude_degrees + UT_minus_12UT_in_seconds/240.0

    IF ( ssa >= 360.0 ) ssa = ssa - 360.0
 
    rlt = 180.0 + ssa
    IF ( rlt > 360.0 ) rlt = rlt - 360.0
    rlt = rlt*DTR

    !DO  n = 2 , 14
    !           vx2(n,nyy-1) = Wind_southwards_ms1(n,nyy-1,l)
    !           vy2(n,nyy-1) = Wind_eastwards_ms1(n,nyy-1,l)
    !           eps2(n,nyy-1) = eps(n,nyy-1,l)
    !END DO
    vx2(2:14,nyy-1) = Wind_southwards_ms1(2:14,nyy-1,l)  ! array assignment lrm20130731
    vy2(2:14,nyy-1) = Wind_eastwards_ms1(2:14,nyy-1,l)
    eps2(2:14,nyy-1) = eps(2:14,nyy-1,l)


    !***********************
    !    latitude loop m   *
    !***********************
    !DO 1720 m = nyy , nxx
    latloop1720 : DO m = nyy , nxx
       !print *,'line 424 : GT_thermosphere : m = ',m

       !- 
       ! set height independent parameters
       !-  
       mn = m + 1
       ms = m - 1
       colatitude_degrees = aa - (m-1.0)*bb/2.0
       colatitude_radians = colatitude_degrees*DTR

       !- check what hemisphere
       IF ( m <= m_equator - 1) THEN
            ihemi = 2
       ELSE
            ihemi = 1
       ENDIF

       splus = SIN(colatitude_radians+deltha/4.0)
       smin = SIN(colatitude_radians-deltha/4.0)
       sth = SIN(colatitude_radians)
       sth2 = sth*sth
       cth = COS(colatitude_radians)
       sx2 = delphi*R0*sth
       sx4 = (sx2*sx2)/4.0
       sx6 = sx4/R0SQ
       sx7 = deltha*sth
       sx8 = R0*sth
  
       latitude_radians = PI2 - colatitude_radians
 
       magnetic_latitude_radians = Magnetic_latitude_degrees(m,l)*dtr
       magnetic_longitude_radians = Magnetic_longitude_degrees(m,l)*dtr
       magnetic_colatitude_degrees = 90. - Magnetic_latitude_degrees(m,l)

       !
       !   calculation of the MLT requires the magnetic coords of the
       !   sub solar point.....
       !
       latitude_subsolar_point_degrees = solar_declination_angle_radians/dtr
       longitude_subsolar_point_degrees = 360. - UT_minus_12UT_in_seconds/240.

       call GEO2MG_apex(Magnetic_longitude_degrees, &
                               Magnetic_latitude_degrees, &
                               longitude_subsolar_point_degrees, &
                               latitude_subsolar_point_degrees, &
                               magnetic_longitude_subsolar_point_degrees, &
                               magnetic_latitude_subsolar_point_degrees)
       !g
       !g  essa is then the MLT in degrees....
       !g
       magnetic_local_time_in_degrees = magnetic_longitude_degrees(m,l) - &
                                        magnetic_longitude_subsolar_point_degrees

       IF ( magnetic_local_time_in_degrees >= 360.0 ) &
                   magnetic_local_time_in_degrees = magnetic_local_time_in_degrees - 360.0
       IF ( magnetic_local_time_in_degrees < 0.0 ) &
                   magnetic_local_time_in_degrees = magnetic_local_time_in_degrees + 360.
 
       rm2 = ixx - magnetic_colatitude_degrees/bb*2.

       !IF ( rm2 <= mj1 ) fp2 = -1.
       !IF ( rm2 <= mj1 ) rj1 = rm2
       IF ( rm2 <= mj1 ) then  ! lrm20130731
            fp2 = -1.
            rj1 = rm2
       END IF

       !IF ( rm2 >= mj2 ) fp2 = 1.
       !IF ( rm2 >= mj2 ) rj1 = ixx - rm2
       IF ( rm2 >= mj2 ) then  ! lrm20130731
            fp2 = 1.
            rj1 = ixx - rm2
       END IF

       !--------------------------------------------------------
       !-- set rj1 to point back into the electric field table
       !-- ( for expanded elec. field runs)
       !--------------------------------------------------------
       rj1 = rj1/fden

       B_magnitude_Tesla = B_magnitude_nT(m,l)*1.E-09
       Dip_angle_radians = Dip_angle_degrees(m,l)*DTR
       Sine_Dip_angle = sin(Dip_angle_radians)
       Cosine_Dip_angle = cos(Dip_angle_radians)

       ! angdif is the angle between mag and geo. frames
 
       angdif = (magnetic_local_time_in_degrees-ssa)*DTR
       cosdif = COS(angdif)
       sindif = SIN(angdif)

        brad   = -B_magnitude_Tesla*Sine_Dip_angle
        btheta = -B_magnitude_Tesla*Cosine_Dip_angle*cosdif
        bphi   = -B_magnitude_Tesla*Cosine_Dip_angle*sindif

        sza = ACOS(-COS(latitude_radians)*COS(solar_declination_angle_radians)*COS(rlt) &
                        +SIN(latitude_radians) &
                        *SIN(solar_declination_angle_radians))
        csza = COS(sza)

        elecx = 0.0
        elecy = 0.0

        ! new lower boundary insert:
        ! Tides at level 1. Lower boundary winds and heights are
        ! calculated here (once for each time step). 
        ! IMW August 1995

        htold(m,l) = ht(1,m,l)

        CALL TIDES(nn, m, l, hough11, hough22, hough23, &
                        hough24, hough25, &
                        amplitude, &
                        ht, &
                        Wind_southwards_ms1, Wind_eastwards_ms1, &
                        realTimeStep, rmt, temp, &
                        tidalPhase, &
                        temp0, vy0, temp0av, dh0)

        CALL SPECIFIC_HEAT(psao(1,m,l), psmo(1,m,l), psmn(1,m,l), cp)


        levelsloop1625 : DO n = 1 , n_levels

           nu = n + 1
           nd = n - 1

           !IF ( n <= 14 .AND. n > 1 ) THEN
           IF ( n <= (n_levels-1) .AND. n > 1 ) THEN
                temp(n) = (-(Wind_southwards_ms1(n,m,l)**2 + &
                             Wind_eastwards_ms1(n,m,l)**2)/2.0 + eps( n,m,l))/cp(n)
           ELSE
                temp(15) = temp(14)
                eps(15,m,l) = cp(15)*temp(15) &
                              + (Wind_southwards_ms1(15,m,l)**2+ &
                              Wind_eastwards_ms1(15,m,l)**2)/2.0
           ENDIF

           !print *,'GT_thermosphere : Wind_southwards_ms1(n,m,l) = ',Wind_southwards_ms1(n,m,l)
           !print *,'GT_thermosphere : Wind_eastwards_ms1(n,m,l) = ', Wind_eastwards_ms1(n,m,l)
           !print *,'GT_thermosphere : eps(n,m,l), cp(n) = ', eps(n,m,l), cp(n)


           ! ****************************************************************
           ! ** modification of lower boundary temp inserted on 9 Feb 1995 **
           ! ** this is the place in the code where lower boundary temps   **
           ! ** need to be changed - if it is done during setup only the   **
           ! ** information gets lost when reading in eps from previous    **
           ! ** run. Since lat loop doesn't automatically go from 1 to 91  **
           ! ** need to set the pole values extra, otherwise it uses those **
           ! ** commented all level 1 temperature settings out since they  **
           ! ** are now oscillated in the TIDES routine                  **
           ! ** IMW, OCTOBER 1995                                    **
           ! ****************************************************************
           IF (n == 1) THEN
               eps(1,m,l) = temp(1) * cp(1) + &
                            (Wind_southwards_ms1(1,m,l)**2 + &
                             Wind_eastwards_ms1(1,m,l)**2)/2.0

               ! the following pole settings are necessary
               IF (m == 3) THEN
                   eps(1,1,l) = eps(1,m,l)
                   eps(1,2,l) = eps(1,m,l)
               ENDIF

               eps(1,90,l) = eps(1,89,l)
               eps(1,91,l) = eps(1,89,l)
           ENDIF

           ! ** end of modification ***

           !print *,'GT_thermosphere : n, pressure(n) = ',n, pressure(n)
           !print *,'GT_thermosphere : rmt(n,m,l) = ', rmt(n,m,l)
           !print *,'GT_thermosphere : GSCON, temp(n) = ',GSCON, temp(n)

           neutral_density_1d(n) = (pressure(n)*rmt(n,m,l))/(GSCON*temp(n))


           neutral_density_3d(n,m,l) = neutral_density_1d(n)
           scht(n) = (GSCON*temp(n))/(rmt(n,m,l)*GRAV)
           IF ( m /= yy ) vxmu = Wind_southwards_ms1(n,m,l) + vx2(n,ms)

           vxmd = Wind_southwards_ms1(n,m,l) &
                  + Wind_southwards_ms1(n,mn,l)

           IF ( m  ==  2 ) vxmu = Wind_southwards_ms1(n,m,l) &
                                  + Wind_southwards_ms1(n,ms,l)

           IF ( l  ==  1 ) THEN
                dvydy = (Wind_eastwards_ms1(n,m,2) &
                         -Wind_eastwards_ms1(n,m,20))/delphi
           ELSEIF ( l  ==  20 ) THEN
                        dvydy = (vyl(n,m)-vy2(n,m))/delphi
           ELSE
                        dvydy = (Wind_eastwards_ms1(n,m,le) &
                                -vy2(n,m))/delphi
           ENDIF


           ! div is horizontal divergence of velocity

           IF ( m  ==  1 ) THEN
                div(n) = ((-2.0*vxmd*smin)/deltha+dvydy)/R0
           ELSE
                div(n) = ((vxmu*splus-vxmd*smin)/deltha+dvydy) &
                                 /(R0*sth)
           ENDIF

     END DO levelsloop1625

     !-------------------------------------------------------------------------
     !  calculation of w(=dp/dt) (called om here) by expansion using Taylor
     !  series, see Tim's thesis, p.108
     !  w is the vertical velocity relative to the pressure level. Later,
     !  the total vertical velocity wvz is calculated using w (ie.om)
     !  and the velocity of the pressure level itself.
     !-------------------------------------------------------------------------
      !DO 1630 n = 1 , 13
      !levelsloop1630 : DO n = 1, 13
      levelsloop1630 : DO nb = 14, 2, -1  ! lrm20130731

         !nb = 15 - n  ! lrm20130731

         ib = nb + 1
         jb = nb - 1
         om(15) = 0.0
         pdn = pressure(nb)*div(nb)
         pdnu = pressure(ib)*div(ib)
         pdnd = pressure(jb)*div(jb)
         ddddom = 0.0
         dddom = 0.0

         IF ( nb /= 2 .AND. nb /= 14 ) THEN
              kb = nb + 2
              lb = nb - 2
              pdnuu = pressure(kb)*div(kb)
              pdndd = pressure(lb)*div(lb)
              ddddom = (pdnuu-2.0*pdnu+2.0*pdnd-pdndd)/2.0
              dddom = (pdnu+pdnd-2.0*pdn)
         ENDIF

         dom = pdn
         ddom = (pdnu-pdnd)/2.0
         om(nb) = om(ib) - dom - ddom/2.0 - dddom/6.0 - &
                  ddddom/24.0
         om1(nb,m,l) = om(nb)

         IF ( nb == 2 ) THEN
              om(1) = om(2) - pressure(1)*div(1)
              om1(1,m,l) = om(1)
         ENDIF

    END DO levelsloop1630

    c8 = 0.0

    !DO n = 1 , 15
    !   edep(n) = 0.0
    !END DO
    edep = 0.0  ! do whole array assignment to 0 lrm20130731
            

    !
    !c  use shaun switch for foster electric field in gtms simulations
    !c
    !c           IF ( rm2 <= mj1 .OR. rm2 >= mj2 ) THEN

    !   expanded electric field

    j1 = INT(rj1)

    IF ( j1 < 1 ) j1 = 1
    IF ( j1 > mjk ) j1 = mjk
    j2 = j1 + 1
    IF ( j1  >=  mjk ) j2 = mjk
    fak = magnetic_local_time_in_degrees/18.0
    i1 = INT(fak) + 1
    IF ( i1 > 20 ) i1 = 20
    i2 = i1 + 1
    IF ( i2 == 21 ) i2 = i2 - 20
    fak = fak - INT(fak)
    fakj = rj1 - j1


    !-----------------------------------------
    !  use a2 ; b2 electric field models
    !-----------------------------------------
    elecx = exns(ihemi,j1,i1)*(1.-fak)*(1.-fakj) &
            + exns(ihemi,j1,i2)*(1.-fakj) &
            *fak + exns(ihemi,j2,i2)*fak*fakj + &
            exns(ihemi,j2,i1)*fakj*(1.-fak)

    elecy = eyns(ihemi,j1,i1)*(1.-fak)*(1.-fakj) &
            + eyns(ihemi,j1,i2)*(1.-fakj) &
            *fak + eyns(ihemi,j2,i2)*fak*fakj + &
            eyns(ihemi,j2,i1)*fakj*(1.-fak)

    elecz = ezns(ihemi,j1,i1)*(1.-fak)*(1.-fakj) &
            + ezns(ihemi,j1,i2)*(1.-fakj) &
            *fak + ezns(ihemi,j2,i2)*fak*fakj + &
            ezns(ihemi,j2,i1)*fakj*(1.-fak)

    !c  ** ------------------------------------------------------
    !c  the potential patterns are perp. to magnetic field.
    !c  the horizontal component elecx should *sin(dip)
    !c  original separations will also have factor sin(dip)
    !c  therefore they cancel.
    !g ----------------------------------------------------------
    exd = (elecx*cosdif-elecy*sindif)*fp2
    eyd = (elecx*sindif+elecy*cosdif)
    ezd = elecz

    !c  **
    !c  scale up foster/holt electric fields
    !c  this should be done in foster for potl also.
    !c               elecx=exd*1.3
    !c               elecy=eyd*1.3
    !
    !  expanded electric field
    !
    elecx = exd*examp
    elecy = eyd*examp
    elecz = ezd*examp

    !            ENDIF 
    !
    !  end of expanded electric field inserts
    !
    elx(m,l) = elecx
    ely(m,l) = elecy
    elz(m,l) = elecz

    elecx = elx(m,l)
    elecy = ely(m,l)

    !-----------------------------
    ! call tiros if required
    !-----------------------------
     qiont = 0.0 ! whole array assignment to 0 lrm20130731

     level = newl

     IF ( ABS(magnetic_latitude_degrees(m,l)) > 50. ) THEN

          CALL TIROS(magnetic_latitude_degrees(m,l), magnetic_local_time_in_degrees, &
                           level, qiont, neutral_density_1d, gw, &
                           emaps, cmaps, profile, dmsp)
     ENDIF

     !DO n = 1 , 15
     !   qiont(n) = qiont(n) + bq(n)
     !   qion3d(n,m,l) = qiont(n)
     !END DO
     qiont = qiont + bq  ! whole array assignment lrm20130731
     qion3d(:,m,l) = qiont(:)

     ex2d(m,l) = elecx
     ey2d(m,l) = elecy
     ez2d(m,l) = elecz
!g

     !do n = 1 , 15
     !   ti1(n,m,l) = Ti_Oplus_K(n,m,l)
     !   !telec(n,m,l) = Te_K(n,m,l)  NOT USED FOR ANYTHING
     !enddo 
     ti1(:,m,l) = Ti_Oplus_K(:,m,l)  ! array assignment lrm20130731

     !g
     !g  New insert to allow for electric field variation - taken from CMAT model
     !g  originally from Mihail Codrescu....
     !g
     ! calculate small scale electric field variation, add on for Joule heating
     ! later, but not on E-field as mean acceleration shouldn't be affected
     elxr = 0.0
     elyr = 0.0

     !------------------------------------
     ! if using random number generator
     !------------------------------------
     if ((efield_var == 2) .and. (nnloop < 5 .or. mod(nnloop, efieldfreq) == 0)) then

 4416   continue

        !-------------------------------------------------------------------------
        ! in future the standard deviation (here 1.e-2)  will be a function of
        ! location, TIROS level, and possibly even F10.7.
        !                   if(abs(elecx) > 0.001) elxr=g05ddf(0., 1.e-2)
        !                   if(abs(elecy) > 0.001) elyr=g05ddf(0., 1.e-2)
        !-------------------------------------------------------------------------
        write(6,2177) nnloop, m, l, elecx, elxr, elecy, elyr
 2177   format('random field ', 3i4, 4e12.4)

        ! restrict to 2 sigma
        if (abs(elxr) > 2.0e-2) then
            write(6,*) 'elxr too large', elxr
            goto 4416
        endif

        if (abs(elyr) > 2.0e-2) then
            write(6,*) 'elyr too large', elyr
            goto 4416
        endif

      endif ! ((efield_var == 2) .and. (nnloop < 5 .or. mod(nnloop, efieldfreq) == 0))


      if (efield_var == 3) then
          if(abs(elecx) > 0.001) elxr=1.e-2
          if(abs(elecy) > 0.001) elyr=1.e-2
      endif

      !
      ! set energy deposition
      !
      !DO n = 1 , 15
      !   edep(n) = qion3d(n,m,l)*effqia(n)*5.6E-18/neutral_density_1d(n)  lrm20130731
      !END DO
      edep(:) = qion3d(:,m,l)*effqia(:)*5.6E-18/neutral_density_1d(:)
 

      !print *,'GT_thermosphere : line 806 : b/f SOLAR_EUV, neutral_density_1d = ',neutral_density_1d
      CALL SOLAR_EUV(f107, rmt(1,m,l), scht, neutral_density_1d, sza, teuv)


      CALL INFRARED_COOLING(temp, ht(1,m,l), neutral_density_1d,QIR)


      !DO n = 2 , 14
      !    aeuv(n,m,l) = teuv(n)
      !END DO
      aeuv(2:14,m,l) = teuv(2:14)

      !DO n = 2 , 14
      !   qeuv(n) = aeuv(n,m,l)
      !END DO
      qeuv(2:14) = aeuv(2:14,m,l)


      !DO n = 1 , 15
      !   Vx_1d_copy(n) = Wind_southwards_ms1(n,m,l)
      !   Vy_1d_copy(n) = Wind_eastwards_ms1(n,m,l)
      !   Eps_1d_copy(n) = eps(n,m,l)
      !END DO
      Vx_1d_copy(:) = Wind_southwards_ms1(:,m,l)
      Vy_1d_copy(:) = Wind_eastwards_ms1(:,m,l)
      Eps_1d_copy(:) = eps(:,m,l)


      !DO n = 1 , 14
      DO n = 1 , n_levels-1
         IF ( l == 20 ) THEN
              vxe(n) = vxl(n,m)
              vye(n) = vyl(n,m)
              epse(n) = epsl(n,m)
              hte(n) = ht(n,m,1)
         ELSE
              vxe(n) = Wind_southwards_ms1(n,m,le)
              vye(n) = Wind_eastwards_ms1(n,m,le)
              epse(n) = eps(n,m,le)
              hte(n) = ht(n,m,le)
         ENDIF

         vxs(n) = vx2(n,ms)
         vys(n) = vy2(n,ms)
         epss(n) = eps2(n,ms)
         hts(n) = ht(n,ms,l)

         IF ( l == 1 ) THEN
            vxw(n) = Wind_southwards_ms1(n,m,20)
            vyw(n) = Wind_eastwards_ms1(n,m,20)
            epsw(n) = eps(n,m,20)
            htw(n) = ht(n,m,20)
         ELSE
            vxw(n) = vx2(n,m)
            vyw(n) = vy2(n,m)
            epsw(n) = eps2(n,m)
            htw(n) = ht(n,m,lw)
         ENDIF

      END DO ! n = 1 , n_levels-1

      !---
      !---
      !DO n = 2 , 14
      DO n = 2, n_levels-1
         tepso = eps(n,m,l) - (Wind_southwards_ms1(n,m,l)**2 &
                             + Wind_eastwards_ms1(n,m,l)**2) &
                             /2.0
         tepss = epss(n) - (vxs(n)**2+vys(n)**2)/2.0
         tepsn = eps(n,mn,l) - (Wind_southwards_ms1(n,mn,l) &
                             **2+Wind_eastwards_ms1(n,mn,l)**2) &
                             /2.0

         tepse = epse(n) - (vxe(n)**2+vye(n)**2)/2.0
         tepsw = epsw(n) - (vxw(n)**2+vyw(n)**2)/2.0

         x8 = (vxs(n)-2.0*Wind_southwards_ms1(n,m,l) &
                           +Wind_southwards_ms1(n,mn,l))/sx5
         x9 = cth*(vxs(n)-Wind_southwards_ms1(n,mn,l))/sx7
         x10 = (vxe(n)-2.0*Wind_southwards_ms1(n,m,l) &
                           +vxw(n))/sx6
         x11 = -Wind_southwards_ms1(n,m,l)/sth2
         x12 = -2.*cth*(vye(n)-vyw(n))/(sth2*delphi)
         x13 = 2.*Wind_southwards_ms1(n,m,l)

         a8(n) = umut(n)*scht(n)*(x8+x9+x10+x11+x12+x13) &
                             /(R0SQ*neutral_density_1d(n))
         y8 = (vys(n)-2.0*Wind_eastwards_ms1(n,m,l) &
                          +Wind_eastwards_ms1(n,mn,l))/sx5
         y9 = cth*(vys(n)-Wind_eastwards_ms1(n,mn,l))/sx7
         y10 = (vye(n)-2.0*Wind_eastwards_ms1(n,m,l) &
                           +vyw(n))/sx6
         y11 = -Wind_eastwards_ms1(n,m,l)/sth2
         y12 = 2.*cth*(vxe(n)-vxw(n))/(sth2*delphi)
         y13 = 2.*Wind_eastwards_ms1(n,m,l)
         b8(n) = umut(n)*scht(n)*(y8+y9+y10+y11+y12+y13) &
                             /(R0SQ*neutral_density_1d(n))
         z8 = (tepss-2.0*tepso+tepsn)/sx5
         z9 = cth*(tepss-tepsn)/sx7
         z10 = (tepse-2.0*tepso+tepsw)/sx6
         c10(n) = (km(n)+kt(n))*scht(n)*(z8+z9+z10) &
                              /(R0SQ*cp(n)*neutral_density_1d(n))

        END DO ! n = 2, 14


        !write(6,*) 'GT_thermosphere : Temperature_K = ',Temperature_K  

        ! THIS IS THE PROPER WAY  *********************
        teff(:) = (Temperature_K(:,m,l) + Ti1(:,m,l))/2.0

        ! FOR DEBUGGING ONLY ******* :
        !teff(n) = Temperature_K(n,m,l)

        O_plus_1d(:) = O_plus_density_m3(:,m,l)
        NO_plus_1d(:)= NO_plus_density_m3(:,m,l)
        O2_plus_1d(:)= O2_plus_density_m3(:,m,l)

        !DO n = 1 , 15
        DO n = 1 , n_levels
           !! THIS IS THE PROPER WAY  *********************
           !teff(n) = (Temperature_K(n,m,l) + Ti1(n,m,l))/2.0

           !! FOR DEBUGGING ONLY ******* :
           !!teff(n) = Temperature_K(n,m,l)

           rnumden = pressure(n)/(BZ*Temperature_K(n,m,l))
           p1(n) = rnumden*psao(n,m,l)*rmt(n,m,l)/16.
           p2(n) = rnumden*psmo(n,m,l)*rmt(n,m,l)/32.
           p33(n) = rnumden*psmn(n,m,l)*rmt(n,m,l)/28.
        END DO ! n = 1, 15


        !print *,'GT_thermosphere : p1 = ',p1
        !print *,'GT_thermosphere : p2 = ',p2
        !print *,'GT_thermosphere : teff = ',teff
        !print *,'GT_thermosphere : m, l = ',m, l

        CALL IONNEUT(p1, p2, p33, O_plus_1d, NO_plus_1d, &
                               O2_plus_1d, teff, &
                               rvin, ramin)
        !g
        !DO 1695 n = 2 , 14 ! ------------------------------------------------------------------
        looplevels1695 : DO n = 2 , n_levels-1 ! ------------------------------------------------------------------

           nu = n + 1
           nd = n - 1

           r_variable = (ramin(n)*rvin(n))/(Electron_charge_Coulombs*B_magnitude_Tesla)
           sigped = (Electron_density_m3(n,m,l)*Electron_charge_Coulombs*r_variable)/ &
                           (B_magnitude_Tesla*(1.0+r_variable**2))
           ped(n,m,l) = sigped
           sped(n) = sigped
           sighal = sigped*r_variable
           hall(n,m,l) = sighal
           SIGPAR=(Electron_charge_Coulombs*Electron_density_m3(n,m,l))/(B_magnitude_Tesla*r_variable)

           !c  **
           !c  insert full layer conductivities for low latitude electrodynamics
           !c  **

           jth(n) = sigped*(elecx+Wind_eastwards_ms1(n,m,l)*brad) &
                    + sighal*(elecy-Wind_southwards_ms1(n,m,l)*brad)*Sine_Dip_angle
           jphi(n)= &
                  +sigped*(elecy-brad*Wind_southwards_ms1(n,m,l)) &
                  -sighal*(elecx+brad*Wind_eastwards_ms1(n,m,l))/Sine_Dip_angle
           jrad(n) = sigped*(elecz-Wind_eastwards_ms1(n,m,l)*btheta &
                    +Wind_southwards_ms1(n,m,l)*bphi) &
                    -sighal*(elecy-Wind_southwards_ms1(n,m,l)*brad)*Cosine_Dip_angle

           jth_3d(n,m,l) = jth(n)
           jphi_3d(n,m,l) = jphi(n)
           jrad_3d(n,m,l) = jrad(n)


           !g
           ! electric fields with variability used in Joule heating
           if (efield_var == 2 .or. efield_var == 3) then

               if (efield_var == 3) then  ! if not using random (not
                       if(elecy*elyr == -1) elyr=-1.*elyr  ! strictly valid
                       if(elecx*elxr == -1) elxr=-1.*elxr
               endif

               rjth(n) = sigped*((elecx+elxr)+ &
                              Wind_eastwards_ms1(n,m,l)*brad)/Sine_Dip_angle &
                              **2 - sighal*(Wind_southwards_ms1(n,m,l) &
                                    *brad-(elecy &
                                    +elyr)) &
                              /Sine_Dip_angle

               rjphi(n) = -sigped*(Wind_southwards_ms1(n,m,l) &
                                 *brad-(elecy+elyr)) &
                               - sighal*((elecx+elxr)+ &
                               Wind_eastwards_ms1(n,m,l)*brad)/Sine_Dip_angle
           else
               rjth(n)=jth(n)
               rjphi(n)=jphi(n)
           endif

           !g
           !c  **
           !c  **

           vionx(n,m,l) = Wind_southwards_ms1(n,m,l) &
                   +(jphi(n)*brad-bphi*jrad(n)) &
                   /(Electron_density_m3(n,m,l)*rvin(n)*ramin(n))

           viony(n,m,l) = Wind_eastwards_ms1(n,m,l) &
                   +(btheta*jrad(n)-jth(n)*brad) &
                   /(Electron_density_m3(n,m,l)*rvin(n)*ramin(n))

           dvx(n) = (Wind_southwards_ms1(nu,m,l)-Vx_1d_copy(nd))/2.0
           dvy(n) = (Wind_eastwards_ms1(nu,m,l)-Vy_1d_copy(nd))/2.0
           dt = (temp(nu)-temp(nd))/2.0
           dumut(n) = (umut(nu)-umut(nd))/2.0
           dkmkt = (km(nu)+kt(nu)-km(nd)-kt(nd))/2.0
           ddt = temp(nu) + temp(nd) - 2.0*temp(n)
           ddvx(n) = Wind_southwards_ms1(nu,m,l) + Vx_1d_copy(nd) &
                               - 2.0*Wind_southwards_ms1(n,m,l)
           ddvy(n) = Wind_eastwards_ms1(nu,m,l) + Vy_1d_copy(nd) - &
                          2.0*Wind_eastwards_ms1(n,m,l)

           !---------------------------------------------------------------------------
           !       energy equation
           !
           !       c1 and c2 form horizontal advection of eps
           !       c5 is vertical advection of eps
           !       c6 is solar heating and IR cooling
           !       c3 (first line) is vertical molecular and turbulent heat conduction
           !       c3 (3rd line) is vertical turbulent heat conduction due to
           !          adiabatic lapse rate
           !       c4 is vertical viscous drag term
           !       c10 (calculated above) is horizontal molecular and turbulent
           !          heat conduction
           !       c7 is ion drag and Joule heating term
           !---------------------------------------------------------------------------
           c1 = -Wind_southwards_ms1(n,m,l) &
                          *((epss(n)+GRAV*hts(n))-(eps(n,mn,l)+ &
                          GRAV*ht(n,mn,l)))/sx1
           c2 = -Wind_eastwards_ms1(n,m,l) &
                          *((epse(n)+GRAV*hte(n))-(epsw(n)+GRAV*htw(n))) &
                          /sx2
           c3 = GRAV*((km(n)+kt(n))*ddt+dt*dkmkt)/pressure(n)
           c3 = c3*temp(n)**0.71*tin(n)/(tin(n)**0.71*temp(n))
           c3 = c3 + (kt(nu)*scht(nu)-kt(nd)*scht(nd)) &
                *GRAV**2/(pressure(n)*cp(n)*2.0)

           c4 = GRAV*(Wind_southwards_ms1(n,m,l)*umut(n) &
                          *ddvx(n)+Wind_southwards_ms1(n,m,l) &
                          *dvx(n)*dumut(n)+umut(n)*dvx(n)**2 &
                          +Wind_eastwards_ms1(n,m,l) &
                          *umut(n)*ddvy(n) &
                          +Wind_eastwards_ms1(n,m,l)*dvy(n)*dumut(n) &
                          +umut(n)*dvy(n)**2)/pressure(n)
           c51 = om(n)*(eps(nu,m,l)-Eps_1d_copy(nd))/ &
                         (2.0*pressure(n))
           c52 = +om(n)/neutral_density_1d(n)
           c5 = c51 + c52

           !g------------------------------------------------------------------------------
           !g  new definitions for c77 to take into account the E field variations......
           !g------------------------------------------------------------------------------
           c7(n) = (rjth(n)*(elecx+elxr)+ &
                    rjphi(n)*(elecy+elyr))/neutral_density_1d(n)

           c77(n) = (rjth(n)*((elecx+elxr) &
                              +Wind_eastwards_ms1(n,m,l)*brad)+ &
                              rjphi(n)*((elecy+elyr)- &
                              Wind_southwards_ms1(n,m,l)*brad))/neutral_density_1d(n)


           if (efield_var == 1) then
              c7(n)=2.*c7(n)
              c77(n)=2.*c77(n)
           endif


           c7_3d(n,m,l) = c7(n)
           c77_3d(n,m,l) = c77(n)

           c8 = edep(n)

           !c  **
           !c  evaluate recombination heating
           !c  **
           source1 = 0.0
           source2 = 0.0

           if (n >= 7) then ! only higher altitudes

              ! this could be more efficient ***
              tr = temp(n)*1.2
              k8 = 2.82e-17-7.74e-18*(tr/300.0)+1.073e-18*(tr &
                /300.0)**2-5.17e-20*(tr/300.0)**3+9.65e-22*(tr &
                /300.0)**4

              !if (tr < 1700.0) k3 = 1.533e-18-5.92e-19*( &   lrm20130731
              !                   tr/300.0)+8.6e-20*(tr/300.0)**2
              !if (tr >= 1700.0) k3 = 2.73e-18-1.155e-18*( &
              !                    tr/300.0)+1.483e-19*(tr/300.0)**2

              if (tr < 1700.0) then
                  k3 = 1.533e-18-5.92e-19*(tr/300.0)+8.6e-20*(tr/300.0)**2
              else
                  k3 = 2.73e-18-1.155e-18*( &
                                  tr/300.0)+1.483e-19*(tr/300.0)**2
              endif

              !------------------------------------------------------
              ! evaluate energy from O+ recombination with N2
              ! source1 term for heating
              !------------------------------------------------------
              !factor = 13.0*1.66E-19/neutral_density_1d(n)  ! lrm20130722

              !source1 = O_plus_density_m3(n,m,l)*p33(n)*k3*13.0*1.66E-19/neutral_density_1d(n)  ! original
              !source1 = O_plus_density_m3(n,m,l)*p33(n)*k3*factor  ! lrm20130722
              source1 = O_plus_density_m3(n,m,l)*p33(n)*k3*4.0*1.66E-19/neutral_density_1d(n)  ! ctipe formula lrm20130807


              ! *****************lrm20130705*******************************************************
              !source1 = 0.0 !  ***************** ONLY FOR TESTING *********************************
              !source1 = source1 * .25 !  **** ONLY FOR TESTING *********************
              !**************************************************************************************

              !-----------------------------------------------
              ! evaluate energy from O+ recombination with O2
              !-----------------------------------------------
              !source2 = O_plus_density_m3(n,m,l)*p2(n)*k8*13.0*1.66E-19/neutral_density_1d(n) ! original
              !source2 = O_plus_density_m3(n,m,l)*p2(n)*k8*factor   ! lrm20130722
              source2 = O_plus_density_m3(n,m,l)*p2(n)*k8*8.55*1.66E-19/neutral_density_1d(n)   ! ctipe formula lrm20130807

           endif ! if n>= 7

           !--------------------------------------------------------
           !  replace w/ neutral heating rates read in from ipe 
           !--------------------------------------------------------
           if (useIPEheatingRates) then
                ! keep qir(n) and add to neutral heating rate (add all 7) & need to add NO cooling
                c6 = qir(n) + neutralHeatingRates(n,m,l)/2.0 ! Joules/kg/sec  lrm20131112
            else

                c6 = qir(n) + qeuv(n) + source1 + source2    

           endif
                                                       
           sum3(n) = c1 + c2 + c3 + c4 + c5 + c6 + c7(n) &
                               + c8 + c10(n)
           volume = (ht(nu,m,l)-ht(nd,m,l)) &
                              *R0**2*deltha*sth*delphi/8.0
           stfac(n) = volume*neutral_density_1d(n)
           c77(n) = c77(n)*stfac(n)


           !--------------------------------------------------------
           ! For checking values of source1, source2  lrm20130610
           !--------------------------------------------------------
           if (checkHeatingSources) then
              source1Array(n,m,l) = source1 
              source2Array(n,m,l) = source2
              !k3Array(n,m,l) = k3
              !c7Array(n,m,l) = c7(n)
           endif

        END DO looplevels1695 ! ------------------------------------------------------------



        !c  ************************************************************
        !
        !       momentum equation
        !
        !       a: in x-direction
        !       b: in y-direction
        !       sum1 and sum2 are in units of [m/s^2]
        !       a1/b1 are horizontal advection terms of velocity
        !       a6/b6 are vertical advection terms of velocity
        !       a4/b4 are geopotential terms
        !       a2/b2 are coriolis terms
        !       a3/b3 are vertical viscous drag terms
        !       a8/b8 are horizontal viscous drag terms
        !       a5/b5 are ion drag terms
        !
        !c  ************************************************************
        !DO n = 2 , 14
        DO n = 2 , n_levels-1
           nu = n + 1
           nd = n - 1
           a11 = -Wind_southwards_ms1(n,m,l)*(vxs(n) &
                           -Wind_southwards_ms1(n,mn,l))/sx1
           a12 = -Wind_eastwards_ms1(n,m,l) &
                           *(vxe(n)-vxw(n))/sx2
           a1 = a11 + a12
           a2 = (2.0*ANROT+Wind_eastwards_ms1(n,m,l)/sx8) &
                           *Wind_eastwards_ms1(n,m,l)*cth
           a3 = GRAV*(umut(n)*ddvx(n)+dumut(n)*dvx(n))/pressure(n)
           a4 = -GRAV*(hts(n)-ht(n,mn,l))/sx1
           a5 = (jphi(n)*brad-bphi*jrad(n))/neutral_density_1d(n)
           a6 = +om(n)*(Wind_southwards_ms1(nu,m,l) &
                          -Vx_1d_copy(nd))/(2.0*pressure(n))
           sum1(n) = a1 + a2 + a3 + a4 + a5 + a6 + a8(n)

           !--------------------------------------------------------
           ! For checking values of source1, source2  lrm20130717
           !--------------------------------------------------------
           if (checkHeatingSources) then
              a5Array(n,m,l) = a5
           endif

        END DO ! n = 2, n_levels-1


        !DO n = 2 , 14
        DO n = 2 , n_levels-1

           nu = n + 1
           nd = n - 1
           b11 = -Wind_southwards_ms1(n,m,l) &
                           *(vys(n)-Wind_eastwards_ms1(n,mn,l))/sx1
           b12 = -Wind_eastwards_ms1(n,m,l) &
                           *(vye(n)-vyw(n))/sx2
           b1 = b11 + b12
           b2 = -(2.0*ANROT+Wind_eastwards_ms1(n,m,l)/sx8) &
                          *Wind_southwards_ms1(n,m,l)*cth
           b3 = GRAV*(umut(n)*ddvy(n)+dumut(n)*dvy(n))/pressure(n)
           b4 = -GRAV*(hte(n)-htw(n))/sx2
           B5 = -(JTH(N)*BRAD-JRAD(N)*BTHETA)/neutral_density_1d(n)
           b6 = +om(n)*(Wind_eastwards_ms1(nu,m,l) &
                          -Vy_1d_copy(nd))/(2.0*pressure(n))
           sum2(n) = b1 + b2 + b3 + b4 + b5 + b6 + b8(n)

           !--------------------------------------------------------
           ! For checking values of source1, source2  lrm20130717
           !--------------------------------------------------------
            if (checkHeatingSources) then
               b5Array(n,m,l) = b5
            endif

         END DO ! 


         !DO n = 2 , 14
         !            Wind_southwards_ms1(n,m,l) = realTimeStep*sum1(n) &
         !                                  + Wind_southwards_ms1(n,m,l)
         !            Wind_eastwards_ms1(n,m,l) = realTimeStep*sum2(n) &
         !                                  + Wind_eastwards_ms1(n,m,l)
         !            eps(n,m,l) = realTimeStep*sum3(n) + eps(n,m,l)
         !END DO
         Wind_southwards_ms1(2:n_levels-1,m,l) = realTimeStep*sum1(2:n_levels-1) &
                                           + Wind_southwards_ms1(2:n_levels-1,m,l)
         Wind_eastwards_ms1(2:n_levels-1,m,l) = realTimeStep*sum2(2:n_levels-1) &
                                           + Wind_eastwards_ms1(2:n_levels-1,m,l)
         eps(2:n_levels-1,m,l) = realTimeStep*sum3(2:n_levels-1) + eps(2:n_levels-1,m,l)



         Wind_southwards_ms1(15,m,l) = Wind_southwards_ms1(14,m,l)
         Wind_eastwards_ms1(15,m,l) = Wind_eastwards_ms1(14,m,l)
         t14 = (-(Wind_southwards_ms1(14,m,l)**2 &
                 + Wind_eastwards_ms1(14,m,l)**2)/2.0+eps(14,m,l)) &
                 /cp(14)

         !t15 = t14
         !eps(15,m,l) = cp(15) &
         !              *t15 + (Wind_southwards_ms1(15,m,l)**2 &
         !              + Wind_eastwards_ms1(15,m,l)**2)/2.0
         eps(15,m,l) = cp(15) &
                       *t14 + (Wind_southwards_ms1(15,m,l)**2 &
                       + Wind_eastwards_ms1(15,m,l)**2)/2.0

         vx2(:,m) = Vx_1d_copy(:)
         vy2(:,m) = Vy_1d_copy(:)
         eps2(:,m) = Eps_1d_copy(:)

         IF ( l == 1 ) THEN
            vxl(:,m) = Vx_1d_copy(:)
            vyl(:,m) = Vy_1d_copy(:)
            epsl(:,m) = Eps_1d_copy(:)
         ENDIF


!g
!g  End of the latitude, longitude loops.....
!g
   END DO latloop1720
 END DO lonloop1740


 !---------------------------------------
 ! Calculate values at the poles.......
 !-------------------------------------------

 zz = m_North_Pole - 1

 !DO 1780 m = 1 , m_North_Pole , zz
 latpole1780 : DO m = 1 , m_North_Pole , zz
    fp1 = 1.0
    IF ( bb == 4. ) fp1 = 2.0

    IF ( m == 1 ) THEN
         mm = 1
         fp = -1.0
    ENDIF

    IF ( m == m_North_Pole ) THEN
         mm = 2
         fp = 1.0
    ENDIF

    mp = m - (fp*fp1) + 0.1

    !DO 1750 n = 1 , 15
    levelpole1750 : DO n = 1 , n_levels

       vpx(n,mm) = 0.0
       vpy(n,mm) = 0.0
       epsp(n,mm) = 0.0
       om1(n,m,1) = 0.0
       Electron_density_m3(n,m,1) = 0.0
       neutral_density_3d(n,m,1) = 0.0

       !DO 1745 l = 1 , 20
       lonpole1745 : DO l = 1 , n_lons
          phir = (FLOAT(l)-1.0)*PI/10.0
          cphi = COS(phir)
          sphi = SIN(phir)
          vpx(n,mm) = vpx(n,mm) + Wind_southwards_ms1(n,mp,l) &
                                 *cphi - Wind_eastwards_ms1(n,mp,l)*sphi*fp
          vpy(n,mm) = vpy(n,mm) + Wind_southwards_ms1(n,mp,l)*sphi*fp + &
                                 Wind_eastwards_ms1(n,mp,l)*cphi
          epsp(n,mm) = epsp(n,mm) + eps(n,mp,l)
          om1(n,m,1) = om1(n,m,1) + om1(n,mp,l)
          Electron_density_m3(n,m,1) = Electron_density_m3(n,m,1) + Electron_density_m3(n,mp,l)
          neutral_density_3d(n,m,1) = neutral_density_3d(n,m,1) + neutral_density_3d(n,mp,l)

          !IF ( l == 20 ) THEN
          IF ( l == n_lons ) THEN
             vpx(n,mm) = vpx(n,mm)/20.0
             vpy(n,mm) = vpy(n,mm)/20.0
             epsp(n,mm) = epsp(n,mm)/20.0
             om1(n,m,1) = om1(n,m,1)/20.
             Electron_density_m3(n,m,1) = Electron_density_m3(n,m,1)/20.
             neutral_density_3d(n,m,1) = neutral_density_3d(n,m,1)/20.
          ENDIF


 !1745             CONTINUE
       END DO lonpole1745
 !1750          CONTINUE
    END DO levelpole1750

    !DO 1760 l = 1 , 20
    lonpole1760 : DO l = 1 , n_lons
       phir = (FLOAT(l)-1.0)*PI/10.0
                  cphi = COS(phir)
                  sphi = SIN(phir)

                  !DO 1755 n = 2 , 15
                  DO n = 1 , n_levels ! lm20110922  Start at 1
                     Wind_southwards_ms1(n,m,l) = vpx(n,mm)*cphi + vpy(n,mm)*sphi*fp
                     Wind_eastwards_ms1(n,m,l) = vpy(n,mm)*cphi - vpx(n,mm)*sphi*fp
                     eps(n,m,l) = epsp(n,mm)
                     om1(n,m,l) = om1(n,m,1)
                     Electron_density_m3(n,m,l) = Electron_density_m3(n,m,1)
                     neutral_density_3d(n,m,l) = neutral_density_3d(n,m,1)
 !1755             CONTINUE
                  END DO

 !1760          CONTINUE
    END DO lonpole1760

 !1780       CONTINUE
 END DO latpole1780

 IF ( bb == 4.0 ) THEN

      !DO 1800 l = 1 , 20
      lonbb1800 : DO l = 1 , n_lons
         m = ixx
         DO 1785 n = 1 , 15
                     Wind_southwards_ms1(n,m,l) = (Wind_southwards_ms1(n,m+1,l)+Wind_southwards_ms1(n,m-1,l))/2.0
                     Wind_eastwards_ms1(n,m,l) = (Wind_eastwards_ms1(n,m+1,l)+Wind_eastwards_ms1(n,m-1,l))/2.0
                     eps(n,m,l) = (eps(n,m+1,l)+eps(n,m-1,l))/2.0
                     om1(n,m,l) = (om1(n,m+1,l)+om1(n,m-1,l))/2.0
                     Electron_density_m3(n,m,l) = (Electron_density_m3(n,m+1,l)+Electron_density_m3(n,m-1,l))/2.0
                     neutral_density_3d(n,m,l) = (neutral_density_3d(n,m+1,l)+neutral_density_3d(n,m-1,l))/2.0
 1785             CONTINUE
                  m = 2
                  DO 1790 n = 1 , 15
                     Wind_southwards_ms1(n,m,l) = (Wind_southwards_ms1(n,m+1,l)+Wind_southwards_ms1(n,m-1,l))/2.0
                     Wind_eastwards_ms1(n,m,l) = (Wind_eastwards_ms1(n,m+1,l)+Wind_eastwards_ms1(n,m-1,l))/2.0
                     eps(n,m,l) = (eps(n,m+1,l)+eps(n,m-1,l))/2.0
                     om1(n,m,l) = (om1(n,m+1,l)+om1(n,m-1,l))/2.0
                     Electron_density_m3(n,m,l) = (Electron_density_m3(n,m+1,l)+Electron_density_m3(n,m-1,l))/2.0  !lm20110922
                     neutral_density_3d(n,m,l) = (neutral_density_3d(n,m+1,l)+neutral_density_3d(n,m-1,l))/2.0
 1790             CONTINUE
 1800          CONTINUE
         END DO lonbb1800
 ENDIF

!
! check if smoothing needed this step
!
            sw_smoothing_this_time_step = .false.

            !print *,'GT_thermosphere : nn_smoothing_counter, i_smoothing_frequency = ',nn_smoothing_counter, i_smoothing_frequency

            IF ( nn_smoothing_counter == i_smoothing_frequency ) THEN

                 print *,'GT_thermosphere : SMOOTHING .....................................'

                 sw_smoothing_this_time_step = .true.

                 nn_smoothing_counter = 0

                 DO 1805 l = 1 , 20
                     DO 1804 m = 1 , m_North_Pole

                        CALL SPECIFIC_HEAT(psao(1,m,l), psmo(1,m,l), psmn(1,m,l), cp)

                        DO n = 2 , 15

                           nd = n - 1
                           Temperature_K(n,m,l) = (-(Wind_southwards_ms1(n,m,l)**2+Wind_eastwards_ms1(n,m,l)**2)/ &
                                        2.0+eps(n,m,l))/cp(n)
                           scht(n) = (GSCON*Temperature_K(n,m,l)) &
                                     /(rmt(n,m,l)*GRAV)
                           htt = ht(n,m,l)
                           ndd = n - 2
                           IF ( n <= 2 ) THEN
                              ht(n,m,l) = ht(nd,m,l) + scht(nd)
                           ELSE
                              ht(n,m,l) = ht(nd,m,l) + scht(nd) &
                                 + (scht(n)-scht(ndd)) &
                                 /4.0 + (scht(n)+scht(ndd)-2.0*scht(nd)) &
                                 /6.0
                           ENDIF
                           wvz(n,m,l) = (ht(n,m,l)-htt) &
                                        /realTimeStep - om1(n,m,l)/neutral_density_3d(n,m,l)/GRAV

                       END DO ! n = 2 , 15



                  !------------------------------
                  ! Now calculate wvz for  n = 1
                  !------------------------------
                  !print *,'GT_thermosphere : m, l, ht(1,m,l) = ', m, l, ht(1,m,l)
                  !print *,'GT_thermosphere : htold(m,l), realTimeStep = ', htold(m,l), realTimeStep
                  !print *,'GT_thermosphere : om1(1,m,l), neutral_density_3d(1,m,l) = ', om1(1,m,l), neutral_density_3d(1,m,l)

                  IF (m > 84 .OR. m < 8) THEN 
                      wvz(1,m,l) = 0.
                  ELSE
                      wvz(1,m,l) = (ht(1,m,l)-htold(m,l))/realTimeStep - &
                                    om1(1,m,l)/neutral_density_3d(1,m,l)/GRAV
                  ENDIF

                  !wvz(1,m,l) = (ht(1,m,l)-htold(m,l))/realTimeStep - &
                  !              om1(1,m,l)/neutral_density_3d(1,m,l)/GRAV

                  !IF (m > 84 .OR. m < 8) wvz(1,m,l)=0.  up higher now to get rid of divide by 0


 1804                CONTINUE
 1805             CONTINUE

                  j = m_North_Pole
                  nmin = 2
                  nmax = 15
                  mmin = 1
                  mmax = m_North_Pole

                  s = 0.25
                  CALL Smooth_VnX_in_X(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_VnX_in_X(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_VnX_in_Y(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_VnX_in_Y(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_VnY_in_X(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_VnY_in_X(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_VnY_in_Y(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_VnY_in_Y(j,Wind_southwards_ms1,Wind_eastwards_ms1,deltha,delphi,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_Tn_in_X(j,Temperature_K,deltha,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_Tn_in_X(j,Temperature_K,deltha,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_Tn_in_Y(j,Temperature_K,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_Tn_in_Y(j,Temperature_K,s,nmin,nmax,mmin,mmax)
!
! new major species composition smoothing added 7th Feb 2008 by TJFR
!
                  s = 0.25
                  CALL Smooth_Tn_in_X(j,psao,deltha,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_Tn_in_X(j,psao,deltha,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_Tn_in_Y(j,psao,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_Tn_in_Y(j,psao,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_Tn_in_X(j,psmo,deltha,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_Tn_in_X(j,psmo,deltha,s,nmin,nmax,mmin,mmax)

                  s = 0.25
                  CALL Smooth_Tn_in_Y(j,psmo,s,nmin,nmax,mmin,mmax)

                  s = -0.25
                  CALL Smooth_Tn_in_Y(j,psmo,s,nmin,nmax,mmin,mmax)

                  !print *,'GT_thermosphere : after all smoothing function calls....'
!
! end of inserts TJFR
!
        DO 1810 l = 1 , 20
           DO 1808 m = 1 , m_North_Pole

              !CALL SPECIFIC_HEAT(psao(1,m,l),psmo(1,m,l),psmn(1,m,l),rmt(1,m,l),cp)
              CALL SPECIFIC_HEAT(psao(1,m,l), psmo(1,m,l), psmn(1,m,l), cp)

              DO 1806 n = 2 , 15
                 eps(n,m,l) = cp(n)*Temperature_K(n,m,l) + &
                              (Wind_southwards_ms1(n,m,l)**2+Wind_eastwards_ms1(n,m,l)**2) / 2.0
 1806         CONTINUE

 1808   CONTINUE
 1810   CONTINUE

        ENDIF    ! smoothing endif

!
! call the composition equation if required
!
       IF ( nn_composition_counter  ==  i_neutral_composition_calling_frequency ) THEN

         CALL neutral_composition_equation( &
                   m_North_Pole, pressure, i_neutral_composition_calling_frequency, &
                   f107, Wind_southwards_ms1, Wind_eastwards_ms1, ht, rmt, om1, Temperature_K, &
                   psao, psmo, psmn)

         nn_composition_counter = 0

       ENDIF




            DO 1840 l = 1 , 20
               DO 1820 m = 1 , m_North_Pole
                    CALL SPECIFIC_HEAT(psao(1,m,l), psmo(1,m,l), psmn(1,m,l), cp)
                    DO 1815 n = 2 , 15
                       nd = n - 1
                       temp(n) = (-(Wind_southwards_ms1(n,m,l)**2 + Wind_eastwards_ms1(n,m,l)**2)/2.0+eps(n,m &
                               ,l))/cp(n)
                       IF ( temp(n) < 0. ) WRITE(6,*) &
                          'Vx, Vy, Eps, cp ' , Wind_southwards_ms1(n,m,l) , Wind_eastwards_ms1(n,m,l) , &
                          eps(n,m,l) , cp(n)
                       IF ( temp(n) < 0. ) WRITE(6,*) 'n,m,l = ' , n , &
                          m , l
                       Temperature_K(n,m,l) = temp(n)
                       neutral_density_1d(n) = (pressure(n)*rmt(n,m,l))/(GSCON*temp(n))
                       scht(n) = (GSCON*temp(n))/(rmt(n,m,l)*GRAV)
!
! check for negative scht
!
                       IF ( scht(n) < 0.0 ) then
                          WRITE(6,99019) scht(n) , rmt(n,m,l) , temp(n)
                          stop
                       endif
99019 FORMAT (' negative scht:',e13.6,'  rmt: ',f9.3,' temp: ',f15.2)

                     htt = ht(n,m,l)
                     ndd = n - 2
                     IF ( n <= 2 ) THEN
                        ht(n,m,l) = ht(nd,m,l) + scht(nd)
                     ELSE
                        ht(n,m,l) = ht(nd,m,l) + scht(nd) &
                                    + (scht(n)-scht(ndd)) &
                                    /4.0 + (scht(n)+scht(ndd) &
                                    -2.0*scht(nd))/6.0
                     ENDIF

                     !print *,'GT_thermosphere : realTimeStep = ',realTimeStep
                     !print *,'GT_thermosphere : om1(n,m,l) = ',om1(n,m,l)
                     !print *,'GT_thermosphere : neutral_density_1d(n) = ',neutral_density_1d(n)


                     IF ( .not. sw_smoothing_this_time_step) &
                          wvz(n,m,l) = (ht(n,m,l)-htt) &
                          /realTimeStep - om1(n,m,l)/neutral_density_1d(n)/GRAV

 1815             CONTINUE

            IF ( .not. sw_smoothing_this_time_step) &
                              wvz(1,m,l) = (ht(1,m,l)-htold(m,l))/ &
                              realTimeStep - om1(1,m,l)/neutral_density_1d(1)/GRAV

            IF (m > 84.OR.m < 8) wvz(1,m,l)=0.

 1820          CONTINUE
 1840       CONTINUE


         !---------------------------------
         ! Write out the GT netcdf file?
         !---------------------------------
         if (idump_GT  ==  1) then

             print *,'GT_thermo : writing out netcdf for ',universal_time_seconds
             !call write_gt_netcdf_history(GT_output_dataset)
             ! Pass in the lats, lons
             call write_gt_netcdf_history(GT_output_dataset)

         endif


        !----------------------------------------------
        ! Write out the graphics data information?
        !----------------------------------------------
        !IF ( i_no_day >= i_graphics_out_start ) THEN

           !Universal_Time_hours = Universal_Time_seconds / 3600.
           Universal_Time_hours = DMOD(Universal_Time_seconds, secondsInDay)

        !  if (iout_high == ipint_high) then
        
            iout_high = iout_high + 1

             !---------------------------------------------
             ! Write out to graphics file every 15 minutes
             !---------------------------------------------
             ggg = MOD(iout_high,15)
             if (ggg == 0) then 


                irec_number = irec_number + 1


                !-----------------------------
	        ! Write to the graphics file
                !-----------------------------
                call write_gt_netcdf_graphics(graphics_file, irec_number, &
                                              electron_density_m3, Universal_Time_hours)

                !----------------------------------------------------------
                ! For write values of aeuv, source1, source2 for checking lrm20130610
                ! or for writing psmn, rmt, neutral density, oplus lrm20130624
                !----------------------------------------------------------
                if (checkHeatingSources) then
                    !print *,'gt_thermo : writing out aeuv, source1, source2...............'
                    !write(aeuvUnit,*) aeuv
                    write(source1Unit,*) source1Array
                    write(source2Unit,*) source2Array
                    !write(psmnUnit,*) psmn
                    !write(rmtUnit,*) rmt
                    !write(neutral_density_3dUnit,*) neutral_density_3d
                    write(OplusUnit,*) O_plus_density_m3
                    !write(a5Unit,*) a5Array
                    !write(b5Unit,*) b5Array
                    !write(c7Unit,*) c7Array
                    !write(k3Unit,*) k3Array
                endif

             endif  ! ggg == 0


         !       iout_high = 0

         !   endif ! (iout_high == ipint_high)

        !ENDIF ! ( i_no_day >= i_graphics_out_start )
        !-----------------------------------------------------------------------------------------


       wind_southwards_ms1_copy(:,:,:) = wind_southwards_ms1(:,:,:)
       wind_eastwards_ms1_copy(:,:,:) = wind_eastwards_ms1(:,:,:)
       wvz_copy(:,:,:) = wvz(:,:,:)
       rmt_copy(:,:,:) = rmt(:,:,:)
       Temperature_K_copy(:,:,:) = Temperature_K(:,:,:)
       ht_copy(:,:,:) = ht(:,:,:)

       do n = 1 , 15
       do m = 1 , 91
       do l = 1 , 20
          prsall = 1.0376*EXP(1.-FLOAT(n))
          rnt = prsall/(BZ*temperature_K(n,m,l))
          O_density_copy(n,m,l) = psao(n,m,l)*rnt*rmt(n,m,l)/16.
          O2_density_copy(n,m,l) = psmo(n,m,l)*rnt*rmt(n,m,l)/32.
          N2_density_copy(n,m,l) = psmn(n,m,l)*rnt*rmt(n,m,l)/28.
       enddo
       enddo
       enddo


RETURN

END SUBROUTINE GT_thermosphere


!==============================================================================


SUBROUTINE GT_thermosphere_INIT( &
                 GT_input_dataset, &  ! startup file
                 GT_output_dataset, & ! output file
                 nday, &
                 Universal_Time_seconds, &
                 solar_declination_angle_radians, &
                 hough11 , hough22 , hough23 , hough24 , hough25, &
                 amplitude, &
                 tidalPhase, &
                 wind_southwards_ms1_copy, &  ! read from startup file
                 wind_eastwards_ms1_copy, & ! read from startup file
                 wvz_copy, &  ! upward nuetral wind
                 rmt_copy, & ! composition 
                 Temperature_K_copy, & ! nuetral temperature in 3D pressure coords
                 ht_copy, & ! height at each grid point in 3D
                 Electron_density_m3, & ! ionospheric electron density interpolated to pressure heights
                 O_plus_density_m3, &   ! ionospheric O+ density interpolated to pressure heights
                 NO_plus_density_m3, &  ! ionospheric NO+ density interpolated to pressure heights
                 O2_plus_density_m3, &  ! ionospheric O2+ density interpolated to pressure heights
                 Ti_Oplus_K )           ! ionospheric Ti O+ temperature interpolated to pressure heights

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: GT_input_dataset, GT_output_dataset
INTEGER, INTENT(IN) :: nday
REAL(kind=8), INTENT(IN) :: Universal_Time_seconds
REAL*8, INTENT(OUT) :: solar_declination_angle_radians
REAL*8, INTENT(OUT) :: hough11(181) , hough22(181) , hough23(181) , hough24(181), hough25(181)

TYPE(amplitudeType), INTENT(IN) :: amplitude
TYPE(tidalPhaseType), INTENT(IN) :: tidalPhase

REAL(kind=8), INTENT(OUT) :: wind_southwards_ms1_copy(n_levels, n_lats, n_lons)
REAL(kind=8), INTENT(OUT) :: wind_eastwards_ms1_copy(n_levels, n_lats, n_lons)
REAL(kind=8), INTENT(OUT) :: wvz_copy(n_levels, n_lats, n_lons)
REAL(kind=8), INTENT(OUT) :: rmt_copy(n_levels, n_lats, n_lons)
REAL(kind=8), INTENT(OUT) :: Temperature_K_copy(n_levels, n_lats, n_lons)
REAL(kind=8), INTENT(OUT) :: ht_copy(n_levels, n_lats, n_lons)

REAL*8, INTENT(OUT) :: Ti_Oplus_K(n_levels, n_lats, n_lons)
REAL*8, INTENT(OUT) :: O2_plus_density_m3(n_levels, n_lats, n_lons)
REAL*8, INTENT(OUT) :: Electron_density_m3(n_levels, n_lats, n_lons)
REAL*8, INTENT(OUT) :: O_plus_density_m3(n_levels, n_lats, n_lons)
REAL*8, INTENT(OUT) :: NO_plus_density_m3(n_levels, n_lats, n_lons)

             

! Local variables -----------------------------------------------------



!-------------
! Parameters 
!-------------
REAL*8, parameter :: GRAV = 9.5
REAL*8, parameter :: GSCON = 8.3141E+03
REAL*8, parameter :: p0 = 1.0376
REAL*8, parameter :: R0 = 6.370E06
!REAL*8, parameter  :: PI = 3.14159



REAL*8 :: zn
REAL*8 :: ty
REAL*8 :: cp(n_levels)

INTEGER :: l , m
INTEGER :: n , nd , ndd
INTEGER :: nnstrt

REAL*8 :: temp(n_levels) , scht(n_levels)

REAL*8 :: exsize
REAL*8 :: realTimeStep

LOGICAL, parameter :: debugInit = .FALSE.


! BEGIN CODE ========================================================================
!------------------------
! Initialize  variables - lrm
!------------------------
!ilenp1 = 0 not needed lrm20121108
!filer = 'xxxxxx'  not needed lrm20121108
zn = 0.
ty = 0.
cp = 0.
l = 0
m = 0
n = 0
nd = 0 
ndd = 0
nnstrt = 0  ! THIS SHOULD REALLY BE PASSED IN???**********

temp = 0.

scht = 0.  ! THIS HAS TO BE DONE TO GET THE SAME SCALE HEIGHTS AS THE 
                 ! CALCULATED ON RES-THERMO-01 W/ gtONLY code *********** lrm

exsize = 0.

realTimeStep = 10.   !!!! this needs sorting - set here for initial Tides
                           !!!! but the tides all need looking at.....


wind_southwards_ms1_copy = 0.
wind_eastwards_ms1_copy = 0.
wvz_copy = 0.
rmt_copy = 0.
Temperature_K_copy = 0.
ht_copy = 0.


m_North_Pole = n_lats
m_equator = (m_North_Pole + 1) / 2

!exsize = 0.0  - initializing it above
examp = 1.
fden = 1. + exsize
ixx = 90
mj1 = 45
mj2 = ixx - mj1
mjk = mj1
mj1 = mj1*fden
mj2 = ixx - mj1

!----------------
! input datasets
!----------------
!ilenp1 = INDEX(GT_input_dataset,' ')
WRITE (6,*) 'GT_thermosphere_INIT :     '
!filer = GT_input_dataset(1:ilenp1-1)
WRITE (6,*) 'GT_thermosphere_INIT : INPUT DATA SET :  ', GT_input_dataset

!----------------
! output datasets
!----------------
!ilenp1 = INDEX(GT_output_dataset,' ')
WRITE (6,*) 'GT_thermosphere_INIT :   '

!filer = GT_output_dataset(1:ilenp1-1)
WRITE (6,*) 'GT_thermosphere_INIT : OUTPUT DATA SET :  ', GT_output_dataset

!-------------------------------
! Counters for graphics output
!-------------------------------
irec_number = 0  !  for writing to the graphics file
iout_high = 0 ! counter for # of thermospheric calls

!-------------------------------
! Graphics output file name 
! Changed graphics_file name creation lrm20121108
!-------------------------------
graphics_file = trim(GT_output_dataset)//'.graphics.nc'
call create_gt_netcdf_graphics(graphics_file)

!-----------------------------------------------
!- here we need code to calculate ty from nday
!-----------------------------------------------
ty = (nday+15.5)*12./365.
IF ( ty > 12.0 ) ty = ty - 12.0

solar_declination_angle_radians = ATAN(0.434*SIN(PI/6.0*(ty-3.17)))

!-----------------------------------------------
! Read in the GT_thermosphere startup data.....
!-----------------------------------------------
!call read_gt_netcdf_history(GT_input_dataset)


!-----------------------------------------------
! Read in the GT and GIP startup data .....
!-----------------------------------------------
call read_gtgip_netcdf_history(GT_input_dataset, & 
                               Electron_density_m3, O_plus_density_m3, &
                               NO_plus_density_m3, O2_plus_density_m3, Ti_Oplus_K )




do n = 1 , n_levels
   zn = n - 1
   pressure(n) = P0*EXP(-zn)
enddo


       kt(1) =  3.310000000000000E-005
       kt(2) =  2.190000000000000E-005
       kt(3) =  1.030000000000000E-005
       kt(4) =  8.320000000000000E-006
       kt(5) =  4.760000000000000E-006
       kt(6) =  2.730000000000000E-006
       kt(7) =  4.960000000000000E-007
       kt(8) =  4.920000000000000E-008
       kt(9) =  0.000000000000000E+000
       kt(10) =  0.000000000000000E+000
       kt(11) =  0.000000000000000E+000
       kt(12) =  0.000000000000000E+000
       kt(13) =  0.000000000000000E+000
       kt(14) =  0.000000000000000E+000
       kt(15) =  0.000000000000000E+000


       km(1) =  3.169999991996519E-006
       km(2) =  3.059999992274243E-006
       km(3) =  2.979999992476223E-006
       km(4) =  3.019999992375233E-006
       km(5) =  2.979999992476223E-006
       km(6) =  2.839999992829689E-006
       km(7) =  2.579999993486126E-006
       km(8) =  2.199999994445534E-006
       km(9) =  1.879999995253456E-006
       km(10) =  1.759999995556427E-006
       km(11) =  1.689999995733160E-006
       km(12) =  1.679999995758408E-006
       km(13) =  1.679999995758408E-006
       km(14) =  1.679999995758408E-006
       km(15) =  1.679999995758408E-006

       umut(1)=  6.342924037227632E-008
       umut(2)=  4.286065225067180E-008
       umut(3)=  2.141756411649758E-008
       umut(4)=  1.779472704331060E-008
       umut(5)=  1.113431621574501E-008
       umut(6)=  7.227299305048900E-009
       umut(7)=  2.767096029057651E-009
       umut(8)=  1.657960851574864E-009
       umut(9)=  1.353052747977312E-009
       umut(10)=  1.200138269414521E-009
       umut(11)=  1.080667216327646E-009
       umut(12)=  9.863601239729971E-010
       umut(13)=  9.176341987897271E-010
       umut(14)=  8.732396318664023E-010
       umut(15)=  8.478973454830954E-010

       tin(1)=   192.000000000000
       tin(2)=   181.000000000000
       tin(3)=   183.066541580778
       tin(4)=   187.454914479433
       tin(5)=   206.389014854052
       tin(6)=   248.248965115296
       tin(7)=   339.223960473209
       tin(8)=   501.833599617990
       tin(9)=   692.035982522604
       tin(10)=   832.622041865981
       tin(11)=   917.233854439876
       tin(12)=   961.894411934840
       tin(13)=   982.465096934422
       tin(14)=   991.512738873483
       tin(15)=   995.601858538608

call read_in_and_normalise_Hough_modes(hough11, &
                                              hough22, &
                                              hough23, &
                                              hough24, &
                                              hough25)
DO m = 1, n_lats
   temp0(m) = 192.
   vy0(m) = 0.
   dh0(m) = 0.
ENDDO

temp0av = 210.


;print *,'GT_thermosphere_INIT : nnstrt = ',nnstrt

DO 1600 l = 1 , n_lons

   DO 1500 m = 1 , m_North_Pole

             !--------------------------------------------------------------------------
             ! new lower boundary insert:
             ! since tides are at level 1 need to do these calculations before.
             ! IMW August 1995
             !--------------------------------------------------------------------------
             CALL TIDES(nnstrt, m, l, hough11, hough22, hough23, &
                        hough24, hough25, &
                        amplitude, &
                        ht, &
                        wind_southwards_ms1, wind_eastwards_ms1, &
                        realTimeStep, rmt, temp, &
                        tidalPhase, &
                        temp0, vy0, temp0av, dh0)

             CALL SPECIFIC_HEAT(psao(1,m,l), psmo(1,m,l), psmn(1,m,l), cp)

             DO n = 2 , n_levels
                nd = n - 1
                temp(n) = (-(wind_southwards_ms1(n,m,l)**2+wind_eastwards_ms1(n,m,l)**2)/2.0+eps(n,m,l))/cp(n)
                scht(n) = (GSCON*temp(n))/(rmt(n,m,l)*GRAV)
                ndd = n - 2
                IF ( n <= 2 ) THEN
                   ht(n,m,l) = ht(nd,m,l) + scht(nd)
                ELSE
                   ht(n,m,l) = ht(nd,m,l) + scht(nd) &
                              + (scht(n)-scht(ndd)) &
                              /4.0 + (scht(n)+scht(ndd)-2.0*scht(nd))/6.0
                ENDIF

             END DO ! n = 2, n_levels

             !if (m == 51) print *,'GT_thermosphere_INIT : scht(1), cp(1) = ', scht(1), cp(1)

 1500    CONTINUE

 1600 CONTINUE


 if (debugINIT) then 
       print *,'GT_thermosphere_INIT : ****************************'
       print *,'GT_thermosphere_INIT : scht = ',scht
       print *,'GT_thermosphere_INIT : rmt(:,51,11) = ', rmt(:,51,11)
       print *,'GT_thermosphere_INIT : eps(:,51,11) = ', eps(:,51,11)
       print *,'GT_thermosphere_INIT : '
       print *,'GT_thermosphere_INIT : wind_southwards_ms1(:,51,11) = ',wind_southwards_ms1(:,51,11)
       print *,'GT_thermosphere_INIT : wind_eastwards_ms1(:,51,11) = ',wind_eastwards_ms1(:,51,11)
       print *,'GT_thermosphere_INIT : '
       print *,'GT_thermosphere_INIT : wind_southwards_ms1(1,51,:) = ',wind_southwards_ms1(1,51,:)
       print *,'GT_thermosphere_INIT : wind_eastwards_ms1(1,51,:) = ',wind_eastwards_ms1(1,51,:)
       print *,'GT_thermosphere_INIT : ht(:,1,1) = ', ht(:,1,1)
       print *,'GT_thermosphere_INIT : ht(:,51,11) = ', ht(:,51,11)
       print *,'GT_thermosphere_INIT : ht(:,91,20) = ', ht(:,91,20)
       print *,'GT_thermosphere_INIT : ****************************'
  endif


       wind_southwards_ms1_copy(:,:,:) = wind_southwards_ms1(:,:,:)
       wind_eastwards_ms1_copy(:,:,:) = wind_eastwards_ms1(:,:,:)
       wvz_copy(:,:,:) = wvz(:,:,:)
       rmt_copy(:,:,:) = rmt(:,:,:)
       Temperature_K_copy(:,:,:) = Temperature_K(:,:,:)
       ht_copy(:,:,:) = ht(:,:,:)


RETURN

END SUBROUTINE GT_thermosphere_INIT


! ==============================================================================

subroutine GEO2MG_apex(Magnetic_longitude_degrees, &
                       Magnetic_latitude_degrees, &
                       input_geo_longitude_degrees, &
                       input_geo_latitude_degrees, &
                       output_mag_longitude_degrees, &
                       output_mag_latitude_degrees)

      implicit none

      REAL*8 Magnetic_latitude_degrees(91,20)
      REAL*8 Magnetic_longitude_degrees(91,20)
      REAL*8 input_geo_longitude_degrees
      REAL*8 input_geo_latitude_degrees
      REAL*8 output_mag_longitude_degrees
      REAL*8 output_mag_latitude_degrees

      !----------------------------
      ! Local variables
      !----------------------------
      INTEGER ilon , ilat
      INTEGER ilon_west , ilon_east , ispecial
      INTEGER ilat_north , ilat_south

      REAL*8 geo_lon_deg(20)
      REAL*8 geo_lat_deg(91)
      REAL*8 x2_minus_x1
      REAL*8 x2_minus_x
      REAL*8 x_minus_x1
      REAL*8 y2_minus_y1
      REAL*8 y2_minus_y
      REAL*8 y_minus_y1
      REAL*8 x2x1_times_y2y1
      REAL*8 x2x_times_y2y
      REAL*8 xx1_times_y2y
      REAL*8 x2x_times_yy1
      REAL*8 xx1_times_yy1
      REAL*8 Q11 , Q12 , Q21 , Q22

! Begin code ======================================================

!------------------------
! Initialize variables
!------------------------

ilon = 0
ilat = 0
ilon_west  = 0
ilon_east  = 0
ispecial = 0
ilat_north  = 0
ilat_south = 0


geo_lon_deg = 0
geo_lat_deg = 0
x2_minus_x1 = 0
x2_minus_x = 0
x_minus_x1 = 0
y2_minus_y1 = 0
y2_minus_y = 0
y_minus_y1 = 0
x2x1_times_y2y1 = 0
x2x_times_y2y = 0
xx1_times_y2y = 0
x2x_times_yy1 = 0
xx1_times_yy1 = 0
Q11 = 0
Q12 = 0 
Q21 = 0 
Q22 = 0



      do ilon = 1 , 20
        geo_lon_deg(ilon) = (ilon - 1) * 18.
        if (geo_lon_deg(ilon) > input_geo_longitude_degrees) then
          ilon_east = ilon
          ilon_west = ilon - 1
          ispecial = 0
          if(ilon_west  ==  0) then
          ilon_west = 20
          ispecial = 2
          endif
          goto 2000
        endif
      enddo

        ilon_west = 20
        ilon_east = 1
        ispecial = 1

 2000 continue

      if (ispecial  ==  1) geo_lon_deg(ilon_east) = geo_lon_deg(ilon_east) + 360.

      do ilat = 1 , 91
         geo_lat_deg(ilat) = (ilat - 46) * 2.
         if (geo_lat_deg(ilat) > input_geo_latitude_degrees) then
             ilat_north = ilat
             ilat_south = ilat - 1
             goto 3000
         endif
      enddo

 3000 continue

!     write(6,*) input_geo_longitude_degrees , input_geo_latitude_degrees , &
!                geo_lon_deg(ilon_west) , geo_lon_deg(ilon_east), &
!                geo_lat_deg(ilat_north) , geo_lat_deg(ilat_south)

      x2_minus_x1 = geo_lon_deg(ilon_east) - geo_lon_deg(ilon_west)
!     if ( x2_minus_x1 < 0.0) x2_minus_x1 = x2_minus_x1 + 360.

      x2_minus_x = geo_lon_deg(ilon_east) - input_geo_longitude_degrees
!     if ( x2_minus_x < 0.0) x2_minus_x = x2_minus_x + 360.

      x_minus_x1 = input_geo_longitude_degrees - geo_lon_deg(ilon_west)
!     if ( x_minus_x1 < 0.0) x_minus_x1 = x_minus_x1 + 360.

      y2_minus_y1 = geo_lat_deg(ilat_north) - geo_lat_deg(ilat_south)
      y2_minus_y = geo_lat_deg(ilat_north) - input_geo_latitude_degrees
      y_minus_y1 = input_geo_latitude_degrees - geo_lat_deg(ilat_south)

      x2x1_times_y2y1 = x2_minus_x1 * y2_minus_y1
      x2x_times_y2y = x2_minus_x*y2_minus_y
      xx1_times_y2y = x_minus_x1*y2_minus_y
      x2x_times_yy1 = x2_minus_x*y_minus_y1
      xx1_times_yy1 = x_minus_x1*y_minus_y1

      Q11 = Magnetic_longitude_degrees(ilat_south,ilon_west)
      Q12 = Magnetic_longitude_degrees(ilat_north,ilon_west)
      Q21 = Magnetic_longitude_degrees(ilat_south,ilon_east)
      Q22 = Magnetic_longitude_degrees(ilat_north,ilon_east)

      output_mag_longitude_degrees = Q11*x2x_times_y2y/x2x1_times_y2y1 + &
                                     Q21*xx1_times_y2y/x2x1_times_y2y1 + &
                                     Q12*x2x_times_yy1/x2x1_times_y2y1 + &
                                     Q22*xx1_times_yy1/x2x1_times_y2y1

!     write(6,*) Q11, Q12, Q21, Q22
!     write(6,*) 'wobba ', output_mag_longitude_degrees 

      Q11 = Magnetic_latitude_degrees(ilat_south,ilon_west)
      Q12 = Magnetic_latitude_degrees(ilat_north,ilon_west)
      Q21 = Magnetic_latitude_degrees(ilat_south,ilon_east)
      Q22 = Magnetic_latitude_degrees(ilat_north,ilon_east)

      output_mag_latitude_degrees = Q11*x2x_times_y2y/x2x1_times_y2y1 + &
                                    Q21*xx1_times_y2y/x2x1_times_y2y1 + &
                                    Q12*x2x_times_yy1/x2x1_times_y2y1 + &
                                    Q22*xx1_times_yy1/x2x1_times_y2y1

!     write(6,*) Q11, Q12, Q21, Q22
!     write(6,*) 'wobba2 ',output_mag_latitude_degrees

      

RETURN

END SUBROUTINE GEO2MG_apex





subroutine create_gt_netcdf_graphics(filename)

      include 'netcdf.inc'

      character(len=*), intent(in) :: filename
      integer :: ncid, istat
      ! integer :: id NOT USED
      integer :: id_n_levels , id_n_lats , id_n_lons
      integer :: id_n_time 

      integer :: ids4(4) ! vectors of dim id's
      ! character(len=120) :: long_name NOT USED
!
! Create new data set (overwrite any pre-existing file):
!
      istat = nf_create(filename,NF_CLOBBER,ncid) 
      if (istat /= NF_NOERR)  &
          call handle_ncerr_gt(istat,'Error from nf_create',1)
      write(6,"(/,'write_hist: Created netcdf graphics file ',a)") trim(filename) 

!
! Define dimensions:
!
      istat = nf_def_dim(ncid,"n_levels",n_levels,id_n_levels)
      istat = nf_def_dim(ncid,"n_lats",n_lats,id_n_lats)
      istat = nf_def_dim(ncid,"n_lons",n_lons,id_n_lons)
      istat = nf_def_dim(ncid,"n_time",NF_UNLIMITED,id_n_time)
!
! Define variables:
!
! variables are defined as 4d with the 4th dimension being the unlimited Time:
      ids4 = (/id_n_levels,id_n_lats,id_n_lons,id_n_time/)

      istat = nf_def_var(ncid,"wind_southwards_ms1",NF_REAL,4,ids4,idv_vsouth2)
      istat = nf_def_var(ncid,"wind_eastwards_ms1",NF_REAL,4,ids4,idv_veast2)
      istat = nf_def_var(ncid,"wvz",NF_REAL,4,ids4,idv_wvz2)
      istat = nf_def_var(ncid,"rmt",NF_REAL,4,ids4,idv_rmt2)
      istat = nf_def_var(ncid,"temperature_K",NF_REAL,4,ids4,idv_tn2)
      istat = nf_def_var(ncid,"ht",NF_REAL,4,ids4,idv_ht2)
      istat = nf_def_var(ncid,"psao",NF_REAL,4,ids4,idv_psao2)
      istat = nf_def_var(ncid,"psmo",NF_REAL,4,ids4,idv_psmo2)
      istat = nf_def_var(ncid,"Ne_m3",NF_REAL,4,ids4,idv_Ne2)
      istat = nf_def_var(ncid,"UT",NF_REAL,1,id_n_time,idv_UT2)
!
! Take out of define mode:
      istat = nf_enddef(ncid)
!
! Close data set:
!
      istat = nf_close(ncid)
      if (istat /= NF_NOERR) &
        call handle_ncerr_gt(istat,'Error from nf_close',0)



end subroutine create_gt_netcdf_graphics







subroutine write_gt_netcdf_graphics(filename, it, &
                                    electron_density_m3, Universal_Time_hours)

      include 'netcdf.inc'

      character(len=*), intent(in) :: filename
      integer :: it
      real(kind=8) :: electron_density_m3(n_levels,n_lats,n_lons)
      REAL(kind=8) :: Universal_Time_hours

!-----------------
! Local variables
!-----------------
      integer :: ncid, istat
      REAL(kind=4) :: Universal_Time_hours_sngl


      integer :: data_start(4)
      integer :: data_count(4)

      real(kind=4) :: wind_southwards_ms1_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: wind_eastwards_ms1_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: wvz_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: rmt_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: temperature_K_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: ht_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: psao_sngl(n_levels,n_lats,n_lons)
      real(kind=4) :: psmo_sngl(n_levels,n_lats,n_lons)

      real(kind=4) :: electron_density_m3_sngl(n_levels,n_lats,n_lons)

! Begin code ========================================================================

!------------------------------
! Open the existing data set:
!------------------------------
      istat = nf_open(filename,NF_WRITE,ncid) 
      if (istat /= NF_NOERR)  &
          call handle_ncerr_gt(istat,'Error from nf_open',1)

!-------------------------------
! Write variables to the file:
!-------------------------------
      data_count(1) = n_levels
      data_count(2) = n_lats
      data_count(3) = n_lons
      data_count(4) = 1
      data_start(1) = 1
      data_start(2) = 1
      data_start(3) = 1
      data_start(4) = it

      wind_southwards_ms1_sngl(:,:,:) = real(wind_southwards_ms1(:,:,:))
      wind_eastwards_ms1_sngl(:,:,:)  = real(wind_eastwards_ms1(:,:,:))
      wvz_sngl(:,:,:)                 = real(wvz(:,:,:))
      rmt_sngl(:,:,:)                 = real(rmt(:,:,:))
      temperature_K_sngl(:,:,:)       = real(temperature_K(:,:,:))
      ht_sngl(:,:,:)                  = real(ht(:,:,:))
      psao_sngl(:,:,:)                = real(psao(:,:,:))
      psmo_sngl(:,:,:)                = real(psmo(:,:,:))
      electron_density_m3_sngl(:,:,:) = real(electron_density_m3(:,:,:))
      Universal_Time_hours_sngl       = real(Universal_Time_hours)
 

      istat = nf_put_vara_real(ncid, idv_vsouth2, data_start,data_count, wind_southwards_ms1_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var vnx', 0)

      istat = nf_put_vara_real(ncid, idv_veast2, data_start, data_count, wind_eastwards_ms1_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var vny', 0)

      istat = nf_put_vara_real(ncid, idv_wvz2, data_start, data_count, wvz_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var wvz', 0)

      istat = nf_put_vara_real(ncid, idv_rmt2, data_start, data_count, rmt_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var rmt', 0)

      istat = nf_put_vara_real(ncid, idv_tn2, data_start, data_count, temperature_K_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var tn', 0)

      istat = nf_put_vara_real(ncid, idv_ht2, data_start, data_count, ht_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var ht', 0)

      istat = nf_put_vara_real(ncid, idv_psao2, data_start, data_count, psao_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var psao', 0)

      istat = nf_put_vara_real(ncid, idv_psmo2, data_start, data_count, psmo_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var psmo', 0)

      istat = nf_put_vara_real(ncid, idv_Ne2, data_start, data_count, electron_density_m3_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var Electron Density', 0)

      istat = nf_put_vara_real(ncid, idv_UT2, it, 1, Universal_Time_hours_sngl)
      if (istat /= NF_NOERR) call handle_ncerr_gt(istat, 'Error writing var UT', 0)

!-----------------
! Close data set:
!-----------------
      istat = nf_close(ncid)
      if (istat /= NF_NOERR) &
        call handle_ncerr_gt(istat, 'Error from nf_close', 0)




end subroutine write_gt_netcdf_graphics





!========================================================================================================
! subroutine to write current values of thermosphere 
! 
! Updated w/ fortran 90 netcdf calls - lrm20130529
!----------------------------------------------------------

subroutine write_gt_netcdf_history(filename)

use netcdf
IMPLICIT NONE

 character(len=*), intent(in) :: filename

!REAL*8, intent(in) :: latitudes(n_lats)
! REAL*8, intent(in) :: longitudes(n_lons)
!real :: lats(NLATS), lons(NLONS) from netcdf example

! Local variables ----------------------------------

integer :: levels(n_levels)

integer :: ncid 
integer :: id_n_levels , id_n_lats , id_n_lons
integer :: lat_varid, lon_varid, levels_varid

integer :: idv_vsouth , idv_veast , idv_wvz
integer :: idv_eps , idv_rmt , idv_tn, idv_ht
integer :: idv_psao , idv_psmo , idv_psmn

integer :: ids3(3) ! vectors of dim id's
integer :: I

 character(len=120) :: long_name

 character (len = *), parameter :: PRESSURE_UNITS = "mb"
 character (len = *), parameter :: LAT_UNITS = "degrees_north"
 character (len = *), parameter :: LON_UNITS = "degrees_east"


! BEGIN CODE -----------------------------------------------

! Always check the return code of every netCDF function call. In
! this example program, wrapping netCDF calls with "call check()"
! makes sure that any return which is not equal to nf90_noerr (0)
! will print a netCDF error message and exit.

! Need to create pressure, lat, lon, pressures
 levels = (/ (I, I = 1, n_levels) /)
! xlhs_n = (/ (I*grd_res_n, I = 0, i_new-1) /)

!-----------------------------------------------------------------------
! Create the netCDF file. The nf90_clobber parameter tells netCDF to
! overwrite this file, if it already exists.
!-----------------------------------------------------------------------
 call check( nf90_create(filename, NF90_CLOBBER, ncid) )

 !write(6, "(/,'write_gt_netcdf_history: Created netcdf file ',a)") trim(filename) 

 !--------------------------------------------------------------
 ! Define the dimensions. NetCDF will hand back an ID for each.
 !--------------------------------------------------------------
 call check( nf90_def_dim(ncid, "n_levels", n_levels, id_n_levels) )
 call check( nf90_def_dim(ncid, "n_lats", n_lats, id_n_lats) )
 call check( nf90_def_dim(ncid, "n_lons", n_lons, id_n_lons) )

 !   write(6, "('Defined dimensions on file ',a)") trim(filename)

 ! Define the coordinate variables. They will hold the coordinate
 ! information, that is, the latitudes and longitudes. A varid is
 ! returned for each.
 call check( nf90_def_var(ncid, "n_levels", NF90_REAL, id_n_levels, levels_varid) )
 call check( nf90_def_var(ncid, "n_lats", NF90_REAL, id_n_lats, lat_varid) )
 call check( nf90_def_var(ncid, "n_lons", NF90_REAL, id_n_lons, lon_varid) )




!------------------------
! Define variables:
!------------------------
! 3d variables:
 ids3 = (/id_n_levels, id_n_lats, id_n_lons/)

 call check( nf90_def_var(ncid, "wind_southwards_ms1", NF90_DOUBLE, ids3, idv_vsouth) )

 long_name = "Southwards Wind"
 call check(nf90_put_att(ncid, idv_vsouth, "long_name", trim(long_name)) )

 call check( nf90_put_att(ncid, idv_vsouth, "units", "ms-1" ) )

 call check( nf90_def_var(ncid, "wind_eastwards_ms1", NF90_DOUBLE, ids3, idv_veast) )

 long_name = "Eastwards Wind"
 call check(nf90_put_att(ncid, idv_veast, "long_name", trim(long_name)) )

 call check( nf90_put_att(ncid, idv_veast, "units", "ms-1" ) )

 call check( nf90_def_var(ncid, "wvz", NF90_DOUBLE, ids3, idv_wvz) )

 call check( nf90_def_var(ncid, "eps", NF90_DOUBLE, ids3, idv_eps) )

 call check( nf90_def_var(ncid, "rmt", NF90_DOUBLE, ids3, idv_rmt) )

 call check( nf90_def_var(ncid, "temperature_K", NF90_DOUBLE, ids3, idv_tn) )

 call check( nf90_def_var(ncid, "ht", NF90_DOUBLE, ids3, idv_ht) )

 call check( nf90_def_var(ncid, "psao", NF90_DOUBLE, ids3, idv_psao) )

 call check( nf90_def_var(ncid, "psmo", NF90_DOUBLE, ids3, idv_psmo) )

 call check( nf90_def_var(ncid, "psmn", NF90_DOUBLE, ids3, idv_psmn) )

 !--------------------------------------------------------------------
 ! End define mode. This tells netCDF we are done defining metadata.
 !--------------------------------------------------------------------
 call check( nf90_enddef(ncid) )

!-------------------------------------
! Write data variables to the file. 
! Although netCDF supports
! reading and writing subsets of data, in this case we write all the
! data in one operation.
!-------------------------------------
! write(6, "('Writing variables to file ',a,'..')") trim(filename)

 ! Write the coordinate variable data. This will put the latitudes
 ! and longitudes of our data grid into the netCDF file.
 call check( nf90_put_var(ncid, levels_varid, levels) )
 call check( nf90_put_var(ncid, lat_varid, therm_model_geo_lat_deg) )
 call check( nf90_put_var(ncid, lon_varid, therm_model_geo_long_deg) )

 call check( nf90_put_var(ncid, idv_vsouth, wind_southwards_ms1) )

 call check( nf90_put_var(ncid, idv_veast, wind_eastwards_ms1) )

 call check( nf90_put_var(ncid, idv_wvz, wvz) )

 call check( nf90_put_var(ncid, idv_eps, eps) )

 call check( nf90_put_var(ncid, idv_rmt, rmt) )

 call check( nf90_put_var(ncid, idv_tn, temperature_K) )

 call check( nf90_put_var(ncid, idv_ht, ht) )

 call check( nf90_put_var(ncid, idv_psao, psao) )

 call check( nf90_put_var(ncid, idv_psmo, psmo) )

 call check( nf90_put_var(ncid, idv_psmn, psmn) )

 !------------------
 ! Close the file.
 !------------------
 call check( nf90_close(ncid) )
       

end subroutine write_gt_netcdf_history

! ===============================================================================================================

subroutine check(status)

 use netcdf
 IMPLICIT NONE

 integer, intent ( in) :: status
     
 if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
 end if
end subroutine check

! ===============================================================================================================
!-----------------------------------------------------------------------



subroutine read_gt_netcdf_history(filename)

      include 'netcdf.inc'

      character(len=*), intent(in) :: filename
      integer :: ncid, istat, id


!
! Open file for reading:
!
      istat = nf_open(filename, NF_NOWRITE, ncid)
      if (istat /= NF_NOERR) &
          call handle_ncerr_gt(istat, 'Error from nf_open', 1)
      write(6, "(/,'read_gt_netcdf_history : Opened file ',a,' for reading.')") &
            trim(filename)

      istat = nf_inq_varid(ncid, "wind_southwards_ms1", id)
      istat = nf_get_var_double(ncid, id, wind_southwards_ms1)

      istat = nf_inq_varid(ncid, "wind_eastwards_ms1", id)
      istat = nf_get_var_double(ncid, id, wind_eastwards_ms1)

      istat = nf_inq_varid(ncid, "wvz", id)
      istat = nf_get_var_double(ncid, id, wvz)

      istat = nf_inq_varid(ncid, "eps", id)
      istat = nf_get_var_double(ncid, id, eps)

      istat = nf_inq_varid(ncid, "rmt", id)
      istat = nf_get_var_double(ncid, id, rmt)

      istat = nf_inq_varid(ncid, "temperature_K", id)
      istat = nf_get_var_double(ncid, id, temperature_K)

      istat = nf_inq_varid(ncid, "ht", id)
      istat = nf_get_var_double(ncid, id, ht)

      istat = nf_inq_varid(ncid, "psao", id)
      istat = nf_get_var_double(ncid, id, psao)

      istat = nf_inq_varid(ncid, "psmo", id)
      istat = nf_get_var_double(ncid, id, psmo)

      istat = nf_inq_varid(ncid, "psmn", id)
      istat = nf_get_var_double(ncid, id, psmn)

      istat = nf_close(ncid)
      write(6, "('Completed read of file ',a)") trim(filename)



end subroutine read_gt_netcdf_history


! ------------------------------------------------------------------------



subroutine read_gtgip_netcdf_history(filename, & !nm20110913
                                  Electron_density_m3, O_plus_density_m3, &
                                  NO_plus_density_m3, O2_plus_density_m3, Ti_Oplus_K )

include 'netcdf.inc'

character(len=*),intent(in) :: filename
integer :: ncid,istat,id

integer :: idv_vsouth , idv_veast , idv_wvz
integer :: idv_eps , idv_rmt , idv_tn, idv_ht
integer :: idv_psao , idv_psmo , idv_psmn

!nm20110913
real(kind=8), DIMENSION(n_levels,n_lats,n_lons), INTENT(IN):: &
     Electron_density_m3, O_plus_density_m3, NO_plus_density_m3, O2_plus_density_m3, Ti_Oplus_K

integer :: idv_ne,idv_oplus,idv_noplus,idv_o2plus,idv_tio
logical, parameter :: sw_add_ion_param=.true.

integer :: ids3(3) ! vectors of dim id's
character(len=120) :: long_name

! begin code ==============================================


!
! Open file for reading:
      istat = nf_open(filename,NF_NOWRITE,ncid)
      if (istat /= NF_NOERR) &
        call handle_ncerr_gt(istat,'Error from nf_open',1)
      write(6,"(/,'read_hist: Opened file ',a,' for reading.')") &
        trim(filename)

      istat = nf_inq_varid(ncid,"wind_southwards_ms1",id)
      istat = nf_get_var_double(ncid,id,wind_southwards_ms1)
!     write(6,"('read_hist: wind_southwards_ms1 min,max=',2e12.4)") minval(wind_southwards_ms1), &
!             maxval(wind_southwards_ms1)

      istat = nf_inq_varid(ncid,"wind_eastwards_ms1",id)
      istat = nf_get_var_double(ncid,id,wind_eastwards_ms1)
!     write(6,"('read_hist: wind_southwards_ms1 min,max=',2e12.4)") minval(wind_eastwards_ms1), &
!             maxval(wind_eastwards_ms1)

      istat = nf_inq_varid(ncid,"wvz",id)
      istat = nf_get_var_double(ncid,id,wvz)
!     write(6,"('read_hist: wvz min,max=',2e12.4)") minval(wvz), maxval(wvz)

      istat = nf_inq_varid(ncid,"eps",id)
      istat = nf_get_var_double(ncid,id,eps)
!     write(6,"('read_hist: eps min,max=',2e12.4)") minval(eps),maxval(eps)

      istat = nf_inq_varid(ncid,"rmt",id)
      istat = nf_get_var_double(ncid,id,rmt)
!     write(6,"('read_hist: rmt min,max=',2e12.4)") minval(rmt),maxval(rmt)

      istat = nf_inq_varid(ncid,"temperature_K",id)
      istat = nf_get_var_double(ncid,id,temperature_K)
!     write(6,"('read_hist: temperature_K min,max=',2e12.4)") minval(temperature_K), &
!             maxval(temperature_K)

      istat = nf_inq_varid(ncid,"ht",id)
      istat = nf_get_var_double(ncid,id,ht)
!     write(6,"('read_hist: ht min,max=',2e12.4)") minval(ht),maxval(ht)

      istat = nf_inq_varid(ncid,"psao",id)
      istat = nf_get_var_double(ncid,id,psao)
!     write(6,"('read_hist: psao min,max=',2e12.4)") minval(psao),maxval(psao)

      istat = nf_inq_varid(ncid,"psmo",id)
      istat = nf_get_var_double(ncid,id,psmo)
!     write(6,"('read_hist: psmo min,max=',2e12.4)") minval(psmo),maxval(psmo)

      istat = nf_inq_varid(ncid,"psmn",id)
      istat = nf_get_var_double(ncid,id,psmn)
!     write(6,"('read_hist: psmn min,max=',2e12.4)") minval(psmn),maxval(psmn)

!commented out until proper input files is ready!
      if ( sw_add_ion_param ) then
!nm20110913
! electron density [m-3]
      istat = nf_inq_varid(ncid,"ne",id)
      istat = nf_get_var_double(ncid,id,Electron_density_m3)
     write(6,"('read_hist: Ne[m-3] min,max=',2e12.4)") minval(Electron_density_m3),maxval(Electron_density_m3)
! O+ density [m-3]
      istat = nf_inq_varid(ncid,"oplus",id)
      istat = nf_get_var_double(ncid,id,O_plus_density_m3)
     write(6,"('read_hist: O+[m-3] min,max=',2e12.4)") minval(O_plus_density_m3),maxval(O_plus_density_m3)
! NO+ density [m-3]
      istat = nf_inq_varid(ncid,"noplus",id)
      istat = nf_get_var_double(ncid,id,NO_plus_density_m3)
     write(6,"('read_hist: NO+[m-3] min,max=',2e12.4)") minval(NO_plus_density_m3),maxval(NO_plus_density_m3)
! O2+ density [m-3]
      istat = nf_inq_varid(ncid,"o2plus",id)
      istat = nf_get_var_double(ncid,id,O2_plus_density_m3)
     write(6,"('read_hist: O2+[m-3] min,max=',2e12.4)") minval(O2_plus_density_m3),maxval(O2_plus_density_m3)
! O+ temperature [K]
      istat = nf_inq_varid(ncid,"tio",id)
      istat = nf_get_var_double(ncid,id,Ti_Oplus_K)
     write(6,"('read_hist: To+[K] min,max=',2e12.4)") minval(Ti_Oplus_K),maxval(Ti_Oplus_K)
     end if !( sw_add_ion_param ) then

!nm20110913
      istat = nf_close(ncid)
      write(6,"('Completed read of file ',a)") trim(filename)



end subroutine read_gtgip_netcdf_history



! -----------------------------------------------------------------------



!-----------------------------------------------------------------------

subroutine handle_ncerr_gt(istat, msg, ifatal)

      include 'netcdf.inc'
!
! Handle a netcdf lib error:
!
      integer, intent(in) :: istat, ifatal
      character(len=*), intent(in) :: msg

      write(6,"(/72('-'))")
      write(6,"('>>> Error from netcdf library:')")
      write(6,"(a)") trim(msg)
      write(6,"('istat=',i5)") istat
      write(6,"(a)") nf_strerror(istat)
      write(6,"(72('-')/)")
      if (ifatal > 0) stop 'handle_ncerr'

return
end subroutine handle_ncerr_gt

!-----------------------------------------------------------------------






subroutine read_in_and_normalise_Hough_modes(hough11, &
                                             hough22, &
                                             hough23, &
                                             hough24, &
                                             hough25)

      IMPLICIT NONE

      INTEGER i , j , k
      REAL*8    hough11(181)
      REAL*8    hough22(181)
      REAL*8    hough23(181)
      REAL*8    hough24(181)
      REAL*8    hough25(181)

         READ(36,99013) (hough11(k), k=1,181)
99013 FORMAT (1x,6E12.4)
         READ(33,99013) (hough22(k), k=1,181)
         READ(33,99013) (hough23(k), k=1,181)
         READ(33,99013) (hough24(k), k=1,181)
         READ(40,99013) (hough25(k), k=1,181)
         CLOSE (33)
         CLOSE (36)
         CLOSE (40)
!
! ** in the following have multiplied hough functions by
! ** factors to normalize their amplitudes to 1.0
! ** IMW, May 1996
!
         i = 1
         j = 1
         DO 1330 k = 1 , 181
            IF ( k == i ) THEN
               hough11(j) = hough11(k)
               hough11(j) = 1./1.650*hough11(j)
               j = j + 1
               i = (j-1)*2 + 1
            ENDIF
 1330    CONTINUE

         i = 1
         j = 1
         DO 1350 k = 1 , 181
            IF ( k == i ) THEN
               hough22(j) = hough22(k)
               hough22(j) = 1./1.165*hough22(j)
               j = j + 1
               i = (j-1)*2 + 1
            ENDIF
 1350    CONTINUE

         i = 1
         j = 1
         DO 1400 k = 1 , 181
            IF ( k == i ) THEN
               hough23(j) = hough23(k)
             hough23(j) = -1./1.093*hough23(j)
               j = j + 1
               i = (j-1)*2 + 1
            ENDIF
 1400    CONTINUE

         i = 1
         j = 1
         DO 1450 k = 1 , 181
            IF ( k == i ) THEN
               hough24(j) = hough24(k)
             hough24(j) = 1./1.078*hough24(j)
               j = j + 1
               i = (j-1)*2 + 1
            ENDIF
 1450    CONTINUE

         i = 1
         j = 1
         DO 1455 k = 1 , 181
            IF ( k == i ) THEN
               hough25(j) = hough25(k)
             hough25(j) = 1./1.075*hough25(j)
               j = j + 1
               i = (j-1)*2 + 1
            ENDIF
 1455    CONTINUE

      return

end subroutine read_in_and_normalise_Hough_modes


!==========================================================================
!==========================================================================


Subroutine calculate_magnetic_parameters_using_apex(inFileName, &
                 B_magnitude_apex_nT, B_dip_angle_apex_degrees, &
                 Magnetic_latitude_degrees, Magnetic_longitude_degrees)

  IMPLICIT NONE

  character(*), intent(in) :: inFileName

  REAL*8 :: B_magnitude_apex_nT(91,20)
  REAL*8 :: B_dip_angle_apex_degrees(91,20)
  REAL*8 :: Magnetic_latitude_degrees(91,20)
  REAL*8 :: Magnetic_longitude_degrees(91,20)

  !------------------
  ! Local variables
  !------------------
  INTEGER :: ilat, ilon

  REAL*8 :: B_east_nT , B_north_nT , B_up_nT
  REAL*8 :: GLON2 , GLAT2 , alon , qdlat

  ! REAL*8  B_declination_apex_degrees(91,20) NOT USED
  REAL*8 :: B_magnitude_horizontal_apex_nT

  !REAL*8 :: PI , DTR
  !PARAMETER (PI=3.14159, DTR=PI/180.0)

  !REAL*8 :: PI
  !PARAMETER (PI=3.14159)

  INTEGER, PARAMETER  :: fileUnit = 117

  !===========================================================


  ! Changed July 19, 2011 LRM 
  !------------------------------------------------------------------------
  ! Open the static file Geographic_to_Magnetic_91_20.2000.0.format
  ! LRM fixed this to use specified static_file_location from the input file
  ! --------------------------------------------------------------------------
  open(fileUnit, file= inFileName, &  ! Geographic_to_Magnetic formatted file
       form='formatted', status='old')


  do ilon = 1 , 20
     do ilat = 1 , 91

            !read(117,4123) GLON2 , GLAT2 , alon , qdlat , B_east_nT , B_north_nT , B_up_nT
 !4123  format(2f6.0,2f9.2,3e14.6)
 !READ(*,'(2I5,F10.2)')

        read(fileUnit,'(2f6.0,2f9.2,3e14.6)') GLON2 , GLAT2 , alon , qdlat , B_east_nT , B_north_nT , B_up_nT

        B_magnitude_apex_nT(ilat,ilon) = sqrt((B_east_nT*B_east_nT) + (B_north_nT*B_north_nT) + (B_up_nT*B_up_nT))

        B_magnitude_horizontal_apex_nT = sqrt((B_east_nT*B_east_nT) + (B_north_nT*B_north_nT))
        B_dip_angle_apex_degrees(ilat,ilon) = ( atan ((0.0 - B_up_nT) / B_magnitude_horizontal_apex_nT) ) / dtr

        Magnetic_latitude_degrees(ilat,ilon) = qdlat
        Magnetic_longitude_degrees(ilat,ilon) = alon

        if (Magnetic_latitude_degrees(ilat,ilon) > 87.0) then
            write(6,*) ' mag lats ', ilon , ilat , Magnetic_latitude_degrees(ilat,ilon)
        endif

     enddo ! ilat = 1 , 91
  enddo ! ilon = 1 , 20

  !print *,'calculate_magnetic_parameters_using_apex : Magnetic_latitude_degrees = ', &
  !         Magnetic_latitude_degrees
  !print *,'calculate_magnetic_parameters_using_apex : '
  !print *,'calculate_magnetic_parameters_using_apex : Magnetic_longitude_degrees = ', &
  !         Magnetic_longitude_degrees


  close(fileUnit)

  return

end Subroutine calculate_magnetic_parameters_using_apex


!=========================================================================
!=========================================================================




SUBROUTINE neutral_composition_equation(XX, PREs, I11, F107, VX, VY, HT, ARMt, OM, T, PSAo, PSMo, &
                       PSMn)
      IMPLICIT NONE

      REAL*8 &
           a1 , a11e , a11n , a11o , a11s , a11w , a12e , a12n , a12o , &
           a12s , a12w , a2 , a21e , a21n , a21o , a21s , a21w , a22e , &
           a22n , a22o
      REAL*8 &
           a22s , a22w , a3 , a4 , a5 , a6 , amu , ARMt , b1 , b2 , b3 , &
           b4 , b5 , b6 , beta , bot , bz , cotth
      REAL*8 &
           d1 , d12 , d13 , d2 , d23 , d3 , d4 , da11 , da11x , da11y , &
           da12 , da12x , da12y , da21 , da21x , da21y , da22 , da22x
      REAL*8 da22y , dbx , dby , dbz , ddhtx , ddhty , ddmw , ddmwx , &
             ddmwy , ddp1 , ddp1x , ddp1y , ddp2 , ddp2x , ddp2y , &
             delphi , delta , deltha , den , den2
      REAL*8 &
           deno , dhtx , dhty , dknt , dkntx , dknty , dl1 , dl2 , dn , &
           dn2 , domo , dp1 , dp1x , dp1y , dp2 , dp2x , dp2y , dphi , &
           dr
      REAL*8 &
           dr1 , dr2 , dsht , dshtx , dshty , dt , dth , dvxx , dvyy , &
           dz , f1 , F107 , f2 , f3 , f4 , fac , fp1
      REAL*8 g2 , gamma , grav , gscon , ho2 , HT , hte , htn , hto , &
             hts , htw , oh , OM , omo , omo2 , p1 , p1e , p1n , p1s , p1t
      !REAL*8 p1w , p2 , p2e , p2n , p2s , p2t , p2w , p3 , pi , pr2 , &
      REAL*8 p1w , p2 , p2e , p2n , p2s , p2t , p2w , p3 , pr2 , &
             PREs , ps , PSAo , psaop , PSMn , psmnp , PSMo , psmop , q1
      REAL*8 q2 , r1x , r1y , r1z , r2x , r2y , r2z , R0 , rlat , rmte , &
             rmtn , rmtp , rmts , rmtw , ro , ro2 , s1 , s2 , sht
      REAL*8 &
           shte , shtn , shto , shts , shtw , sign , sth , sum1 , sum2 , &
           sx1 , sx10 , sx15 , sx2 , sx3 , sx4 , sx9
      REAL*8 T , th , theta , tot , ts , VX , vxn , vxo , vxo2 , vxs , &
             VY , vye , vyo , vyo2 , vyw
      REAL*8 &
           x1 , x10 , x2 , x3 , x4 , x5 , x6 , x7 , x8 , x9 , y1 , y10 , &
           y2 , y3 , y4 , y5 , y6 , y7 , y8 , y9
      INTEGER I11 , I3 , ifp , K3 , l , L3 , le , LM , lw , m , MM , &
              mmm , mn , mp , ms , n , nb , nd , NM
      INTEGER nu
!*** End of declarations inserted by SPAG

      PARAMETER (I3=15, K3=91, L3=20, NM=15, MM=91, LM=20)

      DIMENSION VX(NM,MM,LM) , VY(NM,MM,LM) , HT(NM,MM,LM) , &
                OM(NM,MM,LM) , PSMn(NM,MM,LM) , ARMt(NM,MM,LM) , &
                PSAo(NM,MM,LM) , PSMo(NM,MM,LM) , T(I3,K3,L3) , sht(I3) &
                , den(I3) , p1t(I3,K3,L3) , p2t(I3,K3,L3) , oh(I3) , &
                ho2(I3) , dp1(I3) , dp2(I3) , dknt(I3) , da11(I3) , &
                da12(I3) , da21(I3) , da22(I3) , dsht(I3) , domo(I3) , &
                ddmw(I3) , ddp1(I3) , ddp2(I3) , dp1x(I3) , dp2x(I3) , &
                dkntx(I3) , da11x(I3) , da12x(I3) , da21x(I3) , &
                da22x(I3) , dshtx(I3) , dvxx(I3) , ddp1x(I3) , ddp2x(I3) &
                , ddmwx(I3) , dhtx(I3) , ddhtx(I3) , dp1y(I3) , dp2y(I3) &
                , dknty(I3) , da11y(I3) , da12y(I3) , da21y(I3) , &
                da22y(I3) , dshty(I3) , ddp1y(I3) , ddp2y(I3) , &
                ddmwy(I3) , PREs(I3) , dhty(I3) , ddhty(I3) , dvyy(I3) , &
                vyo(I3)

      DIMENSION p1(I3) , p2(I3) , p3(I3) , s1(I3) , s2(I3) , a3(I3) , &
                a4(I3) , a5(I3) , b3(I3) , b4(I3) , b5(I3) , dr1(I3) , &
                dr2(I3) , a11o(I3) , a12o(I3) , a21o(I3) , a22o(I3) , &
                rmtn(I3) , p1n(I3) , p2n(I3) , a11n(I3) , a12n(I3) , &
                a21n(I3) , a22n(I3) , rmts(I3) , p1s(I3) , p2s(I3) , &
                a11s(I3) , a12s(I3) , a21s(I3) , a22s(I3) , shto(I3) , &
                shtn(I3) , shts(I3) , hto(I3) , htn(I3) , hts(I3) , &
                deno(I3) , omo(I3) , vxo(I3) , vxs(I3) , vxn(I3) , &
                rmte(I3) , p1e(I3) , a11e(I3) , a12e(I3) , a21e(I3) , &
                a22e(I3) , rmtw(I3) , p1w(I3) , p2w(I3) , a11w(I3) , &
                a12w(I3) , a21w(I3) , a22w(I3) , shte(I3) , shtw(I3) , &
                hte(I3) , htw(I3) , p2e(I3) , vye(I3) , vyw(I3) , &
                rmtp(15,2) , psaop(15,2) , psmop(15,2) , psmnp(15,2)
      INTEGER XX , zz
      REAL*8 nt(I3) , m1 , m2 , m3 , jp(I3) , k(I3) , l1 , l2
      REAL*8 alpha , rmto(I3) , dmw(I3) , dmwx(I3) , dmwy(I3)
      PARAMETER (R0=6.370E06)

      dth = 2.0
      dphi = 18.0
      dt = 60.
      dz = 1.0
      CALL OHHO2(oh,ho2)
      gscon = 8.3141E+03
      amu = 1.66E-27
      grav = 9.5
      bz = 1.38E-23
      !pi = 3.14159
      dr = pi/180.
      m1 = 16.*amu
      m2 = 32.*amu
      m3 = 28.*amu
      ts = 273.
      ps = 1.01325E+05

      CALL JPROFIL(jp,F107)
!C  **
!C  NOTE :- IN THIS SUBROUTINE ARMT ARE THE VALUES OF MOLECULAR WEIGHT
!C                 IN AMU'S, AND RMT ARE THE VALUES IN KGS.
!c  **
!c  start longitude scan
!c  **
      DO 300 l = 1 , 20
         le = l + 1
         lw = l - 1
         IF ( le == 21 ) le = 1
         IF ( lw == 0 ) lw = 20
!C  **
!C  START LATITUDE SCAN
!C  **
         DO 250 m = 3 , 89
!C  SET CO-LATITUDE(DEGREES)
            theta = 180. - (m-1)*dth
            sign = 1.0
!C  SET DELphi TO TWICE RESOLUTION (DEGREES)
            delphi = dphi*dr*2.
!C  SET DELTHA TO TWICE RESOLUTION (DEGREES)
            deltha = dth*dr*2.
!C  SET LATITUDE (DEGREES)
            rlat = 90.0 - theta
!C  SET CO-LATITUDE (RADIANS)
            th = theta*dr
!C  SET LATITUDE (RADIANS)
            rlat = rlat*dr
            dn = dz*2.
            DO 20 n = 1 , NM
               nt(n) = PREs(n)/bz/T(n,m,l)
               den(n) = nt(n)*ARMt(n,m,l)*amu
               sht(n) = gscon*T(n,m,l)/ARMt(n,m,l)/grav
 20         CONTINUE
            CALL EDDY(NM,k,HT(1,m,l))
!C  **
!C  SET NB TO LEVEL ABOVE WHICH ATMOSPHERE IN DIFFUSIVE EQUILIBRIUM
!C  **
            nb = 13
!C  **
!C  EVALUATE SOURCES AND SINKS
!C  **
            DO 40 n = 1 , nb
               alpha = 2.76D-46*EXP(710./T(n,m,l))
               gamma = 1.E-26
               beta = 4.2E-17
               delta = 3.5E-17
               p1(n) = PSAo(n,m,l)
               p2(n) = PSMo(n,m,l)
               tot = nt(n)
               p3(n) = PSMn(n,m,l)
               ro = p1(n)*den(n)/m1
               ro2 = p2(n)*den(n)/m2
               q1 = 2.*ro2*jp(n)
!C      Q1=2.*RO2*JP(N)*PI*COSANG/2.
!C      IF (COSANG <= 0.0) Q1=0.0
               l1 = 2.*alpha*ro*tot + 2.*gamma*ro + oh(n)*beta + ho2(n) &
                    *delta
               q2 = alpha*ro**2*tot + gamma*ro**2 + &
                    ro*(oh(n)*beta+ho2(n)*delta)*0.5
               l2 = jp(n)
               s1(n) = m1*(q1-l1*ro)/den(n)
               s2(n) = m2*(q2-l2*ro2)/den(n)
 40         CONTINUE
!C  **
!C  EVALUATE A FEW COMMON DENOMINATORS
!C  **
            sth = SIN(th)
            cotth = 1.0/TAN(th)
            sx1 = deltha*R0
            sx3 = (sx1**2)/4.
            sx9 = R0*sth*delphi
            sx10 = (sx9**2)/4.
            ms = m - 1
            mn = m + 1
!C  **
!C  SET UP ONE-DIMENSIONAL ARRAYS
!C  **
            DO 60 n = 1 , NM
!C  AT CURRENT LATITUDE
!C  **
               p1(n) = PSAo(n,m,l)
               p2(n) = PSMo(n,m,l)
!C  **
               rmto(n) = ARMt(n,m,l)*amu
               tot = nt(n)
               fac = (T(n,m,l)/ts)**1.75*(ps/PREs(n))
               d12 = fac*0.26E-04
               d13 = fac*0.26E-04
               d23 = fac*0.181E-04
               f1 = 1./(m3*d13)
               f2 = 1./(m2*d12)
               f3 = 1./(m1*d12)
               f4 = 1./(m3*d23)
               a11o(n) = -(f1+p2(n)*(f2-f1))/tot
               a12o(n) = (f2-f1)*p1(n)/tot
               a21o(n) = (f3-f4)*p2(n)/tot
               a22o(n) = -(f4+p1(n)*(f3-f4))/tot
               shto(n) = sht(n)
               hto(n) = HT(n,m,l)
               deno(n) = den(n)
               omo(n) = OM(n,m,l)
               vxo(n) = VX(n,m,l)
               vyo(n) = VY(n,m,l)
!C  SOUTH POINTS
               p1s(n) = PSAo(n,ms,l)
               p2s(n) = PSMo(n,ms,l)
               rmts(n) = ARMt(n,ms,l)*amu
               tot = PREs(n)/bz/T(n,ms,l)
               fac = (T(n,ms,l)/ts)**1.75*(ps/PREs(n))
               d12 = fac*0.26E-04
               d13 = fac*0.26E-04
               d23 = fac*0.181E-04
               f1 = 1./(m3*d13)
               f2 = 1./(m2*d12)
               f3 = 1./(m1*d12)
               f4 = 1./(m3*d23)
               a11s(n) = -(f1+p2(n)*(f2-f1))/tot
               a12s(n) = (f2-f1)*p1(n)/tot
               a21s(n) = (f3-f4)*p2(n)/tot
               a22s(n) = -(f4+p1(n)*(f3-f4))/tot
               shts(n) = gscon*T(n,ms,l)/ARMt(n,ms,l)/grav
               hts(n) = HT(n,ms,l)
               vxs(n) = VX(n,ms,l)
!C  NORTH POINTS
               p1n(n) = PSAo(n,mn,l)
               p2n(n) = PSMo(n,mn,l)
               rmtn(n) = ARMt(n,mn,l)*amu
               tot = PREs(n)/bz/T(n,mn,l)
               fac = (T(n,mn,l)/ts)**1.75*(ps/PREs(n))
               d12 = fac*0.26E-04
               d13 = fac*0.26E-04
               d23 = fac*0.181E-04
               f1 = 1./(m3*d13)
               f2 = 1./(m2*d12)
               f3 = 1./(m1*d12)
               f4 = 1./(m3*d23)
               a11n(n) = -(f1+p2(n)*(f2-f1))/tot
               a12n(n) = (f2-f1)*p1(n)/tot
               a21n(n) = (f3-f4)*p2(n)/tot
               a22n(n) = -(f4+p1(n)*(f3-f4))/tot
               shtn(n) = gscon*T(n,mn,l)/ARMt(n,mn,l)/grav
               htn(n) = HT(n,mn,l)
               vxn(n) = VX(n,mn,l)


!C  NORTH POINTS

               p1e(n) = PSAo(n,m,le)
               p2e(n) = PSMo(n,m,le)
               rmte(n) = ARMt(n,m,le)*amu
               tot = PREs(n)/bz/T(n,m,le)

               ! Make sure not dividing by 0.
               !if (ts == 0. .OR. pres(n) == 0.) then

                  ! lrm 09162011
               !   print *,'neutral_composition : n,m,le, T(n,m,le) = ',n,m,le, T(n,m,le)
               !   print *,'neutral_composition : ts, ps, pres(n) = ',ts, ps, pres(n)

               !endif

               ! lrm20120903  Separated out calculations for debugging
               !fac = (T(n,m,le)/ts)**1.75*(ps/PREs(n))
               fac = (T(n,m,le)/ts)
               !print *,'n, m, le(n,m,le), fac = ', n, m, le, T(n,m,le), fac
               fac = fac**1.75   ! crashed for 3 day run
               fac = fac*(ps/PREs(n))


               d12 = fac*0.26E-04
               d13 = fac*0.26E-04
               d23 = fac*0.181E-04
               f1 = 1./(m3*d13)
               f2 = 1./(m2*d12)
               f3 = 1./(m1*d12)
               f4 = 1./(m3*d23)
               a11e(n) = -(f1+p2(n)*(f2-f1))/tot
               a12e(n) = (f2-f1)*p1(n)/tot
               a21e(n) = (f3-f4)*p2(n)/tot
               a22e(n) = -(f4+p1(n)*(f3-f4))/tot
               shte(n) = gscon*T(n,m,le)/ARMt(n,m,le)/grav
               hte(n) = HT(n,m,le)
               vye(n) = VY(n,m,le)
!C  NORTH POINTS
               p1w(n) = PSAo(n,m,lw)
               p2w(n) = PSMo(n,m,lw)
               rmtw(n) = ARMt(n,m,lw)*amu
               tot = PREs(n)/bz/T(n,m,lw)
               fac = (T(n,m,lw)/ts)**1.75*(ps/PREs(n))
               d12 = fac*0.26E-04
               d13 = fac*0.26E-04
               d23 = fac*0.181E-04
               f1 = 1./(m3*d13)
               f2 = 1./(m2*d12)
               f3 = 1./(m1*d12)
               f4 = 1./(m3*d23)
               a11w(n) = -(f1+p2(n)*(f2-f1))/tot
               a12w(n) = (f2-f1)*p1(n)/tot
               a21w(n) = (f3-f4)*p2(n)/tot
               a22w(n) = -(f4+p1(n)*(f3-f4))/tot
               shtw(n) = gscon*T(n,m,lw)/ARMt(n,m,lw)/grav
               htw(n) = HT(n,m,lw)
               vyw(n) = VY(n,m,lw)
 60         CONTINUE
!C  **
!C  SET SOME ARRAYS OF GRADIENTS
!C  **
            DO 80 n = 2 , nb
               dn2 = (dn/2.)**2
               pr2 = PREs(n)**2
               sx2 = -dn*PREs(n)
               sx4 = dn2*pr2
               sx15 = PREs(n)/bz
               nu = n + 1
               nd = n - 1
               dmw(n) = (rmto(nu)-rmto(nd))/sx2
               dp1(n) = (p1(nu)-p1(nd))/sx2
               dp2(n) = (p2(nu)-p2(nd))/sx2
               dknt(n) = (k(nu)*nt(nu)/shto(nu)-k(nd)*nt(nd)/shto(nd)) &
                         /sx2
               da11(n) = (a11o(nu)-a11o(nd))/sx2
               da12(n) = (a12o(nu)-a12o(nd))/sx2
               da21(n) = (a21o(nu)-a21o(nd))/sx2
               da22(n) = (a22o(nu)-a22o(nd))/sx2
               dsht(n) = (shto(nu)-shto(nd))/sx2
               domo(n) = (omo(nu)-omo(nd))/sx2
!C
               ddmw(n) = (rmto(nu)-2.*rmto(n)+rmto(nd))/sx4
               ddp1(n) = (p1(nu)-2.*p1(n)+p1(nd))/sx4
               ddp2(n) = (p2(nu)-2.*p2(n)+p2(nd))/sx4
!C
               dmwx(n) = (rmts(n)-rmtn(n))/sx1
               dp1x(n) = (p1s(n)-p1n(n))/sx1

!              if(m == 42.and.l == 14) then
!              write(6,243) m,n,p1s(n),p1n(n),p1s(n)-p1n(n)
!243           format(2i4,7e12.4)
!              endif
!              if(m == 43.and.l == 14) then
!              write(6,243) m,n,p1s(n),p1n(n),p1s(n)-p1n(n)
!              endif

               dp2x(n) = (p2s(n)-p2n(n))/sx1
               dkntx(n) = k(n)*(sx15/T(n,ms,l)-sx15/T(n,mn,l))/sx1
               da11x(n) = (a11s(n)-a11n(n))/sx1
               da12x(n) = (a12s(n)-a12n(n))/sx1
               da21x(n) = (a21s(n)-a21n(n))/sx1
               da22x(n) = (a22s(n)-a22n(n))/sx1
               dshtx(n) = (shts(n)-shtn(n))/sx1
               dhtx(n) = (hts(n)-htn(n))/sx1
               dvxx(n) = (vxs(n)-vxn(n))/sx1

               ddp1x(n) = (p1s(n)-2.*p1(n)+p1n(n))/sx3
               ddp2x(n) = (p2s(n)-2.*p2(n)+p2n(n))/sx3
               ddmwx(n) = (rmts(n)-2.*rmto(n)+rmtn(n))/sx3
               ddhtx(n) = (hts(n)-2.*hto(n)+htn(n))/sx3

               dmwy(n) = (rmte(n)-rmtw(n))/sx9
               dp1y(n) = (p1e(n)-p1w(n))/sx9
               dp2y(n) = (p2e(n)-p2w(n))/sx9
               dknty(n) = k(n)*(sx15/T(n,m,le)-sx15/T(n,m,lw))/sx9
               da11y(n) = (a11e(n)-a11w(n))/sx9
               da12y(n) = (a12e(n)-a12w(n))/sx9
               da21y(n) = (a21e(n)-a21w(n))/sx9
               da22y(n) = (a22e(n)-a22w(n))/sx9
               dshty(n) = (shte(n)-shtw(n))/sx9
               dhty(n) = (hte(n)-htw(n))/sx9
               dvyy(n) = (vye(n)-vyw(n))/sx9

               ddp1y(n) = (p1e(n)-2.*p1(n)+p1w(n))/sx10
               ddp2y(n) = (p2e(n)-2.*p2(n)+p2w(n))/sx10
               ddmwy(n) = (rmte(n)-2.*rmto(n)+rmtw(n))/sx10
               ddhty(n) = (hte(n)-2.*hto(n)+htw(n))/sx10
 80         CONTINUE


!C  **
!C  EVALUATE TURBULENT DIFFUSION TERMS
!C  **
            DO 100 n = 2 , nb
               den2 = deno(n)**2
               g2 = grav**2
               x1 = ddp1x(n) + cotth*dp1x(n)*sign/R0 + den2*g2*ddp1(n) &
                    + ddp1y(n)
               x2 = dmwx(n)*dp1x(n) + den2*g2*dmw(n)*dp1(n) + dmwy(n) &
                    *dp1y(n)
               x3 = ddmwx(n) + cotth*dmwx(n)*sign/R0 + den2*g2*ddmw(n) &
                    + ddmwy(n)
               x4 = dp1x(n)*dkntx(n) + den2*g2*dp1(n)*dknt(n) + dp1y(n) &
                    *dknty(n)
               x5 = dmwx(n)*dkntx(n) + den2*g2*dmw(n)*dknt(n) + dmwy(n) &
                    *dknty(n)
               a5(n) = k(n)/rmto(n)*(rmto(n)*x1+2.*x2+p1(n)*x3) &
                       + (rmto(n)*x4+p1(n)*x5)*shto(n)/deno(n)
               y1 = ddp2x(n) + cotth*dp2x(n)*sign/R0 + den2*g2*ddp2(n) &
                    + ddp2y(n)
               y2 = dmwx(n)*dp2x(n) + den2*g2*dmw(n)*dp2(n) + dmwy(n) &
                    *dp2y(n)
               y3 = x3
               y4 = dp2x(n)*dkntx(n) + den2*g2*dp2(n)*dknt(n) + dp2y(n) &
                    *dknty(n)
               y5 = x5
               b5(n) = k(n)/rmto(n)*(rmto(n)*y1+2.*y2+p2(n)*y3) &
                       + (rmto(n)*y4+p2(n)*y5)*shto(n)/deno(n)
 100        CONTINUE


!C  **
!C  EVALUATE SMOOTHING TERMS
!C  **
            DO 120 n = 2 , nb
               vxo2 = vxo(n)**2
               vyo2 = vyo(n)**2
               omo2 = omo(n)**2
               x1 = vxo(n)*dvxx(n)*dp1x(n)
               x2 = vxo2*ddp1x(n)
               x3 = omo(n)*domo(n)*dp1(n)
               x4 = omo2*(ddp1(n)-dp1(n)/PREs(n))
               x5 = vyo(n)*dvyy(n)*dp1y(n)
               x6 = vyo2*ddp1y(n)
               a4(n) = (x1+x2+x3+x4+x5+x6)*dt*I11/2.
               y1 = vxo(n)*dvxx(n)*dp2x(n)
               y2 = vxo2*ddp2x(n)
               y3 = omo(n)*domo(n)*dp2(n)
               y4 = omo2*(ddp2(n)-dp2(n)/PREs(n))
               y5 = vyo(n)*dvyy(n)*dp2y(n)
               y6 = vyo2*ddp2y(n)
               b4(n) = (y1+y2+y3+y4+y5+y6)*dt*I11/2.

!     if((l == 15).and.(m == 62))then
!     write(6,198)n,x1,x2,x3,x4,x5,x6
!     write(6,198)n,y1,y2,y3,y4,y5,y6
!198  format(1x,i3,1p6e12.2)
!     endif

 120        CONTINUE


!C  **
!C  EVALUATE MOLECULAR DIFFUSION TERMS
!C  **
            DO 140 n = 2 , nb
               den2 = deno(n)**2
               g2 = grav**2
               x1 = ddp1x(n) + cotth*dp1x(n)*sign/R0 + den2*g2*ddp1(n) &
                    + ddp1y(n)
               x2 = (ddmwx(n)+cotth*dmwx(n)*sign/R0+den2*g2*ddmw(n) &
                    +ddmwy(n))/rmto(n)
               x3 = (dmwx(n)*dp1x(n)+den2*g2*dmw(n)*dp1(n)+dmwy(n) &
                    *dp1y(n))/rmto(n)
               x4 = (dmwx(n)**2+dmwy(n)**2-den2*g2*dmw(n)**2)/rmto(n)**2
               x5 = (ddhtx(n)+ddhty(n)+cotth*dhtx(n)*sign/R0)/shto(n) &
                    *deno(n)*grav
               x6 = dhtx(n) &
                    *(p1(n)*m1*dmwx(n)/shto(n)/rmto(n)**2+(1.-m1/rmto(n) &
                    )*(dp1x(n)/shto(n)-p1(n)*dshtx(n)/shto(n)**2)) &
                    *deno(n)*grav
               x10 = dhty(n) &
                     *(p1(n)*m1*dmwy(n)/shto(n)/rmto(n)**2+(1.-m1/rmto &
                     (n))*(dp1y(n)/shto(n)-p1(n)*dshty(n)/shto(n)**2)) &
                     *deno(n)*grav
               x7 = deno(n)*grav*dmw(n)/shto(n)/rmto(n)**2
               x8 = deno(n)*grav*(1.-m1/rmto(n)) &
                    *(dp1(n)/shto(n)-p1(n)*dsht(n)/shto(n)**2)
               x9 = -den2*g2*dsht(n)*(dp1(n)+p1(n)*dmw(n)/rmto(n)) &
                    /shto(n)
               dr1(n) = x1 + x3 + x6 + x8 + x9 + x10 + p1(n) &
                        *(x2+x4+x5*(1.-m1/rmto(n))+x7*m1)
               y1 = ddp2x(n) + cotth*dp2x(n)*sign/R0 + den2*g2*ddp2(n) &
                    + ddp2y(n)
               y2 = x2
               y3 = (dmwx(n)*dp2x(n)+den2*g2*dmw(n)*dp2(n)+dmwy(n) &
                    *dp2y(n))/rmto(n)
               y4 = x4
               y5 = x5
               y6 = dhtx(n) &
                    *(p2(n)*m2*dmwx(n)/shto(n)/rmto(n)**2+(1.-m2/rmto(n) &
                    )*(dp2x(n)/shto(n)-p2(n)*dshtx(n)/shto(n)**2)) &
                    *deno(n)*grav
               y10 = dhty(n) &
                     *(p2(n)*m2*dmwy(n)/shto(n)/rmto(n)**2+(1.-m2/rmto &
                     (n))*(dp2y(n)/shto(n)-p2(n)*dshty(n)/shto(n)**2)) &
                     *deno(n)*grav
               y7 = x7
               y8 = deno(n)*grav*(1.-m2/rmto(n)) &
                    *(dp2(n)/shto(n)-p2(n)*dsht(n)/shto(n)**2)
               y9 = -den2*g2*dsht(n)*(dp2(n)+p2(n)*dmw(n)/rmto(n)) &
                    /shto(n)
               dr2(n) = y1 + y3 + y6 + y8 + y9 + y10 + p2(n) &
                        *(y2+y4+y5*(1.-m2/rmto(n))+y7*m2)
 140        CONTINUE


            DO 160 n = 2 , nb
               r1x = dp1x(n) + p1(n)*dmwx(n)/rmto(n) + p1(n) &
                     *(1.-m1/rmto(n))*dhtx(n)/shto(n)
               r1y = dp1y(n) + p1(n)*dmwy(n)/rmto(n) + p1(n) &
                     *(1.-m1/rmto(n))*dhty(n)/shto(n)
               r1z = -deno(n)*grav*(dp1(n)+p1(n)*dmw(n)/rmto(n)) - p1(n) &
                     *(1.-m1/rmto(n))/shto(n)
               r2x = dp2x(n) + p2(n)*dmwx(n)/rmto(n) + p2(n) &
                     *(1.-m2/rmto(n))*dhtx(n)/shto(n)
               r2y = dp2y(n) + p2(n)*dmwy(n)/rmto(n) + p2(n) &
                     *(1.-m2/rmto(n))*dhty(n)/shto(n)
               r2z = -deno(n)*grav*(dp2(n)+p2(n)*dmw(n)/rmto(n)) - p2(n) &
                     *(1.-m2/rmto(n))/shto(n)
               bot = a12o(n)*a21o(n) - a22o(n)*a11o(n)
               dbx = a12o(n)*da21x(n) + a21o(n)*da12x(n) - a22o(n) &
                     *da11x(n) - a11o(n)*da22x(n)
               dby = a12o(n)*da21y(n) + a21o(n)*da12y(n) - a22o(n) &
                     *da11y(n) - a11o(n)*da22y(n)
               dbz = -deno(n) &
                     *grav*(a12o(n)*da21(n)+a21o(n)*da12(n)-a22o(n) &
                     *da11(n)-a11o(n)*da22(n))
               dl2 = (a21o(n)*dr1(n)+r1x*da21x(n)-r1z*deno(n) &
                     *grav*da21(n)-a11o(n)*dr2(n)-r2x*da11x(n) &
                     +r2z*deno(n)*grav*da11(n)+r1y*da21y(n)-r2y*da11y(n) &
                     )/bot - (a21o(n)*(r1x*dbx+dbz*r1z+r1y*dby)-a11o(n) &
                     *(r2x*dbx+dbz*r2z+r2y*dby))/bot**2
               dl1 = -(a22o(n)*dr1(n)+r1x*da22x(n)-r1z*deno(n) &
                     *grav*da22(n)-a12o(n)*dr2(n)-r2x*da12x(n) &
                     +r2z*deno(n)*grav*da12(n)+r1y*da22y(n)-r2y*da12y(n) &
                     )/bot + (a22o(n)*(r1x*dbx+dbz*r1z+r1y*dby)-a12o(n) &
                     *(r2x*dbx+dbz*r2z+r2y*dby))/bot**2
               a3(n) = -dl1/deno(n)
               b3(n) = -dl2/deno(n)
 160        CONTINUE


            DO 180 n = 2 , nb
               a1 = -vxo(n)*dp1x(n)
               a6 = -vyo(n)*dp1y(n)
               a2 = -omo(n)*dp1(n)
               b1 = -vxo(n)*dp2x(n)
               b6 = -vyo(n)*dp2y(n)
               b2 = -omo(n)*dp2(n)
               sum1 = a1 + a2 + a3(n) + a4(n) + a5(n) + a6 + s1(n)
               sum2 = b1 + b2 + b3(n) + b4(n) + b5(n) + b6 + s2(n)
               p1t(n,m,l) = sum1*dt*I11 + p1(n)
               p2t(n,m,l) = sum2*dt*I11 + p2(n)
 180        CONTINUE


!C  **
!C  EVALUATE VALUES ABOVE NB ASSUMIMG DIFFUSIVE EQUILIBRIUM
!C  **
            DO 200 n = nb , NM - 1
               nu = n + 1
               nd = n - 1

               dp1(n) = (p1(nu)-p1(nd))/dn
               dp2(n) = (p2(nu)-p2(nd))/dn
               dmw(n) = (rmto(nu)-rmto(nd))/dn
               ddmw(n) = (rmto(nu)-2.*rmto(n)+rmto(nd))/dn2
               d1 = -p1(n)*((dmw(n)+m1)/rmto(n)-1.)
               d2 = p1(n)*(-rmto(n)*ddmw(n)+dmw(n)**2+m1*dmw(n))/rmto(n) &
                    **2 + dp1(n)*(rmto(n)-m1-dmw(n))/rmto(n)
               d3 = -p2(n)*((dmw(n)+m2)/rmto(n)-1.)
               d4 = p2(n)*(-rmto(n)*ddmw(n)+dmw(n)**2+m2*dmw(n))/rmto(n) &
                    **2 + dp2(n)*(rmto(n)-m2-dmw(n))/rmto(n)

!C    X1=-(RMTO(NU)-RMTO(ND))/(DN*RMTO(N))
!C    X2=P1T(N,M)*(1.-M1/RMTO(N))
!C    Y1=X1
!C    Y2=P2T(N,M)*(1.-M2/RMTO(N))
!C    SUM1=P1T(N,M)*X1+X2
!C    SUM2=P2T(N,M)*Y1+Y2
!C    P1T(NU,M)=P1T(N,M)+SUM1*DZ(J)
!C    P2T(NU,M)=P2T(N,M)+SUM2*DZ(J)

               p1t(nu,m,l) = p1t(n,m,l) + d1*dz + d2*dz**2/2.
               p2t(nu,m,l) = p2t(n,m,l) + d3*dz + d4*dz**2/2.
               IF ( p1t(nu,m,l) > .9999 ) p1t(nu,m,l) = .9999
               IF ( p1t(nu,m,l) > .9999 ) p2t(nu,m,l) = 0.000099
 200        CONTINUE

            x1 = (oh(1)*beta+ho2(1)*delta)*m2/(jp(1)*2.*m1)
            x2 = 1. - p3(1)
            p1t(1,m,l) = x2/(1.+x1)
            p2t(1,m,l) = 1. - p1t(1,m,l) - p3(1)

!nm020708: modified to avoid the atomic oxygen getting negative by using the assumption of chemical equilibrium
            if (p1t(2,m,l) <= 0.0) then
                x1 = (oh(2)*beta+ho2(2)*delta)*m2/(jp(2)*2.*m1)
                x2 = 1. - p3(2)
                p1t(2,m,l) = x2/(1.+x1)
                p2t(2,m,l) = 1. - p1t(2,m,l) - p3(2)
            endif

            if (p1t(3,m,l) <= 0.0) then
                x1 = (oh(3)*beta+ho2(3)*delta)*m2/(jp(3)*2.*m1)
                x2 = 1. - p3(3)
                p1t(3,m,l) = x2/(1.+x1)
                p2t(3,m,l) = 1. - p1t(3,m,l) - p3(3)
            endif


            if (p1t(4,m,l) <= 0.0) then
                x1 = (oh(4)*beta+ho2(4)*delta)*m2/(jp(4)*2.*m1)
                x2 = 1. - p3(4)
                p1t(4,m,l) = x2/(1.+x1)
                p2t(4,m,l) = 1. - p1t(4,m,l) - p3(4)
            endif

            if (p1t(5,m,l) <= 0.0) then
                x1 = (oh(5)*beta+ho2(5)*delta)*m2/(jp(5)*2.*m1)
                x2 = 1. - p3(5)
                p1t(5,m,l) = x2/(1.+x1)
                p2t(5,m,l) = 1. - p1t(5,m,l) - p3(5)
            endif

            if (p1t(6,m,l) <= 0.0) then
                x1 = (oh(6)*beta+ho2(6)*delta)*m2/(jp(6)*2.*m1)
                x2 = 1. - p3(6)
                p1t(6,m,l) = x2/(1.+x1)
                p2t(6,m,l) = 1. - p1t(6,m,l) - p3(6)
            endif

            if (p1t(7,m,l) <= 0.0) then
                x1 = (oh(7)*beta+ho2(7)*delta)*m2/(jp(7)*2.*m1)
                x2 = 1. - p3(7)
                p1t(7,m,l) = x2/(1.+x1)
                p2t(7,m,l) = 1. - p1t(7,m,l) - p3(7)
            endif


!C  **
!C  END  LATITUDE SCAN
!C  **
 250     CONTINUE

!C  **
!C  END  longitude SCAN
!C  **
 300  CONTINUE


!C  **
!C  WRITE VALUES ACROSS BEFORE RETURN
!C  **
      DO 400 l = 1 , 20
         DO 350 m = 3 , 89
            DO 320 n = 1 , NM
               PSAo(n,m,l) = p1t(n,m,l)
               PSMo(n,m,l) = p2t(n,m,l)
               PSMn(n,m,l) = 1. - p1t(n,m,l) - p2t(n,m,l)
               IF ( PSMn(n,m,l) < 1.E-07 ) PSMn(n,m,l) = 1.E-07
               ARMt(n,m,l) = 1./(PSAo(n,m,l)/m1+PSMo(n,m,l)/m2+PSMn(n,m, &
                             l)/m3)/amu
 320        CONTINUE
 350     CONTINUE
 400  CONTINUE


!C  **
!C  SET POLE VALUES IF NECESSARY
!C  **
      zz = XX - 1
      DO 600 m = 1 , XX , zz
!C  **
         fp1 = 2.0
!C  **
         IF ( m  ==  1 ) THEN
            mmm = 1
            ifp = -1
         ENDIF
         IF ( m  ==  XX ) THEN
            mmm = 2
            ifp = 1
         ENDIF
         mp = m - (ifp*fp1) + 0.1
         DO 450 n = 1 , 15
            rmtp(n,mmm) = 0.
            psaop(n,mmm) = 0.
            psmop(n,mmm) = 0.
            psmnp(n,mmm) = 0.

            DO 420 l = 1 , 20
               rmtp(n,mmm) = rmtp(n,mmm) + ARMt(n,mp,l)
               psaop(n,mmm) = psaop(n,mmm) + PSAo(n,mp,l)
               psmop(n,mmm) = psmop(n,mmm) + PSMo(n,mp,l)
               psmnp(n,mmm) = psmnp(n,mmm) + PSMn(n,mp,l)
 420        CONTINUE

            rmtp(n,mmm) = rmtp(n,mmm)/20.
            psaop(n,mmm) = psaop(n,mmm)/20.
            psmop(n,mmm) = psmop(n,mmm)/20.
            psmnp(n,mmm) = psmnp(n,mmm)/20.
 450     CONTINUE

         DO 500 l = 1 , 20
            DO 460 n = 1 , 15
               ARMt(n,m,l) = rmtp(n,mmm)
               PSAo(n,m,l) = psaop(n,mmm)
               PSMo(n,m,l) = psmop(n,mmm)
               PSMn(n,m,l) = psmnp(n,mmm)
 460        CONTINUE
 500     CONTINUE
 600  CONTINUE


      DO 800 l = 1 , 20
         m = 90
         DO 650 n = 1 , 15
            ARMt(n,m,l) = (ARMt(n,m-1,l)+ARMt(n,m+1,l))/2.0
            PSAo(n,m,l) = (PSAo(n,m-1,l)+PSAo(n,m+1,l))/2.0
            PSMo(n,m,l) = (PSMo(n,m-1,l)+PSMo(n,m+1,l))/2.0
            PSMn(n,m,l) = (PSMn(n,m-1,l)+PSMn(n,m+1,l))/2.0
 650     CONTINUE

         m = 2
         DO 700 n = 1 , 15
            ARMt(n,m,l) = (ARMt(n,m+1,l)+ARMt(n,m-1,l))/2.0
            PSAo(n,m,l) = (PSAo(n,m+1,l)+PSAo(n,m-1,l))/2.0
            PSMo(n,m,l) = (PSMo(n,m+1,l)+PSMo(n,m-1,l))/2.0
            PSMn(n,m,l) = (PSMn(n,m+1,l)+PSMn(n,m-1,l))/2.0
 700     CONTINUE
 800  CONTINUE

      
RETURN
END SUBROUTINE neutral_composition_equation




SUBROUTINE OHHO2(OH,HO2)

      IMPLICIT NONE

      REAL*8 HO2 , ho2i , OH , ohi
      INTEGER I4 , n , ns
      PARAMETER (I4=65)
      DIMENSION OH(15) , HO2(15) , ohi(I4) , ho2i(I4)
      DATA ohi/9.5E+12 , 1.3E+13 , 1.6E+13 , 1.8E+13 , 2.0E+13 , &
           2.1E+13 , 2.2E+13 , 1.9E+13 , 1.3E+13 , 6.0E+12 , 2.1E+12 , &
           1.0E+12 , 3.0E+11 , 1.0E+11 , 4.0E+10 , 1.6E+10 , 7.0E+09 , &
           3.2E+09 , 1.2E+09 , 5.0E+08 , 2.0E+08 , 44*0.0/
      DATA ho2i/7.0E+12 , 9.0E+12 , 1.2E+13 , 1.3E+13 , 1.6E+13 , &
           1.7E+13 , 1.7E+13 , 1.3E+13 , 7.0E+12 , 2.5E+12 , 7.0E+11 , &
           1.5E+11 , 3.5E+10 , 8.0E+09 , 2.5E+09 , 9.0E+08 , 3.0E+08 , &
           1.4E+08 , 6.0E+07 , 2.5E+07 , 1.0E+07 , 44*0.0/

      DO 100 n = 1 , 15
         ns = n*4 + 3
         OH(n) = (ohi(ns)+ohi(ns+1))/2.
         HO2(n) = (ho2i(ns)+ho2i(ns+1))/2.0
 100  CONTINUE

RETURN
END SUBROUTINE OHHO2





SUBROUTINE JPROFIL(J,F107)

      IMPLICIT NONE

      REAL*8 a , b , c , F107 , fht
      INTEGER I3 , I4 , n , ns , nss
      PARAMETER (I4=65,I3=17)
      DIMENSION J(15) , ji(I4) , jii(I4) , fht(I3) , c(I3)
      REAL*8 J , ji , jii
      DATA jii/.4 , .47 , .53 , .66 , .78 , .92 , 1.1 , 1.3 , 1.5 , &
           1.7 , 2. , 2.4 , 3.0 , 3.7 , 4.5 , 5.5 , 6.8 , 8.2 , 10. , &
           12.5 , 15. , 19. , 23. , 28. , .34 , .42 , .51 , .62 , .77 , &
           .84 , .92 , 1. , 1.07 , 1.15 , 1.21 , 1.27 , 1.35 , 1.40 , &
           1.47 , 1.52 , 1.60 , 1.65 , 1.70 , 1.75 , 1.81 , 1.86 , &
           1.95 , 2.0 , 2.05 , 2.10 , 2.15 , 2.20 , 2.23 , 2.26 , 2.3 , &
           2.33 , 2.36 , 2.40 , 2.43 , 2.46 , 2.50 , 2.52 , 2.54 , &
           2.56 , 2.57/
      DATA fht/8*1.2 , 1.85 , 2.5 , 3.15 , 6*3.80/
      DATA c/8*0.900 , 0.680 , 0.43 , 0.18 , 6* - 0.066/

      DO 100 n = 1 , I4
         IF ( n <= 24 ) ji(n) = jii(n)*1.E-08
         IF ( n > 24 ) ji(n) = jii(n)*1.E-06
 100  CONTINUE

      DO 200 n = 1 , 15
         ns = n*4 + 3
         J(n) = (ji(ns)+ji(ns+1))/2.0
         nss = n + 1
         a = fht(nss+1)*0.6 + fht(nss)*0.4
         b = c(nss+1)*0.6 + c(nss)*0.4
         J(n) = J(n)*((a-1.0)*F107/176.+b)
 200  CONTINUE

RETURN
END SUBROUTINE JPROFIL






SUBROUTINE EDDY(NM,DK,HT)

      IMPLICIT NONE

      REAL*8 darg , dh , DK , dk0 , dkm , HT , ht0 , s1 , s2 , s3
      INTEGER n , NM
      DIMENSION DK(NM) , HT(NM)
      dk0 = 1.0E+02
      s1 = 0.030E-03
      s2 = 0.03E-03
      s3 = 0.05E-03
      dkm = 1.5E+02
      ht0 = 105.E+03
      DO 100 n = 1 , NM
         dh = HT(n) - ht0
         darg = -s3*dh**2*1.E-03
         IF ( darg < -600. ) darg = -600.
         IF ( dh <= 0.0 ) DK(n) = dk0*EXP(s1*dh) + (dkm-dk0) &
                                  *EXP(-s2*dh**2*1.E-03)
         IF ( dh > 0.0 ) DK(n) = dkm*EXP(darg)
         IF ( HT(n) > 150.E+03 ) DK(n) = 0.0
 100  CONTINUE

RETURN
END SUBROUTINE EDDY







!SUBROUTINE SPECIFIC_HEAT(P1, P2, P3, RMT, CP) ! RMT NOT USED
SUBROUTINE SPECIFIC_HEAT(P1, P2, P3, CP)

      IMPLICIT NONE
      REAL*8 :: c , CP , P1 , P2 , P3 ! , RMT NOT USED
      INTEGER::  n
      DIMENSION P1(15) , P2(15) , P3(15) , CP(15) , c(3)
      !DIMENSION RMT(15) ! RMT NOT USED
      REAL*8 :: m1 , m2 , m3
      DATA c/1298.9 , 909.26 , 1039.2/
      DATA m1 , m2 , m3/16. , 32. , 28./

      DO n = 1 , 15
         CP(N) = P1(N)*C(1)+P2(N)*C(2)+P3(N)*C(3)
      END DO

RETURN
END SUBROUTINE SPECIFIC_HEAT






SUBROUTINE TIROS(THMagd, ESSa, LL, QT, DEN, GW, EMAps, CMAps, PROfile, &
                       DMSp)

      IMPLICIT NONE

      REAL*8 ch , chi , CMAps , DEN , dfac , diff , DMSp , dprof , ed , &
           eflux , EMAps , ESSa , GW , PROfile , qdmsp , QT , ratio , &
           ri , rj , th
      REAL*8 THMagd
      INTEGER i , i1 , i2 , j1 , j2 , k , kk , l , ld , LL , n , nn
!*** End of declarations inserted by SPAG

!
! routine to get ionisation profile on 15 press levels from
! form data in common block ( generated form tiros data)
!  t. fuller-rowell may 85
!
      REAL*8 ionchr
!
      DIMENSION EMAps(21,20,7) , CMAps(21,20,7) , PROfile(15,21) , &
                ed(5) , QT(15) , ionchr(21) , DEN(15) , DMSp(21,20,5,5) &
                , dprof(4,5) , qdmsp(15)
!
!
      DATA ionchr/.378 , .458 , .616 , .773 , .913 , 1.088 , 1.403 , &
           1.718 , 2.033 , 2.349 , 2.979 , 3.610 , 4.250 , 4.780 , &
           6.130 , 7.392 , 8.653 , 9.914 , 12.436 , 14.957 , 17.479/
      DATA dprof/4.23E19 , 5.62E19 , 5.77E19 , 5.70E19 , 1.04E19 , &
           1.03E20 , 1.22E20 , 1.23E20 , 0.00E19 , 8.36E19 , 2.37E20 , &
           2.61E20 , 0.00E19 , 0.00E18 , 3.07E20 , 5.26E20 , 0.00E19 , &
           0.00E18 , 0.00E18 , 8.57E20/
      DATA qdmsp/15*0.0/
!c  **
!c  **
!c  set up dmsp data level 1 to 5
!c  **
      dfac = 1.0
      IF ( LL <= 3 ) ld = 1
      IF ( LL == 4 .OR. LL == 5 ) ld = 2
      IF ( LL == 6 .OR. LL == 7 ) ld = 3
      IF ( LL == 8 ) ld = 4
      IF ( LL == 9 ) ld = 5
      IF ( LL == 10 ) ld = 5
      IF ( LL == 10 ) dfac = GW/96.
      l = LL - 2
      IF ( l < 1 ) l = 1
      IF ( l > 7 ) l = 7

      ri = ESSa/18.0 + 11.
      i1 = ri
      ri = ri - i1
      IF ( i1 > 20 ) i1 = i1 - 20
      i2 = i1 + 1
      IF ( i2 > 20 ) i2 = i2 - 20
      th = ABS(THMagd) - 50.
      rj = th/2. + 1.
      j1 = rj
      rj = rj - j1
      j2 = j1 + 1

      eflux = rj*ri*EMAps(j2,i2,l) + (1.-rj)*ri*EMAps(j1,i2,l) &
              + rj*(1.-ri)*EMAps(j2,i1,l) + (1.-rj)*(1.-ri) &
              *EMAps(j1,i1,l)

      DO 100 i = 1 , 5
         ed(i) = rj*ri*DMSp(j2,i2,i,ld) + (1.-rj)*ri*DMSp(j1,i2,i,ld) &
                 + rj*(1.-ri)*DMSp(j2,i1,i,ld) + (1.-rj)*(1.-ri) &
                 *DMSp(j1,i1,i,ld)
 100  CONTINUE

      DO 200 nn = 1 , 4
         n = nn + 11
         qdmsp(n) = 0.0
         DO 150 i = 1 , 5
            qdmsp(n) = qdmsp(n) + ed(i)*dprof(nn,i)*DEN(n)
 150     CONTINUE

 200  CONTINUE


      ch = rj*ri*CMAps(j2,i2,l) + (1.-rj)*ri*CMAps(j1,i2,l) + rj*(1.-ri) &
           *CMAps(j2,i1,l) + (1.-rj)*(1.-ri)*CMAps(j1,i1,l)

      IF ( ch < 0.378 ) ch = 0.379
      IF ( ch > 17.479 ) WRITE (6,99001) ch
      eflux = 10.**(eflux)/1000.

      DO 300 kk = 2 , 21
         IF ( ch <= ionchr(kk) ) THEN
            k = kk - 1
            GOTO 400
         ENDIF
 300  CONTINUE

 400  chi = ch - ionchr(k)
      diff = ionchr(kk) - ionchr(k)
      ratio = chi/diff

      DO 500 n = 1 , 15
         QT(n) = (ratio*PROfile(n,kk)+(1.-ratio)*PROfile(n,k)) &
                 *eflux*DEN(n)*dfac
 500  CONTINUE

      DO 600 n = 12 , 15
         QT(n) = QT(n) + qdmsp(n)*dfac
 600  CONTINUE
!
      RETURN
99001 FORMAT ('  ch value outof bound in tiros',f10.6)

END SUBROUTINE TIROS





SUBROUTINE SOLAR_EUV(f107, RMT, SCHt, DEN, SZA, QEUv)

      IMPLICIT NONE

      REAL*8 b , cdo2 , cdrho , CHLGI , CSZa , dens
      REAL*8 rhi , rmwt , sdis , sht , hlg , sion , SZA
      INTEGER i , n 
      REAL*8, INTENT(IN) :: DEN(15) 
      REAL*8 :: SCHt(15) , QEUv(15) 
      REAL*8 fo2(15) , co2(68) , qsr(68) , crho(68) , qion(68)
      REAL*8 rmt(15)

      ! REAL*8 rn_explain(15) NOT USED

      REAL*8 f1(15) , f2(15)
      real*8 f107
      real*8 EUV_factor
      real*8 EUV_eff(15)
      real*8 SZA_factor
 
      DATA qsr/ - 39.206 , -39.218 , -39.242 , -39.271 , -39.302 , &
           -39.344 , -39.400 , -39.458 , -39.528 , -39.599 , -39.685 , &
           -39.795 , -39.917 , -40.057 , -40.211 , -40.365 , -40.531 , &
           -40.700 , -40.870 , -41.044 , -41.206 , -41.384 , -41.556 , &
           -41.725 , -41.898 , -42.077 , -42.259 , -42.440 , -42.604 , &
           -42.788 , -42.959 , -43.155 , -43.352 , -43.544 , -43.733 , &
           -43.893 , -44.043 , -44.211 , -44.355 , -44.509 , -44.675 , &
           -44.836 , -44.985 , -45.150 , -45.330 , -45.492 , -45.660 , &
           -45.827 , -45.994 , -46.160 , -46.322 , -46.480 , -46.642 , &
           -46.784 , -46.935 , -47.080 , -47.204 , -47.320 , -47.447 , &
           -47.559 , -47.667 , -47.742 , -47.813 , -47.879 , -47.936 , &
           -47.991 , -48.033 , -48.085/

      DATA qion/11.891 , 11.891 , 11.889 , 11.887 , 11.886 , 11.885 , &
                11.883 , 11.879 , 11.877 , 11.873 , 11.868 , 11.862 , &
                11.854 , 11.844 , 11.833 , 11.819 , 11.803 , 11.782 , &
                11.758 , 11.727 , 11.690 , 11.646 , 11.593 , 11.529 , &
                11.429 , 11.353 , 11.239 , 11.103 , 10.941 , 10.752 , &
                10.532 , 10.276 , 9.987 , 9.670 , 9.329 , 8.971 , 8.603 , &
                 8.228 , 7.846 , 7.459 , 7.067 , 6.678 , 6.294 , 5.918 , &
                 5.548 , 5.179 , 4.787 , 4.346 , 3.830 , 3.210 , 2.458 , &
                 1.545 , 0.439 , -0.905 , -2.592 , -4.653 , -7.161 , -10.214 , &
               -13.927 , -18.435 , -23.899 , -30.629 , -38.988 , -49.179 , &
               -61.590 , -76.690 , -95.037 , -117.3/

      DATA fo2/0.2096 , 0.2095 , 0.2090 , 0.1997 , 0.1758 , 0.1567 , &
           0.1326 , 0.1085 , 0.08394 , 0.05991 , 0.03922 , 0.02250 , &
           0.01100 , 0.004809 , 0.001963/
  
      DATA EUV_eff/6*0.3, 0.26, 0.45, 0.54, 0.37, 0.09, 0.11, 0.13, 0.16, 0.16/

      !===============================================================================
  
      if (f107 > 175.) then 
         EUV_factor = 0.00572*f107 + 1.376
      else
         EUV_factor = 0.01377*f107 - 0.03277
      endif

      co2(1) = 37.0
      crho(1) = -18.2

      !DO 100 n = 2 , 68
      DO n = 2 , 68
         co2(n) = co2(n-1) + 0.2
         crho(n) = crho(n-1) + 0.2
! 100  CONTINUE
      END DO

      CSZa = COS(SZA)
      SZA_factor = 1.0
      IF ( ABS(CSZa) <= 0.01 ) THEN
         SZA_factor = 0.5
         CSZa = CSZa + 0.015
      ENDIF
      SZA = ACOS(CSZa)

      DO 200 n = 2 , 14

         rmwt = rmt(n)*1.661E-24

         sht = SCHt(n)*1.0E+02
         dens = DEN(n)*1.0E-03

         !print *,'SOLAR_EUV : dens = ',dens

         IF ( CSZa < 0 ) THEN
            QEUv(n) = 0.0
         ELSE
            hlg = LOG(sht)
            rhi = LOG(dens) + hlg
            IF ( CSZa <= 0.5 ) THEN
!old function  cdrho = rhi + CHLGI(SZA,hlg)
               call function_chlgi(SZA,hlg,chlgi)
               cdrho = rhi + CHLGI
               b = hlg + LOG(rmwt/5.31E-23)
!old function  cdo2 = rhi + LOG(fo2(n)) + CHLGI(SZA,b) + 51.2899
               call function_chlgi(SZA,b,chlgi)
               cdo2 = rhi + LOG(fo2(n)) + CHLGI + 51.2899

            ELSE
               cdrho = rhi - LOG(CSZa)
               cdo2 = cdrho + LOG(fo2(n)) + 51.2899
            ENDIF

            !print *,'SOLAR_EUV : rhi, CSZa, cdrho = ', rhi, CSZa, cdrho
            IF ( cdrho <= crho(1) ) THEN
               F1(n) = EXP(qion(1))
            ELSEIF ( cdrho > crho(68) ) THEN
               F1(n) = 0.0
            ELSE
               i = INT(5.0*cdrho+92.0)
               !print *,'SOLAR_EUV : i, qion(i) = ', i, qion(i)
               sion = qion(i) + 5.0*(qion(i+1)-qion(i))*(cdrho-crho(i))
               F1(n) = EXP(sion)
            ENDIF

            IF ( cdo2 <= co2(1) ) THEN
               F2(n) = EXP(qsr(1))
            ELSEIF ( cdo2 > co2(68) ) THEN
               F2(n) = 0.0
            ELSE
               i = INT(5.0*cdo2-184.0)
               sdis = qsr(i) + 5.0*(qsr(i+1)-qsr(i))*(cdo2-co2(i))
               F2(n) = EXP(sdis)
            ENDIF

            QEUv(n) = SZA_factor*1.0E-04*(F1(n)*EUV_factor*EUV_eff(n)/0.3+fo2(n) &
                      *F2(n)/rmwt)
         ENDIF
 200  CONTINUE

RETURN

END SUBROUTINE SOLAR_EUV






SUBROUTINE INFRARED_COOLING(TEMP, HT, DEN, QIR)

      IMPLICIT NONE

      INTEGER n 
      REAL*8 TEMP(15) , DEN(15) , HT(15) , QIR(15)
      REAL*8 hgt , psi
      REAL*8 rn_explain(15)
 
      DATA rn_explain/10.794 , 11.160 , 11.222 , 11.396 , 11.668 , 11.541 , &
                      11.138 , 10.680 , 10.248 , 9.843 , 9.448 , 9.048 , &
                      8.634 , 8.206 , 7.763/

      DO 200 n = 2 , 14

         QIR(n) = -1.67E-18*(10**RN_explain(n))*EXP(-228.0/TEMp(n)) &
                  /(1.0+0.6*EXP(-228.0/TEMp(n))+0.2*EXP(-325.0/TEMp(n)))

         hgt = Ht(n)/1000.0

         QIR(n) = QIR(n)/(10.0*DEN(n))
         psi = EXP(-(hgt-80.0)/20.0)
         IF ( n <= 5 ) QIR(n) = QIR(n)*(1.0-psi)**4
 200  CONTINUE

RETURN

END SUBROUTINE INFRARED_COOLING








SUBROUTINE FUNCTION_CHLGI(XI, X, CHLGI)

      IMPLICIT NONE

      REAL*8 alf , ch , CHLGI , cs , css , P , sq , u , w , X , XI

      PARAMETER (P=1.5707963)

      u = XI - P
      IF ( u >= 0 ) THEN
         sq = 2.569E+04*EXP(-0.5*X)
         w = 0.707107*sq*TAN(u)
         IF ( w <= 3.0 ) THEN
            ch = 1.25331*sq*(2.0-EXP(-w*w)/(1.0+1.38622*w))
            CHLGI = LOG(ch)
            RETURN
         ENDIF
      ELSE
         cs = XI*XI
         css = cs*cs
         alf = 0.235296 - 0.81057/(20.7594-X)
         CHLGI = cs/(2.-0.23*cs-alf*css)
         RETURN
      ENDIF

      CHLGI = LOG(2.50664*sq)

RETURN
END SUBROUTINE FUNCTION_CHLGI






! REAL*8 P1(15) , P2(15) , P3(15)
! REAL*8 PI1(15) , PI2(15) , PI3(15)
! REAL*8 T(15), VIN(15), AMIn(15)
SUBROUTINE IONNEUT(P1, P2, P3, PI1, PI2, PI3, T, VIN, AMIn)

      IMPLICIT NONE

      REAL*8 a , AMIn , amu , b , factor , P1 , P2 , P3 , PI1 , PI2 , &
           sum , summol , T , v1 , v2 , VIN , PI3
      INTEGER n
      DIMENSION P1(15) , P2(15) , P3(15) , T(15) , VIN(15) , AMIn(15) , &
                a(3) , b(3) , PI1(15) , PI2(15) , PI3(15)
      REAL*8 mi1 , mi2 , mi3

      DATA mi1 , mi2 , mi3/16. , 30. , 32./
      DATA a/3.42E-11 , 6.66E-10 , 6.82E-10/
      DATA b/2.44E-10 , 4.28E-10 , 4.34E-10/

      !==========================================================================

      amu = 1.66E-27

      factor = 1.0

      ! For debugging :
      !print *,'IONNEUT : amu, factor = ',amu, factor
      !print *,'IONNEUT : p1 = ',p1
      !print *,'IONNEUT : p2 = ',p2
      !print *,'IONNEUT : p3 = ',p3
      !print *,'IONNEUT : pi1 = ',pi1
      !print *,'IONNEUT : pi2 = ',pi2
      !print *,'IONNEUT : pi3 = ',pi3
      !print *,'IONNEUT : t = ',t
      !print *,'IONNEUT :  '
 
      DO 100 n = 1 , 15
         summol = PI2(n) + PI3(n)
         sum = PI1(n) + PI2(n) + PI3(n)
         v2 = b(1)*P1(n) + b(2)*P2(n) + b(3)*P3(n)
         v1 = a(3)*P3(n) + a(2)*P2(n) + a(1)*P1(n)*factor*SQRT(T(n)) &
              *(1.08-0.139*LOG10(T(n))+4.51E-03*LOG10(T(n))**2)
         VIN(n) = (v1*PI1(n)+v2*summol)*1.E-06/sum
         AMIn(n) = (PI1(n)*mi1+PI2(n)*mi2+PI3(n)*mi3)*amu/sum
 100  CONTINUE

RETURN
END SUBROUTINE IONNEUT






SUBROUTINE FOSTER(EXNs, EYNs, newl, GW, AMP, potential)

      IMPLICIT NONE

        
      INTEGER LONSIZ , LATSIZ , POTSIZ , ACTIVL
      PARAMETER (LONSIZ=20, LATSIZ=22, POTSIZ=23, ACTIVL=7)


      REAL*8 :: EXNs(2,45,LONSIZ) , EYNs(2,45,LONSIZ) , GW
      INTEGER :: newl
      REAL*8 :: AMP 
      REAL*8 potential(POTSIZ,LONSIZ,ACTIVL)

      !-------------------
      ! Local variables
      !-------------------
      !REAL*8 DTR
      !PARAMETER (DTR=3.14159/180.)

      INTEGER :: newll 

      REAL*8 :: agwamp , tag , dag , gwamp , disty

      INTEGER :: i , j , k , m
      
! 
! newl now passed in to subroutine
! newll is reduced by 2 to 
! put power index in range      7th Feb 2008 TJFR
!
      newll=newl-2
!
! end change                    7th Feb 2008 TJFR
!
      gwamp = 1.0
      IF ( newll  ==  8 ) THEN
         newll = 7
         IF ( GW > 96. ) gwamp = GW/96.
      ENDIF
!
      agwamp = AMP*gwamp
!
!   ADF Dec 92: equation below.  220kV overall potential??
!         the other factor of 2 comes from averaging the potential
!               values
!
      tag = 440./agwamp

      DO 100 i = 2 , POTSIZ
         disty = 2001.5*SIN((i-1)*2.*DTR)
         dag = 2.*disty/agwamp
         m = i - 1
         DO 50 j = 1 , LONSIZ
            k = j - 1
            IF ( j == 1 ) k = LONSIZ
            EYNs(1,m,j) = -(potential(m,j,newll)+potential(i,j,newll) &
                        -potential(m,k,newll)-potential(i,k,newll))/dag
            EYNs(2,m,j) = EYNs(1,m,j)
            EXNs(1,m,j) = -(potential(i,j,newll)+potential(i,k,newll) &
                        -potential(m,j,newll)-potential(m,k,newll))/tag
            EXNs(2,m,j) = EXNs(1,m,j)
 50      CONTINUE
 100  CONTINUE

END SUBROUTINE FOSTER




!============================================================================



SUBROUTINE TIDES (nn, m, l, hough11, hough22, hough23, hough24, hough25, & ! input
                  amplitude, & ! input
                  ht, vx, vy, &  ! output 
                  DTIME, rmt, &  ! input
                  temp, &        ! output
                  tidalPhase, &  ! input
                  temp0, vy0, temp0av, ht0) ! input

!
! This subroutine calculates horizontal winds from the geopotential
! height variations which themselves are described by HOUGH modes.
! Geopotential height perturbations are introduced at level 1
! and recalculated for each time step. Horizontal winds form the
! model's lower boundary. Equations used for the winds correspond
! to those from Chapman & Lindzen (1970), p.109, Eq.18 and 19. They
! are also reproduced in Helen Parish's thesis, p.90.
!
! IMW, August 1995
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nn, m, l
      REAL*8, INTENT(IN) :: hough11(181), hough22(181), hough23(181), hough24(181), hough25(181)

      TYPE(amplitudeType), INTENT(IN) :: amplitude

      REAL*8, INTENT(INOUT) :: ht(15, 91, 20)
      REAL*8, INTENT(INOUT) :: vx(15, 91, 20)
      REAL*8, INTENT(INOUT) :: vy(15, 91, 20)

      ! If this gets changed to INTENT(IN), then the results are changed
      REAL*8, INTENT(INOUT) :: DTIME  ! the problem

      REAL*8, INTENT(IN) :: rmt(15,91,20)
      REAL*8, INTENT(INOUT) :: temp(15)

      TYPE(tidalPhaseType), INTENT(IN) :: tidalPhase 


      REAL*8, INTENT(IN) :: temp0(91)
      REAL*8, INTENT(IN) :: vy0(91)
      REAL*8, INTENT(IN) :: temp0av , ht0(91)

! Local variables -------------------------------------------------

      REAL*8 ssa, ssad
      REAL*8 ssa22, ssa23, ssa24, ssa11, ssa25
      REAL*8 lat, dlat
      !REAL*8 PI, R, OM, G, w1, w2, f, der22, der23, der24, derht, pht, ipht
      REAL*8 :: R, OM, G, w1, w2, f, der22, der23, der24, derht, pht, ipht
      REAL*8 phtd, iphtd, derhtd, der11, merhtd, dmerhtd, facd
      REAL*8 vxd(15,91,20), vyd(15,91,20), lambda11, shift11
      REAL*8 merht, dmerht, fac

      REAL*8 GSCON, H 

      REAL*8 lambda22 , lambda23 , lambda24 , a1 , a2 , a3 , a4 , a5

      REAL*8 shift22 , shift23 , shift24 , shift25
      REAL*8 der25 , lambda25

      !PARAMETER (PI=3.14159, R=6.371E+06, OM=7.292E-05, G=9.56)
      PARAMETER (R=6.371E+06, OM=7.292E-05, G=9.56)
      PARAMETER (GSCON=8.3141E+03)

! BEGIN CODE ============================================================

!-------------------------------
! lrm
! Initialize variables to 0.
!-------------------------------
vxd = 0.
vyd = 0.




!-------------------------------------------------------------------------------
! **       w1,w2       are tidal frequencies in 1/min (diurnal, semidiurnal)
! **      H      is scale height (considered constant here)
! **       dlat       is latitudinal step in rad
! **       f       is Coriolis factor
! **      pht      is perturbation height
! **      ipht      is perturbation height times i (needed for temperature)
! **            and is also multiplied by lambda
! **      ht      is height of level 1
! **      ssa      is default phase of tides modes.
!            for all modes and set to local midday, ie.max is there
! **      lambda      is vertical wave number (alpha in most literature)
! **      H is in km - but that is ok since the equivalent depths are
! **      used in km as well
! **    value of H calculated, using globally averaged temperature, temp0av
!-------------------------------------------------------------------------------
      H = 8.3141*temp0av/G/28.8
      !w1 = 1./86400.*2*PI*60
      !w2 = 2./86400.*2*PI*60  20110927 lrm
      w1 = 1./secondsInDay*2*PI*60
      w2 = 2./secondsInDay*2*PI*60


!------------------------------------------------------------------------------
! **      equivalent depth values are taken from Chapman&Lindzen(1970)
! **      pp 128,129,139. Values are in km (as is H).
!
! **      build in traps to eliminate cases where
! **      energy propagates downwards since we are
! **      not interested in those
!------------------------------------------------------------------------------
      a1 = (2./7.)*H/0.6909-(1./4.)
      a2 = (2./7.)*H/7.8519-(1./4.)
      a3 = (2./7.)*H/3.6665-(1./4.)
      a4 = (2./7.)*H/2.1098-(1./4.)
      a5 = (2./7.)*H/1.3671-(1./4.)

      !print *,'TIDES : a1 - a5 = ', a1, a2, a3, a4, a5
      

      IF (a1 < 0.) THEN
          lambda11 = 0.
          GOTO 40
      ENDIF

      lambda11 = SQRT(a1)
 40   IF (a2 < 0.) THEN
          lambda22 = 0.
          GOTO 50
      ENDIF

       lambda22 = SQRT(a2)
 50    IF (a3 < 0.) THEN
           lambda23 = 0.
           GOTO 60
      ENDIF

       lambda23 = SQRT(a3)
 60    IF (a4 < 0.) THEN
           lambda24 = 0.
           GOTO 70
      ENDIF

       lambda24 = SQRT(a4)
 70    IF (a5 < 0.) THEN
           lambda25 = 0.
           GOTO 80
      ENDIF

      lambda25 = SQRT(a5)
! 80   dlat = 4.* PI/180.
 80   dlat = 4.* DTR
      lat = (m-1)*2 - 90.
      !lat = lat * PI/180.
      lat = lat * DTR

      f = 2*OM*SIN(lat)*60

      !ssa = ((FLOAT(l-1)*18.0) + FLOAT(nn)*DTIME/240)*2.0
      ssad = ( (FLOAT(l-1)*18.0) + FLOAT(nn)*DTIME/240)*1.0
      ssa = ssad*2.0  ! lrm20130822

      !IF (ssa >= 360.0) ssa = ssa - 360
      !IF (ssa < 0.0) ssa = ssa + 360
      IF (ssa >= 360.0) then ! lrm20130822
          ssa = ssa - 360.
      ELSE IF (ssa < 0.0) THEN
           ssa = ssa + 360.
      ENDIF

      !IF (ssad >= 360.0) ssad = ssad - 360
      !IF (ssad < 0.0) ssad = ssad + 360
      IF (ssad >= 360.0) then ! lrm20130822
          ssad = ssad - 360.
      ELSE IF (ssad < 0.0) THEN
           ssad = ssad + 360.
      ENDIF



     !write(6,*) 'TIDES : nn, dtime, ssa ssad  ', nn, dtime, ssa , ssad
     !stop

! **
! **      ssa11,ssa22,ssa23,ssa24,ssa25 determine the phases of the various modes.
! **       In general, for a phase shift of 1 minute BACKwards (from 12.0h)
! **      need to ADD (!) the term {1*DTIME/240*2.0} to the original
! **      ssa value.
! **      e.g.: if I want phase at 11.2h instead of 12.0h need to
! **      determine difference (here: 11.2h = 11h+12 min = 12.0h - 48 min,
! **      therefore: use ssa+{48*DTIME/240*2.0} as phase term.
! **      set individual mode phases here from input.
! **      lt values from instream specify the local time
! **      of maximum and are given in the form eg.11.2,
! **      corresponding to 11 hours and 0.2*60= 12 minutes.
! **      IMW, Nov.1995
! **

      shift11 = (12. - tidalPhase%t_11)*60.
      shift22 = (12. - tidalPhase%t_22)*60.
      shift23 = (12. - tidalPhase%t_23)*60.
      shift24 = (12. - tidalPhase%t_24)*60.
      shift25 = (12. - tidalPhase%t_25)*60.

      ssa11 = ssad + shift11*DTIME/240*1.0
      ssa22 = ssa + shift22*DTIME/240*2.0
      ssa23 = ssa + shift23*DTIME/240*2.0
      ssa24 = ssa + shift24*DTIME/240*2.0
      ssa25 = ssa + shift25*DTIME/240*2.0

      !IF (ssa11 >= 360.0) ssa11 = ssa11 - 360
      !IF (ssa11 < 0.0) ssa11 = ssa11 + 360
      IF (ssa11 >= 360.0) ssa11 = ssa11 - 360
      IF (ssa11 < 0.0) ssa11 = ssa11 + 360

      IF (ssa22 >= 360.0) ssa22 = ssa22 - 360
      IF (ssa22 < 0.0) ssa22 = ssa22 + 360
      IF (ssa23 >= 360.0) ssa23 = ssa23 - 360
      IF (ssa23 < 0.0) ssa23 = ssa23 + 360
      IF (ssa24 >= 360.0) ssa24 = ssa24 - 360
      IF (ssa24 < 0.0) ssa24 = ssa24 + 360
      IF (ssa25 >= 360.0) ssa25 = ssa25 - 360
      IF (ssa25 < 0.0) ssa25 = ssa25 + 360

      !ssa11 = ssa11*PI/180.
      !ssa22 = ssa22*PI/180.
      !ssa23 = ssa23*PI/180.
      !ssa24 = ssa24*PI/180.
      !ssa25 = ssa25*PI/180.

      ssa11 = ssa11*DTR
      ssa22 = ssa22*DTR
      ssa23 = ssa23*DTR
      ssa24 = ssa24*DTR
      ssa25 = ssa25*DTR



      pht = hough22(m)*(COS(ssa22))*amplitude%a_22+hough24(m)*(COS(ssa24)) &
           *amplitude%a_24+hough23(m)*(COS(ssa23))*amplitude%a_23+hough25(m)* &
           (COS(ssa25))*amplitude%a_25

      phtd = hough11(m)*(COS(ssa11))*amplitude%a_11

      ipht = hough22(m)*(-SIN(ssa22))*amplitude%a_22*lambda22 + &
               hough23(m)*(-SIN(ssa23))*amplitude%a_23*lambda23 + &
               hough24(m)*(-SIN(ssa24))*amplitude%a_24*lambda24 + &
               hough25(m)*(-SIN(ssa25))*amplitude%a_25*lambda25

      iphtd = hough11(m)*(-SIN(ssa11))*amplitude%a_11*lambda11

      ht(1,m,l) = 80.E+3 + (pht+phtd) + ht0(m)

      !IF (m  ==  1 .OR. m  ==  91) GOTO 100
      IF (m  ==  1 .OR. m  ==  n_lats) GOTO 100

      der11 = (hough11(m+1)-hough11(m-1))/dlat
      der22 = (hough22(m+1)-hough22(m-1))/dlat
      der23 = (hough23(m+1)-hough23(m-1))/dlat
      der24 = (hough24(m+1)-hough24(m-1))/dlat
      der25 = (hough25(m+1)-hough25(m-1))/dlat

      derht = der22*(COS(ssa22))*amplitude%a_22+der23*(COS(ssa23))*amplitude%a_23 + &
              der24*(COS(ssa24))*amplitude%a_24+der25*(COS(ssa25))*amplitude%a_25
      derhtd =  der11*(COS(ssa11))*amplitude%a_11

! ** merht is the same as ipht but without lambda **
      merht = hough22(m)*(-SIN(ssa22))*amplitude%a_22+hough23(m)*(-SIN(ssa23))* &
              amplitude%a_23+hough24(m)*(-SIN(ssa24))*amplitude%a_24+hough25(m)* &
              (-SIN(ssa25))*amplitude%a_25
      merhtd =  hough11(m)*(-SIN(ssa11))*amplitude%a_11
      dmerht = der22*(-SIN(ssa22))*amplitude%a_22+der23*(-SIN(ssa23))* &
                  amplitude%a_23 + der24*(-SIN(ssa24))*amplitude%a_24 + &
                  der25*(-SIN(ssa25))*amplitude%a_25
      dmerhtd = der11*(-SIN(ssa11))*amplitude%a_11

! ** in fac a factor of 60 is necessary to get the correct conversion
! ** of velocity units into m/s
      fac = G*60/(w2*R*(1-f**2/w2**2))
      facd = G*60/(w1*R*(1-f**2/w1**2))
!
! ** calculate zonal wind **
! ** perturbation then added to the background wind value,
! ** vy0(m) (calculated in back.f)
       vy(1,m,l) = fac*(f/w2*derht-2/COS(lat)*pht)
      vyd(1,m,l) = facd*(f/w1*derhtd-1/COS(lat)*phtd)
       vy(1,m,l) = vy(1,m,l) + vyd(1,m,l)
       vy(1,m,l) = vy(1,m,l) + vy0(m)
!
! ** calculate meridional wind **
       vx(1,m,l) = fac*(2*f/(w2*COS(lat))*merht-dmerht)
      vxd(1,m,l) = facd*(1*f/(w1*COS(lat))*merhtd-dmerhtd)
       vx(1,m,l) = vx(1,m,l) + vxd(1,m,l)
!
! ** calculate temperature
!
! **      equation is
! **      taken from Fesen et.al, JGR, 1991, with some corrections of
! **      mistakes made in the paper (as confirmed by C.Fesen in 1995)
! **    perturbation is then added to the background temp value,
! **    temp0(m) (calculated in back.f)
!
!
      temp(1) = rmt(1,m,l)/GSCON * G * (ipht+iphtd)
      temp(1) = temp(1) + temp0(m)

      GOTO 200  ! return

 100  vx(1,1,l)  = 0.
      vx(1,91,l) = 0.
      vy(1,1,l)  = 0.
      vy(1,91,l) = 0.

      IF (m  ==  1)  temp(1) = temp0(2)
      IF (m  ==  91) temp(1) = temp0(90)


 200      RETURN
END SUBROUTINE TIDES





SUBROUTINE low_lat_efield(exns, eyns, ezns, plvu, zonal, f107, nday)

      IMPLICIT NONE

      ! Outputs
      REAL*8, intent(inout) :: exns(2,45,20), eyns(2,45,20), ezns(2,45,20)
      REAL*8, intent(inout) :: plvu(49), zonal(49)

      ! Inputs
      REAL*8, intent(in) :: f107
      INTEGER, intent(in) :: nday

      !------------------------------------------
      ! Local variables
      !------------------------------------------
      REAL*8 &
       !dtr, thmag, cmag, b, dip, essa, rll, &
       thmag, cmag, b, dip, essa, rll, &
       fac , factor, time , &
       plvu_dec_max(49), plvu_dec_min(49), &
       plvu_eq_max(49), plvu_eq_min(49), &
       plvu_jun_max(49), plvu_jun_min(49), &
       plvu_max(49), plvu_min(49), &
       zonal_dec_max(49), zonal_dec_min(49), &
       zonal_eq_max(49), zonal_eq_min(49), &
       zonal_jun_max(49), zonal_jun_min(49), &
       zonal_max(49), zonal_min(49)

      INTEGER ll, ll1, ll2, m, l, i

!g
!g  Below are American sector/December solstice/solar maximum drifts
!g
      data plvu_dec_max / &
      -15.,-25.,-28.,-28.,-35.,-35.,-35.,-35.,-32.,-20.,-10., -5., &
        0.,  5.,  5.,  8., 10., 12., 15., 16., 15., 12., 12., 13., &
       10., 10., 12., 12., 11., 10., 10.,  9.,  9.,  8.,  8.,  8., &
       12., 20., 30., 38., 30., 14., 10.,  8., -5.,-13.,-13.,-14., &
      -15./

!g
!g  Below are American sector/December solstice/solar minimum drifts.
!g
      data plvu_dec_min / &
      -17.,-13.,-10.,-10.,-10., -8., -8., -7., -8., -9., -8., -5., &
        0.,  3.,  4.,  5.,  9., 11., 14., 17., 19., 19., 17., 13., &
       10.,  7.,  5.,  4.,  3.,  0.,  0.,  0.,  0.,  2.,  3.,  3., &
        3.,  4.,  5.,  3.,  0., -2., -4., -5., -8.,-10.,-12.,-16., &
      -17./
!g
!g  Below are American sector/equinox/solar maximum drifts.
!g
      data plvu_eq_max / &
       -30.,-30.,-30.,-25.,-22.,-22.,-20.,-20.,-20.,-25.,-25.,-25., &
       -15.,-10.,  0.,  4.,  9., 10., 13., 15., 19., 21., 20., 18., &
        16., 15., 14., 12., 10., 10., 10., 12., 13., 15., 18., 21., &
        30., 45., 45., 15.,-15.,-30.,-30.,-30.,-28.,-25.,-30.,-33., &
       -30./
!g
!g  Below are American sector/equinox/solar minimum drifts.
!g
      data plvu_eq_min / &
        -18.,-18.,-18.,-18.,-15.,-15.,-14.,-10.,-10.,-10.,-10., -8., &
         -6., -2.,  0.,  5., 12., 15., 20., 20., 20., 21., 22., 22., &
         21., 20., 16., 12., 10.,  8.,  6.,  5.,  4.,  3.,  3.,  5., &
          5.,  5.,  3., -7.,-11.,-11.,-12.,-15.,-16.,-16.,-17.,-18., &
        -18. /
!g
!g  Below are American sector/June solstice/solar maximum drifts
!g
      data plvu_jun_max / &
        -25.,-23.,-22.,-21.,-20.,-20.,-20.,-20.,-18.,-16.,-12.,-10., &
          0.,  5.,  6.,  8., 10., 11., 12., 14., 15., 16., 17., 18., &
         18., 17., 17., 15., 15., 13., 12., 10., 10., 10., 11., 12., &
         17., 12., -3.,-15.,-20.,-25.,-28.,-29.,-30.,-28.,-28.,-27., &
        -25. /
!g
!g  Below are American sector/June solstice/solar minimum drifts
!g
      data plvu_jun_min / &
      -14.,-12.,-10.,-10.,-10., -8., -8., -8.,-8., -8., -8., -10., &
       -8., -4.,  0.,  2., 10., 12., 15., 15., 17., 18., 19., 20., &
       20., 19., 18., 17., 15., 13., 12., 10.,  9.,  7.,  4.,  0., &
       -3., -5.,-10.,-15.,-16.,-17.,-17.,-17.,-17.,-16.,-17.,-16., &
      -14./

!c
!c  **
!c  zonal ion drifts from bele fejer, from 0.25 LT every 0.5 hours
!c  first set December high solar activity.
!c
!     data zonal /88.3,83.3,75.6,70.9,65.3,69.1,60.9,57.2,37.8,
!    1 21.0,2.1,-9.1,-16.5,-19.2,-20.4,-24.6,-32.1,-37.7,-41.9,-44.5,
!    2 -49.1,-51.7,-55.9,-54.0,-50.6,-46.5,-42.2,-40.1,-34.1,-29.7,
!    3 -20.5,-13.4,-2.1,5.6,24.9,38.6,66.6,84.4,112.0,128.1,140.7,
!    4 144.5,140.1,135.8,125.0,113.9,102.6,92.3,88.3/

     do i = 1, 49
        read(41,*) time, zonal_jun_min(i)
     enddo

     do i = 1, 49
        read(41,*) time, zonal_jun_max(i)
     enddo

     do i = 1, 49
        read(41,*) time, zonal_dec_min(i)
     enddo

     do i = 1, 49
       read(41,*) time, zonal_dec_max(i)
     enddo

     do i = 1, 49
       read(41,*) time, zonal_eq_min(i)
     enddo

     do i = 1, 49
        read(41,*) time, zonal_eq_max(i)
     enddo

!g
!g  This subroutine has been edited so that the electric fields are
!g  automatically interpolated with respect to day number and f10.7
!g
!g  First interpolate in day no....
!g
      if (nday >= 355 .or. nday <= 82) then
          if (nday <= 82) then
             factor = float(nday+10)/92.
          else
             factor = float(nday-355)/92.
          endif

          do i = 1, 49
             plvu_max(i) = factor * (plvu_eq_max(i) - plvu_dec_max(i)) &
                        + plvu_dec_max(i)
             plvu_min(i) = factor * (plvu_eq_min(i) - plvu_dec_min(i)) &
                        + plvu_dec_min(i)
             zonal_max(i) = factor * (zonal_eq_max(i) - zonal_dec_max(i)) &
                         + zonal_dec_max(i)
             zonal_min(i) = factor * (zonal_eq_min(i) - zonal_dec_min(i)) &
                         + zonal_dec_min(i)
          enddo

      elseif (nday > 82 .and. nday <= 172) then

          factor = float(nday - 82) / 90.

          do i = 1, 49
             plvu_max(i) = factor * (plvu_jun_max(i) - plvu_eq_max(i)) &
                           + plvu_eq_max(i)
             plvu_min(i) = factor * (plvu_jun_min(i) - plvu_eq_min(i)) &
                           + plvu_eq_min(i)
             zonal_max(i) = factor * (zonal_jun_max(i) - zonal_eq_max(i)) &
                            + zonal_eq_max(i)
             zonal_min(i) = factor * (zonal_jun_min(i) - zonal_eq_min(i)) &
                            + zonal_eq_min(i)
          enddo

      elseif (nday > 172 .and. nday <= 266) then

          factor = float(nday - 172) / 94.

          do i = 1, 49
             plvu_max(i) = factor * (plvu_eq_max(i) - plvu_jun_max(i)) &
                           + plvu_jun_max(i)
             plvu_min(i) = factor * (plvu_eq_min(i) - plvu_jun_min(i)) &
                           + plvu_jun_min(i)
             zonal_max(i) = factor * (zonal_eq_max(i) - zonal_jun_max(i)) &
                            + zonal_jun_max(i)
             zonal_min(i) = factor * (zonal_eq_min(i) - zonal_jun_min(i)) &
                            + zonal_jun_min(i)
          enddo

      elseif (nday > 266 .and. nday < 355) then

          factor = float(nday - 266) / 89.

          do i = 1, 49
             plvu_max(i) = factor * (plvu_dec_max(i) - plvu_eq_max(i)) &
                           + plvu_eq_max(i)
             plvu_min(i) = factor * (plvu_dec_min(i) - plvu_eq_min(i)) &
                           + plvu_eq_min(i)
             zonal_max(i) = factor * (zonal_dec_max(i) - zonal_eq_max(i)) &
                            + zonal_eq_max(i)
             zonal_min(i) = factor * (zonal_dec_min(i) - zonal_eq_min(i)) &
                            + zonal_eq_min(i)
          enddo

      endif
!g
!g  Next interpolate in F10.7....
!g
      if (f107 >= 180.) then

          do i = 1, 49
             plvu(i) = plvu_max(i)
             zonal(i) = zonal_max(i)
          enddo

      elseif (f107 <= 70.) then

          do i = 1, 49
             plvu(i) = plvu_min(i)
             zonal(i) = zonal_min(i)
          enddo

      else

          factor = (f107 - 70.)/110.
          do i = 1, 49
             plvu(i) = factor * (plvu_max(i) - plvu_min(i)) &
                       + plvu_min(i)
             zonal(i) = factor * (zonal_max(i) - zonal_min(i)) &
                        + zonal_min(i)
          enddo

      endif

!c  **
!c  fill the arrays from 0 - 16 degrees latitude with values from
!c  Jicamarca

     !dtr = 3.14159/180.

!c  **
!c  m values are every 2 degrees in magnetic latitude from 2 degrees from
!c  magnetic pole.  Values run from 1 to 45 (mag equator).  e.g. 37 is 16
!c  degrees mag latitude, m=35 is 20 degrees mag latitude.
!c  **
     do 10 m = 37, 45
        thmag = (45.-m)*2.
        cmag = 90. - thmag
        b = 2.7e-5*sqrt(1.+3.*sin(thmag*dtr)**2)
        dip = atan(2.*cos(cmag*dtr)/sin(cmag*dtr))
        do 20 l = 1,20
           essa = (l-1.)*18.
           rll = essa/7.5+25.
           if (rll >= 49.0) rll = rll - 48.0
           ll = rll
           fac = rll-ll
           ll1 = ll
           ll2 = ll1+1
           eyns(1,m,l) = (plvu(ll1)*(1.-fac)+plvu(ll2)*fac)*b
           eyns(2,m,l) = eyns(1,m,l)
           rll = essa/7.5 + 24.5
           if (rll >= 49.0) rll = rll - 48.0
           ll = rll
           fac = rll - ll
           ll1 = ll
           ll2 = ll1 + 1
           exns(1,m,l) = (zonal(ll1)*(1.-fac) + zonal(ll2)*fac)*b*sin(dip)
           exns(2,m,l) = exns(1,m,l)
           ezns(1,m,l) = -(zonal(ll1)*(1.-fac) + zonal(ll2)*fac)*b*cos(dip)
           ezns(2,m,l) = ezns(1,m,l)
   20 continue
   10 continue


!c  **
!c  decrease e field to zero from 16 degrees mag lat to 30 degrees
!c  **
       do 30 l = 1, 20
       do 35 m = 23, 30
          exns(1,m,l) = 0.0
          exns(2,m,l) = 0.0
          eyns(1,m,l) = 0.0
          eyns(2,m,l) = 0.0
          ezns(1,m,l) = 0.0
          ezns(2,m,l) = 0.0
 35    continue

       do 40 m = 31, 36
          exns(1,m,l) = (m-30.)/7.*exns(1,37,l)
          exns(2,m,l) = exns(1,m,l)
          eyns(1,m,l) = (m-30.)/7.*eyns(1,37,l)
          eyns(2,m,l) = eyns(1,m,l)
          ezns(1,m,l) = (m-30.)/7.*ezns(1,37,l)
          ezns(2,m,l) = ezns(1,m,l)
   40 continue

   30 continue

RETURN
END SUBROUTINE low_lat_efield






SUBROUTINE Smooth_Tn_in_X(J, T, DELtha, S, NMIn, NMAx, MMIn, MMAx)

      IMPLICIT NONE
      !REAL*8 a1 , a2 , DELtha , dth , PI , S , T , th , tt
      REAL*8 a1 , a2 , DELtha , dth , S , T , th , tt

      INTEGER J , J1 , K1 , l , m , MMAx , MMIn , n , NMAx , NMIn
      !PARAMETER (J1=91, K1=22, PI=3.14159)
      PARAMETER (J1=91, K1=22)

      DIMENSION T(15,J,20) , tt(J1,K1)
      dth = DELtha/2.

      DO 200 n = NMIn , NMAx

          DO 50 l = 2 , 21
             DO 20 m = MMIn , MMAx
                tt(m,l) = T(n,m,l-1)
 20         CONTINUE
 50      CONTINUE

         DO 100 m = MMIn , MMAx
            tt(m,1) = T(n,m,20)
            tt(m,22) = T(n,m,1)
 100     CONTINUE

         DO 150 l = 2 , 21
            DO 120 m = MMIn + 1 , MMAx - 1
               th = PI - (m-1.)*dth
               a1 = COS(th)/SIN(th)
               a1 = a1*dth*(tt(m-1,l)-tt(m+1,l))/2.0
               a2 = tt(m-1,l) - 2.*tt(m,l) + tt(m+1,l)
               T(n,m,l-1) = tt(m,l) + S*(a1+a2)
 120        CONTINUE
 150     CONTINUE

 200  CONTINUE

RETURN
END SUBROUTINE Smooth_Tn_in_X






SUBROUTINE Smooth_Tn_in_Y(J, T, S, NMIn, NMAx, MMIn, MMAx)

      IMPLICIT NONE

      REAL*8 a1 , S , T , tt
      INTEGER J , J1 , K1 , l , m , MMAx , MMIn , n , NMAx , NMIn
      PARAMETER (J1=91, K1=22)
      DIMENSION T(15,J,20) , tt(J1,K1)


      DO 200 n = NMIn , NMAx
!
         DO 50 l = 2 , 21
            DO 20 m = MMIn , MMAx
               tt(m,l) = T(n,m,l-1)
 20         CONTINUE
 50      CONTINUE

         DO 100 m = MMIn , MMAx
            tt(m,1) = T(n,m,20)
            tt(m,22) = T(n,m,1)
 100     CONTINUE
         DO 150 l = 2 , 21
            DO 120 m = MMIn + 1 , MMAx - 1
               a1 = tt(m,l+1) - 2.*tt(m,l) + tt(m,l-1)
               T(n,m,l-1) = tt(m,l) + S*(a1)
 120        CONTINUE
 150     CONTINUE

 200  CONTINUE

      
RETURN
      
END SUBROUTINE Smooth_Tn_in_Y







SUBROUTINE Smooth_VnX_in_X(J, VX, VY, DELtha, DELphi, S, NMIn, NMAx, MMIn, MMAx)


      IMPLICIT NONE

      REAL*8 a1 , a2 , a3 , a4 , a5 , DELphi , DELtha , dphi , dth , &
             !dth2 , PI , S , sin2 , th , VX , vxt , VY , vyt
             dth2 , S , sin2 , th , VX , vxt , VY , vyt

      INTEGER J , J1 , K1 , l , m , MMAx , MMIn , n , NMAx , NMIn

      !PARAMETER (J1=91,K1=22,PI=3.14159)
      PARAMETER (J1=91,K1=22)

      DIMENSION VX(15,J,20) , VY(15,J,20) , vxt(J1,K1) , vyt(J1,K1)

      dth = DELtha/2.
      dphi = DELphi/2.

      DO 200 n = NMIn , NMAx

         DO 50 l = 2 , 21
            DO 20 m = MMIn , MMAx
               vxt(m,l) = VX(n,m,l-1)
               vyt(m,l) = VY(n,m,l-1)
 20         CONTINUE
 50      CONTINUE

         DO 100 m = MMIn , MMAx
            vxt(m,1) = VX(n,m,20)
            vyt(m,1) = VY(n,m,20)
            vxt(m,22) = VX(n,m,1)
            vyt(m,22) = VY(n,m,1)
 100     CONTINUE

         DO 150 l = 2 , 21
            DO 120 m = MMIn + 1 , MMAx - 1
               th = PI - (m-1.)*dth
               dth2 = dth*dth
               sin2 = SIN(th)**2
!cc      d2=dphi*dphi
               a1 = COS(th)/SIN(th)
               a1 = a1*dth*(vxt(m-1,l)-vxt(m+1,l))/2.0
               a2 = vxt(m-1,l) - 2.*vxt(m,l) + vxt(m+1,l)
               a3 = -vxt(m,l)*dth2/sin2
               a4 = -dth2*COS(th)*(vyt(m,l+1)-vyt(m,l-1))/sin2/dphi/2.0
               a5 = vxt(m,l)*dth2
               VX(n,m,l-1) = vxt(m,l) + S*(a1+a2+a3+a4+a5)
 120        CONTINUE
 150     CONTINUE

 200  CONTINUE

RETURN
END SUBROUTINE Smooth_VnX_in_X







SUBROUTINE Smooth_VnX_in_Y(J, VX, VY, DELtha, DELphi, S, NMIn, NMAx, MMIn, MMAx)

      IMPLICIT NONE

      REAL*8 &
       !a1 , a2 , a3 , d2 , DELphi , DELtha , dphi , dth , PI , S , &
       a1 , a2 , a3 , d2 , DELphi , DELtha , dphi , dth , S , &
           sin2 , th , VX , vxt , VY , vyt

      INTEGER J , J1 , K1 , l , m , MMAx , MMIn , n , NMAx , NMIn
      !PARAMETER (J1=91, K1=22, PI=3.14159)
      PARAMETER (J1=91, K1=22)

      DIMENSION VX(15,J,20) , VY(15,J,20) , vxt(J1,K1) , vyt(J1,K1)

      dth = DELtha/2.
      dphi = DELphi/2.

      DO 200 n = NMIn , NMAx

         DO 50 l = 2 , 21
            DO 20 m = MMIn , MMAx
               vxt(m,l) = VX(n,m,l-1)
               vyt(m,l) = VY(n,m,l-1)
 20         CONTINUE
 50      CONTINUE

         DO 100 m = MMIn , MMAx
            vxt(m,1) = VX(n,m,20)
            vyt(m,1) = VY(n,m,20)
            vxt(m,22) = VX(n,m,1)
            vyt(m,22) = VY(n,m,1)
 100     CONTINUE

         DO 150 l = 2 , 21
            DO 120 m = MMIn + 1 , MMAx - 1
               th = PI - (m-1.)*dth
               sin2 = SIN(th)**2
               d2 = dphi*dphi
               a1 = vxt(m,l+1) - 2.*vxt(m,l) + vxt(m,l-1)
               a2 = -COS(th)*(vyt(m,l+1)-vyt(m,l-1))*dphi/2.0
               a3 = vxt(m,l)*d2*sin2
               VX(n,m,l-1) = vxt(m,l) + S*(a1+a2+a3)
 120        CONTINUE
 150     CONTINUE

 200  CONTINUE

RETURN
END SUBROUTINE Smooth_VnX_in_Y







SUBROUTINE Smooth_VnY_in_X(J, VX, VY, DELtha, DELphi, S, NMIn, NMAx, MMIn, MMAx)

      IMPLICIT NONE

      REAL*8 &
       a1 , a2 , a3 , a4 , a5 , cth , DELphi , DELtha , dphi , dth , &
           !dth2 , PI , S , sin2 , th , VX , vxt , VY , vyt
           dth2 , S , sin2 , th , VX , vxt , VY , vyt

      INTEGER J , J1 , K1 , l , m , MMAx , MMIn , n , NMAx , NMIn

      !PARAMETER (J1=91,K1=22,PI=3.14159)
      PARAMETER (J1=91,K1=22)


      DIMENSION VX(15,J,20) , VY(15,J,20) , vxt(J1,K1) , vyt(J1,K1)

      !=====================================================================

      dth = DELtha/2.
      dphi = DELphi/2.

      DO 200 n = NMIn , NMAx

         DO 50 l = 2 , 21
            DO 20 m = MMIn , MMAx
               vxt(m,l) = VX(n,m,l-1)
               vyt(m,l) = VY(n,m,l-1)
 20         CONTINUE
 50      CONTINUE

         DO 100 m = MMIn , MMAx
            vxt(m,1) = VX(n,m,20)
            vyt(m,1) = VY(n,m,20)
            vxt(m,22) = VX(n,m,1)
            vyt(m,22) = VY(n,m,1)
 100     CONTINUE

         DO 150 l = 2 , 21
            DO 120 m = MMIn + 1 , MMAx - 1
               th = PI - (m-1.)*dth
               dth2 = dth*dth
               sin2 = SIN(th)**2
!c      d2=dphi*dphi
               cth = COS(th)
               a1 = cth*(vyt(m-1,l)-vyt(m+1,l))*dth/2./SIN(th)
               a2 = vyt(m-1,l) - 2.*vyt(m,l) + vyt(m+1,l)
               a3 = -vyt(m,l)*dth2/sin2
               a4 = cth*dth2*(vxt(m,l+1)-vxt(m,l-1))/2./dphi/sin2
               a5 = 2.*vyt(m,l)*dth2
               VY(n,m,l-1) = vyt(m,l) + S*(a1+a2+a3+a4+a5)
 120        CONTINUE
 150     CONTINUE

 200  CONTINUE

RETURN
END SUBROUTINE Smooth_VnY_in_X









SUBROUTINE Smooth_VnY_in_Y(J, VX, VY, DELtha, DELphi, S, NMIn, NMAx, MMIn, MMAx)

      IMPLICIT NONE

      REAL*8 &
           !a1 , a2 , cth , DELphi , DELtha , dphi , dth , PI , S , th , &
           a1 , a2 , cth , DELphi , DELtha , dphi , dth  , S , th , &
           VX , vxt , VY , vyt

      INTEGER J , J1 , K1 , l , m , MMAx , MMIn , n , NMAx , NMIn

      !PARAMETER (J1=91, K1=22, PI=3.14159)
      PARAMETER (J1=91, K1=22)

      DIMENSION VX(15,J,20) , VY(15,J,20) , vxt(J1,K1) , vyt(J1,K1)

      !====================================================================

      dth = DELtha/2.
      dphi = DELphi/2.

      DO 200 n = NMIn , NMAx
!
         DO 50 l = 2 , 21
            DO 20 m = MMIn , MMAx
               vxt(m,l) = VX(n,m,l-1)
               vyt(m,l) = VY(n,m,l-1)
 20         CONTINUE
 50      CONTINUE

         DO 100 m = MMIn , MMAx
            vxt(m,1) = VX(n,m,20)
            vyt(m,1) = VY(n,m,20)
            vxt(m,22) = VX(n,m,1)
            vyt(m,22) = VY(n,m,1)
 100     CONTINUE

         DO 150 l = 2 , 21
            DO 120 m = MMIn + 1 , MMAx - 1
               th = PI - (m-1.)*dth
               cth = COS(th)
               a1 = vyt(m,l+1) - 2.*vyt(m,l) + vyt(m,l-1)
               a2 = cth*dphi*(vxt(m,l+1)-vxt(m,l-1))/2.
               VY(n,m,l-1) = vyt(m,l) + S*(a1+a2)
 120        CONTINUE
 150     CONTINUE

 200  CONTINUE

RETURN
END SUBROUTINE Smooth_VnY_in_Y

!------------------------------------------------------


SUBROUTINE setThermoLatLons()

IMPLICIT NONE

 integer :: I
 real, parameter :: grid_res_lon = 18.
 real, parameter :: grid_res_lat = 2.
!=====================================================


therm_model_geo_long_deg = (/( (I-1)*grid_res_lon, I = 1, gt_lon_dim)  /)

therm_model_geo_lat_deg =  (/( (I-46)*grid_res_lat, I = 1, gt_lat_dim)  /)


      !-----------------------------------------------------------
      ! Set up lat, longs, used for interpolation to the IPE grid
      !-----------------------------------------------------------
      !do ii = 1 , GT_lon_dim
      !   therm_model_geo_long_deg(ii) = (ii - 1) * 18.
      !enddo

      !do ii = 1 , GT_lat_dim
      !   therm_model_geo_lat_deg(ii) = (ii - 46) * 2.
      !enddo


RETURN
END SUBROUTINE setThermoLatLons


END MODULE moduleTHERMOSPHERE
