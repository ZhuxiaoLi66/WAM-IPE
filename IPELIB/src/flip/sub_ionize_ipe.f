!      program IONIZE_IPE
      subroutine IONIZE_IPE ( 
     & flipdim,z,gr,on,o2n,n2n,hn,hen,tn
     &,gm_lat,mlt
     &,tiros_activity_level,GW
     &,qiont_O,qiont_O2,qiont_N2)
      implicit none
!
      CHARACTER (LEN=*), parameter :: filename0='N'
      CHARACTER (LEN=*), parameter :: filename1='ipe_nh1_data'
      CHARACTER (LEN=*), parameter :: filename2='ipe_nh2_data'
      INTEGER, parameter :: UNIT1=101
      INTEGER, parameter :: UNIT2=102
      CHARACTER (LEN=100) :: string_dum
      INTEGER :: istat, i
      INTEGER,INTENT(IN) :: tiros_activity_level
      REAL,INTENT(IN) :: gw
!nm      INTEGER, parameter :: FLIPDIM=228
      INTEGER, INTENT(IN) :: FLIPDIM
      INTEGER            :: imax1 !NH=6; SH=4
      INTEGER, parameter :: imax2=2
      REAL*8,dimension(flipdim), INTENT(IN) ::Z,GR,ON,O2N,N2N,HN,HEN,TN
      REAL, INTENT(IN) :: gm_lat,mlt
      REAL ::  SL,GL(flipdim),BM,SZA,N4S
      REAL :: UN,NNO,EHT,TI,TE,OP,HP,MINP,HEP
     &,PHION,PRODOP,NP,EQN2D,NPLSPRD
      REAL*8,dimension(flipdim) :: qiont_total
      REAL*8,dimension(flipdim),INTENT(OUT) :: qiont_O,qiont_O2,qiont_N2

! before the call to tiros_ionize read in the TIROS energy influx EMAPS
! characteristic energy CMAPS, and TIROS spectra DJSPECTRA
! read in emaps and cmaps and spectra
!
!      call tiros_init_ipe()
!
! tiros_ionize returns ionization rates for O, O2, and N2 for a given 
! geomagnetic latitude GL and magnetic local time MLT based on 
! TIROS/NOAA statistical maps of energy influx and characteristic energy
! the ionization rates assumes an observed spectrum based on the 
! characteristic energy CH at the location, the TIROS spectrum is 
! between 300eV - 100keV, below 300eV assumes a Maxwellian where 
! the average energy of the distribution is CH
! the model atmosphere should be provided by calling program (e.g., IPE)
! a sample model profile is read in from ipe_nh1_data and ipe_nh2_data
!
! input 
! sample IPE profiles ipe_nh1_data and ipe_nh2_data
! geomagnetic latitude GL in radians
! use the geomagnetic latitude of the foot point 
!nm       gm_lat = GL(1)
! magnetic local time MLT in hours
! set the power index 1 to 10
!nm       tiros_activity_level = 7
! set the number of gigawatts of auroral power (used if LL = 10)
!       GW = 140.
! set the magnetic local time as in sample neutral atmosphere profile
!nm       mlt = 0.168238
! output all units number/m3/s
! qiont_total: total ionization rate from aurora
! qiont_0 total atomic oxygen ionization rate from aurora
! qiont_02 total molecular oxygen ionization rate from aurora
! qiont_N2 total molecular nitrogen ionization rate from aurora
!
      call tiros_ionize_ipe(flipdim,z,gr,on,o2n,n2n,hn,hen,tn
     &,gm_lat,mlt,
     &tiros_activity_level,GW,qiont_total,qiont_O,qiont_O2,qiont_N2)
!
      return 
      end subroutine IONIZE_IPE

