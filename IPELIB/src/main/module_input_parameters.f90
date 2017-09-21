!file1: IPEOptization; file2: wam-ipe
!note:20120207: v36: used only activating the perp.transport gradually...
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
      MODULE module_input_parameters
      USE module_precision
      USE module_IPE_dimension,ONLY: NLP,NMP,NPTS2D
      IMPLICIT NONE

!--- IPE wide run parameters
      INTEGER (KIND=int_prec), PUBLIC   :: utime           !UT[sec] IPE internal time management
      INTEGER (KIND=int_prec), PUBLIC   :: nTimeStep=1     !internal number of time steps
      INTEGER (KIND=int_prec), PUBLIC   :: start_time      !=0  !UT[sec]
      INTEGER (KIND=int_prec), PUBLIC   :: stop_time       !=60 !UT[sec]
      INTEGER (KIND=int_prec), PUBLIC   :: time_step=300       !=60 ![sec]
      REAL (KIND=real_prec), PUBLIC     :: dumpFrequency=3600   ! [sec]
      INTEGER (KIND=int_prec), PUBLIC   :: nprocs=1        !Number of processors
      INTEGER (KIND=int_prec), PUBLIC   :: mype=0          !Processor number
      INTEGER (KIND=int_prec), PUBLIC   :: lps,lpe,mps,mpe !Per processor start and stop indexes for lp,mp
      INTEGER (KIND=int_prec), PUBLIC   :: lpHaloSize=99   !lp halo size (big number=NOP for serial)
      INTEGER (KIND=int_prec), PUBLIC   :: mpHaloSize=99   !mp halo size (big number=NOP for serial)
      INTEGER (KIND=int_prec), PUBLIC   :: MaxLpHaloUsed=0 !Max lp halo size used for the entire run
      INTEGER (KIND=int_prec), PUBLIC   :: MaxMpHaloUsed=0 !Max mp halo size used for the entire run

      REAL (KIND=real_prec), PUBLIC :: F107D   !.. Daily F10.7
      REAL (KIND=real_prec), PUBLIC :: F107AV  !.. 81 day average F10.7
!
      INTEGER (KIND=int_prec), PUBLIC :: NYEAR ! year
      INTEGER (KIND=int_prec), PUBLIC :: NDAY  ! day number

      INTEGER (KIND=int_prec), PUBLIC :: internalTimeLoopMax=1  ![times]internal time loop: default 1
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_output=900  ![sec] must be multiple of time_step: default 15m
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_msis=180    !frequency[sec] to call MSIS/HWM: default 3m
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_plasma=60   !frequency[sec] to call plasma: default 1m
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_paraTrans=300!frequency[sec] to call FLIP: default 5min
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_eldyn=180   !frequency[sec] to call eldyn: default 3m(for quiet climatology),60s for storm

      LOGICAL                , PUBLIC :: parallelBuild=.false.

!--- FLIP specific input parameters
      REAL (KIND=real_prec), PUBLIC :: DTMIN_flip=1.0  !.. Minimum time step allowed (&=10 secs?)
      INTEGER (KIND=int_prec),PUBLIC :: sw_INNO=-1  !.. switch to turn on FLIP NO calculation if <0
      REAL (KIND=real_prec), PUBLIC :: FPAS_flip=0.0   !.. Pitch angle scattering fraction
      REAL (KIND=real_prec), PUBLIC :: HPEQ_flip=0.0   !.. Sets initial equatorial [H+][cm-3] if positive
      REAL (KIND=real_prec), PUBLIC :: HEPRAT_flip=0.09 !.. Initial He+/H+ ratio (.01 to 1.0)
      REAL (KIND=real_prec), PUBLIC :: COLFAC_flip=1.7 !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)

      INTEGER (KIND=int_prec),PUBLIC :: sw_PE2S=1    !.. switches photoelectron solutions ON if > 0  !dbg20141210
      INTEGER (KIND=int_prec),PUBLIC :: sw_TEI=1    !.. switches Te/Ti solutions ON if > 0
      INTEGER (KIND=int_prec),PUBLIC :: sw_OHPLS=1  !.. switches O+/H+ solutions ON if > 0

      INTEGER (KIND=int_prec),PUBLIC :: sw_IHEPLS=1 !.. switches He+ diffusive solutions on if > 0
      INTEGER (KIND=int_prec),PUBLIC :: sw_INPLS=1  !.. switches N+ diffusive solutions on if > 0
      INTEGER (KIND=int_prec),PUBLIC :: sw_wind_flip=1  !.. switch for neutral wind input to FLIP:
!:switch ON: default HWM93 
! +1 :add constant value(fac_wind_flip) for test
! 0  :multiply factor(fac_wind_flip): default is ZERO wind to switch OFF wind
      REAL (KIND=real_prec), PUBLIC :: fac_wind_flip = 0.00000
! 2  :assign constant value(fac_wind_flip) for the entire field aligned wind (UNX in flux_tube_solver.f90) 
!-1  :assign the constant value in Vn_ms1 instead of HWM93 in module_neutral.f90
!     note: FLIP assumes positive SOUTHWARD along a field line
      INTEGER (KIND=int_prec),PUBLIC :: sw_depleted_flip=0  !.. switch for depleted flux tube in FLIP: 1:ON; 0:OFF 
      INTEGER (KIND=int_prec),PUBLIC :: start_time_depleted !.. time UT to start to deplete the flux tube
      INTEGER,PUBLIC :: sw_neutral_heating_flip  !.. switch for neutral heating calculation in FLIP: 1:ON; 0:OFF
      REAL (KIND=real_prec), PUBLIC :: init_Te_max !.. max Te[K] in the initial profile
      INTEGER, PUBLIC :: sw_DEBUG_flip           !.. switch to turn on debug writes:0=off; 1=on for solver
      INTEGER, PUBLIC :: sw_ERSTOP_flip          !.. switch to turn on STOP in SUB-WRITE_EFLAG
!dbg20120301: N+ BAND solver issue
      LOGICAL, PUBLIC :: sw_LCE  !local chemical equilibruium below ht_LCE[km]
      REAL (KIND=real_prec), PUBLIC :: ht_LCE !.. max ht[km] for LCE
      REAL (KIND=real_prec), PUBLIC :: ZLBNP_inp !.. ZLBNP
!dbg20120304:
      REAL (KIND=real_prec), PUBLIC :: FNFAC_flip !.. FNFAC in RSPRIM.FOR
!dbg20121129
      LOGICAL, PUBLIC :: sw_optw_flip !=F  !chemical routine is called before He+ solution for too inflated o+ density due to exb drift
!dbg20121130
      LOGICAL, PUBLIC :: sw_init_guess_flip !=F  !this might help in finding a solution for convergence error???
      INTEGER (KIND=int_prec), PUBLIC :: dt_init_guess_flip=60 !max DT for changing init_guess 
!dbg20121130
      REAL (KIND=real_prec), PUBLIC :: ZLBDY_flip=120.  !Lower boundary altitude

!--- MSIS/HWM specific input parameters
      REAL (KIND=real_prec), DIMENSION(7), PUBLIC :: AP   ! magnetic index(daily)
!.. MSIS: or when sw(9)=-1. :                 
!           - array containing:                                         
!             (1) daily ap                                              
!             (2) 3 hr ap index for current time                        
!             (3) 3 hr ap index for 3 hrs before current time           
!             (4) 3 hr ap index for 6 hrs before current time           
!             (5) 3 hr ap index for 9 hrs before current time           
!             (6) average of eight 3 hr ap indicies from 12 to 33 hrs pr
!                    to current time                                    
!             (7) average of eight 3 hr ap indicies from 36 to 57 hrs pr
!                    to current time 
!.. HWM:        ap - two element array with                                    
!             ap(1) = magnetic index(daily) (use 4 in lower atmos.)     
!             ap(2)=current 3hr ap index (used only when sw(9)=-1.) 
!
!--- ELDYN specific input parameters
      REAL (KIND=real_prec), PUBLIC :: kp_eld      =  1.0    ! geomagnetic index
      LOGICAL, PUBLIC :: sw_bnd_wei=.false.
      REAL (KIND=real_prec), PUBLIC :: bnd_wei_eld = 44.   ! weimer boundary setting
      REAL (KIND=real_prec), PUBLIC :: lat_sft_eld = 54.   ! weimer boundary setting
! for now inputs available every minute for 24hrs
!     INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: nLevPI = 1440   !=60min/hr*24hr/dy
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: nLevPI = 6000
      INTEGER (KIND=int_prec),            PUBLIC :: LPI    =  1     ! time(minute) index for magnetic indices
! weimer inputs:(0) default setting: bz/y
      REAL (KIND=real_prec), DIMENSION(nLevPI), PUBLIC :: bz_eld!      = -0.3163809  ! geomagnetic index
      REAL (KIND=real_prec), DIMENSION(nLevPI), PUBLIC :: by_eld!      =  0.000  ! geomagnetic index

! weimer inputs:(1) solar wind parameters from ctip input when sw_ctip_input is ON:
      LOGICAL, PUBLIC :: sw_ctip_input=.false.
      INTEGER (KIND=int_prec), PUBLIC   :: utime0LPI=518400 !start time UT[sec] of the ctip input parameters  
      REAL (KIND=real_prec), DIMENSION(nLevPI), PUBLIC :: swbt! 4.98  !solar wind total magnetic field
      REAL (KIND=real_prec), DIMENSION(nLevPI), PUBLIC :: swangle! 80.88 !solar wind angle
      REAL (KIND=real_prec), DIMENSION(nLevPI), PUBLIC :: swvel! 449.47 !solar wind velocity

! tiros inputs [Foster et al 1986grl]
      INTEGER (KIND=int_prec),  DIMENSION(nLevPI), PUBLIC :: LevPI     ! Power Index 
      REAL    (KIND=real_prec), DIMENSION(nLevPI), PUBLIC :: GWatts!  1.7   ! Power Input [GW]


!--- all the SWITCHes either integer or logical or character
      LOGICAL, PUBLIC :: sw_debug=.false.
      LOGICAL, PUBLIC :: sw_debug_mpi=.false.
      LOGICAL, PUBLIC :: sw_output_fort167=.false.
      LOGICAL, PUBLIC :: sw_output_wind    =.false. !unit=6000,6001
      LOGICAL, PUBLIC :: barriersOn=.false. !true means turn on barriers.
      LOGICAL, PUBLIC :: sw_use_wam_fields_for_restart=.true. !unit=5000,5001

      INTEGER(KIND=int_prec), PUBLIC :: peFort167=0 !default mype=0
      INTEGER(KIND=int_prec), PUBLIC :: mpfort167=10
      INTEGER(KIND=int_prec), PUBLIC :: lpfort167=14
      INTEGER(KIND=int_prec), DIMENSION(2), PUBLIC :: iout
      INTEGER(KIND=int_prec), PUBLIC :: mpstop=80
      INTEGER(KIND=int_prec), PUBLIC :: sw_neutral=1    
!0: WAM debug: use which ever ESMF fields are coming across for debugging purpose
!1: WAM default science mode: specify ESMF fields you wish to use
!2: GT
!3: MSIS(default)
!4: read in files
      LOGICAL, dimension(7), PUBLIC :: swNeuPar!     =.true. !f:OFF (from MSIS); t:ON (from WAM)
!determines which neutral parameters to derive from WAM only when sw_neutral=0/1? 
!1:tn; 2:un1(east); 3:un2(north); 4:un3(up); 5:[O]; 6:[O2]; 7:[N2]
      LOGICAL, PUBLIC :: swEsmfTime =.false.
      INTEGER(KIND=int_prec), PUBLIC :: sw_eldyn
!0:self-consistent eldyn solver; 1:WACCM efield ;2:  ;3: read in external efield
      INTEGER(KIND=int_prec), PUBLIC :: sw_aurora=1
!0:no aurora; 1:tiros; 2:read-in; 3: mhd
      INTEGER(KIND=int_prec), PUBLIC :: sw_pcp        !0:heelis; 1:weimer
      INTEGER(KIND=int_prec), PUBLIC :: sw_grid       !0:APEX; 1:FLIP
! if sw_grid=1 
!dbg20120304: nolonger used
!nm20120304      REAL (KIND=real_prec), PUBLIC :: PCO_flip  
!nm20120304      REAL (KIND=real_prec), PUBLIC :: BLON_flip 
      LOGICAL, PUBLIC :: sw_output_plasma_grid
!JFM  LOGICAL, PUBLIC :: sw_rw_sw_perp_trans
      LOGICAL, PUBLIC :: sw_dbg_perp_trans
      INTEGER(KIND=int_prec), PUBLIC :: sw_perp_transport 
!0:WITHOUT perpendicular transport
!1:THETA only transport included
!2:both THETA&PHI:transport included, NH/SH flux tubes are moving together with the same ExB drift
!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
! if sw_perp_tr=>1
      INTEGER (KIND=int_prec), PUBLIC :: lpmin_perp_trans !=15 :mlatN=78deg???
      INTEGER (KIND=int_prec), PUBLIC :: lpmax_perp_trans !=151:mlatN=5.64deg
      INTEGER (KIND=int_prec), PUBLIC :: sw_th_or_r
!0:th method (ctipe/shawn)
!1:R method (gip)
      INTEGER (KIND=int_prec), PUBLIC :: record_number_plasma_start
      INTEGER (KIND=int_prec), PUBLIC :: sw_record_number
!nm20160329: used only when HPEQ_flip=0.5
      INTEGER (KIND=int_prec), PUBLIC :: ut_start_perp_trans=432000
      INTEGER (KIND=int_prec), PUBLIC :: duration=86400 !used when sw_record_n=1
      INTEGER (KIND=int_prec), PUBLIC :: sw_exb_up
! (0) self consistent electrodynamics
! (1) WACCM E empirical model
! (option) multiply a factor: default is 1.0
      REAL(KIND=real_prec), PUBLIC :: fac_exb_up = 1.0
! (2) GIP empirical model
! (3) SUPIM empirical model
! (4) zero
      INTEGER(KIND=int_prec), PUBLIC :: sw_para_transport 
!0:WITHOUT parallel transport (no calling to flux tube solver)
!1:parallel transport included
      INTEGER(KIND=int_prec), PUBLIC :: sw_ksi
!0: ksi_factor=1.0---no compressional effect/no adiabatic heating
!1: ksi_factor from richards thesis---including compressional effect/adiabatic heating ,,,used until 20120314 with interpolate_ft.v16.f90
!2:20120330: new way of calculating the ksi_factor with interpolate_ft.v17.f90
      INTEGER(KIND=int_prec), PUBLIC :: sw_divv
!0: div * V//=0
!1: div * V// included in the Te/i solver
!dbg20120313 
      REAL   (KIND=real_prec), PUBLIC :: fac_BM
      INTEGER(KIND=int_prec) , PUBLIC :: SMScomm,sendCount,NumPolevalProcs
      INTEGER, PUBLIC :: MPI_COMM_IPE        

! --- NUOPC cap input parameters
      INTEGER,                 PARAMETER :: str_len_max     = 256            ! max string length for VTK base name
      REAL(KIND=real_prec),       PUBLIC :: mesh_height_min =   0._real_prec !  min mesh height (km)
      REAL(KIND=real_prec),       PUBLIC :: mesh_height_max = 782._real_prec !  max mesh height (km)
      INTEGER,                    PUBLIC :: mesh_write      = 0              ! write mesh to VTK file(s): 1=yes, 0=no
      CHARACTER(LEN=str_len_max), PUBLIC :: mesh_write_file = 'ipemesh'      !  default base name for VTK file(s)
!
!---

      NAMELIST/IPEDIMS/NLP,NMP,NPTS2D
      NAMELIST/NMIPE/start_time &
     &,stop_time &
     &,time_step &
     &,dumpFrequency &
     &,F107D   &
     &,F107AV  &
     &,NYEAR  &
     &,NDAY   &
     &,internalTimeLoopMax &
     &,ip_freq_eldyn &
     &,ip_freq_output &
     &,ip_freq_msis &
     &,ip_freq_plasma &
     &,ip_freq_paraTrans
      NAMELIST/NMFLIP/DTMIN_flip  & 
     &,sw_INNO   & 
     &,FPAS_flip   & 
     &,HPEQ_flip   & 
     &,HEPRAT_flip & 
     &,COLFAC_flip & 
     &,sw_PE2S &
     &,sw_TEI &
     &,sw_OHPLS &
     &,sw_IHEPLS &
     &,sw_INPLS  &
     &,sw_wind_flip &
     &,fac_wind_flip &
     &,sw_depleted_flip &
     &,start_time_depleted &
     &,sw_neutral_heating_flip &
     &,init_Te_max &
     &,sw_DEBUG_flip &
     &,sw_ERSTOP_flip &
     &,sw_LCE &
     &,ht_LCE &
     &,ZLBNP_inp &
     &,FNFAC_flip &
     &,sw_optw_flip &
     &,sw_init_guess_flip &
     &,dt_init_guess_flip &
     &,ZLBDY_flip 
      NAMELIST/NMSWITCH/&
           &  sw_neutral     &
           &, swNeuPar       &
           &, swEsmfTime     &
           &, sw_eldyn       &
           &, sw_aurora      &
           &, sw_ctip_input  &
           &, utime0LPI      &
           &, sw_pcp         &
           &, sw_grid        &
           &, sw_output_plasma_grid        &
!JFM       &, sw_rw_sw_perp_trans &
           &, sw_dbg_perp_trans &
           &, sw_perp_transport &
           &, lpmin_perp_trans &
           &, lpmax_perp_trans &
           &, sw_th_or_r &
           &, sw_exb_up &
           &, fac_exb_up &
           &, sw_para_transport &
           &, sw_ksi &
           &, sw_divv &
           &, mpstop  &
           &, sw_debug       &
           &, sw_debug_mpi   &
           &, sw_output_fort167   &
           &, sw_output_wind   &
           &, sw_use_wam_fields_for_restart   & !nm20170728temporary commented out
           &, mpfort167   &
           &, lpfort167   &
           &, peFort167   &
           &, record_number_plasma_start   &
           &, sw_record_number   &
           &, ut_start_perp_trans   &
           &, duration   &
           &, fac_BM   &
           &, iout     &
           &, barriersOn
!nm20120304           &, PCO_flip       &
!nm20120304           &, BLON_flip      &
      NAMELIST/NMMSIS/AP  
      NAMELIST/NMWEIM/  &
     & bnd_wei_eld &
     &,by_eld &
     &,bz_eld &
     &,kp_eld &
     &,lat_sft_eld &
     &,swangle &
     &,swbt &
     &,swvel
      NAMELIST/NMTIROS/  &
     & GWatts &
     &,LevPI
     NAMELIST/IPECAP/ &
              mesh_height_min, &
              mesh_height_max, &
              mesh_write,      &
              mesh_write_file


      PRIVATE
      PUBLIC :: read_input_parameters


      CONTAINS
!---------------------------
! initialise plasma grids
        SUBROUTINE read_input_parameters ( )
        USE module_IPE_dimension,ONLY: NLP,NMP,NPTS2D
        IMPLICIT NONE
!MPI requirement 
      ! Joe : July 19, 2017 : This causes an error if serial compilation
      ! is desired. Should have preprocessing flags around it.
!SMS$INSERT         include "mpif.h"
!---------
        INTEGER(KIND=int_prec),PARAMETER :: LUN_nmlt=1,LUN_nmlt2=2
        CHARACTER(LEN=*),PARAMETER :: INPTNMLT='IPE.inp'
        CHARACTER(LEN=*),PARAMETER :: INPTNMLT2='IPEsw.inp'
        INTEGER(KIND=int_prec) :: IOST_OP=0
        INTEGER(KIND=int_prec) :: IOST_RD=0
        INTEGER (KIND=int_prec), PARAMETER :: LUN_LOG0=10  !output4input parameters only
        CHARACTER (LEN=*), PARAMETER :: filename='logfile_input_params.log'
        INTEGER (KIND=int_prec) :: istat        
!dbg20160408 sms debug
        INTEGER (KIND=int_prec) :: nElements,ierr
        INTEGER (KIND=int_prec) :: mycore !Processor to which mype is assigned
        !MPI communicator to be passed to SMS
     !   INTEGER (KIND=int_prec) :: MPI_COMM_IPE        

!SMS$IGNORE BEGIN
        OPEN(LUN_nmlt, FILE=INPTNMLT,ERR=222,IOSTAT=IOST_OP,STATUS='OLD')
        OPEN(LUN_nmlt2,FILE=INPTNMLT2,ERR=222,IOSTAT=IOST_OP,STATUS='OLD')
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=IPEDIMS,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMIPE ,ERR=222,IOSTAT=IOST_RD)
!SMS$IGNORE END


!SMS$INSERT lpHaloSize=1
!SMS$INSERT mpHaloSize=2

!nm20170906: i must find a way to consolidate these two options to get it to work for wam-ipe.
!>>>>1.IPEOptimization version:
!SMS$INSERT call MPI_COMM_RANK(MPI_COMM_WORLD,mype,istat)
!SMS$ignore begin
!SMS$INSERT print*,'mype,mod(mype,48)',mype,mod(mype,48)
!SMS$ignore end
!!SMS$INSERT call set_affinity (mod(mype,48)) !Pin MPI rank mype to core mod(mype,48)

!!!>>>>2.WAM-IPE version:
!!!set up MPI communicator for SMS
!(1) when NEMS is used, MPI_COMM_IPE=mpiCommunicator has been assigned already in sub-myIPE_Init
!(2) when IPE Standalone is used, MPI_COMM_IPE=MPI_COMM_WORLD is set in driver_ipe
!SMS$SET_COMMUNICATOR( MPI_COMM_IPE )


!SMS$CREATE_DECOMP(dh,<NLP,NMP>,<lpHaloSize,mpHaloSize>: <NONPERIODIC, PERIODIC>)

!SMS$SERIAL(<IOST_RD,istat,OUT>) BEGIN
        IOST_RD = 0
        istat   = 0
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMFLIP,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMSWITCH,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMMSIS,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt2
        READ(LUN_nmlt2,NML=NMWEIM,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt2
        READ(LUN_nmlt2,NML=NMTIROS,ERR=222,IOSTAT=IOST_RD)
        IF (sw_neutral == 1) THEN
          REWIND LUN_nmlt
          READ(LUN_nmlt,NML=IPECAP   ,ERR=222,IOSTAT=IOST_RD)
        END IF

        OPEN(UNIT=LUN_LOG0,FILE=filename,STATUS='unknown',FORM='formatted',IOSTAT=istat)
        IF ( istat /= 0 ) THEN
          WRITE( UNIT=6, FMT=*)'ERROR OPENING FILE',filename
          STOP
        END IF
        WRITE(UNIT=LUN_LOG0, NML=NMIPE)
        WRITE(UNIT=LUN_LOG0, NML=NMFLIP)
        WRITE(UNIT=LUN_LOG0, NML=NMSWITCH)
        WRITE(UNIT=LUN_LOG0, NML=NMMSIS)
        WRITE(UNIT=LUN_LOG0, NML=NMWEIM)
        WRITE(UNIT=LUN_LOG0, NML=NMTIROS)
        WRITE(UNIT=LUN_LOG0,FMT=*)'NMP=',NMP,' NLP=',NLP,' NPTS2D=',NPTS2D
        WRITE(UNIT=LUN_LOG0,FMT=*)'real_prec=',real_prec,' int_prec=',int_prec
        CLOSE(LUN_LOG0)
!SMS$SERIAL END
        CLOSE(LUN_nmlt)
        CLOSE(LUN_nmlt2)
222     IF ( IOST_OP /= 0 ) THEN
          WRITE(UNIT=LUN_nmlt, FMT=*) "OPEN NAMELIST FAILED!", IOST_OP
          STOP
        ELSEIF ( IOST_RD /= 0 ) THEN
          WRITE(UNIT=LUN_nmlt, FMT=*) "READ NAMELIST FAILED!", IOST_RD
          STOP
        ENDIF


stop_time=start_time+duration

WRITE(*,*)" DATE: 08 September, 2011"
WRITE(*,*)"********************************************"
WRITE(*,*)"***      Copyright 2011 NAOMI MARUYAMA   ***"
WRITE(*,*)"***      ALL RIGHTS RESERVED             ***"
WRITE(*,*)"********************************************"
WRITE(*,*)" LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model"
WRITE(*,*)" DEVELOPER: Dr. Naomi Maruyama"
WRITE(*,*)" CONTACT INFORMATION:"
WRITE(*,*)" E-MAIL : Naomi.Maruyama@noaa.gov"
WRITE(*,*)" PHONE  : 303-497-4857"
WRITE(*,*)" ADDRESS: 325 Broadway, Boulder, CO 80305"
WRITE(*,*)"                                            "

!SMS$insert call NNT_NPROCS(nprocs)
!SMS$insert call NNT_ME    (mype  )
!SMS$TO_LOCAL(dh:<1,lps:lbound>,<1,NLP:ubound>) BEGIN
lps = 1
lpe = NLP
!SMS$TO_LOCAL END
!SMS$TO_LOCAL(dh:<2,mps:lbound>,<2,NMP:ubound>) BEGIN
mps = 1
mpe = NMP
!SMS$TO_LOCAL END
print *,'finished reading namelist:',filename
print *,' '
print"(' NLP:                 ',I6)",NLP
print"(' NMP:                 ',I6)",NMP
print"(' mpstop:              ',I6)",mpstop
print"(' stop_time            ',I6)",stop_time
print"(' Number of Processors:',I6)",nprocs
print"(' lpHaloSize:          ',I6)",lpHaloSize
print"(' mpHaloSize:          ',I6)",mpHaloSize
print *,' '
print *,' '
!mycore = mod(mype,NMP/2)
!!SMS$ignore begin
!print*,'mype,mycore=',mype,mycore
!!SMS$ignore end
if(mype == 39) then
!!SMS$INSERT call   set_affinity (24) !Pin MPI rank mype to core mycore
endif
!!SMS$INSERT call   set_affinity (mype) !Pin MPI rank mype to core mype
!!SMS$INSERT call print_affinity (mype)

!dbg20120509        IF ( sw_rw_sw_perp_trans )  CALL setup_sw_perp_transport ()
!note:20120207: v36: used only activating the perp.transport gradually...

if(parallelBuild)then

!dbg20160408 broadcast solar wind parameters to other proccessors
   nElements = size(swbt)
!dbg20160408 sms PPP_BCAST  debug: comment out these lines
!!SMS$INSERT   call MPI_BCAST(swbt   ,nElements,MPI_REAL,0, SMScomm, ierr)
!!SMS$INSERT   call MPI_BARRIER(SMScomm,ierr)
   !
!!SMS$INSERT   call MPI_BCAST(swangle,nElements,MPI_REAL,0, SMScomm, ierr)
!!SMS$INSERT   call MPI_BARRIER(SMScomm,ierr)
   !
!!SMS$INSERT   call MPI_BCAST(swvel  ,nElements,MPI_REAL,0, SMScomm, ierr)
!!SMS$INSERT   call MPI_BARRIER(SMScomm,ierr)
   !
!!SMS$INSERT   call MPI_BCAST(gwatts ,nElements,MPI_REAL,0, SMScomm, ierr)
!!SMS$INSERT   call MPI_BARRIER(SMScomm,ierr)
   !
!!SMS$INSERT   call MPI_BCAST(levpi  ,nElements,MPI_INTEGER,0, SMScomm, ierr)
!!SMS$INSERT   call MPI_BARRIER(SMScomm,ierr)
   !
endif !parallelB
!dbg

!dbg20160711
!SMS$IGNORE begin
print*,mype,'sub-read_input: swNeuPar',swNeuPar
!SMS$IGNORE end
 
        END SUBROUTINE read_input_parameters

END MODULE module_input_parameters
