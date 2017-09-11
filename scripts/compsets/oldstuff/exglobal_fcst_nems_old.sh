#!/bin/ksh
################################################################################
####  UNIX Script Documentation Block
#                      .                                             .
# Script name:         exglobal_fcst.sh.ecf
# Script description:  Runs a global spectral model forecast
#
# Author:        Mark Iredell       Org: NP23         Date: 1999-05-01
#
# Abstract: This script runs a single or an ensemble of global spectral model
#           forecasts. The initial conditions and run parameters are either
#           passed in the argument list or imported.
#
# Script history log:
# 1999-05-01  Mark Iredell
# 2005-01-03  Sarah Lu : add namelist SOIL_VEG; set FSMCL(2:4)=FSMCL2; 
#                      : add FNVMNC,FNVMXC,FNSLPC,FNABSC
# 2006-2012   Shrinivas Moorthi 
#                      : Modified to run ESMF - Stand Alone version of ATM - Only filestyle "L" allowed 
#                      : Added a ESMF config file. The script can run up to 21 ENS members concurrently.
#                      : Added default PE$n values to 0
#                      : Added stochastic ensembles, semi-Lagrangian high frequency output, 
#                      : G3D outputs and many other upgrades related to model changes
#                      : Upgraded for the new physics and nst model
#                      : and rewrote some parts and added ensemble generality
# 2010-01     Weiyu Yang: modified for the ensemble GEFS.
# 2009-2012   Sarah Lu : ESMF_State_Namelist modified
#                      : Added GOCART_CLIM and GOCART_LUTS tracer added; (q, oz, cld) removed
#                      : modify phy_namelist.rc (set p_export/dp_export to 1)
#                      : add passive_tracer to atm_namelist.rc 
#                      : use wildcard to copy files from GOCART_CLIM to DATA
#                      : add AER, modify filename_base, file_io_form, file_io
#                      : add thermodyn_id and sfcpress_id to nam_dyn
#                      : add WRITE_DOPOST, GOCART_POSTOUT, POST_GRIBVERSION
#                      : change GOCART_POSTOUT to GOCART_AER2POST
#                      : modify how LATG is specified
# 2009-2012  Jun Wang  : add write grid component option
#                      : activate reduced grid option and digital filter option
#                      : link atm_namelist.rc to configure_file
#                      : Add restart
#                      : add option to output filtered 3hr output
#                      : add copy MAPL/CHEM config files,set default FILE_IO_FORM
#                      : add grib2 option for POST
#                      : add POST_NCEPGRB2TBL for post
#                      : set dp_import to 1 for NDSL
#                      : set sigio option
# 2010-2012 Henry Juang: add ndslfv, process_split, and mass_dp for NDSL
#                      : add JCAPG for NDSL, add option of IDEA, remove JCAPG
# 2013-07    Sarah Lu  : specify wgrib and nemsioget for multi-platform
# 2013-11  Xingren Wu  : add A2OI for Atm/Ocn/Ice coupling
# 2014-04  Xingren Wu  : add CPLFLX for Atm/Ocn/Ice coupling
# 2014-12  Kate Howard : Rework to run NEMS with GFS scripts
# 2014-2015 S. Moorthi : Clean up, fix slg option, gaea specific  etc
#                      : Clean up, unify gloabl and ngac scripts remove scheduler etc
#                      : Added THEIA option ; turned off ESMF Compliance check etc
#                      : added MICRO_PHY_DATA
# 2015-10 Fanglin Yang : debug and update to be able to run both fcst1 and fcst2 for NEMS GFS
#
# Usage:  exglobal_fcst.sh.ecf SIGI/GRDI SFCI SIGO FLXO FHOUT FHMAX IGEN D3DO NSTI NSTO FHOUT_HF FHMAX_HF
#
#   Input script positional parameters:
#     1             Input grd file 1
#                   defaults to $SIGI or $GRDI; one or the other is required
#     2             Input surface file
#                   defaults to $SFCI; one or the other is required
#     3             Output sigma file with embedded forecast hour '${FH}'
#                   defaults to $SIGO, then to ${COMOUT}/${SIGOSUF}f'${FH}'$SUFOUT
#     4             Output flux file with embedded forecast hour '${FH}'
#                   defaults to $FLXO, then to ${COMOUT}/${FLXOSUF}f'${FH}'$SUFOUT
#     5             Output frequency in hours
#                   defaults to $FHOUT, then to 3
#     6             Length of forecast in hours
#                   defaults to $FHMAX; otherwise FHSEG is required to be set
#     7             Output generating code
#                   defaults to $IGEN, defaults to 0
#     8             Output flux file with embedded forecast hour '${FH}'
#                   defaults to $D3DO, then to ${COMOUT}/d3df'${FH}'$SUFOUT
#     9             Input NST file
#     10            Output NST file
#     11            High frequency output interval ; default 1 hour
#     12            Maximum hour of high frequecny output
#
#   Imported Shell Variables:
#     SIGI/GRDI     Input sigma file
#                   overridden by $1; one or the other is required
#     SFCI          Input surface file
#                   overridden by $2; one or the other is required
#     SIGO          Output sigma file with embedded forecast hour '${FH}'
#                   overridden by $3; defaults to ${COMOUT}/${SIGOSUF}f'${FH}'$SUFOUT
#     FLXO          Output flux file with embedded forecast hour '${FH}'
#                   overridden by $4; defaults to ${COMOUT}/${FLXOSUF}f'${FH}'$SUFOUT
#     D3DO          Output d3d file with embedded forecast hour '${FH}'
#                   overridden by $4; defaults to ${COMOUT}/d3df'${FH}'$SUFOUT
#     NSTO          Output nst file with embedded forecast hour '${FH}'
#                   overridden by $4; defaults to ${COMOUT}/nstf'${FH}'$SUFOUT
#     FHOUT         Output frequency in hours
#                   overridden by $5; defaults to 3
#     FHMAX         Length of forecast in hours
#                   overridden by $6; either FHMAX or FHSEG must be set
#     IGEN          Output generating code
#                   overridden by $7; defaults to 0
#     FIXGLOBAL     Directory for global fixed files
#                   defaults to /nwprod/gsm.v12.0.0/fix/fix_am
#     EXECGLOBAL    Directory for global executables
#                   defaults to /nwprod/gsm.v12.0.0/exec
#     DATA          working directory
#                   (if nonexistent will be made, used and deleted)
#                   defaults to current working directory
#     COMOUT        output directory
#                   (if nonexistent will be made)
#                   defaults to current working directory
#     XC            Suffix to add to executables
#                   defaults to none
#     SUFOUT        Suffix to add to output filenames
#                   defaults to none
#     NCP           Copy command
#                   defaults to cp
#     SIGHDR        Command to read sigma header
#                   (required if JCAP, LEVS, or FHINI are not specified)
#                   defaults to ${EXECGLOBAL}/global_sighdr$XC
#     JCAP          Spectral truncation for model wave
#     LEVS          Number of levels
#                   defaults to the value in the input sigma file header
#     LEVR          Number of levels over which radiation is computed
#                   defaults to LEVS
#     FCSTEXEC      Forecast executable
#                   defaults to ${EXECGLOBAL}/global_fcst$XC
#     SIGI2         Second time level sigma restart file
#                   defaults to NULL
#     CO2CON        Input CO2 radiation (vertical resolution dependent)
#                   defaults to ${FIXGLOBAL}/global_co2con.l${LEVS}.f77
#     MTNVAR        Input mountain variance (horizontal resolution dependent)
#                   defaults to ${FIXGLOBAL}/global_mtnvar.t${JCAP}.f77
#     MTNRSL        A string representing topography resolution
#                   defaults to $JCAP
#     MTNRSLUF      A string representing unfiltered topography resolution
#                   defaults to $MTNRSL
#     O3FORC        Input ozone forcing (production/loss) climatology
#                   defaults to ${FIXGLOBAL}/global_o3prdlos.f77
#     O3CLIM        Input ozone climatology
#                   defaults to ${FIXGLOBAL}/global_o3clim.txt
#     FNGLAC        Input glacier climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_glacier.2x2.grb
#     FNMXIC        Input maximum sea ice climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_maxice.2x2.grb
#     FNTSFC        Input SST climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_sstclim.2x2.grb
#     FNSNOC        Input snow climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_snoclim.1.875.grb
#     FNZORC        Input roughness climatology GRIB file
#                   defaults to 'sib' (From sib vegetation-based lookup table.
#                   FNVETC must be set to sib file: ${FIXGLOBAL}/global_vegtype.1x1.grb)
#     FNALBC        Input albedo climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_albedo4.1x1.grb
#     FNAISC        Input sea ice climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_iceclim.2x2.grb
#     FNTG3C        Input deep soil temperature climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_tg3clim.2.6x1.5.grb
#     FNVEGC        Input vegetation fraction climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_vegfrac.1x1.grb
#     FNVETC        Input vegetation type climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_vegtype.1x1.grb
#     FNSOTC        Input soil type climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_soiltype.1x1.grb
#     FNSMCC        Input soil moisture climatology GRIB file
#                   defaults to ${FIXGLOBAL}/global_soilmcpc.1x1.grb
#     FNVMNC        Input min veg frac climatology GRIB file    
#                   defaults to ${FIXGLOBAL}/global_shdmin.0.144x0.144.grb
#     FNVMXC        Input max veg frac climatology GRIB file    
#                   defaults to ${FIXGLOBAL}/global_shdmax.0.144x0.144.grb
#     FNSLPC        Input slope type climatology GRIB file    
#                   defaults to ${FIXGLOBAL}/global_slope.1x1.grb
#     FNABSC        Input max snow albedo climatology GRIB file    
#                   defaults to ${FIXGLOBAL}/global_snoalb.1x1.grb
#     OROGRAPHY     Input orography GRIB file (horiz resolution dependent)
#                   defaults to ${FIXGLOBAL}/global_orography.t$JCAP.grb
#     OROGRAPHY_UF  Input unfiltered orography GRIB file (resolution dependent)
#                   defaults to ${FIXGLOBAL}/global_orography_uf.t$JCAP.grb
#     LONSPERLAT    Input txt file containing reduced grid information
#                   defaults to ${FIXGLOBAL}/global_lonsperlat.t$MTNRSL.txt}
#     FNMSKH        Input high resolution land mask GRIB file
#                   defaults to ${FIXGLOBAL}/seaice_newland.grb
#     FNTSFA        Input SST analysis GRIB file
#                   defaults to none
#     FNACNA        Input sea ice analysis GRIB file
#                   defaults to none
#     FNSNOA        Input snow analysis GRIB file
#                   defaults to none
#     AERODIR       Input aersol climatology directory
#                   defaults to ${FIXGLOBAL}
##########################################################################
#     FIX_RAD       Directory for global fixed files
#                   Defaults to $${FIXGLOBAL}
#     EMISDIR       Input earth's surface emissivity data directory
#                   defaults to ${FIX_RAD} - export IEMS=1 to activate
#     SOLCDIR       11 year cycle Solar constat data directory
#                   defaults to ${FIX_RAD} - export ISOL=1,2,3,4, or 10 to activate
#     VOLCDIR       Volcanic aerosol  data directory
#                   defaults to ${FIX_RAD} - export IAER=100,101, or 110 to activate
#     CO2DIR        Historical CO2 data directory
#                   defaults to ${FIX_RAD} - export ICO2=1 or 2 to activate
#                   ICO2=1 gives annual mean and ICO2=2 uses monthly 2D data
##########################################################################
#     GOCART_CLIM   Directory for gocart climo files
#                   Defaults to $${FIXGLOBAL}
#     GOCARTC_LUTS  Directory for gocart luts files
#                   Defaults to $${FIXGLOBAL}
#     SIGR1         Output first time level sigma restart file
#                   defaults to ${DATA}/sigr1 which is deleted
#     SIGR2         Output second time level sigma restart file
#                   defaults to ${DATA}/sigr2 which is deleted
#     SFCR          Output surface restart file
#                   defaults to ${DATA}/sfcr which is deleted
#     NSTR          Output nst restart file
#                   defaults to ${DATA}/nstr which is deleted
#     SFCO          Output surface file with embedded forecast hour '${FH}'
#                   defaults to ${COMOUT}/${SFCOSUF}f'${FH}'$SUFOUT
#     LOGO          Output log file with embedded forecast hour '${FH}'
#                   defaults to ${COMOUT}/logf'${FH}'$SUFOUT
#     INISCRIPT     Preprocessing script
#                   defaults to none
#     LOGSCRIPT     Log posting script
#                   defaults to none
#     ERRSCRIPT     Error processing script
#                   defaults to 'eval [[ $err = 0 ]]'
#     ENDSCRIPT     Postprocessing script
#                   defaults to none
#     FHINI         Starting forecast hour
#                   defaults to the value in the input sigma file header
#     FHSEG         Number of hours to integrate
#                   (only required if FHMAX is not specified)
#                   defaults to 0
#     DELTIM        Timestep in seconds
#                   defaults to 3600/($JCAP/20)
#     FHRES         Restart frequency in hours
#                   defaults to 24
#     FHZER         Zeroing frequency in hours
#                   defaults to 6
#     FHLWR         Longwave radiation frequency in seconds
#                   defaults to 3600
#     FHSWR         Shortwave radiation frequency in seconds
#                   defaults to 3600
#     FHROT         Forecast hour to Read One Time level
#                   defaults to 0
#     FHDFI         Half number of hours of digital filter initialization
#                   defaults to 0
#     FHCYC         Surface cycling frequency in hours
#                   defaults to 0 for no cycling
#     IDVC          Integer ID of the vertical coordinate type
#                   defaults to that in the header for the input upperair
#                   file. IDVC=1 for sigma; IDVC=2 for pressure/sigma hybrid
#     TFILTC        Time filter coefficient
#                   defaults to 0.85
#     DYNVARS       Other namelist inputs to the dynamics executable
#                   defaults to none set
#     PHYVARS       Other namelist inputs to the physics executable
#                   defaults to none set
#     TRACERVARS    Other namelist inputs to the forecast executable
#                   defaults to none set
#     FSMCL2        Scale in days to relax to soil moisture climatology
#                   defaults to 99999 for no relaxation
#     FTSFS         Scale in days to relax to SST anomaly to zero
#                   defaults to 90
#     FAISS         Scale in days to relax to sea ice to climatology
#                   defaults to 99999
#     FSNOL         Scale in days to relax to snow to climatology
#                   defaults to 99999
#     FSICL         Scale in days to relax to sea ice to climatology
#                   defaults to 99999
#     FZORL         Scale in days to relax to roughness climatology.
#                   defaults to 99999 because the 'sib' option sets
#                   roughness from a lookup table and is static.
#     CYCLEVARS     Other namelist inputs to the surface cycling
#                   defaults to none set
#     NTHREADS      Number of threads
#                   defaults to 1
#     SPECTRAL_LOOP Number of spectral loops
#                   defaults to 2
#     NTHSTACK      Size of stack per thread
#                   defaults to 64000000
#     FILESTYLE     File management style flag
#                   ('L' for symbolic links in $DATA is the only allowed style),
#     PGMOUT        Executable standard output
#                   defaults to $pgmout, then to '&1'
#     PGMERR        Executable standard error
#                   defaults to $pgmerr, then to '&1'
#     pgmout        Executable standard output default
#     pgmerr        Executable standard error default
#     REDOUT        standard output redirect ('1>' or '1>>')
#                   defaults to '1>', or to '1>>' to append if $PGMOUT is a file
#     REDERR        standard error redirect ('2>' or '2>>')
#                   defaults to '2>', or to '2>>' to append if $PGMERR is a file
#     VERBOSE       Verbose flag (YES or NO)
#                   defaults to NO
#
#   Exported Shell Variables:
#     PGM           Current program name
#     pgm
#     ERR           Last return code
#     err
#
#   Modules and files referenced:
#     scripts    : $INISCRIPT
#                  $LOGSCRIPT
#                  $ERRSCRIPT
#                  $ENDSCRIPT
#
#     programs   : $FCSTEXEC
#
#     input data : $1 or $SIGI
#                  $2 or $SFCI
#                  $SIGI2
#                  $FNTSFA
#                  $FNACNA
#                  $FNSNOA
#
#     fixed data : $CO2CON
#                  $MTNVAR
#                  $O3FORC
#                  $O3CLIM
#                  $FNGLAC
#                  $FNMXIC
#                  $FNTSFC
#                  $FNSNOC
#                  $FNZORC
#                  $FNALBC
#                  $FNAISC
#                  $FNTG3C
#                  $FNVEGC
#                  $FNVETC
#                  $FNSOTC
#                  $FNSMCC
#                  $FNVMNC
#                  $FNVMXC
#                  $FNSLPC
#                  $FNABSC
#                  $FNMSKH
#                  $OROGRAPHY
#                  $OROGRAPHY_UF
#                  $LONSPERLAT
#
#     output data: $3 or $SIGO
#                  $4 or $FLXO
#                  $SFCO
#                  $LOGO
#                  $SIGR1
#                  $SIGR2
#                  $SFCR
#                  $NSTR
#                  $PGMOUT
#                  $PGMERR
#
#     scratch    : ${DATA}/fort.11
#                  ${DATA}/fort.12
#                  ${DATA}/fort.14
#                  ${DATA}/fort.15
#                  ${DATA}/fort.24
#                  ${DATA}/fort.28
#                  ${DATA}/fort.29
#                  ${DATA}/fort.48
#                  ${DATA}/fort.51
#                  ${DATA}/fort.52
#                  ${DATA}/fort.53
#                  ${DATA}/SIG.F*
#                  ${DATA}/SFC.F*
#                  ${DATA}/FLX.F*
#                  ${DATA}/LOG.F*
#                  ${DATA}/D3D.F*
#                  ${DATA}/G3D.F*
#                  ${DATA}/NST.F*
#                  ${DATA}/sigr1
#                  ${DATA}/sigr2
#                  ${DATA}/sfcr
#                  ${DATA}/nstr
#                  ${DATA}/NULL
#
# Remarks:
#
#   Condition codes
#      0 - no problem encountered
#     >0 - some problem encountered
#
#  Control variable resolution priority
#    1 Command line argument.
#    2 Environment variable.
#    3 Inline default.
#
# Attributes:
#   Language: POSIX shell
#   Machine: WCOSS, GAEA, THEIA
#
####
################################################################################
#  Set environment.
export VERBOSE=${VERBOSE:-"NO"}
if [[ $VERBOSE = YES ]] ; then
  echo $(date) EXECUTING $0 $* >&2
  set -x
fi
export COMPLIANCECHECK=${COMPLIANCECHECK:-OFF}
export ESMF_RUNTIME_COMPLIANCECHECK=$COMPLIANCECHECK:depth=4
#
export machine=${machine:-WCOSS}
export machine=$(echo $machine|tr '[a-z]' '[A-Z]')
export FCST_LAUNCHER=${FCST_LAUNCHER:-${APRUN:-""}}
if [ $machine = THEIA ]; then
  export MPICH_FAST_MEMCPY=${MPICH_FAST_MEMCPY:-"ENABLE"}
  export MPI_BUFS_PER_PROC=${MPI_BUFS_PER_PROC:-2048}
  export MPI_BUFS_PER_HOST=${MPI_BUFS_PER_HOST:-2048}
  export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}


#. /apps/lmod/5.8/init/ksh
# module load intel/14.0.2
# module load  impi/4.1.3.048
# module load intel/15.1.133 impi/5.0.3.048 netcdf/4.3.0 szip/2.1

. /apps/lmod/lmod/init/sh
#source $BASE_GSM/sorc/global_fcst.fd/NEMS/src/conf/modules.nems
 source $BASE_NEMS/src/conf/modules.nems


elif [ $machine = GAEA ]; then
  export MPICH_FAST_MEMCPY=${MPICH_FAST_MEMCPY:-"ENABLE"}
  export MPICH_MAX_SHORT_MSG_SIZE=${MPICH_MAX_SHORT_MSG_SIZE:-4096}
  export MPICH_UNEX_BUFFER_SIZE=${MPICH_UNEX_BUFFER_SIZE:-1024000000}
  export MPICH_PTL_UNEX_EVENTS=${MPICH_PTL_UNEX_EVENTS:-400000}
  export MPICH_PTL_OTHER_EVENTS=${MPICH_PTL_OTHER_EVENTS:-100000}
  export MPMD_PROC=${MPMD_PROC:-NO}
  export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
elif [ $machine = WCOSS ] ; then
  if [ ${LOADICS:-YES} = YES ] ; then
#   . /usrx/local/Modules/3.2.9/init/ksh
#   module unload ics
#   export ICS_VERSION=${ICS_VERSION:-15.0.1}
#   module load ics/$ICS_VERSION


#.  /usrx/local/Modules/default/init/bash
.   /usrx/local/Modules/default/init/sh
    source $BASE_NEMS/src/conf/modules.nems

    export MKL_CBWR=${MKL_CBWR:-AVX}          # Needed for bit reproducibility with mkl
    export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
    export SAVE_ALL_TASKS=${SAVE_ALL_TASKS:-yes}
    export PROFILE_BY_CALL_SITE=${PROFILE_BY_CALL_SITE:-yes}

    if [ ${USEBULKXFER:-NO} = YES ] ; then
      module unload ibmpe
      module load ibmpe/1.3.0.8p
      export MP_USE_BULK_XFER=yes
      export MP_BULK_MIN_MSG_SIZE=512K
      export MP_RC_USE_LMC=yes
    fi
  fi
# export MP_EAGER_LIMIT=${MP_EAGER_LIMIT:-64K}
# export FORT_BUFFERED=${FORT_BUFFERED:-true}
  export MP_EUIDEVICE=${MP_EUIDEVICE:-min}

  export MP_EUILIB=${MP_EUILIB:-us}

# export MP_TASK_AFFINITY=${MP_TASK_AFFINITY:-"cpu:$NTHREADS"}


#test
  export MPICH_ALLTOALL_THROTTLE=${MPICH_ALLTOALL_THROTTLE:-0}
  export MP_SINGLE_THREAD=${MP_SINGLE_THREAD:-yes}

#set -x
#ulimit -s unlimited
#ulimit -c unlimited

#export OMP_NUM_THREADS=1
#export MP_EUILIB=us
#export MP_MPILIB=mpich2


   
  export VPROF_PROFILE=${VPROF_PROFILE:-no}
  export MP_COREFILE_FORMAT=${MP_COREFILE_FORMAT:-"lite"}
# export MP_USE_TOKEN_FLOW_ CONTROL=${MP_USE_TOKEN_FLOW_ CONTROL:-yes}
# export MP_S_ENABLE_ERR_PRINT=yes
fi

export APRUN=${APRUN:-""}
export model=${model:-global}
export NEMSIO_IN=${NEMSIO_IN:-".true."}
export NEMSIO_OUT=${NEMSIO_OUT:-".true."}
export ENS_NUM=${ENS_NUM:-1}
export FM=${FM}

#  Command line arguments.
if [ $NEMSIO_IN = .false. ] ; then
 export SIGI=${1:-${SIGI:-?}}
else
 export GRDI=${1:-${GRDI:-?}}
 export SIGI=${1:-${GRDI:-?}}
fi

export SFCI=${2:-${SFCI:-?}}
export SIGO=${3:-${SIGO}}
export FLXO=${4:-${FLXO}}
export FHOUT=${5:-${FHOUT:-3}}
export FHMAX=${6:-${FHMAX:-0}}
export IGEN=${7:-${IGEN:-0}}
export D3DO=${8:-${D3DO}}
export NSTI=${9:-${NSTI:-?}}
export NSTO=${10:-${NSTO}}
export FHOUT_HF=${11:-${FHOUT_HF:-0}}
export FHMAX_HF=${12:-${FHMAX_HF:-0}}
export AERO=${13:-${AERO}}

# DHOU 02/28/2008 Modified for general case
# DHOU 01/07/2008 Added two input for the GEFS_Cpl module
# FHM_FST is the FHMAX for the integration before the first stop
# FH_INC is the FHMAX_increase for the integration before next stop
export FH_INC=${FH_INC:-100000000}
export ENS_SPS=${ENS_SPS:-.false.}
export ADVANCECOUNT_SETUP=${ADVANCECOUNT_SETUP:-0}
export HOUTASPS=${HOUTASPS:-10000}

export SPS_PARM1=${SPS_PARM1:-"0.005 10.0 0.005 10.0 0.0 0.0 0.0 0.0 0.0 0.0"}
export SPS_PARM2=${SPS_PARM2:-"0.105 0.03 0.12 42.0 0.0 0.0 0.0 0.0 0.0 0.0"}
export SPS_PARM3=${SPS_PARM3:-"0.2 0.34 -0.34 3.0 0.0 0.0 0.0 0.0 0.0 0.0"}

[[ $ENS_NUM -lt 2 ]]&&ENS_SPS=.false.
if [ $ENS_SPS = .false. ] ; then export FH_INC=$FHMAX ; fi

#  Directories.
export gsm_ver=${gsm_ver:-v14.0.0}
export BASEDIR=${BASEDIR:-/nwprod}
export NWPROD=${NWPROD:-$BASEDIR}
export FIXSUBDA=${FIXSUBDA:-fix/fix_am}
export FIXGLOBAL=${FIXGLOBAL:-$NWPROD/gsm.$gsm_ver/$FIXSUBDA}
export FIX_RAD=${FIX_RAD:-$FIXGLOBAL}
export FIX_IDEA=${FIX_IDEA:-$FIXGLOBAL}
export FIX_NGAC=${FIX_NGAC:-$NWPROD/fix/fix_ngac}
export PARMSUBDA=${PARMSUBDA:-parm/parm_am}
export PARMGLOBAL=${PARMGLOBAL:-$NWPROD/gsm.$gsm_ver/$PARMSUBDA}
export PARM_NGAC=${PARM_NGAC:-$NWPROD/parm/parm_ngac}
export EXECGLOBAL=${EXECGLOBAL:-$NWPROD/exec}
export DATA=${DATA:-$(pwd)}
export COMOUT=${COMOUT:-$(pwd)}

#  Filenames.
MN=${MN:-""}
export XC=${XC}
export SUFOUT=${SUFOUT}
export NCP=${NCP:-"/bin/cp -p"}

export nemsioget=${nemsioget:-/nwprod/ngac.v1.0.0/exec/nemsio_get}
if [ $NEMSIO_IN = .true. ]; then
 export JCAP=${JCAP:-$($nemsioget ${GRDI}$FM jcap |grep -i "jcap" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export LEVS=${LEVS:-$($nemsioget ${GRDI}$FM levs|grep -i "levs" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export LEVR=${LEVR:-$LEVS}
 export LONF=${LONF:-$($nemsioget ${GRDI}$FM LONF|grep -i "lonf" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 if [[ ${RESTART:-".false."} = .true. ]] ; then
  export LATG=${LATG:-$($nemsioget ${GRDI}$FM LATF|grep -i "latf" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 else
  export LATG=${LATG:-$($nemsioget ${GRDI}$FM LATG|grep -i "latg" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 fi
 export LONR=${LONR:-$LONF}
 export LATR=${LATR:-$LATG}
 export NTRAC=${NTRAC:-$($nemsioget ${GRDI}$FM NTRAC|grep -i "NTRAC" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export IDVC=${IDVC:-$($nemsioget ${GRDI}$FM IDVC |grep -i "IDVC" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export IDVM=${IDVM:-$($nemsioget ${GRDI}$FM IDVM |grep -i "IDVM" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export FHINI=${FHINI:-$($nemsioget ${GRDI}$FM NFHOUR |grep -i "NFHOUR" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
else


### AMK for sfca03 IAU ###
 export SFCHDR=${SFCHDR:-${EXECGLOBAL}/global_sfchdr$XC}
 export CHGSFCFHREXEC=${CHGSFCFHREXEC:-/swpc/save/swpc.spacepara/util/chgsfcfhr/chgsfcfhr}
 ##########################


 export SIGHDR=${SIGHDR:-${EXECGLOBAL}/global_sighdr$XC}
 export JCAP=${JCAP:-$(echo jcap|$SIGHDR ${SIGI}$FM)}
 export LEVS=${LEVS:-$(echo levs|$SIGHDR ${SIGI}$FM)}
 export LEVR=${LEVR:-$LEVS}
 export LONR=${LONR:-$(echo lonr|$SIGHDR ${SIGI}$FM)}
 export LATR=${LATR:-$(echo latr|$SIGHDR ${SIGI}$FM)}
 export LONF=${LONF:-$(echo lonf|$SIGHDR ${SIGI}$FM)}
 export LATG=${LATG:-$(echo latf|$SIGHDR ${SIGI}$FM)}
 export NTRAC=${NTRAC:-$(echo ntrac|$SIGHDR ${SIGI}$FM)}
 export IDVC=${IDVC:-$(echo idvc|$SIGHDR ${SIGI}$FM)}
 export IDVM=${IDVM:-$(echo idvm|$SIGHDR ${SIGI}$FM)}
 export FHINI=${FHINI:-$(echo ifhr|$SIGHDR ${SIGI}$FM)}
fi

export LONB=${LONB:-$LONF}
export LATB=${LATB:-$LATG}
export THERMODYN_ID=${THERMODYN_ID:-$((IDVM/10))}
export SFCPRESS_ID=${SFCPRESS_ID:-$((IDVM-(IDVM/10)*10))}
export NMTVR=${NMTVR:-14}
export LSOIL=${LSOIL:-4}
export NTOZ=${NTOZ:-2}
export NTCW=${NTCW:-3}
export NCLD=${NCLD:-1}
export NGPTC=${NGPTC:-30}

export ADIAB=${ADIAB:-.false.}
export nsout=${nsout:-0}
export LDFIFLTO=${LDFIFLTO:-.false.}
export LDFI_GRD=${LDFI_GRD:-.false.}
export DFILEVS=${DFILEVS:-$LEVS}
export liope=${liope:-.false.}
export NUM_FILE=${NUM_FILE:-3}
export QUILTING=${QUILTING:-.true.}
export REDUCED_GRID=${REDUCED_GRID:-.true.}
export PASSIVE_TRACER=${PASSIVE_TRACER:-.false.}
export NST_FCST=${NST_FCST:-0}
export IAER=${IAER:-0}
export GOCART=${GOCART:-0}
#export NGRID_A2OI=${NGRID_A2OI:-20}
#export A2OI_OUT=${A2OI_OUT:-.false.}
#export CPLFLX=${CPLFLX:-.false.}
export NDSLFV=${NDSLFV:-.false.}
export EXPLICIT=${EXPLICIT:-.false.}
export MASS_DP=${MASS_DP:-.false.}
export PROCESS_SPLIT=${PROCESS_SPLIT:-.false.}
export dp_import=${dp_import:-1}
export p_import=${p_import:-1}
export dpdt_import=${dpdt_import:-0}
if [ $NDSLFV = .true. ] ; then
 export MASS_DP=.true.
 export PROCESS_SPLIT=.false.
 export dp_import=1
fi
export ZFLXTVD=${ZFLXTVD:-.false.}
export SEMI_IMPLICIT_TEMP_PROFILE=${SEMI_IMPLICIT_TEMP_PROFILE:-.false.}
#
export FCSTEXEC=${FCSTEXEC:-${EXECGLOBAL}/${model}_fcst$XC}
export GRDI2=${GRDI2:-NULL}
export SIGI2=${SIGI2:-NULL}
export CO2CON=${CO2CON:-${FIXGLOBAL}/global_co2con.l${LEVS}.f77}
export MTNRSL=${MTNRSL:-$JCAP}
export MTNRSLUF=${MTNRSLUF:-$MTNRSL}
export MTNVAR=${MTNVAR:-${FIXGLOBAL}/global_mtnvar.t$MTNRSL.f77}
export O3FORC=${O3FORC:-${FIXGLOBAL}/global_o3prdlos.f77}
export O3CLIM=${O3CLIM:-${FIXGLOBAL}/global_o3clim.txt}
export FNGLAC=${FNGLAC:-${FIXGLOBAL}/global_glacier.2x2.grb}
export FNMXIC=${FNMXIC:-${FIXGLOBAL}/global_maxice.2x2.grb}
#export FNTSFC=${FNTSFC:-${FIXGLOBAL}/cfs_oi2sst1x1monclim19822001.grb}
export FNTSFC=${FNTSFC:-${FIXGLOBAL}/RTGSST.1982.2012.monthly.clim.grb}
export FNSNOC=${FNSNOC:-${FIXGLOBAL}/global_snoclim.1.875.grb}
#export FNZORC=${FNZORC:-${FIXGLOBAL}/global_zorclim.1x1.grb}
export FNZORC=${FNZORC:-sib}
export FNALBC=${FNALBC:-${FIXGLOBAL}/global_albedo4.1x1.grb}
#export FNAISC=${FNAISC:-${FIXGLOBAL}/cfs_ice1x1monclim19822001.grb}
export FNAISC=${FNAISC:-${FIXGLOBAL}/CFSR.SEAICE.1982.2012.monthly.clim.grb}
export FNTG3C=${FNTG3C:-${FIXGLOBAL}/global_tg3clim.2.6x1.5.grb}
export FNVEGC=${FNVEGC:-${FIXGLOBAL}/global_vegfrac.0.144.decpercent.grb}
export FNVETC=${FNVETC:-${FIXGLOBAL}/global_vegtype.1x1.grb}
export FNSOTC=${FNSOTC:-${FIXGLOBAL}/global_soiltype.1x1.grb}
#export FNSMCC=${FNSMCC:-${FIXGLOBAL}/global_soilmcpc.1x1.grb}
export FNSMCC=${FNSMCC:-${FIXGLOBAL}/global_soilmgldas.t${JCAP}.${LONR}.${LATR}.grb}
export FNVMNC=${FNVMNC:-${FIXGLOBAL}/global_shdmin.0.144x0.144.grb}
export FNVMXC=${FNVMXC:-${FIXGLOBAL}/global_shdmax.0.144x0.144.grb}
export FNSLPC=${FNSLPC:-${FIXGLOBAL}/global_slope.1x1.grb}
export FNABSC=${FNABSC:-${FIXGLOBAL}/global_snoalb.1x1.grb}
export FNMSKH=${FNMSKH:-${FIXGLOBAL}/seaice_newland.grb}
export OROGRAPHY=${OROGRAPHY:-${FIXGLOBAL}/global_orography.t$MTNRSL.grb}
export OROGRAPHY_UF=${OROGRAPHY_UF:-${FIXGLOBAL}/global_orography_uf.t$MTNRSLUF.grb}
export LONSPERLAT=${LONSPERLAT:-${FIXGLOBAL}/global_lonsperlat.t${JCAP_TMP}.$LONB_TMP.$LATB_TMP.txt}
export LONSPERLAR=${LONSPERLAR:-$LONSPERLAT}
export FNTSFA=${FNTSFA}
export FNACNA=${FNACNA}
export FNSNOA=${FNSNOA}
#
export AERODIR=${AERODIR:-${FIX_RAD}}
export EMISDIR=${EMISDIR:-${FIX_RAD}}
export SOLCDIR=${SOLCDIR:-${FIX_RAD}}
export VOLCDIR=${VOLCDIR:-${FIX_RAD}}
export CO2DIR=${CO2DIR:-${FIX_RAD}}
export GOCART_CLIM=${GOCART_CLIM:-${FIX_RAD}}
export GOCART_LUTS=${GOCART_LUTS:-${FIX_RAD}}
#export ALBDIR=${ALBDIR:-${FIX_RAD}}
export IEMS=${IEMS:-0}
export ISOL=${ISOL:-0}
export IAER=${IAER:-0}
export ICO2=${ICO2:-0}
#
LOCD=${LOCD:-""}
export COMENS=$COMOUT'$LOCD'
export GRDR1=${GRDR1:-${COMENS}/grdr1}
export GRDR2=${GRDR2:-${COMENS}/grdr2}
export SIGR1=${SIGR1:-${COMENS}/sigr1}
export SIGR2=${SIGR2:-${COMENS}/sigr2}
export SFCR=${SFCR:-${COMENS}/sfcr}
export NSTR=${NSTR:-${COMENS}/nstr}

export SIGS1=${SIGS1:-${COMENS}/sigs1}
export SIGS2=${SIGS2:-${COMENS}/sigs2}
export SFCS=${SFCS:-${COMENS}/sfcs}
export NSTS=${NSTS:-${COMENS}/nsts}

## History Files
export SIGO=${SIGO:-${COMENS}/${SIGOSUF}f'${FHIAU}''${MN}'$SUFOUT}
export SFCO=${SFCO:-${COMENS}/${SFCOSUF}f'${FHIAU}''${MN}'$SUFOUT}
export FLXO=${FLXO:-${COMENS}/${FLXOSUF}f'${FHIAU}''${MN}'$SUFOUT}
export LOGO=${LOGO:-${COMENS}/logf'${FHIAU}''${MN}'$SUFOUT}
export D3DO=${D3DO:-${COMENS}/d3df'${FHIAU}''${MN}'$SUFOUT}
export NSTO=${NSTO:-${COMENS}/${NSTOSUF}f'${FHIAU}''${MN}'$SUFOUT}
export AERO=${AERO:-${COMOUT}/aerf'${FH}''${MN}'$SUFOUT}

export INISCRIPT=${INISCRIPT}
export ERRSCRIPT=${ERRSCRIPT:-'eval [[ $err = 0 ]]'}
export LOGSCRIPT=${LOGSCRIPT}
export ENDSCRIPT=${ENDSCRIPT}

#  Other variables.
export FHSEG=${FHSEG:-0}
export FHMAX=${FHMAX:-$((10#$FHINI+10#$FHSEG))}
export DELTIM=${DELTIM:-$((3600/(JCAP/20)))}
export DTPHYS=${DTPHYS:-$((DELTIM/2))}                     
export FHRES=${FHRES:-24}
export FHZER=${FHZER:-6}
export FHLWR=${FHLWR:-3600}
export FHSWR=${FHSWR:-3600}
export FHROT=${FHROT:-0}
export FHDFI=${FHDFI:-1}

export FHCYC=${FHCYC:-0}

export nhours_dfini=${nhours_dfini:-$FHDFI}
export GB=${GB:-0}
export gfsio_in=${gfsio_in:-.false.}
if [ $gfsio_in = .true. ] ; then export GB=1 ; fi

#        WAM related namelist variables
#        ------------------------------
export IDEA=${IDEA:-.false.}
export WAM_IPE_COUPLING=${WAM_IPE_COUPLING:-.false.}
export HEIGHT_DEPENDENT_G=${HEIGHT_DEPENDENT_G:-.false.}
export F107_KP_SIZE=${F107_KP_SIZE:-56}
export F107_KP_DATA_SIZE=${F107_KP_DATA_SIZE:-56}
export F107_KP_SKIP_SIZE=${F107_KP_SKIP_SIZE:-0}
export F107_KP_INTERVAL=${F107_KP_INTERVAL:-10800}


## for post
export WRITE_DOPOST=${WRITE_DOPOST:-.false.}
export GOCART_AER2POST=${GOCART_AER2POST:-.false.}
export POST_GRIBVERSION=${POST_GRIBVERSION:-grib1}
export POSTCTLFILE=${POSTCTLFILE:-$PARM_NGAC/ngac_postcntrl.parm}
export POST_PARM=${POST_PARM:-$PARM_NGAC/ngac_postcntrl.xml}
export POST_AVBLFLDSXML=${POST_AVBLFLDSXML:-$PARM_NGAC/ngac_post_avblflds.xml}
export POST_NCEPGRB2TBL=${POST_NCEPGRB2TBL:-$NWPROD/lib/sorc/g2tmpl/params_grib2_tbl_new}

## copy/link post related files
if [[ $WRITE_DOPOST = .true. ]] ; then
 if [[ $POST_GRIBVERSION = grib1 ]] ; then
   #ln -sf ${POSTCTLFILE} fort.14
   ${NCP} ${POSTCTLFILE} fort.14
 elif [[ $POST_GRIBVERSION = grib2 ]] ; then
   ${NCP} ${POST_PARM}        postcntrl.xml
   ${NCP} ${POST_AVBLFLDSXML} post_avblflds.xml
   ${NCP} ${POST_NCEPGRB2TBL} params_grib2_tbl_new
 fi
 ln -sf griddef.out fort.110
 MICRO_PHYS_DATA=${MICRO_PHYS_DATA:-${POST_LUTDAT:-$NWPROD/$PARMSUBDA/nam_micro_lookup.dat}}
 ${NCP} $MICRO_PHYS_DATA ./eta_micro_lookup.dat
fi
#
# Total pe = WRT_GROUP*WRTPE_PER_GROUP + fcst pes
#
export WRT_GROUP=${WRT_GROUP:-1}
export WRTPE_PER_GROUP=${WRTPE_PER_GROUP:-1}
export WRITE_NEMSIOFLAG=${WRITE_NEMSIOFLAG:-.true.}
export QUILTING=${QUILTING:-.true.}
export GOCART_AER2POST=${GOCART_AER2POST:-.false.}
if [ $NEMSIO_IN = .true. ]; then
  export FILE_IO_FORM=${FILE_IO_FORM:-"'bin4' 'bin4' 'bin4'"}
else
  export FILE_IO_FORM=${FILE_IO_FORM:-"'grib' 'grib' 'grib'"}
fi

#
export LWRTGRDCMP=${LWRTGRDCMP:-".true."}
if [ $NEMSIO_OUT = .false. -a $WRITE_DOPOST = .false. ] ; then
  export LWRTGRDCMP=.false.
fi
# number of output files, default =3, for adiab num_file=1
ioform_sig=${ioform_sig:-bin4}
ioform_sfc=${ioform_sfc:-bin4}
ioform_flx=${ioform_flx:-bin4}
if [[ $ADIAB = .true. ]] ; then
  export NUM_FILE=1 ;
  export FILENAME_BASE="'SIG.F'"
  export FILE_IO_FORM="'bin4'"
else
  export FILENAME_BASE="'SIG.F' 'SFC.F' 'FLX.F'"
  export FILE_IO_FORM=${FILE_IO_FORM:-"'bin4' 'bin4' 'bin4'"}
  export NUM_FILE=3

  if [ $NST_FCST -gt 0 ] ; then
    export FILENAME_BASE=${FILENAME_BASE}" 'NST.F'"
    export FILE_IO_FORM=${FILE_IO_FORM}" 'bin4'"
    export NUM_FILE=$((NUM_FILE+1))
    if [ $NST_FCST -eq 1 ]; then
      NST_SPINUP=1
    fi
  fi
  if [ $GOCART == 1 ] ; then
    export FILENAME_BASE=${FILENAME_BASE}" 'AER.F'"
    export FILE_IO_FORM=${FILE_IO_FORM}" 'grib'"
    export NUM_FILE=$((NUM_FILE+1))
  fi
  echo "NUM_FILE=$NUM_FILE,GOCART=$GOCART,NST_FCST=$NST_FCST,FILENAME_BASE=$FILENAME_BASE"
fi
export NST_SPINUP=${NST_SPINUP:-0}


#wanghj
export NST_FCST=${NST_FCST:-0}
export NST_SPINUP=${NST_SPINUP:-0}
export NST_RESERVED=${NST_RESERVED:-0}
export ZSEA1=${ZSEA1:-0}
export ZSEA2=${ZSEA2:-0}

export nstf_name="$NST_FCST,$NST_SPINUP,$NST_RESERVED,$ZSEA1,$ZSEA2"
export NST_ANL=${NST_ANL:-.false.}



#
if [ $IDVC = 1 ] ; then
 export HYBRID=.false.
 export GEN_COORD_HYBRID=.false.
elif [ $IDVC = 2 ] ; then
 export HYBRID=.true.
 export GEN_COORD_HYBRID=.false.
elif [ $IDVC = 3 ] ; then
 export HYBRID=.false.
 export GEN_COORD_HYBRID=.true.
fi
export TFILTC=${TFILTC:-0.85}
export DYNVARS=${DYNVARS:-""}
export PHYVARS=${PHYVARS:-""}
export TRACERVARS=${TRACERVARS:-""}
export FSMCL2=${FSMCL2:-99999}
export FTSFS=${FTSFS:-90}
export FAISS=${FAISS:-99999}
export FSNOL=${FSNOL:-99999}
export FSICL=${FSICL:-99999}
export CYCLVARS=${CYCLVARS}
export POSTGPVARS=${POSTGPVARS}
export NTHREADS=${NTHREADS:-1}
export SEMILAG=${SEMILAG:-${semilag:-.false.}}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-${NTHREADS:-1}}
export SPECTRAL_LOOP=${SPECTRAL_LOOP:-2}
export FILESTYLE=${FILESTYLE:-'L'}
export PGMOUT=${PGMOUT:-${pgmout:-'&1'}}
export PGMERR=${PGMERR:-${pgmerr:-'&2'}}
export MEMBER_NAMES=${MEMBER_NAMES:-''}

export REDOUT=${REDOUT:-'1>'}
export REDERR=${REDERR:-'2>'}
export print_esmf=${print_esmf:-.false.}

################################################################################
#  Preprocessing
$INISCRIPT
pwd=$(pwd)
if [[ -d $DATA ]] ; then
   mkdata=NO
else
   mkdir -p $DATA
   mkdata=YES
fi
cd $DATA||exit 99
if [ $IDEA = .true. ] ; then ${NCP} $FIX_IDEA/global_idea* . ; fi
[[ -d $COMOUT ]]||mkdir -p $COMOUT
################################################################################
#  Make forecast
export PGM='$FCST_LAUNCHER $DATA/$(basename $FCSTEXEC)'
export pgm=$PGM
$LOGSCRIPT
${NCP:-cp} $FCSTEXEC $DATA

#------------------------------------------------------------
if [ $FHROT -gt 0 ] ; then export RESTART=.true. ; fi
export RESTART=${RESTART:-.false.}
if [ $RESTART = .false. ] ; then # when restarting should not remove - Weiyu
  rm -f NULL
fi

FH=$((10#$FHINI))
[[ $FH -lt 10 ]]&&FH=0$FH
if [[ $FHINI -gt 0 ]] ; then
   if [ $FHOUT_HF -ne $FHOUT -a $FH -lt $FHMAX_HF ] ; then
    FH=$((10#$FHINI+10#$FHOUT_HF))
   else
    FH=$((10#$FHINI+10#$FHOUT))
   fi
   [[ $FH -lt 10 ]]&&FH=0$FH
fi
while [[ 10#$FH -le $FHMAX ]] ; do
   if [[ $FH -le $HOUTA ]] ; then
     FNSUB=$NMSUB
   else
     FNSUB=""
   fi
   if [ $DOIAU = YES ]; then
     if [ 10#$FH -lt 10#6 ]; then
       FHIAU=$((10#6-10#$FH))
       FHIAU=m$FHIAU
     else
       FHIAU=$((10#$FH-10#6))
       [[ $FHIAU -lt 10 ]]&&FHIAU=0$FHIAU
     fi
   else
     FHIAU=$FH
   fi
   eval rm -f ${LOGO}${FNSUB}
   if [ $FHOUT_HF -ne $FHOUT -a $FH -lt $FHMAX_HF ] ; then
     ((FH=10#$FH+10#$FHOUT_HF))
   else
     ((FH=10#$FH+10#$FHOUT))
   fi
   [[ $FH -lt 10 ]]&&FH=0$FH
done
if [[ $FILESTYLE = "L" ]] ; then
   #ln -fs $CO2CON fort.15
   #ln -fs $MTNVAR fort.24
   #ln -fs $O3FORC fort.28
   #ln -fs $O3CLIM fort.48

   ${NCP} $CO2CON fort.15
   ${NCP} $MTNVAR fort.24
   ${NCP} $O3FORC fort.28
   ${NCP} $O3CLIM fort.48
else
  echo 'FILESTYLE' $FILESTYLE 'NOT SUPPORTED'
  exit 222
fi

#for m in 01 02 03 04 05 06 07 08 09 10 11 12
#do
# ln -fs $AERODIR/global_aeropac3a.m$m.txt aeropac3a.m$m
#done

AEROSOL_FILE=${AEROSOL_FILE:-global_climaeropac_global.txt}
EMMISSIVITY_FILE=${EMMISSIVITY_FILE:-global_sfc_emissivity_idx.txt}

#ln -fs $AERODIR/$AEROSOL_FILE     aerosol.dat
#ln -fs $EMISDIR/$EMMISSIVITY_FILE sfc_emissivity_idx.txt
#ln -fs $OROGRAPHY                 orography
#ln -fs $OROGRAPHY_UF              orography_uf
#ln -fs $LONSPERLAT                lonsperlat.dat
#ln -fs $LONSPERLAR                lonsperlar.dat
${NCP} $AERODIR/$AEROSOL_FILE     aerosol.dat
${NCP} $EMISDIR/$EMMISSIVITY_FILE sfc_emissivity_idx.txt
${NCP} $OROGRAPHY                 orography
${NCP} $OROGRAPHY_UF              orography_uf
${NCP} $LONSPERLAT                lonsperlat.dat
${NCP} $LONSPERLAR                lonsperlar.dat


if [ $IEMS -gt 0 ] ; then
 EMMISSIVITY_FILE=${EMMISSIVITY_FILE:-global_sfc_emissivity_idx.txt}
 #ln -fs $EMISDIR/$EMMISSIVITY_FILE sfc_emissivity_idx.txt
  ${NCP} $EMISDIR/$EMMISSIVITY_FILE sfc_emissivity_idx.txt
fi
if [ $ISOL -gt 0 ] ; then
 cd $SOLCDIR
 for file in `ls | grep solarconstant` ; do
  ${NCP} $file $DATA/$(echo $file |sed -e "s/global_//g")
 done
fi
if [ $IAER -gt 0 ] ; then
 cd $VOLCDIR
 for file in `ls | grep volcanic_aerosols` ; do
  ${NCP:-cp} $file $DATA/$(echo $file |sed -e "s/global_//g")
 done
 cd $DATA
 #${NCP:-cp} $GOCART_CLIM/* $DATA
 ${NCP:-cp} $GOCART_LUTS/NCEP_AEROSOL.bin $DATA
fi
if [ $ICO2 -gt 0 ] ; then
 cd $CO2DIR
 for file in `ls | grep co2historicaldata` ; do
  ${NCP:-cp} $file $DATA/$(echo $file |sed -e "s/global_//g")
 done
 CO2_seasonal_cycle=${CO2_seasonal_cycle:-global_co2monthlycyc1976_2006.txt}
 ${NCP} $CO2_seasonal_cycle $DATA/co2monthlycyc.txt
fi
cd $DATA
export PHYVARS="IEMS=$IEMS,ISOL=$ISOL,IAER=$IAER,ICO2=$ICO2,$PHYVARS"

#
#     For one member case i.e. control
#     --------------------------------
mins=$((DELTIM/60))
secs=$((DELTIM-(DELTIM/60)*60))
[[ $mins -lt 10 ]] &&mins=0$mins
[[ $secs -lt 10 ]] &&secs=0$secs
export FHINI=$((FHINI+0))
export FHROT=$((FHROT+0))

if [[ $ENS_NUM -le 1 ]] ; then
  FH=$((10#$FHINI))
  [[ $FH -lt 10 ]]&&FH=0$FH
  if [[ $FHINI -gt 0 ]] ; then
    if [ $FHOUT_HF -ne $FHOUT -a $FH -lt $FHMAX_HF ] ; then
     FH=$((10#$FHINI+10#$FHOUT_HF))
    else
     FH=$((10#$FHINI+10#$FHOUT))
    fi
    [[ $FH -lt 10 ]]&&FH=0$FH
  fi
#        For Initial Conditions
#        ----------------------
  if [ $FHINI -eq  $FHROT ]; then
    if [ $NEMSIO_IN = .true. ]; then
      ln -fs $GRDI  grid_ini
      ln -fs $SIGI  sig_ini
    else
      ln -fs $SIGI  sig_ini
    fi
    ln -fs $SFCI  sfc_ini
    ln -fs $NSTI  nst_ini
    if [ $FHROT -gt 0 ] ; then
      export RESTART=.true.
    else
      export RESTART=.false.
    fi
  else
    ln -fs $GRDI  grid_ini
    ln -fs $GRDI2 grid_ini2
    ln -fs $SIGI  sig_ini
    ln -fs $SIGI2 sig_ini2
    ln -fs $SFCI  sfc_ini
    ln -fs $NSTI  nst_ini
    ln -fs $PLASI ipe_grid_plasma_params
    ln -fs $NEUTI ipe_grid_neutral_params
    export RESTART=.true.
  fi
#        For output
#        ----------
  while [[ 10#$FH -le $FHMAX ]] ; do
    if [[ $FH -le $HOUTA ]] ; then
      FNSUB=$NMSUB
    else
      FNSUB=""
    fi
    if [ $DOIAU = YES ]; then
      if [ 10#$FH -lt 10#6 ]; then
        FHIAU=$((10#6-10#$FH))
        FHIAU=m$FHIAU
      else
        FHIAU=$((10#$FH-10#6))
        [[ $FHIAU -lt 10 ]]&&FHIAU=0$FHIAU
      fi
    else
      FHIAU=$FH
    fi
    if [ $FH -eq 00 ] ; then
      SUF2=:${mins}:${secs}
    else
      SUF2=""
    fi
    eval ln -fs ${SIGO}$FNSUB SIG.F${FH}$SUF2
    eval ln -fs ${SFCO}$FNSUB SFC.F${FH}$SUF2
    eval ln -fs ${FLXO}$FNSUB FLX.F${FH}$SUF2
    eval ln -fs ${LOGO}$FNSUB LOG.F${FH}$SUF2
    eval ln -fs ${D3DO}$FNSUB D3D.F${FH}$SUF2
    eval ln -fs ${NSTO}$FNSUB NST.F${FH}$SUF2
    eval ln -fs ${AERO}$FNSUB AER.F${FH}$SUF2
    eval ln -fs ${

    if [ $FHOUT_HF -ne $FHOUT -a $FH -lt $FHMAX_HF ] ; then
     ((FH=10#$FH+10#$FHOUT_HF))
    else
     ((FH=10#$FH+10#$FHOUT))
    fi
    [[ $FH -lt 10 ]]&&FH=0$FH

  done
  eval ln -fs $GRDR1 GRDR1
  eval ln -fs $GRDR2 GRDR2
  eval ln -fs $SIGR1 SIGR1
  eval ln -fs $SIGR2 SIGR2
  eval ln -fs $SFCR  SFCR
  eval ln -fs $NSTR  NSTR
else
#
#   For Ensemble runs (members > 1)
#   -------------------------------
  for MN in $MEMBER_NAMES ; do
    IMN=`echo $MN|cut -c2-3`
    IMN=$((IMN+0))
    if [ $IMN -eq 0 ] ; then IMN=$ENS_NUM ; fi
    if [ $IMN -lt 10 ] ; then IMN=0$IMN ; fi
    echo 'IMN=' $IMN
    if [ ${USESUBDIR:-NO} = YES ] ; then
     if [ $IMN -eq $ENS_NUM ] ; then LOCD=/c00 ; else LOCD=/p$IMN ; fi
    fi
    mkdir -p `eval echo \$COMENS`

#      This is just faking the ensemble ICs.
#   ${NCP:-cp} $SIGI  ${SIGI}${MN}
#   ${NCP:-cp} $SFCI  ${SFCI}${MN}
#   ${NCP:-cp} $SIGI2 ${SIGI2}${MN}
#        For Initial Conditions
#        ----------------------
    eval ln -fs ${GRDI}${MN}  grid_ini_${IMN}
    eval ln -fs ${GRDI2}${MN} grid_ini2_${IMN}
    eval ln -fs ${SIGI}${MN}  sig_ini_${IMN}
    eval ln -fs ${SIGI2}${MN} sig_ini2_${IMN}
    eval ln -fs ${SFCI}${MN}  sfc_ini_${IMN}
    eval ln -fs ${NSTI}${MN}  nst_ini_${IMN}

    if [ $FHINI -eq  $FHROT ]; then
      if [ $FHROT -gt 0 ] ; then
        export RESTART=.true.
      else
        export RESTART=.false.
      fi
    else
      export RESTART=.true.
    fi

#        For output
#        ----------
    FH=$((10#$FHINI))
    [[ $FH -lt 10 ]]&&FH=0$FH
    if [[ $FHINI -gt 0 ]] ; then
      if [ $FHOUT_HF -ne $FHOUT -a $FH -lt $FHMAX_HF ] ; then
       FH=$((10#$FHINI+10#$FHOUT_HF))
      else
       FH=$((10#$FHINI+10#$FHOUT))
      fi
      [[ $FH -lt 10 ]]&&FH=0$FH
    fi
    while [[ $FH -le $FHMAX ]] ; do
      if [[ $FH -le $HOUTA ]] ; then
        FNSUB=$NMSUB
      else
        FNSUB=""
      fi
      if [ $DOIAU = YES ]; then
        if [ 10#$FH -lt 10#6 ]; then
          FHIAU=$((10#6-10#$FH))
          FHIAU=m$FHIAU
        else
          FHIAU=$((10#$FH-10#6))
          [[ $FHIAU -lt 10 ]]&&FHIAU=0$FHIAU
        fi
      else
        FHIAU=$FH
      fi

      eval ln -fs ${SIGO}$FNSUB SIG.F${FH}${SUF2}_${IMN}
      eval ln -fs ${SFCO}$FNSUB SFC.F${FH}${SUF2}_${IMN}
      eval ln -fs ${FLXO}$FNSUB FLX.F${FH}${SUF2}_${IMN}
      eval ln -fs ${LOGO}$FNSUB LOG.F${FH}${SUF2}_${IMN}
      eval ln -fs ${D3DO}$FNSUB D3D.F${FH}${SUF2}_${IMN}
      eval ln -fs ${NSTO}$FNSUB NST.F${FH}${SUF2}_${IMN}
      eval ln -fs ${AERO}$FNSUB AER.F${FH}${SUF2}_${IMN}

      if [ $FHOUT_HF -ne $FHOUT -a $FH -lt $FHMAX_HF ] ; then
       ((FH=10#$FH+10#$FHOUT_HF))
      else
       ((FH=10#$FH+10#$FHOUT))
      fi
      [[ $FH -lt 10 ]]&&FH=0$FH
    done

# 02/29/2008 DHOU,  added new files for the output after SPS
#     FH=$FHMAX
# 09/09/2008 DHOU,  changed the time of output after SPS, fro end to earlier
#             for the digital filtering after resolution change
      FH=$HOUTASPS
    if [[ $FH -lt 10000 ]] ; then
      [[ $FH -lt 10 ]]&&FH=0$FH
      eval ln -fs $SIGS SIG.S${FH}_${IMN}
      eval ln -fs $SFBS SFB.S${FH}_${IMN}
      eval ln -fs $FLXS FLX.S${FH}_${IMN}
    fi

    eval ln -fs ${SIGR1}_${IMN} SIGR1_${IMN}
    eval ln -fs ${SIGR2}_${IMN} SIGR2_${IMN}
    eval ln -fs ${SFCR}_{IMN}   SFCR_${IMN}
    eval ln -fs ${NSTR}_${IMN}  NSTR_${IMN}
# 02/29/2008 DHOU,  added new files for the re-start-files after SP
    if [ $ENS_NUM -gt 2 ] ; then
      eval ln -fs ${SIGS1}_${IMN} SIGS1_${IMN}
      eval ln -fs ${SIGS2}_${IMN} SIGS2_${IMN}
      eval ln -fs ${SFCS}_${IMN}  SFCS_${IMN}
      eval ln -fs ${NSTS}_${IMN}  NSTS_${IMN}
    fi
  done
fi

#
# Create Configure file (i.e. .rc file) here
# PE$n are to be imported from outside.  If PE$n are not set from outside, the
# model would give equal processors for all ensembel members.
#
c=1
while [ $c -le $ENS_NUM ] ; do
 eval export PE$c=\${PE$c:-0}
 c=$((c+1))
done

export wgrib=${wgrib:-$NWPROD/util/exec/wgrib}

if [ $FHINI -eq 0 ]; then
  if [[ $ENS_NUM -le 1 ]] ; then
    if [ $NEMSIO_IN = .true. ]; then
     if [ $ioform_sig = 'grib' ] ; then
      export CDATE_NEMS=$($wgrib -4yr $GRDI | grep -i hgt |awk -F: '{print $3}' |awk -F= '{print $2}')
     else
      export CDATE_NEMS=$($nemsioget $GRDI idate |grep -i "idate" |awk -F= '{print $2}')
     fi
    else
      export CDATE_SIG=$(echo idate|$SIGHDR ${SIGI})
    fi
  else
    MN=c00
    if [ $NEMSIO_IN = .true. ]; then
     if [ $ioform_sig = 'grib' ] ; then
      export CDATE_NEMS=$($wgrib -4yr ${GRDI}${MN} | grep -i hgt |awk -F: '{print $3}' |awk -F= '{print $2}')
     else
      export CDATE_NEMS=$($nemsioget $GRDI idate |grep -i "idate" |awk -F= '{print $2}')
     fi
    else
      export CDATE_SIG=$(echo idate|$SIGHDR ${SIGI}i${MN})
    fi
  fi
else
  if [[ $ENS_NUM -le 1 ]] ; then
    if [ $NEMSIO_IN = .true. ]; then
     if [ $ioform_sig = 'grib' ] ; then
      export CDATE_NEMS=$($wgrib -4yr $GRDI | grep -i hgt |awk -F: '{print $3}' |awk -F= '{print $2}')
     else
      export CDATE_NEMS=$($nemsioget $GRDI idate | grep -i "idate" |awk -F= '{print $2}')
     fi
    else
      export CDATE_SIG=$(echo idate|$SIGHDR ${SIGI})
    fi
  else
    MN=c00
    if [ $NEMSIO_IN = .true. ]; then
     if [ $ioform_sig = 'grib' ] ; then
      export CDATE_NEMS=$($wgrib -4yr ${GRDI}${MN} | grep -i hgt |awk -F: '{print $3}' |awk -F= '{print $2}')
     else
      export CDATE_NEMS=$($nemsioget ${GRDI}${MN} idate | grep -i "idate" |awk -F= '{print $2}')
     fi
    else
      export CDATE_SIG=$(echo idate|$SIGHDR ${SIGI}${MN})
    fi
  fi
fi

if [ $NEMSIO_IN = .true.  ]; then
    INI_YEAR=$(echo $CDATE_NEMS | awk -F" " '{print $1}')
    echo "cdate=$CDATE_NEMS, ini_year=${INI_YEAR}"
    INI_MONTH=$(echo $CDATE_NEMS | awk -F" " '{print $2}')
    INI_DAY=$(echo $CDATE_NEMS | awk -F" " '{print $3}')
    INI_HOUR=$(echo $CDATE_NEMS | awk -F" " '{print $4}')
    yyyy=$INI_YEAR
    mm=$INI_MONTH; if [ $mm -lt 10 ]; then mm=0$mm ; fi
    dd=$INI_DAY; if [ $dd -lt 10 ]; then dd=0$dd ; fi
    hh=`expr $INI_HOUR + 0 `; if [ $hh -lt 10 ]; then hh=0$hh ; fi
    export CDATE=${yyyy}${mm}${dd}${hh}
else
    INI_YEAR=$(echo $CDATE_SIG | cut -c1-4)
    echo "cdate=$CDATE_SIG, ini_year=${INI_YEAR}"
    INI_MONTH=$(echo $CDATE_SIG | cut -c5-6)
    INI_DAY=$(echo $CDATE_SIG | cut -c7-8)
    INI_HOUR=$(echo $CDATE_SIG | cut -c9-10)
    export CDATE=$CDATE_SIG


#wanghj
 NWPROD=/nwprod
 NDATE=$NWPROD/util/exec/ndate
 CDATE0=`$NDATE +06 $CDATE`

 INI_YEAR0=$(echo $CDATE0 | cut -c1-4)
 INI_MONTH0=$(echo $CDATE0 | cut -c5-6)
 INI_DAY0=$(echo $CDATE0 | cut -c7-8)
 INI_HOUR0=$(echo $CDATE0 | cut -c9-10)

 CDATE3=`$NDATE +03 $CDATE`
 INI_YEAR3=$(echo $CDATE3 | cut -c1-4)
 INI_MONTH3=$(echo $CDATE3 | cut -c5-6)
 INI_DAY3=$(echo $CDATE3 | cut -c7-8)
 INI_HOUR3=$(echo $CDATE3 | cut -c9-10)

 START_UT_SEC=`expr ${INI_HOUR3} \* 3600`


fi

## copy configure files needed for NEMS GFS
${NCP} ${MAPL:-$PARM_NGAC/MAPL.rc}                    MAPL.rc
${NCP} ${CHEM_REGISTRY:-$PARM_NGAC/Chem_Registry.rc}  Chem_Registry.rc

## copy configure files and fixed files needed for GOCART
if [ $GOCART == 1 ] ; then
 ${NCP} ${CONFIG_DU:-$PARM_NGAC/DU_GridComp.rc}             DU_GridComp.rc
 ${NCP} ${CONFIG_SU:-$PARM_NGAC/SU_GridComp.rc}             SU_GridComp.rc
 ${NCP} ${CONFIG_OC:-$PARM_NGAC/OC_GridComp.rc}             OC_GridComp.rc
 ${NCP} ${CONFIG_OCx:-$PARM_NGAC/OC_GridComp---full.rc}     OC_GridComp---full.rc
 ${NCP} ${CONFIG_BC:-$PARM_NGAC/BC_GridComp.rc}             BC_GridComp.rc
 ${NCP} ${CONFIG_SS:-$PARM_NGAC/SS_GridComp.rc}             SS_GridComp.rc
 ${NCP} ${AOD_REGISTRY:-$PARM_NGAC/Aod-550nm_Registry.rc}   $DATA/Aod_Registry.rc
#jw  ${NCP} ${AOD_REGISTRY:-$PARM_NGAC/Aod-550nm_Registry.rc}   $DATA/Aod-550nm_Registry.rc
 ${NCP} $PARM_NGAC/AEROSOL_LUTS.dat                         $DATA/

 ln -sf $FIX_NGAC  ngac_fix
fi

# Add time dependent variables F10.7 and kp into the WAM model.
#--------------------------------------------------------------
if [ $IDEA = .true. ]; then
# ${NCP} $WAMINDIR/wam_input-${INI_YEAR}${INI_MONTH}${INI_DAY}T${INI_HOUR}15.xml $DATA/wam_input2.xsd
  ${NCP} $WAMINDIR/wam_input-${INI_YEAR0}${INI_MONTH0}${INI_DAY0}T${INI_HOUR0}15.xml $DATA/wam_input2.xsd

# ${NCP} $WAMINDIR/../wanghj/wam_input-${INI_YEAR}${INI_MONTH}${INI_DAY}T${INI_HOUR}15.xml $DATA/wam_input2.xsd
  ${NCP} $WAMINDIR/../wanghj/wam_input-${INI_YEAR0}${INI_MONTH0}${INI_DAY0}T${INI_HOUR0}15.xml $DATA/wam_input2.xsd


  if [ $machine = WCOSS ] ; then
    PERL=/u/Weiyu.Yang/bin/perl/bin/perl-static
  elif [ $machine = THEIA ] ; then
    PERL=perl
  fi
  $PERL $BASE_GSM/ush/WAM_XML_to_ASCII.pl
  ${NCP} $DATA/wam_input.asc $DATA/wam_input_f107_kp.txt


# wanghj - copy some other data
${NCP} ${EXPDIR}/WAM_INPUT/* $DATA


#IAU stuff
if [ $DOIAU = YES ]; then

   export RESTART=.false.

  #ln -sf ${ROTDIR}/sfcf03.gdas.${CDATE} ${DATA}/sfc_ini

   export FHRES=3
   export FHOUT=1
   export FHZER=3

fi



export DATADIR=/nems/noscrub/emc.nemspara/RT/WAM-IPE/WAM-IPE_NEMS201606-20170131/data
/bin/cp -f ${DATADIR}/MED_SPACEWX/gsm%wam%T62_ipe%80x170/wam3dgridnew2.nc ${DATA}


if [ $WAM_IPE_COUPLING = .true. ]; then

   /bin/cp -f ${DATADIR}/MED_SPACEWX/gsm%wam%T62_ipe%80x170/ipe3dgrid2.nc ${DATA}/.

   ln -sf ${DATADIR}/MED_SPACEWX/gsm%wam%T62_ipe%80x170/wam3dgridnew_20160427.nc ${DATA}/.

   ln -sf ${DATADIR}/IPE/GIP_apex_coords_global_lowres_new20120705 ${DATA}/ipe_grid


#wanghj - need to update this
#  ln -sf ${DATADIR}/IPE/cases/20090115_1day_spacewx_80x170/stup* ${DATA}/.
#  ln -sf ${COMOUT}/ipe${CDATE}/stup* ${DATA}/.

   export ARCION=${NOSCRUB}/${LOGNAME}/archive/ipe201707

#  ln -sf ${ARCION}/ipe${CDATE}/stup* ${DATA}/.

ln -sf ${ARCION}/ipe${CDATE}/plasma00 ${DATA}/stup00
ln -sf ${ARCION}/ipe${CDATE}/plasma01 ${DATA}/stup01
ln -sf ${ARCION}/ipe${CDATE}/plasma02 ${DATA}/stup02
ln -sf ${ARCION}/ipe${CDATE}/plasma03 ${DATA}/stup03
ln -sf ${ARCION}/ipe${CDATE}/plasma04 ${DATA}/stup04
ln -sf ${ARCION}/ipe${CDATE}/plasma05 ${DATA}/stup05
ln -sf ${ARCION}/ipe${CDATE}/plasma06 ${DATA}/stup06
ln -sf ${ARCION}/ipe${CDATE}/plasma07 ${DATA}/stup07
ln -sf ${ARCION}/ipe${CDATE}/plasma08 ${DATA}/stup08
ln -sf ${ARCION}/ipe${CDATE}/plasma09 ${DATA}/stup09
ln -sf ${ARCION}/ipe${CDATE}/plasma10 ${DATA}/stup10
ln -sf ${ARCION}/ipe${CDATE}/plasma11 ${DATA}/stup11
ln -sf ${ARCION}/ipe${CDATE}/plasma12 ${DATA}/stup12
ln -sf ${ARCION}/ipe${CDATE}/plasma13 ${DATA}/stup13
ln -sf ${ARCION}/ipe${CDATE}/plasma14 ${DATA}/stup14
ln -sf ${ARCION}/ipe${CDATE}/plasma15 ${DATA}/stup15
ln -sf ${ARCION}/ipe${CDATE}/plasma16 ${DATA}/stup16

ln -sf ${ARCION}/ipe${CDATE}/ut_rec   ${DATA}/stup_ut_rec


#for neutral
#ln -sf ${ARCION}/ipe${CDATE}/ipe_grid_neut_params_ut ${DATA}/stup_ipe_grid_neut_params_ut

ln -sf ${ARCION}/ipe${CDATE}/ipe_grid_neut_N2_den ${DATA}/stup_ipe_grid_neut_N2_den
ln -sf ${ARCION}/ipe${CDATE}/ipe_grid_neut_O2_den ${DATA}/stup_ipe_grid_neut_O2_den
ln -sf ${ARCION}/ipe${CDATE}/ipe_grid_neut_O_den  ${DATA}/stup_ipe_grid_neut_O_den
ln -sf ${ARCION}/ipe${CDATE}/ipe_grid_neut_temp   ${DATA}/stup_ipe_grid_neut_temp  
ln -sf ${ARCION}/ipe${CDATE}/ipe_grid_neut_wind   ${DATA}/stup_ipe_grid_neut_wind  


    export FHINI_SFC=${CDATE_SFC:-$(echo fhour|$SFCHDR ${SFCI}$FM)}

   /bin/cp -f ${EXPDIR}/IPE_INPUT/SMSnamelist ${DATA}

   /bin/cp -f ${EXPDIR}/IPE_INPUT/coeff* ${DATA}
   /bin/cp -f ${EXPDIR}/IPE_INPUT/wei96* ${DATA}

#  /bin/cp -f /swpc/save/$LOGNAME/para_gfs/prwam/IPE_INPUT/IPE.inp ${DATA}


#wanghj
export coupling_interval_fast_sec=180.0
export PE1=64
export FHDFI=0
export nhours_dfini=$FHDFI
export WRT_GROUP=1
export WRTPE_PER_GROUP=3
export FILE_IO_FORM="'grib' 'bin4' 'grib'"

#export ESMF_RUNTIME_COMPLIANCECHECK=ON:depth=4
#export TASKSIZE=64


#this is for previous day, need to change to reflect real run DOY
#export DOY=`date -d ${INI_MONTH}/${INI_DAY}/${INI_YEAR} +%j`
 export DOY=`date -d ${INI_MONTH3}/${INI_DAY3}/${INI_YEAR3} +%j`


#extract f107 & kp for ionosphere model from f107/kp data file

F107AVG=`grep 'F10 81 Day Avg' $DATA/wam_input_f107_kp.txt | awk '{print $5}'`

F107DAY=`grep ${INI_YEAR3}-${INI_MONTH3}-${INI_DAY3}T${INI_HOUR3} $DATA/wam_input_f107_kp.txt | awk '{print $2}'`

KP3HR=`grep ${INI_YEAR3}-${INI_MONTH3}-${INI_DAY3}T${INI_HOUR3} $DATA/wam_input_f107_kp.txt | awk '{print $3}'`

#need to convert from kp to ap

#just some approximation first
AP3HR=`expr ${KP3HR} \* 4`

echo 'F107AVG, F107DAY, KPX, APX = ' ${F107AVG}, ${F107DAY}, ${KP3HR}, ${AP3HR}


cat  > IPE.inp <<EOF
&ipedims
  nlp=170
  nmp=80
  npts2d=44514
/
&nmflip
  colfac_flip=1.7
  dt_init_guess_flip=60
  dtmin_flip=1.0
  fnfac_flip=1.0
  fpas_flip=0.0
  heprat_flip=0.09
  hpeq_flip=0.0
  ht_lce=200.0
  init_te_max=3000.0
  start_time_depleted=0
  sw_debug_flip=0
  sw_depleted_flip=0
  sw_erstop_flip=0
  sw_ihepls=1
  sw_init_guess_flip=f
  sw_inno=-1
  sw_inpls=1
  sw_lce=f
  sw_neutral_heating_flip=0
  sw_ohpls=1
  sw_optw_flip=t
  sw_tei=1
  sw_wind_flip=1
  zlbdy_flip=120.00
  zlbnp_inp=115.00
/
&nmipe
  f107av=${F107AVG}
  f107d=${F107DAY}
  internalTimeLoopMax=1
  ip_freq_eldyn=180
  ip_freq_msis=180
  ip_freq_plasma=180
  ip_freq_output=3600
  nday=$DOY
  nyear=2000
  start_time=${START_UT_SEC}
  time_step=180
/
  &nmmsis
  ap(1)=${AP3HR}
  ap(2)=${AP3HR}
  ap(3)=${AP3HR}
  ap(4)=${AP3HR}
  ap(5)=${AP3HR}
  ap(6)=${AP3HR}
  ap(7)=${AP3HR}
  kp_eld=${KP3HR}
/
&nmswitch
  duration=43200
  fac_bm=1.00
  iout(1)=1
  iout(2)=60
  lpFort167=57
  lpmax_perp_trans=151
  lpmin_perp_trans=15
  mpFort167=71
  mpstop=80
  peFort167=56
  record_number_plasma_start=0
  sw_dbg_perp_trans=f
  sw_debug=f
  sw_debug_mpi=f
  sw_divv=0
  sw_eldyn=1
  sw_exb_up=1
  sw_grid=0
  sw_ksi=2
  sw_neutral=1
  swNeuPar(1)=t
  swNeuPar(2)=t
  swNeuPar(3)=t
  swNeuPar(4)=t
  swNeuPar(5)=f
  swNeuPar(6)=f
  swNeuPar(7)=f
  swEsmfTime=t
  sw_output_fort167=f
  sw_output_wind=t
  sw_output_plasma_grid=f
  sw_use_wam_fields_for_restart=t
  sw_para_transport=1
  sw_pcp=0
  sw_perp_transport=1
  sw_record_number=2
  sw_th_or_r=1
  ut_start_perp_trans=${START_UT_SEC}
/
EOF

fi


fi


#
# jw: generate configure file
#


if [ $WAM_IPE_COUPLING = .true. ]; then

cat << EOF > $DATA/nems.configure
#############################################
####  NEMS Run-Time Configuration File  #####
#############################################

# EARTH #
EARTH_component_list: MED ATM IPM
EARTH_attributes::
  Verbosity = max
::

# MED #
MED_model:                      spaceweather
MED_petlist_bounds:             16 23
MED_attributes::
  Verbosity = max
  DumpFields = false
  DumpRHs = false
::

# ATM #
ATM_model:                      gsm
ATM_petlist_bounds:             0 15
ATM_attributes::
  Verbosity = max
::

# IPM #
IPM_model:                      ipe
IPM_petlist_bounds:             24 63
IPM_attributes::
  Verbosity = max
::

# Run Sequence #
runSeq::
  @180.0
    ATM -> MED :remapMethod=redist
    MED
    MED -> IPM :remapMethod=redist
    ATM
    IPM
  @
::
EOF


else


cat << EOF > $DATA/nems.configure
EARTH_component_list:       ATM
ATM_model:                  ${atm_model:-gsm}
runSeq::
  ATM
::
EOF

fi


core=${core:-gfs}
#cat << EOF > $DATA/atmos.configure
#core: $core
#EOF


if [ $WAM_IPE_COUPLING = .true. ]; then

cat << EOF > $DATA/atmos.configure
core: $core
atm_model:                  ${atm_model:-gsm}
atm_coupling_interval_sec:  ${coupling_interval_fast_sec:-600}
EOF

else

cat << EOF > $DATA/atmos.configure
core: $core
atm_model:                  ${atm_model:-gsm}
atm_coupling_interval_sec:
EOF

fi 


cat << EOF > atm_namelist.rc

core: $core
print_esmf:     ${print_esmf}

nhours_dfini=0

#nam_atm +++++++++++++++++++++++++++
nlunit:                  35
deltim:                  ${DELTIM}.0
fhrot:                   $FHROT
namelist:                atm_namelist
total_member:            $ENS_NUM
grib_input:              $GB
PE_MEMBER01:             $PE1
PE_MEMBER02:             $PE2
PE_MEMBER03:             $PE3
PE_MEMBER04:             $PE4
PE_MEMBER05:             $PE5
PE_MEMBER06:             $PE6
PE_MEMBER07:             $PE7
PE_MEMBER08:             $PE8
PE_MEMBER09:             $PE9
PE_MEMBER10:             $PE10
PE_MEMBER11:             $PE11
PE_MEMBER12:             $PE12
PE_MEMBER13:             $PE13
PE_MEMBER14:             $PE14
PE_MEMBER15:             $PE15
PE_MEMBER16:             $PE16
PE_MEMBER17:             $PE17
PE_MEMBER18:             $PE18
PE_MEMBER19:             $PE19
PE_MEMBER20:             $PE20
PE_MEMBER21:             $PE21

# For stachastic purturbed runs -  added by Dhou and Wyang
  --------------------------------------------------------
#  ENS_SPS, logical control for application of stochastic perturbation scheme
#  HH_START, start hour of forecast, and modified ADVANCECOUNT_SETUP
#  HH_INCREASE and HH_FINAL are fcst hour increment and end hour of forecast
#  ADVANCECOUNT_SETUP is an integer indicating the number of time steps between
#  integrtion_start and the time when model state is saved for the _ini of the
#  GEFS_Coupling, currently is 0h.

HH_INCREASE:             $FH_INC
HH_FINAL:                $FHMAX
HH_START:                $FHINI
ADVANCECOUNT_SETUP:      $ADVANCECOUNT_SETUP

ENS_SPS:                 $ENS_SPS
HOUTASPS:                $HOUTASPS

#ESMF_State_Namelist +++++++++++++++

RUN_CONTINUE:            .false.

#
dt_int:                  $DELTIM
dt_num:                  0
dt_den:                  1
start_year:              $INI_YEAR
start_month:             $INI_MONTH
start_day:               $INI_DAY
start_hour:              $INI_HOUR
start_minute:            0
start_second:            0
nhours_fcst:             $FHMAX
restart:                 $RESTART
nhours_fcst1:            $FHMAX
im:                      $LONB
jm:                      $LATB
global:                  .true.
nhours_dfini:            $nhours_dfini
adiabatic:               $ADIAB
lsoil:                   $LSOIL
passive_tracer:          $PASSIVE_TRACER
dfilevs:                 $DFILEVS
ldfiflto:                $LDFIFLTO
ldfi_grd:                $LDFI_GRD
num_tracers:             $NTRAC
lwrtgrdcmp:              $LWRTGRDCMP
nemsio_in:               $NEMSIO_IN

#jwstart added quilt
###############################
#### Specify the I/O tasks ####
###############################


quilting:                $QUILTING   #For asynchronous quilting/history writes
read_groups:             0
read_tasks_per_group:    0
write_groups:            $WRT_GROUP
write_tasks_per_group:   $WRTPE_PER_GROUP

num_file:                $NUM_FILE                   #
filename_base:           $FILENAME_BASE
file_io_form:            $FILE_IO_FORM                     
file_io:                 'DEFERRED' 'DEFERRED' 'DEFERRED' 'DEFERRED'  #
write_dopost:            $WRITE_DOPOST          # True--> run do on quilt
post_gribversion:        $POST_GRIBVERSION      # True--> grib version for post output files
gocart_aer2post:         $GOCART_AER2POST
write_nemsioflag:        $WRITE_NEMSIOFLAG      # True--> Write nemsio run history files
nfhout:                  $FHOUT
nfhout_hf:               $FHOUT_HF
nfhmax_hf:               $FHMAX_HF
nsout:                   $nsout

io_recl:                 100
io_position:             ' '
io_action:               'WRITE'
io_delim:                ' '
io_pad:                  ' '

#jwend

EOF


# addition import/export variables for stochastic physics
export sppt_import=${sppt_import:-0}
export sppt_export=${sppt_export:-0}
export shum_import=${shum_import:-0}
export shum_export=${shum_export:-0}
export skeb_import=${skeb_import:-0}
export skeb_export=${skeb_export:-0}
export vc_import=${vc_import:-0}
export vc_export=${vc_export:-0}


#
cat atm_namelist.rc > dyn_namelist.rc
cat << EOF >> dyn_namelist.rc

SLG_FLAG:                        $SEMILAG

#ESMF_State_Namelist +++++++++++++++
idate1_import:                    1
z_import:                         1
ps_import:                        1
div_import:                       0
vor_import:                       0
u_import:                         1
v_import:                         1
temp_import:                      1
tracer_import:                    1
p_import:                         1
dp_import:                        1
dpdt_import:                      1

idate1_export:                    1
z_export:                         1
ps_export:                        1
div_export:                       0
vor_export:                       0
u_export:                         1
v_export:                         1
temp_export:                      1
tracer_export:                    1
p_export:                         1
dp_export:                        1
dpdt_export:                      1
sppt_wts_export:                  ${sppt_export}
shum_wts_export:                  ${shum_export}
skeb_wts_export:                  ${skeb_export}
vc_wts_export:                    ${vc_export}

EOF

cat atm_namelist.rc > phy_namelist.rc
cat << EOF >> phy_namelist.rc

#Upper_Air_State_Namelist +++++++++++++++
idate1_import:                    1
z_import:                         1
ps_import:                        1
div_import:                       0
vor_import:                       0
u_import:                         1
v_import:                         1
temp_import:                      1
tracer_import:                    1
p_import:                         1
dp_import:                        1
dpdt_import:                      1
sppt_wts_import:                  ${sppt_import}
shum_wts_import:                  ${shum_import}
skeb_wts_import:                  ${skeb_import}
vc_wts_import:                    ${vc_import}

idate1_export:                    1
z_export:                         1
ps_export:                        1
div_export:                       0
vor_export:                       0
u_export:                         1
v_export:                         1
temp_export:                      1
tracer_export:                    1
p_export:                         1
dp_export:                        1
dpdt_export:                      1

# Surface state.
#---------------
orography_import:                 1
t_skin_import:                    1
soil_mois_import:                 1
snow_depth_import:                1
soil_t_import:                    1
deep_soil_t_import:               1
roughness_import:                 1
conv_cloud_cover_import:          1
conv_cloud_base_import:           1
conv_cloud_top_import:            1
albedo_visible_scattered_import:  1
albedo_visible_beam_import:       1
albedo_nearir_scattered_import:   1
albedo_nearir_beam_import:        1
sea_level_ice_mask_import:        1
vegetation_cover_import:          1
canopy_water_import:              1
m10_wind_fraction_import:         1
vegetation_type_import:           1
soil_type_import:                 1
zeneith_angle_facsf_import:       1
zeneith_angle_facwf_import:       1
uustar_import:                    1
ffmm_import:                      1
ffhh_import:                      1
sea_ice_thickness_import:         1
sea_ice_concentration_import:     1
tprcp_import:                     1
srflag_import:                    1
actual_snow_depth_import:         1
liquid_soil_moisture_import:      1
vegetation_cover_min_import:      1
vegetation_cover_max_import:      1
slope_type_import:                1
snow_albedo_max_import:           1

orography_export:                 1
t_skin_export:                    1
soil_mois_export:                 1
snow_depth_export:                1
soil_t_export:                    1
deep_soil_t_export:               1
roughness_export:                 1
conv_cloud_cover_export:          1
conv_cloud_base_export:           1
conv_cloud_top_export:            1
albedo_visible_scattered_export:  1
albedo_visible_beam_export:       1
albedo_nearir_scattered_export:   1
albedo_nearir_beam_export:        1
sea_level_ice_mask_export:        1
vegetation_cover_export:          1
canopy_water_export:              1
m10_wind_fraction_export:         1
vegetation_type_export:           1
soil_type_export:                 1
zeneith_angle_facsf_export:       1
zeneith_angle_facwf_export:       1
uustar_export:                    1
ffmm_export:                      1
ffhh_export:                      1
sea_ice_thickness_export:         1
sea_ice_concentration_export:     1
tprcp_export:                     1
srflag_export:                    1
actual_snow_depth_export:         1
liquid_soil_moisture_export:      1
vegetation_cover_min_export:      1
vegetation_cover_max_export:      1
slope_type_export:                1
snow_albedo_max_export:           1

EOF


# additional namelist parameters for stochastic physics.  Default is off
export SPPT=${SPPT:-"0.0,0.0,0.0,0.0,0.0"}
export ISEED_SPPT=${ISEED_SPPT:-0}
export SPPT_LOGIT=.TRUE.
export SPPT_LOGIT=${SPPT_LOGIT:-.TRUE.}
export SPPT_TAU=${SPPT_TAU:-"21600,2592500,25925000,7776000,31536000"}
export SPPT_LSCALE=${SPPT_LSCALE:-"500000,1000000,2000000,2000000,2000000"}

export SHUM=${SHUM:-"0.0, -999., -999., -999, -999"}
export ISEED_SHUM=${ISEED_SHUM:-0}
export SHUM_TAU=${SHUM_TAU:-"2.16E4, 1.728E5, 6.912E5, 7.776E6, 3.1536E7"}
export SHUM_LSCALE=${SHUM_LSCALE:-"500.E3, 1000.E3, 2000.E3, 2000.E3, 2000.E3"}

export SKEB=${SKEB:-"0.0, -999., -999., -999, -999"}
export ISEED_SKEB=${ISEED_SKEB:-0}
export SKEB_TAU=${SKEB_TAU:-"2.164E4, 1.728E5, 2.592E6, 7.776E6, 3.1536E7"}
export SKEB_LSCALE=${SKEB_LSCALE:="1000.E3, 1000.E3, 2000.E3, 2000.E3, 2000.E3"}
export SKEB_VFILT=${SKEB_VFILT:-40}
export SKEB_DISS_SMOOTH=${SKEB_DISS_SMOOTH:-12}

export VC=${VC:-0.0}
export ISEED_VC=${ISEED_VC:-0}
export VCAMP=${VCAMP:-"0.0, -999., -999., -999, -999"}
export VC_TAU=${VC_TAU:-"4.32E4, 1.728E5, 2.592E6, 7.776E6, 3.1536E7"}
export VC_LSCALE=${VC_LSCALE:-"1000.E3, 1000.E3, 2000.E3, 2000.E3, 2000.E3"}


#
#   WARNING WARNING FILESTYLE "C" will not work for Component Ensembles!!!
#
#eval $PGM <<EOF $REDOUT$PGMOUT $REDERR$PGMERR
#totalview poe -a $PGM <<EOF $REDOUT$PGMOUT $REDERR$PGMERR
#
cat  > atm_namelist <<EOF
 &nam_dyn
  FHOUT=$FHOUT, FHMAX=$FHMAX, IGEN=$IGEN, DELTIM=$DELTIM,
  FHRES=$FHRES, FHROT=$FHROT, FHDFI=$FHDFI, nsout=$nsout,
  nxpt=1, nypt=2, jintmx=2, lonf=$LONF, latg=$LATG,
  jcap=$JCAP, levs=$LEVS,  levr=$LEVR,
  ntrac=$NTRAC, ntoz=$NTOZ, ntcw=$NTCW, ncld=$NCLD,
  ngptc=$NGPTC, hybrid=$HYBRID, tfiltc=$TFILTC,
  gen_coord_hybrid=$GEN_COORD_HYBRID, zflxtvd=$ZFLXTVD,
  spectral_loop=$SPECTRAL_LOOP, explicit=$EXPLICIT,
  ndslfv=$NDSLFV,mass_dp=$MASS_DP,process_split=$PROCESS_SPLIT,
  reduced_grid=$REDUCED_GRID,lsidea=$IDEA,
  wam_ipe_coupling=$WAM_IPE_COUPLING,
  height_dependent_g=$HEIGHT_DEPENDENT_G,
  semi_implicit_temp_profile=$SEMI_IMPLICIT_TEMP_PROFILE,
  thermodyn_id=$THERMODYN_ID, sfcpress_id=$SFCPRESS_ID,
  dfilevs=$DFILEVS, liope=$liope,
  FHOUT_HF=$FHOUT_HF, FHMAX_HF=$FHMAX_HF,
  $DYNVARS /
 &nam_phy
  FHOUT=$FHOUT, FHMAX=$FHMAX, IGEN=$IGEN, DELTIM=$DELTIM,
  DTPHYS=$DTPHYS,
  FHRES=$FHRES, FHROT=$FHROT, FHCYC=$FHCYC, FHDFI=$FHDFI,
  FHZER=$FHZER, FHLWR=$FHLWR, FHSWR=$FHSWR,nsout=$nsout,
  nxpt=1, nypt=2, jintmx=2, lonr=$LONR, latr=$LATR,
  jcap=$JCAP, levs=$LEVS, levr=$LEVR, reduced_grid=$REDUCED_GRID,
  ntrac=$NTRAC, ntoz=$NTOZ, ntcw=$NTCW, ncld=$NCLD,
  lsoil=$LSOIL, nmtvr=$NMTVR, lsidea=$IDEA,
  f107_kp_size=$F107_KP_SIZE,
  f107_kp_interval=$F107_KP_INTERVAL,
  f107_kp_skip_size=$F107_KP_SKIP_SIZE,
  ngptc=$NGPTC, hybrid=$HYBRID, tfiltc=$TFILTC,
  gen_coord_hybrid=$GEN_COORD_HYBRID,
  thermodyn_id=$THERMODYN_ID, sfcpress_id=$SFCPRESS_ID,
  FHOUT_HF=$FHOUT_HF, FHMAX_HF=$FHMAX_HF,
  nstf_name=${nstf_name},liope=$liope,
  $PHYVARS /
 &TRACER_CONSTANT
  $TRACERVARS /
 &SOIL_VEG
  LPARAM = .FALSE./
 &NAMSFC
  FNGLAC="$FNGLAC",
  FNMXIC="$FNMXIC",
  FNTSFC="$FNTSFC",
  FNSNOC="$FNSNOC",
  FNZORC="$FNZORC",
  FNALBC="$FNALBC",
  FNAISC="$FNAISC",
  FNTG3C="$FNTG3C",
  FNVEGC="$FNVEGC",
  FNVETC="$FNVETC",
  FNSOTC="$FNSOTC",
  FNSMCC="$FNSMCC",
  FNMSKH="$FNMSKH",
  FNTSFA="$FNTSFA",
  FNACNA="$FNACNA",
  FNSNOA="$FNSNOA",
  FNVMNC="$FNVMNC",
  FNVMXC="$FNVMXC",
  FNSLPC="$FNSLPC",
  FNABSC="$FNABSC",
  LDEBUG=.false.,
  FSMCL(2)=$FSMCL2,
  FSMCL(3)=$FSMCL2,
  FSMCL(4)=$FSMCL2,
  FTSFS=$FTSFS,
  FAISS=$FAISS,
  FSNOL=$FSNOL,
  FSICL=$FSICL,
  FTSFL=99999,
  FAISL=99999,
  FVETL=99999,
  FSOTL=99999,
  FvmnL=99999,
  FvmxL=99999,
  FSLPL=99999,
  FABSL=99999,
  FSNOS=99999,
  FSICS=99999,

  $CYCLVARS  /
 &NAMPGB
  $POSTGPVARS /
EOF

ln -sf atm_namelist.rc ./model_configure


ulimit -s unlimited
ulimit -c unlimited


### AMK affect the sfc file so it matches the sigf initial date and forecast hour
 if [ $DOIAU = YES ]; then
    export CDATE_SFC=${CDATE_SFC:-$(echo idate|$SFCHDR ${SFCI}$FM)}
   #export FHINI_SFC=${CDATE_SFC:-$(echo fhour|$SFCHDR ${SFCI}$FM)}
    export FHINI_SFC=${FHINI_SFC:-$(echo fhour|$SFCHDR ${SFCI}$FM)}
    eval $CHGSFCFHREXEC $SFCI $CDATE_SIG $FHINI
 fi


eval $FCSTENV $PGM $REDOUT$PGMOUT $REDERR$PGMERR


### AMK change sfc file back
 if [ $DOIAU = YES ]; then
    eval $CHGSFCFHREXEC $SFCI $CDATE_SFC $FHINI_SFC
 fi


if [ $WAM_IPE_COUPLING = .true. ]; then

#wanghj: copy ionosphere data to ROTDIR/COMOUT
mkdir -p ${ARCION}/ipe${CDATE0}
/bin/cp -f ${DATA}/plasma* ${ARCION}/ipe${CDATE0}
/bin/cp -f ${DATA}/ut_rec  ${ARCION}/ipe${CDATE0}

#any other data to copy???
/bin/cp -f ${DATA}/IPE.inp ${ARCION}/ipe${CDATE0}
/bin/cp -f -p ${DATA}/FLIP_ERR ${ARCION}/ipe${CDATE0}
/bin/cp -f -p ${DATA}/logfile_input_params.log ${ARCION}/ipe${CDATE0}

#also copy neutral data here too
/bin/cp -f -p ${DATA}/ipe_grid_neut* ${ARCION}/ipe${CDATE0}


#also get ready for the next cycle? this could be done more compactly

#ln -sf ${ARCION}/ipe${CDATE0}/plasma00 ${ARCION}/ipe${CDATE0}/stup00
#ln -sf ${ARCION}/ipe${CDATE0}/plasma01 ${ARCION}/ipe${CDATE0}/stup01
#ln -sf ${ARCION}/ipe${CDATE0}/plasma02 ${ARCION}/ipe${CDATE0}/stup02
#ln -sf ${ARCION}/ipe${CDATE0}/plasma03 ${ARCION}/ipe${CDATE0}/stup03
#ln -sf ${ARCION}/ipe${CDATE0}/plasma04 ${ARCION}/ipe${CDATE0}/stup04
#ln -sf ${ARCION}/ipe${CDATE0}/plasma05 ${ARCION}/ipe${CDATE0}/stup05
#ln -sf ${ARCION}/ipe${CDATE0}/plasma06 ${ARCION}/ipe${CDATE0}/stup06
#ln -sf ${ARCION}/ipe${CDATE0}/plasma07 ${ARCION}/ipe${CDATE0}/stup07
#ln -sf ${ARCION}/ipe${CDATE0}/plasma08 ${ARCION}/ipe${CDATE0}/stup08
#ln -sf ${ARCION}/ipe${CDATE0}/plasma09 ${ARCION}/ipe${CDATE0}/stup09
#ln -sf ${ARCION}/ipe${CDATE0}/plasma10 ${ARCION}/ipe${CDATE0}/stup10
#ln -sf ${ARCION}/ipe${CDATE0}/plasma11 ${ARCION}/ipe${CDATE0}/stup11
#ln -sf ${ARCION}/ipe${CDATE0}/plasma12 ${ARCION}/ipe${CDATE0}/stup12
#ln -sf ${ARCION}/ipe${CDATE0}/plasma13 ${ARCION}/ipe${CDATE0}/stup13
#ln -sf ${ARCION}/ipe${CDATE0}/plasma14 ${ARCION}/ipe${CDATE0}/stup14
#ln -sf ${ARCION}/ipe${CDATE0}/plasma15 ${ARCION}/ipe${CDATE0}/stup15
#ln -sf ${ARCION}/ipe${CDATE0}/plasma16 ${ARCION}/ipe${CDATE0}/stup16

#ln -sf ${ARCION}/ipe${CDATE0}/ut_rec   ${ARCION}/ipe${CDATE0}/stup_ut_rec


#copy ionosphere files to ROTDIR/COMOUT to prepare for hpss archiving

/bin/cp -f -p ${DATA}/plasma00  ${COMOUT}/plasma00.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma01  ${COMOUT}/plasma01.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma02  ${COMOUT}/plasma02.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma03  ${COMOUT}/plasma03.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma04  ${COMOUT}/plasma04.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma05  ${COMOUT}/plasma05.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma06  ${COMOUT}/plasma06.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma07  ${COMOUT}/plasma07.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma08  ${COMOUT}/plasma08.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma09  ${COMOUT}/plasma09.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma10  ${COMOUT}/plasma10.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma11  ${COMOUT}/plasma11.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma12  ${COMOUT}/plasma12.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma13  ${COMOUT}/plasma13.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma14  ${COMOUT}/plasma14.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma15  ${COMOUT}/plasma15.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/plasma16  ${COMOUT}/plasma16.${CDUMP}.${CDATE0}

/bin/cp -f -p ${DATA}/FLIP_ERR  ${COMOUT}/FLIP_ERR.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/IPE.inp   ${COMOUT}/IPE.inp.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/ut_rec    ${COMOUT}/ut_rec.${CDUMP}.${CDATE0}

/bin/cp -f -p ${DATA}/logfile_input_params.log  ${COMOUT}/logfile_input_params.log.${CDUMP}.${CDATE0}

/bin/cp -f -p ${DATA}/ipe_grid_neut_N2_den    ${COMOUT}/ipe_grid_neut_N2_den.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/ipe_grid_neut_O2_den    ${COMOUT}/ipe_grid_neut_O2_den.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/ipe_grid_neut_O_den     ${COMOUT}/ipe_grid_neut_O_den.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/ipe_grid_neut_temp      ${COMOUT}/ipe_grid_neut_temp.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/ipe_grid_neut_wind      ${COMOUT}/ipe_grid_neut_wind.${CDUMP}.${CDATE0}
/bin/cp -f -p ${DATA}/ipe_grid_neut_params_ut ${COMOUT}/ipe_grid_neut_params_ut.${CDUMP}.${CDATE0}



fi



export ERR=$?
export err=$ERR
$ERRSCRIPT||exit 2


rm -f NULL
rm -f fort.11 fort.12 fort.14
rm -f fort.15 fort.24 fort.28 fort.29 fort.48
rm -f orography
rm -f orography_uf
rm -f fort.51 fort.52 fort.53
#rm -f SIG.F* SFC.F* FLX.F* LOG.F* D3D.F AER.F*
#rm -f sigr1 sigr2 sfcr nstr
################################################################################
#  Postprocessing
cd $pwd
[[ $mkdata = YES ]]&&rmdir $DATA
$ENDSCRIPT
set +x
if [[ "$VERBOSE" = "YES" ]] ; then
   echo $(date) EXITING $0 with return code $err >&2
fi
exit $err
