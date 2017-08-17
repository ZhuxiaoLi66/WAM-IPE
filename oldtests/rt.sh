#!/bin/bash
set -eu

hostname

die() { echo "$@" >&2; exit 1; }
usage() {
  set +x
  echo
  echo "Usage: $0 -c <model> | -f | -s | -l <file> | -m | -h"
  echo
  echo "  -c  create new baseline results for <model>"
  echo "  -f  run full suite of regression tests"
  echo "  -s  run standard suite of regression tests"
  echo "  -l  runs test specified in <file>"
  echo "  -m  compare against new baseline results"
  echo "  -h  display this help"
  echo "  -F  run on an unsupported platform"
  echo "Only Phase 1 Tide and Gyre, and Theia are supported."
  echo "You must use the -F option to run elsewhere."
  echo
  exit 1
}

[[ $# -eq 0 ]] && usage

source detect_machine.sh

export dprefix=""
export MACHINE_ID=${MACHINE_ID:-wcoss}
if [ $MACHINE_ID = wcoss ]; then
  source /usrx/local/Modules/default/init/sh
  export RTPWD=/nems/noscrub/emc.nemspara/RT/NEMSGSM/trunk-2017011116/data
  export RTBAS=/nems/noscrub/emc.nemspara/RT/NEMSGSM/trunk-2017011116/template
  # export pex=1           # for wcoss phase1
  # export pex=${pex:-2}   # default - phase2
  export QUEUE=debug # dev
  export ACCNR=dev
  if [ $pex -eq 2 ] ; then
   export QUEUE=debug$pex
   export ACCNR=dev$pex
  fi
# export STMP=/stmp$pex
  export STMP=/ptmpp$pex
  export PTMP=/ptmpp$pex
  export SCHEDULER=lsf
  export SIGHDR=/global/save/Shrinivas.Moorthi/para/sorc/global_sighdr.fd/global_sighdr
  cp gfs_fcst_run.IN_IBM gfs_fcst_run.IN
  cp gfs_bsub.IN_wcoss gfs_bsub.IN
elif [ $MACHINE_ID = gaea ]; then
  export RTPWD=/lustre/f1/unswept/ncep/Ratko.Vasic/wx20rv/REGRESSION_TEST    
  export RTBAS=/lustre/f1/unswept/ncep/Ratko.Vasic/wx20rv/REGRESSION_TEST_baselines
  export STMP=/lustre/f1/ncep
  export PTMP=/lustre/f1/ncep
  export SCHEDULER=moab
  cp gfs_fcst_run.IN_Linux gfs_fcst_run.IN
elif [ $MACHINE_ID = theia ]; then
  source /apps/lmod/lmod/init/sh
  export ACCNR
  export dprefix=/scratch4/NCEPDEV
  export RTPWD=/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/RT/WAM-IPE/uncapped-2017011116/data
  export RTBAS=/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/RT/WAM-IPE/uncapped-2017011116/template
  export STMP=$dprefix/stmp4
  export PTMP=$dprefix/stmp3
  export SCHEDULER=pbs
  export MPIEXEC=mpirun
  export SIGHDR=$dprefix/global/save/Shrinivas.Moorthi/para/sorc/global_sighdr.fd/global_sighdr
  export SLG=.false.
  cp gfs_fcst_run.IN_Linux gfs_fcst_run.IN
elif [ $MACHINE_ID = yellowstone ]; then
  export ACCNR=P35071400
  # export ACCNR=UCUB0024
  export QUEUE=small
  export STMP=/glade/scratch
  export PTMP=/glade/scratch
  export SCHEDULER=lsf
  export MPIEXEC=mpirun
  export SIGHDR=/glade/p/work/theurich/para/global_sighdr.fd/global_sighdr
  export SLG=.false.
  cp gfs_fcst_run.IN_Linux gfs_fcst_run.IN
  cp gfs_bsub.IN_yellowstone gfs_bsub.IN
else
  die "Unknown machine ID, please edit detect_machine.sh file"
fi
export pex=${pex:-""}


################################################################
# RTPWD - Path to previously stored regression test answers
# RTBAS - Path to input directory to use for baseline generation
################################################################




export CREATE_BASELINE=false
CB_arg=''
TESTS_FILE='rt.conf'
SET_ID='standard'
ALLOW_UNSUPPORTED=false
while getopts ":c:fsl:mhF" opt; do
  case $opt in
    F)
      ALLOW_UNSUPPORTED=true
      ;;
    c)
      export CREATE_BASELINE=true
      CB_arg=$OPTARG
      SET_ID=' '
      ;;
    f)
      SET_ID=' '
      ;;
    s)
      SET_ID='standard'
      ;;
    l)
      TESTS_FILE=$OPTARG
      SET_ID=' '
      ;;
    m)
      export RTPWD=${STMP}/${USER}/REGRESSION_TEST
      ;;
    h)
      usage
      ;;
    \?)
      usage
      die "Invalid option: -$OPTARG"
      ;;
    :)
      usage
      die "Option -$OPTARG requires an argument."
      ;;
  esac
done

###################################
# PATHRT - Path to regression test
###################################

export PATHRT=$( pwd -P ) # Path to oldtests with symlinks removed (bash-specific)
# Path to nems trunk:
export PATHTR=$( dirname ${PATHRT} )
if [[ -s "$PATHTR/NEMS/NEMSAppBuilder" ]] ; then
    export PATHTR=$PATHTR/NEMS
fi

###################################
# Verify regtest data
###################################

if [[ "$CREATE_BASELINE" == false ]] ; then
    baseline_fingerprint="$RTBAS/REGTEST-FINGERPRINT.md"
else
    baseline_fingerprint="$RTPWD/REGTEST-FINGERPRINT.md"
fi

if ( ! cmp "$baseline_fingerprint" "$PATHRT/REGTEST-FINGERPRINT.md" ) ; then
    echo "Baseline fingerprint does not match repo fingerprint" 1>&2
    echo "Baseline fingerprint: $RTBAS/REGTEST-FINGERPRINT.md" 1>&2
    echo "Repo fingerprint: $PATHRT/REGTEST-FINGERPRINT.md" 1>&2
    echo "You may have the wrong data directory." 1>&2
    exit 1
else
    echo "Regression test fingerprint matches."
    echo "Baseline fingerprint: $RTBAS/REGTEST-FINGERPRINT.md"
    echo "Repo fingerprint: $PATHRT/REGTEST-FINGERPRINT.md"
    echo "Rejoice."
fi

if [[ ( "$MACHINE_ID" != wcoss && "$MACHINE_ID" != theia ) \
      || "$pex" == 2 ]] ; then
    echo "Warning: using unsupported system $MACHINE_ID $pex" 1>&2
    if [[ "$ALLOW_UNSUPPORTED" != true ]] ; then
        echo "Aborting.  Use the -F option (F = \"Force\") to run on this system." 1>&2
        exit 1
    fi
fi

shift $((OPTIND-1))
[[ $# -gt 0 ]] && usage

mkdir -p $STMP/$USER
mkdir -p $PTMP/$USER

if [[ $CREATE_BASELINE == true ]]; then
  #
  # prepare new regression test directory
  #
  export RTPWD_U=${STMP}/${USER}/REGRESSION_TEST
  rm -rf ${RTPWD_U}
  echo "copy REGRESSION_TEST_baselines"
  mkdir -p ${STMP}/${USER}
  cp -r "${RTBAS:?}" "${RTPWD_U:?}"

  if [[ $CB_arg != gfs ]]; then
    echo "copy gfs"
    cp ${RTPWD}/GFS_EULERIAN/*             ${RTPWD_U}/GFS_EULERIAN/.
    cp ${RTPWD}/WAM_gh_l150/*              ${RTPWD_U}/WAM_gh_l150/.
    cp ${RTPWD}/WAM_gh_l150_nemsio/*       ${RTPWD_U}/WAM_gh_l150_nemsio/.
    cp ${RTPWD}/GFS_GOCART_NEMSIO/*        ${RTPWD_U}/GFS_GOCART_NEMSIO/.
    cp ${RTPWD}/GFS_SLG_ADIABATIC/*        ${RTPWD_U}/GFS_SLG_ADIABATIC/.
    cp ${RTPWD}/GFS_SLG/*                  ${RTPWD_U}/GFS_SLG/.
    cp ${RTPWD}/GFS_SLG_RSTHST/*           ${RTPWD_U}/GFS_SLG_RSTHST/.
    cp ${RTPWD}/GFS_SLG_48PE/*             ${RTPWD_U}/GFS_SLG_48PE/.
    cp ${RTPWD}/GFS_SLG_T574/*             ${RTPWD_U}/GFS_SLG_T574/.
    cp ${RTPWD}/GFS_SLG_NSST/*             ${RTPWD_U}/GFS_SLG_NSST/.
    cp ${RTPWD}/GFS_SLG_STOCHY/*           ${RTPWD_U}/GFS_SLG_STOCHY/.
    cp ${RTPWD}/GFS_SLG_LAND/*             ${RTPWD_U}/GFS_SLG_LAND/.
  fi
  if [[ $CB_arg != post ]]; then
    echo "copy post"
#    cp -r ${RTPWD}/GFS_GOCART_POST/*       ${RTPWD_U}/GFS_GOCART_POST/.
  fi
fi

export REGRESSIONTEST_LOG=${PATHRT}/RegressionTests_$MACHINE_ID.log
COMPILE_LOG=${PATHRT}/Compile_$MACHINE_ID.log

date > "$COMPILE_LOG"
echo "Start Regression test" >> "$COMPILE_LOG"

date > ${REGRESSIONTEST_LOG}
echo "Start Regression test" >> ${REGRESSIONTEST_LOG}
(echo;echo;echo)             >> ${REGRESSIONTEST_LOG}

export RUNDIR_ROOT=${PTMP}/${USER}/rt_$$
mkdir -p ${RUNDIR_ROOT}

source default_vars.sh

export TEST_NR=0
rm -f fail_test
export TEST_NAME
cat $TESTS_FILE | while read line; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line == \#* ]] && continue

  if [[ $line == APPBUILD* ]] ; then
      APPBUILD_OPTS=`echo $line | cut -d'|' -f2`
      SET=`     echo $line | cut -d'|' -f3`
      MACHINES=`echo $line | cut -d'|' -f4`
      [[ $SET_ID != ' ' && $SET != *${SET_ID}* ]] && continue
      [[ $MACHINES != ' ' && $MACHINES != *${MACHINE_ID}* ]] && continue

      echo "Removing old exe/NEMS.x"
      echo "Removing old exe/NEMS.x" >> "$COMPILE_LOG"
      echo "Removing old exe/NEMS.x" >> "$REGRESSIONTEST_LOG"
      rm -f "$PATHTR/exe/NEMS.x" "$PATHTR/src/conf/modules.nems"

      echo "NEMSAppBuilder $APPBUILD_OPTS"
      echo "NEMSAppBuilder $APPBUILD_OPTS" >> "$COMPILE_LOG"
      echo "NEMSAppBuilder $APPBUILD_OPTS" >> "$REGRESSIONTEST_LOG"
      cd $PATHTR/..
      pwd -P
      pwd
      ./NEMS/NEMSAppBuilder rebuild $APPBUILD_OPTS >> "$COMPILE_LOG" 2>&1
      cd $PATHRT
      if [[ ! -s "$PATHTR/exe/NEMS.x" || ! -x "$PATHTR/exe/NEMS.x" || \
            ! -s "$PATHTR/src/conf/modules.nems" ]] ; then
          echo "NEMSAppBuilder failed.  Abort." >> "$REGRESSIONTEST_LOG"
          echo "No executable: $PATHTR/exe/NEMS.x" >> "$REGRESSIONTEST_LOG"
          echo "For details, look in \"$COMPILE_LOG\"" >> "$REGRESSIONTEST_LOG"
          exit 1
      fi
    continue
  elif [[ $line == COMPILE* ]] ; then
      NEMS_VER=`echo $line | cut -d'|' -f2`
      SET=`     echo $line | cut -d'|' -f3`
      MACHINES=`echo $line | cut -d'|' -f4`
      ESMF_VER=`echo $line | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//'`
      [[ $SET_ID != ' ' && $SET != *${SET_ID}* ]] && continue
      [[ $MACHINES != ' ' && $MACHINES != *${MACHINE_ID}* ]] && continue

      echo "Removing old exe/NEMS.x"
      echo "Removing old exe/NEMS.x" >> "$COMPILE_LOG"
      echo "Removing old exe/NEMS.x" >> "$REGRESSIONTEST_LOG"
      rm -f "$PATHTR/exe/NEMS.x" "$PATHTR/src/conf/modules.nems"

      echo "Compiling $NEMS_VER $ESMF_VER"
      echo "Compiling $NEMS_VER $ESMF_VER" >> "$COMPILE_LOG"
      echo "Compiling $NEMS_VER $ESMF_VER" >> "$REGRESSIONTEST_LOG"
      cd $PATHTR/src
      ./configure ${ESMF_VER}_${MACHINE_ID} >> "$COMPILE_LOG" 2>&1
      if [[ ! -s "$PATHTR/src/conf/modules.nems" ]] ; then
          echo "Configure failed.  Abort." >> "$REGRESSIONTEST_LOG"
          echo "File missing: $PATHTR/src/conf/modules.nems" >> "$REGRESSIONTEST_LOG"
          echo "For details, look in \"$COMPILE_LOG\"" >> "$REGRESSIONTEST_LOG"
          echo "Configuring ${ESMF_VER}_${MACHINE_ID}" > fail_test
          exit 1
      fi
      (
          set -e
          source ./modules.nems.sh
          module list
#          source conf/modules.nems              >> "$COMPILE_LOG" 2>&1
          module list                           >> "$COMPILE_LOG" 2>&1
          gmake clean                           >> "$COMPILE_LOG" 2>&1
          gmake ${NEMS_VER} J=-j2               >> "$COMPILE_LOG" 2>&1
      )
      cd $PATHRT
      if [[ ! -s "$PATHTR/exe/NEMS.x" || ! -x "$PATHTR/exe/NEMS.x" ]] ; then
          echo "Compilation failed.  Abort." >> "$REGRESSIONTEST_LOG"
          echo "No executable: $PATHTR/exe/NEMS.x" >> "$REGRESSIONTEST_LOG"
          echo "For details, look in \"$COMPILE_LOG\"" >> "$REGRESSIONTEST_LOG"
          echo "Compiling $NEMS_VER $ESMF_VER" > fail_test
          exit 1
      fi
    continue
  elif [[ $line == RUN* ]] ; then
    cd "$PATHRT"
    TEST_NAME=`echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//'`
    SET=`      echo $line | cut -d'|' -f3`
    MACHINES=` echo $line | cut -d'|' -f4`
    CB=`       echo $line | cut -d'|' -f5`
    [[ -e "tests/$TEST_NAME" ]] || die "run test file tests/$TEST_NAME does not exist"
    [[ $SET_ID != ' ' && $SET != *${SET_ID}* ]] && continue
    [[ $MACHINES != ' ' && $MACHINES != *${MACHINE_ID}* ]] && continue
    [[ $CREATE_BASELINE == true && $CB != *${CB_arg}* && 'all' != *${CB_arg}* ]] && continue

    (( TEST_NR += 1 ))
    (
      export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}
      source tests/$TEST_NAME
      export JBNME=`basename $RUNDIR_ROOT`_${TEST_NR}
      echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}" >> ${REGRESSIONTEST_LOG}
      echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}"
      ./${RUN_SCRIPT} || die "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR} failed"
    )
    continue
  else
    die "Unknown command $line"
  fi
done

if [ -e fail_test ]; then
  echo "FAILED TESTS: "
  echo "FAILED TESTS: " >> ${REGRESSIONTEST_LOG}
  for failed_test_name in `cat fail_test`
  do
    echo "Test " ${failed_test_name} " failed "
    echo "Test " ${failed_test_name} " failed " >> ${REGRESSIONTEST_LOG}
  done
  echo ; echo REGRESSION TEST FAILED
  (echo ; echo REGRESSION TEST FAILED) >> ${REGRESSIONTEST_LOG}
else
  echo ; echo REGRESSION TEST WAS SUCCESSFUL
  (echo ; echo REGRESSION TEST WAS SUCCESSFUL) >> ${REGRESSIONTEST_LOG}
fi

# Finalize, Clenaup
rm -f err out gfs_fcst_run \
nems.configure gfs_qsub gfs_fcst_run.IN ngac_qsub ngac_bsub gfs_bsub \
configure_file_01 configure_file_02 configure_file_03 configure_file_04 \
atmos.configure fail_test
# atmos.configure fail_test

date >> ${REGRESSIONTEST_LOG}

exit
