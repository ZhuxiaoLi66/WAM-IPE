#!/bin/bash



setup(){
  HASHID=$(git log | grep commit | head -1 | awk -F " " '{print substr($2,1,8)}')
  echo 'Model Hash ID :' ${HASHID}
  STMPDIR="/scratch4/NCEPDEV/stmp4/${USER}/${HASHID}"
  PTMPDIR="/scratch4/NCEPDEV/stmp3/${USER}/${HASHID}"
}

run_model(){

# Modify the config file so that the JOBNAME corresponds to the hash id
sed -i '/JOBNAME/c\export JOBNAME='${HASHID} ${CONFIG_FILE}
# Submit the job and pipe the output to a temporary file
./submit.sh ${CONFIG_FILE} 1 ${NCYCLES}

}

model_plots(){

  # Copy the IPE.inp file from STMPDIR to PTMPDIR -- this file is needed by the Convert_mpi.batch script
  cp ${STMPDIR}/IPE.inp ${PTMPDIR}
  
  # Set the directory where IPE output is located for Convert_mpi.batch to make plots
  sed -i '/IPE_RUNDIR=/c\IPE_RUNDIR="'${PTMPDIR}/'"' ../../IPELIB/scripts/Convert_mpi.batch
  
  mkdir -p /scratch3/NCEPDEV/swpc/noscrub/wam-ipe_regression-plots/${HASHID}

  CONFIG=$( echo ${CONFIG_FILE} | awk -F "." '{print $1}' )
  
  # Do the polar plots for IPE
  sed -i '/PLOTDIR=/c\PLOTDIR="/scratch3/NCEPDEV/swpc/noscrub/wam-ipe_regression-plots/'${HASHID}'/'${CONFIG}'/ipe/polar_plots/"' ./Convert_mpi.batch
  sed -i '/PLOTTYPE=/c\PLOTTYPE="-p"' ./Convert_mpi.batch
  ./Convert_mpi.batch
  
#  # Do the mercator plots for IPE
#  sed -i '/PLOTDIR=/c\PLOTDIR="/scratch3/NCEPDEV/swpc/noscrub/wam-ipe_regression-plots/'${HASHID}'/'${CONFIG}'/ipe/mercator_plots/"' ./Convert_mpi.batch
#  sed -i '/PLOTTYPE=/c\PLOTTYPE=""' ./Convert_mpi.batch
#  ./Convert_mpi.batch

}

model_runtimes(){

  nProc=$(ls ${STMPDIR}/PET*.ESMF_LogFile | wc -l)
  
  for petfile in $(ls ${STMPDIR}/PET350.ESMF_LogFile)
  do
    datestr=''
    grep sub-update_IPE\ finished ${petfile} | while read -r line ; do

      if [[ -n $datestr ]]; then
        ./timepet.py "${line:0:15}" "$datestr" >> timepet.out
      fi
      datestr=${line:0:15}

    done

  done
}

# ----- Parse through command line options ----- #
HELP="no"
PLOT="no"
POSITIONAL=()

if [ $# -eq 0 ]; then
    HELP="yes"
fi

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r|--regression-config)
    CONFIG_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
    HELP="yes"
    shift # past argument
    shift # past value
    ;;
    -n|--ncycles)
    NCYCLES="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--plot)
    PLOT="yes"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo ${CONFIG_FILE}
if [ -z "$CONFIG_FILE" ]; then
  HELP="yes"
fi

if [ -z "$NCYCLES" ]; then
  NCYCLES=1
fi 

if [ "${HELP}" = "yes" ]; then
  echo 'Usage : ./handle_regressions.sh [options]'
  echo '  Options :'
  echo '    -h   | --help '
  echo '        Display this help message '
  echo ' '
  echo '    -r  <config-file>   | --regression-config <config-file>'
  echo '        Sets which base regression test config file to use'   
  echo '        This option is required'
  echo ' '
  echo '    -n  <number of cycles>   | --ncycles <number of cycles>'
  echo '        Sets the number of forecast cycles to run for the test'   
  echo ' '
  echo '    -p  | --plot'
  echo '        Only does the plotting of model output. Can only be run'
  echo '        after forecast simulations have been completed'   
  echo ' '
  echo ' --------------------------------------------------------------------- '
  echo ' '
  echo ' Suggested Usage : '
  echo '   (1) Run the coupled_20130316 test case for five cycles. This attempts'
  echo '       5 days of model forecasts through 1 day forecasts with restarts. '
  echo ' '
  echo '         ./handle_regressions.sh -r coupled_20130316.config -n 5        '
  echo ' '
  echo '   (2) Once the simulations are complete, you can generate the standard '
  echo '       set of plots for later inspection'
  echo ' '
  echo '         ./handle_regressions.sh -r coupled_20130316.config -p          '
  echo ' '
  echo ' --------------------------------------------------------------------- '

fi

# ----- Parse through command line options ----- #


if [ "${HELP}" = "no" ]; then

  setup
  if [ "${PLOT}" = "no" ]; then
    run_model
  fi

  if [ "${PLOT}" = "yes" ]; then
    model_plots
    model_runtimes
  fi
   
fi


