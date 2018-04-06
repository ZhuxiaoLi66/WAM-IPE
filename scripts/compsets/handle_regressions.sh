#!/bin/bash


export CONFIG_FILE='coupled_20130316.config'
export NCYCLES=3



#


setup(){
  HASHID=$(git log | grep commit | head -1 | awk -F " " '{print substr($2,1,8)}')
  echo 'Model Hash ID :' ${HASHID}
  STMPDIR='/scratch4/NCEPDEV/stmp4/${USER}/${HASHID}'
  PTMPDIR='/scratch4/NCEPDEV/stmp3/${USER}/${HASHID}'
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
  sed -i '/IPE_RUNDIR/c\IPE_RUNDIR="'${PTMPDIR}'"' ../../IPELIB/scripts/Convert_mpi.batch
  
  mkdir -p /scratch3/NCEPDEV/swpc/noscrub/wam-ipe_regression-plots/${HASHID}
  
  # Do the polar plots for IPE
  sed -i '"/PLOTDIR/c\PLOTDIR=/scratch3/NCEPDEV/swpc/noscrub/wam-ipe_regression-plots/'${HASHID}'/polar_plots/"' ../../IPELIB/scripts/Convert_mpi.batch
  sed -i '/PLOTTYPE/c\PLOTTYPE="-p"' ../../IPELIB/scripts/Convert_mpi.batch
  ../../IPELIB/scripts/Convert_mpi.batch
  
  # Do the mercator plots for IPE
  sed -i '"/PLOTDIR/c\PLOTDIR=/scratch3/NCEPDEV/swpc/noscrub/wam-ipe_regression-plots/'${HASHID}'/mercator_plots/"' ../../IPELIB/scripts/Convert_mpi.batch
  sed -i '/PLOTTYPE/c\PLOTTYPE=""' ../../IPELIB/scripts/Convert_mpi.batch
  ../../IPELIB/scripts/Convert_mpi.batch

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
    LIBPATH="$2"
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

if [ "${HELP}" = "yes" ]; then
  echo 'Usage : ./handle_regressions.sh [options]'
  echo '  Options :'
  echo '    -h   | --help '
  echo '        Display this help message '
  echo ' '
  echo '    -r  <config-file>   | --regression-config <config-file>'
  echo '        Sets which base regression test config file to use'   
  echo ' '
  echo '    -n  <number of cycles>   | --ncycles <number of cycles>'
  echo '        Sets the number of forecast cycles to run for the test'   
  echo ' '
  echo '    -p  | --plot'
  echo '        Only does the plotting of model output. Can only be run'
  echo '        after forecast simulations have been completed'   
  echo ' '
fi

# ----- Parse through command line options ----- #


if [ "${HELP}" = "no" ]; then

  setup
  if [ "${PLOT}" = "no" ]; then
    run_model
  fi

  if [ "${PLOT}" = "yes" ]; then
    model_plots
  fi
   
fi


