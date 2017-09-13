#!/bin/bash -l
#
#PBS -N ipe_convert
#PBS -A swpc
#PBS -q bigmem
#PBS -j oe
#PBS -l procs=16
#PBS -l walltime=00:06:00
#PBS -W umask=027

set -aeux

### setup
## IPE_CONV is the directory containing ipe_conv_slice.py
cd $PBS_O_WORKDIR
IPE_CONV=`pwd`
#IPE_CONV=/scratch3/NCEPDEV/swpc/save/Adam.Kubaryk/util/ipe_conv

## INPUT_DIR should contain plasma*, ipe_grid_neut_* files, ut_rec, and IPE.inp
INPUT_DIR=/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/rt_123253/swpc%20130316_1hr_newspacewx_gsm%wam%T62_ipe%80x170
## IPE.inp is used to get the timestep currently... if it cannot find it, the program will default to a 30-minute output timestep

## note, STUP="-s" does not work because stup13-15 don't exist?
# if INPUT_DIR has stup* instead of plasma*, uncomment the following line:
#STUP="-s"

## TIMESTAMP format YYYYMMDDHHmm corresponding to the initial date, i.e. first read/write will be at (TIMESTAMP + DELTIM)
TIMESTAMP=201303160000

## OUTPUT_DIR will contain your ipe_grid_*_params.T* files
OUTPUT_DIR=/scratch4/NCEPDEV/stmp4/$USER/ipe_converted_output

### initialize
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

### modules
. /apps/lmod/lmod/init/bash
module purge
module use -a /contrib/modulefiles
module load intel anaconda/anaconda3-4.4.0

### run
PYTHON=${PYTHON:-"python"}
IPE_CONV_EXEC=${IPE_CONV_EXEC:-"ipe_conv_slice.py"}
STUP=${STUP:-""}
$PYTHON $IPE_CONV/$IPE_CONV_EXEC -i $INPUT_DIR -o $OUTPUT_DIR -t $TIMESTAMP $STUP
