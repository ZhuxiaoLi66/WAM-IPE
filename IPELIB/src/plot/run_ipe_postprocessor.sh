#!/bin/bash

# cd to the postprocessor code directory #

cd /scratch3/NCEPDEV/swpc/noscrub/George.Millward/convert_ipe 

# compile the postprocessor

module load intel
ifort -o postprocess_ipe interface_to_fixed_height_single_file.f90

# define the results directory and the names of input and output files

export RESULTS_DIRECTORY="/scratch3/NCEPDEV/swpc/noscrub/George.Millward/ipe_results/" 
export INPUT_PLASMA_FILE="ipe_grid_plasma_params.20170716T1954"
export OUTPUT_PLASMA_FILE="ipe_TEC_nmf2_hmf2.20170716T1954"

# run the postprocessor

./postprocess_ipe 

exit 0
