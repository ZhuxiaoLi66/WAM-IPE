#!/bin/csh


# Set the full path to the restart directory with the grid files can be found.
setenv RESDIR /scratch3/NCEPDEV/swpc/noscrub/ipe_initial_conditions/2013160300/grid_80x170/


# Set the full path to the i2hg executable
setenv I2HG_EXE /scratch3/NCEPDEV/swpc/noscrub/Joseph.Schoonover/Software/WAM-IPE/IPELIB/bin/i2hg

# Set the path to the directory where your IPE output can be found
setenv IPE_RUNDIR ipe_theia_intel_parallel_16/





# Do not modify below this point unless you really know what you are doing


${I2HG_EXE} --plasma-file ${IPE_RUNDIR}ipe_grid_plasma_params.iter_00432000 \
            --output-dir ${IPE_RUNDIR} \
            --grid-file ${RESDIR}GIP_Fixed_GEO_grid_lowres_corrected.bin
