#!/bin/csh

setenv GRID /scratch3/NCEPDEV/swpc/noscrub/Joseph.Schoonover/ipe_data/GIP_Fixed_GEO_grid_lowres_corrected.bin

setenv I2HG_EXE /scratch3/NCEPDEV/swpc/noscrub/Joseph.Schoonover/Software/WAM-IPE/IPELIB/bin/i2hg

setenv IPE_RUNDIR ipe_theia_intel_parallel_16/



# Do not modify below this point unless you really know what you are doing
${I2HG_EXE} -i ${IPE_RUNDIR}ipe_grid_plasma_params.iter_00432000 -o ${IPE_RUNDIR} -g ${GRID}
