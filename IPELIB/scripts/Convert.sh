#!/bin/sh


# Set the full path to the restart directory with the grid files can be found.
#setenv RESDIR /scratch3/NCEPDEV/swpc/noscrub/ipe_initial_conditions/2013160300/grid_80x170/
export RESDIR=/scratch3/NCEPDEV/swpc/noscrub/ipe_initial_conditions/2013160300/grid_80x170/


# Set the full path to the i2hg executable
#setenv I2HG_EXE /scratch3/NCEPDEV/swpc/noscrub/Joseph.Schoonover/Software/WAM-IPE/IPELIB/bin/i2hg
export I2HG_EXE=/scratch3/NCEPDEV/swpc/noscrub/Joseph.Schoonover/Software/WAM-IPE/IPELIB/bin/i2hg

# Set the path to the directory where your IPE output can be found
#setenv IPE_RUNDIR ipe_theia_intel_parallel_16/
export IPE_RUNDIR=1504711534_ipe_theia_intel_parallel_32/



# Do not modify below this point unless you really know what you are doing

#Bring in the IPE grid
cp ${RESDIR}ipe_grid ./

for file in ${IPE_RUNDIR}ipe_grid_plasma_params.*
do
  echo $file
  timestamp=$( echo $file | awk -F'.' '{print $2}')
  echo ${IPE_RUNDIR}ipe_grid_neutral_params.$timestamp
  neutralfile=ipe_grid_neutral_params.$timestamp

${I2HG_EXE} --plasma-file $file \
            --output-dir ${IPE_RUNDIR} \
            --neutral-file ${IPE_RUNDIR}${neutralfile}\
            --grid-file ${RESDIR}GIP_Fixed_GEO_grid_lowres_corrected.bin

break
done

