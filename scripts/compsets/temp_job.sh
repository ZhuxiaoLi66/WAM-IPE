#!/bin/bash
#PBS -N newscripts
#PBS -A swpc
#PBS -l walltime=00:30:00
#PBS -o /scratch4/NCEPDEV/stmp3/Adam.Kubaryk/newscripts/fcst.%J
#PBS -e /scratch4/NCEPDEV/stmp3/Adam.Kubaryk/newscripts/fcst.%J
#PBS -l procs=40
#PBS -q debug
#PBS -W umask=027
set -ax

cd /scratch3/NCEPDEV/swpc/save/Adam.Kubaryk/nems/mdate_compatible_timestamp/WAM-IPE/scripts/compsets

##-------------------------------------------------------
## source config file
##-------------------------------------------------------

. /scratch3/NCEPDEV/swpc/save/Adam.Kubaryk/nems/mdate_compatible_timestamp/WAM-IPE/scripts/compsets/config/workflow.sh /scratch3/NCEPDEV/swpc/save/Adam.Kubaryk/nems/mdate_compatible_timestamp/WAM-IPE/scripts/compsets/user.config

##-------------------------------------------------------
## execute forecast
##-------------------------------------------------------

mkdir -p /scratch4/NCEPDEV/stmp4/Adam.Kubaryk/newscripts
cd /scratch4/NCEPDEV/stmp4/Adam.Kubaryk/newscripts

export VERBOSE=YES

. /scratch3/NCEPDEV/swpc/save/Adam.Kubaryk/nems/mdate_compatible_timestamp/WAM-IPE/scripts/compsets/exglobal/exglobal_fcst_nems.sh
if [ 0 != 0 ]; then echo "forecast failed, exit"; exit; fi
echo "fcst done"

exit 
