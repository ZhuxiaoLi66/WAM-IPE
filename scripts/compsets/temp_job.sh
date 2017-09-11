#!/bin/bash
#BSUB -J newscripts
#BSUB -P SPACE-T2O
#BSUB -W 30
#BSUB -o /ptmpp2/Adam.Kubaryk/newscripts/fcst.%J
#BSUB -e /ptmpp2/Adam.Kubaryk/newscripts/fcst.%J
#BSUB -n 40
#BSUB -q debug2
#BSUB -a poe
#BSUB -x
#BSUB -network type=sn_all:mode=US
#BSUB -R span[ptile=24]
set -ax

cd /global/save/Adam.Kubaryk/nems/WAM-IPE/scripts/compsets

##-------------------------------------------------------
## source config file
##-------------------------------------------------------

. fcst_config

##-------------------------------------------------------
## execute forecast
##-------------------------------------------------------

mkdir -p /stmpp2/Adam.Kubaryk/newscripts
cd /stmpp2/Adam.Kubaryk/newscripts

export VERBOSE=YES

. /global/save/Adam.Kubaryk/nems/WAM-IPE/scripts/compsets/exglobal_fcst_nems.sh
if [ 0 != 0 ]; then echo "forecast failed, exit"; exit; fi
echo "fcst done"

exit 
