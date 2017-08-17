#! /bin/bash
#! DATE: 08 September, 2011
#!********************************************
#!***      Copyright 2011 NAOMI MARUYAMA   ***
#!***      ALL RIGHTS RESERVED             ***
#!********************************************
#! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
#! DEVELOPER: Dr. Naomi Maruyama
#! CONTACT INFORMATION:
#! E-MAIL : Naomi.Maruyama@noaa.gov
#! PHONE  : 303-497-4857
#! ADDRESS: 325 Broadway, Boulder, CO 80305
#!--------------------------------------------  
# cmprs2hpss.bash
#purpose of this script:
#(1) create tgz file
#(2) put the tgz file onto the HPSS system 
#NOTE: two steps are required because htar does not have the compress option...
HOME_DIR=/home/Naomi.Maruyama
#PROJECT_DIR="$HOME_DIR"/ipe
RUN_DIR=/scratch2/portfolios/BMC/idea/Naomi.Maruyama
#---which dir you want to create tgz?
RUNID=v52.bkup20120301dbg
#---when was this run performed?
RUN_DATE=20120229
#---temporary working dir where the tgz file is created before putting to hpss
WORK_DIR=/scratch2/portfolios/NCEPDEV/ptmp/Naomi.Maruyama/tmp
#---which dir you want to put on hpss?
HPSS_DIR=/BMC/idea/Naomi.Maruyama/ipe/runs
#
sw_cmprs="T"
sw_hpss="T"
sw_del="T"
#
cd "$RUN_DIR"
pwd
flnm="$RUNID"."$RUN_DATE".tgz
if [ "$sw_cmprs" = "T" ]; then
    echo "compressing"
    echo "$RUNID"
    tar cvzf "$WORK_DIR"/"$flnm"  ./"$RUNID"
    echo "compressing completed"
fi
#
if [ "$sw_hpss" = "T" ]; then
    cd  "$WORK_DIR"
    pwd
    module load hpss
    echo "putting to hpss dir:"
    echo "$HPSS_DIR"
#hsi put /home/Naomi.Maruyama/ptmp/tmp/v52.20120310.tgz : /BMC/idea/Naomi.Maruyama/ipe/runs/v52.20120310.tgz
    hsi put "$WORK_DIR"/"$flnm" : "$HPSS_DIR"/"$flnm"
    hsi ls -tl "$HPSS_DIR"
    ls -tl 
    echo "putting to hpss completed"
    rm  "$flnm"
fi
#
if [ "$sw_del" = "T" ]; then
    cd "$RUN_DIR"
    pwd
    echo "start deleting the dir"
    echo "$RUN_DIR"
#interactively delete the files
    mv "$RUNID" "$WORK_DIR"/
fi