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
#
HOME_DIR=/home/Naomi.Maruyama
PROJECT_DIR="$HOME_DIR"/ipe
###RUN_DIR="$PROJECT_DIR"/run
#---storage directory
STR_DIR=/scratch2/portfolios/BMC/idea/Naomi.Maruyama
#---run directory
RUN_DIR=/scratch2/portfolios/NCEPDEV/ptmp/Naomi.Maruyama
RUN_DATE=20120402
TEST0=3d.trans
RUNID=v57
###NEW_RUN_DIR="$RUN_DATE"."$TEST0"."$RUNID"
NEW_RUN_DIR="$RUNID"
###for ipe_grid
#0. APEX original
#1. APEX CORRECTED
#2. APEX CORRECTED new Q
#3. tilted dipole
#4. tilted dipole : new Q
sw_new_dir="T"
sw_grid="4"
sw_nmltinp="T"
sw_startup="T"
sw_rscrpt="T"
sw_efield="T"
###for IPE.inp
if [ "$sw_nmltinp" = "T" ]; then
###    RUN_DIR1="$RUNID0"/ut"$UT0"
###    RUN_DIR1="$RUNID0"/
    RUN_DIR1=v33.bk20120313
fi
#
###for startup
if [ "$sw_startup" = "T" ]; then
    RUNID0=v33.jet
    UT0=248306
    RUN_DIR2="$RUNID0"/but"$UT0"
###    RUN_DIR2="$RUNID0"
    sw_startup_special="F" 
fi
#
if [ "$sw_new_dir" = "T" ]; then
     mkdir "$RUN_DIR"/"$NEW_RUN_DIR"
fi
cd    "$RUN_DIR"/"$NEW_RUN_DIR"
pwd
#
if [ "$sw_nmltinp" = "T" ]; then
    cp "$STR_DIR"/"$RUN_DIR1"/IPE.inp  "$PWD"/
fi
#
if [ "$sw_startup" = "T" ]; then
###i should use the loop!!!
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma00  ./stup00
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma01  ./stup01
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma02  ./stup02
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma03  ./stup03
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma04  ./stup04
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma05  ./stup05
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma06  ./stup06
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma07  ./stup07
    ln -s "$STR_DIR"/"$RUN_DIR2"/plasma08  ./stup08
    ln -s "$STR_DIR"/"$RUN_DIR2"/ut_rec    ./stup_ut_rec

  if [ "$sw_startup_special" = "F" ]; then
      ln -s "$STR_DIR"/"$RUN_DIR2"/plasma09  ./stup09
      ln -s "$STR_DIR"/"$RUN_DIR2"/plasma10  ./stup10
      ln -s "$STR_DIR"/"$RUN_DIR2"/plasma11  ./stup11
  else
      ln -s "$STR_DIR"/"$RUN_DIR2"/plasma13  ./stup09
      ln -s "$STR_DIR"/"$RUN_DIR2"/plasma14  ./stup10
      ln -s "$STR_DIR"/"$RUN_DIR2"/plasma15  ./stup11
  fi

fi
GRID_DR="$STR_DIR"/field_line_grid
fnAPEX0="$GRID_DR"/20110419lowres_global/GIP_apex_coords_global_lowres
fnAPEX1="$fnAPEX0"CORRECTED20120112
fntilt="$GRID_DR"/20111025dipole/GIP_tilt_coords_global_lowresCORRECTED20120112
FN0=ipe_grid
# APEX
if [ "$sw_grid" = "0" ]; then
    ln -s "$fnAPEX0"  ./"$FN0"

elif [ "$sw_grid" = "1" ]; then
    ln -s $fnAPEX1  ./"$FN0"

elif [ "$sw_grid" = "2" ]; then
# APEX new Q
    ln -s "$fnAPEX1".nQ  ./"$FN0"

elif [ "$sw_grid" = "3" ]; then
# tilted dipole 
    ln -s "$fntilt"  ./"$FN0"

elif [ "$sw_grid" = "4" ]; then
# tilted dipole : new Q
    ln -s "$fntilt".nQ  ./"$FN0"
else
    echo  "!STOP! INVALID sw_grid"
    echo  "$sw_grid"
fi
#
if [ "$sw_rscrpt" = "T" ]; then
#	RUN_DIR4="$RUNID0"
	RUN_DIR4="$RUN_DIR1"
	FN0=r.ipe
	ls -tl  "$STR_DIR"/"$RUN_DIR4"/"$FN0"
	cp  "$STR_DIR"/"$RUN_DIR4"/"$FN0"  ./
#FN="$FN0"."$RUNID"
####note: i do not know how to readin the file name into a variable
#    cat filename.log
#    FN=`cat filename.log`
#    echo 'FN' "$FN"
#####modify the file name from v45 to v48
#####modify what's in the file from v45 to v48: cat ...| sed    
fi
#
if [ "$sw_efield" = "T" ]; then

    RUN_DIR3="$HOME_DIR"/sandbox/efield
    ln -s "$RUN_DIR3"/coeff_lflux.dat  ./
    ln -s "$RUN_DIR3"/coeff_hflux.dat  ./
    ln -s "$RUN_DIR3"/wei96.cofcnts    ./
fi
