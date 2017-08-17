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
### mote_output
### ---PURPOSE---move output files on the run directory before re=running the new run
mv  ../plasma0[0-9]*  ./
mv  ../plasma1[0-9]*  ./
mv  ../fort.*            ./
mv  ../logfile*.log  ./
mv  ../FLIP_ERR  ./
mv  ../ut_rec ./
mv  ../IPE.inp  ./
mv  ../error ./
mv  ../output ./
mv  ../ipe*.log ./
#mv  ../STDIN.*  ./
#mv  ../gmon.out  ./
