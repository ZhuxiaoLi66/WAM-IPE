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
### delete_output
### ---PURPOSE---delete output files on the run directory before re=running the new run
HOME_DIR=/lfs0/projects/idea/maruyama/sandbox/ipe/run
rm  plasma0[0-9]*  plasma1[0-7]*
rm  fort.*
rm  logfile*.log  FLIP_*  ut_rec*
rm  error output
rm  gmon.out
#######rm  r.ipe.*.e*  r.ipe.*.o* r.ipe.*.pe*  r.ipe.*.po* 
