#!/bin/bash

export SCRIPTSDIR=`pwd`
export CONFIGDIR=`pwd`/config

# run detect_machine, will exit if machine is unknown or config file doesn't exist for machine
. $CONFIGDIR/detect_machine.sh

# load user setup file
. $1

# load all machine configuration logic
. $CONFIGDIR/$FMID.config

# load computational logic
. $CONFIGDIR/compute.config

# load non-machine-specific general configuration: dependent directories, executables, etc.
. $CONFIGDIR/general.config

if [ $RESTART = .true. ] ; then
. $CONFIGDIR/restart.config
else
. $CONFIGDIR/coldstart.config
fi

# load ESMF variables
. $CONFIGDIR/esmf.config

# load WAM-specific configuration
. $CONFIGDIR/wam.config

if [ $WAM_IPE_COUPLING = .true. ] ; then
# load IPE-specific configuration
. $CONFIGDIR/ipe.config

# load the coupled configuration
. $CONFIGDIR/coupled.config

# load the namelist options
if [ $NAMELIST_VER = 'compset' ] ; then
. $CONFIGDIR/wam-ipe_ocr_dpnamelist.config
else
. $CONFIGDIR/wam-ipe_dpnamelist.config
fi

else
if [ $NAMELIST_VER = 'compset' ] ; then
. $CONFIGDIR/wam_ocr_dpnamelist.config
else
. $CONFIGDIR/wam_ocr_dpnamelist.config
fi

fi

# run our checks to make sure we're not walking into any walls before we try to run a job
. $CONFIGDIR/checks.sh
