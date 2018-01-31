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

# setup input I/O. note: some values are overwritten if RESTART=.true.
. $CONFIGDIR/coldstart.config

# overwrite the input files and adjust FH* settings for restart if in use
if [ $RESTART = .true. ] ; then
. $CONFIGDIR/restart.config
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
. $CONFIGDIR/wam-ipe_dpnamelist.config


else # standalone
. $CONFIGDIR/wam_dpnamelist.config
fi

# run our checks to make sure we're not walking into any walls before we try to run a job
. $CONFIGDIR/checks.sh
