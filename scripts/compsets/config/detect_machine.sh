#!/bin/bash
### basically the idea is to just check if directories exist that are specific to certain machines

if [[ -e /gpfs/hps && -e /usrx ]] ; then
	export MACHINE=wcoss
	if [[ -d /global ]] ; then
		# We are on WCOSS Phase 1 or 2.
		if ( ! cat /proc/cpuinfo | grep 'processor.*32' > /dev/null ) ; then
			# Fewer than 32 fake (hyperthreading) cpus, so Phase 1.
			export FMID=wcoss.phase1
			export PEX=1
		else
			export FMID=wcoss.phase2
			export PEX=2
		fi
	else
		# WCOSS Cray
		export FMID=wcoss.cray
		export PEX=c
	fi
elif [[ -e /scratch4 && -e /scratch3 ]] ; then
	export FMID=theia
	export MACHINE=theia
elif [[ -e /glade ]] ; then
	export FMID=yellowstone
	export MACHINE=yellowstone
elif [[ -e /pan2 && -e /lfs3 ]] ; then
	export FMID=jet
	export MACHINE=jet
elif ( hostname | grep -i gaea ) ; then
	export FMID=gaea
	export MACHINE=gaea
elif [[ -e /nobackupp8 ]] ; then
	export FMID=pleiades
	export MACHINE=pleiades
else
	# cannot ID machine
	echo "cannot identify current machine... check config/detect_machine"
	exit 1
fi

# now check if config file exists
if [ ! -e $CONFIGDIR/$FMID.config ] ; then
	echo "machine configuration file not found for FMID: $FMID"
	exit 2
fi
