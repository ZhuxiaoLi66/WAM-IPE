#!/bin/bash

pwd=$(pwd)

export FHMAX=1

## set restart
cycle=${2:-1}
if [[ $cycle == 1 ]] ; then
	export RESTART=.false.
else
	export RESTART=.true.
fi

export HOUR_START=$((${cycle} * ${FHMAX}-${FHMAX}))
echo "HOUR "${HOUR_START}

