#!/bin/bash
re='^[0-9]+$' # used to check for integers

## should probably compartmentalize this file
echo "running checks for $FMID"
echo $TASKS
# basic directory/executable checks
for var in SIGHDR SFCHDR NEMSIOGET APRUN NDATE MDATE FIXGLOBAL BASE_NEMS FCSTEXEC EXGLOBALFCSTSH MTNVAR O3FORC O3CLIM \
	   OROGRAPHY OROGRAPHY_UF LONSPERLAT LONSPERLAR ; do
	eval val=\$$var	# get the value ($val) of the variable that $var references
	echo "checking for $var: $val"
	if [ -z $val ] || [ ! -e $val ] ; then # makes sure the variable has been set and the path/file exists
		echo "   $var not found! exiting."
		if [ $var = 'APRUN' ] ; then
			echo "   common issue: mpi library isn't loaded, or the full path hasn't been specified"
		fi
		exit 1
	fi
done

# frequency checks
echo "checking properties of the frequency options (FHOUT, FHRES, FHDFI, FHCYC, FHZER)"
if [ -z $FHOUT ] || [ $FHOUT -gt $FHMAX ] ; then
	echo "   FHOUT is invalid, setting to $FHMAX"
	export FHOUT=$FHMAX
fi
# rule 1: positive DELTIM
if [ -z $DELTIM ] || [ $DELTIM -le 0 ] ; then
	echo "   invalid DELTIM (must be positive)! exiting." ; exit 1
fi
# rule 2: FHOUT positive, multiple of DELTIM to within $tol
if [ -z $FHOUT ] || [ $FHOUT -le 0 ] || [ $(((FHOUT*3600%DELTIM)/DELTIM)) != 0 ] ; then
	echo "   FHOUT invalid: must be positive and a multiple of DELTIM! exiting." ; exit 1
fi
# rule 3: FHSWR positive, multiple of DELTIM to within $tol
if [ -z $FHSWR ] || [ $FHSWR -le 0 ] || [ $(((FHSWR%DELTIM)/DELTIM)) != 0 ] ; then
	echo "   FHSWR invalid: must be positive and a multiple of DELTIM! exiting." ; exit 1
fi
# rule 4: FHLWR, same as rule 3, but also a multiple of FHSWR
if [ -z $FHLWR ] || [ $FHLWR -le 0 ] || [ $(((FHLWR%DELTIM)/DELTIM)) != 0 ] || [ $((FHLWR%FHSWR)) != 0 ] ; then
	echo "   FHLWR invalid: must be positive, a multiple of DELTIM, and a multiple of FHSWR! exiting." ; exit 1
fi
# rule 5: FHZER, same as rule 3, but also a multiple of FHOUT
if [ -z $FHZER ] || [ $FHZER -le 0 ] || [ $(((FHZER*3600%DELTIM)/DELTIM)) != 0 ] || [ $((FHZER%FHOUT)) != 0 ] ; then
	echo "   FHZER invalid: must be positive, a multiple of DELTIM, and a multiple of FHOUT! exiting." ; exit 1
fi
# rule 6: FHRES, same as rule 3, but also a multiple of FHLWR and FHZER
if [ -z $FHRES ] || [ $FHRES -le 0 ] || [ $(((FHRES*3600%DELTIM)/DELTIM)) != 0 ] || [ $((FHRES*3600%FHLWR)) != 0 ] || [ $((FHRES%FHZER)) != 0 ] ; then
	echo "   FHRES invalid: must be positive, a multiple of DELTIM, FHLWR, and FHZER! exiting." ; exit 1
fi
# rule 7: FHDFI, same as rule 3, but also a multiple of FHLWR and <= FHRES
if [ -z $FHDFI ] || [ $FHDFI -lt 0 ] || [ $(((FHDFI*3600%DELTIM)/DELTIM)) != 0 ] || [ $((FHDFI*3600%FHLWR)) != 0 ] || [ $FHDFI -gt $FHRES ] ; then
	echo "   FHDFI invalid: must be non-negative, a multiple of DELTIM, FHLWR, and <= FHRES! exiting." ; exit 1
fi
# rule 8: FHCYC, same as rule 3, but also a multiple of FHLWR
if [ -z $FHCYC ] || [ $FHCYC -lt 0 ] || [ $(((FHCYC%DELTIM)/DELTIM)) != 0 ] || [ $((FHCYC*3600%FHLWR)) != 0 ] ; then
	echo "   FHCYC invalid: must be non-negative, a multiple of DELTIM, and a multiple of FHLWR! exiting." ; exit 1
fi

# wam.config checks
if [ $IDEA = .true. ] ; then
	# WAM-specific directory checks
	for var in FIX_IDEA GRIDSDIR RT_WAM DATADIR WAMINDIR ; do
		eval val=\$$var
		echo "checking for $var: $val"
		if [ -z $val ] || [ ! -e $val ] ; then
			echo "   $var not found! exiting." ; exit 1
		fi
	done

	echo "checking to make sure CDATE is a valid length: $CDATE"
	if [ ${#CDATE} != 10 ] ; then # y10k problem!
		echo "   $CDATE is ${#CDATE} characters long, needs to be 10 (YYYYMMDDHH)! exiting." ; exit 1
	fi

	# ipe.config checks
	if [ $WAM_IPE_COUPLING = .true. ] ; then
		echo "checking to make sure NPROCIPE is valid: $NPROCIPE"
		valid=0
		for i in "${VALID_IPE_PET[@]}" ; do
			if [ $NPROCIPE = $i ] ; then valid=1 ; fi
		done
		[ $valid = 0 ] && echo "   could not match NPROCIPE against the list of allowed configurations, check ipe.config! exiting." && exit 1

		echo "checking IPEDECOMP against NPROCIPE: (${IPEDECOMPARR[0]}x${IPEDECOMPARR[1]})"
		if [[ ${#IPEDECOMPARR[*]} -ne 2 ]] || ! [[ ${IPEDECOMPARR[0]} =~ $re ]] || ! [[ ${IPEDECOMPARR[1]} =~ $re ]] ; then
			echo "   bad IPEDECOMP, must be two space-separated integers! exiting." ; exit 1
		else
			if [[ $((${IPEDECOMPARR[0]}*${IPEDECOMPARR[1]})) != $NPROCIPE ]] ; then
				echo "IPEDECOMP (${IPEDECOMPARR[0]}x${IPEDECOMPARR[1]}) inconsistent with NPROCIPE: ${NPROCIPE}! exiting." ; exit 1
			fi
		fi
		if [ ! -z $IPEDATETIME ] ; then
			echo "checking to make sure IPEDATETIME is a valid length: $IPEDATETIME"
			if [ ${#IPEDATETIME} != 12 ] ; then # also y10k problem!
				echo "   $IPEDATETIME is ${#IPEDATETIME} characters long, needs to be 13 (YYYYMMDDHHmm)! exiting." ; exit 1
			fi
		fi
		echo "checking for IPEGRID: $IPEGRID"
		if [ ! -e $IPEGRID ] ; then
			echo "   IPEGRID not found! exiting." ; exit 1
		fi
		echo "checking to make sure IPE output frequency is valid (IPEFREQ): $IPEFREQ"
		if [ $IPEFREQ -gt $((FHMAX*3600)) ] ; then # IPEFREQ > model run time
			echo "   IPEFREQ too high... setting IPEFREQ to $((FHMAX*3600))"
			export IPEFREQ=$((FHMAX*3600))
		elif [ $IPEFREQ -lt $DELTIM ] ; then # IPEFREQ < model integration step time
			echo "   setting IPEFREQ to DELTIM"
			export IPEFREQ=$DELTIM
		elif [ ! $(($IPEFREQ % $DELTIM)) = '0' ] ; then # IPEFREQ not a multiple of DELTIM
			echo "   $IPEFREQ is not a multiple of DELTIM: $DELTIM"
			echo "   setting IPEFREQ to $((FHMAX*3600))"
			export IPEFREQ=$((FHMAX*3600))
		fi
	fi
fi

# check for ROTDIR
echo "checking for ROTDIR: $ROTDIR"
if [ -z $ROTDIR  ] || [ ! -d $ROTDIR ] ; then
	"   ROTDIR not found. creating $ROTDIR"
	mkdir -p $ROTDIR
fi
# check if ICs are in place
if [ $RESTART = .false. ] ; then # cold start
	echo "checking for atmospheric/surface initial conditions in ROTDIR"
	if [ $NEMSIO_IN = .true. ] ; then # nemsio initial conditions
		IC_SEARCH=$NEMSIOSEARCH
		echo "   assuming nemsio format for NEMSIO_IN=.true., searching for $IC_SEARCH$CDATE"
	else
		IC_SEARCH=$SIGSFCIOSEARCH
		echo "   assuming sigio/sfcio format for NEMSIO_IN=.false., searching for $IC_SEARCH$CDATE"
	fi
	# this [ -n "$(ls $VAR)" ] format is the best way I found to reliably search for wildcard matches
	# could probably be improved upon?
	if ! [[ -n "$(ls $ROTDIR/$IC_SEARCH$CDATE)" ]] ; then
		echo "   ICs not found in ROTDIR. checking for IC_DIR"
		if [ ! -z $IC_DIR ] ; then
			echo "   IC_DIR has been set to $IC_DIR"
			if ! [[ -n "$(ls $IC_DIR/$IC_SEARCH$CDATE)" ]] ; then
				echo "   but ICs for $CDATE are not found! exiting." ; exit 1
			else
				echo "   found ICs! copying over."
				$NCP $IC_DIR/$IC_SEARCH$CDATE $ROTDIR/.
			fi
		else # IC_DIR is unset, we don't know where to look
			echo "   IC_DIR has not been set: cannot find initial conditions! exiting." ; exit 1
		fi
	fi
	# now we check to see that the surface idate&fhour match the atmospheric idate&fhour
	echo "making sure our ICs match idate and fhour"
	if [ $NEMSIO_IN = .true. ] ; then
		export ATMIN=`ls $ROTDIR/${NEMSIOATM}*$CDATE | head -1`
		export SFCIN=`ls $ROTDIR/${NEMSIOSFC}*$CDATE | head -1`
		if [ $($NEMSIOGET $ATMIN fhour) != $($NEMSIOGET $SFCIN fhour) ] || \
		   [ $($NEMSIOGET $ATMIN idate) != $($NEMSIOGET $SFCIN idate) ] ; then
			echo "   $ATMIN and $SFCIN do not have matching fhour and idate! exiting." ; exit 1
		fi
	else
		export ATMIN=`ls $ROTDIR/${SIGIOATM}*$CDATE | head -1`
		export SFCIN=`ls $ROTDIR/${SFCIOSFC}*$CDATE | head -1`
		if [ $($SIGHDR $ATMIN fhour) != $($SFCHDR $SFCIN fhour) ] || \
		   [ $($SIGHDR $ATMIN idate) != $($SFCHDR $SFCIN idate) ] ; then
			echo "   $ATMIN and $SFCIN do not have matching fhour and idate! exiting." ; exit 1
		fi
	fi
else # restart conditions
	echo "checking for atmospheric/surface restart files in RESTARTDIR $RESTARTDIR"
	NFHOUR_ARR=()
	IDATE_ARR=()
	for file in $SIGR1 $SIGR2 $SFCR $GRDR1 $GRDR2 $FORT1051 ; do
		if ! [[ -n "$(ls $file)" ]] ; then # file not found
			echo "   restart file $file not found! exiting." ; exit 1
		else # file found
			# it's good that we found the file, but we also want to be sure these files are compatible.
			# in other words, we want to have idate and fhour that all match.
			if [ $file != $FORT1051 ] ; then # $FORT1051 is binary format
			# we pull the nemsio nfhour and idate into arrays:
			nfhour=`$NEMSIOGET ${file} nfhour | tr -s ' ' | cut -d' ' -f 3   | sed -e 's/\s//g'`
			idate=` $NEMSIOGET ${file} idate  | tr -s ' ' | cut -d' ' -f 3-7 | sed -e 's/\s//g'`
			# likely unnecessary to separate out the array addition, but we do it anyway
	                NFHOUR_ARR+=("$nfhour")
	                IDATE_ARR+=("$idate")
			fi
	        fi
	done

	# count the number of unique entries in our arrays
	uniq_fhour=($(echo "${NFHOUR_ARR[@]}" | tr ' ' '\n' | sort -u | wc -l))
	uniq_idate=($(echo "${IDATE_ARR[@]}" | tr ' ' '\n' | sort -u | wc -l))
	# if they are not one, we have disagreeing restart files
	if [ ! $uniq_fhour = 1 ] || [ ! $uniq_idate = 1 ] ; then
		echo "   there's an issue with your restart files; check the idate and nfhour to make sure they match! exiting." ; exit 1
	fi
fi # end restart block

# now do the IPE initial conditions
if [ $WAM_IPE_COUPLING = .true. ] ; then
	if [ -z $IPEDATETIME ] ; then # user has not defined IPEDATETIME
		# if RESTART=.false. then CDATE should be valid and we can skip this nonsense
		if [ $RESTART = .false. ] ; then
			CIPEDATE=$CDATE${IPE_MINUTES:-00}
		else # if it is true, then we need to figure out where we actually are in time
			# can use NFHOUR_ARR from above to our advantage, but we need to pad the IDATE out with zeros, as nemsio_get does not
			idate=` $NEMSIOGET $SIGR1 idate  | tr -s ' ' | cut -d' ' -f 3-7`
			iyear=` echo $idate | cut -d' ' -f 1`
			imonth=`printf "%02d" $(echo $idate | cut -d' ' -f 2)`
			iday=`  printf "%02d" $(echo $idate | cut -d' ' -f 3)`
			ihour=` printf "%02d" $(echo $idate | cut -d' ' -f 4)`
			idate=${iyear}${imonth}${iday}${ihour}
			CIPEDATE=`$NDATE $NFHOUR_ARR $idate`
			# NDATE adds NFHOUR_ARR[0] to idate properly
		fi
	fi
	IPESEARCH=$IPEBASESEARCH$CIPEDATE
	echo "searching ROTDIR, then RESTARTDIR, then IPE_IC_DIR for $IPESEARCH"
	# then we search ROTDIR, then RESTARTDIR, then IPE_IC_DIR
	if [[ -n "$(ls $ROTDIR/$IPESEARCH)" ]] ; then
		echo "   found in ROTDIR"
	elif [[ -n "$(ls $RESTARTDIR/$IPESEARCH)" ]] ; then
		echo "   found in RESTARTDIR, copying to ROTDIR"
		$NCP $RESTARTDIR/$IPESEARCH $ROTDIR
	# IPE_IC_DIR has been defaulted to IC_DIR if IC_DIR is defined
	elif [ ! -z ${IPE_IC_DIR} ] && [[ -n "$(ls $IPE_IC_DIR/$IPESEARCH)" ]] ; then
		echo "   found in IPE_IC_DIR, copying to ROTDIR"
		$NCP $IPE_IC_DIR/$IPESEARCH $ROTDIR
	else # can't find any matching IPE initial conditions
		echo "   can't find your IPE files... check the filename convention! exiting." ; exit 1
	fi
	export PLASI=${PLASI:-$ROTDIR/ipe_grid_plasma_params.$CIPEDATE}
	export NEUTI=${NEUTI:-$ROTDIR/ipe_grid_neutral_params.$CIPEDATE}
fi

echo "our enviroment seems to be good, moving to submit the job"

return 0
