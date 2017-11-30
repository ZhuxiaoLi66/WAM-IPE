#!/bin/bash

pwd=$(pwd)

## set restart
cycle=${2:-1}
if [[ $cycle == 1 ]] ; then
	export RESTART=.false.
else
	export RESTART=.true.
fi

## source config
. $pwd/config/workflow.sh $pwd/$1

## create job file
tmp=temp_job.sh
rm -rf $tmp
touch $tmp
chmod +x $tmp

## set up PBS/LSF/whatever
echo "#!/bin/bash" > $tmp
# the below is a little hacky but it sure works
d=0
while
        d=$((d+1))
        eval SUBFLAG=\$SUBFLAG$d
        [[ -n "$SUBFLAG" ]]
do
        echo `eval echo $SUBFLAG` >> $tmp
done
## and now back to the regularly scheduled program
cat >> $tmp << EOF
set -ax

cd $SCRIPTSDIR

## set restart
if [[ $cycle == 1 ]] ; then
	export RESTART=.false.
else
	export RESTART=.true.
fi

##-------------------------------------------------------
## source config file
##-------------------------------------------------------

. $pwd/config/workflow.sh $pwd/$1

##-------------------------------------------------------
## execute forecast
##-------------------------------------------------------

rm -rf $RUNDIR
mkdir -p $RUNDIR
cd $RUNDIR

export VERBOSE=YES

. $EXGLOBALFCSTSH
if [ $? != 0 ]; then echo "forecast failed, exit"; exit; fi
echo "fcst done"

if [[ $((cycle+1)) -le $3 ]] ; then
echo "resubmitting $1 for cycle $((cycle+1)) out of $3"
cd $SCRIPTSDIR
. $pwd/submit.sh $1 $((cycle+1)) $3
else
echo "cycle $((cycle+1)) > $3, done!"
fi

exit $status
EOF

$SCHEDULER_SUB < $tmp
rm -rf $tmp
