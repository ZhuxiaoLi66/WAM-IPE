#!/bin/bash

pwd=$(pwd)

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

##-------------------------------------------------------
## source config file
##-------------------------------------------------------

. $pwd/config/workflow.sh $pwd/$1

##-------------------------------------------------------
## execute forecast
##-------------------------------------------------------

mkdir -p $RUNDIR
cd $RUNDIR

export VERBOSE=YES

. $EXGLOBALFCSTSH
if [ $? != 0 ]; then echo "forecast failed, exit"; exit; fi
echo "fcst done"

exit $status
EOF

$SCHEDULER_SUB < $tmp
rm -rf $temp
