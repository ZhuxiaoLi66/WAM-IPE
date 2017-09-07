#!/bin/bash
set -ax

CONFIG=$1

## source config
. $CONFIG

## create job file
tmp=temp_job.sh
rm -rf $tmp
touch $tmp
chmod +x $tmp

cat >> $tmp << EOF
#!/bin/bash
`eval echo "$NAMEFLAG"`
`eval echo "$PROJFLAG"`
`eval echo "$WALLFLAG"`
`eval echo "$STDOUTFLAG"`
`eval echo "$STDERRFLAG"`
`eval echo "$SELECTFLAG"`
`eval echo "$QUEUEFLAG"`
`eval echo "$SPANFLAG"`
`eval echo "$EXTRAFLAG1"`
`eval echo "$EXTRAFLAG2"`
`eval echo "$EXTRAFLAG3"`

set -ax

cd $SCRIPTSDIR

##-------------------------------------------------------
## source config file
##-------------------------------------------------------

. $CONFIG

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
