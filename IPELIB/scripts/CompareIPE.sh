#!/bin/sh/

#
#


REFFILE=$1
COMFILE=$2
OUTDIR=$3

# Extract the file name without the path
filename=$(echo $1 | awk -F / '{print $NF}')


# Make the output directory if needed
[ ! -d $OUTDIR ] && mkdir -p $OUTDIR


# Take the difference of the two netcdf files passed in
# and store the output in a temporary file (diff.nc)
ncdiff $REFFILE $COMFILE diff.nc

# Use the ncap script to compute absolute max and rms difference
ncap2 -v -O -S ncapOperations.nco diff.nc $OUTDIR/$filename


