#!/bin/bash
#set -ax
# ********* Settings that you should modify ********* #


# Reference output directory
REFDIR='/scratch3/NCEPDEV/swpc/scrub/Joseph.Schoonover/ptmp/Joseph.Schoonover/march2013_old_mediator_1x80_5day/'

# Comparison directory
COMPDIR='/scratch3/NCEPDEV/swpc/scrub/Raffaele.Montuoro/ptmp/Raffaele.Montuoro/march2013_updated_mediator_1x80_5d/'



# ----------------------------------------------------------------------- #
 

# Create a list of neutral files that are not bitwise identical

for file in ${REFDIR}ipe_grid_neutral_params.*
do
  filename=$(echo $file | awk -F / '{print $NF}')
  if [ -e $COMPDIR/$filename ]
  then

     if ! cmp ${REFDIR}${filename} ${COMPDIR}${filename} >/dev/null 2>&1
     then

       echo $filename
       echo $filename >> difference_list.txt

       echo $filename >> comparisons.txt
       echo '-----------------------' >> comparisons.txt
       ./comp ${REFDIR}${filename} ${COMPDIR}${filename}  > raw.txt
       grep "DIFF REP" raw.txt >> comparisons.txt
       echo '-----------------------' >> comparisons.txt
 
     fi
 


  fi

done


#diff -q ${REFDIR}'ipe_grid_neutral_params.*' ${COMPDIR}'ipe_grid_neutral_params.*'


#more neutral_files.txt

