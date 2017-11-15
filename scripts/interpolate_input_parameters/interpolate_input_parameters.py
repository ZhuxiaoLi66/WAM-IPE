#!/usr/bin/env python
import numpy as np
from os import path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import datetime
from itertools import chain

def failure(fail_point):
  # fail: string
  print('error during '+fail_point)
  exit(0)

def flatten(l):
  # l: 2D list to convert to 1D
  return list(chain.from_iterable(l))

def interpolate(arr,mins_per_segment):
  output = []
  for pair in zip(arr,arr[1:]):
    output.append(numpy.linspace(pair[0],pair[1],mins_per_segment+1)[:-1])
  output.append([arr[-1]])
  return output

def get_kp_f107(start,end):
  dates = [start]
  next = new_timestamp(dates[-1],24)
  while next != end:
    dates.append(next)
    next = new_timestamp(dates[-1],24)
  dates.append(end)
  print dates

  kp = []
  f107 = []

  try:
    for cdate in dates:              # for each date to pull info for
      with open(cdate[:4]) as file:  # open yearly database
        for line in file:            # for each file
          if line[:6] == cdate[2:8]: # match date
            # append to lists
            kp.append([real(line[i:i+2])/10 for i in range(12,28,2)])
            f107 = float(line[65:71])

  except:
    failure('yearly kp_ap database read')

  return flatten(interpolate(flatten(kp))),flatten(interpolate(flatten(f107)))

def new_timestamp(init,diff):
  # init: YYYYMMDDHH string
  # diff: integer hours
  dt = datetime.datetime.strptime(init,'%Y%m%d%H')
  return dt+datetime.timedelta(0,60*60*diff)

def parse(start_date,end_date):
  # start_date: YYYYMMDDHH string
  # end_date:   YYYYMMDDHH string
  starting_min = float(start_date[-2:])*60
  ending_min   = float(end_date[-2:]  )*60

  # first determine which dates we need to pull data for
  if float(starting_min) / min_per_kp_segment > 0.5: # start interpolation from current day
    min_f107 = start_date
  else:                                              # start interpolation from prior day
    min_f107 = new_timestamp(start_date, 24)
  
  if float(ending_min) / min_per_kp_segment > 0.5:   # end interpolation at next day
    max_f107 = new_timestamp(end_date, 24)
  else:                                              # end interpolation at current day
    max_f107 = end_date

  # now get it for the days described
  kp,f107 = get_kp_f107(min_f107,max_f107)

  # and select the right ones
  
  return kp[

def output(kp,f107,power):

## parsing options
parser = ArgumentParser(description='Parse KP_AP files into binned data', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--interval',  help='interval length (minutes) (default=1)', type=int, default=1)
parser.add_argument('-d', '--duration',  help='duration of run (hours) (default=24)', type=int, default=24)
parser.add_argument('-s', '--startdate', help='starting date of run (YYYYMMDDhh)', required=True, type=str)
args = parser.parse_args()

## globals
min_per_kp_segment   = 3*60
min_per_f107_segment = 24*60
end_date = new_timestamp(args.startdate,args.duration)
