#!/usr/bin/python

import subprocess

bucket = 'https://storage.googleapis.com/ipe-data/'

files = ['dwm07b104i.dat',
         'gd2qd.dat',
         'global_idea_coeff_hflux.dat',
         'global_idea_wei96.cofcnts',
         'hwm123114.bin',
         'ionprof',
         'IPE_Grid.nc',
         'IPE_State.apex.201303160000.nc',
         'IPE.inp',
         'ipe.pbs',
         'tiros_spectra']

files = [bucket + f  for f in files]
for f in files:
  print(f)

subprocess.call(['wget'] + files)
                                  
