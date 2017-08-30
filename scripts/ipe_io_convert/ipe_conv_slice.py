#!/usr/bin/env python
## Python 3 now
import numpy as np
import math
from subprocess import check_output
from os import path, stat
from sys import exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from scipy.io import FortranFile
from multiprocessing import Pool,Process

def failure(fail_point):
	print('error reading '+fail_point)
	exit(0)

def get_grid_size(var_file):
	# read header from unformatted fortran output file to get record size
	try:
		f = open(var_file,'rb')
	except:
		failure(var_file)
	return read_header(f)

def get_deltim():
	# using the frequency of the plasma update as 
	try:
		file = path.join(args.input_directory,'IPE.inp')
		deltim = int(check_output(['/bin/grep','ip_freq_plasma',file]).strip().decode('utf-8').split('=')[1])//60
	except:
		print('Cannot find IPE.inp, defaulting to 3 minute timesteps')
		deltim = 3
	return deltim

def get_timestep(i):
	# pull time at step i from ut_rec (or stup_ut_rec for stup files)
	try:
		ut_rec = 'stup_ut_rec' if var_basename == 'stup' else 'ut_rec'
		f = open(path.join(args.input_directory,ut_rec))
		for tmp in range(i):
			f.readline()
		minutes = int(f.readline().split(' ')[-1])//60
	except:
		failure('ut_rec')
	return minutes

def get_filesize(var_file):
	return stat(var_file).st_size

def get_total_timesteps(var_file):
	# if the plasma file is 46+TB in size, this method won't work. good until then, though!
	try:
		grid_size = get_grid_size(var_file)
		timesteps = int(get_filesize(var_file)/grid_size/byte_size)
	except:
		failure(var_file)
	return timesteps

def new_timestamp(mdate,init,diff):
	# uses mdate to add the time to the time valid for the IPE output
	try:
		newstamp = check_output([mdate,str(diff),init]).strip().decode("utf-8")
	except:
		failure('mdate?')
	return newstamp[:8] + 'T' + newstamp[8:]

def read_header(var_file):
	return int(np.fromfile(var_file,dtype=np.dtype('Int32'),count=1))//byte_size

def read_record(var_file,record_length):
	return np.fromfile(var_file,dtype=np.dtype('Float32'),count=record_length)

def var_parse(var_file):
	# get the actual filename from global variables
	var_file = path.join(args.input_directory,var_basename+var_ext[var_file])
	print('reading '+var_file)

	# get dimensions for output array
	timesteps=slice[1]-slice[0]+1
	grid_size=get_grid_size(var_file)

	# open file handle
	f = open(var_file,'rb')

	# initialize output array
	var_arr = np.ndarray((timesteps,grid_size),dtype=np.dtype('Float32'))

	# fill it
	for i in range(slice[0]): # seek past old data if slice[0] != 0
		record_length = read_header(f)
		f.seek(record_length*byte_size)
		record_length = read_header(f)
	for i in range(slice[1]-slice[0]+1): # now fill our output array
		record_length = read_header(f)
		var_arr[i,:] = read_record(f,record_length)
		record_length = read_header(f)

	return var_arr

def get_slices(var_file,max_slice_size = 500000000):
	# need to do 2GB+ files in slices because multiprocessing uses an older version of Pickle, and
	# it turns out that slicing is considerably faster anyway. we use 500MB slices by default, runs 
	# about the same as 750MB, and faster than 1.5GB and 2GB... haven't tested beyond that.
	num_slices = int(math.ceil(float(get_filesize(var_file))/max_slice_size))
	total_timesteps = get_total_timesteps(var_file)
	slices = []
	for i in range(num_slices):
		if i+1 < num_slices:
			slices.append([int(i*total_timesteps/num_slices),int((i+1)*total_timesteps/num_slices)-1])
		else: # handle edge case if total_timesteps/num_slices not round
			slices.append([int(i*total_timesteps/num_slices),total_timesteps-1])
	return slices

def write_output(output,i,type):
	# find difference between the start of the run and the current timestep,
	# and add the model timestep because the output routine used to be buggy
	minute_diff = get_timestep(i)-get_timestep(0)+get_deltim()
	file = path.join(args.output_directory,
			 'ipe_grid_'+type+'_params.'+new_timestamp(mdate,args.timestamp,minute_diff))
	print('writing out '+file)
	# scipy FortranFile is nicer code than manually writing headers with numpy
	fout = FortranFile(file,'w')
	for var in range(len(var_ext)):
		fout.write_record(output[var,:]) 

def main(type):
	global slice
	test_file = path.join(args.input_directory,var_basename+var_ext[0])
	slices = get_slices(test_file)
	print (slices)
	for slice in slices:
		print('loading '+type+'...')

		## parallelized input
		p = Pool(len(var_ext))
		output = np.transpose(p.map(var_parse,range(len(var_ext))),(1,0,2))
		
		print('loaded '+type+' files, beginning write sequence')
		
		## parallelized output
		processes = [Process(target=write_output, args=(output[i-slice[0]],i,type,)) for i in range(slice[0],slice[1]+1)]
		for r in processes:
			r.start()
		for r in processes:
			r.join()

### parsing options
parser = ArgumentParser(description='convert old IPE IO to new IPE IO', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input_directory',  help='directory where IPE files are stored', type=str, required=True)
parser.add_argument('-o', '--output_directory', help='directory where output files are stored', type=str, required=True)
parser.add_argument('-t', '--timestamp',        help='starting timestamp in format YYYYMMDDHHmm', type=str, required=True)
parser.add_argument('-s', '--stup',             help='process stup* files instead of plasma*', dest='do_stup', default=False, action='store_true')
args = parser.parse_args()

### global path/file definitions
mdate='/scratch3/NCEPDEV/swpc/save/Adam.Kubaryk/util/mdate.fd/mdate'
plasma_max = 16 # number of plasmaXX files, assumed to be 00 to 15, not sure what the deal is with plasma16?
byte_size = 4

### PLASMA
if args.do_stup:
	var_basename = 'plasma' # 'stup' stup13-15 don't exist?
else:
	var_basename = 'plasma'

var_ext = []
for i in range(plasma_max):
	var_ext.append("{:02d}".format(i))

main('plasma')

### NEUTRAL
var_basename = 'ipe_grid_neut_'
var_ext = ['temp','wind','O_den','N2_den','O2_den']

main('neutral')
