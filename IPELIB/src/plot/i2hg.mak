

FC  = ifort
OPT = -O3     

i2hg : interface_to_fixed_height_single_file.f90	
	${FC} ${OPT} interface_to_fixed_height_single_file.f90 -o $@
	mv $@ ../../bin/
