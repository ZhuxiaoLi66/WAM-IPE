

FC  = ifort
OPT = -O0 -g -traceback -check all -debug all -fpe0 -fpp
LIB = -L/apps/netcdf/4.3.0-intel/lib -lnetcdff -lnetcdf
INC = -I/apps/netcdf/4.3.0-intel/include

OBJS = module_precision.o \
       module_IPE_dimension.o \
       module_input_parameters.o \
       module_FIELD_LINE_GRID_MKS.o \
       IPEToHeightGrid.o

i2hg : ${OBJS}	
	${FC} ${OPT} ${OBJS} ${LIB} ${INC} -o $@
	mv $@ ../../bin/

IPEToHeightGrid.o : IPEToHeightGrid.f90 module_precision.o module_IPE_dimension.o module_input_parameters.o module_FIELD_LINE_GRID_MKS.o
	${FC} ${OPT} -c IPEToHeightGrid.f90 ${LIB} ${INC} -o $@

module_input_parameters.o : ../main/module_input_parameters.f90 module_precision.o module_IPE_dimension.o
	${FC} ${OPT} -c ../main/module_input_parameters.f90 -o $@

module_IPE_dimension.o: ../main/module_IPE_dimension.f90 module_precision.o
	${FC} ${OPT} -c ../main/module_IPE_dimension.f90 -o $@

module_FIELD_LINE_GRID_MKS.o : ../main/module_field_line_grid.f90 module_precision.o module_IPE_dimension.o
	${FC} ${OPT} -c ../main/module_field_line_grid.f90 -o $@

module_precision.o : ../main/module_precision.f90
	${FC} ${OPT} -c ../main/module_precision.f90 -o $@
