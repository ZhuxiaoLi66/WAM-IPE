

FC  = mpiifort
OPT = -O3 -fpp 
#OPT = -O0 -g -traceback -fpp
LIB = -L${NETCDF}/lib -lnetcdff -lnetcdf
INC = -I${NETCDF}/include

OBJS = module_precision.o \
       module_IPE_dimension.o \
       FLIP_GRID.o \
       module_physical_constants.o \
       module_input_parameters.o \
       get_pvalue_dipole.o \
       module_FIELD_LINE_GRID_MKS.o \
       module_init_plasma_grid.o \
       module_read_plasma_grid_global.o \
       module_io.o \
       module_open_file.o \
       allocate_arrays.o \
       get_flip_grid.o \
       get_sinim.o \
       module_unit_conversion.o \
       IPEToHeightGrid_mpi.o

i2hg_mpi : ${OBJS}	
	${FC} ${OPT} ${OBJS} ${LIB} ${INC} -o $@
	mv $@ ../../bin/

IPEToHeightGrid_mpi.o : IPEToHeightGrid_mpi.f90 module_precision.o module_IPE_dimension.o module_input_parameters.o \
                    module_FIELD_LINE_GRID_MKS.o module_init_plasma_grid.o
	${FC} ${OPT} -c IPEToHeightGrid_mpi.f90 ${LIB} ${INC} -o $@

allocate_arrays.o : ../main/allocate_arrays.f90 \
                    module_precision.o module_IPE_dimension.o module_FIELD_LINE_GRID_MKS.o \
                    module_input_parameters.o
	${FC} ${OPT} -c ../main/allocate_arrays.f90 -o $@

FLIP_GRID.o : ../flip/FLIP_GRID.f
	${FC} ${OPT} -c ../flip/FLIP_GRID.f -o $@

get_flip_grid.o : ../main/get_flip_grid.f90 FLIP_GRID.o \
                  module_precision.o module_FIELD_LINE_GRID_MKS.o module_input_parameters.o \
                  module_physical_constants.o module_unit_conversion.o
	${FC} ${OPT} -c ../main/get_flip_grid.f90 -o $@

get_sinim.o : ../plasma/get_sinim.f90 \
              module_precision.o
	${FC} ${OPT} -c ../plasma/get_sinim.f90 -o $@

module_unit_conversion.o : ../main/module_unit_conversion.f90 module_precision.o
	${FC} ${OPT} -c ../main/module_unit_conversion.f90 -o $@

get_pvalue_dipole.o : ../plasma/get_pvalue_dipole.f90 module_precision.o
	${FC} ${OPT} -c ../plasma/get_pvalue_dipole.f90 -o $@

module_init_plasma_grid.o : ../main/module_init_plasma_grid.f90 \
                            module_read_plasma_grid_global.o module_IPE_dimension.o \
                            module_precision.o module_physical_constants.o \
                            module_input_parameters.o module_FIELD_LINE_GRID_MKS.o
	${FC} ${OPT} -c ../main/module_init_plasma_grid.f90 -o $@

module_read_plasma_grid_global.o : ../main/module_read_plasma_grid_global.f90 \
                                    module_precision.o module_IPE_dimension.o module_physical_constants.o \
                                   module_input_parameters.o module_io.o module_FIELD_LINE_GRID_MKS.o \
                                   module_open_file.o
	${FC} ${OPT} -c ../main/module_read_plasma_grid_global.f90 -o $@

module_physical_constants.o : ../main/module_physical_constants.f90 \
                              module_precision.o
	${FC} ${OPT} -c ../main/module_physical_constants.f90 -o $@

module_io.o : ../main/module_io.f90 module_precision.o module_IPE_dimension.o
	${FC} ${OPT} -c ../main/module_io.f90 -o $@

module_open_file.o : ../main/module_open_file.f90 module_precision.o module_IPE_dimension.o \
                     module_input_parameters.o
	${FC} ${OPT} -c ../main/module_open_file.f90 -o $@

module_input_parameters.o : ../main/module_input_parameters.f90 module_precision.o module_IPE_dimension.o
	${FC} ${OPT} -c ../main/module_input_parameters.f90 -o $@

module_IPE_dimension.o: ../main/module_IPE_dimension.f90 module_precision.o
	${FC} ${OPT} -c ../main/module_IPE_dimension.f90 -o $@

module_FIELD_LINE_GRID_MKS.o : ../main/module_field_line_grid.f90 module_precision.o module_IPE_dimension.o
	${FC} ${OPT} -c ../main/module_field_line_grid.f90 -o $@

module_precision.o : ../main/module_precision.f90
	${FC} ${OPT} -c ../main/module_precision.f90 -o $@
