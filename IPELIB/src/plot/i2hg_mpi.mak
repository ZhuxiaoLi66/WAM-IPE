FC  = mpiifort
OPT = -O3 -fpp 
#OPT = -O0 -g -traceback -fpe0
LIB = -L${NETCDF}/lib -lnetcdff -lnetcdf
INC = -I${NETCDF}/include

OBJS = module_precision.o \
       IPEToHeightGrid_mpi.o

i2hg_mpi : ${OBJS}	
	${FC} ${OPT} ${OBJS} ${LIB} ${INC} -o $@
	mv $@ ../../bin/

IPEToHeightGrid_mpi.o : IPEToHeightGrid_mpi.f90
	${FC} ${OPT} -c IPEToHeightGrid_mpi.f90 ${LIB} ${INC} -o $@

module_precision.o : ../main/module_precision.f90
	${FC} ${OPT} -c ../main/module_precision.f90 -o $@

clean :
	rm -rf *.mod *.o
