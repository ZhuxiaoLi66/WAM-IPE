

ifeq (${FC},ifort)

  CPPFLAGS=-fpp

  ifeq (${DEBUG},yes)
 
    OPT=-g -O0 -ftz -traceback -check all -debug all -warn all -fpe0

  endif

endif

ifeq (${FC},gfortran)
  CPPFLAGS=-cpp
endif

ifeq (${FC},pgfortran)
  CPPFLAGS=-Mpreprocess
endif


# IPE Model CPP-Flags
ifeq (${COUPLED},yes)

  CPPFLAGS+= -DCOUPLED

  ifeq (${COUPLED_TO_WAM},yes)
    CPPFLAGS+= -DCOUPLED_TO_WAM
  endif

endif


LIB=${NCLIB}
INC=${NCINC}

# MSISE00 Library
libnrlmsis.a:
	make --directory=./msise00/ libnrlmsis.a

# FLIP Library


# IPE Objects
IPE_Precision.o : IPE_Precision.F90
	${FC} ${CPPFLAGS} -c IPE_Precision.F90 -o $@

IPE_Common_Routines.o : IPE_Precision.o IPE_Common_Routines.F90
	${FC} ${CPPFLAGS} -c IPE_Common_Routines.F90 ${LIB} ${INC} -o $@

IPE_Constants_Dictionary.o : IPE_Precision.o IPE_Constants_Dictionary.F90 
	${FC} ${CPPFLAGS} -c IPE_Constants_Dictionary.F90 -o $@

IPE_Model_Parameters_Class.o : IPE_Precision.o IPE_Constants_Dictionary.o IPE_Common_Routines.o IPE_Model_Parameters_Class.F90
	${FC} ${CPPFLAGS} -c IPE_Model_Parameters_Class.F90 ${LIB} ${INC} -o $@

IPE_Grid_Class.o : IPE_Precision.o IPE_Constants_Dictionary.o IPE_Common_Routines.o IPE_Grid_Class.F90
	${FC} ${CPPFLAGS} ${OPT} -c IPE_Grid_Class.F90 ${LIB} ${INC} -o $@

IPE_Neutrals_Class.o : IPE_Precision.o IPE_Constants_Dictionary.o IPE_Grid_Class.o IPE_Neutrals_Class.F90
	${FC} ${CPPFLAGS} ${OPT} -c IPE_Neutrals_Class.F90 -o $@



# Test Programs

OBJ = IPE_Precision.o \
      IPE_Common_Routines.o \
      IPE_Constants_Dictionary.o \
      IPE_Model_Parameters_Class.o \
      IPE_Grid_Class.o \
      Convert_LegacyGrid_to_NetCDF.o 

legacy2netcdf : ${OBJ} 
	${FC} ${OPT} ${OBJ} ${LIB} ${INC} -o $@
	mv $@ ../bin/

Convert_LegacyGrid_to_NetCDF.o : IPE_Grid_Class.o IPE_Model_Parameters_Class.o test/Convert_LegacyGrid_to_NetCDF.F90
	${FC} ${CPPFLAGS} ${OPT} -c test/Convert_LegacyGrid_to_NetCDF.F90 ${LIB} ${INC} -o $@ 

.PHONY : clean

clean :
	rm *.o *.mod
