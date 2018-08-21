bin_PROGRAMS = ipe legacy2netcdf

ipe_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90 IPE_Driver.F90

legacy2netcdf_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90 Convert_LegacyPlasma_to_NetCDF.F90

AM_CPPFLAGS = ${NETCDF4_FFLAGS} -I flip/ -I msise00/ -I empirical_efield/
ipe_LDADD = flip/libflip.a msise00/libmsise.a empirical_efield/libefield.a
legacy2netcdf_LDADD = flip/libflip.a msise00/libmsise.a empirical_efield/libefield.a
