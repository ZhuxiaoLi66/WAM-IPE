bin_PROGRAMS = ipe legacy2netcdf grid2netcdf eregrid
lib_LIBRARIES = libipe.a

include_HEADERS = ipe_precision.$(FC_MODEXT) ipe_constants_dictionary.$(FC_MODEXT) ipe_common_routines.$(FC_MODEXT) ipe_forcing_class.$(FC_MODEXT) ipe_time_class.$(FC_MODEXT) ipe_grid_class.$(FC_MODEXT) ipe_model_parameters_class.$(FC_MODEXT) ipe_electrodynamics_class.$(FC_MODEXT)  ipe_neutrals_class.$(FC_MODEXT) ipe_plasma_class.$(FC_MODEXT) ipe_model_class.$(FC_MODEXT)

libipe_a_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90

ipe_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90 IPE_Driver.F90

legacy2netcdf_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90 Convert_LegacyPlasma_to_NetCDF.F90

grid2netcdf_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90 Convert_LegacyGrid_to_NetCDF.F90

eregrid_SOURCES = IPE_Precision.F90 IPE_Constants_Dictionary.F90 IPE_Common_Routines.F90 IPE_Forcing_Class.F90 IPE_Time_Class.F90 IPE_Grid_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Electrodynamics_Class.F90  IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90  IPE_Model_Class.F90 Eldyn_Converter.F90

AM_CPPFLAGS = ${NETCDF4_FFLAGS} -I flip/ -I msise00/ -I empirical_efield/
ipe_LDADD = flip/libflip.a msise00/libmsise.a empirical_efield/libefield.a
legacy2netcdf_LDADD = flip/libflip.a msise00/libmsise.a empirical_efield/libefield.a
grid2netcdf_LDADD = flip/libflip.a msise00/libmsise.a empirical_efield/libefield.a
eregrid_LDADD = flip/libflip.a msise00/libmsise.a empirical_efield/libefield.a
