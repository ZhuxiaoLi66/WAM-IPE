cmake_minimum_required (VERSION 2.8)
project (IPE)
set (IPE_VERSION_MAJOR 1)
set (IPE_VERSION_MINOR 0)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# External
set (NETCDF_F90 "YES")
find_package (NETCDF REQUIRED)
target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})


# Libraries with source code shipped with IPE

# FLIP
add_subdirectory(flip)
add_library(flip CMINOR.f CTIPE-int.f ELECXS.f FLIP_GRID.f INIT-PROFILES.f KEMPRN.f MINORA.f Neut_Heating.f Photoel-Freqs.f RSDENA_EXB.f RSJACA.f RSLPSD.f RSLPST.f RSPE2B.f RSPRIM.f RSTEMC_EXB.f Rates.f )

# MSIS
add_subdirectory(msise00)
add_library(nrlmsis physics_constants.f utils_spline.f physics_msis.f utils_constants.f hwm14.f)


# IPE Modules and Classes
add_library(ipe IPE_Common_Routines.F90 IPE_Constants_Dictionary.F90 IPE_FileIO.F90 IPE_Forcing_Class.F90 IPE_Grid_Class.F90 IPE_Model_Class.F90 IPE_Model_Parameters_Class.F90 IPE_Neutrals_Class.F90 IPE_Plasma_Class.F90 IPE_Precision.F90 IPE_Time_Class.F90 )

