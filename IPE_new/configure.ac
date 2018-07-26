AC_INIT([ipe], [1.0], [ipe_devteam@noaa.gov], [https://github.com/SWPC-IPE/WAM-IPE])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIRS([m4])

# Detect Fortran compilers
AC_PROG_CC
AC_PROG_FC
AC_PROG_F77

AM_PROG_AR
AC_PROG_RANLIB

# NetCDF Tests
AX_LIB_NETCDF4
if test "$with_netcdf4" = "yes"; then
  CC="$NETCDF4_CC"
else
  AC_MSG_ERROR([Unable to find NetCDF4, we need NetCDF4.])
fi

if test "$with_netcdf4_fortran" = "yes"; then
  LDFLAGS+="$NETCDF4_LDFLAGS" 
  LDFLAGS+="$NETCDF4_FLIBS" 
  LDFLAGS+="$NETCDF4_CLIBS" 
else
  AC_MSG_ERROR([Unable to find NetCDF4-Fortran, we need NetCDF4-Fortran.])
fi



AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/flip/Makefile
 src/msise00/Makefile
])
AC_OUTPUT
