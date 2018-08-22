AC_PREREQ([2.69])
AC_INIT([ipe], [1.0], [ipe_devteam@noaa.gov], [https://github.com/SWPC-IPE/WAM-IPE])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIRS([m4])

# Detect Fortran compilers
AC_PROG_CC
AC_PROG_FC([ifort gfortran pgfortran])
AC_PROG_F77([ifort gfortran pgfortran])

FCFLAGS="-O2"

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

AC_LANG_PUSH([Fortran])

# Add the free-line-length-none flag for GCC compilers
AX_CHECK_COMPILE_FLAG([-ffree-line-length-none], [FCFLAGS+=" -ffree-line-length-none"])

# Check for debug request
# Set --enable-debug yes turns on the -DHAVE_CUDA CPP-flag if the PGI compilers are used.
AC_ARG_ENABLE( debug, AS_HELP_STRING([--enable-debug, default:no]))
AS_IF([test "x$enable_debug" = "xyes"],[

  FCFLAGS="-O0 -g"
  # GCC
  AX_CHECK_COMPILE_FLAG([-fbacktrace], [FCFLAGS+=" -fbacktrace"])
  # Intel
  AX_CHECK_COMPILE_FLAG([-traceback], [FCFLAGS+=" -traceback"])
  AX_CHECK_COMPILE_FLAG([-ftrapuv], [FCFLAGS+=" -ftrapuv"])

])         

FFLAGS="$FCFLAGS"

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/flip/Makefile
 src/msise00/Makefile
 src/empirical_efield/Makefile
])
AC_OUTPUT
