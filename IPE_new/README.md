# Ionosphere Plasmasphere Electrodynamics (IPE) Model

# Reporting Issues
Issues should be reported to

https://github.com/SWPC-IPE/WAM-IPE/issues


# Getting started

## Dependencies
The IPE code depends on 
NetCDF and NetCDF-Fortran.

## Building the code
The IPE build system uses autotools to generate makefiles to compile 
source code to a library( libipe.a ) and executables ( ipe ).

To build the code, from the head directory do

```
./build_ipe_fromscratch.sh
```

These instructions create the following directories : 
```
$ IPE_new/install/bin
$ IPE_new/install/share
$ IPE_new/run
```
and sets up the following symbolic links
```
IPE_new/run/ipe => IPE_new/install/bin/ipe
IPE_new/run/eregrid => IPE_new/install/bin/eregrid
IPE_new/run/legacy2netcdf => IPE_new/install/bin/legacy2netcdf
```
If you modify code and want to update the executables, simply run
`make install` from the `IPE_new` directory.

### Executables
`ipe`
The main executable that forward steps the IPE model based on parameters within
the IPE.inp file


`eregrid`
This executable is used to convert electric potential output from other
models to the IPE grid. Currently, this executable simply calls the
Update routine for IPE_Electrodynamics. Tools for working with OpenGGCM
and Geospace are currently in production.


`legacy2netcdf`
This executable is used to convert legacy IPE pickup files to the currently used
NetCDF format.


## Running the model
If you're on a system that can access https://storage.googleapis.com, downloading the
input decks for running the demo version of IPE, simply run
```
python scripts/download_ipe-data.py
```
If you need input-decks, e-mail Joseph.Schoonover@noaa.gov.

Once you download the data, place the data in a directory with the ipe executable. The
executable can be found under
```
</path/to/install>/bin/ipe
```

On systems with a PBS batch scheduler, you can use the sripts/ipe.pbs job submission
script as a template.

Otherwise, you can run simply by doing
```
./ipe
```
in the directory where your input data is located.

## Linking IPE to other code
If you want to use IPE's API within your own program, you'll need to link to libipe.a,
libmsise.a, libflip.a, and libefield.a and include the path to the .mod files
