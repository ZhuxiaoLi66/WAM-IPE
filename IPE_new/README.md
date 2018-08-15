# Ionosphere Plasmasphere Electrodynamics (IPE) Model

# Getting started

## Dependencies
The IPE code depends on 
NetCDF and NetCDF-Fortran.

## Building the code
The IPE build system uses autotools to generate makefiles to compile 
source code to a library( libipe.a ) and executables ( ipe ).

To build the code, from the head directory do

```
autoreconf --install
./configure --prefix=</path/to/install>
make
make install
```

These instructions create the following directories : 
```
</path/to/install>/bin
</path/to/install>/share
```

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
