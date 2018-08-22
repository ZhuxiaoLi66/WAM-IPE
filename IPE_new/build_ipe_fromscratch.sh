#!/bin/sh

IPE_PATH=$(pwd)

autoreconf --install
./configure --prefix="${IPE_PATH}/install"
make
make install

# Make the run directory and add symbolic link of executable from install directory
if [ ! -d run ];then
  mkdir run
fi

if [ ! -f ${IPE_PATH}/run/ipe ];then
  ln -s ${IPE_PATH}/install/bin/ipe ${IPE_PATH}/run/ipe
fi

if [ ! -f ${IPE_PATH}/run/legacy2netcdf ];then
  ln -s ${IPE_PATH}/install/bin/legacy2netcdf ${IPE_PATH}/run/legacy2netcdf
fi

if [ ! -f ${IPE_PATH}/run/eregrid ];then
  ln -s ${IPE_PATH}/install/bin/eregrid ${IPE_PATH}/run/eregrid
fi

