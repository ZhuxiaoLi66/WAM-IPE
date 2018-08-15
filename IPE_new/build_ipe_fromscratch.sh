#!/bin/sh

IPE_PATH=$(pwd)

autoreconf --install
./configure --prefix="${IPE_PATH}/install"
make
make install

# Make the run directory and add symbolic link of executable from install directory
mkdir run
ln -s ${IPE_PATH}/install/bin/ipe ${IPE_PATH}/run/ipe

