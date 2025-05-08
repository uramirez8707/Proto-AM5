#!/bin/bash

source $MODULESHOME/init/bash

module unload cray-netcdf cray-hdf5 fre
module unload PrgEnv-pgi PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module load PrgEnv-intel/8.5.0
module unload intel intel-classic intel-oneapi
module load intel-classic/2023.2.0
module load cray-hdf5/1.12.2.11
module load cray-netcdf/4.9.0.11
module load libyaml/0.2.5

ulimit -s unlimited

module list
