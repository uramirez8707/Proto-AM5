#!/bin/bash

export BASEDIR=`pwd`
echo "BASEDIR = ${BASEDIR}"
cd ./exec

echo "Current directory:"
pwd
echo " "

echo "Sourcing setup_environment.sh"
source env.sh

#make clean
#make -j 4

