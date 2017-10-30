#!/bin/bash

if [ -e petsc-build ]; then
    echo "Directory petsc-build exists, please remove it."
    exit 1;
fi

mkdir petsc-build; cd petsc-build;

echo Downloading petsc-3.7.7.tar.gz package
wget -q http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.7.tar.gz

echo Unpacking petsc-3.7.7.tar.gz package
tar -xf petsc-3.7.7.tar.gz
cd petsc-3.7.7

echo Configuring petsc-3.7.7
PETSC_DIR=$(pwd) ./configure --prefix=${HOME}/.local/opt/petsc/petsc-3.7.7/ > /dev/null

echo Building petsc-3.7.7
PETSC_DIR=$(pwd) make > /dev/null

echo Installing petsc-3.7.7
PETSC_DIR=$(pwd) make install > /dev/null

echo Cleaning
cd ../../
rm -rf petsc-build

echo Done. You can set PETSC_DIR to ${HOME}/.local/opt/petsc/petsc-3.7.7/
