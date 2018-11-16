#!/bin/bash

if [ ! -d "openmm-7.3.0" ]
then
  if [ ! -f openmm-7.3.0-py27_cuda92_1.tar.bz2 ]
  then
    wget "https://anaconda.org/omnia/openmm/7.3.0/download/linux-64/openmm-7.3.0-py27_cuda92_rc_1.tar.bz2"
  fi
  mkdir openmm-7.3.0
  cd openmm-7.3.0
  tar xf ../openmm-7.3.0-py27_cuda92_rc_1.tar.bz2
  cd ..
fi
export OMM_INSTALL_DIR=$PWD/openmm-7.3.0

if [ ! -d "LuaJIT-2.0.5" ]
then
  if [ ! -f LuaJIT-2.0.5.tar.gz ]
  then
    wget "https://luajit.org/download/LuaJIT-2.0.5.tar.gz"
  fi
  tar xf LuaJIT-2.0.5.tar.gz
  cd LuaJIT-2.0.5
  make PREFIX=$PWD/build -j 8 install
  cd ..
fi
export LUAJIT_INSTALL_DIR=$PWD/LuaJIT-2.0.5/build

CC=gcc CXX=g++ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$OMM_INSTALL_DIR;$LUAJIT_INSTALL_DIR" ..

