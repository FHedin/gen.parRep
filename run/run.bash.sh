#!/bin/bash

ORIG_DIR=$PWD
TMP_DIR=$(mktemp --tmpdir=$ORIG_DIR -u)
mkdir -p $TMP_DIR && echo "Running in tmp dir $TMP_DIR"

cd $TMP_DIR
ln -s $ORIG_DIR/mol

#ln -s $ORIG_DIR/mol/ala2vac/input_parRep.lua input.lua
#ln -s $ORIG_DIR/mol/ala2vac/input_generalized_parRep.lua input.lua

#ln -s $ORIG_DIR/mol/ala2vac_transient/input_parRep.lua input.lua
#ln -s $ORIG_DIR/mol/ala2vac_transient/input_generalized_parRep.lua input.lua

ln -s $ORIG_DIR/mol/1d7h_i/input_generalized_parRep.lua input.lua

ln -s $ORIG_DIR/parRep

export OPENMM_CPU_THREADS=1
export OPENMM_PLUGIN_DIR=$ORIG_DIR/../build/openmm-7.3.0/lib/plugins

mpirun -np 1  ./parRep -i input.lua -log dbg : \
       -np 7  ./parRep -i input.lua -log dbg -o out.txt -e err.txt

