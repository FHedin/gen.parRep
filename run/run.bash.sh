#!/bin/bash

# --------------------------------------------------------------
# create a local tmp dir for the run 

ORIG_DIR=$PWD
TMP_DIR=$(mktemp --tmpdir=$ORIG_DIR -u)
mkdir -p $TMP_DIR && echo "Running in tmp dir $TMP_DIR"

cd $TMP_DIR
ln -s $ORIG_DIR/../mol

# --------------------------------------------------------------

# use one of the following input files

#ln -s mol/ala2vac/input_parRep.lua input.lua
#ln -s mol/ala2vac/input_generalized_parRep.lua input.lua

#ln -s mol/ala2vac_transient/input_parRep.lua input.lua
ln -s mol/ala2vac_transient/input_generalized_parRep.lua input.lua

#ln -s mol/1d7h_i/input_generalized_parRep.lua input.lua

ln -s $ORIG_DIR/../build/parRep

# When using the OMM CPU platform this decides how many CPU threads per replica to use
export OPENMM_CPU_THREADS=1
# if the following is not present only the slow Reference OMM platform will be available
export OPENMM_PLUGIN_DIR=$ORIG_DIR/../build/openmm-7.3.0/lib/plugins

# to save random seeds to a file for later use
mpirun -np 1  ./parRep -i input.lua -log dbg --out-seeds seeds.bin : \
       -np 7  ./parRep -i input.lua -log dbg -o out.txt -e err.txt --out-seeds seeds.bin

# to run using random seeds from a previous execution for reproducibility (not guaranteed to work on all platforms or if the software environnment is not exactly the same)
# mpirun -np 1  ./parRep -i input.lua -log dbg --inp-seeds seeds.bin : \
#        -np 7  ./parRep -i input.lua -log dbg -o out.txt -e err.txt --inp-seeds seeds.bin
