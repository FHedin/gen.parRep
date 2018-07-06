#!/bin/bash

# set this variable to the number of cpu cores you want to allow for each replica
export OPENMM_CPU_THREADS=1
# required for the executable to locate OpenMM plugins
export OPENMM_PLUGIN_DIR=$HOME/bin/openmm/lib/plugins

export VGRIND=/usr/local/bin/valgrind

export LD_PRELOAD=/usr/local/lib/valgrind/libmpiwrap-amd64-linux.so

mpirun -np 32 \
$VGRIND \
parRep -i input_file.lua -log dbg -o out.txt -e err.txt

