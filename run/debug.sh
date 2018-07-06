#!/bin/bash

# set this variable to the number of cpu cores you want to allow for each replica
export OPENMM_CPU_THREADS=1
# required for the executable to locate OpenMM plugins
export OPENMM_PLUGIN_DIR=$HOME/bin/openmm/lib/plugins

mpirun -x OPENMM_CPU_THREADS -x OPENMM_PLUGIN_DIR -np 6 xterm -e gdb --ex run --args parRep -i input_file.lua -log dbg

