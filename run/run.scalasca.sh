#!/bin/bash

# set this variable to the number of cpu cores you want to allow for each replica
export OPENMM_CPU_THREADS=1
# required for the executable to locate OpenMM plugins
export OPENMM_PLUGIN_DIR=$HOME/bin/openmm/lib/plugins

scalasca -analyze      -f parrep.filt mpirun -np 8 parRep -i input_file.lua -log info -o out.txt -e err.txt 
#scalasca -analyze -q -t -f parrep.filt mpirun -np 8 parRep -i input_file.lua -log info -o out.txt -e err.txt
