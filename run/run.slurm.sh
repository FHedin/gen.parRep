#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH -J test_ala2vac

ORIG_DIR=$PWD
TMP_DIR=$(mktemp --tmpdir=$ORIG_DIR -u)
mkdir -p $TMP_DIR && echo "Running in tmp dir $TMP_DIR"

cd $TMP_DIR
ln -s $ORIG_DIR/mol
ln -s $ORIG_DIR/mol/ala2vac/input_generalized_parRep.lua input.lua
ln -s $ORIG_DIR/parRep

# set this variable to the number of cpu cores you want to allow for each replica
export OPENMM_CPU_THREADS=1
# required for the executable to locate OpenMM plugins
export OPENMM_PLUGIN_DIR=$HOME/bin/openmm/lib/plugins

# running with 8 replicas (much much more expected for a production run !!)
#  the first rep prints everything to the terminal
#  the 7 following to a text file : out.txt and err.txt are template names, rep 1 will write to out.1.txt, rep 2 to out.2.txt ...
srun --mpi=pmi2 ./parRep -i input.lua -log dbg -o out.txt -e err.txt

