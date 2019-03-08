#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH -J test_ala2

ORIG_DIR=$PWD
TMP_DIR=$(mktemp --tmpdir=$ORIG_DIR -u)
mkdir -p $TMP_DIR && echo "Running in tmp dir $TMP_DIR"

cd $TMP_DIR
ln -s $ORIG_DIR/../mol
ln -s mol/ala2vac_transient/input_generalized_parRep.lua input.lua
ln -s $ORIG_DIR/../build/parRep

export OPENMM_CPU_THREADS=1
export OPENMM_PLUGIN_DIR=$ORIG_DIR/../build/openmm-7.3.0/lib/plugins

srun --mpi=pmi2 ./parRep -i input.lua -log warn -o out.txt -e err.txt

