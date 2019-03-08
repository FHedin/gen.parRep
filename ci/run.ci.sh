#!/bin/bash

rm -f *.db *.log *.txt

export OPENMM_CPU_THREADS=1
export OPENMM_PLUGIN_DIR=$PWD/../build/openmm-7.3.0/lib/plugins

mpirun -np 8 $PWD/../build/parRep -i input_generalized_parRep.lua -log warn -o out.txt -e err.txt --inp-seeds ./ref.sim/seeds.bin

MPI_RET_CODE=$?

# compare the db file from the run and the reference, they should be the same
sqldiff lj7.db ref.sim/lj7.db > diff.txt

SQL_RET_CODE=$?

if [ -s "diff.txt" ]
then 
  TEST_RET_CODE=-1
else
  TEST_RET_CODE=0
fi

cat out.*.txt

# success only if there was no error code from any of the previous commands
exit $(($MPI_RET_CODE|$SQL_RET_CODE|$TEST_RET_CODE))
