#!/bin/bash -l
#PBS -l mem=32gb
#PBS -l walltime=24:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=2:ppn=20
#PBS -M jakub.krajniak@cs.kuleuven.be
#PBS -A lp_polymer_goa_project

module purge

cd $PBS_O_WORKDIR

module load espressopp/chem
# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)
mpdboot -f $PBS_NODEFILE -n $n_node -r ssh -v

for n in 1 2 4 8 16 20 24 32 40; do
  mpirun -n $n python -u lennard_jones.py &> run_${n}.log
done
