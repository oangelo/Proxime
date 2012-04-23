#!/bin/sh
#PBS -N Chaos
#PBS -M angelo@if.uff.br
#PBS -m bae
#PBS -l nice=16,nodes=6:ppn=4,walltime=1068:00:00 
cd $PBS_O_WORKDIR

date
export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
mpirun -machinefile $PBS_NODEFILE -np $NPROCS /nfs/fiscomp/angelo/Prog_Caos/Caos/chaos
date
