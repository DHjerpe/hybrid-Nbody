#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 00:20
#SBATCH -N 2
#SBATCH -p node

# module load gcc openmpi
export OMP_NUM_THREADS=$2


echo 100 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 100 200 1e-3 0 -x OMP_NUM_THREADS

echo 200 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 200 200 1e-3 0 -x OMP_NUM_THREADS

echo 300 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 300 200 1e-3 0 -x OMP_NUM_THREADS

echo 400 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 400 200 1e-3 0 -x OMP_NUM_THREADS

echo 500 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 500 200 1e-3 0 -x OMP_NUM_THREADS

echo 600 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 600 200 1e-3 0 -x OMP_NUM_THREADS

echo 1000 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 1000 200 1e-3 0 -x OMP_NUM_THREADS

echo 2000 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 2000 200 1e-3 0 -x OMP_NUM_THREADS

echo 4000 particles, $1 processes, $2 threads 
mpirun -np 2 ./nBody 4000 200 1e-3 0 -x OMP_NUM_THREADS


