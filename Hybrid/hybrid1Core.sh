#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 15:00
#SBATCH -N 1 # nodes to reserve 
#SBATCH -p node

# Number of MPI tasks
#SBATCH -n 1  # 8 -> 4 nodes
#
# Number of cores per task (proccess)
#SBATCH -c 1  # 10 cores in each processor

module load gcc openmpi
export OMP_NUM_THREADS=1


echo 160 particles, $2 processes, $3 threads 
mpirun ./nBody 160 200 1e-3 0 -x OMP_NUM_THREADS

echo 320 particles, $2 processes, $3 threads 
mpirun ./nBody 320 200 1e-3 0 -x OMP_NUM_THREADS

echo 640 particles, $2 processes, $3 threads 
mpirun ./nBody 640 200 1e-3 0 -x OMP_NUM_THREADS

echo 1280 particles, $2 processes, $3 threads 
mpirun ./nBody 1280 200 1e-3 0 -x OMP_NUM_THREADS

echo 2560 particles, $2 processes, $3 threads 
mpirun ./nBody 2560 200 1e-3 0 -x OMP_NUM_THREADS

echo 5120 particles, $2 processes, $3 threads 
mpirun ./nBody 5120 200 1e-3 0 -x OMP_NUM_THREADS

echo 10240 particles, $2 processes, $3 threads 
mpirun ./nBody 10240 200 1e-3 0 -x OMP_NUM_THREADS
