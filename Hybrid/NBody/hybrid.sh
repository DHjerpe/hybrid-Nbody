#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 00:20
#SBATCH -N 2
#SBATCH -p node

module load gcc openmpi

# echo $1 particles

echo 2000 particles
echo "two nodes, 1 threads on each node"
mpirun -np 2 ./galsim 2000 input_data/ellipse_N_02000.gal 200 1e-5 0 1


echo 2000 particles
echo "two nodes, 20 threads on each node" 
mpirun -np 2 ./galsim 2000 input_data/ellipse_N_02000.gal 200 1e-5 0 20


echo 2000 particles
echo "two nodes, 40 processes total"
mpirun -np 40 ./galsim 2000 input_data/ellipse_N_02000.gal 200 1e-5 0 1


# echo 4 processors
# mpirun -np 4 --map-by core ./wave $1

# echo 9 processors
# mpirun -np 9 --map-by core ./wave $1

# echo 12 processors
# mpirun -np 12 --map-by core ./wave $1

# echo 16 processors
# mpirun -np 16 --map-by core ./wave $1

# echo 25 processors
# mpirun -np 25 --map-by core ./wave $1

# echo 36 processors
# mpirun -np 36 --map-by core ./wave $1

# echo 40 processors
# mpirun -np 40 --map-by core ./wave $1
