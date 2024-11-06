#!/bin/bash
#SBATCH --job-name=Tensile-Gr-PMMA
#SBATCH --time=8-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=50G
#SBATCH --partition=regular

# Load the module for LAMMPS
module load LAMMPS/23Jun2022-foss-2021b-kokkos

# Run LAMMPS
mpirun -np 10 lmp -in tensile.in
