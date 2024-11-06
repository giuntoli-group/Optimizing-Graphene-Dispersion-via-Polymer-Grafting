# Optimizing-Graphene-Dispersion-via-Polymer-Grafting

This repository contains LAMMPS codes used in our paper to simulate the graphene-p(MMA) nanocomposites with varying grafted chain lengths and grafting densities. 

You can read our paper [here](https://pubs.acs.org/doi/************). 

Here, we have three folders:

1. Grafted-Chain-Length (7 subfolders with different grafted chain lengths.) 
2. Grafting-Density (7 subfolders with different grafting densities.)
3. Tensile (1 subfolder for tensile test in x direction of graphene/p(MMA) system with grafted chain length of 20.)
4. Python-script (DispersionState.py, ClusterAnalysis.py, LAMMPS_Data_Morphologies.py, ElectricalConductivity.py)

To run the scripts and obtain the equilibrated system, execute:

lmp -in Graphene-PMMA.in

To conduct the tensile deformation, execute the following command based on Final.data generated in the previous simulation:

lmp -in tensile.in
