  # -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 11:16:14 2023

@author: Lenovo
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import MDAnalysis as mda

G = [0.000,0.025,0.050,0.075,0.100,0.125,0.150]
for g in range(len(G)):
    u = mda.Universe(r"E:\Ethan\Graphene-Grafting-PMMA\MD\Model.data",
                     r"E:\Ethan\Graphene-Grafting-PMMA\MD\GraftingDensity\{:.3f}\graphene-wrap-1.dump".format(G[g]),format='LAMMPSDUMP')  # 后缀
    print(u.trajectory)
    
    mass = u.select_atoms('type 1').masses[0]
    
    ids = 0 
    for ts in u.trajectory[2990:2991:1]:
        print(ts.frame)
        lx = u.dimensions[0]
        ly = u.dimensions[1]
        lz = u.dimensions[2]
        pos = u.select_atoms('type 1').positions
        len(pos)
        data = pos.reshape(50, 370, 3)
        A = []
        I = []
        S = []
        MA = []
        MI = []
        MS = []
        index = 1
        for i in range(len(data)):
            group_a = data[i]
            group_b = np.delete(data, i, axis=0).reshape((50-1)*370, 3)
            
            # Calculate distances between each coordinate in group_a and group_b
            distances = np.sqrt(np.sum((group_a[:, np.newaxis] - group_b) ** 2, axis=2))
            
            # Initialize lists to store coordinates based on conditions
            aggregated_coordinates = []
            intercalated_coordinates = []
            separation_coordinates = []
            am = []
            im = []
            sm = []
            for m, row in enumerate(distances):
                if np.any(row <= 8):
                    aggregated_coordinates.append(group_a[m])
                    am.append(index)
                
                elif np.any((row > 8) & (row <= 18)):
                    intercalated_coordinates.append(group_a[m])
                    im.append(index)
                    
                elif np.any(row > 18):
                    separation_coordinates.append(group_a[m])
                    sm.append(index)
            index += 1        
            # Append the lists of coordinates to their respective lists for each frame
            A.append(aggregated_coordinates)
            I.append(intercalated_coordinates)
            S.append(separation_coordinates)
            MA.append(am)
            MI.append(im)
            MS.append(sm)
        # ... (your previous code)
        ids += 1
        with open(r'E:\Ethan\Graphene-Grafting-PMMA\MD\GraftingDensity\DispersedData\dispersion-G-{:.3f}-{}.data'.format(G[g],ts.frame),'w')as LAMMPS:
            # First line is a comment line 
            LAMMPS.write('The random Graphene grafted PI CG system from Python\n\n')
            #----------------Header Line----------------#
            LAMMPS.write('{} atoms\n\n'.format(50*int(len(aggregated_coordinates) + len(intercalated_coordinates) + len(separation_coordinates))))
            #-----------------Types defination--------------#
            LAMMPS.write('{} atom types\n\n'.format(3))
            #--------------Specify Masses and dimensions------------------#
            LAMMPS.write('{} {} xlo xhi\n'.format(0, lx)) 
            LAMMPS.write('{} {} ylo yhi\n'.format(0, ly))
            LAMMPS.write('{} {} zlo zhi\n'.format(0, lz))
            LAMMPS.write('\nMasses\n\n')
            LAMMPS.write('{} {} \n'.format(1, mass))
            LAMMPS.write('{} {} \n'.format(2, mass))
            LAMMPS.write('{} {} \n'.format(3, mass))
        
            #-------------- Specify Atoms information ------------------#
            # Atoms section
            LAMMPS.write('\nAtoms  # full \n\n')
            # Atom_style: full----atom-Id; molecule-ID; atom-type; q; x; y; z;
        	# Write Atoms Section
            id = 1
            for l in range (len(A)):
                for n in range (len(A[l])):
                    LAMMPS.write('{} {} {} {} {} {} {} {} {} {}\n'.format(id, MA[l][n], 1, 0,
                                                                 A[l][n][0],
                                                                 A[l][n][1],
                                                                 A[l][n][2],0,0,0))   
                    id = id + 1
            for l in range (len(I)):
                for n in range (len(I[l])):
                    LAMMPS.write('{} {} {} {} {} {} {} {} {} {}\n'.format(id, MI[l][n], 2, 0,
                                                                 I[l][n][0],
                                                                 I[l][n][1],
                                                                 I[l][n][2],0,0,0))   
                    id = id + 1
            for l in range (len(S)):
                for n in range (len(S[l])):
                    LAMMPS.write('{} {} {} {} {} {} {} {} {} {}\n'.format(id, MS[l][n], 3, 0,
                                                                 S[l][n][0],
                                                                 S[l][n][1],
                                                                 S[l][n][2],0,0,0))   
                    id = id + 1
               