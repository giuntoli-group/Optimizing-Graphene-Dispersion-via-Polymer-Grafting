# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 11:50:25 2023

@author: Yang Wang

To count the number of graphene CG beads in Aggregated, Intercalated, and Unbound morphologies.

"""

import numpy as np
import math
import matplotlib.pyplot as plt
import MDAnalysis as mda

u = mda.Universe(r"G:\Ethan\Graphene-Grafting-PMMA\MD\Model.data",
                 r"G:\Ethan\Graphene-Grafting-PMMA\MD\GraftingDensity\0.000\Equilibration-graphene.lammpstrj",format='LAMMPSDUMP')
print(u.trajectory)

interval  = 30

ids = 0 
for ts in u.trajectory[0:len(u.trajectory):50]:
    print(ts.frame)
    pos = u.select_atoms('type 1').positions
    len(pos)
    data = pos.reshape(50, 370, 3)
    aggregated = []
    intercalated = []
    seperation = []
    for i in range(len(data)):
        group_a = data[i]
        group_b = np.delete(data, i, axis=0).reshape((50-1)*370, 3)
        count1 = 0
        count2 = 0
        count3 = 0
        distances = np.sqrt(np.sum((group_a[:, np.newaxis] - group_b) ** 2, axis=2))
        for row in distances:
            if np.any(row <= 8):
                count1 += 1
                continue 
            if np.any((row > 8) & (row <= 18)):
                count2 += 1
                continue 
            if np.any(row > 18):
                count3 += 1
        aggregated.append(count1/370)
        intercalated.append(count2/370)
        seperation.append(count3/370)

    with open (r'G:\Ethan\Graphene-Grafting-PMMA\MD\dispersion-{}.txt'.format(ids),'w') as file:
        for s in range(len(aggregated)):    
            file.write('{} {} {}\n'.format(aggregated[s],intercalated[s],seperation[s]))
    ids += 1
    
u = mda.Universe(r"G:\Ethan\Graphene-Grafting-PMMA\MD\GraftingDensity\0.000\Equil.data",
                 r"G:\Ethan\Graphene-Grafting-PMMA\MD\GraftingDensity\0.000\Annealing.lammpstrj",format='LAMMPSDUMP')
print(u.trajectory)

interval  = 30

for ts in u.trajectory[0:len(u.trajectory):50]:
    print(ts.frame)
    pos = u.select_atoms('type 1').positions
    len(pos)
    data = pos.reshape(50, 370, 3)
    aggregated = []
    intercalated = []
    seperation = []
    for i in range(len(data)):
        group_a = data[i]
        group_b = np.delete(data, i, axis=0).reshape((50-1)*370, 3)
        count1 = 0 
        count2 = 0 
        count3 = 0 
        distances = np.sqrt(np.sum((group_a[:, np.newaxis] - group_b) ** 2, axis=2))
        for row in distances:
            if np.any(row <= 8):
                count1 += 1
                continue
            if np.any((row > 8) & (row <= 18)):
                count2 += 1
                continue
            if np.any(row > 18):
                count3 += 1
        aggregated.append(count1/370)
        intercalated.append(count2/370)
        seperation.append(count3/370)

    with open (r'G:\Ethan\Graphene-Grafting-PMMA\MD\dispersion-{}.txt'.format(ids),'w') as file:
        for s in range(len(aggregated)):    
            file.write('{} {} {}\n'.format(aggregated[s],intercalated[s],seperation[s]))
    ids += 1













