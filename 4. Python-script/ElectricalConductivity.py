# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 22:39:40 2024

@author: Yang Wang
"""

#%%
import numpy as np
import math
import matplotlib.pyplot as plt
import MDAnalysis as mda
from scipy.spatial import distance_matrix 
from scipy.spatial.distance import cdist
import networkx as nx
from scipy.spatial import distance
from mpl_toolkits.mplot3d import Axes3D
import scienceplots
from scipy.linalg import pinv
from collections import defaultdict, deque
from itertools import combinations

#plt.style.use('ieee')
plt.rcParams['font.family'] = 'Arial'
#cmap = plt.get_cmap('Blues')

#%%
GD = [GraftingDensity]
distance_threshold = 20
num_graphene       = 50

font_axis = 34
font_label = 38
bead_size = 0.4
com_size = 10

#%%
sigma_ave = []
sigma_error = []
edges_ave = []
edges_error = []
sheet   = []
for g in range(len(GD)):
    u1 = mda.Universe(r"E:\Ethan\Graphene-Grafting-PMMA\MD\GraftingLength\L={}\Final.data".format(GD[g]),
                      r"E:\Ethan\Graphene-Grafting-PMMA\MD\GraftingLength\L={}\wrap.dcd".format(GD[g]))
    u2 = mda.Universe(r"E:\Ethan\Graphene-Grafting-PMMA\MD\GraftingLength\L={}\Final.data".format(GD[g]),
                      r"E:\Ethan\Graphene-Grafting-PMMA\MD\GraftingLength\L={}\unwrap.dcd".format(GD[g]))
    print(u1.trajectory)
    
    sigma_all = [] 
    edge_all = []
    sheet_all = []
    for ts in u1.trajectory[0:500:100]:
        data1 = u1.select_atoms('type 1').positions # Wrap
        data2 = u2.select_atoms('type 1').positions # Wrap
        left_wall_threshold  = (min(u1.atoms.positions[:,0])+min(u1.atoms.positions[:,1])+min(u1.atoms.positions[:,2]))/3
        right_wall_threshold = (max(u1.atoms.positions[:,0])+max(u1.atoms.positions[:,1])+max(u1.atoms.positions[:,2]))/3
        box_size_threshold   = right_wall_threshold - left_wall_threshold - 2

        box_dims = np.array([left_wall_threshold, right_wall_threshold,
                             left_wall_threshold, right_wall_threshold,
                             left_wall_threshold, right_wall_threshold])

        graphene_data1 = np.array(data1.reshape(50,370,3))  # Wrap
        graphene_data2 = np.array(data2.reshape(50,370,3))  # Unwrap
        
        box_size = right_wall_threshold - left_wall_threshold
  
        G = nx.Graph()
        for i in range(num_graphene):
            for j in range(num_graphene):
                if np.any(cdist(graphene_data1[i], graphene_data1[j]) <= distance_threshold):
                    G.add_edge(i, j)
                    
        components = list(nx.connected_components(G))
        print(len(components),len(components[0]))
        
        def check_through_box(coords):
            coords = np.vstack(coords)
            max_coords = coords.max(axis=0)
            min_coords = coords.min(axis=0)
            coord_diff_x = max_coords[0] - min_coords[0]
            coord_diff_y = max_coords[1] - min_coords[1]
            coord_diff_z = max_coords[2] - min_coords[2]
            return coord_diff_x > box_size_threshold and coord_diff_y > box_size_threshold and coord_diff_z > box_size_threshold
        
        results = []
        for component in components:
            component_coords1 = [graphene_data1[i] for i in component] 
            component_coords2 = [graphene_data2[i] for i in component]
            if check_through_box(component_coords1) and check_through_box(component_coords2):
                through_box = check_through_box(component_coords1) and check_through_box(component_coords2)
                results.append({
                    "graphene_ids": list(component),
                    "through_box": through_box})
            else:
                results.append({
                    "graphene_ids": list(component),
                    "through_box": False })
        excluded_numbers = []
        for result in results:
            print(f"Graphene IDs: {result['graphene_ids']}, Through Box: {result['through_box']}")
            if result['through_box'] == False:
                excluded_numbers.extend(result['graphene_ids'])
        #excluded_numbers.extend([2, 36, 12, 48, 24])
        print(excluded_numbers)

        #%%
        # Selected graphene sheet IDs for analysis
        selected_ids = [i for i in range(0, 50) if i not in excluded_numbers]
        
        # Calculate the center of mass for each graphene sheet
        center_of_mass = np.mean(graphene_data2, axis=1)
        
        # Filter the graphene data and center of mass based on the selected IDs
        selected_graphene_data = graphene_data1[selected_ids]
        selected_center_of_mass = center_of_mass[selected_ids]

        # Adjust center of mass to fit within the box dimensions
        for i in range(3):  # Loop over x, y, z dimensions
            if box_dims[2*i] != box_dims[2*i+1]:  # Check if the boundaries are different (not fixed)
                selected_center_of_mass[:, i] = np.where(selected_center_of_mass[:, i] < box_dims[2*i], 
                                                selected_center_of_mass[:, i] + (box_dims[2*i+1] - box_dims[2*i]),
                                                selected_center_of_mass[:, i])
                selected_center_of_mass[:, i] = np.where(selected_center_of_mass[:, i] > box_dims[2*i+1], 
                                                selected_center_of_mass[:, i] - (box_dims[2*i+1] - box_dims[2*i]),
                                                selected_center_of_mass[:, i])

        # Function to calculate the minimum distance between beads of two graphene sheets
        def min_distance_between_sheets(sheet1, sheet2):
            return np.min(distance.cdist(sheet1, sheet2, 'euclidean'))
        
        # Create an array to store the pairs of sheets that are connected
        connected_pairs = []
        # Iterate over all pairs of graphene sheets to check distances
        for i in range(len(selected_graphene_data)):
            for j in range(i + 1, len(selected_graphene_data)):
                if min_distance_between_sheets(selected_graphene_data[i], selected_graphene_data[j]) < distance_threshold:
                    min_dist = min_distance_between_sheets(selected_graphene_data[i], selected_graphene_data[j])
                    connected_pairs.append((i, j, min_dist))
        #connected_pairs.append((32,49,19.8))
        #connected_pairs.append((21,34,19.9))
        #connected_pairs.append((25,39,19.9))
        
        len(connected_pairs)

        #%%
        # Adjust nodes to start from 0 for Python indexing and extract distances
        connected_edges = [(u, v) for u, v, _ in connected_pairs]  # Remove the third element 'd'
        # Create a graph
        G = nx.Graph()
        G.add_edges_from(connected_edges)

        # Calculate adjacency matrix A
        nodelist = sorted(G.nodes())
        A = nx.adjacency_matrix(G, nodelist=nodelist).toarray()
        
        # Calculate degree matrix D
        degrees = dict(G.degree())
        D = np.diag([degrees[node] for node in nodelist])
        
        # Calculate Laplacian matrix L
        L = nx.laplacian_matrix(G, nodelist=nodelist).toarray()
        
        # Calculate the pseudoinverse of the Laplacian matrix
        L_pinv = pinv(L)
        
        # Function to calculate effective resistance between specified nodes
        def effective_resistance(L_pinv, start_node, end_node):
            R_start_end = L_pinv[start_node, start_node] + L_pinv[end_node, end_node] - 2 * L_pinv[start_node, end_node]
            return R_start_end

        # Loop over all unique pairs of nodes and calculate conductivities
        conductivities = []
        for start_node, end_node in combinations(range(len(nodelist)), 2):
            resistance = effective_resistance(L_pinv, start_node, end_node)
            conductivity = 1 / resistance
            conductivities.append(conductivity)
            #print(f"Conductivity between node {nodelist[start_node]} and node {nodelist[end_node]}: {conductivity:.4f} siemens")
            
        print("Overall Electrical Conductivity of the network: {}".format(sum(conductivities)))
        sigma_all.append(sum(conductivities))
        edge_all.append(len(connected_pairs))
        sheet_all.append(len(selected_ids))
    sigma_ave.append(np.mean(sigma_all))
    sigma_error.append(np.std(sigma_all))
    edges_ave.append(np.mean(edge_all))
    edges_error.append(np.std(edge_all))
    sheet.append(np.mean(sheet_all))
