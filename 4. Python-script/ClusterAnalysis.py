# -*- coding: utf-8 -*-
"""
Created on Mon May  8 11:07:51 2023

@author: Yang Wang

For obtaining the clusters of graphene

"""

#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from ovito.io import import_file, export_file
from ovito.modifiers import ClusterAnalysisModifier, ExpressionSelectionModifier

matplotlib.use('Agg')

plt.rcParams["font.family"] = "Times New Roman" 
cmap = plt.get_cmap('hot')
cutoff = 12.11

gd = [GraftDensity]
divide = 1000000
ts = 4 # fs
for i in range(len(gd)):  
    color = cmap(i/len(gd))
    #%%
    pipeline1 = import_file("/scratch/p309582/Graphene-grafting-PMMA/MD/GraftingDensity/{:.3f}/All.dump".format(gd[i]))
    pipeline1.modifiers.append(ExpressionSelectionModifier(expression = 'ParticleType==1'))
    step1       = []
    cluster1    = []
    for f in range(pipeline1.source.num_frames):
        pipeline1.modifiers.append(ClusterAnalysisModifier(
            cutoff=cutoff,
            only_selected=True,
            sort_by_size=True, 
            compute_com=True,
            compute_gyration=True,
            unwrap_particles=True))
        export_file(pipeline1, "/scratch/p309582/Graphene-grafting-PMMA/MD/GraftingDensity/Cluster/Python/cluster-{:.3f}-stage1-{}.txt".format(gd[i],f), 'txt/table', key='clusters',
                    multiple_frames=True, start_frame=f, end_frame=f)
        step1.append(50000*ts*f)
        data1 = pipeline1.compute(frame=f)
        cluster_table = data1.tables['clusters']
        cluster1.append(max(cluster_table['Cluster Identifier'][...]))

    #%%
    fig = plt.figure(figsize=(6, 4.5))
    step    = np.concatenate([step1])
    cluster = np.concatenate([cluster1])
    plt.xlabel('Time (ns)',fontdict={'size':20})
    plt.ylabel('Number of Clusters',fontdict={'size':20})
    plt.xticks(size=20)   
    plt.yticks(size=20)   
    plt.rcParams['xtick.direction']='in'
    plt.rcParams['ytick.direction']='in'
    plt.ylim(0,51)
    plt.plot(step/divide,cluster, color=color)
    plt.scatter(step/divide,cluster, color=color,label='gd = {:.3f}'.format(gd[i]))
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}
    legend = plt.legend(prop=font,frameon=False)
    plt.savefig("/scratch/p309582/Graphene-grafting-PMMA/MD/GraftingDensity/Cluster/Python/ClusterAnalysis-{:.3f}gd.jpg".format(gd[i]),bbox_inches='tight',dpi=800)
    with open  ('/scratch/p309582/Graphene-grafting-PMMA/MD/GraftingDensity/Cluster/Python/ClusterAnalysis-{:.3f}gd.txt'.format(gd[i]),'w')as Myfile:
        for x in range(len(step)):
            Myfile.write('{} {}\n'.format(step[x]/divide,cluster[x]))


