#VERSO_STEP_2.py
#Script to run VERSO STEP #2. Please refer to the readme at this link: https://github.com/BIMIB-DISCo/VERSO/tree/VERSO/VERSO_STEP_2
import os
os.chdir(os.getcwd())
os.getcwd()

import dill
import pickle
import gc
gc.collect()
pickle.HIGHEST_PROTOCOL

import scanpy as sc
import leidenalg
import scanpy.external as sce
import numpy as np
import scipy as sp
from scipy.stats import norm
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.axes as axxx
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import anndata as an
from scipy import sparse
from pybiomart import Dataset

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=400, dpi_save=500)

import matplotlib.cm as cm

def colormap(fig,nax):
    ax1=fig.get_axes()[nax]
    clim=ax1.get_children()[0].get_clim()
    cnorm=colors.Normalize(clim[0],clim[1])
    fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap=None), ax=ax1)
    fig.delaxes(fig.get_axes()[-2])


# Please LOAD the INPUT FILES. For the specifications refer to the readme at this link: 
#https://github.com/BIMIB-DISCo/VERSO/tree/VERSO/VERSO_STEP_2

#(A)
input_data_matrix_filename = "VF_matrix.csv"

#(B)
samples_info_matrix_filename = "samples_info_matrix.txt"

#(C)
variants_info_matrix_filename = "variants_info_matrix.txt"

#(D)
options_file = "configuration_VERSO.txt"

options = pd.read_csv(options_file, header=0, sep='\t')
options.set_index("Parameter", drop=True, inplace=True)

input_sample_filter = str(options[options.index=='sample_filter']['Value'][0])
input_variant_filter = str(options[options.index=='variant_filter']['Value'][0])
input_knn = int(options[options.index=='knn']['Value'])
input_n_pca = int(options[options.index=='n_pca']['Value'])
input_distance_metric = str(options[options.index=='distance_metric']['Value'][0])
input_leiden_res = float(options[options.index=='leiden_res']['Value'])
input_a_UMAP = float(options[options.index=='a_UMAP']['Value'])
input_b_UMAP = float(options[options.index=='b_UMAP']['Value'])
input_spread_UMAP = float(options[options.index=='spread_UMAP']['Value'])
input_min_dist_UMAP = float(options[options.index=='min_dist_UMAP']['Value'])

#defining the data structures
data_all = pd.read_csv(input_data_matrix_filename, header=0, sep=',')
assignments = pd.read_csv(samples_info_matrix_filename, header=0, sep='\t')
metadata = pd.read_csv(variants_info_matrix_filename, header=0, sep='\t')

data_all.set_index("Unnamed: 0", drop=True, inplace=True)
assignments.set_index("Unnamed: 0", drop=True, inplace=True)
metadata.set_index("Unnamed: 0", drop=True, inplace=True)

adata_all=an.AnnData(data_all.T)
adata_all.var_names_make_unique()
adata_all.obs_names_make_unique()

for i in range(len(assignments.columns)):
    temp = assignments.columns[i]
    adata_all.obs[temp] = 'NA'
    adata_all.obs[temp][:] = assignments[temp]

#selecting the variants
if input_variant_filter == 'yes':
    list_selected_variants = metadata[metadata['Selected']==1].index
    adata_all_subset_variants = adata_all[:, list_selected_variants]
if input_variant_filter != 'yes':
    adata_all_subset_variants = adata_all

#selecting the samples
if input_sample_filter == 'yes':
    adata_all_subset_variants_samples = adata_all_subset_variants[adata_all_subset_variants.obs['Selected']==1].copy()
if input_sample_filter != 'yes':
    adata_all_subset_variants_samples = adata_all_subset_variants

#performing VERSO STEP #2 on all clonal genotypes
genotype_list = (np.unique(adata_all_subset_variants_samples.obs['Genotype']))
how_many_clusters = len(genotype_list)

for x in range(how_many_clusters):
    print('G'+str(genotype_list[x]))
    temp = ('G'+str(genotype_list[x]))
    temp_title = ('Genotype '+str(genotype_list[x]))
        
    adata_all_cluster = adata_all_subset_variants_samples[adata_all_subset_variants_samples.obs['Genotype'] == genotype_list[x]].copy()   
    how_many_samples = len(adata_all_cluster.obs)
    how_many_variables = len(adata_all_cluster.var)

    temp_title = temp_title + ' (n = '+ str(how_many_samples) +')'
       
    if how_many_samples > 20:
        
        sc.pp.pca(adata_all_cluster, n_comps=10,zero_center = 'False',random_state=1)       
        sc.pp.neighbors(adata_all_cluster, knn = 'TRUE', n_pcs = input_n_pca, 
                        n_neighbors=input_knn,random_state=1,
                        metric = input_distance_metric, method = "umap")
        sc.tl.leiden(adata_all_cluster, key_added='Leiden_clusters', random_state=1,resolution = input_leiden_res)

        sc.tl.umap(adata_all_cluster, random_state=1, a = input_a_UMAP, b = input_b_UMAP, 
                   spread = input_spread_UMAP, min_dist = input_min_dist_UMAP)

        temp_file = temp+'.svg'
        temp_output = temp+'_OUTPUT.csv'
        temp_file_distances = temp+'_distances.txt'

        if how_many_samples > 500:
            size_nodes = 50
        if how_many_samples <= 500:
            size_nodes = 100
            
            
        for i in range(len(assignments.columns)):
            if assignments.columns[i] != 'Genotype':
            	if assignments.columns[i] != 'Selected':
		            temp_file = temp+'_'+assignments.columns[i]+'.svg'
		            temp_title_fin = temp_title + ' ' + assignments.columns[i]
		            sc.pl.umap(adata_all_cluster, color=([str(assignments.columns[i])]), 
		                       legend_loc='right margin', edges = 'TRUE', 
		                       edges_width = 0.2, size = size_nodes,  
		                       title = temp_title_fin , show=False)
		            plt.savefig(temp_file,bbox_inches='tight')
        
        temp_file = temp + '_Leiden_clusters.svg'
        temp_title_fin = temp_title + ' Leiden_clusters'
        sc.pl.umap(adata_all_cluster, color=('Leiden_clusters'),                    
                   legend_loc='right margin', edges = 'TRUE', 
                   edges_width = 0.2, size = size_nodes,                     
                   title = temp_title_fin , show=False)
        plt.savefig(temp_file,bbox_inches='tight')    
        
        adata_all_cluster.write_csvs(temp_output)
        
        #saving the distances of the knn graph
        c = sp.sparse.csc_matrix.todense(adata_all_cluster.obsp['distances'])      
        np.savetxt(temp_file_distances, c)