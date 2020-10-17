# No actions are required at this step, please just run the notebook 
#from IPython.core.interactiveshell import InteractiveShell
#InteractiveShell.ast_node_interactivity = "all"

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


 # Please LOAD the INPUT FILES with the following specifications: 

# A) Variant Frequency (VF) Matrix, with n rows (variants) and m columns (samples), each entry includes the VF 
#    for that variant in that sample. Pre-processing step might be required to impute NA values (e.g., due to
#    low coverage). In the example analyses, NA's are imputed with 0's. 
#    file format -> comma separated value (csv) file: input_data_matrix_filename = "xxx.csv"

# B) Samples info matrix, with m rows (samples) and k columns (attributes), each entry includes the value of a 
#    specific attribute, e.g., country, date, clade, etc., for any given sample. 
#    file format -> comma separated value (csv) file: samples_info_matrix_filename = "yyy.csv"

# C) Variants info matrix, with n rows (variants) and z columns (attributes), each entry includes the value of a 
#    specific attribute for any variant to be considered in the analysis. 
#    file format -> comma separated value (csv) file: variants_info_matrix_filename = "zzz.csv"

#(A)
input_data_matrix_filename = "VF_matrix.csv"

#(B)
samples_info_matrix_filename = "samples_info_matrix.txt"

#(C)
variants_info_matrix_filename = "variants_info_matrix.txt"



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

adata_all.obs['Genotype'] = 'NA'
adata_all.obs['Outlier'] = 'NA'
adata_all.obs['Contact_group'] = 'NA'
adata_all.obs['Week'] = 'NA'


adata_all.obs['Genotype'][:] = assignments['Genotype']
adata_all.obs['Outlier'][:] = assignments['outlier']
adata_all.obs['Contact_group'][:] = assignments['Epi']
adata_all.obs['Week'][:] = assignments['week']

adata_all


#selecting only the variants present in the dataset and those observed in at most 1 clonal genotype, 
#excluding the clonal variants employed in the reconstruction of the phylogenetic tree via VERSO STEP #1
list_variants_no_homoplasy = metadata[metadata['nCluster']<2].index
list_variants_no_clonal = metadata[metadata['inference']==0].index
index_final_variants = list_variants_no_homoplasy.intersection(list_variants_no_clonal)

adata_all_subset_variants = adata_all[:, index_final_variants]
adata_all_subset_variants

#selecting only the samples that are NOT outliers (with <100 minor variants)
adata_all_subset_variants_samples = adata_all_subset_variants[adata_all_subset_variants.obs['Outlier']==0].copy()
adata_all_subset_variants_samples


#performing VERSO STEP #2 on all clonal genotypes
how_many_clusters = len(np.unique(adata_all_subset_variants_samples.obs['Genotype']))
how_many_clusters
for x in range(how_many_clusters):

    y = x + 1
    if y<10:
        print('C'+ str(0)+str(y))
        temp = ('C'+ str(0)+str(y))
        temp_title = ('Dataset #1 - Genotype '+ str(0)+str(y))
    else: 
        print('C'+str(y))
        temp = ('C'+str(y))
        temp_title = ('Dataset #1 - Genotype '+str(y))
        
    adata_all_cluster = adata_all_subset_variants_samples[adata_all_subset_variants_samples.obs['Genotype'] == temp].copy()
    
    how_many_samples = len(adata_all_cluster.obs)
    temp_title = temp_title + ' (n = '+ str(how_many_samples) +')'
    
    how_many_variables = len(adata_all_cluster.var)
    
    if how_many_samples > 20:
        number_of_neighbours_to_be_used = 10
        number_of_PCA = 10
        
        sc.pp.pca(adata_all_cluster, n_comps=10,zero_center = 'False',random_state=1)       
        sc.pp.neighbors(adata_all_cluster, knn = 'TRUE', n_pcs = number_of_PCA, 
                        n_neighbors=number_of_neighbours_to_be_used,random_state=1,
                        metric = 'braycurtis', method = "umap")
        sc.tl.leiden(adata_all_cluster, key_added='groups', random_state=1,resolution = 1)

        sc.tl.umap(adata_all_cluster, random_state=1, a = 0.01, b = 1, spread = 1, min_dist = 2)

        temp_file = temp+'.svg'
        temp_output = 'output'+temp+'.csv'
        temp_file_distances = temp+'_distances.txt'

        if how_many_samples > 500:
            size_nodes = 40
        if how_many_samples <= 500:
            size_nodes = 100
        sc.pl.umap(adata_all_cluster, color=['groups'], legend_loc='right margin', edges = 'TRUE', edges_width = 0.2, size = size_nodes, title = temp_title , show=False)       
        #to visualize the contact information on the UMAP plots, uncomment the following line
        # sc.pl.umap(adata_all_cluster, color=['Contact_group'], legend_loc='right margin', edges = 'TRUE', edges_width = 0.2, size = 50,  title = temp ,show=False)
        
        plt.savefig(temp_file,bbox_inches='tight')
        adata_all_cluster.write_csvs(temp_output)
        
        #saving the distances of the knn graph
        c = sp.sparse.csc_matrix.todense(adata_all_cluster.obsp['distances'])      
        np.savetxt(temp_file_distances, c)

