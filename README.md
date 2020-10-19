Viral Evolution ReconStructiOn (VERSO)
================

[![Actions Status](https://github.com/BIMIB-DISCo/VERSO/workflows/check-master/badge.svg)](https://github.com/BIMIB-DISCo/VERSO/actions?query=workflow%3Acheck-master)
[![Actions Status](https://github.com/BIMIB-DISCo/VERSO/workflows/check-development/badge.svg)](https://github.com/BIMIB-DISCo/VERSO/actions?query=workflow%3Acheck-development)

In this repository we provide an implementation of the *VERSO* (Viral Evolution ReconStructiOn) framework. 

*VERSO* is an algorithmic framework that processes variants profiles from viral samples, to produce phylogenetic models of viral evolution from clonal variants and to subsequently quantify the intra-host genomic diversity of samples. VERSO includes two separate and subsequent steps, which are perfomed via the scripts provided in this repository. 

The manuscript of the method is available as a preprint on bioRxiv at: https://www.biorxiv.org/content/10.1101/2020.04.22.044404v2

![Imgur Graphical Abstract](https://github.com/BIMIB-DISCo/VERSO/blob/VERSO/graphical_abstract.png) 

# VERSO STEP #1: robust phylogenomic inference from clonal variant profiles
VERSO STEP #1 employs a probabilistic maximum-likelihood framework for the reconstruction of robust phylogenetic trees from binarized mutational profiles of clonal variants. 

VERSO STEP #1 takes as input:  a n (samples) × m (variants) binary mutational profile matrix, as defined on clonal SNVs only.  In this case, an entry in a given sample is equal to 1 (present) if the VF is larger than a certain threshold (in our analyses, equal to 90%), it is equal to 0 if lower than a distinct threshold (in our analyses, equal to 5%), and is considered as missing (NA) in the other cases.

VERSO outputs: (i) a phylogenetic model where each leaf correspond to a sample, whereas internal nodes correspond to accumulating clonal variants (in Newick format), (ii) the corrected clonal genotype of each sample.

VERSO is provided as an R open source tool. 

## RUNNING
To run the tools for VERSO STEP #1 please go to the folder: https://github.com/BIMIB-DISCo/VERSO/tree/VERSO/VERSO_STEP_1 and refer to the related readme, which include instructions and guidelines. 


# VERSO STEP #2: characterization of intra-host genomic diversity
In the second step, VERSO takes into account the variant frequency profiles of groups of samples with the same clonal genotype (identified via VERSO STEP #1), in order to characterize their intra-host genomic diversity and visualize it on a low-dimensional space. 

VERSO STEP #2 takes as input a n (samples) × m (variants) variant frequency (VF) profile matrix, in which each entry includes the VF ∈ (0,1) of a given mutation in a certain sample, after filtering out: (i) the clonal variants employed in STEP #1 and (ii) the minor variants possibly involved in homoplasies.

As output, VERSO STEP #2 delivers both the partitioning of samples in homogeneous clusters and the visualizationin a low-dimensional space, also allowing to label samples according to other covariates, such as, e.g., collection date or geographical location. 

VERSO STEP #2 is provided as a Python script which employs the SCANPY suite of tools: https://scanpy.readthedocs.io/en/stable/

## RUNNING
To run the tools for VERSO STEP #2 please go to the folder: https://github.com/BIMIB-DISCo/VERSO/tree/VERSO/VERSO_STEP_2 and refer to the related readme, which include instructions and guidelines. 


Please feel free to contact us if you have problems running our tool at daniele.ramazzotti1@gmail.com or d.maspero@campus.unimib.it. 
