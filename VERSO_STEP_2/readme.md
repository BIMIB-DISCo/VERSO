# SUMMARY
VERSO STEP #2 processes the variant frequency (VF) profiles of groups of samples with the same clonal genotype (identified via VERSO STEP #1), in order to characterize their intra-host genomic composition and visualize it on a low-dimensional space. 
This step requires VF profiles generated from raw sequencing data via variant calling and the prior execution of STEP #1.

VERSO STEP #2 is provided as Python script (filename: "VERSO_STEP_2.py") and employs standard libraries included in the SCANPY suite.  

# REQUIREMENTS
Please install the following libraries.

* [Python 3.6.x] Follow the instructions at: https://www.python.org/downloads/

* [SCANPY 1.6] Run this command via PIP to install the proper version:

	<pre><code>pip install scanpy==1.6</code></pre>
Further details at: https://scanpy.readthedocs.io/en/stable/installation.html
# INPUTS
  
VERSO STEP #2 requires **4** input files which must be positioned in the folder of the script and detailed as follows.

## FILE A
* **File name: _"VF_matrix.csv"_**

* File type: csv (comma separated value)

* Content: a csv file including the variant Frequency (VF) Matrix, with:
	* _n_ (variants) + 1 rows  
	* _m_ (samples) + 1 columns 
	
* The entry in position 1,1 (1st row, 1st column) must be left empty. 

	The first row (from position 2 to position _m_+1) must include sample IDs.

	The first column (from position 2 to position _n_+1) must include variant IDs.

	Each remaining entry in position _i,j_ includes the VF of variant i in sample j. 


**Warning** : VERSO STEP #2 does not process or impute missing values (e.g., _NA_ or _NaN_).

Therefore, a pre-processing step must be executed by the user to impute possible missing values in the _"VF_matrix.csv file"_.

## FILE B
* **File name: _"samples_info_matrix.txt"_**

* File type: txt (tab-delimited text file)

* Content: a txt file including the key information about samples, with: 

	* _m_ (samples) rows 
	
	* _k_ (attributes) columns

* The entry in position 1,1 (1st row, 1st column) must be left empty. 

	The first column (from position 2 to position _m_+1) must include sample IDs.

	The first row (from position 2 to position _k_+1) must include the attribute labels.
		The first attribute must be: "Genotype" (clonal genotype of each sample, retrieved from VERSO STEP #1, numbered from 1). 
		The second attribute must be: "Selected" (if the sample is selected this value has to be **1**)

	Each remaining entry in position _i,j_ includes the value of a specific attribute _j_ in sample _i_.  

**Warning** : an attribute named "Genotype" (i.e., the number of the clonal genotype of each sample returned by VERSO STEP #1) must necessarily included to perform STEP #2.

## FILE C
* **File name: _"variants_info_matrix.txt"_**

* File type: txt (tab-delimited text file)

* Content: Variants info matrix, with _n_ rows (variants) and _z_ columns (attributes), each entry includes the value of a specific attribute for any variant to be considered in the analysis. 
 
* The entry in position 1,1 (1st row, 1st column) must be left empty. 

 	The first column (from position 2 to position _n_+1) must include variant IDs.

	The first row (from position 2 to position _z_+1) must include the attribute labels.
		The first attribute must be: "Selected" (if the variant is included in the analysis this value has to be **1**). 

	Each remaining entry in position _i,j_ includes the value of a specific attribute _j_ in variant _i_.  

## FILE D
* **File name: _"configuration_VERSO.txt"_**

* File type: txt (tab-delimited text file)

* Content: file including the optional parameters of VERSO STEP #2. Im detail: 

* The first row must include the headers: 
	* Parameter 
	* Value

* From the second row on, please include the value of the following parameters: 
	* [sample_filter] yes (use only samples with value 1 in "Selected" attribute) / no (use all samples)
	* [variant_filter] yes (use only variants with value 1 in "Selected" attribute) / no (use all variants)
	* [knn]: number of k nearest neighbours
	* [n_pca]: number of principal components to be used 
	* [distance_metric]: metric to be employed in the knn (e.g., 'braycurtis') (see SCANPY documentation: https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.neighbors.html)
	* [leiden_res]: resolution of the Leiden clustering algorithm on the knn graph
	* [a_UMAP]: parameter a of the UMAP plot
	* [b_UMAP]: parameter b of the UMAP plot
	* [spread_UMAP]: parameter spread of the UMAP plot
	* [min_dist_UMAP]: parameter min_dist of the UMAP plot

# RUNNING
Launch the Python script from the terminal, with the following command: 

	python VERSO_STEP_2.py

# OUTPUTS
VERSO STEP #2 returns as output:

* 1) the SVG images including the UMAP plots related to the distinct clonal genotypes included in the datasets. 
	** A single SVG is produced for each sample attribute (excluded "Genotype" and "Selected"). 
	
	The file names are numbered according to the clonal genotype ID and the Attribute: G01_Attribute_1.svg,  G02_Attribute_1.svg, etc. 
* 2) the distance among samples, numbered according to the clonal genotype ID: G01_distances.txt, G01_distances.txt, etc.  

* 3) the metadata for each clonal genotype in folders names as: G01_OUPUTS, G02_OUPUTS, etc.txt

# EXAMPLE FILES
As an example, the GitHub folder contains sample files from Dataset #1 of the original publication: https://www.biorxiv.org/content/10.1101/2020.04.22.044404v2
