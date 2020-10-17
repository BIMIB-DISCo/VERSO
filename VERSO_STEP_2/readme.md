# SUMMARY
VERSO STEP #2 processes the variant frequency (VF) profiles of groups of samples with the same clonal genotype (identified via VERSO STEP #1), in order to characterize their intra-host genomic composition and visualize it on a low-dimensional space. 
This step requires VF profiles generated from raw sequencing data via variant calling and the prior execution of STEP #1.

VERSO STEP #2 is provided as Python script (filename: "VERSO_STEP_2.py") and employs standard libraries included in the SCANPY suite. 

# REQUIREMENTS
Please install the following libraries.

* [Python 3.6.x] Follow the instructions at: https://www.python.org/downloads/

* [Scanpy 1.6] Follow the instructions at: https://scanpy.readthedocs.io/en/stable/installation.html


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
* **File name: _"samples_info_matrix.csv"_**

* File type: csv (comma separated value)

* Content: a csv file including the key information about samples, with: 

	* _m_ (samples) rows 
	
	* _k_ (attributes) columns

* The entry in position 1,1 (1st row, 1st column) must be left empty. 

	The first column (from position 2 to position _m_+1) must include sample IDs.

	The first row (from position 2 to position _k_+1) must include the attribute name.

	Each remaining entry in position _i,j_ includes the value of a specific attribute _j_ in sample _i_.  


**Warning** : an attribute named "Genotype" (i.e., the number of the clonal genotype of each sample returned by VERSO STEP #1) must necessarily included to perform STEP #2.


## FILE C
* **File name: _"variants_info_matrix.csv"_**

* File type: csv (comma separated value)

* Content: Variants info matrix, with _n_ rows (variants) and _z_ columns (attributes), each entry includes the value of a specific attribute for any variant to be considered in the analysis. 
 
* The entry in position 1,1 (1st row, 1st column) must be left empty. 

 	The first column (from position 2 to position _n_+1) must include variant IDs.

	The first row (from position 2 to position _z_+1) must include the attribute name.

	Each remaining entry in position _i,j_ includes the value of a specific attribute _j_ in variant _i_.  

## FILE D
* **File name: _"configuration_VERSO.txt"_**

* File type: txt (textual)

* Content: file including the optional parameters of VERSO STEP #2. Im detail: 

* The first row must include the headers: 
	* Parameter 
	* Value

* From the second row on, please include the value of the following parameters: 
	* [knn]: number of k nearest neighbours (default = 10)
	* [n_pca]: number of principal components (default = 10)
	* [distance_metric]: metrics to be employed (see SCANPY documentation) (default = 'braycurtis')
	* [leiden_res]: resolution of the Leiden clustering algorighm (default = 1)
	* [a_UMAP]: parameter a of the UMAP plot (default = 0.01)
	* [b_UMAP]: parameter b of the UMAP plot (default = 1)
	* [spread_UMAP]: parameter spread of the UMAP plot (default = 1)
	* [min_dist_UMAP]: parameter min_dist of the UMAP plot (default = 2)
	



# RUNNING
Launch the Python script from the terminal, with the following command: 

	python VERSO_STEP_2.py

# OUTPUTS
VERSO STEP #2 returns as output:

* 1) the SVG images including the UMAP plots related to the distinct clonal genotypes included in the datasets. The file names are numbered according to the clonal genotype ID: _C01.svg_, _C02.svg_, etc. 

* 2) the metadata for each clonal genotype in folders names as: OUTPUT_C01, OUTPUT_C02, etc. 
