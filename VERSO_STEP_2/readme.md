# SUMMARY
VERSO STEP #2 processes the variant frequency (VF) profiles of groups of samples with the same clonal genotype (identified via VERSO STEP #1), in order to characterize their intra-host genomic composition and visualize it on a low-dimensional space. 
This step requires VF profiles generated from raw sequencing data via variant calling and the prior execution of STEP $\#1$.

VERSO STEP #2 is provided as Python script (filename: "VERSO_STEP_2.py") and employs standard libraries included in the SCANPY suite. 

# REQUIRIMENTS
Please install the following libraries.

(1) [Python 3.6.x] Follow the instructions at: https://www.python.org/downloads/

(2) [Scanpy 1.6] https://scanpy.readthedocs.io/en/stable/installation.html


# INPUTS
  
VERSO STEP #2 requires *4* input files which must be positioned in the folder of the script and detailed as follows.

**FILE A)**
File name: "VF_matrix.csv"

File type: csv (comma separated value)

Content: a csv file including the variant Frequency (VF) Matrix, with:

	n (variants) + 1 rows  
	
	m (samples) + 1 columns 
	

The entry in position 1,1 (1st row, 1st column) must be left empty. 

The first row (from position 2 to position m+1) must include sample IDs.

The first column (from position 2 to position n+1) must include variant IDs.

Each remaining entry in position i,j includes the VF of variant i in sample j. 


Warning: VERSO STEP #2 does not process or impute missing values (e.g., NA or NaN).

Therefore, a pre-processing step must be executed by the user to impute possible missing values in the VF_matrix.csv file.

####################

**FILE B)**
File name: "samples_info_matrix_filename.csv"

File type: csv (comma separated value)

Content: a csv file including the key information about samples, with: 

	m (samples) rows 
	
	k (attributes) columns

The entry in position 1,1 (1st row, 1st column) must be left empty. 

The first column (from position 2 to position m+1) must include sample IDs.

The first row (from position 2 to position k+1) must include the attribute name.

Each remaining entry in position i,j includes the value of a specific attribute j in sample i.  


Warning: an attribute named "Genotype" (i.e., the number of the clonal genotype of each sample returned by VERSO STEP #1) must necessarily included to perform STEP #2.

####################

**FILE C)**
File name: "variants_info_matrix_filename.csv"

File type: csv (comma separated value)

Content: Variants info matrix, with n rows (variants) and z columns (attributes), each entry includes the value of a 
    specific attribute for any variant to be considered in the analysis. 
        

**FILE D)**
File name: configuration_VERSO.txt"

File type: txt (textual)

Content: file including the optional parameters of VERSO STEP #2


# RUNNING
Launch the Python script from the terminal, with the following command: python VERSO_STEP_2.py

# OUTPUTS
VERSO STEP #2 returns as output:

1) the SVG images including the UMAP plots related to the distinct clonal genotypes included in the datasets. The file names are numbered according to the clonal genotype ID: C01.svg, C02.svg, etc. 

2) the metadata for each clonal genotype in folders names as: OUTPUT_C01, OUTPUT_C02, etc. 
