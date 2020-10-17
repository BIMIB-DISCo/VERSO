# SUMMARY
VERSO STEP #1 processes binarized variant profiles (0 variant is absent, 1 present) for clonal variants obtained via variant calling or consensus sequences, in order to infer the maximum log-likelihood phylogenetic tree, where variants are inner nodes and samples are tips.

VERSO STEP #1 is provided as a set of R scripts (directory R) and a demo che be executed by running the script main.R. 

# REQUIREMENTS
Please install the following R libraries.

* [ape, see https://cran.r-project.org/web/packages/ape/index.html] Command:
<pre><code>if (!require("ape")) install.packages("ape")</code></pre>

* [Rfast, see https://cran.r-project.org/web/packages/Rfast/index.html] Command:
<pre><code>if (!require("Rfast")) install.packages("Rfast")</code></pre>

# INPUTS
  
VERSO STEP #1 requires as input an R matrix reporting binarized variant profiles for a set of clonal variants. 

# RUNNING
VERSO STEP #1 R demo script can be executed either by R GUI or from terminal, with the following command: 

	Rscript main.R

# OUTPUTS
VERSO STEP #1 returns as output and R list reporting the inferred maximum log-likelihood phylogenetic tree. 

# EXAMPLE FILES
As an example, the GitHub folder contains a demo RData named variants.RData which includes variants for a selected set of 15 SARS-CoV-2 samples obtained by variant calling from raw data available from NCBI BioProject PRJNA610428. 
