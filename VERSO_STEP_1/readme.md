# SUMMARY
VERSO STEP #1 processes binarized variant profiles (0 variant is absent, 1 present) for clonal variants obtained via variant calling or consensus sequences, in order to infer the maximum log-likelihood phylogenetic tree, where variants are inner nodes and samples are tips.

VERSO STEP #1 is provided as a set of R scripts (directory R) and a demo che be executed by running the script main.R. 

# REQUIREMENTS
Please install the following R libraries.

* [ape, see https://cran.r-project.org/web/packages/ape/index.html] with the command:
<pre><code>if (!require("ape")) install.packages("ape")</code></pre>

* [Rfast, see https://cran.r-project.org/web/packages/Rfast/index.html] with the command:
<pre><code>if (!require("Rfast")) install.packages("Rfast")</code></pre>

# INPUTS
  
VERSO STEP #1 requires as input an R matrix reporting binarized variant profiles for a set of clonal variants. 

This can be generated either by a preprocessed file (see data/variants.txt) or from aligned consensus sequences (see as an example the files in the directory example_alignment_consensus_sequences and specifically the file gisaid_example_aligned.fasta. 

In the first case (preprocessed file), data can be directly imported in R by reading the file. In the second case, raw fasta files first needs to be aligned; this can be done with many tools and here we provide as an example a bash script to do this task using augur from the Nextstrain pipeline (https://github.com/nextstrain/augur). To replicate our analysis, run the bash script run_alignment.sh; configuration of Nextstrain/augur is left to the user. 

We finally provide the R code to import VERSO STEP #1 input data both from preprocessed file and from consensus sequences in the script main_part1_data_import.R. It can be executed either by R GUI or from terminal, with the following command: 

	Rscript main_part1_data_import.R

# RUNNING
VERSO STEP #1 R demo script can be executed either by R GUI or from terminal, with the following command: 

	Rscript main_part2_inference.R

# OUTPUTS
VERSO STEP #1 returns as output an R list reporting the inferred maximum log-likelihood phylogenetic tree. 

# EXAMPLE FILES
As an example, the GitHub folder contains a demo RData named variants.RData which includes variants for a selected set of 15 SARS-CoV-2 samples obtained by variant calling from raw data available from NCBI BioProject PRJNA610428. 
