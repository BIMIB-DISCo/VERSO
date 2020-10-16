# install VERSO dependencies if needed
if (!require("ape")) install.packages("ape")
if (!require("Rfast")) install.packages("Rfast")

# source VERSO R scripts
source("R/frontend.R")
source("R/inference.R")
source("R/utils.R")

# load example (toy) data
load("data/variants.RData")

# set the grid search (alpha and beta)
alpha = c(0.01,0.05)
beta = c(0.01,0.05)

# perform the inference
inference = VERSO(D = variants,
                  alpha = alpha,
                  beta = beta,
                  check_indistinguishable = TRUE,
                  num_rs = 5,
                  num_iter = 100,
                  n_try_bs = 50,
                  num_processes = 1,
                  seed = 12345,
                  verbose = TRUE)

# plot the inferred phylogenetic tree
plot(inference$phylogenetic_tree)
