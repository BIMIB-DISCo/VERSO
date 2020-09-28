#' Perform the inference of the maximum log-likelihood VERSO phylogenetic tree.
#' @title VERSO
#'
#' @examples
#' data(variants)
#' inference = VERSO(D = variants, 
#'                   alpha = c(0.01,0.05), 
#'                   beta = c(0.01,0.05), 
#'                   check_indistinguishable = TRUE, 
#'                   num_rs = 10, 
#'                   num_iter = 100, 
#'                   n_try_bs = 50, 
#'                   num_processes = 1, 
#'                   seed = 12345, 
#'                   verbose = FALSE)
#'
#' @param D Input data for the inference reporting presence (as 1), absense (as 0) or missing information (as NA) for a set of variants.
#' @param alpha False positive error rate provided as a verctor; if a vector of alpha (and beta) is provided, the inference is performed 
#' for multiple values and the solution at maximum-likelihood is returned.
#' @param beta False negative error rate provided as a verctor; if a vector of beta (and alpha) is provided, the inference is performed 
#' for multiple values and the solution at maximum-likelihood is returned.
#' @param check_indistinguishable Boolean. Shall I remove any indistinguishable variant from input data prior inference?
#' @param num_rs Number of restarts during MCMC inference.
#' @param num_iter Maximum number of MCMC steps to be performed during the inference.
#' @param n_try_bs Number of steps without changes in likelihood of best solution after which to stop the MCMC.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either 1, NA or NULL.
#' @param seed Seed for reproducibility.
#' @param verbose Boolean. Shall I print to screen information messages during the execution?
#' @return A list of 8 elements: B, C, phylogenetic_tree, corrected_genotypes, genotypes_prevalence, genotypes_summary, log_likelihood and error_rates. 
#' Here, B returns the maximum likelihood variants tree (inner nodes of the phylogenetic tree), C the attachment of patients to genotypes and phylogenetic_tree 
#' VERSO phylogenetic tree, including both variants tree and patients attachments to variants; corrected_genotypes is the corrected genotypes, which corrects 
#' D given VERSO phylogenetic tree, genotypes_prevalence the number of patients and observed prevalence of each genotype and genotypes_summary provide a 
#' summary of association of mutations to genotypes; finally log_likelihood and error_rates return the likelihood of the inferred phylogenetic moldel and 
#' best values of alpha and beta as estimated by VERSO.
#' @export VERSO
#' @import ape
#' @import parallel
#' @importFrom phangorn keep.tip
#' @importFrom Rfast rowMaxs
#' @importFrom stats runif dist
#'
"VERSO" <- function( D, alpha = NULL, beta = NULL, check_indistinguishable = TRUE, num_rs = 10, num_iter = 10000, n_try_bs = 1000, num_processes = Inf, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    # set storage mode to integer
    storage.mode(D) <- "integer"
    
    # remove any indistinguishable variant from input data prior inference
    if(check_indistinguishable) {
        D <- check.indistinguishable(D)
        storage.mode(D) <- "integer"
    }

    # initialize error rates alpha and beta if not provided as inputs
    if(is.null(alpha)) {
        alpha = c(0.01,0.05,0.10)
        beta = c(0.01,0.05,0.10)
    }
    
    if(verbose) {
        cat(paste0("Performing inference for a total of ",length(alpha)," different values of alpha and beta.","\n"))
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.null(parallel)&&length(alpha)>1) {
        if(is.na(num_processes)||is.null(num_processes)||num_processes==1) {
            parallel <- NULL
        }
        else if(num_processes==Inf) {
            cores <- as.integer((detectCores()-1))
            if(cores < 2) {
                parallel <- NULL
            }
            else {
                num_processes <- min(num_processes,length(alpha))
                parallel <- makeCluster(num_processes,outfile="")
                close_parallel <- TRUE
            }
        }
        else {
            num_processes <- min(num_processes,length(alpha))
            parallel <- makeCluster(num_processes,outfile="")
            close_parallel <- TRUE
        }
        if(verbose && !is.null(parallel)) {
            cat("Executing",num_processes,"processes via parallel...","\n")
        }
    }

    # now start the inference
    if(is.null(parallel)) {

        # sequential computation
        inference <- list()
        for(i in 1:length(alpha)) {
            
            if(verbose) {
                cat('Performing inference for alpha =',paste0(alpha[i],collapse=" | "),'and beta =',paste0(beta[i],collapse=" | "),'\n')
            }

            inference[[i]] <- learn.VERSO.phylogenetic.tree(D = D, 
                                                            alpha = alpha[i], 
                                                            beta = beta[i], 
                                                            initialization = initialization, 
                                                            num_rs = num_rs, 
                                                            num_iter = num_iter, 
                                                            n_try_bs = n_try_bs, 
                                                            seed = round(runif(1)*10000), 
                                                            verbose = verbose)
            
        }

    }
    else {

        # parallel computation
        res_clusterEvalQ <- clusterEvalQ(parallel,library("Rfast"))
        clusterExport(parallel,varlist=c("D","alpha","beta","initialization","num_rs","num_iter","n_try_bs","verbose"),envir=environment())
        clusterExport(parallel,c("learn.VERSO.phylogenetic.tree","initialize.B","move.B","compute.C"),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        inference <- parLapply(parallel,1:length(alpha),function(x) {
            
            if(verbose) {
                cat('Performing inference for alpha =',paste0(alpha[x],collapse=" | "),'and beta =',paste0(beta[x],collapse=" | "),'\n')
            }
            
            inference <- learn.VERSO.phylogenetic.tree(D = D, 
                                                       alpha = alpha[x], 
                                                       beta = beta[x], 
                                                       initialization = initialization, 
                                                       num_rs = num_rs, 
                                                       num_iter = num_iter, 
                                                       n_try_bs = n_try_bs, 
                                                       seed = round(runif(1)*10000), 
                                                       verbose = verbose)

        })
        
    }
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }

    # return the solution at maximum log-likelihood among the inferrend ones
    lik <- NULL
    for(i in 1:length(inference)) {
        lik <- c(lik,inference[[i]][["log_likelihood"]])
    }
    best <- which(lik==max(lik))[1]
    inference <- inference[[best]]
    B <- inference$B
    storage.mode(B) <- "integer"
    C <- inference$C
    log_likelihood <- inference[["log_likelihood"]]
    error_rates <- list(alpha=alpha[best],beta=beta[best])

    # compute corrected genotypes
    C_matrix <- array(0,c(nrow(D),ncol(B)))
    rownames(C_matrix) <- rownames(D)
    colnames(C_matrix) <- colnames(B)
    for(i in 1:nrow(C_matrix)) {
        C_matrix[i,which(rownames(B)==C[rownames(C_matrix)[i],"Genotype"])] <- 1
    }
    corrected_genotypes <- C_matrix %*% B
    corrected_genotypes <- corrected_genotypes[,colnames(D)]
    storage.mode(corrected_genotypes) <- "integer"

    # compute genotypes prevalence
    genotypes_prevalence <- array(NA,c(nrow(B),2))
    rownames(genotypes_prevalence) <- rownames(B)
    colnames(genotypes_prevalence) <- c("#Patients","Prevalence")
    patients <- table(C[,"Genotype"])
    prevalence <- patients/nrow(C)
    genotypes_prevalence[names(patients),"#Patients"] <- as.numeric(patients)
    genotypes_prevalence[names(prevalence),"Prevalence"] <- as.numeric(prevalence)

    # compute genotypes summary
    genotypes_summary <- list()
    i = 1
    for(genotypes in rownames(B)) {
        mut_list <- colnames(B)[B[i,]==1]
        genotypes_summary[[genotypes]] <- mut_list
        i = i + 1
    }

    # finally build VERSO phylogenetic tree
    phylogenetic_tree <- get.phylo(adjacency_matrix=as.adj.matrix(B),valid_genotypes=B,samples_attachments=C)

    # return results of the inference
    results <- list(B=B,C=C,phylogenetic_tree=phylogenetic_tree,corrected_genotypes=corrected_genotypes,genotypes_prevalence=genotypes_prevalence,genotypes_summary=genotypes_summary,log_likelihood=log_likelihood,error_rates=error_rates)
    return(results)

}
