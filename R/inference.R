# learn VERSO phylogenetic tree from data
"learn.VERSO.phylogenetic.tree" <- function( D, alpha = 10^-3, beta = 10^-3, initialization = NULL, num_rs = 10, num_iter = 10000, n_try_bs = 1000, seed = NULL, verbose = TRUE ) {
    
    # set the seed
    set.seed(seed)

    # initialize global variables (best solution among all restarts)
    B_global <- NULL
    C_global <- NULL
    lik_global <- NULL
    
    # we remove any NA from data
    if(any(is.na(D))) {
        D[which(is.na(D),arr.ind=TRUE)] <- -3
    }
    storage.mode(D) <- "integer"
    
    # perform a number of num_rs restarts
    for(i in 1:num_rs) {
        
        if(verbose) {
            message("Performing restart number ",i," out of ",num_rs)
        }
        
        # initialize B
        if(i==1) {
            if(is.null(initialization)) {
                B <- initialize.B(D)
            }
            else {
                B <- initialization
                storage.mode(B) <- "integer"
            }
        }
        else {
            # perform random shuffling of nodes/variants ordering
            B <- B_global
            colnames(B) <- c("r",sample(1:(ncol(B)-1)))
        }
        
        # compute C given B
        res <- compute.C(B,D,alpha,beta)
        C <- res$C
        lik <- res$lik
        
        # initialize result variables (best solution for the current restart)
        B_best <- B
        C_best <- C
        lik_best <- lik
        count_lik_best_cons <- 0
        
        # repeat MCMC moves until num_iter number of iterations is performed
        for(j in 1:num_iter) {
            
            if(verbose&&(j%%100)==0) {
                message("Performing iteration number ",j," out of ",num_iter," | Current best log-likelihood ",lik_best)
            }
            
            # perform a move on B
            B_tmp <- move.B(B)
            
            # compute C at maximun likelihood given B_tmp and returns its likelihood
            res <- compute.C(B_tmp,D,alpha,beta)
            C_tmp <- res$C
            lik_tmp <- res$lik
            
            # if likelihood at current step is better than best likelihood, replace best model with current
            if(lik_tmp>lik_best) {
                B_best <- B_tmp
                C_best <- C_tmp
                lik_best <- lik_tmp
                count_lik_best_cons <- 0
                B <- B_best
                C <- C_best
                lik <- lik_best
            }
            else {
                # update count
                count_lik_best_cons <- count_lik_best_cons + 1
                if(count_lik_best_cons>n_try_bs) {
                    # print a message
                    if(verbose) {
                        message("Not improving likelihood of best solution after ",n_try_bs," iterations. Skipping to next restart")
                    }
                    break;
                }
                # take the current state with a probability proportional to the ratio of the two likelihoods
                rho <- min(exp((lik_tmp-lik)),1)
                if(runif(1)<=rho) {
                    B <- B_tmp
                    C <- C_tmp
                    lik <- lik_tmp
                }
            }
            
        }

        # compare current best solution with best global solution
        if(is.null(lik_global)||(lik_best>lik_global)) {
            B_global <- B_best
            C_global <- C_best
            lik_global <- lik_best
        }
        
    }

    # renaming
    rownames(B_global) <- paste0("G",1:nrow(B_global))
    colnames(B_global) <- c("Reference",colnames(D)[as.numeric(colnames(B_global)[2:ncol(B_global)])])
    C_global <- matrix(paste0("G",C_global),ncol=1)
    rownames(C_global) <- rownames(D)
    colnames(C_global) <- "Genotype"

    
    # return maximum likelihood solution
    return(list(B=B_global,C=C_global,log_likelihood=lik_global))
    
}

# randomly initialize B
"initialize.B" <- function( D ) {
    
    # data structure where to save B
    B <- array(0L,c((ncol(D)+1),(ncol(D)+1)))
    rownames(B) <- c("r",1:ncol(D))
    colnames(B) <- c("r",sample(1:ncol(D)))
    diag(B) <- 1L
    B[,1] <- 1L
    
    # add arcs with probability 50% to obtain a random tree topology
    p <- 0.50
    for(i in 2:(nrow(B)-1)) {
        if(runif(1)<p) {
            B[(i+1),] <- B[i,] + B[(i+1),]
        }
    }
    B[which(B>1)] <- 1L
    
    # return B
    return(B)
    
}

# performing either relabeling or edge changing moves on B
"move.B" <- function( B ) {
    
    # sample a random probability of choosing a move
    p <- runif(1)

    # perform pairwise relabeling with 55% probability
    if(p<0.55) {

        # nodes relabeling
        chosen <- sample(2:ncol(B),2,replace=FALSE)
        tmp <- colnames(B)[chosen[1]]
        colnames(B)[chosen[1]] <- colnames(B)[chosen[2]]
        colnames(B)[chosen[2]] <- tmp

    }
    # perform structural moves with 40% probability
    else if(p>=0.55&&p<0.95) {

        # change one arch
        is_not_valid <- TRUE
        while(is_not_valid) {
            ch_1 <- sample(3:nrow(B),1)
            ch_2 <- sample(1:(ch_1-1),1)
            # a pair of two nodes are a valid set if the nodes are not already directly connected
            if(!(all(B[ch_1,1:ch_2]==B[ch_2,1:ch_2])&&sum(B[ch_1,])==(sum(B[ch_2,])+1))) {
                is_not_valid <- FALSE
            }
        }

        # performing move on ch_1
        ch_1_bkp <- B[ch_1,1:ch_1]
        B[ch_1,1:(ch_1-1)] <- c(1L,rep(0L,(ch_1-2)))
        B[ch_1,] <- B[ch_1,] + B[ch_2,]
        
        # performing move on children of ch_1
        if(ch_1 != nrow(B)) {
            for(i in (ch_1+1):nrow(B)) {
                if(all(ch_1_bkp==B[i,1:ch_1])) {
                    B[i,1:(ch_1-1)] <- c(1L,rep(0L,(ch_1-2)))
                    B[i,] <- B[i,] + B[ch_2,]
                }
            }
        }
        B[which(B>1)] <- 1L
        
    }
    # perform full relabeling with 5% probability
    else if(p>=0.95) {

        # random relabeling of all clones
        colnames(B) <- c("r",sample(1:(ncol(B)-1)))

    }
    
    # return B
    return(B)

}

# compute attachments matrix C at maximum likelihood given B and D
"compute.C" <- function( B, D, alpha = 10^-3, beta = 10^-3 ) {
    
    # determine indeces to order D such that it matches B
    idx_srt <- as.integer(colnames(B)[2:ncol(B)])
    
    # go through each patient and compute likelihood for all possible attachments
    lik_matrix <- array(0L,c(nrow(D),ncol(B)))
    curr_D <- cbind(rep(1,nrow(D[,idx_srt,drop=FALSE])),D[,idx_srt,drop=FALSE])
    for(k in 1:nrow(B)) {
        curr_C = matrix(rep(0L,nrow(B)),nrow=1)
        curr_C[1,k] <- 1L
        r_D_tilde <- (curr_C%*%B)*2
        sum_cell <- as.matrix(sweep(curr_D,MARGIN=2,r_D_tilde,"+"))
        lik_matrix[,k] <- (beta^Rfast::rowSums(sum_cell==2)) * ((1-beta)^Rfast::rowSums(sum_cell==0)) * ((alpha)^Rfast::rowSums(sum_cell==1)) * ((1-alpha)^Rfast::rowSums(sum_cell==3))
    }

    # compute maximum likelihood attachments
    C <- rowMaxs(lik_matrix,value=FALSE)
    storage.mode(C) <- "integer"
    lik <- sum(log(rowMaxs(lik_matrix,value=TRUE)))
    
    # return maximum likelihood attachments
    return(list(C=C,lik=lik))
    
}
