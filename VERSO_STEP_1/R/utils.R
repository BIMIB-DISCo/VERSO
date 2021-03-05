# remove any indistinguishable variant from input data
"check.indistinguishable" <- function( data ) {
    
    # check for indistinguishable events
    indistinguishable <- as.numeric(which(duplicated(t(data))))
    if(length(indistinguishable)>0) {
        # merge names of indistinguishable events
        valid_colnames <- colnames(data)[-indistinguishable]
        new_colnames <- valid_colnames
        invalid_colnames <- colnames(data)[indistinguishable]
        for(i in invalid_colnames) {
            for(j in valid_colnames) {
                if(all(is.na(data[,i])==is.na(data[,j]))) { # if NAs in i and j are the same
                    if(all(data[,i]==data[,j],na.rm=TRUE)) { # if not-NAs entries in i and j are the same
                        new_colnames[which(valid_colnames==j)] <- paste0(new_colnames[which(valid_colnames==j)],"|",invalid_colnames[which(invalid_colnames==i)])
                        next;
                    }
                }
            }
        }
        # remove indistinguishable events from data
        data <- data[,valid_colnames,drop=FALSE]
        colnames(data) <- new_colnames
    }
    
    # return data
    return(data)
    
}

# build a phylogenetic tree from a variants tree (adjacency_matrix) and related samples attachments (samples_attachments)
"get.phylo" <- function( adjacency_matrix, valid_genotypes, samples_attachments ) {
    
    # compute Manhattan distance among valid genotypes
    distance_genotypes <- as.matrix(dist(valid_genotypes,method="manhattan"))

    # consider the variants tree (adjacency_matrix) to build inner nodes structure of VERSO phylogenetic tree
    edges <- which(adjacency_matrix==1,arr.ind=TRUE)
    parents_list <- rownames(adjacency_matrix)[as.numeric(edges[,"row"])]
    children_list <- rownames(adjacency_matrix)[as.numeric(edges[,"col"])]
    edges <- cbind(parents_list,children_list)
    edges_weights <- NULL
    for(i in 1:nrow(edges)) {
        curr_p <- rownames(valid_genotypes)[which(colnames(valid_genotypes)==edges[i,1])]
        curr_c <- rownames(valid_genotypes)[which(colnames(valid_genotypes)==edges[i,2])]
        edges_weights <- c(edges_weights,distance_genotypes[curr_p,curr_c])
    }
    nodes_list <- unique(as.vector(edges))
    Nnode <- length(nodes_list)
    parents_list_numeric <- parents_list
    children_list_numeric <- children_list
    for(nlist in 1:Nnode) {
        curr_nl <- which(parents_list==nodes_list[nlist])
        if(length(curr_nl)>0) {
            parents_list_numeric[curr_nl] <- nlist
        }
        curr_nl <- which(children_list==nodes_list[nlist])
        if(length(curr_nl)>0) {
            children_list_numeric[curr_nl] <- nlist
        }
    }
    edges <- cbind(as.numeric(parents_list_numeric),as.numeric(children_list_numeric))
    colnames(edges) <- c("Parent","Child")
    
    # now consider samples attachments
    attachments <- samples_attachments[,"Genotype"]
    tip.label <- c(names(attachments),"Reference_Genotype")
    overhead <- length(tip.label)
    edges <- edges + overhead
    for(i in 1:length(attachments)) {
        curr_att <- which(nodes_list==colnames(valid_genotypes)[which(rownames(valid_genotypes)==attachments[i])])+overhead
        edges <- rbind(edges,c(curr_att,i))
        edges_weights <- c(edges_weights,0)
    }
    edges <- rbind(edges,c(as.numeric(edges[1,"Parent"]),(i+1)))
    edges_weights <- c(edges_weights,0)
    rownames(edges) <- 1:nrow(edges)

    # build VERSO phylogenetic tree
    phylogenetic_tree <- list(edge=edges,tip.label=tip.label,Nnode=Nnode,edge.length=edges_weights,node.label=nodes_list,root.edge=0)
    class(phylogenetic_tree) <- "phylo"
    phylogenetic_tree <- keep.tip(phylogenetic_tree,tip.label)

    # return VERSO phylogenetic tree
    return(phylogenetic_tree)

}

# convert B to an adjacency matrix
"as.adj.matrix" <- function( B ) {

    # create the data structure where to save the adjacency matrix obtained from B
    adj_matrix <- array(0L,dim(B))
    rownames(adj_matrix) <- colnames(B)
    colnames(adj_matrix) <- colnames(B)

    # set arcs in the adjacency matrix
    for(i in 1:(nrow(B)-1)) {
        for(j in ((i+1):nrow(B))) {
            if(all(B[i,1:i]==B[j,1:i])&&(sum(B[i,])==(sum(B[j,])-1))) {
                adj_matrix[i,j] <- 1
            }
        }
    }

    # return the adjacency matrix obtained from B
    return(adj_matrix)

}

# build B from an adjacency matrix where we assume genotypes and mutations to be both ordered
"as.B" <- function( adj_matrix, D ) {
    
    # build data structure to save results
    n_clones <- nrow(adj_matrix)
    B <- diag(n_clones)
    
    # build B
    for(k in 1:n_clones) {
        idx_child <- which(adj_matrix[k,]==1,arr.ind=TRUE)
        if(length(idx_child)==1) {
            B[idx_child,] <- B[k,] + B[idx_child,]
        }
        else if(length(idx_child)>1) {
            B[idx_child,] <- sweep(B[idx_child,],2,B[k,],"+")
        }
    }
    rownames(B) <- c("r",1:(nrow(B)-1))
    mycolnames <- "r"
    for(i in 2:nrow(adj_matrix)) {
        mycolnames <- c(mycolnames,as.character(which(rownames(adj_matrix)[i]==colnames(D))))
    }
    colnames(B) <- mycolnames
    
    # return B
    return(B)
    
}

# perform node relabeling
"relabeling" <- function(B) {
    
    # relabeling
    chosen <- sample(2:ncol(B),2,replace=FALSE)
    tmp <- colnames(B)[chosen[1]]
    colnames(B)[chosen[1]] <- colnames(B)[chosen[2]]
    colnames(B)[chosen[2]] <- tmp
    return(B)
    
}

# perform prune and reattach
"prune.and.reattach" <- function(B) {
    
    # change one arch
    is_not_valid <- TRUE
    while(is_not_valid) {
        
        #Select source node
        ch_1 <- sample(x=2:nrow(B),size=1)
        ch_1_gen <- B[ch_1,1:ch_1]
        remaing_node <- as.numeric(which(apply(B[,1:ch_1], c(1), FUN = function(x){!all(x == ch_1_gen)})))
        
        # chose the target node from the nodes not included in the subtree where ch_1 is the root
        if(length(remaing_node) > 1) {
            
            ch_2 <- sample(x=remaing_node,size=1)
            
        } else if(length(remaing_node)==1) {
            
            ch_2 <- remaing_node
            
        } else {
            # if there aren't any nodes, select a different source
            next;
        }

        # a pair of two nodes is valid if the nodes are not already directly connected
        if(!(all(B[ch_1,1:ch_2]==B[ch_2,1:ch_2])&sum(B[ch_1,])==(sum(B[ch_2,])+1))) {
            is_not_valid <- FALSE
        }
        
    }
    
    descendent_nodes <- setdiff(1:nrow(B), remaing_node)
    
    # extract descendent node submatrix
    rem_B <- B[remaing_node,remaing_node,drop=FALSE]
    new_B <- matrix(data = 0L, nrow = nrow(rem_B), ncol = ncol(B))
    colnames(new_B) <- c(colnames(B)[remaing_node], colnames(B)[descendent_nodes])
    new_B[1:length(remaing_node),1:length(remaing_node)] <- rem_B
    gen_ch_2 <- new_B[which(colnames(new_B)==colnames(B)[ch_2]),1:length(remaing_node)]
    desc_B <- cbind(matrix(rep(gen_ch_2,each=length(descendent_nodes)),  nrow = length(descendent_nodes)), 
                    B[descendent_nodes,descendent_nodes,drop=FALSE])
    new_B <- rbind(new_B,desc_B)
    
    return(new_B)
    
}
