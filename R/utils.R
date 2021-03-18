# remove any indistinguishable variant from input data
check.indistinguishable <- function( data ) {
    
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
get.phylo <- function( adjacency_matrix, valid_genotypes, samples_attachments ) {
    
    # compute Manhattan distance among valid genotypes
    distance_genotypes <- as.matrix(dist(valid_genotypes,method="manhattan"))

    # consider the variants tree (adjacency_matrix) to build inner nodes structure of VERSO phylogenetic tree
    edges <- which(adjacency_matrix==1,arr.ind=TRUE)
    parents_list <- rownames(adjacency_matrix)[as.numeric(edges[,"row"])]
    children_list <- rownames(adjacency_matrix)[as.numeric(edges[,"col"])]
    edges <- cbind(parents_list,children_list)
    edges_weights <- NULL
    for(i in seq_len(nrow(edges))) {
        curr_p <- rownames(valid_genotypes)[which(colnames(valid_genotypes)==edges[i,1])]
        curr_c <- rownames(valid_genotypes)[which(colnames(valid_genotypes)==edges[i,2])]
        edges_weights <- c(edges_weights,distance_genotypes[curr_p,curr_c])
    }
    nodes_list <- unique(as.vector(edges))
    Nnode <- length(nodes_list)
    parents_list_numeric <- parents_list
    children_list_numeric <- children_list
    for(nlist in seq_len(Nnode)) {
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
    for(i in seq_len(length(attachments))) {
        curr_att <- which(nodes_list==colnames(valid_genotypes)[which(rownames(valid_genotypes)==attachments[i])])+overhead
        edges <- rbind(edges,c(curr_att,i))
        edges_weights <- c(edges_weights,0)
    }
    edges <- rbind(edges,c(as.numeric(edges[1,"Parent"]),(i+1)))
    edges_weights <- c(edges_weights,0)
    rownames(edges) <- seq_len(nrow(edges))

    # build VERSO phylogenetic tree
    phylogenetic_tree <- list(edge=edges,tip.label=tip.label,Nnode=Nnode,edge.length=edges_weights,node.label=nodes_list,root.edge=0)
    class(phylogenetic_tree) <- "phylo"
    phylogenetic_tree <- ape::keep.tip(phylogenetic_tree,tip.label)

    # return VERSO phylogenetic tree
    return(phylogenetic_tree)

}

# convert B to an adjacency matrix
as.adj.matrix <- function( B ) {

    # create the data structure where to save the adjacency matrix obtained from B
    adj_matrix <- array(0L,dim(B))
    rownames(adj_matrix) <- colnames(B)
    colnames(adj_matrix) <- colnames(B)

    # set arcs in the adjacency matrix
    for(i in seq_len((nrow(B)-1))) {
        for(j in ((i+1):nrow(B))) {
            if(all(B[i,seq_len(i)]==B[j,seq_len(i)])&&(sum(B[i,])==(sum(B[j,])-1))) {
                adj_matrix[i,j] <- 1
            }
        }
    }

    # return the adjacency matrix obtained from B
    return(adj_matrix)

}

# build B from an adjacency matrix where we assume genotypes and mutations to be both ordered
as.B <- function( adj_matrix, D ) {
    
    # build data structure to save results
    n_clones <- nrow(adj_matrix)
    B <- diag(n_clones)
    
    # build B
    for(k in seq_len(n_clones)) {
        idx_child <- which(adj_matrix[k,]==1,arr.ind=TRUE)
        if(length(idx_child)==1) {
            B[idx_child,] <- B[k,] + B[idx_child,]
        }
        else if(length(idx_child)>1) {
            B[idx_child,] <- sweep(B[idx_child,],2,B[k,],"+")
        }
    }
    rownames(B) <- c("r",seq_len((nrow(B)-1)))
    mycolnames <- "r"
    for(i in 2:nrow(adj_matrix)) {
        mycolnames <- c(mycolnames,as.character(which(rownames(adj_matrix)[i]==colnames(D))))
    }
    colnames(B) <- mycolnames
    
    # return B
    return(B)
    
}

# draw B
draw.B <- function( B, mut_label = colnames(B)[-1], last_mut_node_label =  TRUE) {
    
    Broot <- Node$new('r')
    Broot$mut <- B[1,]
    
    nClone <- nrow(B)
    Clones <- list(Broot)
    
    for(rP in 1:(nrow(B)-1)) {
        
        for(rC in ((rP+1):nrow(B))) {
            
            if(all(Clones[[rP]]$mut[1:rP]==B[rC,1:rP])&&(sum(Clones[[rP]]$mut)==(sum(B[rC,])-1))) {
                if(last_mut_node_label) {
                    mutName <- tail(mut_label[which(B[rC,-1]==1)], 1)
                } else {
                    mutName <- paste(mut_label[which(B[rC,-1]==1)],collapse="")  
                }
                
                Clones[[rC]] <- Clones[[rP]]$AddChild(mutName)
                Clones[[rC]]$mut <- B[rC,]
                
            }
            
        }
        
    }
    
    return(Broot)
    
}
