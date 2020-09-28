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
    }
    
    # return data
    return(data)
    
}
