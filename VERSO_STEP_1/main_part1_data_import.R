# reading binarized clonal profiles from preprocessed file
variants_preprocessed_file = read.table("data/variants.txt",header=TRUE,sep=",",row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
positions_of_interest = NULL
for(i in strsplit(colnames(variants_preprocessed_file),"_")) {
    positions_of_interest = c(positions_of_interest,i[[1]])    
}
positions_of_interest = as.numeric(positions_of_interest)

# reading example raw variants obtained from aligned consensus sequences
variants_raw_aligned_consensus_sequences = read.table(file="example_alignment_consensus_sequences/gisaid_example_aligned.fasta",header=FALSE,sep="\n",check.names=FALSE,stringsAsFactors=FALSE)

# get samples headers
start = NULL
for(i in 1:nrow(variants_raw_aligned_consensus_sequences)) {
    if(length(grep(">",variants_raw_aligned_consensus_sequences[i,]))>=1) {
        start = c(start,i)
    }
}
samples = list()
for(i in 1:length(start)) {
    if(i!=length(start)) {
        curr_s = start[i]
        curr_e = start[(i+1)]-1
    }
    else {
        curr_s = start[i]
        curr_e = nrow(variants_raw_aligned_consensus_sequences)
    }
    samples[[i]] = c(curr_s,curr_e)
}

# make data matrix
i = 1
reference = paste0(variants_raw_aligned_consensus_sequences[(samples[[i]][1]+1):samples[[i]][2],],collapse="")
data_matrix_seq = array(NA,c((length(samples)-1),nchar(reference)))
colnames(data_matrix_seq) = strsplit(paste0(variants_raw_aligned_consensus_sequences[(samples[[1]][1]+1):samples[[1]][2],],collapse=""),"")[[1]]
rownames_data = NULL
for(i in 2:length(samples)) {
    curr_data = paste0(variants_raw_aligned_consensus_sequences[(samples[[i]][1]+1):samples[[i]][2],],collapse="")
    data_matrix_seq[(i-1),] = strsplit(curr_data,"")[[1]]
    rownames_data = c(rownames_data,gsub(">","",variants_raw_aligned_consensus_sequences[(samples[[i]][1]),]))
    cat(i/length(samples),"\n")
}
rownames(data_matrix_seq) = rownames_data
variants_raw_aligned_consensus_sequences = data_matrix_seq

# make final binarized matrix
variants_binarized_aligned_consensus_sequences = array(0,dim(variants_raw_aligned_consensus_sequences[,positions_of_interest]))
rownames(variants_binarized_aligned_consensus_sequences) = rownames(variants_raw_aligned_consensus_sequences[,positions_of_interest])
colnames(variants_binarized_aligned_consensus_sequences) = colnames(variants_preprocessed_file)
for(i in 1:ncol(variants_binarized_aligned_consensus_sequences)) {
    curr_mutated = names(which(variants_raw_aligned_consensus_sequences[,positions_of_interest[i],drop=TRUE]!=colnames(variants_raw_aligned_consensus_sequences[,positions_of_interest[i],drop=FALSE])))
    variants_binarized_aligned_consensus_sequences[curr_mutated,i] = 1    
}
variants_binarized_aligned_consensus_sequences = variants_binarized_aligned_consensus_sequences[sort(rownames(variants_binarized_aligned_consensus_sequences)),]

# print the results for comparison
print("variants_preprocessed_file")
print(variants_preprocessed_file)
print("variants_binarized_aligned_consensus_sequences")
print(variants_binarized_aligned_consensus_sequences)
