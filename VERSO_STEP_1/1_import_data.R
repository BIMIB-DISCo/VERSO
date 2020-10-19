# set working directory
setwd("/Users/daniele/Desktop/GISAID/alignment/post_processing_results")

# read data
data = read.table(file="results/alignment.fasta",header=FALSE,sep="\n",check.names=FALSE,stringsAsFactors=FALSE)

# get samples headers
start = NULL
for(i in 1:nrow(data)) {
    if(length(grep(">",data[i,]))>=1) {
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
        curr_e = nrow(data)
    }
    samples[[i]] = c(curr_s,curr_e)
}

# make data matrix
i = 1
reference = paste0(data[(samples[[i]][1]+1):samples[[i]][2],],collapse="")
data_matrix_seq = array(NA,c((length(samples)-1),nchar(reference)))
colnames(data_matrix_seq) = strsplit(paste0(data[(samples[[1]][1]+1):samples[[1]][2],],collapse=""),"")[[1]]
rownames_data = NULL
for(i in 2:length(samples)) {
    curr_data = paste0(data[(samples[[i]][1]+1):samples[[i]][2],],collapse="")
    data_matrix_seq[(i-1),] = strsplit(curr_data,"")[[1]]
    rownames_data = c(rownames_data,gsub(">","",data[(samples[[i]][1]),]))
    cat(i/length(samples),"\n")
}
rownames(data_matrix_seq) = rownames_data

# save results
save(data_matrix_seq,file="RData/data_matrix_seq.RData")
