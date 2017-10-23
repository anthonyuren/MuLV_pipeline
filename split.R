# Read in the matrix file
setwd("/csc/analysis/Cscbioinf/mdore/Uren_Pooled4_runs/Merged/MergedBams/SortedPaired/Mate2/SplitBySample/")
m <- read.table("matrix.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)


# Calculate a total for each row and append to the matrix
row_total <- rowSums(m[,5:ncol(m)])
m <- cbind(m[,1:4], row_total, m[,5:ncol(m)])

## Extract only the count values into a matrix and then into an array
m_count <- m[,5:ncol(m)]
m_meta <- m[,1:4]
rm(m)

## Break the matrix down by sample
for(i in 1:ncol(m_count)) {
	rowNums <- m_count[,i]>0
	m_meta_new <- m_meta[rowNums,]
	m_count_new <- m_count[rowNums,i]
	sampleName <- colnames(m_count)[i]
	sampleName <- gsub("X", "", sampleName)
	#sampleName <- gsub("_", ".", sampleName)
	#sampleParts <- unlist(strsplit(sampleName, "\\."))
	#sampleName <- sampleParts[9]
	m_new <- cbind(m_meta_new, m_count_new)
	cat(paste(i, "\n"))
	write.table(m_new, file=paste("m_split_", sampleName, ".txt", sep=""), sep="\t", row.names=FALSE)
	write.table(m_meta_new, file=paste("m_meta_split_", sampleName, ".txt", sep=""), sep="\t", row.names=FALSE)
	write.table(m_count_new, file=paste("m_count_split_", sampleName, ".txt", sep=""), sep="\t", row.names=FALSE)
}

