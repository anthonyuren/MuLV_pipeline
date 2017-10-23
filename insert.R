#Scripts for clustering of insertions

args <- commandArgs(trailingOnly = FALSE)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
cat(myargument)
cat(paste(myargument, ".txt", sep=""))

# Read in the matrix file
sample <- read.table( paste(myargument, ".txt", sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=NULL)

colnames(sample)<-c("chr","LTRpos","ori","fragments","reads","startPos","endPos","Frequency")

finalinsert<-matrix(nrow=0, ncol=11)
                      
colnames(finalinsert)<-c("run_id","sample_id","chr","ori","fragments","reads","bestbase","medianbase","allbase","startpos","endpos")

for (i in 1:nrow(sample)){
  base<-""
  count<-""
  medianList<-NA
  for (el in unlist(strsplit(sample[i,8], ", "))){
   base<-c(base, unlist(strsplit(el, ":"))[1])
   thisbase<-unlist(strsplit(el, ":"))[1]
   count<-c(count, unlist(strsplit(el, ":"))[2])
   thiscount<-unlist(strsplit(el, ":"))[2]
   medianList<-na.omit(c(medianList, rep(as.numeric(thisbase), as.numeric(thiscount))))
  }
    
  medianbase<-median(sort(medianList))
  bestbase<-base[which(count == max(na.omit(as.numeric(count))))]
  finalinsert<-rbind(finalinsert, c(1,gsub("[^\\d]+", "", myargument, perl=TRUE), #run_id, sample_id
                                    sample[i,1],sample[i,3],sample[i,4],  #chr, ori, fragments
                                    sample[i,5], toString(bestbase), toString(medianbase), toString(sample[i,8]), #reads, bestbase, medianbase, all bases
                                    sample[i,6],sample[i,7]) ) #start pos, end pos
   
  
}

finalinsert<-finalinsert[order(finalinsert[,3],finalinsert[,4],as.numeric(finalinsert[,8])),]

write.table(finalinsert,file="inserts.txt", sep="\t", col.names = F, row.names = F, append = TRUE, quote = FALSE)
