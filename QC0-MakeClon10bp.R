# make clon10bp.txt file
# a large file listing inserts by their chromosome orientation, position, sample and mouse with many empty fields that will be filled in later
# inserts for each sample are assigned a unique id in the last column

#
setwd(".")
# load sample list
mouselist <- read.table("sampleTable.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=NULL)
colnames(mouselist) <- c("sample_id","mouse_id")

setwd("TOTs")
insertfiles <- dir(pattern = "TOT")
insertnums <-insertfiles
all_inserts <- matrix(nrow = 0, ncol = 5)

for (i in insertfiles) {
  insertnum <- i
  num<-unlist(strsplit(i, "TOT"))
  num<-unlist(strsplit(num,"_"))[1]
  mouse <- mouselist$mouse_id[which(mouselist$sample_id == num)]
  print(paste(num,mouse,sep=" "))
  table <- read.table(paste(insertnum), sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=NULL)
  table <- cbind(table,num)
  table <- cbind(table,mouse)
  all_inserts <- rbind(all_inserts,table)
}

clon10bp <- as.data.frame(matrix(data = "0", nrow =nrow(all_inserts), ncol = 38))

clon10bp[,1] <- gsub("chr","",all_inserts[,1]) #chr
clon10bp[,1] <- gsub("X","20",clon10bp[,1])
clon10bp[,1] <- gsub("Y","21",clon10bp[,1])
clon10bp[,2] <- all_inserts[,6] #startbase
clon10bp[,3] <- all_inserts[,7] #endbase
clon10bp[,5] <- all_inserts[,3] #ori
clon10bp[,12] <- all_inserts[,9] #sample
clon10bp[,13] <- all_inserts[,10] #mouse
clon10bp[,18] <- all_inserts[,4] #fragments
clon10bp[,19] <- all_inserts[,5] #reads
clon10bp[,20] <- all_inserts[,2] #bestbase
clon10bp[,38] <- 1:nrow(all_inserts) #give each insert a unique ID

write.table(clon10bp, file = "clon10bp.txt" , sep="\t", col.names = F, row.names = F, append = FALSE, quote=FALSE) 

