setwd(".")
library(plyr)

#Script that take the mouse ID and sample ID for each replicated insertion

# run this section once to break the insert list up by chromosomes to speed up processing of next steps

#read input file, a large matrix but only a handful of fields will be used
m <- read.table("clon10bp.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE)

#dilution/mixed samples are removed from the entire dataset prior to filtering
m<-m[(m[,13]!=1208),]
m<-m[(m[,13]!=1189),]

# #order by orientation, chromosome and then startpos
x<-m[order(m[,5], m[,1] ,as.numeric(m[,2])),]

#for each chromosome create an insert file
for (i in 1:max(x[,1])) {
 write.table( x[(x[1]==i),] , file=paste(i,"_chr.txt", sep=""), sep="\t", col.names = F, row.names = F, quote=FALSE)
}

##############################################################################

# args <- commandArgs(trailingOnly = FALSE)

chr_files <- dir(pattern = "chr.txt") #list all the relevant files
for (file_num in 1:length(chr_files)) { # for each file i.e. for each chromosome forward file and reverse file
  
myargument <- chr_files[file_num]
myargument2 <- sub(".txt","",myargument)
cat(myargument, "\n")
cat(myargument2, "\n")
# cat(myargument)
# cat(paste(myargument, ".txt", sep=""))

# Read in the matrix file
x <- read.table(myargument, sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=NULL)

row.names(x)<-NULL

#matrix that contains only orientation, best base, sample id, insert id, mouse id, chromosome, pos, end pos
contam<-x[,c(5,20,12,38,13,1,18,2,3)]
colnames(contam)<-c("ori", "bestbase", "sample", "insert_id", "mouse","chr", "fragments","startpos","endpos")

#windowSize<-9  ###
duplicate<-matrix(nrow=0, ncol=10)
colnames(duplicate)<-c("ori", "bestbase", "sample", "insert_id", "mouse","chr", "fragments","startpos","endpos", "dup")

for (i in 1:nrow(contam)){
  
  # Determine the current position and criteria
  currentChr <- contam[i,6]
  currentOri <- contam[i,1]
  currentBestBase <- contam[i,2]
  
  gap <- 0
  if(i < nrow(contam) && i > 1 ) { # if i is not the first or last line of inserts_table
    as.numeric(unlist(strsplit(currentBestBase, ", ")))
    
    gap1 <- contam$endpos[i+1] - contam$startpos[i] #gap ahead
    gap2 <- contam$startpos[i] - contam$endpos[i-1]  #gap behind
    gap <- min(c(gap1,gap2))
  }
  
  if (gap <1000) { #if the next base or previous base is close (1000 is overkill, 1 should do it)
    
    if (grepl(",",currentBestBase)){
      
      density<-0
      
       for (el in unlist(strsplit(currentBestBase, ", "))) {
 
        # Count how many inserts appear within the window and meet the qualifying criteria
        density <- max(density, nrow(contam[(contam[6]==currentChr & grepl(paste("\\<",el,"\\>", sep=""), contam[,2]) & contam[1]==currentOri),]))
      }
      
      duplicateV<-cbind(contam[i,],density)
      duplicate<-rbind(duplicate, duplicateV)
      
    } else{
      
      # Count how many inserts appear within the window and meet the qualifying criteria
      density <- nrow(contam[(contam[6]==currentChr & grepl(paste("\\<",currentBestBase,"\\>", sep=""), contam[,2]) & contam[1]==currentOri),])
      duplicateV<-cbind(contam[i,],density)
      duplicate<-rbind(duplicate, duplicateV)
    }  
    } else { # if the next base and previous base are more than "gap" away
      density <- 1
      duplicateV<-cbind(contam[i,],density)
      duplicate<-rbind(duplicate, duplicateV)
    }
  
  if (i%%1000==0){
    cat(i)
  }
  
}
duplicate <- duplicate[,c(1,2,3,4,5,6,7,10)]

write.table(duplicate, file=paste(myargument2, "_dup.txt", sep=""), sep="\t", col.names = F, row.names = F, quote=FALSE) 

}
