setwd(".")

#read input file
m <- read.table("clon10bp.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE)

# remove dilution/mixed samples
m<-m[(m[,13]!=1208),]
m<-m[(m[,13]!=1189),]

#order by orientation and then start pos (because best base is a string)
x<-m[order(m[,5], as.numeric(m[,2])),]
# x<-x[-33241,]

row.names(x)<-NULL

#matrix that contains only orientation, best base, sample id, insert id, mouse id, chromosome
contam<-x[,c(5,20,12,38,13,1,18)]
colnames(contam)<-c("ori", "bestbase", "sample", "insert_id", "mouse","chr", "fragments")

#add a column with a unique string (chr-ori-base) so it's easier to check later   ########################
contam<-cbind(contam, paste("chr",contam$chr, " ori", contam$ori,  " ", contam$bestbase, sep=" "))
colnames(contam)<-c("ori", "bestbase", "sample", "insert_id", "mouse","chr", "fragments", "string")

# prior to running this script merge these individual chr files
#cat 1_chr_dup.txt 2_chr_dup.txt 3_chr_dup.txt 4_chr_dup.txt 5_chr_dup.txt 6_chr_dup.txt 7_chr_dup.txt 8_chr_dup.txt 9_chr_dup.txt 10_chr_dup.txt 11_chr_dup.txt 12_chr_dup.txt 13_chr_dup.txt 14_chr_dup.txt 15_chr_dup.txt 16_chr_dup.txt 17_chr_dup.txt 18_chr_dup.txt 19_chr_dup.txt 20_chr_dup.txt 21_chr_dup.txt > tot.txt

duplicate <- read.table("tot.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=NULL)

#take only the duplicates (i.e. scripts with counter > 1)
duplicates_only<-duplicate[(duplicate[,8]>1),]

MouseListCat<-matrix(nrow=0,ncol=4)

# deal with only one chromosome at a time to reduce search time for each base
for (chromosome in 1:21) {
  duplicateLess<-duplicates_only[(duplicates_only[,6]==chromosome),]
  contamination<-contam[(contam[,6]==chromosome),]
#for each string (chr-ori-base) take the list of mice that have that insert
MouseList<-matrix(nrow=nrow(duplicateLess),ncol=4)

for (i in 1:nrow(duplicateLess)) {
  MouseList[i,1]<-toString(paste("chr",duplicateLess[i,6], " ori", duplicateLess[i,1],  " ", duplicateLess[i,2], sep=" "))
  
  if (grepl(",",duplicateLess[i,2])) { # if there is more than one bestbase
    
    for (el in seq(as.numeric(min(unlist(strsplit(duplicateLess[i,2], ", ")))), as.numeric(max(unlist(strsplit(duplicateLess[i,2], ", ")))), 1)) {
      
      if (nrow(contamination[(grepl(paste("\\<",el,"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]),])!=0) { #if at least one row is found to have the same bestbase, orientation and chromosome
            if (is.na(MouseList[i,2])) { # if there are no values in this cell yet i.e. this is the first base to be processed
            
            MouseList[i,2] <- toString(contamination[(grepl(paste("\\<",el,"\\>", sep="" ), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 5]) #get all mouse id that match this base/ori/chr
            MouseList[i,3] <- toString(contamination[(grepl(paste("\\<",el,"\\>", sep="" ), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 7]) #get all fragments that match this base/ori/chr
            MouseList[i,4] <- toString(contamination[(grepl(paste("\\<",el,"\\>", sep="" ), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 4]) #get all insert id that match this base/ori/chr
            
            } else { # if this is not the first base number to be processed for this row, append the details for the next base
           
            MouseList[i,2] <- paste(MouseList[i,2], toString(contamination[(grepl(paste("\\<",el,"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 5]), sep=", ") #append all mouse id that match this base/ori/chr
            MouseList[i,3] <- paste(MouseList[i,3], toString(contamination[(grepl(paste("\\<",el,"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 7]), sep=", ") #append all fragments that match this base/ori/chr
            MouseList[i,4] <- paste(MouseList[i,4], toString(contamination[(grepl(paste("\\<",el,"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 4]), sep=", ") #append all insert id that match this base/ori/chr
            }
      }
    }
  } else {
    MouseList[i,2] <- toString(contamination[(grepl(paste("\\<",duplicateLess[i,2],"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 5])
    MouseList[i,3] <- toString(contamination[(grepl(paste("\\<",duplicateLess[i,2],"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 7])
    MouseList[i,4] <- toString(contamination[(grepl(paste("\\<",duplicateLess[i,2],"\\>", sep=""), contamination[,2])) & (contamination[6]==duplicateLess[i,6]) & (contamination[1]==duplicateLess[i,1]), 4])

  } 
  if (i%%1000==0){
	cat(paste("chr",chromosome,i,sep=" ")) #progress counter
}
}
MouseListCat <- rbind(MouseListCat, MouseList[!is.na(MouseList[,1]),] ) # append all the rows that have values i.e. column 1 is not NA
write.table(unique(MouseListCat), file=paste(chromosome,"countsInsertsMOUSE.txt",sep=""), row.names=FALSE, sep="\t", quote=FALSE)
}

# #for each string (chr-ori-base) take the list of samples
# SampleList<-matrix(nrow=nrow(duplicateLess),ncol=3)
# for (i in 1:nrow(duplicateLess)){
#   SampleList[i,1]<-toString(paste("chr",duplicateLess[i,6], " ori", duplicateLess[i,1],  " ", duplicateLess[i,2], sep=" "))
#   SampleList[i,2]<-toString(contam[(contam[,7]==SampleList[i,1]), 3])
#   SampleList[i,3]<-toString(contam[(contam[,7]==SampleList[i,1]), 7])
#   
# }

write.table(unique(MouseListCat), file="countsInsertsMOUSE.txt", row.names=FALSE, sep="\t", quote=FALSE)
#write.table(unique(SampleList), file="countsInsertsSAMPLE.txt", row.names=FALSE, sep="\t", quote=FALSE)
