#Scripts for clustering of insertions

args <- commandArgs(trailingOnly = FALSE)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
cat(myargument)
cat(paste(myargument, ".txt", sep=""))

# Read in the matrix file
sortP <- read.table( paste(myargument, ".txt", sep=""), sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=NULL, skip=1)

#sortP <-read.table( "100_TCTCTCGACT-CTACAGTGCT_sorted_paired.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=NULL, skip=1)

#order by orientation
sortP<-sortP[order(sortP[,4]),]

#take only pos orientation
sortPos<-sortP[(sortP[,4]=="+"),]
#order by chr and LTR pos
sortPos<-sortPos[order(sortPos[,2], as.numeric(sortPos[,3])),]
#take only neg orientation
sortNeg<-sortP[(sortP[,4]=="-"),]
#order by chr and LTR pos
sortNeg<-sortNeg[order(sortNeg[,2], as.numeric(sortNeg[,3])),]

colnames(sortP)<-c("Id","chr","LTRpos","ori","BAMreadLength","pos","pos2","lengthSequence", "Combination Exists")
colnames(sortPos)<-c("Id","chr","LTRpos","ori","BAMreadLength","pos","pos2","lengthSequence", "Combination Exists")
colnames(sortNeg)<-c("Id","chr","LTRpos","ori","BAMreadLength","pos","pos2","lengthSequence", "Combination Exists")

chrlist<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
           "chr21","chr22","chrX","chrY")

totPos<-matrix(nrow=0,ncol=9)
colnames(totPos)<-c("chr","LTRpos","ori","fragments","reads", "startPos", "endPos","pos","pos2")
totNeg<-matrix(nrow=0,ncol=9)
colnames(totNeg)<-c("chr","LTRpos","ori","fragments","reads", "startPos", "endPos","pos","pos2")


#for each chromosome, cluster insertions
for (i in unique(sortPos[,2])){
  
  if (i %in% chrlist){
    
    #####CLUSTERING
    m<-NA
    #distance matrix
    m <- dist(unique(sortPos[(sortPos[,2]==i),3]), method = "euclidean")
    
    if (length(m)==0){
      el<-unique(sortPos[(sortPos[,2]==i),3])
      totPos<-rbind(totPos, c(i, #chr
                              el, #LTR grouped vector
                              "+", #ori
                              nrow(unique(sortPos[(el == sortPos[,3] & sortPos[,2]==i),c(2,3,4,6,7)])), #fragments
                              nrow(sortPos[(el == sortPos[,3] & sortPos[,2]==i),]), #reads
                              el, #startPos
                              el, #endPos
                              toString(unique(sortPos[(el == sortPos[,3] & sortPos[,2]==i),6])), #pos
                              toString(unique(sortPos[(el == sortPos[,3] & sortPos[,2]==i),7])))) #pos2  
    }
    
    else{
      #hierarchical clustering
      hcm <- hclust(m, method="single")
      
       #cut tree for grouping by 10bp
      k<-cutree(hcm,h = 10)
    
    
    
    #####GROUPING
    
    for (jj in unique(k)){
      label<-unique(sortPos[(sortPos[,2]==i),3])
      vecLTR<-label[which(k==jj)]
      
       for(el in vecLTR){
         
        totPos<-rbind(totPos, c(i, #chr
                                toString(vecLTR), #LTR grouped vector
                                "+", #ori
                                nrow(unique(sortPos[(el == sortPos[,3] & sortPos[,2]==i),c(2,3,4,6,7)])), #fragments
                                nrow(sortPos[(el == sortPos[,3] & sortPos[,2]==i),]), #reads
                                min(vecLTR), #startPos
                                max(vecLTR), #endPos
                                toString(unique(sortPos[(el == sortPos[,3] & sortPos[,2]==i),6])), #pos
                                toString(unique(sortPos[(el == sortPos[,3] & sortPos[,2]==i),7])))) #pos2  
        
       
       }
     }
    }
  }
  
}

#same for negative orientation
for (i in unique(sortNeg[,2])){
  
  if (i %in% chrlist){
    #####CLUSTERING
    m<-NA
    #distance matrix
    m <- dist(unique(sortNeg[(sortNeg[,2]==i),3]), method = "euclidean")
    
    if (length(m)==0){
      el<-unique(sortNeg[(sortNeg[,2]==i),3])
      totNeg<-rbind(totNeg, c(i, #chr
                              el, #LTR grouped vector
                              "-", #ori
                              nrow(unique(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),c(2,3,4,6,7)])), #fragments
                              nrow(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),]), #reads
                              el, #startPos
                              el, #endPos
                              toString(unique(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),6])), #pos
                              toString(unique(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),7])))) #pos2  
    }
    
    else{
      #hierarchical clustering
      hcm <- hclust(m, method="single")
      
      #cut tree for grouping by 10bp
      k<-cutree(hcm,h = 10)
      
      
      
      #####GROUPING
      
      for (jj in unique(k)){
        label<-unique(sortNeg[(sortNeg[,2]==i),3])
        vecLTR<-label[which(k==jj)]
        
        for(el in vecLTR){
          
          totNeg<-rbind(totNeg, c(i, #chr
                                  toString(vecLTR), #LTR grouped vector
                                  "-", #ori
                                  nrow(unique(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),c(2,3,4,6,7)])), #fragments
                                  nrow(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),]), #reads
                                  min(vecLTR), #startPos
                                  max(vecLTR), #endPos
                                  toString(unique(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),6])), #pos
                                  toString(unique(sortNeg[(el == sortNeg[,3] & sortNeg[,2]==i),7])))) #pos2  
          
          
        }
      }
    }
  }
  
}


######group the tables totNeg and totPos with (merging fragments (in the same cluster, unique pos2), sum of reads, add number of fragments to each LTR position)

totPos2<-matrix(nrow=0,ncol=ncol(totPos)+1)
colnames(totPos2)<-c(colnames(totPos), "Frequency")


for (i in 1:nrow(totPos)){
  if (grepl(",", totPos[i,2])){
    vec2<-NA
    #vec2<-strsplit(unlist(totPos[i,2]), ", ")
    vec2<-totPos[i,2]
    totPos2<-rbind(totPos2, c(totPos[i,1], #chr
                              totPos[i,2], #LTR grouped vector
                              totPos[i,3], #ori
                              #length(unique(totPos[i,9])), #fragments
                              length(unique(unlist(strsplit(toString(totPos[(totPos[,2]==vec2),9]), ", ")))),
                              sum(as.numeric(totPos[(totPos[,2]==vec2),5])), #reads
                              totPos[i,6], #startPos
                              totPos[i,7], #endPos
                              totPos[i,8],
                              totPos[i,9],
                              toString(paste(totPos[(totPos[,2]==vec2),8], totPos[(totPos[,2]==vec2),4], sep=":"))
                              )) #pos2
  }
  else{
    vec2<-totPos[i,2]
    totPos2<-rbind(totPos2, c(totPos[i,],paste(totPos[(totPos[,2]==vec2),8], totPos[(totPos[,2]==vec2),4], sep=":"))) 
  }
  
}

totPos3<-unique(totPos2[,c(1,2,3,4,5,6,7,10)])
####

totNeg2<-matrix(nrow=0,ncol=ncol(totNeg)+1)
colnames(totNeg2)<-c(colnames(totNeg), "Frequency")


for (i in 1:nrow(totNeg)){
  if (grepl(",", totNeg[i,2])){
    vec2<-NA
    #vec2<-strsplit(unlist(totNeg[i,2]), ", ")
    vec2<-totNeg[i,2]
    totNeg2<-rbind(totNeg2, c(totNeg[i,1], #chr
                              totNeg[i,2], #LTR grouped vector
                              totNeg[i,3], #ori
                              #length(unique(totNeg[i,9])), #fragments
                              length(unique(unlist(strsplit(toString(totNeg[(totNeg[,2]==vec2),9]), ", ")))),
                              sum(as.numeric(totNeg[(totNeg[,2]==vec2),5])), #reads
                              totNeg[i,6], 
                              totNeg[i,7], 
                              totNeg[i,8],
                              totNeg[i,9],
                              toString(paste(totNeg[(totNeg[,2]==vec2),8], totNeg[(totNeg[,2]==vec2),4], sep=":"))
                              ))
  }
  else{
    vec2<-totNeg[i,2]
    totNeg2<-rbind(totNeg2, c(totNeg[i,], paste(totNeg[(totNeg[,2]==vec2),8], totNeg[(totNeg[,2]==vec2),4], sep=":"))) 
  }
  
}

totNeg3<-unique(totNeg2[,c(1,2,3,4,5,6,7,10)])

tot<-matrix(nrow=0, ncol=ncol(totNeg3))

tot<-rbind(tot, totNeg3)
tot<-rbind(tot, totPos3)

write.table(tot, file=paste("TOT",myargument,".txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")
