setwd(".")

#read input file
m <- read.table("countsInsertsMOUSE.txt", sep="\t", row.names=NULL, header=TRUE, stringsAsFactors=FALSE)
m <- cbind(m, NA)
m <- cbind(m, NA)

# take unique mouse and fragments
for (i in 1:nrow(m)){
  vecUniqueMouse<-as.numeric(c(unique(unlist(strsplit(m[i,2], ", ")))))
  vecMouse<-as.numeric(c(unlist(strsplit(m[i,2], ", "))))
  vecFragments<-as.numeric(c(unlist(strsplit(m[i,3], ", "))))
  vecUniqueFragments<-c()
  flag=FALSE
  
  # take the fragments of unique mice
  for (el in vecUniqueMouse){
    # if there's more than one mouse

    if ((length(which(el==vecMouse)) > 1)) {  
       vecUniqueFragments<-c(vecUniqueFragments, sum(vecFragments[which(el==vecMouse)])) 

    } else {
      if (length(which(el==vecMouse)) == 1){
       vecUniqueFragments<-c(vecUniqueFragments, vecFragments[which(el==vecMouse)]) 
     }
   }
  }
  
  m[i,6]<-toString(vecUniqueFragments)
  m[i,5]<-toString(vecUniqueMouse)
}

m2<-m[,c(1,5,6,4)]
colnames(m2)<-c("Insert","Mouse", "Fragments", "InsertID")

write.table(m2, file="KeepingInserts.txt", quote=FALSE, sep="\t", row.names=FALSE)
