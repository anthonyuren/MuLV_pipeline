# remove duplicate inserts keeping only the earliest cloned insert
IP <- read.table("InsertsAndPlates-2018.txt", sep="\t", row.names=NULL, header=TRUE, stringsAsFactors=FALSE)
#all the inserts in the db
Allin <- read.table("clon10bp.txt", sep="\t", row.names=NULL, header=FALSE, stringsAsFactors=FALSE)

# create a new matrix with inserts to remove
# find all columns from InsertsAndPlates that did not have an insert to keep 
# i.e. the plate where that insert first appeared did not have one insert that is 10 fold more abundant than the others.
# grep for anything lacking a digit i.e. !(grepl("\\d") in the column IP$EarliestMouse

itr<-IP[(!(grepl("\\d", IP[,7]))),]


#take only "Insert" "PlatesSamplesFragments" "AllmouseInfo" (last should be empty since these rows are not keepers)
itr<-itr[,c(1,4,10)]

itr1<-matrix(nrow=0, ncol=ncol(itr)-1)

# list these mice as one line each i.e. separate the joint lines with multiple inserts at the same position indo individual lines.

for (i in 1:nrow(itr)){
  for (el in unlist(strsplit(itr[i,2], " / "))){
    if (!(el == itr[i,3])){
      itr1<-rbind(itr1, c(itr[i,1],el)) 
    }
  }
}

itr1<-cbind(itr1, NA)
itr1<-cbind(itr1, NA)
itr1<-cbind(itr1, NA)
itr1<-cbind(itr1, NA)
itr1<-cbind(itr1, NA)

# split the insert details into chrom ori and bestbase individual positions

for (i in 1:nrow(itr1)){
  #chr
  itr1[i,3]<-(gsub("[^0-9]", "", unlist(strsplit(itr1[i,1], "   ")), ""))[1]
  #ori
  itr1[i,4]<-(gsub("ori ", "", unlist(strsplit(itr1[i,1], "  ")), ""))[2]
  
  if (grepl(",", itr1[i,1])){
    #insert
    itr1[i,5]<-toString(unlist(strsplit(itr1[i,1], "   "))[2])
  }
  else{
    #insert
    itr1[i,5]<-(gsub("[^0-9]", "", unlist(strsplit(itr1[i,1], "   ")), ""))[2]
  }
  
  if (grepl(",", itr1[i,2])){
    for (mm in unlist(strsplit(itr1[i,2], ", "))){
      if(is.na(itr1[i,6])){
        #mouse
        itr1[i,6]<-(gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[1]
        #sample
        itr1[i,7]<-(gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[2]
      }
      else{
        #mouse
        itr1[i,6]<-paste(itr1[i,6], (gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[1], sep=", ")
        #sample
        itr1[i,7]<-paste(itr1[i,7], (gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[2], sep=", ")
      }
    }
  }
  else{
    #mouse
    itr1[i,6]<-(gsub("[^0-9]", "", unlist(strsplit(itr1[i,2], " : ")), ""))[1]
    #sample
    itr1[i,7]<-(gsub("[^0-9]", "", unlist(strsplit(itr1[i,2], " : ")), ""))[2]
  }
}

####

# find all columns from InsertsAndPlates that have an insert to keep i.e. the plate where that insert first appeared had one insert that is 10 fold more abundant than all the others on that plate.
# grep for anything lacking a digit i.e. !(grepl("\\d") in the column IP$EarliestMouse

#mouse to keep
itk<-IP[(grepl("\\d", IP[,7])),]

#take only insert info and mouse to keep
itk<-itk[,c(1,4,10)]

####

itr2<-matrix(nrow=0, ncol=ncol(itk)-1)

# list these mice as one line each i.e. separate the joint lines with multiple inserts at the same position indo individual lines.
# remove all the insert details except the insert/s from the mouse that will be kept

for (i in 1:nrow(itk)){
  for (el in unlist(strsplit(itk[i,2], " / "))){
   if (!(el == itk[i,3])){
    itr2<-rbind(itr2, c(itk[i,1],el)) 
   }
  }
}

#read the table to extract sample-mouse-insert to delete

itr2<-cbind(itr2, NA)
itr2<-cbind(itr2, NA)
itr2<-cbind(itr2, NA)
itr2<-cbind(itr2, NA)
itr2<-cbind(itr2, NA)

for (i in 1:nrow(itr2)){
  #chr
  itr2[i,3]<-(gsub("[^0-9]", "", unlist(strsplit(itr2[i,1], "   ")), ""))[1]
  #ori
  itr2[i,4]<-(gsub("ori ", "", unlist(strsplit(itr2[i,1], "  ")), ""))[2]
  
  if (grepl(",", itr2[i,1])){
    #insert
    itr2[i,5]<-toString(unlist(strsplit(itr2[i,1], "   "))[2])
  }
  else{
    #insert
    itr2[i,5]<-(gsub("[^0-9]", "", unlist(strsplit(itr2[i,1], "   ")), ""))[2]
  }
  
  if (grepl(",", itr2[i,2])){
    for (mm in unlist(strsplit(itr2[i,2], ", "))){
      if(is.na(itr2[i,6])){
        #mouse
        itr2[i,6]<-(gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[1]
        #sample
        itr2[i,7]<-(gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[2]
      }
      else{
        #mouse
        itr2[i,6]<-paste(itr2[i,6], (gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[1], sep=", ")
        #sample
        itr2[i,7]<-paste(itr2[i,7], (gsub("[^0-9]", "", unlist(strsplit(mm, " : ")), ""))[2], sep=", ")
      }
    }
  }
  else{
    #mouse
    itr2[i,6]<-(gsub("[^0-9]", "", unlist(strsplit(itr2[i,2], " : ")), ""))[1]
    #sample
    itr2[i,7]<-(gsub("[^0-9]", "", unlist(strsplit(itr2[i,2], " : ")), ""))[2]
  }
}

itrtot<-rbind(itr1,itr2)

itrtot<-itrtot[,c(3:7)]

itrtot1<-matrix(nrow=0, ncol=ncol(itrtot))

# split rows with a comma up so there is one row per position

for (i in 1:nrow(itrtot)){
  count<-0
  for (el in unlist(strsplit(itrtot[i,4], ", "))){
    count<-count+1
    itrtot1<-rbind(itrtot1, c(itrtot[i,c(1:3)], el, unlist(strsplit(itrtot[i,5], ", "))[count])) 
  }
}

Allin<-Allin[order(Allin[,1], Allin[,5], Allin[,20]),]

save(itrtot1, Allin, file = "QC5.RData")
