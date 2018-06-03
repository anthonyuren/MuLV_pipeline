
load("QC5.RData")
finalt2<-matrix(nrow=0, ncol=ncol(Allin))

# rbind for large tables gets very slow so best to break it up into groups of 100,000 lines
# the loops can be run sequentially or concurrently and the tables can be concatenated at the end

for (i in 1:100000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finalta.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 100001:200000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finaltb.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 200001:300000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finaltc.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 300001:400000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finaltd.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 400001:500000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finalte.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 500001:600000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finaltf.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 600001:700000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finaltg.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 700001:800000){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finalth.txt", row.names=FALSE, sep="\t", quote=FALSE)

finalt2<-matrix(nrow=0, ncol=ncol(Allin))
for (i in 800001:nrow(Allin)){ ##divided in a-b-c-d
  if (length(itrtot1[(itrtot1[,1] == Allin[i,1]) & (itrtot1[,2]==Allin[i,5]) & (itrtot1[,3]==Allin[i,20]) &
                     (itrtot1[,4]==Allin[i,13]) & (itrtot1[,5]==Allin[i,12]), ])==0){
    finalt2<-rbind(finalt2, Allin[i,])
  }
}  

write.table(finalt2, file="finalti.txt", row.names=FALSE, sep="\t", quote=FALSE)
