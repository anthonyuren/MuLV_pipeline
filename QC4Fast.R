
setwd(".")

m3<-read.table("KeepingInserts.txt", sep="\t", row.names=NULL, header=TRUE, stringsAsFactors=FALSE, check.names=F)

#read all inserts file
mCLONfull <- read.table("clon10bp.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE)
mCLONfull <- mCLONfull[order(mCLONfull[,1] , mCLONfull[,5]),]
#read sample table (for plates and wells)
sampletable<-read.table("sampleTable.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE, check.names=F)
mousetable<-matrix(nrow=nrow(m3),ncol=11)

#list of unique mice
MouseList<-sort(unique(unlist(strsplit(m3[,2], ", "))))
chr <- 0 #need an initializer
ori <- "neither" #need an initializer

#for each insert present in more than one mouse i.e. for every row in m3
for (rowM in 1:nrow(m3)){
  
  #for each mouse that contains that insert
  for (mouse in (unlist(strsplit(m3[rowM,2], ", ")))) {
    
    #add a row to mousetable for that insert for that mouse, keep same columns
    mousetable[rowM,1]<-m3[rowM,1]
    mousetable[rowM,2]<-m3[rowM,2]
    mousetable[rowM,3]<-m3[rowM,3]
    mousetable[rowM,4]<-NA
    
    #sample_id contains ALL the samples related to this mouse
    sample_id<-sampletable[(sampletable[3]==mouse),2]
    #vector of plates
    plateV<-c()
    #get chr for this duplicated insert
    chrold <- chr
    chr <- gsub("[^\\d]+", "",unlist(strsplit(m3[rowM,1], "  "))[1], perl=TRUE)
      
    #get ori for this duplicated insert
    oriold <- ori
    if (grepl( "\\-", (unlist(strsplit(m3[rowM,1], "  "))[2]))){
      ori<-"-"
    }
    else if (grepl( "\\+", (unlist(strsplit(m3[rowM,1], "  "))[2]))){
      ori<-"+"
    }
    
    if (chr != chrold | ori != oriold) { # refresh the mCLON list if the chromosome and/or orientation of this row differs from the previous row
      mCLON <- mCLONfull[which(mCLONfull[,1]==chr & mCLONfull[,5]==ori),]
    }
    
    #get bestbase(s) for this duplicated insert
    bestbase<-gsub(" ", "",unlist(strsplit(unlist(strsplit(m3[rowM,1], "  "))[3], ", ")))
    #for each bestbase take the samples and fragment numbers for that insert in this mouse
    sampleR<-c()
    fragments<-c()
  
    for (currentBestBase in seq(as.numeric(min(bestbase)), as.numeric(max(bestbase)), 1)){ # for every base between the minimum best base and maximum best base
      if (length(mCLON[(mCLON[,13]==mouse & mCLON[,12] %in% sample_id & mCLON[,1]==chr & mCLON[,5]==ori & 
                                                        grepl(paste("\\<",currentBestBase,"\\>", sep=""), mCLON[,20])),12]) > 0) { #if a row for this currentBestBase exists (it may not since it increments by 1 between the min and max value)
        
        sampleR<-c(sampleR, toString(mCLON[(mCLON[,13]==mouse & mCLON[,12] %in% sample_id & mCLON[,1]==chr & mCLON[,5]==ori & 
                                                              grepl(paste("\\<",currentBestBase,"\\>", sep=""), mCLON[,20])),12])) # get the sample numbers
        
        fragments<-c(fragments, toString(mCLON[(mCLON[,13]==mouse & mCLON[,12] %in% sample_id & mCLON[,1]==chr & mCLON[,5]==ori & 
                                                              grepl(paste("\\<",currentBestBase,"\\>", sep=""), mCLON[,20])),18])) # get the fragment numbers
      }
    }
    
    #assign a plate/s number based on the sample number, store as list plateV, each plate is 96 samples therefore can get the plate number using "ceiling(as.numeric(sampleR)/96)"
    
    #if only one sample
    if (length(unlist(strsplit(sampleR, ", ")))==1){
      plate <- ceiling(as.numeric(sampleR)/96)
      #well<-sampletable[(sampletable[3]==mouse),13]
      plateV<-plate
    }
    else{
      #if more samples, take plates of each sample
      for (i in as.numeric(unlist(strsplit(sampleR, ", ")))){
      plate <- ceiling(i/96)
      #well<-sampletable[(sampletable[3]==mouse),13]
      plateV<-c(plateV, plate) 
      }
    }
    
    # keep a string for each mouse-sample-plate-fragment
    if (is.na(mousetable[rowM,5])){
      mousetable[rowM,5]<- toString(paste("M(",mouse,") : S(",(unlist(strsplit(sampleR, ", "))),") : P(",plateV,") : F(",(unlist(strsplit(fragments, ", "))),")"), sep="")
    }
    else{
      mousetable[rowM,5]<- paste(mousetable[rowM,5], "/", toString(paste("M(",mouse,") : S(",(unlist(strsplit(sampleR, ", "))),") : P(",plateV,") : F(",(unlist(strsplit(fragments, ", "))),")"), sep=""))
    }  
    
    #list of fragments
    if (is.na(mousetable[rowM,6])){
      mousetable[rowM,6]<- toString(unlist(strsplit(fragments, ", ")))
    }
    else{
      mousetable[rowM,6]<- paste(mousetable[rowM,6], "/", toString(unlist(strsplit(fragments, ", "))))
    } 
    
    #list of plates
    if (is.na(mousetable[rowM,7])){
      mousetable[rowM,7]<- toString(plateV)
    }
    else{
      mousetable[rowM,7]<- paste(mousetable[rowM,7], "/", toString(plateV))
    }
  } 
  
#print(mousetable[rowM,])  # e.g. 
  # [1] "chr 1  ori -   4623109"                                                         
  # [2] "3197, 3032"                                                                     
  # [3] "1, 1"                                                                           
  # [4] NA                                                                               
  # [5] "M( 3197 ) : S( 313 ) : P( 4 ) : F( 1 ) / M( 3032 ) : S( 505 ) : P( 6 ) : F( 1 )"
  # [6] "1 / 1"                                                                          
  # [7] "4 / 6"
  # if any mouse has two samples with the same fragment these are separated by comma e.g.
  # M( 3230 ) : S( 615 ) : P( 7 ) : F( 1 ), M( 3230 ) : S( 638 ) : P( 7 ) : F( 1 )
  # after row is created
  # delete duplicates
  
  temp1<-""
  temp2<-""
  temp3<-""
  
  # if the insert is only present in one mouse i.e. if there are no "/" in that row 
  if (!(grepl ("/", mousetable[rowM,5]))) {
    mousetable[rowM,5]<- toString(unique(unlist(strsplit(mousetable[rowM,5], ", "))))
    
    im<-which((unlist(strsplit(mousetable[rowM,5], ", "))) == (unique(unlist(strsplit(mousetable[rowM,5], ", ")))))
    mousetable[rowM,6]<- toString(unlist(strsplit(mousetable[rowM,6], ", "))[im]) # take the values of column 6 that are multiple inserts withing the same mouse and split them
    mousetable[rowM,7]<- toString(unlist(strsplit(mousetable[rowM,7], ", "))[im]) # take the values of column 7 that are multiple inserts withing the same mouse and split them
    print(c("a row representing a single mouse",mousetable[rowM,5]))
  }
  
  # else the insert is present in more than one mouse, take each string for each mouse-sample-plate-fragment
  
  else { 
    for (ju in unlist(strsplit(mousetable[rowM,5], " / "))) { 
      # if more samples per mouse
      if (grepl (",", ju)) {
        im1<-which((unlist(strsplit(mousetable[rowM,5], " / "))) == ju)
        ju2<-unique(unlist(strsplit(ju, ", ")))
        
        Fr2<-c()
        Pl2<-c()
        # for each index of lists
        for (el in ju2){
          im2<-which(unlist(strsplit(ju, ", ")) == el)
          Fr2<-c(Fr2, unique(as.numeric(unlist(strsplit(unlist(strsplit(mousetable[rowM,6], " / "))[im1], ", ")))[im2]))
          Pl2<-c(Pl2, unique(as.numeric(unlist(strsplit(unlist(strsplit(mousetable[rowM,7], " / "))[im1], ", ")))[im2]))
        }
                
        # if string is empty
        if (temp1==""){
          temp1<- toString(unique(ju2))
          temp2<- toString(Fr2)
          temp3<- toString(Pl2)
        }
        else{
          # concat
          temp1<- paste(temp1, toString(unique(ju2)), sep=" / ")
          temp2<- paste(temp2, toString(Fr2), sep=" / ")
          temp3<- paste(temp3, toString(Pl2), sep=" / ")
        }
      }
      # if one sample per mouse
      else{
        # take the index of 
        im3<-which((unlist(strsplit(mousetable[rowM,5], " / "))) == ju)
        # if string is empty
        if (temp1==""){
          temp1<- toString(ju)
          temp2<- toString(unlist(strsplit(mousetable[rowM,6], " / "))[im3])
          temp3<- toString(unlist(strsplit(mousetable[rowM,7], " / "))[im3])
        }
        
        else{
          # concat
          temp1<- paste(temp1, toString(ju), sep=" / ")
          temp2<- paste(temp2, toString(unlist(strsplit(mousetable[rowM,6], " / "))[im3]), sep=" / ")
          temp3<- paste(temp3, toString(unlist(strsplit(mousetable[rowM,7], " / "))[im3]), sep=" / ")
        }
      }
    }
    # replace the row
    mousetable[rowM,5]<- toString(temp1)
    mousetable[rowM,6]<- toString(temp2)
    mousetable[rowM,7]<- toString(temp3)
  }
  
  print(c(mousetable[rowM,]))   # e.g. 
  
  ###################################                                                 
  # vector of fragments used
  frV<-c()
  # take each plate (FIRST CHECK, PLATES: FIND THE EARLIEST)
  plateUnl<-unlist(strsplit(mousetable[rowM,7], " / "))
  
  # if more sample per mouse
  if (grepl(",",mousetable[rowM,7])){
    
    #index of mouse/mice with multiple samples
    indStrComPl<-which(grepl(",",unlist(strsplit(mousetable[rowM,7], " / "))))
    
    # split each list and take the minimum
    
    # if only one mouse has more samples
    if (length(indStrComPl)==1){
      
      # minimum plate is the minimum of the list of plates
      PlateMin<-min(as.numeric(unlist(strsplit(unlist(strsplit(mousetable[rowM,7], " / "))[indStrComPl], ", "))))
      #new vector with new plates (one per mouse)
      newVP<-unlist(strsplit(mousetable[rowM,7], " / "))
      #replace at right place the list with the minimum of the list
      newVP[indStrComPl]<-PlateMin
      #take the minimum of the new vector
      PlateMin2<-min(as.numeric(newVP))
      #take the index (that will correspond to the mouse index)
      mouseEarlIndex<-which(PlateMin2 == newVP)
    }
    
    # otherwise split more mice and take all the minimums
    else{
      # for each index of lists
      for (el in indStrComPl){
        #minimum plate is the minimum of the list of plates
        PlateMin<-min(as.numeric(unlist(strsplit(unlist(strsplit(mousetable[rowM,7], " / "))[el], ", "))))
        #new vector with new plates (one per mouse)
        newVP<-unlist(strsplit(mousetable[rowM,7], " / "))
        #replace at right place the list with the minimum of the list
        newVP[el]<-PlateMin
      }
      #take the minimum of the new vector
      PlateMin2<-min(as.numeric(newVP))
      #take the indexes (that will correspond to the mouse indexes)
      mouseEarlIndex<-which(PlateMin2 == newVP)
    
    }
  }
  
  # if only one sample per mouse, take the index of the mouse with minimum plate
  else{
    mouseEarlIndex<-which(min(as.numeric(unlist(strsplit(mousetable[rowM,7], " / ")))) == (as.numeric(unlist(strsplit(mousetable[rowM,7], " / ")))))
  }
  
  ##TAKE THE EARLIEST MOUSE (index/indexes ready)
  
  # if there is only one insert/mouse on the earliest plate choose that one
  if (length(mouseEarlIndex)==1){
    #update row with mouse number
    mousetable[rowM,8]<-(as.numeric(unlist(strsplit(mousetable[rowM,2], ", "))))[mouseEarlIndex]
    
  }
  # if the earliest plate has more than one mouse with the same insert compare the relative abundance of inserts
  else if (length(mouseEarlIndex)>1){
    
    # default is "-"
    mousetable[rowM,8]<-toString("-")
    # for each mouse in the same plate
    for (mID in mouseEarlIndex){
      
      # if the plate contains more fragments (samples)
      if (grepl(",",mousetable[rowM,6])){
        
        indStrCom<-which(grepl(",",unlist(strsplit(mousetable[rowM,6], " / "))))
        
        # if only one mouse has more samples
        if (length(indStrCom)==1){
          
          # maximum fragment is the maximum of the list of fragments
          FragMax<-max(as.numeric(unlist(strsplit(unlist(strsplit(mousetable[rowM,6], " / "))[indStrCom], ", "))))
          # new vector with new fragments (one per mouse)
          newVF<-unlist(strsplit(mousetable[rowM,6], " / "))
          # replace at right place the list with the maximum of the list
          newVF[indStrCom]<-FragMax
        }
        
        # otherwise split more mice and take all the maximums
        else{
          # for each index of lists
          for (el in indStrCom){
            # maximum fragment is the maximum of the list of fragments
            FragMax<-max(as.numeric(unlist(strsplit(unlist(strsplit(mousetable[rowM,6], " / "))[el], ", "))))
            # new vector with new fragments (one per mouse)
            newVF<-unlist(strsplit(mousetable[rowM,6], " / "))
            # replace at right place the list with the maximum of the list
            newVF[el]<-c(newVF, FragMax)
          }
        }
        frV<-newVF
      }
      
      # if one fragment per plate
      else{
        # take vector of fragments of earliest MICE
        frV<-c(frV, (as.numeric(unlist(strsplit(mousetable[rowM,6], " / "))))[mID])
      }
    }
    
    ratio <- 0
    frV <- as.numeric(frV)
    
    # take the 2 max values and then ratio
    el1<-max(frV)
    el2<-max(frV[frV!=max(frV)])
    
    ratio<-as.numeric(el1)/as.numeric(el2)
    
    # if either fragment is 10 fold more abundant than the other (e.g. frag 211 vs frag 19) keep the more abundant one
    if (ratio < 0.1){
      indexfrag<-which(el2==frV) 
      indexR<-mouseEarlIndex[indexfrag]
      mousetable[rowM,8]<-toString(as.numeric(unlist(strsplit(mousetable[rowM,2], ", ")))[indexR])
      mousetable[rowM,9]<-ratio
    } else if (ratio > 10){
      indexfrag<-which(el1==frV) 
      indexR<-mouseEarlIndex[indexfrag]
      mousetable[rowM,8]<-toString(as.numeric(unlist(strsplit(mousetable[rowM,2], ", ")))[indexR])
      mousetable[rowM,9]<-ratio
    }
    else{ # if neither is 10 fold more abundant than the other keep neither
      mousetable[rowM,9]<-"0"
    }
  }
  
  # record the winning fragment information (sample, fragment, plate) in the final column
  FinalIndex<-which(mousetable[rowM,8]==(unlist(strsplit(mousetable[rowM,2], ", ")))) 
  mousetable[rowM,10]<-toString(unlist(strsplit(mousetable[rowM,6], " / "))[FinalIndex])
  mousetable[rowM,11]<-toString(unlist(strsplit(mousetable[rowM,5], " / "))[FinalIndex])
  cat(rowM)
}

colnames(mousetable)<-c("Insert","Mouse","FragmentsByMouse","InsertIdList","PlatesSamplesFragments","FragmentsBySample","Plates","EarliestMouse","Ratio",
                        "MouseInfoFragments", "AllmouseInfo")

write.table(mousetable[,c(1,2,3,5,6,7,8,9,10,11)], file="InsertsAndPlates.txt", quote=FALSE, sep="\t", row.names=FALSE)