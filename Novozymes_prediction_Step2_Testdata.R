####
#Amudha
#12/3/22
#Novozymes Enzyme prediction - Kaggle competion
###

################## install required packages  ########################
install.packages("bioseq")
library(bioseq)
library(dplyr)


################## Test Data enhacement #######################

#Let us create the sequence weightage for prediction

#Step 1 Add amino acid weightage for each sequence

test_prop <- seq_stat_prop(aa(test$protein_sequence))
test_propdf <-  as.data.frame(do.call(rbind, test_prop))
test <- cbind (test,test_propdf )


# Step 1 Add mutation information to testing set:
# Add mutation information to testing set:
test[,c('type','resid','wt','mut')] <- do.call(rbind,lapply(test$protein_sequence,function(seq){
  # case 1 = wild type:
  if(seq==wtseq){ 
    return(c('WT',-1,NaN,NaN))
    # case 2 = substitution:
  } else if(nchar(seq)==nchar(wtseq)){ 
    i <- mapply(function(x,y) which(x!=y)[1], strsplit(seq,""), strsplit(wtseq,""))
    return(c('SUB',i,substr(wtseq,i,i),substr(seq,i,i)))
    # case 3 = deletion:
  } else if(nchar(seq)<nchar(wtseq)){ 
    wtsub <- substr(wtseq,1,nchar(seq))
    i <- mapply(function(x,y) which(x!=y)[1], strsplit(seq,""), strsplit(wtsub,""))
    return(c('DEL',i,substr(wtseq,i,i),NaN))
  }
}))

head(test)

#Step 2 - add the b factor from pdb file to test data

# Read AlphaFold2 result for wild type sequence:
pdb <- unique(test_pdb[test_pdb$V1=='ATOM',c(6,11)])
colnames(pdb) <- c('resid','b')
head(pdb)


# Add B factor to the testing set:
test_data_withbfactor <- merge(test,pdb,all.x=T)
test_data_withbfactor <- test_data_withbfactor[order(test_data_withbfactor$seq_id),]
head(test_data_withbfactor)


# Download blosum matrix and add score to testing set:
download.file('https://home.cc.umanitoba.ca/~psgendb/doc/local/pkg/ugene/data/weight_matrix/blosum100.txt', destfile="BLOSUM100.txt")
blosum <- read.table('BLOSUM100.txt')
test_data_withbfactor$blosum <- apply(test_data_withbfactor,1,function(x){
  if(x['type']=='WT'){
    return(0)
  } else if(x['type']=='DEL'){
    return(-10)
  } else {
    return(blosum[x['wt'],x['mut']])
  }
})
test_data_withbfactor$blosum[test_data_withbfactor$blosum>0] <- 0

head(test_data_withbfactor)

#### All set with Test data ######

