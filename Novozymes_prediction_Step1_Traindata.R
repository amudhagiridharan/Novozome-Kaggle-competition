####
#Amudha
#12/3/22
#Novozymes Enzyme prediction - Kaggle competion
###

################## install required packages  ########################
install.packages("bioseq")
library(bioseq)
library(dplyr)

################## Data load ################## 
# read the train, train update, test and pdb data sets as provided by competition hosts

train <-  read.csv(file="/Users/giridharangovindan/Downloads/train.csv", header=T)
train_update <-  read.csv(file="/Users/giridharangovindan/Downloads/train_updates_20220929.csv", header=T)
test <-  read.csv(file="/Users/giridharangovindan/Downloads/test.csv", header=T)
test_pdb <- read.table("/Users/giridharangovindan/Downloads/wildtype_structure_prediction_af2.pdb",fill=T)


# Wild type sequence provided in the Kaggle "Dataset Description":
wtseq <- 'VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK'


head(train)
str(train) #31390 obs
################## Data cleaning #######################

# Competion host mentions that train update has the fixes for data issues in train dataset
# NAs in train data has to be deleted as per competion host instruction


# let us remove all the rows with data issues as advised by host
train <- train[!(train$seq_id %in% train_update$seq_id),]
str(train) # 28956 obs

#Classify the train update data and identify the rows to be used
train_update_null <-  subset(train_update,is.na(train_update$tm))
head(train_update_null) # we are going to ignore this dataset. using only for validation purpose
train_update_swap <- subset(train_update,!is.na(train_update$tm))
head(train_update_swap)

#use the train_update_swap dataframe to add back the rows with swapped pH & tm
train <- rbind(train, train_update_swap)
str(train) # 28981 obs


# check for duplicates
train[duplicated(train),] # no duplicates

################## Train Data Wrangling #######################

#on protein sequence analysis, we find that it is important to find the sequence count
# and cluster the similar mutants for further modeling

#Let us try to cluster the mutants which seem similar


#Step 1 add sequence char count to train data

train$charcnt <- seq_nchar(aa(train$protein_sequence))

#Step 2: Let us filter out the protein sequences which have less frequency of occurrance
# we might not be use them in modeling as the datapoints wouldn't be sufficient.
#Threshold set is 20

seq_freq_threshold = 20

train_filtered <- train %>% 
  group_by(charcnt) %>% 
  filter(n() >= seq_freq_threshold)

#Step 3 add cluster number to train data

train_filtered$clusternum = -1
train_filtered$wildtype = ''
for (i in unique(train_filtered$charcnt)) {
  train_filtered[train_filtered$charcnt == i, ]$clusternum <- seq_cluster(aa(train_filtered[train_filtered$charcnt == i, ]$protein_sequence))
}



# Step 4 Identify the wild type for train data !! long running query

#for each protein sequence length(charcnt) and cluster number, loopthrough and 
#find the protein sequence consensus (wildtype)

for (i in unique(train_filtered$charcnt)) {
  for (j in unique (train_filtered[train_filtered$charcnt == i, ]$clusternum)){
    train_filtered[(train_filtered$charcnt == i & train_filtered$clusternum == j) , ]$wildtype <- seq_consensus(aa(train_filtered[(train_filtered$charcnt == i & train_filtered$clusternum == j), ]$protein_sequence))
  }
}

# Step 5 Order it by charcnt for easiness

train_filtered <- train_filtered[order(train_filtered$charcnt),]


#Step 6 Add amino acid weightage for each sequence

train_filtered_prop <- seq_stat_prop(aa(train_filtered$protein_sequence))
train_filtered_propdf <-  as.data.frame(do.call(rbind, train_filtered_prop))
train_filtered <- cbind (train_filtered,train_filtered_propdf )

# Step 7 Add the group info to identify potential model training data

#grouping train data by datasource, PH , charcnt and cluster
#note a set of charcnt and cluster belong to one wild type

train_grouped <- train_filtered %>%     # Create ID by group
  group_by(charcnt,clusternum,pH,data_source) %>%
  dplyr::mutate(group = cur_group_id())


head(train_grouped)

# Step 8 - Final step to filter out only the top 25 groups to train our model

train_grouped_top25 <- train_grouped %>% 
  group_by(group) %>% 
  filter(n() > 25)

train_grouped_top25 <- as.data.frame(train_grouped_top25)

train_grouped_top25[,c('type','resid','wt','mut')] = NA
for(i in 1:nrow(train_grouped_top25)){
  if(train_grouped_top25$protein_sequence[i]==train_grouped_top25$wildtype[i]){ 
    train_grouped_top25[i ,c('type','resid','wt','mut')] = as.list(c('WT',-1,NaN,NaN))
    # case 2 = substitution:
  }
  else if(nchar(train_grouped_top25$protein_sequence[i])==nchar(train_grouped_top25$wildtype[i])){ 
    P <- mapply(function(x,y) which(x!=y)[1], strsplit(train_grouped_top25$protein_sequence[i],""), strsplit(train_grouped_top25$wildtype[i],""))
    train_grouped_top25[i ,c('type','resid','wt','mut')]=as.list(c('SUB',P,substr(train_grouped_top25$wildtype[i],P,P),substr(train_grouped_top25$protein_sequence[i],P,P)))
    # case 3 = deletion:
  } else if(nchar(train_grouped_top25$protein_sequence[i])<nchar(train_grouped_top25$wildtype[i])){ 
    wtsub <- substr(train_grouped_top25$wildtype[i],1,nchar(train_grouped_top25$protein_sequence[i]))
    P <- mapply(function(x,y) which(x!=y)[1], strsplit(train_grouped_top25$protein_sequence[i],""), strsplit(wtsub,""))
    train_grouped_top25[i ,c('type','resid','wt','mut')]=as.list(c('DEL',P,substr(train_grouped_top25$wildtype[i],P,P),NaN))
  }
} 



head(train_grouped_top25)


###### All set with Training data ######
