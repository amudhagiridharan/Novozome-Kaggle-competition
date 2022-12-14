####
#Amudha
#12/3/22
#Novozymes Enzyme prediction - Kaggle competion
#Step3 - Model prediction
###
library(xgboost)
#filter the necessary rows for modeling from train dataframe
#traingrp <- data.frame(train_grouped_top25[-c(1,2,4,6,7,8,30)])

set.seed(200)
#filter the necessary rows for modeling from train dataframe
traingrp <- train_grouped_top25 %>% dplyr::select(pH,tm,A,C,G,T,W,S,M,K,R,Y,B,D,H,V,N,E,F,I,L,P,Q,resid,wt,mut)
traingrp$resid = as.integer(traingrp$resid)
# group var won't filter so I'm talking it out here
#traingrp = traingrp[c(-1)]
# makign sure it's a dataframe for the matrix
traingrp = as.data.frame(traingrp)
head(traingrp)

#filter test rows similar to train data
test_data_withbfactor$tm <- 0
testxgb <- test_data_withbfactor %>% dplyr::select(pH,tm,A,C,G,T,W,S,M,K,R,Y,B,D,H,V,N,E,F,I,L,P,Q,resid,wt,mut)
testxgb = as.data.frame(testxgb)
str(testxgb)
xgb_train = xgb.DMatrix(data=(data.matrix(traingrp[,-2])), label=(traingrp[,2] ))
xgb_test = xgb.DMatrix(data=(data.matrix(testxgb[,-2])), label=(testxgb[,2] ))

# Using Cross Validation for best parameter
xgbcv = xgb.cv(data = xgb_train, nfold =25, nrounds =1000, early_stopping_rounds = 40)
opt_iterations = xgbcv$best_iteration

xgb <- xgboost(data = xgb_train, max.depth=5,nrounds=opt_iterations)
#xgb <- xgboost(data = xgb_train, max.depth=5,nrounds=25)
pred_xgb = predict(xgb, xgb_test)
head(pred_xgb)
submission <-  data.frame(seq_id = test$seq_id)
submission$tm <- ((-rank(test_data_withbfactor$b)/length(submission$seq_id))+(rank(pred_xgb))/length(submission$seq_id))+(rank(test_data_withbfactor$blosum)/length(submission$seq_id))

#Deletion ensemble

# Deletion type:
idx <- test_data_withbfactor$type  == 'DEL'
submission[idx,'tm'] <- ((-rank(test_data_withbfactor$b[idx])/length(submission$seq_id))+(rank(test_data_withbfactor$blosum[idx])/length(submission$seq_id)))


# Wild type:
idx <-  test_data_withbfactor$type  == 'WT'
submission[idx,'tm'] <- max(submission$tm) + 1



head(submission)  
write.csv(submission,"/Users/giridharangovindan/Downloads/submissions.csv")

