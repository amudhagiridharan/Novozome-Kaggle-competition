####
#Amudha
#12/3/22
#Novozymes Enzyme prediction - Kaggle competion
#Step3 - Model prediction
###
library(xgboost)
#filter the necessary rows for modeling from train dataframe
traingrp <- data.frame(train_grouped_top25[-c(1,2,4,6,7,8,30)])
head(traingrp)

#filter test rows similar to train data
testxgb <- test_data_withbfactor[-c(1,2,3,5,27,28,29,30,31)]
testxgb$tm <- 0
str(testxgb)
xgb_train = xgb.DMatrix(data=(data.matrix(traingrp[,-2])), label=(traingrp[,2] ))
xgb_test = xgb.DMatrix(data=(data.matrix(testxgb[,-23])), label=(testxgb[,23] ))
xgb <- xgboost(data = xgb_train, max.depth=5,nrounds=25)
pred_xgb = predict(xgb, xgb_test)
head(pred_xgb)
submission <-  data.frame(seq_id = test$seq_id)
submission$tm <- ((-rank(test_data_withbfactor$b)/length(submission$seq_id))+(rank(pred_xgb))/length(submission$seq_id))+(rank(test_data_withbfactor$blosum)/length(submission$seq_id))
head(submission)  
write.csv(submission,"/Users/giridharangovindan/Downloads/submissions_new_file.csv")

