rm(list=ls())
#Variable importance in wavelet filtered data when put to machine learning
#load("./RData/006_FinalEEGdataForMLRelPower_32Hz.rda")
load("./RData/006_FinalEEGdataForMLRelPower_64Hz.rda")
source("./func-room.R")
colnames(Final_EEGData_forML)[72]<-"class"
colnames(Final_EEGData_forML)[72]
colnames(Final_EEGData_forML)[1]<-c("subjectName")
colnames(Final_EEGData_forML)[1]
eeg.df_02<-cbind(Final_EEGData_forML[,-c(1)],Final_EEGData_forML[,c(1)])
head(eeg.df_02[,71])
colnames(eeg.df_02)[72]<-"subjectName"
colnames(eeg.df_02)<-gsub(pattern = " ",x = colnames(eeg.df_02),replacement = "_")

ConfusionMatrix_multiRF <- function(data){
  require(caret)
  data<-eeg.df_02[,-72]  
  train.test <- getTrainTest(
    subject.data = data,
    class.name = 'class',
    seed = 1094
  )
  fold.accuracy <- c()
  test.accuracy <- c()
  rf_model_list<-vector(mode="list")
  for(i in 1:10){
  #i=1
  message('fold ',i)
  train.idx <- unlist(train.test$train$fold10[-i])
  train.dat <- train.test$train$data[train.idx,]
  test.idx <- unlist(train.test$train$fold10[i])
  idx.state <- grep('class',colnames(train.test$train$data))
  test.dat <- train.test$train$data[test.idx,-idx.state]
  test.y <- as.factor(train.test$train$data[test.idx,idx.state])
  train.dat$class <- as.factor(train.dat$class)
  library(randomForest)
  set.seed(8596+i)
  rf.model <- randomForest(
    class~.,
    data=train.dat,
    ntree=2000)
  rf_model_list[[i]]<-rf.model
  
  varImpPlot(x = rf.model)
  }
  predicted <- predict(rf.model, test.dat,type = 'class')  # predicted scores
  library(caret)
  cm <- confusionMatrix(data=predicted,reference = test.y)
  return(rf_model_list)    
  #fold.accuracy[i] <- cm$overall[[1]]
  #predicted.test <- predict(rf.model, train.test$test$features,type = 'class')  # predicted scores
  #test.accuracy[i] <- confusionMatrix(
  #  data=predicted.test,
  #  reference = as.factor(train.test$test$y)
  #)$overall[[1]]
  #}
  #return(list(test = test.accuracy,
  #           fold10 = fold.accuracy))
  }


cm_eeg<-ConfusionMatrix_multiRF( data= eeg.df[,-86])

fold10_cv<-lapply(cm_eeg,varImp)
testing<-plyr::ldply(fold10_cv,.fun = t)
col_median<-apply(testing,2,median)
median_order<-order(col_median,decreasing = TRUE)
ordered_testing<-testing[,median_order]

melted_ordtest<-reshape::melt(ordered_testing)
head(melted_ordtest)

library(ggplot2)
#model.m <- reshape::melt(model.perf)
#model.m$value <- model.m$value*100
colnames(melted_ordtest) <- c("ElectrodeBand","MeanGini")
p <- ggplot(melted_ordtest,aes(x=ElectrodeBand,y=MeanGini)) +
  geom_boxplot(aes(fill=ElectrodeBand)) +
  xlab('Electrode EEGBand') + ylab('Mean Gini') +
  ggtitle('Wavelet PSD RF Variable Importance 10 Fold CV') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
png('./figures/000-VarImpRF.png',
    width = 3500,height = 1200,res = 200)
p
dev.off()
