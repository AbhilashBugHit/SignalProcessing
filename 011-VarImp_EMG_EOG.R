rm(list=ls())

source("./func-room.R")
load("./RData/EMG.df_app-03.rda")

ConfusionMatrix_multiRF <- function(data){
  require(caret)
  #data<-emg.df_02[,-ncol(emg.df_02)]  
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

cm_eeg<-ConfusionMatrix_multiRF( data= emg.df_02[,-ncol(emg.df_02)])

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
colnames(melted_ordtest) <- c("ElectrodeStats","MeanGini")

  p <- ggplot(melted_ordtest,aes(x=ElectrodeStats,y=MeanGini)) +
  geom_boxplot(aes(fill=ElectrodeStats)) +
  xlab('Electrode EEGBand') + ylab('Mean Gini') +
  ggtitle('EMG/EOG interval State Classification Variable Importance 10 Fold CV') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

png('./figures/011-VarImpRF_EMGEOG_timeinterval.png',width = 3500,height = 1200,res = 200)
p
dev.off()
