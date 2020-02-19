rm(list=ls())
#setwd("~/Documents/Sandbox/ceeri-igib-ai-project/")
load('./RData/009_ML_df.rda')
source('./func-room.R')

dir.create("./data")
dir.create('./data/emg-models-app03',showWarnings = TRUE)

ML_df$SubjectName<-as.character(ML_df$SubjectName)
#colnames(Final_EEGData_forML)
Final_EMGData_forML<-ML_df

colnames(Final_EMGData_forML)[1]<-c("subjectName")

emg.df_02<-cbind(Final_EMGData_forML[,-c(1)],Final_EMGData_forML[,c(1)])
colnames(emg.df_02)[57]<-c("class")
#colnames(eeg.df_02)
colnames(emg.df_02)[ncol(emg.df_02)]<-"subjectName"
colnames(emg.df_02)<-gsub(pattern = " ",x = colnames(emg.df_02),replacement = "_")

save(emg.df_02,file = "./RData/EMG.df_app-03.rda",compress=TRUE)

rf.models <- multiRF(data = emg.df_02[,-(ncol(emg.df_02))])
save(
  rf.models,
  file = './data/emg-models-app03/rf.models.RData'
)

svm.radial <- multiSVM(
  data = emg.df_02[,-(ncol(emg.df_02))],
  kernel.type = 'radial'
)

save(
  svm.radial,
  file = './data/emg-models-app03/svm.radial.RData'
)

svm.polynomial <- multiSVM(
  data = emg.df_02[,-(ncol(emg.df_02))],
  kernel.type = 'polynomial'
)
save(
  svm.polynomial,
  file = './data/emg-models-app03/svm.polynomial.RData'
)

svm.sigmoid <- multiSVM(
  data = emg.df_02[,-(ncol(emg.df_02))],
  kernel.type = 'sigmoid'
)
save(
  svm.sigmoid,
  file = './data/emg-models-app03/svm.sigmoid.RData'
)
svm.linear <- multiSVM(
  data = emg.df_02[,-(ncol(emg.df_02))],
  kernel.type = 'linear'
)
save(
  svm.linear,
  file = './data/emg-models-app03/svm.linear.RData'
)
knn.model <- multiKNN(
  data = emg.df_02[,-(ncol(emg.df_02))]
)
save(
  knn.model,
  file = './data/emg-models-app03/knn.model.RData'
)

rm(list=ls())
load('./data/emg-models-app03/rf.models.RData')
load('./data/emg-models-app03/svm.radial.RData')
load('./data/emg-models-app03/svm.sigmoid.RData')
load('./data/emg-models-app03/svm.polynomial.RData')
load('./data/emg-models-app03/svm.linear.RData')
load('./data/emg-models-app03/knn.model.RData')

model.perf <- rbind(
  data.frame(
    cv = rf.models$fold10,
    test.20 = rf.models$test,
    model = 'RF'),
  data.frame(
    cv = knn.model$fold10,
    test.20 = knn.model$test,
    model = 'knn'),
  data.frame(
    cv = svm.radial$fold10,
    test.20 = svm.radial$test,
    model = 'SVM radial'),
  data.frame(
    cv = svm.polynomial$fold10,
    test.20 = svm.polynomial$test,
    model = 'SVM polynomial'),
  data.frame(
    cv = svm.linear$fold10,
    test.20 = svm.linear$test,
    model = 'SVM linear'),
  data.frame(
    cv = svm.sigmoid$fold10,
    test.20 = svm.sigmoid$test,
    model = 'SVM sigmoid')
)

library(ggplot2)
model.m <- reshape::melt(model.perf)
model.m$value <- model.m$value*100
colnames(model.m)[2] <- 'robustness'
p <- ggplot(model.m,aes(x=model,y=value)) +
  geom_boxplot(aes(fill=robustness)) +
  xlab('Model approach') + ylab('Accuracy(%)') +
  ggtitle('Diff ML model acc stats on EOG/EMG time distribution \napproach 03') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

dir.create("./figures")
png('./figures/010_EMG_stats-ClassificationResults.png',
    width = 1000,height = 1200,res = 200)
print(p)
dev.off()
  