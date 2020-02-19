rm(list=ls())
#setwd("~/Documents/Sandbox/ceeri-igib-ai-project/")
load("./RData/014_EpochLabelEEG_PSD_RelPower_64Hz.rda")
source('./func-room.R')

dir.create("./data")
dir.create('./data/eeg-nback-models-app01',showWarnings = TRUE)

#colnames(Final_EEGData_forML)
eeg.df_02<-Final_EEGData_forML
head(eeg.df_02[,71])
colnames(eeg.df_02)
colnames(eeg.df_02)[71]<-"subjectName"
colnames(eeg.df_02)<-gsub(pattern = " ",x = colnames(eeg.df_02),replacement = "_")

sub_ix<-71

rf.models <- multiRF(data = eeg.df_02[,-sub_ix])
save(
  rf.models,
  file = './data/eeg-nback-models-app01/rf.models.RData'
)

svm.radial <- multiSVM(
  data = eeg.df_02[,-sub_ix],
  kernel.type = 'radial'
)

save(
  svm.radial,
  file = './data/eeg-nback-models-app01/svm.radial.RData'
)

svm.polynomial <- multiSVM(
  data = eeg.df_02[,-sub_ix],
  kernel.type = 'polynomial'
)
save(
  svm.polynomial,
  file = './data/eeg-nback-models-app01/svm.polynomial.RData'
)

svm.sigmoid <- multiSVM(
  data = eeg.df_02[,-sub_ix],
  kernel.type = 'sigmoid'
)
save(
  svm.sigmoid,
  file = './data/eeg-nback-models-app01/svm.sigmoid.RData'
)
svm.linear <- multiSVM(
  data = eeg.df_02[,-sub_ix],
  kernel.type = 'linear'
)
save(
  svm.linear,
  file = './data/eeg-nback-models-app01/svm.linear.RData'
)
knn.model <- multiKNN(
  data = eeg.df_02[,-sub_ix]
)
save(
  knn.model,
  file = './data/eeg-nback-models-app01/knn.model.RData'
)

rm(list=ls())
load('./data/eeg-nback-models-app01/rf.models.RData')
load('./data/eeg-nback-models-app01/svm.radial.RData')
load('./data/eeg-nback-models-app01/svm.sigmoid.RData')
load('./data/eeg-nback-models-app01/svm.polynomial.RData')
load('./data/eeg-nback-models-app01/svm.linear.RData')
load('./data/eeg-nback-models-app01/knn.model.RData')

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
  ggtitle('Diff ML model accuracy stats on EEG \nNback 5 min label approach') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

dir.create("./figures")
png('./figures/015_StdProcess_5minNbackLabel_CfxnResults.png',
    width = 1000,height = 1200,res = 200)
print(p)
dev.off()
