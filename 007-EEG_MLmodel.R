rm(list=ls())
#setwd("~/Documents/Sandbox/ceeri-igib-ai-project/")
load('./RData/006_FinalEEGdataForMLRelPower_32Hz.rda')
source('./func-room.R')

dir.create("./data")
dir.create('./data/eeg-models-app03',showWarnings = TRUE)

#colnames(Final_EEGData_forML)
colnames(Final_EEGData_forML)[1]<-c("subjectName")
class(Final_EEGData_forML)

eeg.df_02<-cbind(Final_EEGData_forML[,-c(1)],Final_EEGData_forML[,c(1)])
#head(eeg.df_02[,71])
colnames(eeg.df_02)
colnames(eeg.df_02)[ncol(eeg.df_02)]<-"subjectName"
colnames(eeg.df_02)<-gsub(pattern = " ",x = colnames(eeg.df_02),replacement = "_")


rf.models <- multiRF(data = eeg.df_02[,-c(72)])
save(
  rf.models,
  file = './data/eeg-models-app03/rf.models.RData'
)

svm.radial <- multiSVM(
  data = eeg.df_02[,-c(72)],
  kernel.type = 'radial'
)

save(
  svm.radial,
  file = './data/eeg-models-app03/svm.radial.RData'
)

svm.polynomial <- multiSVM(
  data = eeg.df_02[,-72],
  kernel.type = 'polynomial'
)
save(
  svm.polynomial,
  file = './data/eeg-models-app03/svm.polynomial.RData'
)

svm.sigmoid <- multiSVM(
  data = eeg.df_02[,-72],
  kernel.type = 'sigmoid'
)
save(
  svm.sigmoid,
  file = './data/eeg-models-app03/svm.sigmoid.RData'
)
svm.linear <- multiSVM(
  data = eeg.df_02[,-72],
  kernel.type = 'linear'
)
save(
  svm.linear,
  file = './data/eeg-models-app03/svm.linear.RData'
)
knn.model <- multiKNN(
  data = eeg.df_02[,-72]
)
save(
  knn.model,
  file = './data/eeg-models-app03/knn.model.RData'
)

rm(list=ls())
load('./data/eeg-models-app03/rf.models.RData')
load('./data/eeg-models-app03/svm.radial.RData')
load('./data/eeg-models-app03/svm.sigmoid.RData')
load('./data/eeg-models-app03/svm.polynomial.RData')
load('./data/eeg-models-app03/svm.linear.RData')
load('./data/eeg-models-app03/knn.model.RData')

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
  ggtitle('Diff ML model accuracy stats on EEG \napproach 03, method median of \n 180_Window_10_SlideInterval_windowed FFT-PSD-integration') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

dir.create("./figures")
png('./figures/007_StdProcess-ClassificationResults_64hz_windowedPSDmedian.png',
    width = 1000,height = 1200,res = 200)
print(p)
dev.off()
    