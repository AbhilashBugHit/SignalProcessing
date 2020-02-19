rm(list=ls())
load("./RData/EEG_rawdata_from_Emotiv.rda")

EEG_raw_data_emotiv_labeled[[1]][1:6,1]

# 3:16 are the EEG channels and 18,19 are the gyroscope X and Y
data_column_indices<-c(1,3:16,18,19,26:40)

EEG_channel_data<-lapply(EEG_raw_data_emotiv_labeled,function(x){return(x[,data_column_indices])})

dir.create("3_EEG_GyroXY")

for(i in 1:length(EEG_channel_data))
{
  write.table(EEG_channel_data[[i]],file = paste("./3_EEG_GyroXY/",names(EEG_channel_data)[i],"_EEG_Gyro.csv",sep=""),quote = FALSE,row.names = FALSE,col.names = TRUE,sep=",")
}

