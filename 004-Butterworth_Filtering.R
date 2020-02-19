rm(list=ls())
require("eegkit")
csv_files_RawData<-list.files(path="./4_EEG_GyroRemoved/",pattern = "gyroRemoved")
subjectName<-unlist(lapply(lapply(csv_files_RawData,function(x){strsplit(x,split="_")[[1]]}),function(y){y[1]}))

dir.create("5_ButterworthFiltered",showWarnings = TRUE)

pb<-txtProgressBar(min=1,max = length(csv_files_RawData),initial = 1,style = 3)
for(i in 1:length(csv_files_RawData))
{
  setTxtProgressBar(pb,i)
  if(length(readLines(con=paste("./4_EEG_GyroRemoved/",csv_files_RawData[i],sep="")))>10)
  {
    rawdata_14ch<-read.delim(file = paste("./4_EEG_GyroRemoved/",csv_files_RawData[i],sep=""), header = TRUE,sep = ",",stringsAsFactors = FALSE)
    dm<-data.matrix(rawdata_14ch)
    filtered_14ch<-eegkit::eegfilter(x = dm,Fs = 128,lower = 0.5,upper = 63.5,method = "butter",order = 4, plot = FALSE)
    write.table(filtered_14ch,file = paste("./5_ButterworthFiltered/",subjectName[i],"_BWFiltered",sep=""),row.names = FALSE,quote = FALSE,sep=",")
  }
}
