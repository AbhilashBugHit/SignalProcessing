rm(list=ls())

# Version number indicates 1 (direct EEG PSD with squaring and adding up of fourier co-efficients 
# of whole 10 minute epoch)

# Version number 2 indicates sliding window of EEG PSD of 10 seconds in a 3 minute window
# with simpson integration of AUC for fourier co-efficients corresponding to each EEG band
# to obtain the total power within the EEG band.

version=1;

if(version==2){source("./CustFuncs/EEG_6Band_PSD.R")}
if(version==1){source("./CustFuncs/EEG_5Band_PSD.R")}
dir.create("./RData")
dir.create(path = "./6_EEGPSD")
dir.create(path = "./6_EEGPSD/Normal")
dir.create(path = "./6_EEGPSD/Fatigue")

all_files<-list.files(path="./5_ButterworthFiltered/",pattern = "BWFiltered")
subjectNames<-unlist(lapply(strsplit(all_files,split = "_"),function(x){x[1]}))


# Skip first 5 minutes for headset acclimatization and electrode setting time
# Take first 10 minutes EEG as normal period
# For the full recording, take the last 10 minutes as fatigued period.
# Sampling rate is 128Hz, so 10 mins is 128*10*60

norm_psd_list<-vector(mode="list")
fatig_psd_list<-vector(mode="list")


pb<-txtProgressBar(min=1,max=length(all_files),style = 3)
for(i in 1:length(all_files)){
  setTxtProgressBar(pb,i)
  df<-read.table(file=paste("./5_ButterworthFiltered/",all_files[i],sep=""), header=TRUE,sep=",",stringsAsFactors=FALSE)
  skip_ix = 128*60*5
  norm_start_ix = 1+skip_ix
  norm_interval = (128*60*10)
  norm_end_ix = norm_start_ix + norm_interval 
  norm_ix <- c(norm_start_ix:norm_end_ix)
  fatig_end_ix = nrow(df)
  fatig_interval <- (128*60*10)
  fatig_start_ix = nrow(df) - fatig_interval
  fatig_ix <- c(fatig_start_ix:fatig_end_ix)
  norm_matrix<-df[norm_ix,]
  fatig_matrix<-df[fatig_ix,]
  
  total_epochs_reqd<- skip_ix+norm_interval+fatig_interval
  epochs_avlbl<-nrow(df)
  if(epochs_avlbl>=total_epochs_reqd)
  {
  if(version==1){norm_psd<-EEG_5Band_PSD(multiChEEG = norm_matrix,Fs = 128,lower = 0.5,upper = 63)}
  if(version==2){norm_psd<-EEG_slidWin_PSD(multiChEEG = norm_matrix,Fs = 128,window_size = 180,slide_interval = 10)}
  if(version==1){fatig_psd<-EEG_5Band_PSD(multiChEEG = fatig_matrix,Fs = 128,lower = 0.5,upper = 63)}
  if(version==2){fatig_psd<-EEG_slidWin_PSD(multiChEEG = fatig_matrix,Fs = 128,window_size = 180,slide_interval = 10)}
  write.table(x = norm_psd,file = paste("./6_EEGPSD/Normal/",subjectNames[i],"_PSD",".csv",sep=""), sep=",", quote = FALSE,row.names = TRUE,col.names = TRUE)
  write.table(x = fatig_psd,file = paste("./6_EEGPSD/Fatigue/",subjectNames[i],"_PSD",".csv",sep=""), sep=",", quote = FALSE,row.names = TRUE,col.names = TRUE)
  }
}

