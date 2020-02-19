
# Power spectral decomposition for EEG signal to build feature table
rm(list=ls())
source("./CustFuncs/EEG_5Band_PSD.R")
filenames<-list.files(path = "./012_TimeSegmentation/",full.names = TRUE)
labeled_EEG_time_segments<-lapply(filenames,function(filex){subx<-read.delim(file = filex,header=TRUE,sep = ",",stringsAsFactors = FALSE);
sub_eegpsd<-EEG_5Band_PSD(multiChEEG = subx,Fs = 128,lower = 0.5,upper = 63.5)
return(sub_eegpsd)})
names(labeled_EEG_time_segments)<-gsub(x = filenames, pattern="./012_TimeSegmentation//",replacement = "")

dir.create(path = "./013_EEGFeatureSegmentPSD")

for(i in 1:length(labeled_EEG_time_segments))
{
  outfile_table<-labeled_EEG_time_segments[[i]]
  fullpath<- paste("./013_EEGFeatureSegmentPSD/",names(labeled_EEG_time_segments)[i],"_PSD",sep="")
  write.table(x = outfile_table,file =fullpath ,sep=",",row.names = TRUE,col.names = TRUE,quote = FALSE)
}
    