rm(list=ls())

library(reshape)
library(magrittr)
library(dplyr)
library(plyr)
library(tidyr)

fileNames<-list.files(path = "./6_EEGPSD/Normal/")
subjectNames<-unlist(lapply(strsplit(x = fileNames,split = "_"),function(x){x[[1]][1]}))

norm_files<-list.files(path = "./6_EEGPSD/Normal/",full.names = TRUE)
fatig_files<-list.files(path = "./6_EEGPSD/Fatigue/",full.names = TRUE)


norm_psd_list<-lapply(X = norm_files,FUN = function(filename){read.table(file = filename,header = TRUE,stringsAsFactors = FALSE,sep = ",")})
fatig_psd_list<-lapply(X = fatig_files,FUN = function(filename){read.table(file = filename,header = TRUE,stringsAsFactors = FALSE,sep = ",")})

names(norm_psd_list)<-subjectNames
names(fatig_psd_list)<-subjectNames

npl_test<-lapply(norm_psd_list,function(x){
  y<-cbind(rownames(x),x);
  rownames(y)<-NULL;
  colnames(y)[1]<-"PSD"
  return(y)})

fpl_test<-lapply(fatig_psd_list,function(x){
  y<-cbind(rownames(x),x);
  rownames(y)<-NULL;
  colnames(y)[1]<-"PSD"
  return(y)})

norm_df<-plyr::ldply(npl_test)
fatig_df<-plyr::ldply(fpl_test)  


df_to_ft<-function(norm_df){
melted_norm_df<-reshape2::melt(data = norm_df)
melted_norm_df$PSD<-as.character(melted_norm_df$PSD)
melted_norm_df$variable<-as.character(melted_norm_df$variable)
powband<-apply(melted_norm_df[,c(2,3)],1,function(x){paste(x,sep="_",collapse = "_")})
melted_norm_df$variable<-powband
melted_norm_df<-melted_norm_df[,-c(2)]
norm_feature<-tidyr::spread(data=melted_norm_df,key = variable,value = value)
return(norm_feature)
}

norm_feature<-df_to_ft(norm_df = norm_df)
fatig_feature<-df_to_ft(norm_df = fatig_df)
norm_feature$class<-as.factor("Normal")
fatig_feature$class<-as.factor("Fatigue")
full_ft<-rbind(norm_feature,fatig_feature)
colnames(full_ft)[1]<-c("SubjectName")
Final_EEGData_forML<-full_ft

save(Final_EEGData_forML,file = "./RData/006_FinalEEGdataForMLRelPower_32Hz.rda",compress = TRUE)
#-----------------------------------------------------#