rm(list=ls())

library(reshape)
library(magrittr)
library(dplyr)
library(plyr)
library(tidyr)

fileNames<-list.files(path = "./013_EEGFeatureSegmentPSD/",full.names = TRUE)
fname2<-gsub(x = fileNames,pattern = "./013_EEGFeatureSegmentPSD//",replacement = "")
subjectNames<-unlist(lapply(strsplit(x = fname2,split = "_PSD"),function(x){x[[1]][1]}))

psd_list<-lapply(X = fileNames,FUN = function(filename){read.table(file = filename,header = TRUE,stringsAsFactors = FALSE,sep = ",")})
#fatig_psd_list<-lapply(X = fatig_files,FUN = function(filename){read.table(file = filename,header = TRUE,stringsAsFactors = FALSE,sep = ",")})

names(psd_list)<-subjectNames
#names(fatig_psd_list)<-subjectNames

npl_test<-lapply(psd_list,function(x){
  y<-cbind(rownames(x),x);
  rownames(y)<-NULL;
  colnames(y)[1]<-"PSD"
  return(y)})


PSD_df<-plyr::ldply(npl_test)


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

PSD_feature<-df_to_ft(norm_df = PSD_df)
class_label<-unlist(lapply(strsplit(x = PSD_feature$.id,split = "_"),function(x){x[3]}))
class_label<-gsub(pattern = " ",replacement = "_",x = class_label)
subjectNames<-unlist(lapply(strsplit(x = PSD_feature$.id,split = "_"),function(x){paste(x[1:2],sep="_",collapse = "_")}))
PSD_feature$subjectNames<-subjectNames
PSD_feature$class<-class_label
Final_EEGData_forML<-PSD_feature[,-1]


save(Final_EEGData_forML,file = "./RData/014_EpochLabelEEG_PSD_RelPower_64Hz.rda",compress = TRUE)
#-----------------------------------------------------#