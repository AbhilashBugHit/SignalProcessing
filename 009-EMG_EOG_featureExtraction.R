
# Step 1: segment the data in normal and fatigue class
# step 2: run the SVM classifier on the data
# step 3: modify the ML algorithm to use variance, kurtosis and skew?
# step 3: record the time stamps for each eog event (sample index)
# 
rm(list=ls())

#--------------------------------------------------------------------------------------#
# Build an kNN classifier, classify instances of good and bad signals
# Remove the bad signals using 
# Calculate Variance-Kurtosis-Vpp features for SVM
ML_window_stats<-function(slide_window){
  require(e1071)
  varian<-var(slide_window)
  kurtos<-kurtosis(slide_window)
  #ampliPP<-diff(range(slide_window))
  skew<-skewness(slide_window)
  #ML_features<-c(varian,kurtos,ampliPP)
  ML_features<-c(varian,kurtos,skew)
  #names(ML_features)<-c("Variance","Kurtosis","AmplitudePtP")
  names(ML_features)<-c("Variance","Kurtosis","Skewness")
  return(ML_features)
}
load("./InputRData/badEEGtrainingData.rda")
load("./InputRData/goodEEGtrainingData.rda")
GoodEEG_Features<-t(apply(goodEEGmatrix,1,ML_window_stats))
BadEEG_Features<-t(apply(badEEGmatrix,1,ML_window_stats))

#par(mfrow=c(2,3))
#apply(GoodEEG_Features,2,hist)
#apply(BadEEG_Features,2,hist)
#head(GoodEEG_Features)
#dev.off()

class_label<-as.factor(c(rep("Good",nrow(GoodEEG_Features)),rep("Bad",nrow(BadEEG_Features))))
AllEEG_Features<-rbind(GoodEEG_Features,BadEEG_Features)
AllEEG_df<-data.frame(AllEEG_Features)
#head(AllEEG_df)
#apply(AllEEG_df,2,class)
AllEEG_df$class<-class_label
#colnames(AllEEG_df)<-c("Variance","Kurtosis","Vpp","class")
colnames(AllEEG_df)<-c("Variance","Kurtosis","Skewness","class")
source("~/Desktop/CompleteWorkflow/func-room.R")

# Run classifiers on the variables
#multiKNN(AllEEG_df)
#multiRF(AllEEG_df)
#multiSVM(AllEEG_df)

EEG_KNN <- function(data){
  library(caret)
  knn.model <- knn3(class~., data=data)
  return(knn.model)
}

# K nearest neighbours model is trained on training data to recognize the EOG artifact from
# new signals.
knn.model.eeg<-EEG_KNN(AllEEG_df) 
# So now a model that can call out a good or bad EEG on the basis of statistical features is 
# available
save(knn.model.eeg,file="./RData/EOGPredictor_KNN.rda",compress=TRUE)
save(ML_window_stats,file="./RData/MLwindowStats.rda",compress=TRUE)
rm(list=ls())
#--------------------------------------------------------------------------------------#
#--------------------- CLASS SEGMENTATION into NORMAL and FATIGUE ---------------------#

all_files<-list.files(path="./5_ButterworthFiltered/",pattern = "BWFiltered")
subjectNames<-unlist(lapply(strsplit(all_files,split = "_"),function(x){x[1]}))
dir.create(path = "./9_ClassSegmentation")
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
    write.table(x = norm_matrix,file = paste("./9_ClassSegmentation/",subjectNames[i],"_Norm",".csv",sep=""), sep=",", quote = FALSE,row.names = TRUE,col.names = TRUE)
    write.table(x = fatig_matrix,file = paste("./9_ClassSegmentation/",subjectNames[i],"_Fatig",".csv",sep=""), sep=",", quote = FALSE,row.names = TRUE,col.names = TRUE)  
  }
}
 #--------------------------------------------------------------------------------------#

#
# Get the SVM classifier that you created for identifying EOG/EMG.
# Count the number of EOG EMG events and build a summary of the frequency of events versus time
# For normal and fatigue state, count the distribution of events in 10 minutes each.
# Keep track of time between EOG events as well, in normal and fatigue epochs
# estimate the parameters of the total frequency of distribution as mean, variance, kurtosis, skew.
# treat each minute of the 10 minutes as a bin and look at the number of EOG events in normal and 
# fatigued.

# Pseudocode
# load the normal and fatigue segments from the files
# they will be 128*60*10 samples long
# iterate over them with the sliding window

Norm_segmentFiles<-list.files(path="./9_ClassSegmentation/",pattern = "Norm",full.names = TRUE)
Fatig_segmentFiles<-list.files(path="./9_ClassSegmentation/",pattern = "Fatig",full.names = TRUE)

subjectNames<-gsub(pattern = "_Norm.csv",replacement = "",list.files(path="./9_ClassSegmentation/",pattern = "Norm"))

NormSegments<-lapply(Norm_segmentFiles,function(x){
  read.table(file=x, header=TRUE,sep=",",stringsAsFactors=FALSE)
  })
names(NormSegments)<-subjectNames

FatigSegments<-lapply(Fatig_segmentFiles,function(x){
  read.table(file=x, header=TRUE,sep=",",stringsAsFactors=FALSE)
})
names(FatigSegments)<-subjectNames

save(NormSegments, file = "./InputRData/NormSegments.rda",compress = TRUE)
save(FatigSegments, file = "./InputRData/FatigSegments.rda",compress = TRUE)

rm(list=ls())

load("./RData/EOGPredictor_KNN.rda")
load("./RData/MLwindowStats.rda")
load("./InputRData/NormSegments.rda")
load("./InputRData/FatigSegments.rda")


# loop through each individual, then each spectrum, then within the time series
# then slide windows within the time series and calculate the 3 statistics
# pass the statistics and classify segments using trained SVM
# for Output classes from SVM that are bad, retrieve indices.
# Use indices to retrieve signal from spectrum in a list
# pass the list to autoencoder filter, obtain filtered signal as output
# replace output into input time series at the correct indices.

Fs=128;BlinkSec=0.45;HalfWidthBSec=BlinkSec/2;

pb<-txtProgressBar(min = 0,max = length(NormSegments),style=3)

Subjects_NormChannelList<-vector(mode="list",length = length(NormSegments))
Subjects_FatigChannelList<-vector(mode="list",length = length(NormSegments))
for(sub_ix in 1:length(NormSegments))
{
  setTxtProgressBar(pb,value = sub_ix)
  BadSegments_NormIndexList<-vector(mode="list",length = ncol(NormSegments[[sub_ix]]))
  BadSegments_FatigIndexList<-vector(mode="list",length = ncol(NormSegments[[sub_ix]]))
  
    for(ch_ix in 1:ncol(NormSegments[[sub_ix]]))
      {
      # defining the sliding window interval and sliding window overlap
      SWInterval<-round(BlinkSec*Fs)
      SWOverlap<-round(HalfWidthBSec*Fs)
    
      #rm(BlinkSec,Fs,HalfWidthBSec)
      # generating the window index sequence
      start_win_ix<-seq(from = 1,to = nrow(NormSegments[[sub_ix]]),by=SWOverlap)
      end_win_ix<-start_win_ix+SWInterval
      #create a matrix of 2 rows where each row contains the start and end indices
      # describes the co-ordinates of the sliding window
      index_mat<-rbind(start_win_ix,end_win_ix)
      rm(start_win_ix,end_win_ix,SWOverlap,SWInterval)
    
      ml_stats_norm=matrix(nrow = 1,ncol = 3)
      ml_stats_fatig=matrix(nrow = 1,ncol = 3)
      colnames(ml_stats_norm)<-c("Variance","Kurtosis","Skewness")
      colnames(ml_stats_fatig)<-c("Variance","Kurtosis","Skewness")
    
      # use the indices to create sliding windows and cut the window
      # out into variable data and calculate the 
      # three statistics for each of them, store the value of the three statistics
      # as an rbind operation in the ml_stats matrix
        for(slidWin in 1:ncol(index_mat))
          {
          vert_bars<-c(index_mat[1,slidWin],index_mat[2,slidWin])
          data_norm<-NormSegments[[sub_ix]][c(vert_bars[1]:vert_bars[2]),ch_ix]
          data_fatig<-FatigSegments[[sub_ix]][c(vert_bars[1]:vert_bars[2]),ch_ix]
          ml_stats_norm<-rbind(ml_stats_norm,ML_window_stats(slide_window = data_norm))
          ml_stats_fatig<-rbind(ml_stats_fatig,ML_window_stats(slide_window = data_fatig))
          }
    
      # So ml_stats_* now contains statistics for the sliding window, the co-ordinate system
      # will now refer to the indices of the ML_window statistics. It will initially be 
      # dimensioned the same way as index_mat.
    
      # ml_stats contains 3 feature stats for all the windows and is of same ncol as index_mat
      # we check out which row indices in ml_stats contain NA's
      # we want to track which indices in index_mat got a "bad" tag
      # so we keep a separate index vector, pass it to ml_stats after removing NA ix
      # So now we get a bad tag on indices only present in ml_ix
      # we retrieve the indices of ml_ix that got a bad tag
    
      ml_stats_norm<-data.frame(ml_stats_norm[-1,])
      ml_stats_fatig<-data.frame(ml_stats_fatig[-1,])
    
      #------------Removing NA indices --------------------#
      NA_ix<-unique(unlist(apply(ml_stats_norm,2,function(x){which(is.na(x))})))
      ml_ix_norm<-1:nrow(ml_stats_norm)
      ml_ix_norm<-ml_ix_norm[-c(NA_ix)]
    
      NA_ix<-unique(unlist(apply(ml_stats_fatig,2,function(x){which(is.na(x))})))
      ml_ix_fatig<-1:nrow(ml_stats_fatig)
      ml_ix_fatig<-ml_ix_fatig[-c(NA_ix)]
    
      # To preserve the index ordering, we store the original index and then remove
      # the NA's from that index. We use this indexing to present the data to the 
      # ML stats (which does not tolerate NA's).
      # But now the indexing that we retrieve from the ML classification will
      # correspond to the indexing scheme of the filtered list
      # how to remap it back to the original list?
    
    #-----------------------------------------------------#
      segment_prediction_norm<-predict(object=knn.model.eeg, newdata = ml_stats_norm[ml_ix_norm,],type = 'class')
      segment_prediction_fatig<-predict(object=knn.model.eeg, newdata = ml_stats_fatig[ml_ix_fatig,],type = 'class')
    
      bad_segs_norm<-ml_ix_norm[which(segment_prediction_norm=="Bad")] # this tracks the bad segments indices correctly.
      bad_segs_fatig<-ml_ix_fatig[which(segment_prediction_fatig=="Bad")]
      # So bad_segs_norm/fatig specifically refers to the original indexing of ML
      # window stats sans NA's and now tagged with bad segments and therefore
      # can be remapped to index_mat to obtain original sliding window indices    
      
      #print(table(segment_prediction_norm))
      #print(table(segment_prediction_fatig))
      # index_mat[,bad_segs_norm]
      # index_mat[,bad_segs_fatig]
  
      BadSegments_FatigIndexList[[ch_ix]]<-index_mat[,bad_segs_fatig]  
      BadSegments_NormIndexList[[ch_ix]]<-index_mat[,bad_segs_norm]
      
      }# end of for loop for channel index iteration
    
    Subjects_NormChannelList[[sub_ix]]<-BadSegments_NormIndexList
    Subjects_FatigChannelList[[sub_ix]]<-BadSegments_FatigIndexList
}  
 
save(Subjects_NormChannelList,file = "./RData/EOGMEG_Events_NormalSubject.rda",compress=TRUE)
save(Subjects_FatigChannelList,file = "./RData/EOGMEG_Events_FatigueSubject.rda",compress=TRUE)

# So now we have the EOG events tagged across sliding windows for 25 subjects in their 
# Normal segment and Fatigue segment across 14 channels.

# Can we look at some statistics regarding the how to record the EOG frequency?
# It seems that 10 second windows a lot of EOG events depending on people
# Can you 

# The start point of each window is spaced 29 samples apart because we overlap by
# 0.225 seconds on what is essentially a 0.45 second event

# By looking at the diff of the rows of the EOG events stored and marking all 
# events above 29, we can find the frequency of unique EOG events spaced in time
# events that are spaced at 29 are sequential frames of the same EOG event

# Diffs that are >29 are the timespaces between two unique EOG events
# There was talk of this spacing being informative of normal and fatigue status
# I will get a distribution of this timing for each subject across different channels
# my feature table would be the Subject on the Rows, CH_Moments(timing_distribution)

#------------- ADDING SUBJECT NAMES AND CHANNEL NAMES TO NORMAL AND FATIGUE SEGMENTS------------#
# Extracting sample and channel names which have gone missing in the workflow.  
sampleFile<-list.files(path = "./5_ButterworthFiltered/",pattern = "BWFiltered",full.names = TRUE)
ch_names<-unlist(strsplit(readLines(con = sampleFile[1],n = 1),split = ","))
prefix<-gsub(pattern = "./5_ButterworthFiltered//",replacement = "",sampleFile)
subjectNames<-gsub(prefix,pattern = "_BWFiltered",replacement = "")

names(Subjects_FatigChannelList)<-subjectNames
names(Subjects_NormChannelList)<-subjectNames

Subjects_NormChannelList_named<-lapply(Subjects_NormChannelList,function(sub){
  names(sub)<-ch_names
  return(sub)})

Subjects_FatigChannelList_named<-lapply(Subjects_FatigChannelList,function(sub){
  names(sub)<-ch_names
  return(sub)})
#--------------------------------------------------------------------------------------------#
#============ Function to take EOG x-coordinates and convert to feature table==================#
EEG_matrix_to_FT<-function(Subjects_ChannelList_named){
#---------- Calculating intervals of EOG events---------------------------------#
Sub_EOG_timeDistrib<-lapply(Subjects_ChannelList_named,function(sub){
  lapply(sub,function(ch){
  diff(t(ch))
  })
})

#-----------Subsetting single column for intervals ------------------------------#
Sub_EOG_timeDist2<-lapply(Sub_EOG_timeDistrib,function(sub)
  {
  lapply(sub,function(ch)
    {
    if(length(dim(ch))!=0)
      {
      return(ch[,1])
      }
    })
  })
#-----------------------------------------------------------------------------

EOG_Interval_Stats<-function(EOG_Event_Time_Distribution_Vec){
  EETD<-EOG_Event_Time_Distribution_Vec
  mu<-mean(EETD,na.rm=TRUE)
  sigsq<-var(EETD,na.rm = TRUE)
  kurty<-e1071::kurtosis(EETD,na.rm = TRUE)
  skewy<-e1071::skewness(EETD,na.rm = TRUE)
  EETDstats<-c(mu,sigsq,kurty,skewy)
  names(EETDstats)<-c("Mean","Variance","Kurtosis","Skewness")
  return(EETDstats)
  }

Sub_EOG_timeDist_stats<-lapply(Sub_EOG_timeDist2,function(sub)
  {
  lapply(sub,function(ch)
    {
    if(length(ch)!=0)
      {
      EOG_Interval_Stats(ch)
      }
    })
  }
  )

#NormSegments[[1]][[1]]
#Sub_Norm_EOG_timeDist2[[1]][[1]]

Sub_EOG_timeDist_stats_Ch2df<-lapply(Sub_EOG_timeDist_stats,function(sub){plyr::ldply(sub)})

Sub_EOG_timeDist_stats_Ch2dfnamed<-lapply(Sub_EOG_timeDist_stats_Ch2df,function(x){names(x)[1]<-"Channel"
return(x)
})

Sub_Ch_Stats<-plyr::ldply(Sub_EOG_timeDist_stats_Ch2dfnamed)

melted_Sub_Ch_Stats<-melt(Sub_Ch_Stats)
id_Chstats_value<-tidyr::unite(data=melted_Sub_Ch_Stats,col="ChannelStats",Channel,variable)
feature_table<-tidyr::spread(data=id_Chstats_value,key = ChannelStats,value = value)

return(feature_table)
}
#==============================================================================================#

# Random Forest based imputation of missing data
RFImpute<-function(df){
impute_df<-missForest::missForest(df[,-1])
df2<-cbind(Norm_df[,1],impute_df[[1]])
names(df2)[1]<-"SubjectName"
return(df2)
}

#-------------------------OUTPUT testing section -------------------------------------------#
# predicted scores
Norm_df<-EEG_matrix_to_FT(Subjects_NormChannelList_named)
Fatig_df<-EEG_matrix_to_FT(Subjects_FatigChannelList_named)
Norm_Imp_df<-RFImpute(Norm_df)
Fatig_Imp_df<-RFImpute(Fatig_df)
Norm_Imp_df$Class<-c("Normal")
Fatig_Imp_df$Class<-c("Fatigue")
ML_df<-rbind(Norm_Imp_df,Fatig_Imp_df)

save(ML_df,file = "./RData/009_ML_df.rda",compress = TRUE)

#------------------------------------------------------------------------------------------#