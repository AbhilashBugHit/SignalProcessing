# Takes the multichannel EEG data  for a single user along with channels/columns 18 to 31 to create a stacked bar-plot that shows the percentage distribution of each signal quality metric for every channel
# columns with CQ prefixed to the channel name are the Channel Quality metrics.  

#---------------------------------------------------------------------------------------------#
SignalQuality_Plot<-function(SingleUser_data,UserName){
CQ_score_table<-apply(SingleUser_data[,18:31],2,function(x){table(x)})
#sum(CQ_score_table[[1]])
if(is.vector(CQ_score_table)==TRUE && is.list(CQ_score_table)==FALSE){
  CQ_score_table<-matrix(CQ_score_table,ncol=1)
  label<-unique(unlist(unique(as.vector(SingleUser_data[,18:31]))))
  
  if(length(label)==1){
    colnames(CQ_score_table)<-label 
    }
  message("Only one value of CQ was found; returned a vector")
  }

if(is.list(CQ_score_table)==TRUE){
  CQ_came_as_matrix=FALSE
  CQ_score_norm<-lapply(CQ_score_table,function(x){(x/sum(x,na.rm = TRUE))*100})  
  CQ_range<-as.character(c(0:4))
  
  for(i in 1:length(CQ_score_norm)){
    current_scores<-names(CQ_score_norm[[i]])
    missing<-setdiff(CQ_range,current_scores)
    
    if(length(missing)>0)
    {
      CQ_score_norm[[i]]<-c(CQ_score_norm[[i]],rep(0,length(missing)))
      names(CQ_score_norm[[i]])[which(CQ_score_norm[[i]]==0)]<-missing
      CQ_score_norm[[i]]<-CQ_score_norm[[i]][order(names(CQ_score_norm[[i]]))]
    }
    } # end of for loop
  message("Non uniform: all CQ scores were not observed, added scores not observed at 0%")
  CQ_df<-plyr::ldply(CQ_score_norm)
  # CQ_df$.id<-colnames(SingleUser_data[,18:31])
  }

if(is.matrix(CQ_score_table)==TRUE && is.list(CQ_score_table)==FALSE){
  if(ncol(CQ_score_table)>1){
  message("All CQ Scores were observed,returned matrix")
  CQ_came_as_matrix=TRUE
  CQ_score_norm<-apply(CQ_score_table,2,function(x){(x/sum(x))*100})
  CQ_df<-data.frame(t(CQ_score_norm))
  CQ_df$.id<-colnames(SingleUser_data[,18:31])
  }
  
  if(ncol(CQ_score_table)==1)
    {
    message("Only 1 type of CQ score")
    CQ_score_norm<-apply(CQ_score_table,1,function(x){(x/sum(x))*100})
    CQ_df<-data.frame(CQ_score_norm)
    CQ_df$.id<-colnames(SingleUser_data[,18:31])
    }
  #message()
  }

  
#if(is.matrix(CQ_score_norm)!=TRUE){
#  CQ_df<-plyr::ldply(CQ_score_norm)
 # CQ_df$.id<-colnames(SingleUser_data[,18:31])
  
 # if(is.vector(CQ_score_norm)==TRUE && is.list(CQ_score_norm)==FALSE)
#      {
#      CQ_score_norm2<-matrix(CQ_score_norm,ncol=1)
#      rownames(CQ_score_norm2)<-colnames(SingleUser_data[,18:31])
 #     colnames(CQ_score_norm2)<-as.character(label)
    #CQ_score_norm<-t(CQ_score_norm)
#      CQ_df<-data.frame(CQ_score_norm2)
#      CQ_df$.id<-colnames(SingleUser_data[,18:31])
#      }
#    }

#CQ_score_norm

#CQ_df
CQ_df_melted<-reshape2::melt(CQ_df)
#UserName=SubjectNames[i]
ggp<-ggplot(CQ_df_melted, aes(fill=variable, y=value, x=.id)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle(label = UserName)+
  xlab("Electrode")+
  ylab("Percentage")+
  labs(fill="Signal Quality")

  
print(ggp)
}