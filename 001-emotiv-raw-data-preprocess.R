xl_files<-list.files(path = "./1_RawEmotivData/",pattern=".csv")

files_sans_headers<-lapply(xl_files,function(x){read.csv(file = paste("./1_RawEmotivData/",x,sep=""),header=FALSE,sep = ",",stringsAsFactors = FALSE,skip = 1)})  

first_lines<-unlist(lapply(xl_files,function(x){return(readLines(con = paste("./1_RawEmotivData/",x,sep=""),n = 1))}))

header_stripper<-function(first_lines){
raw_line<-unlist(strsplit(first_lines,split = ","))
process_line<-grep(pattern = "labels",x = raw_line)
mat_header<-strsplit(x = raw_line[process_line],split = " ")[[1]][-1]
return(mat_header)
}

header_list<-lapply(first_lines,header_stripper)
unique(unlist(lapply(header_list,length)))

for(i in 1:length(files_sans_headers))
{
  colnames(files_sans_headers[[i]])<-header_list[[i]]
}

dir.create("2_FormattedEmotivData")
for(i in 1:length(files_sans_headers))
{
  write.table(x = files_sans_headers[[i]],file = paste("./2_FormattedEmotivData/",xl_files[i],"_formatted_",sep=""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
}

EEG_raw_data_emotiv_labeled<-files_sans_headers
rm(files_sans_headers)
subject_names<-unlist(lapply(xl_files,function(x){strsplit(x,split = "_")[[1]][1]}))
names(EEG_raw_data_emotiv_labeled)<-subject_names
dir.create("RData")
save(EEG_raw_data_emotiv_labeled,file = "./RData/EEG_rawdata_from_Emotiv.rda",compress = TRUE)

