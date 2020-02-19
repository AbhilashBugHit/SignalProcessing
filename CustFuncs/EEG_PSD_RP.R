EEG_PSD_RP<-function(x,Fs,lower,upper,relative_power=TRUE)
  {
  require(eegkit)
  xfft<-eegkit::eegfft(x, Fs, lower, upper)
  nc <- ncol(xfft$strength)
  xfft$strength <- xfft$strength^2
  xfft_b<-xfft
  relPow_xfft<-xfft$strength/rowSums(xfft$strength)
  xfft$strength<-relPow_xfft
  imagebar(xfft$frequency, 1:nc, xfft$strength,xlab="Frequency",ylab="Channel",zlab="Relative Power")
 if(relative_power==TRUE){invisible(xfft)}
if(relative_power==FALSE){invisible(xfft_b)}
}

