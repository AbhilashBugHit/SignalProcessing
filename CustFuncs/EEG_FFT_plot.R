EEG_FFT_plot<-function(input1,input2,Fsamp=128,Lfreq=0.5,Hfreq=35,plot_type="l"){
EEG_FFT1<-eegkit::eegfft(x = input1, Fs =Fsamp, lower = Lfreq, upper = Hfreq)
EEG_FFT2<-eegkit::eegfft(x = input2, Fs =Fsamp, lower = Lfreq, upper = Hfreq)
plot(EEG_FFT1$frequency,y = EEG_FFT1$strength,type=plot_type,col="red",main="FFT comparison",lwd=0.7,xlab="Frequency",ylab="Strength")
lines(EEG_FFT2$frequency,y =EEG_FFT2$strength,type=plot_type,lwd=0.7,col="blue")
abline(v=c(0.5,4,8,12,30),lty=3)
text_label_Ycoord<- max(c(EEG_FFT1$strength,EEG_FFT2$strength))*0.60
band_label<-c("Delta","Theta","Alpha","Beta","Gamma")
text(x = c(0.5,4,8,12,30),y = c(rep(text_label_Ycoord,length(band_label))),labels = band_label,pos = 4)
legend(x = 0.7*Hfreq,y = text_label_Ycoord,legend = c("Input 1","Input 2"),col = c("red","blue"),lty=c(1,1),bty = "n")
EEGFFT_List<-list()
EEGFFT_List$fft1<-EEG_FFT1
EEGFFT_List$fft2<-EEG_FFT2
invisible(EEGFFT_List)
}