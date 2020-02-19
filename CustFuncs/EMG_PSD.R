# You retrieved the original signal now do a power spectral decomposition
# for each segment. pred_ind_mat contains indices that span an EOG event.
# So for each segment create a PSD type statistic that can be used for
# concatenating as a vector. So maybe look at HF/LF proportion
# So for an EMG event, there may be high-frequency components and low-frequency components
# Does the ratio between the components alter in mental fatigue? I am not sure about it.
# So we can try this with approach 3 first and then what i think is emg is 1 to 500Hz
# so with a sampling rate of 128 Hz, the maximum I can detect is 64
# so 1Hz-8Hz 9Hz-16Hz 17Hz-24Hz 25Hz-32Hz 33-63Hz
# Now we calculate these as features. Wait, we should combine indices that are consecutive
# and make them as one event. 
#the final output statistic has to be a 4 band vector containing relative power
# estimates for each channel (14) in the 5 minute segment.

# Now the EMG/EOG signals are really non-stationary. Their PSD estimates are best done
# not over the whole signal but in sliding windows that move across the sample
# Also maybe integrate over time to produce a better estimate of the power within the 
# EMG/EOG segment.

#EMG_segment<-ACEMG_UL[[2]]

EMG_PSD<-function(EMG_segment,
                  Fs=128,
                  lower=0.1,
                  upper=63,
                  B1=c(1,7.999),
                  B2=c(8,15.999),
                  B3=c(16,23.999),
                  B4=c(24,31.999),
                  B5=c(32,64),BandIntegrate=FALSE){
  # so 1Hz-8Hz 9Hz-16Hz 17Hz-24Hz 25Hz-32Hz 33-63Hz
  #Fs=128;lower=0.5;upper=63;B1=c(1,7.999);B2=c(8,15.999);B3=c(16,23.999); B4=c(24,31.999);B5=c(32,64);
  require(eegkit)
  xfft <- eegfft(EMG_segment, Fs, lower, upper)
  if(BandIntegrate==FALSE){  
    B1_ix<-which(xfft$frequency>=B1[1] & xfft$frequency<=B1[2])
    B2_ix<-which(xfft$frequency>=B2[1] & xfft$frequency<=B2[2])
    B3_ix<-which(xfft$frequency>=B3[1] & xfft$frequency<=B3[2])
    B4_ix<-which(xfft$frequency>=B4[1] & xfft$frequency<=B4[2])
    B5_ix<-which(xfft$frequency>=B5[1] & xfft$frequency<=B5[2])
    
    # The power of a Discrete Fourier Transform is the MEAN of 
    # the sum of squares of its Fourier co-efficients; ref: Parseval theorem for DFT
    B1_power<-mean(xfft$strength[B1_ix]^2,na.rm = TRUE)
    B2_power<-mean(xfft$strength[B2_ix]^2,na.rm = TRUE)
    B3_power<-mean(xfft$strength[B3_ix]^2,na.rm = TRUE)
    B4_power<-mean(xfft$strength[B4_ix]^2,na.rm = TRUE)
    B5_power<-mean(xfft$strength[B5_ix]^2,na.rm = TRUE)
    
    total_power<-B1_power+
      B2_power+
      B3_power+
      B4_power+
      B5_power
    
    #rm(B1_power,B2_power,B3_power,B4_power,B5_power)
    #rm(B1_ix,B2_ix,B3_ix,B4_ix,B5_ix)   
    
    B1_Power<-B1_power/total_power
    B2_Power<-B2_power/total_power
    B3_Power<-B3_power/total_power
    B4_Power<-B4_power/total_power
    B5_Power<-B5_power/total_power
  
    EMG_5Band<-c(B1_Power,B2_Power,B3_Power,B4_Power,B5_Power)
    }  
  
  # An approach to EEG Band Power calculation mentioned here 
  # https://raphaelvallat.com/bandpower.html. Apparently you have to 
  # integrate the power derived over the frequency interval with simpsons
  if(BandIntegrate==TRUE){
    # Frequency delineation into EEG bands using published cut-offs
    B1_ix<-which(xfft$frequency>=B1[1] & xfft$frequency<=B1[2])
    B2_ix<-which(xfft$frequency>=B2[1] & xfft$frequency<=B2[2])
    B3_ix<-which(xfft$frequency>=B3[1] & xfft$frequency<=B3[2])
    B4_ix<-which(xfft$frequency>=B4[1] & xfft$frequency<=B4[2])
    B5_ix<-which(xfft$frequency>=B5[1] & xfft$frequency<=B5[2])
    
    # Extracting the frequency by indices obtained in previous step
    B1_band<-xfft$frequency[B1_ix]
    B2_band<-xfft$frequency[B2_ix]
    B3_band<-xfft$frequency[B3_ix]
    B4_band<-xfft$frequency[B4_ix]
    B5_band<-xfft$frequency[B5_ix]
    
    # Squaring Fourier co-efficients to get Power in mV^2
    B1_coefs<-xfft$strength[B1_ix]^2
    B2_coefs<-xfft$strength[B2_ix]^2
    B3_coefs<-xfft$strength[B3_ix]^2
    B4_coefs<-xfft$strength[B4_ix]^2
    B5_coefs<-xfft$strength[B5_ix]^2
    
    #Simpson integration module-------------------------#
    simpson <- function(x, y)
    {
      if(length(x) < 5)
        stop("Must have at least 5 values")
      if(length(x) %% 2 == 0)
      {
        x<-x[-length(x)]
        y<-y[-length(y)]
      }
      #stop("Number of values must be odd")
      ord <- order(x)
      x <- x[ord]
      y <- y[ord]
      diffs <- diff(x)
      delta <- mean(diffs)
      if((max(diffs) - min(diffs))/delta > 1e-6)
        stop("X-values must be equally spaced")
      coefs <- c(1, 4, rep(c(2, 4), (length(x) - 3)/2), 1)
      sum(coefs*y)*delta/3
    }
    #---------------------------------------------------#
    
    # Simpson integration over all the coefficients in each EEG band and frequency vector
    B1_power<-simpson(x = B1_band,y = B1_coefs)
    B2_power<-simpson(x = B2_band,y = B2_coefs)
    B3_power<-simpson(x = B3_band,y = B3_coefs)
    B4_power<-simpson(x = B4_band,y = B4_coefs)
    B5_power<-simpson(x = B5_band,y = B5_coefs)
    
    # Calculating the total power of the EEG PSD for relative power calculation downstream
    total_power<-B1_power+
      B2_power+
      B3_power+
      B4_power+
      B5_power
    #Sigma power is excluded in the calculation above because it overlaps with beta band
    # and we are trying to integrate over the whole spectrum piece/band by piece/band
    
    # Calculating relative power from the Simpson integral
    B1_Power<-B1_power/total_power
    B2_Power<-B2_power/total_power
    B3_Power<-B3_power/total_power
    B4_Power<-B4_power/total_power
    B5_Power<-B5_power/total_power
    
    EMG_5Band<-c(B1_Power,
                B2_Power,
                B3_Power,
                B4_Power,
                B5_Power)  
  }
  
  return(EMG_5Band)
}

#int_F<-EMG_PSD(EMG_segment = EMG_segment,BandIntegrate = FALSE)
#int_T<-EMG_PSD(EMG_segment = EMG_segment,BandIntegrate = TRUE)

#Welch_PSD(EMG_segment = EMG_segment,Band_Integrate = FALSE)

Welch_PSD<-function(EMG_segment,
                    Fs=128,
                    window_size=180,
                    slide_interval=10,Band_Integrate=TRUE){
  # take an EEG segment, what ever size it is
  # create windows of 180 seconds and slide of 10 seconds
  # generate data windows and for each window
  # do PowerSpectralDensity into 6 EEG bands
  # Do SimpsonIntegration on each of these 6 EEG bands to calculate power
  # Convert Band Power to relative power by dividing all EEG band powers by total EEG power
  # now you have 15 odd values of relative power of EEG band for each window
  # choose the median value of this distribution of values
  
  #Fs=128;window_size=180;slide_interval=10;
  
  win_size=window_size*Fs
  slide_int=slide_interval*Fs
  start_win_ix<-seq(from=1,
                    to=length(EMG_segment)-win_size,
                    by=slide_int)
  
  end_win_ix<- start_win_ix+win_size  
  indices<-cbind(start_win_ix,end_win_ix)
  windowed_EMG_matrix<-apply(indices,1,function(x){EMG_segment[x[1]:x[2]]})
  windowed_PSD_matrix<-apply(windowed_EMG_matrix,2,function(emg_seg){EMG_PSD(EMG_segment = emg_seg,Fs = Fs,BandIntegrate = Band_Integrate)})
  EMG_Segment_Power<-apply(windowed_PSD_matrix,1,function(x){median(x,na.rm=TRUE)})
  names(EMG_Segment_Power)<-paste("B",1:length(EMG_Segment_Power),sep="")
  return(EMG_Segment_Power)
}
