Welch_PSD<-function(multiChEEG,
                          Fs=128,
                          window_size=180,
                          slide_interval=10){
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
                    to=nrow(multiChEEG)-win_size,
                    by=slide_int)
  
  end_win_ix<- start_win_ix+win_size  
  indices<-cbind(start_win_ix,end_win_ix)
  windowed_EEG_list<-apply(indices,1,function(x){multiChEEG[x[1]:x[2],]})
  windowed_PSD_list<-lapply(windowed_EEG_list,function(eeg_seg){EEG_6Band_PSD(multiChEEG = eeg_seg,Fs = Fs)})
  # EEG_6Band_PSD takes a time by amplitude multi-channel EEG data and calculates the PSD
  # using an FFT followed by simpsons integration over the band region of the fft power
  # we do this for 10 second shifts over 180 second windows and get a list of matrices
  # that have EEG bands on rows and Channel Names on columns, the entries are PSD power
  
  # Now we have to take a median over these matrices and aggregate their values
  # We melt the list of matrices to a list of dataframes with columns 
  # EEG Band, EEG Channel, EEG power (lapply(melt))
  # We convert the list of melted dataframes to a single dataframe using ldply
  # Now each EEG band has a distribution of Band Powervalues calculated from the moving windows
  # we want to take the median of this distribution for each EEG Band and EEG Channel
  # so now we have an aggregated dataframe which has EEG Band, EEG Channel, EEG Power
  # and the EEG Power is the median of the distribution of values is the previous dataframe
  # After this, we convert this aggregated data frame into a matrix with EEG band on rows
  # and EEG channels on columns and drop the EEG band column and add it to rownames
  # transpose this matrix before returning to preserve compatibility with existing workflow.
  meltry<-lapply(windowed_PSD_list,reshape2::melt)
  df_meltry<-plyr::ldply(meltry)
  agg_df_meltry<-aggregate(x = df_meltry$value,by=list(df_meltry$Var1,df_meltry$Var2),FUN = median,na.rm=TRUE,drop=TRUE)
  names(agg_df_meltry)<-c("EEG_Band","EEG_Channel","Band_Power")
  EEGbandPower<-tidyr::spread(data=agg_df_meltry,key = EEG_Band,value = Band_Power)
  rownames(EEGbandPower)<-EEGbandPower[,1]
  EEGbandPower_df<-t(EEGbandPower[,c(-1)])
  return(EEGbandPower_df)
}

EEG_6Band_PSD<-function(multiChEEG,
                        Fs=128,
                        lower=0.5,
                        upper=63,
                        delta=c(0.5,3.999),
                        theta=c(4,7.999),
                        alpha=c(8,11.999),
                        sigma=c(12,15),
                        beta=c(12,29.999),
                        gamma=c(30,64)){
 
  # Simpson integration  
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
  
  # multiChEEG has rows of samples and columns correspond to EEG channels
  require(eegkit)
  xfft <- eegfft(multiChEEG, Fs, lower, upper)
  
  # Frequency delineation into EEG bands using published cut-offs
  delta_ix<-which(xfft$frequency>=delta[1] & xfft$frequency<=delta[2])
  theta_ix<-which(xfft$frequency>=theta[1] & xfft$frequency<=theta[2])
  alpha_ix<-which(xfft$frequency>=alpha[1] & xfft$frequency<=alpha[2])
  sigma_ix<-which(xfft$frequency>=sigma[1] & xfft$frequency<=sigma[2])
  beta_ix<-which(xfft$frequency>=beta[1] & xfft$frequency<=beta[2])
  gamma_ix<-which(xfft$frequency>=gamma[1] & xfft$frequency<=gamma[2])
  
  # Extracting the frequency by indices obtained in previous step
  delta_band<-xfft$frequency[delta_ix]
  theta_band<-xfft$frequency[theta_ix]
  alpha_band<-xfft$frequency[alpha_ix]
  sigma_band<-xfft$frequency[sigma_ix]
  beta_band<-xfft$frequency[beta_ix]
  gamma_band<-xfft$frequency[gamma_ix]
  
  # Squaring Fourier co-efficients to get Power in mV^2
  delta_coefs_AllCH<-xfft$strength[delta_ix,]^2
  theta_coefs_AllCH<-xfft$strength[theta_ix,]^2
  alpha_coefs_AllCH<-xfft$strength[alpha_ix,]^2
  sigma_coefs_AllCH<-xfft$strength[sigma_ix,]^2
  beta_coefs_AllCH<-xfft$strength[beta_ix,]^2
  gamma_coefs_AllCH<-xfft$strength[gamma_ix,]^2
  
  # Simpson integration over all the coefficients in each EEG band and frequency vector
  delta_power_AllCH<-apply(delta_coefs_AllCH,2,function(CH_N){simps_integ=simpson(x = delta_band,y = CH_N);return(simps_integ)})
  theta_power_AllCH<-apply(theta_coefs_AllCH,2,function(CH_N){simps_integ=simpson(x = theta_band,y = CH_N);return(simps_integ)})
  alpha_power_AllCH<-apply(alpha_coefs_AllCH,2,function(CH_N){simps_integ=simpson(x = alpha_band,y = CH_N);return(simps_integ)})
  sigma_power_AllCH<-apply(sigma_coefs_AllCH,2,function(CH_N){simps_integ=simpson(x = sigma_band,y = CH_N);return(simps_integ)})
  beta_power_AllCH<-apply(beta_coefs_AllCH,2,function(CH_N){simps_integ=simpson(x = beta_band,y = CH_N);return(simps_integ)})
  gamma_power_AllCH<-apply(gamma_coefs_AllCH,2,function(CH_N){simps_integ=simpson(x = gamma_band,y = CH_N);return(simps_integ)})
  
  # Calculating the total power of the EEG PSD for relative power calculation downstream
  total_power_perCH<-delta_power_AllCH+
    theta_power_AllCH+
    alpha_power_AllCH+
    beta_power_AllCH+
    gamma_power_AllCH
  #Sigma power is excluded in the calculation above because it overlaps with beta band
  # and we are trying to integrate over the whole spectrum piece/band by piece/band
  
  # Calculating relative power from the Simpson integral
  Delta_Power<-delta_power_AllCH/total_power_perCH
  Theta_Power<-theta_power_AllCH/total_power_perCH
  Alpha_Power<-alpha_power_AllCH/total_power_perCH
  Sigma_Power<-sigma_power_AllCH/total_power_perCH
  Beta_Power<-beta_power_AllCH/total_power_perCH
  Gamma_Power<-gamma_power_AllCH/total_power_perCH
  
 AllChannel_RP_EEGBand<-rbind(Delta_Power,
                              Theta_Power,
                              Alpha_Power,
                              Sigma_Power,
                              Beta_Power,
                              Gamma_Power)
 #rownames(AllChannel_RP_EEGBand)<-c("Delta_Power","Theta_Power","Alpha_Power","Sigma_Power","Beta_Power","Gamma_Power")
  return(AllChannel_RP_EEGBand)
}



