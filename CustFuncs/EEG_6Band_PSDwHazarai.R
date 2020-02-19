EEG_6Band_PSD<-function(multiChEEG,
                        Fs=128,
                        lower=0.5,
                        upper=63,
                        delta=c(0.5,3.999),
                        theta=c(4,7.999),
                        alpha=c(8,11.999),
                        sigma=c(12,15),
                        beta=c(12,29.999),
                        gamma=c(30,64),add_derived=FALSE,BandIntegrate=FALSE){
  
  # multiChEEG has rows of samples and columns correspond to EEG channels
  require(eegkit)
  xfft <- eegfft(multiChEEG, Fs, lower, upper)
  if(BandIntegrate==FALSE){  
  delta_ix<-which(xfft$frequency>=delta[1] & xfft$frequency<=delta[2])
  theta_ix<-which(xfft$frequency>=theta[1] & xfft$frequency<=theta[2])
  alpha_ix<-which(xfft$frequency>=alpha[1] & xfft$frequency<=alpha[2])
  sigma_ix<-which(xfft$frequency>=sigma[1] & xfft$frequency<=sigma[2])
  beta_ix<-which(xfft$frequency>=beta[1] & xfft$frequency<=beta[2])
  gamma_ix<-which(xfft$frequency>=gamma[1] & xfft$frequency<=gamma[2])
  
  # The power of a Discrete Fourier Transform is the MEAN of 
  # the sum of squares of its Fourier co-efficients; ref: Parseval theorem for DFT
  # Band_power_AllCH is a matrix of FFT co-efficients on rows and EEG Channels on columns
  # the operation below gives vectors of the average of squared co-efficients in each EEGband
  # across 14 channels, adding them together gives a vector of total power in 14 channels
  delta_power_AllCH<-colMeans(xfft$strength[delta_ix,]^2,na.rm = TRUE)
  theta_power_AllCH<-colMeans(xfft$strength[theta_ix,]^2,na.rm = TRUE)
  alpha_power_AllCH<-colMeans(xfft$strength[alpha_ix,]^2,na.rm = TRUE)
  sigma_power_AllCH<-colMeans(xfft$strength[sigma_ix,]^2,na.rm = TRUE)
  beta_power_AllCH<-colMeans(xfft$strength[beta_ix,]^2,na.rm = TRUE)
  gamma_power_AllCH<-colMeans(xfft$strength[gamma_ix,]^2,na.rm = TRUE)
  
  # Sigma band is not included in the total power calculations since it already overlaps beta
  total_power_perCH<-delta_power_AllCH+
    theta_power_AllCH+
    alpha_power_AllCH+
    beta_power_AllCH+
    gamma_power_AllCH
  
  Delta_Power<-delta_power_AllCH/total_power_perCH
  Theta_Power<-theta_power_AllCH/total_power_perCH
  Alpha_Power<-alpha_power_AllCH/total_power_perCH
  Sigma_Power<-sigma_power_AllCH/total_power_perCH
  Beta_Power<-beta_power_AllCH/total_power_perCH
  Gamma_Power<-gamma_power_AllCH/total_power_perCH
  }  
  
  # An approach to EEG Band Power calculation mentioned here 
  # https://raphaelvallat.com/bandpower.html. Apparently you have to 
  # integrate the power derived over the frequency interval with simpsons
  if(BandIntegrate==TRUE){
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
  }
  
  AllChannel_RP_EEGBand<-rbind(Delta_Power,Theta_Power,Alpha_Power,Sigma_Power,Beta_Power,Gamma_Power)
  
  # Derived features are simple ratios of EEG Band Powers or ratios of sum of different
  # EEG band powers that can express informative activity in terms of how the power distribution
  # changes within different EEG bands. The module below will add those derived features if activated
  if(add_derived==TRUE){
  Delta_By_Beta_Power<-Delta_Power/Beta_Power
  AlphaTheta_By_Beta_Power<-(Alpha_Power+Theta_Power)/Beta_Power
  Alpha_By_Beta_Power<-Alpha_Power/Beta_Power
  AlphaTheta_By_AlphaBeta_Power<-(Alpha_Power+Theta_Power)/(Alpha_Power+Beta_Power)
  Theta_By_Beta_Power<-Theta_Power/Beta_Power
  DeltaTheta_By_AlphaBetaSigma_Power<-(Delta_Power+Theta_Power)/(Alpha_Power+Beta_Power+Sigma_Power)
  Theta_By_Alpha_Power<-Theta_Power/Alpha_Power
  AllChannel_RP_EEGBand<-rbind(AllChannel_RP_EEGBand,
                               Delta_By_Beta_Power,
                               AlphaTheta_By_Beta_Power,
                               Alpha_By_Beta_Power,
                               AlphaTheta_By_AlphaBeta_Power,
                               Theta_By_Beta_Power,
                               DeltaTheta_By_AlphaBetaSigma_Power,
                               Theta_By_Alpha_Power)
  }
  
  return(AllChannel_RP_EEGBand)
}