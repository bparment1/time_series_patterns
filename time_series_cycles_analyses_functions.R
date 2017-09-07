############### SESYNC Research Support: Animals Trade project ########## 
#### General functions to examine and detect periodic cycles such as seasonality.
## 
## DATE CREATED: 08/17/2017
## DATE MODIFIED: 09/07/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO: - Fourier
##        - windowing function
##        - generate artificial dataset
##        - windowed Fourier
##        - multitaper methods 
##
## COMMIT: adding multitaper methods for spectrum analyses
##

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)                             # Plotting functionalities
library(lubridate)                           # Dates manipulation functionalities
library(dplyr)                               # Data wrangling
library(forecast)                            # ARIMA and other time series methods
library(multitaper)                          # Multitaper estimation of spectrum
library(GeneCycle)                           # Fisher test for harmonics and Time series functionalities
library(xts)                                 # Extension for time series object and analyses
library(zoo)                                 # Time series object and analysis
library(mblm)                                # Theil Sen estimator


###### Functions used in this script

##create an output directory
create_dir_fun <- function(out_dir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

sine_structure_fun <-function(x,T,phase_val,a,b){
  #Create sine for a one dimensional series
  #Note that sine function uses radian unit.
  #a=amplitude
  #b=mean or amplitude 0 of the series
  #T= stands for period definition
  #phase=phase angle (in radian!!)
  
  y <- a*sin((x*2*pi/T)+ phase_val) + b
}


generate_dates_by_step <-function(start_date,end_date,step_date){
  #library(xts) declare out of this function
  #library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(as.numeric(strftime(st,"%Y")),as.numeric(strftime(en,"%Y"))) #extract year
  
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
    #paste(yday(ll,)
  }
  
  #
  dates_modis <-as.Date(unlist((ll_list))) 
  #wiht 001
  dates_DOY_modis <- paste(year(dates_modis),sprintf("%03d", yday(dates_modis)),sep="")
  dates_obj <- list(dates_modis,dates_DOY_modis)
  names(dates_obj) <- c("dates","doy")  
  return(dates_obj)
}


# compute the Fourier Transform

spectrum_analysis_fft_run <- function(x){
  #
  #This functions examines a single time series and generate spectrum summary to detect harmonics 
  #INPUT:
  #x : a time series as a vector
  #OUTPUT:
  # 1) periodogram_obj:
  # 2) spectrum_obj: list of the follwing element
  # - "p"
  # - "freq_df"
  # - "ranked_freq_df",
  # - "n_orig"
  # - "n_used"
  #
  
  #### Begin #####
  
  ### Part I: use periodogram, line spectrum to find harmonics contribution
  
  p <- periodogram(x,fast=F)
  #spectrum_val <- spectrum(as.numeric(x)) #not padding of power 2
  
  n_orig <- length(x)
  freq_df <- data.frame(freq=p$freq, spec=p$spec)
  total_variance <- sum(freq_df$spec)
  freq_df$index <- as.numeric(row.names(freq_df))
  freq_df$period <- 1/freq_df$freq
  freq_df$period_orig <- n_orig/freq_df$index
  freq_df$variance <- freq_df$freq/total_variance*100
  ranked_freq_df <- freq_df[order(-freq_df$spec),]
  
  #barplot(freq_df$variance)
  #p$orig.n* as.numeric(rownames(top2)
  n_used <- p$n.used
  
  peaks_df <- as.data.frame(findpeaks(p$spec))
  names(peaks_df) <- c("val","loc","start","end")
  
  periodogram_obj <- list(p,freq_df,ranked_freq_df,p$orig.n,p$n.used,peaks_df)
  names(periodogram_obj) <- c("p","freq_df","ranked_freq_df","n_orig","n_used","peaks_df")

  ### Part II: use spectrum, power density to find harmonics contribution

  spectrum_val <- spectrum(as.numeric(x),fast=F) #not padding of power 2
  #spectrum_val <- spectrum(as.numeric(x))
  length(spectrum_val$freq)

  #spectrum_val$freq #no option to remove padding
  
  n_orig <- length(x)
  
  freq_df <- data.frame(freq=spectrum_val$freq, spec=spectrum_val$spec)
  
  total_variance <- sum(freq_df$spec)
  freq_df$index <- as.numeric(row.names(freq_df))
  freq_df$period <- 1/freq_df$freq
  freq_df$period_orig <- n_orig/freq_df$index
  freq_df$variance <- freq_df$freq/total_variance*100
  ranked_freq_df <- freq_df[order(-freq_df$spec),]
  
  #barplot(freq_df$variance)
  #p$orig.n* as.numeric(rownames(top2)
  n_used <- p$n.used
  
  peaks_df <- as.data.frame(findpeaks(spectrum_val$spec))
  names(peaks_df) <- c("val","loc","start","end")
  
  spectrum_obj <- list(spectrum_val,freq_df,ranked_freq_df,p$orig.n,p$n.used,peaks_df)
  names(spectrum_obj) <- c("spectrum","freq_df","ranked_freq_df","n_orig","n_used","peaks_df")
  
  ## Prepare return object:
  
  spectrum_analyis_fft_obj <- list(periodogram_obj,spectrum_obj)
  names(spectrum_analyis_fft_obj) <- c("periodogram_obj","spectrum_obj")
  
  return(spectrum_analyis_fft_obj)
}


generate_harmonic <- function(i,coef_fft_df,a0){
  
  a <- coef_fft_df$amp[i] #amplitude in this case
  b <- a0 # this should be the average!!
  #T <- 230
  #T <- 10
  
  T <- coef_fft_df$T[i]
  #T <- 23
  
  phase_val <-coef_fft_df$phase[i] #get the phase from fft?
  
  n_val <- nrow(coef_fft_df)
  #x_input <- 1:230 #index sequence corresponding to timesetps for time series
  x_input <- 1:n_val
  
  ux <- sine_structure_fun(x_input,T,phase_val,a,b)
  #plot(ux,type="b")
  return(ux) #harmonics for frequency 
  
}


extract_harmonic_fft_parameters_run <- function(x){
  #
  ## This functions generate harmonics using fft from x (time) series.
  #INPUTS
  #1)x: series to process with fft
  #2)a0: average of series or reference point
  #2)selected_f: number of harmonics to retain, if NULL generate to 5 or 10 (depending on the inputs)
  #OUTPUTS
  #1) 
  #
  
  ########## BEGIN ############
  
  options(scipen=999)  #remove scientific writing
  
  x_fft <- fft(x,inverse=T) # transformed fft
  #Let j be the frequency
  # aj*cos(wj*t) + bj*sin(wjt)
  #
  #x_trans <- fft(x,inverse=T,fast=F) # transformed fft, no padding to get to power of 2
  
  #
  #pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/(N0 * xfreq)
  #N0 is the total number of steps
  
  ## To get the amplitude, you need to normalize by the number of element
  ## You multiply by two because it is symmetrical and the power is spread
  ## in both negative and positive
  #sqrt(a^2+b^2): this gets the magnitude of the imaginary vector
  mod_val <- sqrt((Im(x_fft)^2 + (Re(x_fft))^2)) #module of the imaginery number
  ## power spectrum:
  mod_val <- (mod_val/length(x_fft))*2 #normalized by length and multiply by 2
  amp_val <- (as.numeric(Mod(x_fft))/length(x_fft))*2
  amplitude <- abs(x_fft)/length(x_fft)*2 #
  
  phase <- Arg(x_fft) # atan(Im(x_fft)/Re(x_fft))
  
  #phase[11]
  #1.296784 
  #> (phase[11]/pi)*180
  #11 
  #74.30026 
  #> 360/23
  #[1] 15.65217
  #> 90-360/23
  #[1] 74.34783
  #> 
  
  #amp <- amp_val
  #amp[1] <- 0 #Set the amplitude to zero for the harmonic zero
  #https://www.mathworks.com/matlabcentral/answers/162846-amplitude-of-signal-after-fft-operation
  
  n <- length(x)
  n_half <- n/2 
  plot(amplitude,type="h")
  plot(1:n_half,amp[1:n_half],"h")
  
  #http://www.mathworks.com/help/matlab/ref/fft.html
  P2 = abs(x_trans/n)
  P1 = P2[(1:n)/2+1]
  P1(2:end-1) = 2*P1(2:end-1);
  
  phase <- as.numeric(Arg(x_trans))
  #phase[11]
  #barplot(phase)
  #amp_scaled <- 1/(amp_val/n)
  #amp_scaled <- n/amp_val
  #amp_scaled <- log(amp_val)
  #amp_scaled[1] <- 0
  #1/(113.070230382360506382611/230)
  
  n_selected <- n_half -1
  period_orig <- 230/(1:n_half-1)
  #period_orig <- 230/(1:n_half)
  harmonic_val <- 1:n_half
  frequency_val <- 1/period_orig
  freq_rad <- 2*pi*frequency_val
  phase_deg <- phase*180/pi
  variance_freq <- (amplitude^2)/2
  variance_perc <- (variance_freq/sum(variance_freq))*100
  
  coef_fft_df <- data.frame(harmonic_val,
                            period_orig,
                            frequency_val,
                            amplitude[1:n_half],
                            phase[1:n_half],
                            phase_deg[1:n_half],
                            variance_freq[1:n_half],
                            variance_perc[1:n_half])
  names(coef_fft_df) <- c("harmonic","period_orig","frequency",
                          "amplitude","phase","phase_deg","variance","var_percent")
  
  #x_in <- 1:230
  #amp[10]*sin(+ phase[10])
  

  return(coef_fft_df)
}


filter_frequency_and_generate_harmonics <- function(x,selected_f=NULL){
  
  ### Generate a sequence from sine
  
  if(is.null(selected_f)){
  #then take the top 10
  #type_spatialstructure[5] <- "periodic_x1"
  
  #debug(harmonic_analysis_fft_run)
  spectrum_analysis_fft_obj <- spectrum_analysis_fft_run(x) 
  
  ranked_freq_df <- spectrum_analysis_fft_obj$spectrum_obj$ranked_freq_df
  selected_f <- ranked_freq_df$freq[1:10]
  }
  ##

  ### Get ranked frequencies

  selected_coef_fft_df <- coef_fft_df[selected_f,]
  #i<-1
  #generate_harmonic(i,coef_fft_df=selected_coef_fft_df,a0)
  a0 <- mean(x)

  #may want to use diff or detrend?

  n_selected <- nrow(selected_coef_fft_df)

  harmonics_fft_list <- lapply(1:n_selected,
                             FUN= generate_harmonic,
                             coef_fft_df=selected_coef_fft_df,
                             a0=a0)
  #Can add them all together?

  harmonic_obj <- list(harmonics_fft_list,harmonic_fft_obj)
  
  return()

}


adding_temporal_structure <- function(list_param){
  #
  #This functions generate different temporal components.
  # 
  
  ## Parse input arguments
  
  nt <- list_param$nt
  phase <- list_param$phase
  temp_periods <- list_param$temp_periods
  amp <- list_param$amp
  temp_period_quadrature <- list_param$temp_period_quadrature
  random_component <- list_param$random_component #mean and sd used in rnorm
  
  ### Start #####
  
  #type_temporalstructure <- character (length=6)
  x <- as.numeric(1:nt)
  
  # Generate temporal signal with periodic signal
  # can have multiple periods!!!
  
  y1_list<-vector("list",length=length(temp_periods))
  
  ### if only one amplitude, then all frequencies have the same strength
  if(length(amp)==1){
    amp <- rep(amp,length(temp_periods))
  }
  
  for (i in 1:length(temp_periods)){
    T<-temp_periods[i]
    amp_val <- amp[i]
    
    #a<- 1 #amplitude in this case
    #b<- 0
    #T<- 12
    #phase_val <- 0
    x_input<- x
    
    a <- amp_val
    b <- 0
    phase_val <- phase
    
    #undebug(sine_structure_fun)
    y1 <- sine_structure_fun(x_input,T,phase_val,a,b)
    #sine_structure_fun <-function(x,T,phase_val,a,b){
    
    y1_list[[i]] <-  y1
  }
  names(y1_list) <- paste("t_period",temp_periods,sep="_")
  # Generate periodic Quadrature:
  
  #T <- temp_period_quadrature
  #y2_list<-vector("list",length=length(temp_periods))
  #y2_list[[1]] <- temp_structure_fun(x,T,pi/2)
  #y2_list[[2]] <- temp_structure_fun(x,T,0)
  
  #y2<-y2_list[[1]] + y2_list[[2]]
  
  # Generate temporal trends
  a <- 0.1
  b <- 0
  y3 <-  a*x + b  
  
  # Generate temporal temporal randomness
  y4 <- runif(nt) 
  
  #y5<- rnorm(nt)
  y5 <- rnorm(nt,random_component[1],random_component[2])
  
  #Prepare return object
  dfrm1 <-do.call(cbind,y1_list)
  #dfrm2 <- do.call(cbind,list(quadrature=y2,trend=y3,unif=y4,norm=y5))
  dfrm2 <- do.call(cbind,list(trend=y3,unif=y4,norm=y5))
  
  temp_pattern_dfrm <- (cbind(dfrm1,dfrm2)) #data.frame containing temporal patterns...
  temp_pattern_dfrm<-as.data.frame(temp_pattern_dfrm)
  #head(temp_pattern_dfrm)
  file_name<-paste("table_temporal_patterns","_",out_suffix,".txt",sep="")
  write.table(temp_pattern_dfrm,file=file_name,sep=",")
  
  return(temp_pattern_dfrm)
}


mtm_spectrum_analysis_fun <- function(){
  #This functions using the multitaper methods to find harmonics in signal.
  #
  #Functions needs to be update
  
  ##### Start ####
  
  #require(multitaper);
  data(willamette);
  resSpec <- spec.mtm(willamette, k=10, nw=5.0, nFFT = "default",
                      centreWithSlepians = TRUE, Ftest = TRUE,
                      jackknife = FALSE, maxAdaptiveIterations = 100,
                      plot = TRUE, na.action = na.fail) 
  
  resSpec <- spec.mtm(x_ts1, k=10, nw=5.0, nFFT = "default",
                      centreWithSlepians = TRUE, Ftest = TRUE,
                      jackknife = FALSE, maxAdaptiveIterations = 100,
                      plot = TRUE, na.action = na.fail) 
  
  resSpec <- spec.mtm(x_ts1_lm, k=10, nw=5.0, nFFT = "default",
                      centreWithSlepians = TRUE, Ftest = TRUE,
                      jackknife = FALSE, maxAdaptiveIterations = 100,
                      plot = TRUE, na.action = na.fail) 
  ### control the padding, request none
  resSpec <- spec.mtm(x_ts1_lm, k=10, nw=5.0, nFFT = 230,
                      centreWithSlepians = TRUE, Ftest = TRUE,
                      jackknife = FALSE, maxAdaptiveIterations = 100,
                      plot = TRUE, na.action = na.fail) 
  which.max(resSpec$spec)#harmonic 9 instead of harmonic 10, the tapering affect the peak
  which.min(resSpec$mtm$Ftest)
  
  length(resSpec$spec)
  
  ### control the padding, request none
  resSpec <- spec.mtm(x_ts1_lm, k=10, nw=9.0, nFFT = 230,
                      centreWithSlepians = TRUE, Ftest = TRUE,
                      jackknife = FALSE, maxAdaptiveIterations = 100,
                      plot = TRUE, na.action = na.fail) 
  which.max(resSpec$spec)#harmonic 9 instead of harmonic 10, the tapering affect the peak
  which.min(resSpec$mtm$Ftest)
  length(resSpec$spec)
  
  return()
  
}

################### End of script ################
