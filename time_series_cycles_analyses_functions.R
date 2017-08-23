############### SESYNC Research Support: Animals Trade project ########## 
#### General functions to examine and detect periodic cycles such as seasonality.
## 
## DATE CREATED: 08/17/2017
## DATE MODIFIED: 08/23/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO: - Fourier
##        - windowed Fourier
##
## COMMIT: adding functions from main script
##

###################################################
#

###### Library used

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate) 

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

harmonic_analysis_fft_run <- function(x){
  #
  #Function to compute 
  
  p <- periodogram(x)
  
  spectrum_val <- spectrum(as.numeric(x))
  n_orig <- length(x)
  freq_df <- data.frame(freq=p$freq, spec=p$spec)
  total_variance <- sum(freq_df$spec)
  freq_df$index <- as.numeric(row.names(freq_df))
  freq_df$period <- 1/freq_df$freq
  freq_df$period_orig <- n_used/freq_df$index
  freq_df$variance <- freq_df$freq/total_variance*100
  ranked_freq_df <- freq_df[order(-freq_df$spec),]
  
  #barplot(freq_df$variance)
  #p$orig.n* as.numeric(rownames(top2)
  n_used <- p$n.used
  
  ## Prepare return object:
  
  harmonic_fft_obj <- list(p,freq_df,ranked_freq_df,p$orig.n,p$n.used)
  names(harmonic_fft_obj) <- c("p","freq_df","ranked_freq_df","n_orig","n_used")
  
  return(harmonic_fft_obj)
}

################### End of script ################
