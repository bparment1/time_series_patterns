############################## Trend methods time seores #################### 
##
## Functions generated through various research projects (STA) and SESYNC research support.
## Performing trend analyses of time series data to with theil sen, OLS and Mann Kendall.
##
## DATE CREATED: 08/11/2017
## DATE MODIFIED: 05/15/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis 
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

####### This script contains the following functions:
#
#1) calculate_theil_sen_time_series

############################
###Loading R library and packages                                                      
#library(gstat) #spatial interpolation and kriging methods
library(sp) # spatial/geographic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
library(readxl) #functionalities to read in excel type data
library(sf) # spatial ojbects simple feature model implementation OGC
#library(gstat)
#library(spacetime)
library(mblm)

###### Functions used in this script and sourced from other files

calculate_theil_sen_trend <- function(y,mod_obj=F){
  #y,n,harmonic_val=NULL,mod_obj=F,figure=F
  #This function generates Theil Sen slope estimate using mblm package.
  #
  #INPUTS
  #1) i: column index for input data.frame, used as y variable (trend variable)
  #2) data_df: input data containing time series
  #4) out_dir: output directory
  #5) out_suffix: output suffix appended to written output on drive
  #
  #OUTPUTS
  #1)
  #
  #TO DO: 
  
  ############# BEGIN FUNCTION #############
  
  #setting up the data input: data.frame with 2 columns
  time_index <- 1:length(y) #

  #stores dates for later processing
  #if(inherits(data_df,"zoo")){
  #  dates_val <- date(df_mblm)
  #}else{
  #  dates_val <- NULL
  #}
  
  dates_val <- NULL # set this later
  
  df_val <- data.frame(y=y,time_index=time_index) #convert to data.frame since it was a zoo df
  #df_mblm$time_index <- time_index
  #names(df_mblm)
  
  if(method="theil_sen"){
    df_val <- na.omit(df_val) # removes NA, might need to think about that later on
    ### automate formula input
    formula_str <- paste0(names(df_val)[1]," ~ ","time_index")
    formula_val <- as.formula(formula_str) #transform object into formula
    
    mod_obj <- mblm(formula_val,df_val)
  }
  
  if(method="ols"){
    #df_val <- na.omit(df_val) # removes NA, might need to think about that later on
    ### automate formula input
    formula_str <- paste0(names(df_val)[1]," ~ ","time_index")
    formula_val <- as.formula(formula_str) #transform object into formula
    
    mod_obj <- lm(formula_val,df_val)
    
  }
  
  ##### Extract information from model object mblm
  
  #slope_theil_sen <- coef(mod_mblm)[2]
  #intercept_theil_sen <- coef(mod_mblm)[1]
  
  slope <- coef(mod_obj)[2]
  intercept <- coef(mod_obj)[1]
  
  #ID_ts <- subset_name
  method_str <- "theil_sen"
  
  if(!is.null(dates_val)){
    nt <- length(dates_val)
    start_date <- dates_val[1]
    end_date <- dates_val[nt]
  }else{
    nt <- NA
    start_date <- NA
    end_date <- NA
  }
  
  duration_ts <- nt
  slope_sign <- sign(slope_theil_sen)
  
  df_trend <- data.frame(ID_ts=ID_ts,
                             intercept=intercept_theil_sen,
                             slope=slope_theil_sen,
                             slope_sign=slope_sign,
                             method=method_str,
                             start_date= start_date,
                             end_date = end_date,
                             duration=nt)
  
  #### Prepare object to return
  obj_theil_sen <- list(mod_mblm,df_trend)
  names(obj_theil_sen) <- c("mod_mblm","df_trend")
  
  ##### save to disk
  
  if(mod_obj==TRUE){
    #obj_theil_sen_filename <- file.path(out_dir,paste("theil_sen_obj_",subset_name,"_",out_suffix,".RData",sep=""))
    #save(obj_theil_sen,file= obj_theil_sen_filename)
  } 
  
  return(obj_theil_sen)
}

########################## End of script #######################################

