############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/23/2017
## DATE MODIFIED: 05/26/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: adding output dir function
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


calculate_theil_sen_time_series <- function(i,data_df,save_opt=FALSE,out_dir=".",out_suffix=""){
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
  time_index <- 1:nrow(data_df) #
  subset_name <- names(data_df)[i]
  df_mblm <- subset(data_df,select=subset_name) #subset using relevant name
  
  #stores dates for later processing
  if(inherits(data_df,"zoo")){
    dates_val <- date(df_mblm)
  }else{
    dates_val <- NULL
  }
  
  df_mblm <- as.data.frame(df_mblm) #convert to data.frame since it was a zoo df
  df_mblm$time_index <- time_index
  #names(df_mblm)
  
  ### automate formula input
  formula_str <- paste0(names(df_mblm)[1]," ~ ","time_index")
  formula_mblm <- as.formula(formula_str) #transform object into formula
  
  mod_mblm <- mblm(formula_mblm,df_mblm)
  
  ##### Extract information from model object mblm

  slope_theil_sen <- coef(mod_mblm)[2]
  intercept_theil_sen <- coef(mod_mblm)[1]
  ID_ts <- subset_name
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
  
  df_theil_sen <- data.frame(ID_ts=ID_ts,
                             intercept=intercept_theil_sen,
                             slope=slope_theil_sen,
                             slope_sign=slope_sign,
                             method=method_str,
                             start_date= start_date,
                             end_date = end_date,
                             duration=nt)
  
  #### Prepare object to return
  obj_theil_sen <- list(mod_mblm,df_theil_sen)
  names(obj_theil_sen) <- c("mod_mblm","df_theil_sen")
  
  ##### save to disk
  
  if(save_opt==TRUE){
    obj_theil_sen_filename <- file.path(out_dir,paste("theil_sen_obj_",subset_name,"_",out_suffix,".RData",sep=""))
    save(obj_theil_sen,file= obj_theil_sen_filename)
  } 
  
  return(obj_theil_sen)
}


################### End of script ################
