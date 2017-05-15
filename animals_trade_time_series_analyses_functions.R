############### SESYNC Research Support: Animals Trade project ########## 
## General purpose function to analyze time series data.
## Includes:
## 1) PCA analyses (both time series and other)
## 2) Theil Sen slope estimator
## 3) Mann Kenall trend test 
##
## DATE CREATED: 05/22/2017
## DATE MODIFIED: 05/22/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 2
## PROJECT: Animals trade
## Issue: #21977
## link: https://base.sesync.org/issues/21977
## 
## TO DO:
##
## COMMIT: animals trade, issue #21977, time series processing pca and more
##

###################################################
#

### Loading R library and packages  

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate)
library(nlme) #mixed models and GLS
library(forecast) #moving average, ARIMA and other tools
library(psych) #pca/eigenvector decomposition functionalities
library(GPArotation)

###### Functions used in this script

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


################### End of script ################