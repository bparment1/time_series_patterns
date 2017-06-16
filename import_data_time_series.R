############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/31/2017
## DATE MODIFIED: 06/16/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: fixing for import of new google data
##

###################################################
#

###### Library used

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate) 

###### Functions used in this script

##### Functions used in this script 

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

functions_processing_data_script <- "processing_data_google_search_time_series_functions_06162017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2
source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/nfs/edaut-data/Time Series MARSS"
#ARGS 2
infile_name <- "vert_sp_gst_original_06122017.csv"
#ARGS 3
start_date="2004-01-01"
#ARGS 4
end_date <- NULL
#ARGS 5
scaling_factor <- 1000
#ARGS 6
out_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/outputs" #parent directory where the new output directory is located
#out_dir <- "/nfs/edaut-data/Time Series MARSS/outputs"
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <- "vert_sp_gst_06162017"
#ARGS 9
n_col_start_date <- 4

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

######### Import data #################

#undebug(import_data_ts)

data_ts_filename <- import_data_ts(infile_name = infile_name,
                          in_dir=in_dir,
                          scaling=scaling,
                          n_col_start_date=4,
                          start_date=start_date,
                          end_date=NULL,
                          out_dir=out_dir, 
                          out_suffix=out_suffix)

data_ts <- read.table(data_ts_filename,sep=",",header=T)
  
######################## END OF SCRIPT ##############################