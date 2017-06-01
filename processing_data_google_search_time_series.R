############### SESYNC Research Support: Animals Trade ########## 
## Processing data from google search on species for the animals-trade project at SESYNC.
## 
## DATE CREATED: 05/31/2017
## DATE MODIFIED: 06/01/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut 
## Version: 1
## PROJECT: Animals trade by Elizabeth Daut
## ISSUE: 
## TO DO:
##
## COMMIT: fix naming problem in data.frame related to dates
##
## Links to investigate:

###################################################
#

###### Library used
###Loading R library and packages                                                      

library(raster)                 # loading the raster package
library(gtools)                 # loading R helper programming tools/functions
library(sp)                     # spatial objects in R
library(gplots)                 # plotting functions such as plotCI
library(rgdal)                  # gdal driver for R
library(RColorBrewer)           # color scheme, palettes used for plotting
library(gdata)                  # read different format (including .xlsx)
library(plotrix)                # plot options and functions 
library(rasterVis)              # raster visualization
library(colorRamps)             # contains matlab.like palette
library(zoo)                    # time series objects and methods
library(xts)                    # extension of time series objects
library(BMS)                    # contains hex2bin and bin2hex
library(bitops)                 # bit operations
library(gtools)                 #
library(maptools)               #
library(rgeos)                  # spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(sphet)                  # spatial analyis, regression eg.contains spreg for gmm estimation
library(forecast)               # arima and other time series methods
library(lubridate)              # date and time handling tools
library(parallel)               # access to parallelization functions
require(RCurl)
require(stringr)
require(XML)

###### Functions used in this script sourced from other files

#function_rainfall_time_series_NEST_analyses <- "rainfall_time_series_NEST_function_12112015.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/NEST/R_NEST" #path to script #PARAM 
#script_path <- "/home/parmentier/Data/rainfall/NEST"
#source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.

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

## Add processing of downloaded files!!
## 1) crop and reproject if needed
## 2) creation of TimeRaster object with summaries?
## 3)
### Other functions ####

function_processing_data <- "processing_data_google_search_time_series_functions_06012017b.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 
source(file.path(script_path,function_processing_data)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data" #local bpy50 , param 1
#in_dir <- "/home/parmentier/Data/rainfall/NEST" #NCEAS, param 
out_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/outputs" #param 2

#Elizabeth inputs
#in_dir <- "/nfs/edaut-data/Time Series MARSS"
#out_dir <- "/nfs/edaut-data/Time Series MARSS/outputs"

n_col_start_date <- 4
scaling <- 1000
infile_name <- "gst40k_original_5262017.csv"

start_date <- "2004-01-01"
end_date <- NULL # param 4

num_cores <- 4 #param 8
create_out_dir_param=TRUE # param 9

NA_value <- -9999 # param 10
NA_flag_val <- NA_value #param 11

out_suffix <-"processing_ts_06012017" #output suffix for the files and ouptu folder #param 12

#download_file <- FALSE #param 14
#unzip_files <- F #param 15

############## START SCRIPT ############################

######### PART 0: Set up the output dir ################

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

### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

#############

debug(import_data_ts)
lf_processed <- import_data_ts(infile_name,
                               in_dir=in_dir,
                               scaling=1,
                               n_col_start_date=n_col_start_date,
                               start_date="2004-01-01",
                               end_date=NULL,
                               out_dir=out_dir, 
                               out_suffix=out_suffix)

df_test <- read.table(lf_processed[1])
df_test <- read.table(lf_processed[1])

df_test <- read.table(lf_processed[1],sep=",",header = T)
names(df)
names_test <- as.character(seq(as.Date(start_date),by = "month",length.out=160))

library(assertive.code)
is_valid_variable_name(names_test)


############################# END OF SCRIPT ##############################