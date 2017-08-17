############### SESYNC Research Support: Spatial demography ########## 
## Performing PCA on animals trade data at SESYNC.
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 08/17/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: Animals Trade, Elizabeth Daut
## ISSUE: 
## TO DO:
##
## COMMIT: initial pca run for animals trade data
##
## Links to investigate:

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
library(ggplot2)
library(lubridate)
library(dplyr)
library(forecast)

###### Functions used in this script and sourced from other files

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

###### Functions used in this script

functions_time_series_analyses_script <- "time_series_functions_08012017.R" #PARAM 1
functions_processing_data_script <- "processing_data_google_search_time_series_functions_07202017.R" #PARAM 1
functions_time_series_cycles_analyses_script <- "time_series_cycles_analyses_functions_08172017.R" #PARAM 1


#script_path <- "C:/Users/edaut/Documents/gst_ts" #path to script #PARAM 2
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2
source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.

function_pca_eof <- "pca_eof_functions_08022017.R" #PARAM 1
source(file.path(script_path,function_pca_eof)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
#in_dir <- "C:/Users/edaut/Documents/gst_ts"
in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data"

#ARGS 2
#infile_name <- "vert_sp_gst_original_08162017.csv"
#infile_name_gst_totals <- "total_monthly_gst_averages.csv"
#use test data:
infile_name <- "dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt"

## Use data with known cycles:

#ARGS 3
start_date <- "2004-01-01"
#ARGS 4
end_date <- NULL
#ARGS 5
scaling_factor <- 100 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
#scaling_factor <- 1000 
#ARGS 6
#out_dir <- NULL #"C:/Users/edaut/Documents/gst_ts/outputs" #parent directory where the new output directory is placed
out_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/outputs"
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <-"cycles_test_08172017" #output suffix for the files and ouptut folder #param 12

#ARGS_9
n_col_start_date <- 4
#ARGS 10
num_cores <- 2

################# START SCRIPT ###############################

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

#data_df <- read.table(file.path(in_dir,infile_names),sep=",",header=T)

############### PART 1: Imported and time series transformation #####
## Step 1: read in the data and generate time stamps

#infile_name <- "vert_sp_gst_original_08162017.csv"
#data_ts_filename <- import_data_ts(infile_name = infile_name,
#                                   in_dir = in_dir,
#                                   scaling = scaling,
#                                   n_col_start_date=4,
#                                   start_date = start_date,
#                                   end_date=NULL,
#                                   out_dir = out_dir,
#                                   out_suffix = out_suffix)

#df_original <- read.table(data_ts_filename,sep=",",fill=T,header=T,check.names = F)

#Set up dates as column names  ###############  
#IS THIS NEEDED b/c ALREADY IN PROCESSING??? ###############

# range_dates <- names(df_original)[n_col_start_date:ncol(df_original)]
# range_dates_str <- as.character(range_dates)
# range_dates<- ymd(range_dates_str)
# 
# #Transform and scale data
# df_subset <- df_original
# df_ts <- t(df_original[n_col_start_date:ncol(df_subset)])
# df_ts <- as.data.frame(df_ts)
# #dim(df_ts)
# df_ts_scaled <- df_ts[,]*scaling_factor  
# 
# #create ts zoo object
# df_ts <- zoo(df_ts_scaled,range_dates)
# #class(df_ts)
# #combine country-species as column names
# names_countries <- as.character(df_subset$country)
# names_species <- as.character(df_subset$sci_name)
# names_species <- sub(" ","_",names_species)
# names_col <- paste(names_countries,names_species,sep="_")
# names(df_ts) <- names_col
# #View(df_ts)

#http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html

infile_name <- file.path(in_dir,infile_name)
data_df <- read.table(infile_name,header=T,sep=",",stringsAsFactors = F)
View(data_df)
names(data_df)
dim(data_df)
data_df <- data_df[,1:230]
df_ts <- (t(data_df))
dim(df_ts)
date_range <- c("2001.01.01","2010.12.31") #PARAM 15, NDVI Katrina
range_dates <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina
class(range_dates)

df_ts <- zoo(df_ts,range_dates)

(df_ts[1:10,])

plot(df_ts[1,])
dim(df_ts)

xs <- seq(-2*pi,2*pi,pi/100)
wave.1 <- sin(3*xs)
wave.2 <- sin(10*xs)
par(mfrow = c(1, 2))
plot(xs,wave.1,type="l",ylim=c(-1,1)); abline(h=0,lty=3)
plot(xs,wave.2,type="l",ylim=c(-1,1)); abline(h=0,lty=3)

#https://www.r-bloggers.com/smoothing-techniques-using-basis-functions-fourier-basis/
#https://stats.stackexchange.com/questions/1207/period-detection-of-a-generic-time-series/1214#1214

findfrequency

find.freq <- function(x)
?spec.ar



{
  n <- length(x)
  spec <- spec.ar(c(x),plot=FALSE)
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        period <- round(1/spec$freq[nextmax])
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  return(period)
}


############################## END OF SCRIPT #############################################

