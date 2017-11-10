
## Functions used in the processing of data from google search on species for the animals-trade project at SESYNC.

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(raster)                              # raster functions and spatial utilities
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(plyr)                                # Various tools including rbind.fill
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots


###### Functions used in this script

import_data_ts <- function(infile_name,in_dir=".",scaling=1,n_col_start_date=3,start_date="2012-11-01",end_date=NULL,out_dir=".", out_suffix=""){  #new data start in Nov. 2012
#import_data_ts <- function(infile_name,in_dir=".",scaling=1,n_col_start_date=4,start_date="2011-01-01",end_date=NULL,out_dir=".", out_suffix=""){
  #Import data sent by google and format columns for future processing and visualization.
  
  ##### BEGIN FUNCTION ####
  df <- read.table(file.path(in_dir,infile_name),sep=",",fill=T,quote="",header=F)
  
  #remove first row that contains names
  #rename columns using:
  n_col <- ncol(df)
  nt <- n_col - n_col_start_date  #we are dropping the last date (column) because it is often incomplete month #### NO LONGER NECESSARY
  #nt <- 60
  
  range_dates <- seq.Date(from=as.Date(start_date),by="month",
                          length.out = nt)
  range_dates_str <- as.character(range_dates)
  #class(range_dates)
  #start_date
  #names_col <- c("g_id","sci_name","country",range_dates_str)
  names_col <- c("sci_name","country",range_dates_str)  # g_id no longer included
  df <- df[-1,-n_col] #remove the first row with incomplete header  
    names(df) <- names_col
  df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])
  

  out_filename_tmp <- sub(extension(infile_name),"",infile_name)#
  out_filename <- file.path(out_dir,paste0(out_filename_tmp,"_",out_suffix,".csv"))
  write.table(df,file=out_filename,sep=",")
  return(out_filename)
}

############################# END OF SCRIPT ########################