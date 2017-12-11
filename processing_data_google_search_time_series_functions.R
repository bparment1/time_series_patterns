############### SESYNC Research Support: Animals Trade ########## 
## Functions used in the processing of data from google search on species for the animals-trade project at SESYNC.
## 
## DATE CREATED: 05/31/2017
## DATE MODIFIED: 12/11/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut 
## Version: 1
## PROJECT: Animals trade by Elizabeth Daut
## ISSUE: 
## TO DO:
##
## COMMIT: documenting and testing function to import data sent by google
##
## Links to investigate:

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(raster)                              # raster functions and spatial utilities
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)

###### Functions used in this script

import_data_ts <- function(infile_name,version=1,in_dir=".",scaling=1,n_col_start_date=4,start_date="2011-01-01",end_date=NULL,out_dir=".", out_suffix=""){
  #
  #Goal: Import data sent by google and format columns for future processing and visualization.
  #
  #INPUTS:
  #1) infile_name: text file containing data in csv format
  #2) version: indicate which version, version 1 contains gid while version does not
  #3) in_dir: input directory
  #4) n_col_start_date: first column containing date information
  #5) start_date: starte date for the dataset, default is "2011-01-01"
  #6) end_date: end data for the dataset, default is NULL
  #7) out_dir: output directory, default is current dir "."
  #8) out_suffix: suffix string to add to outputs
  #OUTPUTS
  #1) out_filename: output filename
  
  ##### BEGIN FUNCTION ####
  
  ### Add quote="" otherwise EOF warning and error in reading
  df <- read.table(file.path(in_dir,infile_name),sep=",",fill=T,quote="",header=F,stringsAsFactors = F)
  
  #remove first row that contains names
  #rename columns using:
  n_col <- ncol(df)
  nt <- n_col - n_col_start_date + 1 #number of time steps
  
  range_dates <- seq.Date(from=as.Date(start_date),by="month",
                          length.out = nt )
  range_dates_str <- as.character(range_dates)

  if(version==1){
    #If version1, we are dropping the last date because it is often incomplete
    names_col <- c("g_id","sci_name","country",range_dates_str)
    df <- df[-1,-n_col] #remove the first row with incomplete header and last column with incomplete data
    names(df) <- names_col
    df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])
  }
  
  if(version==2){
    
    #if version2, we do not drop the last column and the data does not contain gid
    names_col <- c("sci_name","country",range_dates_str)
    df <- df[-1,] #remove the first row with incomplete header, keep the last column
    names(df) <- names_col
    df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])
    df[,n_col] <- sub('"',"",df[,n_col])
    
    #df[,62]
    #[997] "0.2015965068465905\""   "0.21044053064948223\""  "0\""                    "0.73523085537128918\"" 
  }

  ### Apply scaling factor
  
  out_filename_tmp <- sub(extension(infile_name),"",infile_name)#
  out_filename <- file.path(out_dir,paste0(out_filename_tmp,"_",out_suffix,".csv"))
  write.table(df,file=out_filename,sep=",")
  return(out_filename)
}

############################# END OF SCRIPT ########################