############### SESYNC Research Support: Urbanization Impact on Biodiversity ########## 
##
## 
## DATE CREATED: 05/30/2017
## DATE MODIFIED: 05/30/2017
## AUTHORS: Benoit Parmentier 
## Version: 1
## PROJECT: Urbanization impact on biodiversity
## ISSUE: 
## TO DO:
##
## COMMIT: initial commit
##
## Links to investigate:
#https://gis.stackexchange.com/questions/119993/convert-line-shapefile-to-raster-value-total-length-of-lines-within-cell
#https://edzer.github.io/sfr/articles/sfr.html
#https://geographicdatascience.com/2017/01/06/first-impressions-from-sf-the-simple-features-r-package/

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

###### Functions used in this script


import_data_ts <- function(infile_name,in_dir=".",scaling=1,n_col_start_date=4,start_date="2004-01-01",end_date=NULL,out_dir=".", out_suffix=""){
  #Import data sent by google and format columns for future processing and visualization.
  
  ##### BEGIN FUNCTION ####
  
  ### Add quote="" otherwise EOF warning and error in reading
  df <- read.table(file.path(in_dir,infile_name),sep=",",fill=T,quote="",header=F)
  
  #remove first row that contains names
  #rename columns using:
  #n_col_start_date <- 3 #may need to be changed
  n_col <- ncol(df)
  nt <- n_col - n_col_start_date #we are dropping the last date because it is often incomplete
  
  range_dates <- seq.Date(from=as.Date(start_date),by="month",length.out = nt )
  range_dates_str <- as.character(range_dates)
  #class(range_dates)
  #start_date
  names_col <- c("g_id","sci_name","country",range_dates_str)
  df <- df[-1,-n_col] #remove the first row with incomplete header and last column with incomplete data
  names(df) <- names_col
  df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])
  
  out_filename_tmp <- sub(extension(infile_name),"",infile_name)#
  out_filename <- file.path(out_dir,paste0(out_filename_tmp,"_",out_suffix,".csv"))
  write.table(df,filename=out_filename)
  return(out_filename)
}
