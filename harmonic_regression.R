############################## Harmonic regression #################### 
##
## Functions generetaed through various research projects and SESYNC research support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 05/08/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: General use script
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

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

#Benoit setup
script_path <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/scripts"

crop_data_processing_functions <- "harmonic_regression_functions_05082019b.R"
source(file.path(script_path,crop_data_processing_functions))

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/data"

#ARGS 2
infile_name_df <- "dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt" #use this data to test filtering
infile_name_raster <- "reg2_NDVI_katrina.tif"
#ARGS 3
#start_date <- "2004-01-01"
start_date <- "2012-11-01"  #new data starts in November 2012
#ARGS 4
end_date <- NULL
#ARGS 5
out_dir <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/outputs"
#ARGS 6
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 7
out_suffix <-"example_ts_04242019" #output suffix for the files and ouptut folder #param 12
#ARGS 8
num_cores <- 2 # number of cores

#range_window <- c("2012-01-01","2017-01-01")

################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

options(scipen=999)

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #xcan modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#######################################
### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

infile_name_df <- file.path(in_dir,infile_name_df)
#
data_df <- read.table(infile_name_df,header=T,sep=",",stringsAsFactors = F)
#names(data_df)
#start_date <- "2004-01-01"
#start_date <- "2012-11-01"  #new data starts in November 2012

#y ~ A0 + b1 cos(x) + b2* sin(x)
#y ~ b0 + b1*x1 + b2*x2

y_all <- as.numeric(data_df[1400,1:230])
y_all

plot(y_all)
plot(y_all[1:23])
y <- y_all[1:24]
n <- length(y)


#debug(harmonic_regression)
harmonic_results <- harmonic_regression(y,n,
                                        harmonic_val=NULL,
                                        mod_obj=F,
                                        figure=F)

#View(harmonic_results)

####################
#### This is synthetic value

n <- 24
x <- seq(1, 24)
p <- 1 #harmonic 1
omega= 2*pi*p/n

y <- 2*cos(omega*x) + rnorm(n, sd=0.2)
# y_clean <- sin(2*x + 5)
plot(y)

harmonic_results2 <- harmonic_regression(y,n,
                                        harmonic_val=NULL,
                                        mod_obj=T,
                                        figure=F)

harmonic_results2$l_harmonic_obj
mod <- harmonic_results2$l_harmonic_obj[[1]]$mod

summary(mod)
plot(y)
lines(mod$fitted.values)

barplot(Im(fft(y)))
plot(Re(fft(y)))

#### generate function
#### split in annual windows

x <- y_all

#debug(split_sequence)
split_obj <- split_sequence(y_all,n=23)
split_obj$list_y[[9]]
length(split_obj$list_y[[9]])

harmonic_results3 <- harmonic_regression(split_obj$list_y[[9]],n=23+1,
                                         harmonic_val=NULL,
                                         mod_obj=T,figure=F)

mod <- harmonic_results3$l_harmonic_obj[[1]]$mod

summary(mod)
plot(y)
lines(mod$fitted.values)

#### Testing with raster time series and run across multiple time

infile_name_raster <- file.path(in_dir,infile_name_raster)
#
#data_df <- read.table(infile_name,header=T,sep=",",stringsAsFactors = F)
r <- brick(infile_name_raster)
names(r)
plot(r,y=1)

NAvalue(r)
plot(r,y=1,colNA="black")


### harmonic 1 amplitude for first year

x <- 1:230

#debug(split_sequence)

split_obj <- split_sequence(x,n=23)
intervals_df <- split_obj$intervals
#lengh(split_obj$list_y[[1]])
split_obj$list_y[[1]]
n_split <- nrow(intervals_df)

## make this a function later

### for A1, function seems to work for A1 and A2
l_r_A1 <- vector("list",length=n_split)
n_val <- 23
harmonic_val <- 1
var_name <- "A"

for(i in 1:n_split){
  start_val <- intervals_df$start[i]
  end_val <- intervals_df$end[i]
  #p <- calc(subset(r,1:23), fun=harmonic_reg_f1)
  #r_out <- try(calc(subset(r,end_val:start_val), fun=harmonic_reg_f1))
  #multiple arg does not work
  #p <- calc(subset(r,end_val:start_val), fun=harmonic_reg_f1,
  #          n=24,harmonic=1)
  n_val <- end_val - start_val + 1
  #test <- calc(subset(r,end_val:start_val), 
  #          fun=function(y){harmonic_reg_raster(y,n=n_val,harmonic=harmonic_val)})
  #debug(harmonic_reg_raster)
  #harmonic_reg_raster(y,n=n_val,harmonic=harmonic_val,var_name=)
  
  r_out <- try(calc(subset(r,end_val:start_val), 
               fun=function(y){harmonic_reg_raster(y,
                                                   var_name="A",
                                                   n=n_val,
                                                   harmonic=NULL)}))
  names(r_out)
  l_r_A1[[i]] <- r_out
  #plot(r_out)
  #rm(p)
}

r_A1 <- stack(l_r_A1)
plot(r_A1)
plot(l_r_A1[[5]])
plot(l_r_A1[[6]])

l_r_A1
plot(r_A1,y=1:2)

################################### End of script #######################################

