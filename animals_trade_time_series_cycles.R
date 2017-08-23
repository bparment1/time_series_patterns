############### SESYNC Research Support: Spatial demography ########## 
## Performing PCA on animals trade data at SESYNC.
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 08/20/2017
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
functions_time_series_cycles_analyses_script <- "time_series_cycles_analyses_functions_08202017.R" #PARAM 1

#script_path <- "C:/Users/edaut/Documents/gst_ts" #path to script #PARAM 2
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2

source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_cycles_analyses_script)) #source all functions used in this script 1.

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

#?na.omit
data_df <- na.omit(data_df)
dim(data_df)


df_ts <- (t(data_df))
dim(df_ts)
date_range <- c("2001.01.01","2010.12.31") #PARAM 15, NDVI Katrina
range_dates <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina
class(range_dates)

df_ts <- zoo(df_ts,range_dates)

(df_ts[1:10,])

plot(df_ts[,1])
dim(df_ts)


### 
nt <- 230
??fft
vect_z <- df_ts[,1]
test <- fft(vect_z)
test2 <- fft(vect_z,inverse=T) #not normalized

nt <- length(vect_z)

class(test)

test2

FF = abs(fft(vect_z)/sqrt(nt))^2
FF = abs(fft())
P = (4/128)*FF[1:65] # Only need the first (n/2)+1 values of the FFT result.
P = (4/128)*FF[1:65] # Only need the first (n/2)+1 values of the FFT result.

P=(4/nt)*FF[(nt/2)+1]


plot(FF)
f = (0:64)/128 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
nt_half <- nt/2
f = (0:nt_half)/nt
plot(f, P, type="l") # This plots the periodogram; type = “l” creates a line plot.  Note: l is lowercase L, not number 1.

ifft <- function(x) { fft(x, inverse=TRUE ) / length(x) }
tslm

#plot(data_df[1,])
xs <- seq(-2*pi,2*pi,pi/100)
wave.1 <- sin(3*xs)
wave.2 <- sin(10*xs)
par(mfrow = c(1, 2))
plot(xs,wave.1,type="l",ylim=c(-1,1)); abline(h=0,lty=3)
plot(xs,wave.2,type="l",ylim=c(-1,1)); abline(h=0,lty=3)

#https://www.r-bloggers.com/smoothing-techniques-using-basis-functions-fourier-basis/
#https://stats.stackexchange.com/questions/1207/period-detection-of-a-generic-time-series/1214#1214

findfrequency

?spec.ar

find.freq_test <- function(x){
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



#https://anomaly.io/detect-seasonality-using-fourier-transform-r/

# Install and import TSA package
install.packages("TSA")
library(TSA)

# read the Google Analaytics PageView report
raw = read.csv("20131120-20151110-google-analytics.csv")

# compute the Fourier Transform

harmonic_analysis_fft_run <- function(x){
  #
  #Function to compute 
  
  p <- periodogram(x)
  
  spectrum_val <- spectrum(as.numeric(x))
  n_orig <- length(x)
  freq_df <- data.frame(freq=p$freq, spec=p$spec)
  total_variance <- sum(freq_df$spec)
  freq_df$index <- as.numeric(row.names(freq_df))
  freq_df$period <- 1/freq_df$freq
  freq_df$period_orig <- n_used/freq_df$index
  freq_df$variance <- freq_df$freq/total_variance*100
  ranked_freq_df <- freq_df[order(-freq_df$spec),]
  
  #barplot(freq_df$variance)
  #p$orig.n* as.numeric(rownames(top2)
  n_used <- p$n.used
  
  ## Prepare return object:
  
  harmonic_fft_obj <- list(p,freq_df,ranked_freq_df,p$orig.n,p$n.used)
  names(harmonic_fft_obj) <- c("p","freq_df","ranked_freq_df","n_orig","n_used")
  
  return(harmonic_fft_obj)
}

#compute harmonic
#get nth harmonic to remove the period

# display the 2 highest "power" frequencies
top2
str(p)

extract_harmonic_fft_run <- function(x){
  
  x_trans <- fft(x) # transformed fft
  
  amp_val <- as.numeric(Mod(x_trans))
  
  amp <- amp_val
  amp[1] <- 0
  
  n <- length(x)
  
  n_half <- n/2
  plot(1:n_half,amp[1:n_half],"h")
  
  phase <- Arg(x_trans)  
  phase[10]
  barplot(phase)
  
  x_in <- 1:230
  amp[10]*sin(+ phase[10])
  
  
  ### Generate a sequence from sine
  #type_spatialstructure[5] <- "periodic_x1"
  a <- amp[10] #amplitude in this case
  b<- 0 # this should be the average!!
  T <- 230
  T <- 10
  phase <- 0
  x_input<-1:230
  
  ux <- sine_structure_fun(x_input,T,phase,a,b)
  plot(ux)
  
}


#The dominant peak area occurs somewhere around a frequency of 0.05.  Investigation of the periodogram values indicates that the peak occurs at nearly exactly this frequency.  This corresponds to a period of about 1/.05 = 20 time periods.  That’s 10 years, since this is semi-annual data.  
#Thus there appears to be a dominant periodicity of about 10 years in sunspot activity.

#sunspots=scan("sunspots.dat")
#plot(sunspots,type="b")
#x = diff(sunspots)
#I = abs(fft(x)/sqrt(458))^2
#P = (4/458)*I[1:230]
#freq = (0:229)/458
#plot(freq,P,type="l") 

sine_structure_fun <-function(x,T,phase,a,b){
  #Create sine for a one dimensional series
  #Note that sine function uses radian unit.
  #a=amplitude
  #b=mean or amplitude 0 of the series
  #T= stands for period definition
  #phase=phase angle (in radian!!)
  
  y <- a*sin((x*pi/T)+ phase) + b
}

### Generate a sequence from sine
#type_spatialstructure[5] <- "periodic_x1"
a<- 2 #amplitude in this case
b<- 0
T<- 230
phase <- 0
x_input<-1:230

ux <- sine_structure_fun(x_input,T,phase,a,b)
plot(ux)


# convert frequency to time periods
X <- fft(x)
fq <- 2 * pi /n
frq <- 0
FL <- 0
Fl[1] <- X[1]^2 / n*2

for( j in 2:(n/2)){
  FL[j] <- 2 * (X[j] )
  d
  
}
test[1]
test[1]^2 / n*2



#https://anomaly.io/seasonal-trend-decomposition-in-r/
  
############################## END OF SCRIPT #############################################

