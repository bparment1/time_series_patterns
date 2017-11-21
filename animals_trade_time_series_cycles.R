############### SESYNC Research Support: Animals Trade ############# 
## 
## Performing spectral analyses to uncover cycles on animals trade data at SESYNC.
## Different methods are available. Use multiple methods for multiple lines of evidence:
## Here is a quick lists of methods:
## 1) Spectrum: Uses spectrum generation with some detrending, the spectrum is smoothed
## 2) FFT periodogram: Uses raw outputs from FFT and plot squared of amplitude of frequencies
## 3) multitaper: Uses taper to estimate spectrum
## 4) harmonic regression (not implemented yet)
## 5) lag PCA/EOF (not implemented yet)
## 6) windowed Fourier T. (not implemented yet)
## 7) wavelet (not implemented yet)
## 
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 11/20/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: Animals Trade, Elizabeth Daut
## ISSUE: 
## TO DO: - add windowed/Short-Term Fourier transform option
##        - add wavelet option
##        - lag analyis with PCA to remove seasonality?
##        - PCA on spectral density matrix to identify and remove important harmonics
##        - implement harmonic option
##
## COMMIT: testing window generate on time series
##
## Links to investigate:
##
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
library(ggplot2)                             # Plotting functionalities
library(lubridate)                           # Dates manipulation functionalities
library(dplyr)                               # Data wrangling
library(forecast)                            # ARIMA and other time series methods
library(multitaper)                          # Multitaper estimation of spectrum
#library(GeneCycle)                          # Fisher test for harmonics and Time series functionalities: conflicts with periodogram function
library(xts)                                 # Extension for time series object and analyses
library(zoo)                                 # Time series object and analysis
library(mblm)                                # Theil Sen estimator
library(TSA)                                 # Time series analyses functionalities
library(Rwave)                               # Wavelet package R
library(pracma)                              # pracma
library(seewave)                             # time waves, time series functionalities
library(stats)
#detach(package:igraph) #conflict spectrum from this package with the general R stat psecturm

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
functions_time_series_cycles_analyses_script <- "time_series_cycles_analyses_functions_11202017.R" #PARAM 1


#script_path <- "C:/Users/edaut/Documents/gst_ts" #path to script #PARAM 2
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2

source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_cycles_analyses_script)) #source all functions used in this script 1.

#function_pca_eof <- "pca_eof_functions_08022017.R" #PARAM 1
#source(file.path(script_path,function_pca_eof)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
#in_dir <- "C:/Users/edaut/Documents/gst_ts"
in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data"

#ARGS 2
infile_name <- "vert_sp_gst_original_08162017.csv"
#infile_name_gst_totals <- "total_monthly_gst_averages.csv"
#use test data:
#infile_name <- "dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt" #use this data to test filtering

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
out_suffix <-"cycles_test_11202017" #output suffix for the files and ouptut folder #param 12

#ARGS_9
n_col_start_date <- 4
#ARGS 10
num_cores <- 2

#ARGS 11
selected_ts_infile <- "/nfs/bparmentier-data/Data/projects/animals_trade/data/selected_species_with_pattern_112018.csv"
range_window <- c("2012-01-01","2017-01-01")

################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

options(scipen=999)


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

#######################################
### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

#infile_name <- file.path(in_dir,infile_name)
#data_df <- read.table(infile_name,header=T,sep=",",stringsAsFactors = F)
#names(data_df)
#dim(data_df)
#data_df <- data_df[,1:230]

#?na.omit
#data_df <- na.omit(data_df)
#dim(data_df)

############### PART 1: Imported and time series transformation #####
## Step 1: read in the data and generate time stamps

infile_name <- "vert_sp_gst_original_08162017.csv"
data_ts_filename <- import_data_ts(infile_name = infile_name,
                                   in_dir = in_dir,
                                   scaling = scaling,
                                   n_col_start_date=4,
                                   start_date = start_date,
                                   end_date=NULL,
                                   out_dir = out_dir,
                                   out_suffix = out_suffix)

df_original <- read.table(data_ts_filename,sep=",",fill=T,header=T,check.names = F)
#object.size(df_original,units="Gb")

#Set up dates as column names  ###############  
#IS THIS NEEDED b/c ALREADY IN PROCESSING??? ###############

range_dates <- names(df_original)[n_col_start_date:ncol(df_original)]
range_dates_str <- as.character(range_dates)
range_dates<- ymd(range_dates_str)
 
#Transform and scale data

## get subset list in a text file?


if(is.null(selected_ts_infile)){
  df_subset <- df_original
}else{
  df_selected_ts_names <- read.table(selected_ts_infile,sep=",",header=T,check.names = F,stringsAsFactors = F)
  #df_subset <- subset(df_original,df_original$sci_name %in% df_selected_ts_names$sci_name)
  df_subset <- filter(df_original, sci_name %in% df_selected_ts_names$sci_name )
  df_subset <- filter(df_subset, country %in% df_selected_ts_names$country)
}

df_ts <- t(df_subset[n_col_start_date:ncol(df_subset)]) #transpose the data
df_ts <- as.data.frame(df_ts)
# #dim(df_ts)
df_ts_scaled <- df_ts[,]*scaling_factor  
 
#create ts zoo object
df_ts <- zoo(df_ts_scaled,range_dates)
class(df_ts)
dim(df_ts)
#combine country-species as column names
names_countries <- as.character(df_subset$country)
names_species <- as.character(df_subset$sci_name)
names_species <- sub(" ","_",names_species)
names_col <- paste(names_countries,names_species,sep="_")
names(df_ts) <- names_col
#View(df_ts)

#http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html


### Now let's remove the the most important components

#functions_time_series_cycles_analyses_script <- "time_series_cycles_analyses_functions_11072017b.R" #PARAM 1
#source(file.path(script_path,functions_time_series_cycles_analyses_script)) #source all functions used in this script 1.

j<- 10
x_ts1 <- df_ts[,j]
names(df_ts)
plot(x_ts1,type="l",main=names(df_ts)[j])
x_ts1 <- window(x_ts1,start=range_window[1],end=range_window[2])
plot(x_ts1,type="l",main=names(df_ts)[j])

freq_range <- c(11,13)
### Use the new function
test <- filter_freq(x_ts=x_ts1,
                    freq_range=freq_range,
                    w_length=NULL,
                    overlap_w=NULL)

plot(x_ts1,type="l",ylim=c(0,max(x_ts1)))
lines(as.numeric(test[,1]),type="l",col="red")
plot(test,type="l",col="red")

freq_range <- c(9,14)

undebug(filter_freq)
test <- filter_freq(x_ts=as.numeric(x_ts1),
                    freq_range=freq_range,
                    w_length=NULL,
                    overlap_w=80)
test <- filter_freq(x_ts=as.numeric(x_ts1),
                    freq_range=freq_range,
                    w_length=NULL,
                    overlap_w=80)
plot(test)
length(test)
plot(x_ts1)
plot(test)

y_range <- range(c(test,x_ts1))
plot(x_ts1,type="l",ylim=y_range)
lines(test,type="l",col="red")


x_ts_filtered <- filter_frequency_and_generate_harmonics(x_ts=x_ts1,freq_range=freq_range)


selected_period <- 12

undebug(filter_frequency_and_generate_harmonics)


x_ts1_filtered <- filter_frequency_and_generate_harmonics(x_ts1,
                                                         freq_range=NULL,
                                                         selected_period=selected_period,
                                                         variance_threshold=NULL,
                                                         peak_opt=NULL)

plot(x_ts3,type="l")
lines(x_ts_filtered,col="red")

coef_fft_df <- extract_harmonic_fft_parameters_run(as.numeric(x_ts1[1:60]))
View(coef_fft_df)
plot(coef_fft_df$amplitude[-1],type="h",xlab="harmonic") #highest amplitude is 23

test<- fft(x_ts1)
time_val <- 1:length(x_ts1)
mod <- lm(x_ts1~time_val)

vect_z <- residuals(mod)
plot(vectz)
coef_fft_df <- extract_harmonic_fft_parameters_run(as.numeric(vect_z))
View(coef_fft_df)
plot(coef_fft_df$amplitude[-1],type="h",xlab="harmonic") #highest amplitude is 23

freq_range <- c(4,6)
debug(filter_frequency_and_generate_harmonics)
x_ts_filtered <- filter_frequency_and_generate_harmonics(x_ts=x_ts1,freq_range=freq_range)
plot(x_ts1,ylim=c(-10000,100000))


zoo(x_ts1_filtered,
lines(c(0,x_ts_filtered,0),col="red")

plot(x_ts_filtered)

############################## END OF SCRIPT #############################################

#https://www.r-bloggers.com/fir-filter-design-and-digital-signal-processing-in-r/
  
#https://dsp.stackexchange.com/questions/22637/implementing-a-bandpass-filter-in-r-this-codes-logic-confuses-me

#https://math.stackexchange.com/questions/1002/fourier-transform-for-dummies

# https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/findpeaks
# 
# x <- seq(0, 1, len = 1024)
# pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
# hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
# wdt <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
# 
# pSignal <- numeric(length(x))
# 
# for (i in seq(along=pos)) {
#   pSignal <- pSignal + hgt[i]/(1 + abs((x - pos[i])/wdt[i]))^4
# }
# plot(pSignal,type="l")
# 
# findpeaks(pSignal, npeaks=3, threshold=4, sortstr=TRUE)
# findpeaks(pSignal, sortstr=TRUE)
# 
# pSignal[1:798] <- 0
# p
# #https://anomaly.io/seasonal-trend-decomposition-in-r/
# 
# ### Generate function for formal test of presence of periodicity/frequency in 
# ## the data.
# #https://en.wikipedia.org/wiki/Welch%27s_method
# #The periodogram is a biased estimate in most cases so we need
# #to estimate the spectrum with other methods or try to get
# #a better estimate. This often results in lower resolution identification
# #of the frequency but improved estimate.
# #https://stats.stackexchange.com/questions/12164/testing-significance-of-peaks-in-spectral-density
# 
# 
# # convert frequency to time periods
# X <- fft(x_ts1)
# fq <- 2 * pi /nt
# frq <- 0
# FL <- 0
# Fl[1] <- X[1]^2 / nt*2
# 
# for( j in 2:(n/2)){
#   FL[j] <- 2 * (X[j] )
#   d
#   
# }
# test[1]
# test[1]^2 / n*2
# 
# #The dominant peak area occurs somewhere around a frequency of 0.05.  Investigation of the periodogram values indicates that the peak occurs at nearly exactly this frequency.  This corresponds to a period of about 1/.05 = 20 time periods.  That’s 10 years, since this is semi-annual data.  
# #Thus there appears to be a dominant periodicity of about 10 years in sunspot activity.
# 
# #sunspots=scan("sunspots.dat")
# #plot(sunspots,type="b")
# #x = diff(sunspots)
# #I = abs(fft(x)/sqrt(458))^2
# #P = (4/458)*I[1:230]
# #freq = (0:229)/458
# #plot(freq,P,type="l") 
# 
# test2
# 
# FF = abs(fft(vect_z)/sqrt(nt))^2
# FF = abs(fft())
# P = (4/128)*FF[1:65] # Only need the first (n/2)+1 values of the FFT result.
# P = (4/128)*FF[1:65] # Only need the first (n/2)+1 values of the FFT result.
# 
# P=(4/nt)*FF[(nt/2)+1]
# 
# 
# plot(FF)
# f = (0:64)/128 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
# nt_half <- nt/2
# f = (0:nt_half)/nt
# plot(f, P, type="l") # This plots the periodogram; type = “l” creates a line plot.  Note: l is lowercase L, not number 1.
# 
# ifft <- function(x) { fft(x, inverse=TRUE ) / length(x) }
# tslm
# 
# #plot(data_df[1,])
# xs <- seq(-2*pi,2*pi,pi/100)
# wave.1 <- sin(3*xs)
# wave.2 <- sin(10*xs)
# par(mfrow = c(1, 2))
# plot(xs,wave.1,type="l",ylim=c(-1,1)); abline(h=0,lty=3)
# plot(xs,wave.2,type="l",ylim=c(-1,1)); abline(h=0,lty=3)
# 
# #https://www.r-bloggers.com/smoothing-techniques-using-basis-functions-fourier-basis/
# #https://stats.stackexchange.com/questions/1207/period-detection-of-a-generic-time-series/1214#1214
# 
# findfrequency
# 
# ?spec.ar
# 
# find.freq_test <- function(x){
#   n <- length(x)
#   spec <- spec.ar(c(x),plot=FALSE)
#   
#   if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
#   {
#     period <- round(1/spec$freq[which.max(spec$spec)])
#     if(period==Inf) # Find next local maximum
#     {
#       j <- which(diff(spec$spec)>0)
#       if(length(j)>0)
#       {
#         nextmax <- j[1] + which.max(spec$spec[j[1]:500])
#         period <- round(1/spec$freq[nextmax])
#       }
#       else
#         period <- 1
#     }
#   }
#   else
#     period <- 1
#   return(period)
# }
# 
# 
# #https://anomaly.io/detect-seasonality-using-fourier-transform-r/
# 
# # Install and import TSA package
# install.packages("TSA")
# library(TSA)
# 
# # read the Google Analaytics PageView report
# raw = read.csv("20131120-20151110-google-analytics.csv")
# 
# #compute harmonic
# #get nth harmonic to remove the period
# 
# # display the 2 highest "power" frequencies
# top2
# str(p)
# 
# ####
# P=abs(2*fft(x_ts1)/230)^2
# P[10]
# #see TSA book p.179
# 
# y<- x_ts1_lm
# tappercent=.05
# N= length(y)
# fn = 1/(2*dt)
# tapy = spec.taper(y, p=tappercent)
# ##tapy = tapy-mean(tapy)
# plot(tapy,type="l")
# Y = fft(tapy)
# Pyy = (Mod(Y)^2)/(N*N)
# ## Pyy = Y * Conj(Y)
# n = floor(length(Pyy)/2)
# Syy = Pyy[1:n]
# 
# fs = (0:(length(Syy)-1))/length(Syy)
# 
# plot(fs, Syy, type='l', xlab="frequency",
#      ylab="Power Density", log='')
# 
# fs = (0:(length(Syy)-1))*fn/length(Syy)
# plot(fs, Syy, type='l', xlab="frequency",
#      + ylab="Power Density", log='')
# 
# #testing for significant harmonics
# #library(GeneCycle)
# ?fisher.g.test

#Examples


#https://dsp.stackexchange.com/questions/6220/why-is-it-a-bad-idea-to-filter-by-zeroing-out-fft-bins
