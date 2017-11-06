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
## DATE MODIFIED: 11/06/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: Animals Trade, Elizabeth Daut
## ISSUE: 
## TO DO: - add windowed Fourier transform option
##        - add wavelet option
##        - lag analyis with PCA to remove seasonality?
##        - PCA on spectral density matrix to identify and remove important harmonics
##
## COMMIT: testing filtering on real and synthetic data with known cycles
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

#functions_time_series_analyses_script <- "time_series_functions_08012017.R" #PARAM 1
#functions_processing_data_script <- "processing_data_google_search_time_series_functions_07202017.R" #PARAM 1
functions_time_series_cycles_analyses_script <- "time_series_cycles_analyses_functions_11022017.R" #PARAM 1

#script_path <- "C:/Users/edaut/Documents/gst_ts" #path to script #PARAM 2
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2

#source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
#source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_cycles_analyses_script)) #source all functions used in this script 1.

#function_pca_eof <- "pca_eof_functions_08022017.R" #PARAM 1
#source(file.path(script_path,function_pca_eof)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
#in_dir <- "C:/Users/edaut/Documents/gst_ts"
in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data"

#ARGS 2
#infile_name <- "vert_sp_gst_original_08162017.csv"
#infile_name_gst_totals <- "total_monthly_gst_averages.csv"
#use test data:
infile_name <- "dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt" #use this data to test filtering

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
out_suffix <-"cycles_test_11062017" #output suffix for the files and ouptut folder #param 12

#ARGS_9
n_col_start_date <- 4
#ARGS 10
num_cores <- 2

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

infile_name <- file.path(in_dir,infile_name)
data_df <- read.table(infile_name,header=T,sep=",",stringsAsFactors = F)
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
#
(df_ts[1:10,])

plot(df_ts[,1])
dim(df_ts)

### 
nt <- 230

#??fft
vect_z <- df_ts[,1]
test <- fft(vect_z)
plot(as.complex(test),type="p")
plot(Real(test))
plot(Im(test))

class(test) # zoo object
### See part 3 for more in depth analyses and removal of frequencies for 

######################################################################
################ PART 2: Generate Synthetic Data Time Series and Run tests ##########
######### Test 1: using time step from 800

## For the test use 8000
### Generate a sequence from sine with 8000 steps with:
##  - periods of 800  and Amplitude 2
##  - periods of 1600 and Amplitude 1
##  - mean average of 0
##  - no phase delay i.e. 0
##  - random noise 0 0.01

amp<- c(2,1) #amplitude in this case
b<- 0
T<- c(800,1600)
phase_val <- 0

nt <- 8000  #change length for test with seewave ffilter
temp_periods <- c(T)
temp_period_quadrature <- NULL #no quadrature
random_component <- c(0,0.1)

list_param <- list(nt,phase_val,temp_periods,amp,
                   temp_period_quadrature,
                   random_component)

names(list_param) <- c("nt","phase","temp_periods",
                       "amp","temp_period_quadrature",
                       "random_component")

ts_synthetic_8000 <- adding_temporal_structure(list_param)
names(ts_synthetic_8000)

### Using patterns, generate synthetic time series
x_ts1_8000 <- ts_synthetic_8000$t_period_800 + ts_synthetic_8000$trend + ts_synthetic_8000$norm
x_ts2_8000 <- ts_synthetic_8000$t_period_800 + ts_synthetic_8000$t_period_1600 
x_ts3_8000 <- ts_synthetic_8000$t_period_800 + ts_synthetic_8000$t_period_1600 + ts_synthetic_8000$norm
x_ts4_8000 <- ts_synthetic_8000$t_period_800/2 + ts_synthetic_8000$t_period_1600

plot(x_ts2_8000,type="l") #peaks for period 1600 and 800
spectrum(x_ts2_8000)

periodogram(x_ts2_8000)

## Test to filter out periods/frequencies
## Use default filter window: Hanning
#ffilter(as.matrix(x_ts2),f=230,from=40,to=45)
x_ts_filtered <- ffilter(as.matrix(x_ts3_8000),f=8000,from=0.18,to=0.2)
plot(x_ts3_8000,col="red")
lines(x_ts_filtered) #mostly filtered out!!!

x_ts_filtered <- ffilter(as.matrix(x_ts3_8000),f=8000,from=100,to=2000)
#test<- ffilter(as.matrix(x_ts2),f=8000,from=0.18,to=22)

plot(x_ts3_8000,col="red")
lines(x_ts_filtered) #mostly filtered out!!!
periodogram(x_ts_filtered) #still peak but if with noise might not appear
periodogram(x_ts3_8000,xlim=c(0,0.02)) #zoom in
periodogram(x_ts_filtered,xlim=c(0,0.02)) #still peak but if with noise might not appear

## dealt with the conflicts in namespace
spectrum(x_ts_filtered)
spectrum(x_ts3_8000)
spectrum(x_ts3_8000,xlim=c(0,0.2))
spectrum(x_ts_filtered,xlim=c(0,0.2))

spectrum(x_ts3_8000,xlim=c(0,0.2))
spectrum(x_ts_filtered,xlim=c(0,0.2))

################# Test 2 with artificial data #######
### Generate a sequence from sine: 230 steps

## For the test use 230
### Generate a sequence from sine with 8000 steps with:
##  - periods of 23  and Amplitude 2
##  - periods of 46 and Amplitude 1
##  - mean average of 0
##  - no phase delay i.e. 0
##  - random noise 0 0.01


#type_spatialstructure[5] <- "periodic_x1"
amp<- c(2,1) #amplitude in this case
b<- 0
T<- c(23,46) #annual cycle of 23

phase_val <- 0

nt <- 230  #change lenght for test with seewave ffilter
temp_periods <- c(T)
temp_period_quadrature <- NULL
random_component <- c(0,0.1)

list_param <- list(nt,phase_val,temp_periods,amp,
                   temp_period_quadrature,
                   random_component)

names(list_param) <- c("nt","phase","temp_periods",
                       "amp","temp_period_quadrature",
                       "random_component")
#nt <- list_param$nt
#phase <- list_param$phase
#temp_periods <- list_param$temp_periods
#amp <- list_param$amp
#temp_period_quadrature <- list_param$temp_period_quadrature
#random_component <- list_param$random_component #mean and sd used in rnorm

#debug(adding_temporal_structure)
ts_synthetic <- adding_temporal_structure(list_param)
names(ts_synthetic)  

#y <- a*sin((x_input*2*pi/T)+ phase_val) + b
#plot(y)
#ux <- sine_structure_fun(x_input,T,phase_val,a,b)
#plot(ux)

x_ts1 <- ts_synthetic$t_period_23 + ts_synthetic$trend + ts_synthetic$norm
x_ts2 <- ts_synthetic$t_period_23 + ts_synthetic$t_period_46
x_ts3 <- ts_synthetic$t_period_23 + ts_synthetic$t_period_46 + ts_synthetic$norm
x_ts4 <- ts_synthetic$t_period_23/2 + ts_synthetic$t_period_46
x_ts5 <- ts_synthetic$t_period_23

plot(x_ts1,type="l")
plot(x_ts5,type="l")
x_fft <- fft(x_ts5)
plot(Real(x_fft), type="l")
mid_val<- length(Real(x_fft))/2
abline(v=mid_val)
plot(Im(x_fft), type="l")
abline(v=mid_val)

plot(Real(x_fft), type="h")
plot(Im(x_fft), type="h")

test <- Im(x_fft)
test <- Real(x_fft)
plot(Real(x_fft), type="h")

test2<- fft(test)
plot(test2,type="l")
mesh_val<- meshgrid(Real(test2),Im(test2))

plot(Real(test2),type="l")
plot(Im(test2),type="l")

#x_ts1 <- test$t_period_23 + test$trend + test$unif + test$norm

plot(ts_synthetic$t_period_23,type="l")
plot(ts_synthetic$trend,type="l")
plot(ts_synthetic$unif)

#Error in seq.default(1, n - wl, wl - (ovlp * wl/100)) : 
#  wrong sign in 'by' argument

plot(x_ts3,type="l")
lines(x_ts2,type="l",col="red")
lines(x_ts4,type="l",col="green")

### Length of the window must be set otherwise it is set to 1024
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=23) #this works
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=46) #this works
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=w_length) #this works
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=230) #this does not work
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=115) #this works but note that data is cut off!!!!

0.128*8000 # use this ratio?
w_length<- 0.128*230
w_length
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=w_length,ovlp=90) #this works
x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=0.18,to=0.2,wl=w_length,ovlp=0) #this works

length(x_ts_filtered)

#?ffilter
plot(x_ts3,type="l")
lines(x_ts_filtered)
#> x_ts3 - x_ts_filtered
#Error: dims [product 228] do not match the length of object [230]
#In addition: Warning message:
# In x_ts3 - x_ts_filtered :
# longer object length is not a multiple of shorter object length
x_ts_obj <- ts(x_ts3,start=c(2001,1),end=c(2010,23),frequency = 23)
plot(x_ts_obj)

plot(x_ts3,type="l")
test_decomp <- decompose(x_ts_obj,filter=23)
plot(test_decomp$seasonal)
test <- x_ts_obj -test_decomp$seasonal
plot(test)

### Maybe test fftfilt from signal package?
### Also filter

hanning_filt<- hanning.w(23)
hanning_filt<- hanning.w(230)

plot(hanning_filt,type="l")
periodogram(hanning_filt)
spectrum(hanning_filt)
fft_x<-fft(hanning_filt)
plot(fft_x)
plot(Real(fft_x),type="l")
plot(Im(fft_x),type="l")

#undebug(ffilter)
#x_ts_filtered <- ffilter(as.matrix(x_ts3),f=230,from=20,to=50,wl=23) #this does not work

#debug(ffilter)
#test<- ffilter(as.matrix(x_ts2),f=8000,from=0.18,to=22)

plot(x_ts3,col="red",type="l")
lines(x_ts_filtered) #mostly filtered out!!!
periodogram(x_ts_filtered) #still peak but if with noise might not appear
periodogram(x_ts3)
periodogram(x_ts3,xlim=c(0,0.02)) #zoom in
periodogram(x_ts_filtered,xlim=c(0,0.02)) #still peak but if with noise might not appear
spectrum(x_ts_filtered)
spectrum(x_ts3)
spectrum(x_ts3,xlim=c(0,0.2))
spectrum(x_ts_filtered,xlim=c(0,0.2))

periodogram(x_ts2)

#### Now find out if you can see the cycles

nt <- 230
time_steps <- 1:nt
val_df <- data.frame(x_ts1,time_steps)

mod <- lm(x_ts1 ~ time_steps,val_df)
plot(mod$residuals,type="l")
x_ts1_lm <- mod$residuals

periodogram(x_ts1) #need to remove trends or will impact low frequencies
periodogram(x_ts1_lm)
periodogram(x_ts2)

spec_obj <- spectrum(x_ts1,fast=F)
spec_obj <- spectrum(x_ts1_lm,fast=F)
spec_obj <- spectrum(x_ts2,fast=F)

findpeaks(spec_obj$spec)

p <-periodogram(x_ts1_lm,fast=F) #spectral leakage!!

## Fix this
#debug(spectrum_analysis_fft_run)
test <- spectrum_analysis_fft_run(x_ts1_lm)

### Confidence interval 
#The simplest is to compare to the null model of whtie noise. White noise has a horizontal spectrum line because
#variance is not concentrated in particular frequencies.
#However positive autocorrelation can skew amplitdues towards low frequencies so a null test against white noise
#can be problematic.
#other option can be to use another null model comparing to an autoregressive process.

x_ts1_diff <- diff(x_ts1)
plot(x_ts1,type="l")
plot(x_ts1_diff,type="l") #loosing one data point, also note that this affected the amplitude too!!!

#undebug(harmonic_analysis_fft_run)

### find harmonic cycles
undebug(spectrum_analysis_fft_run)
spectrum_analysis_fft_obj_diff <- spectrum_analysis_fft_run(x_ts1_diff)
spectrum_analysis_fft_obj_lm <- spectrum_analysis_fft_run(x_ts1_lm)

spectrum_analysis_fft_obj_ts2 <- spectrum_analysis_fft_run(x_ts2)

### Get fft coef for harmonics:
#undebug(extract_harmonic_fft_parameters_run)
coef_fft_obj_ts2 <- extract_harmonic_fft_parameters_run(x_ts2)
coef_fft_obj_ts3 <- extract_harmonic_fft_parameters_run(x_ts3)

View(coef_fft_obj_ts3)


### Wavelet can detect frequencies that are localized in time.
#Wavelet can then detect changes in dominant frequencies across a time series.
#It can be used for non stationary time series
wavelet_run(x_ts3)
wavelet_run(x_ts1_lm)
wavelet_run(x_ts4) #same power for period 23 and 46

#Check other wavelets options!!
#http://jaysthesisblog.com/R/wavelets_intro.html

### Now let's remove the the most important components

functions_time_series_cycles_analyses_script <- "time_series_cycles_analyses_functions_11062017.R" #PARAM 1
source(file.path(script_path,functions_time_series_cycles_analyses_script)) #source all functions used in this script 1.


plot(x_ts3,type="l")
freq_range <- c(20,50)
### Use the new function
test <- filter_freq(x_ts=x_ts3,
                    freq_range=freq_range,
                    w_length=NULL,
                    overlap_w=NULL)
  
test <- filter_freq(x_ts=x_ts3,
                    freq_range=freq_range,
                    w_length=NULL,
                    overlap_w=90)


debug(filter_frequency_and_generate_harmonics)

test <- filter_frequency_and_generate_harmonics(x_ts=x_ts3,freq_range=freq_range)

  
  
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
