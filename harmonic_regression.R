############################## Harmonic regression #################### 
##
## Research Support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 04/17/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis Managin Hurricanes
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

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

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/data"

#ARGS 2
infile_name <- "dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt" #use this data to test filtering
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
out_suffix <-"example_ts_04162019" #output suffix for the files and ouptut folder #param 12
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
#start_date <- "2004-01-01"
start_date <- "2012-11-01"  #new data starts in November 2012

#y ~ A0 + b1 cos(x) + b2* sin(x)
#y ~ b0 + b1*x1 + b2*x2

y_all <- as.numeric(data_df[1400,1:230])
y_all

plot(y_all)
plot(y_all[1:24])
y <- y_all[1:24]
n <- length(y)

harmonic_regression<- function(y,n,harmonic_val){
  
  ## Default to two first harmonic:
  if(is.null(harmonic_val){
    p <- 1:2
  }
  #harmonic_val = 1 # pr from 1 to n/2
  
  #n<-24
  
  omega = lappply(p,function(p){2*pi*p/n})
  
  t <- 1:n
  cos_val =cos(omega*t)
  sin_val =sin(omega*t)
  plot(cos_val)
  plot(sin_val)
  
  in_df <- data.frame(y=y,cos_val=cos_val,sin_val=sin_val)
  mod <- lm(y~cos_val + sin_val,data=in_df)
  summary(mod)
  
  plot(y)
  lines(mod$fitted.values)
  points(mod$fitted.values)
  
  plot(mod$fitted.values,ylim=c(3000,10000))
  points(y,pch=2)
  
  plot(mod$fitted.values,ylim=c(4000,8000))
  points(y,pch=2)
  
  mod$coefficients[2]
  plot(y)
  return()
}
#first harmonic

p = 1 #from 1 to n/2

n<-24
omega= 2*pi*p/n

t <- 1:n
cos_val =cos(omega*t)
sin_val =sin(omega*t)
plot(cos_val)
plot(sin_val)

in_df <- data.frame(y=y,cos_val=cos_val,sin_val=sin_val)
mod <- lm(y~cos_val + sin_val,data=in_df)
summary(mod)

plot(y)
lines(mod$fitted.values)
points(mod$fitted.values)

plot(mod$fitted.values,ylim=c(3000,10000))
points(y,pch=2)

plot(mod$fitted.values,ylim=c(4000,8000))
points(y,pch=2)

mod$coefficients[2]
plot(y)

#amplitude
A = atan(-mod$coefficients[3]/mod$coefficients[2])

A

### FFT, does not work with missing values
barplot(Im(fft(y)))
plot(Re(fft(y)))

####################
#### This is synthetic value

x <- seq(1, 23)
p <-1
n=23
omega= 2*pi*p/n

y <- 2*cos(omega*x) + rnorm(23, sd=0.2)
# y_clean <- sin(2*x + 5)
plot(y)

t <- 1:n
cos_val =cos(omega*t)
sin_val =sin(omega*t)
plot(cos_val)
plot(sin_val)

in_df <- data.frame(y=y,cos_val=cos_val,sin_val=sin_val)
mod <- lm(y~cos_val + sin_val,data=in_df)
summary(mod)
plot(y)
lines(mod$fitted.values)

barplot(Im(fft(y)))
plot(Re(fft(y)))

#### generate function
#### split in annual windows

################################### End of script #######################################

# infile_name <- "dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt" #use this data to test filtering
# 
# 
# #https://stackoverflow.com/questions/34444193/fitting-harmonic-model-in-r-of-the-form-asinbxc
# #As you dont provide any data, here some simulated data:
# 
# set.seed(1)
# x <- seq(0, 10, length = 100)
# y <- sin(2*x + 5) + rnorm(100, sd=0.5)
# y_clean <- sin(2*x + 5)
# 
# # Plot of the noised data (y) and the noiseless line (y_clean) which you want to approximate through your model:
# plot(x, y, ylim = c(-2,3))
# lines(x,y_clean, col = "green")
# 
# 
# # Model estimation and adding of the fitted values to the previous plot:
# model <- nls(y~a*sin(b*x+c), start = list(a=1,b=1.5, c=1),control = list(maxiter = 500))
# 
# lines(x, fitted(model), col = "red")
# 
# legend("topright", col=c("green", "red"), legend = c("\"true\" curve", "fitted line"), bty = "n", lty = 1)
# 
# library(TSA)
# data(tempdub)
# # first creates the first pair of harmonic functions and then fit the model
# har.=harmonic(tempdub,1)
# model4=lm(tempdub~har.)
# summary(model4)
# plot(model4)
# 
# plot(har.[,1])
# plot(har.[,1],type="l")
# lines(har.[,2],type="l",col="red")
# ?harmonic
# plot(tempdub)
# length(tempdub)
# 
# data(tempdub)
# month.=season(tempdub) # the period sign is included to make the printout from
# # the commands two line below clearer; ditto below.
# model2=lm(tempdub~month.-1) # -1 removes the intercept term
# summary(model2)
# 
# test <- tempdub
# test[c(33,45,102,111,137)] <- NA
# test
# har.=harmonic(test,1)
# plot(har.[,1],type="l")
# 
# test2 <- harmonic(test,m=2)
# test2
# dim(test2)
# class(test2)
# 
# plot(test2[,1],type="l")
# lines(test2[,2],type="l",col="pink")
# lines(test2[,3],type="l",col="green")
# lines(test2[,2],type="l",col="pink")
# 
# 
# plot(test2[,1],type="l")
# lines(test2[,2],type="l",col="pink")
# lines(test2[,3],type="l",col="green")
# lines(test2[,2],type="l",col="pink")
# 
# harmonic
# as.vector(test)
# 
# fft(as.vector(test))
# 
# ###
# #A0/2 + Sum(a* cost(rt + b* sin(rt)))
# 
# #first harmonic
# 
# #https://stat.ethz.ch/pipermail/r-help/2008-May/162879.html
# 
# y ~ a + c*sin(x+b)
# 
# #so the amplitude of the sine wave is adjustable (otherwise, you assume (or
# #                                                                        know) that the amplitude is 1).  Then
# 
# y ~ a + c*sin(b)*cos(x) + c*cos(b)*sin(x)
# y ~ A0 + b1 cos(x) + b2* sin(x)
# 
# #follicles ~ sin(2*pi*Time)+cos(2*pi*Time)
# #
# #or
# 
# y ~ b0 + b1*x1 + b2*x2
# 
# #which is a linear regression form that you can do using the lm function.
# #After you get b0, b1 and b2 you do
# 
# #a = b0
# #b1^2 + b2^2 = c^2*(sin^2(b) + cos^2(b)) = c^2  ====> c = sqrt(b1^2 + b2^2)
# 
# #b1/b2 = tan(b)  ====>  b = arctan(b1/b2)
# 
# 
# x <- 1:12
# phase <- 0
# a<- 1
# b <- 0
# T <- 12 #if
# y <- a*sin((2*x*pi/T)+ phase) + b
# 
# plot(y,type="l")
# plot(y~x)
# 
# x <- 1:12
# phase <- 0
# a<- 1
# b <- 0
# T <- 6 #if
# y <- a*sin((x*2*pi/T)+ phase) + b
# 
# y
# plot(y,type="l")
# plot(y~x)
# 
# #Generate spatial pattern 5:     
# n_col <- ncol(r_init)
# n_row <- nrow(r_init)
# 
# #u <- xFromCol(r_init,col=1:n_col)
# #add padding option later...buffer from a specific distance and tailling of at 0.1
# u <- 1:n_col
# a<- 1 #amplitude in this case
# b<- 0
# T<- n_col
# phase <- 0
# use_cos <- FALSE
# ux <- sine_structure_fun(u,T,phase,a,b,use_cos)
# ux_rep <-rep(ux,time=n_row)  
# r1 <-setValues(r_init,ux_rep)  #note efficient in memory might need to revise this
# #plot(r)
# 
# v <- 1:n_row
# a<- 1 #amplitude in this case
# b<- 0
# T<- n_row
# phase <- 0
# use_cos <- FALSE
# vx <- sine_structure_fun(v,T,phase,a,b,use_cos)
# vx_rep <- unlist((lapply(1:n_row,FUN=function(j){rep(vx[j],time=n_col)})))  
# 
# 
# 
# sine_structure_fun <-function(x,T,phase,a,b,use_cos=FALSE){
#   
#   #Create sine for a one dimensional series
#   #Note that sine function uses radian unit.
#   #a=amplitude
#   #b=mean or amplitude 0 of the series
#   #T= stands for period definition
#   #phase=phase angle (in radian!!)
#   #cos: use cosine instead of sine if TRUE
#   
#   if(use_cos==FALSE){
#     y <- a*sin((x*pi/T)+ phase) + b
#   }else{
#     y <- a*cos((x*pi/T)+ phase) + b
#   }
#   return(y)
# }
