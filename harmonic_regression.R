############################## Harmonic regression #################### 
##
## Research Support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 05/01/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis Managing Hurricanes
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

#Benoit setup
script_path <- "/nfs/bparmentier-data/Data/projects/managing_hurricanes/scripts"

crop_data_processing_functions <- "harmonic_regression_functions_05012019d.R"
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
y <- y_all[1:23]
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
                                         mod_obj=T,figure=F)
harmonic_results2 <- harmonic_regression(y,n,
                                        harmonic_val=NULL,
                                        mod_obj=T,
                                        figure=F)

mod <- harmonic_results2$l_harmonic_obj[[1]]$mod

summary(mod)
plot(y)
lines(mod$fitted.values)

barplot(Im(fft(y)))
plot(Re(fft(y)))

#### generate function
#### split in annual windows

x <- y_all

debug(split_sequence)
split_obj <- split_sequence(y_all,n=23)
split_obj$list_y[[9]]
length(split_obj$list_y[[9]])

harmonic_results3 <- harmonic_regression(test[[1]],n=23,
                                         harmonic_val=NULL,
                                         mod_obj=T,figure=F)

mod <- harmonic_results2$l_harmonic_obj[[1]]$mod

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

#https://matinbrandt.wordpress.com/2013/11/15/pixel-wise-time-series-trend-anaylsis-with-ndvi-gimms-and-r/
  
### need to check the computation of Amplitudes!!!!

harmonic_reg_f1 <- function(y){
  harmonic_results <-harmonic_regression(y,n=23,
                        harmonic_val=NULL,
                        mod_obj=T,figure=F)
  df_in <- subset(harmonic_results$harmonic_df,harmonic==1)
  A <- df_in$A
  
  return(A)
}

harmonic_reg_f1 <- function(y,n=24,harmonic=1){
  harmonic_results <-harmonic_regression(y,n=n,
                                         harmonic_val=NULL,
                                         mod_obj=T,figure=F)
  df_in <- subset(harmonic_results$harmonic_df,harmonic==1)
  A <- df_in$A
  
  return(A)
}

#Does work with calc
#Note that A has two outputs right now and it creates multiple outputs in raster
harmonic_reg_raster <- function(y,n=24,harmonic=1){
  harmonic_results <-harmonic_regression(y,n=n,
                                         harmonic_val=NULL,
                                         mod_obj=T,figure=F)
  df_in <- subset(harmonic_results$harmonic_df,harmonic==harmonic)
  A <- df_in$A
  
  return(A)
}

#fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}

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
  r_out <- try(calc(subset(r,end_val:start_val), 
               fun=function(y){harmonic_reg_raster(y,n=n_val,harmonic=harmonic_val)}))
  
  l_r_A1[[i]] <- r_out
  #plot(p)
  #rm(p)
}

r_A1 <- stack(l_r_A1)
plot(r_A1)
l_r_A1
plot(r_A1,y=1:2)

################################### End of script #######################################

