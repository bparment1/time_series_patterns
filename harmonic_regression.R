############################## Harmonic regression #################### 
##
## Research Support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 04/22/2019
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

crop_data_processing_functions <- "harmonic_regression_functions_04222019.R"
source(file.path(script_path,crop_data_processing_functions))

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
out_suffix <-"example_ts_04222019" #output suffix for the files and ouptut folder #param 12
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


harmonic_results <- harmonic_regression(y,n,harmonic_val=NULL,mod_obj=F,figure=F)


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

mod <- harmonic_results2$l_harmonic_obj[[1]]$mod

summary(mod)
plot(y)
lines(mod$fitted.values)

barplot(Im(fft(y)))
plot(Re(fft(y)))

#### generate function
#### split in annual windows

x <- y_all

split_sequence <- function(x,n,overlap=1){
  if(overlap==0){
    n_splits <- floor(length(x)/n)
    n_modified <- n - overlap
    intervals_val <- seq(1,to=length(x),by=n_modified)
    length(intervals_val)
    
  }
}
################################### End of script #######################################

