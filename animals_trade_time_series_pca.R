############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/22/2017
## DATE MODIFIED: 05/24/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 2
## PROJECT: Animals trade
## Issue: #21977
## link: https://base.sesync.org/issues/21977
## 
## TO DO:
##
## COMMIT: animals trade, issue #21977, time series processing pca and more
##

###################################################
#

### Loading R library and packages  

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate) # parse and manipulate dates and date-time objects in a friendly manner
library(nlme) #mixed models and GLS
library(forecast) #moving average, ARIMA and other tools
library(psych) #pca/eigenvector decomposition functionalities
library(GPArotation) #rotation functions for PCA and EFA
library(plotrix) #additional plotting and drawing options 

###### Functions used in this script

script_path <- "/research-home/bparmentier/Data/projects/animals_trade/scripts"
function_time_series_analyses <- "animals_trade_time_series_analyses_functions_05232017.R" #PARAM 1
functions_pca_script <- "pca_eof_functions_05242017.R" #PARAM 1
source(file.path(script_path,function_time_series_analyses ))
source(file.path(script_path,functions_pca_script)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/research-home/bparmentier/Data/projects/animals_trade"
#ARGS 2
out_dir <- "/research-home/bparmentier/Data/projects/animals_trade/outputs"
#ARGS 3
infile_name <- "Gekko gecko PCA.csv" #input dataset 
#ARGS 4
#date_range <- c("2011.01.01","2017.03.01") #
date_range <- c("2011/1/1","2017/3/1")
#ARGS 5
scaling_factor <- 1000
#ARGS 6
ref_poly_shp_fname <- ""  #country shapefile
#ARGS 7
out_suffix <-"animals_trade_time_series_05242017" #output suffix for the files and ouptu folder #PARAM 8
#ARGS 8
create_out_dir_param <- TRUE 
#ARGS 9
npc <- 6 #number of pca to produce 
#Args 10
produce_scores <- TRUE

################# START SCRIPT ###############################


### PART I: READ AND PREPARE DATA FOR ANALYSES #######


## First create an output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

## Step 1: read in the data and generate time stamps

dates_val <- seq(as.Date(date_range[1]), as.Date(date_range[2]), by="month")
df_dat_animals <- read.table(file.path(in_dir,infile_name),sep=",",header=T)

df_dat_animals <-  df_dat_animals * scaling_factor

dim(df_dat_animals)
#View(df_dat_animals)

### Step 2: Subset and transpose to create
names(df_dat_animals)
dates_years <- year(dates_val)

df_dat <- df_dat_animals[,15:89]
#df_dat <- subset(df_dat_animals,) #improve later
names_col <- as.character(df_dat_animals$country)

###########################
### PART II: PCA analyses #######

## Use revised functions that combines everything:

selected_var1 <- names(df_dat)
npc <- length(selected_var1)

data_df_subset <- subset(df_dat,select=selected_var1)

names(data_df_subset)

## By default we use the correlation matrix in T mode, this will be changed later on.
#debug(run_pca_fun)
out_suffix_str <- paste0("pca_analyses1_cor_",out_suffix)
npc <- length(selected_var1)
var_labels <- selected_var1
pcs_selected <- list(c(1,2),c(2,3),c(3,4),c(4,5))
#undebug(run_pca_analysis)  
run_pca_analysis(data_df=data_df_subset,
                 matrix_val=NULL,
                 npc=npc,
                 pcs_selected=pcs_selected,
                 time_series_loadings=F,
                 save_opt=FALSE,
                 var_labels=var_labels,
                 mode_val=T, 
                 rotation_opt="none",
                 scores_opt=produce_scores,
                 out_dir=out_dir,
                 out_suffix=out_suffix_str)

#functions_pca_script <- "pca_eof_functions_05232017c.R" #PARAM 1
#source(file.path(script_path,functions_pca_script)) #source all functions used in this script 1.

################### End of script ################