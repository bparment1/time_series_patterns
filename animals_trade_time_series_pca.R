############### SESYNC Research Support: Spatial demography ########## 
## Performing PCA on animals trade data at SESYNC.
## 
## DATE CREATED: 06/15/2017
## DATE MODIFIED: 08/16/2017
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
#script_path <- "C:/Users/edaut/Documents/gst_ts" #path to script #PARAM 2
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2
source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.

function_pca_eof <- "pca_eof_functions_08022017.R" #PARAM 1
source(file.path(script_path,function_pca_eof)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
#in_dir <- "C:/Users/edaut/Documents/gst_ts"
in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data"

#ARGS 2
infile_name <- "vert_sp_gst_original_08162017.csv"
infile_name_gst_totals <- "total_monthly_gst_averages.csv"

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
out_suffix <-"vert_sp_pca_08022017" #output suffix for the files and ouptut folder #param 12

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

infile_name <- "vert_sp_gst_original_08162017.csv"
#infile_name <- file.path(in_dir,infile_name)

data_ts_filename <- import_data_ts(infile_name = infile_name,
                                   in_dir = in_dir,
                                   scaling = scaling,
                                   n_col_start_date=4,
                                   start_date = start_date,
                                   end_date=NULL,
                                   out_dir = out_dir,
                                   out_suffix = out_suffix)

df_original <- read.table(data_ts_filename,sep=",",fill=T,header=T,check.names = F)

#Set up dates as column names  ###############  
#IS THIS NEEDED b/c ALREADY IN PROCESSING??? ###############

range_dates <- names(df_original)[n_col_start_date:ncol(df_original)]
range_dates_str <- as.character(range_dates)
range_dates<- ymd(range_dates_str)

#Transform and scale data
df_subset <- df_original
df_ts <- t(df_original[n_col_start_date:ncol(df_subset)])
df_ts <- as.data.frame(df_ts)
#dim(df_ts)
df_ts_scaled <- df_ts[,]*scaling_factor  

#create ts zoo object
df_ts <- zoo(df_ts_scaled,range_dates)
#class(df_ts)
#combine country-species as column names
names_countries <- as.character(df_subset$country)
names_species <- as.character(df_subset$sci_name)
names_species <- sub(" ","_",names_species)
names_col <- paste(names_countries,names_species,sep="_")
names(df_ts) <- names_col
#View(df_ts)


########### REMOVE SPECIES-COUNTRY COMBINATIONS THAT ARE OF NO INTEREST
# HOW TO DO THIS NOW THAT IN COLUMNS ????????????????
# remove_sp <- readLines("remove_sp_co_60K.txt") # remove species in countries that have no interest
# df_ts <- filter(df_ts, !(sci_name %in% remove_sp_co_60K))

#Remove country-species combos with low gst volumes 
# really don't want to filter again for high-gst activity b/c has already been filter for that, 
# so set the threshold really low to get all sp


##### PART 2: Run different variants of PCA analysis

### 1) Run PCA on genes ncis with original data

## By default we use the correlation matrix in T mode,
#this will be changed later on.

data_df <- as.data.frame(t(df_ts))

selected_var1 <- names(data_df)[1:163]
out_suffix_str <- paste0("all_species_countires_analyses1_cor_",out_suffix)
npc <- 10
var_labels <- selected_var1

pcs_selected <- list(c(1,2),c(2,3),c(3,4),c(4,5))

#debug(run_pca_analysis)  
pca_obj1 <- run_pca_analysis(data_df,
                             matrix_val=NULL,
                             npc=npc,
                             pcs_selected=pcs_selected,
                             time_series_loadings=T,
                             var_labels=var_labels,
                             mode_val=T,
                             rotation_opt="none",
                             scores_opt=TRUE,
                             out_dir=".",
                             out_suffix=out_suffix_str )


scores_df <- pca_obj1$principal_pca_obj$scores

loadings_df <- pca_obj1$principal_pca_obj$loadings

dim(loadings_df)

#create ts zoo object
loadings_df_ts <- zoo(loadings_df,range_dates)
plot(loadings_df_ts[,])

dim(scores_df)
names(scores_df)

country_species <- rownames(data_df)
scores_df$country_species <- country_species
rownames(scores_df) <- country_species

scores_ordered_df <- arrange(scores_df,desc(PC1))
scores_ordered_df[1:20,c(11,1)]

### Select highest scores for PC1
selected_country_names <- scores_ordered_df[1:10,11] #[1] "MEX_Ailuropoda_melanoleuca"

rownames(df_ts)

df_ts[,c("MEX_Ailuropoda_melanoleuca")]
df_ts_selected <- df_ts[,c(selected_country_names)]
plot(df_ts_selected)

# ### 2) Run PCA on animals trade with screeing
# 
# 
# ## By default we use the correlation matrix in T mode,
# #this will be changed later on.
# 
# #Scale values to ease computation of eigenvalues
# selected_var2 <- names(genes_ncis_df)[1:1006]
# out_suffix_str <- paste0("ncis_analyses2_scaled_cor_",out_suffix)
# npc <- 10
# var_labels <- selected_var2
# 
# pcs_selected <- list(c(1,2),c(2,3),c(3,4),c(4,5))
# 
# #debug(run_pca_analysis)  
# pca_obj2 <- run_pca_analysis(genes_ncis_rescaled_df,
#                              matrix_val=NULL,
#                              npc=npc,
#                              pcs_selected=pcs_selected,
#                              time_series_loadings=F,
#                              var_labels=var_labels,
#                              mode_val=T,
#                              rotation_opt="none",
#                              scores_opt=FALSE,
#                              out_dir=".",
#                              out_suffix=out_suffix_str )

####
### 3) Run PCA on genes ncis square matrix of data

# ## By default we use the correlation matrix in T mode,
# #this will be changed later on.
# 
# #Scale values to ease computation of eigenvalues
# #genes_ncis_rescaled_df <- genes_ncis_df[,2:1007]*1000
# 
# selected_var3 <- names(genes_ncis_rescaled_df)[1:1006]
# out_suffix_str <- paste0("ncis_analyses3_scaled_square_matrix_",out_suffix)
# npc <- 10
# var_labels <- selected_var3
# 
# pcs_selected <- list(c(1,2),c(2,3),c(3,4),c(4,5))
# 
# #debug(run_pca_analysis)  
# matrix_square <- t(as.matrix(genes_ncis_rescaled_df)) %*% as.matrix(genes_ncis_rescaled_df)
#   
# pca_obj3 <- run_pca_analysis(genes_ncis_rescaled_df,
#                              matrix_val=matrix_square,
#                              npc=npc,
#                              pcs_selected=pcs_selected,
#                              time_series_loadings=F,
#                              var_labels=var_labels,
#                              mode_val=T,
#                              rotation_opt="none",
#                              scores_opt=FALSE,
#                              out_dir=".",
#                              out_suffix=out_suffix_str )
# 
# 
# ## screening should be done with input of "last_month_no"
# df_w_ts_screening <- window(df_ts,start=as.Date("2016-05-01"),end= as.Date("2017-05-01")) # HOW TO DO THIS BY LAST 12 MONTHS AUTOMATICALLY
# 
# index_selected <- lapply(df_w_ts_screening, 
#                          FUN=above_threshold,
#                          threshold_val=0.02)  #### THIS IS SO HIGH BECAUSE DATA HAVE BEEN SCALED
# columns_selected <- unlist(index_selected)
# df_ts_subset <- df_ts[,columns_selected]
# dim(df_ts_subset)  


############################## END OF SCRIPT #############################################

