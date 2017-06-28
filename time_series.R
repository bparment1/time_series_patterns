############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/15/2017
## DATE MODIFIED: 06/28/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: fixing bug in theil Sen function
##

###################################################
#

###### Library used

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate) 

###### Functions used in this script

functions_time_series_analyses_script <- "time_series_functions_06282017b.R" #PARAM 1
functions_processing_data_script <- "processing_data_google_search_time_series_functions_06162017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2
source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data"
#ARGS 2
#infile_name <- ""
infile_name <- "vert_sp_gst_original_06122017.csv"
#ARGS 3
start_date <- "2004-01-01"
#ARGS 4
end_date <- NULL
#ARGS 5
scaling_factor <- 1000 
#ARGS 6
out_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/outputs" #parent directory where the new output directory is located
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <- "60k_time_series_analyses_06282017"
#ARGS_9
n_col_start_date <- 4
#ARGS 10
num_cores <- 8

################# START SCRIPT ###############################


### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

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


############### PART 1: Imported and time series transformation #####
## Step 1: read in the data and generate time stamps

data_ts_filename <- import_data_ts(infile_name = infile_name,
                                   in_dir = in_dir,
                                   scaling = scaling,
                                   n_col_start_date=4,
                                   start_date = start_date,
                                   end_date=NULL,
                                   out_dir = out_dir,
                                   out_suffix = out_suffix)

df <- read.table(data_ts_filename,sep=",",fill=T,header=T,check.names = F)

#dates_val <- seq(as.Date(date_range[1]), as.Date(date_range[2]), by="month")
#test<- df_dat_animals * scaling_factor

dim(df)
#View(df_dat_animals)
range_dates <- names(df)[n_col_start_date:ncol(df)]
range_dates_str <- as.character(range_dates)

range_dates<- ymd(range_dates_str)

#class(range_dates)
#start_date
#names_col <- c("g_id","sci_name","country",range_dates_str)
#df <- df[-1,-n_col] #remove the first row with incomplete header and last column with incomplete data
#names(df) <- names_col
#df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])

df_subset <- df

df_ts <- t(df[n_col_start_date:ncol(df_subset)])
df_ts <- as.data.frame(df_ts)
df_ts <- zoo(df_ts,range_dates)

names_countries <- as.character(df_subset$country)
names_species <- as.character(df_subset$sci_name)
names_species <- sub(" ","_",names_species)

names_col <- paste(names_countries,names_species,sep="_")
names(df_ts) <- names_col
View(df_ts)

########
## Example of windowing by dates

df_w_ts_ref <- window(df_ts,start=as.Date("2016-05-01"),end= as.Date("2017-05-01"))

#df_subset <- mclapply(df_ts[1:16],
#                      thresold_val=0.00004,
#              FUN=above_threshold,
#              mc.preschedule=FALSE,
#              mc.ores=num_cores)

index_selected <- lapply(df_w_ts_ref,
                      FUN=above_threshold,
                      threshold_val=0.04)
columns_selected <- unlist(index_selected)
df_ts_subset <- df_ts[,columns_selected]
dim(df_ts_subset)

#undebug(plot_ts)

test_plot <- plot_ts(df_subset,in_dir=".",scaling=1,
                     n_col_start_date=4,start_date="2004-01-01",
                     end_date=NULL,selected_countries="USA",
                     selected_species="Aegithina tiphia",save_fig = TRUE,out_dir=".", 
                     out_suffix=out_suffix)

############## PART 2: THEIL SEN AND TREND DETECTION ############

### Theil Sen slope slope calculation

#make this a function later: example with the first country
#undebug(calculate_theil_sen_time_series)
mod_mblm_test <- calculate_theil_sen_time_series(i=2,
                                                 data_df=df_ts_subset,
                                                 out_dir=out_dir,
                                                 out_suffix=out_suffix)
mod_mblm_test
#debug(calculate_theil_sen_time_series)  
list_mod_mblm_test <- lapply(1:ncol(df_ts_subset), # input parameter i as a list
                             FUN=calculate_theil_sen_time_series,
                             data_df=df_ts_subset,
                             out_dir=out_dir,
                             out_suffix=out_suffix)


list_mod_mblm_test <- lapply(1:10, # input parameter i as a list
                             FUN=calculate_theil_sen_time_series,
                             data_df=df_ts_subset,
                             out_dir=out_dir,
                             out_suffix=out_suffix)


list_mod_mblm_test <- mclapply(1:16, # input parameter i as a list
                             FUN=calculate_theil_sen_time_series,
                             data_df=df_ts_subset,
                             out_dir=out_dir,
                             out_suffix=out_suffix,
                             mc.preschedule=FALSE,
                             mc.cores=num_cores)


class(list_mod_mblm_test[[1]])
names(list_mod_mblm_test[[1]])
list_df_theil_sen <- lapply(list_mod_mblm_test,FUN=function(x){x$df_theil_sen})

df_theil_sen <- do.call(rbind,list_df_theil_sen)

View(df_theil_sen)

test_df <- arrange(df_theil_sen,desc(slope)) #order data by magnitude of slope (theil sen)
View(test_df)


###############





####
#undebug(calculate_theil_sen_time_series)

test <- trend_pattern_detection(df_ts_subset,range1=NULL,range2=NULL,out_suffix="",out_dir=".")
  
out_filename <- paste0("theil_sen_trend_detection_",out_suffix,".txt")
write.table(test,out_filename,sep=",",row.names = F)
out_filename <- paste0("theil_sen_trend_detection_",out_suffix,".csv")
write.table(test,out_filename,row.names = F)

range1 <- c("2016-05-01","2016-12-01")
range2 <- c("2016-05-01","2017-04-01")
test1 <- trend_pattern_detection(df_ts_subset,range1=range1,range2=range2,out_suffix="",out_dir=".")

range1 <- c("2016-05-01","2016-12-01")
range2 <- c("2016-05-01","2017-04-01")

test2 <- trend_pattern_detection(df_ts_subset,range1=NULL,range2=NULL,out_suffix="",out_dir=".")





#########

################### End of script ################
