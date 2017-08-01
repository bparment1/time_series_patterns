############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/15/2017
## DATE MODIFIED: 07/17/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: using normalized data with common species removed
##

###################################################
#

###### Library used

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate)


###### Functions used in this script

functions_time_series_analyses_script <- "time_series_functions_07272017_rolling_right.R" #PARAM 1
functions_processing_data_script <- "processing_data_google_search_time_series_functions_07202017.R" #PARAM 1
#script_path <- "C:/Users/edaut/Documents/gst_ts" #path to script #PARAM 2
script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2
source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.


#####  Parameters and argument set up ###########

#ARGS 1
#in_dir <- "C:/Users/edaut/Documents/gst_ts"
in_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/data"

#ARGS 2
#infile_name <- "hi_gst60k_species_concatenated2.csv" #hand normalized
infile_name <- "vert_sp_gst_original_06122017.csv"
infile_name_gst_totals <- "total_monthly_gst_averages.csv"

#ARGS 3
start_date <- "2004-01-01"
#ARGS 4
end_date <- NULL
#ARGS 5
scaling_factor <- 100000 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
#scaling_factor <- 1000 
#ARGS 6
#out_dir <- NULL #"C:/Users/edaut/Documents/gst_ts/outputs" #parent directory where the new output directory is placed
out_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/outputs"
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <- "60k_ts_analyses_08012017c"
#ARGS_9
n_col_start_date <- 4
#ARGS 10
num_cores <- 8

################# START SCRIPT ###############################


### PART I. READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}

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

df_original <- read.table(data_ts_filename,sep=",",fill=T,header=T,check.names = F)

df_gst_totals <- read.table(file.path(in_dir,infile_name_gst_totals),
                            sep=",",header=T)

#View(df_gst_totals)
#names(df_original)
#View(df_original)
##### Normalization 

selected_country_names_totals <- names(df_gst_totals)
df_original_subset <- subset(df_original,selected_country_names_totals %in% df_original$country)
df_original_subset <- subset(df_original_subset, df_original_subset$country != c("All"))

#View(df_original_subset)
rownames(df_original_subset) <- NULL


test<- normalize_by_country_totals(df_input=df_original_subset,df_gst_totals,country_name="USA",n_col_start_date)
  
normalize_by_country_totals <- function(df_input,df_gst_totals,country_name,n_col_start_date){
  #This functions normalizes data using monthly totals by countries.
  
  #### BEGIN #####
  
  n_col <- ncol(df_input)
  
  n_select <- n_col_start_date-1

  
  df_input_subset <- df_input[,n_col_start_date:n_col]
  df_input_subset <- subset(df_input_subset,country_name==df_input$country)
  df_input_info <- df_input_subset[,1:n_select]
  
  df_input_subset <- t(df_input_subset)
  

  df_gst_totals_subset <- subset(df_gst_totals,select=country_name)
  monthly_averages <- as.numeric(df_gst_totals_subset[,1])
  
  df_input_norm <- df_input_subset/monthly_averages
  
  df_input_norm <- t(df_input_norm)
  rownames(df_input_norm)<- NULL
  
  df_import_norm <- cbind(df_input_info,df_input_norm)

  return(df_input_norm)
}

#### 

df_norm <- df_original_norm
#dim(df_norm)

#remove "ALL" country category (when using full dataset)
#df_norm <- filter(df_norm, !(country == "All"))  


#Select for increasing trends in last month
names(df_norm) #confirm column names
selected_last <- df_norm[,c("2017-05-01")] > df_norm[,c("2017-04-01")]  
#sum(selected_last) #how many maintained in df
df_last <- df_norm[selected_last,]
df_norm <- df_last
#dim(df_norm)


#Set up dates as column names  ###############  
#IS THIS NEEDED b/c ALREADY IN PROCESSING??? ###############

range_dates <- names(df_norm)[n_col_start_date:ncol(df_norm)]
range_dates_str <- as.character(range_dates)
range_dates<- ymd(range_dates_str)

#Transform and scale data
df_subset <- df_norm
df_ts <- t(df_norm[n_col_start_date:ncol(df_subset)])
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

##### SHOULDN'T USE df_w_ts_ref BECAUSE IT'S USED BELOW FOR SOMETHING ELSE????

## screening should be done with input of "last_month_no"
df_w_ts_screening <- window(df_ts,start=as.Date("2016-05-01"),end= as.Date("2017-05-01")) # HOW TO DO THIS BY LAST 12 MONTHS AUTOMATICALLY
index_selected <- lapply(df_w_ts_screening, 
                      FUN=above_threshold,
                      threshold_val=0.02)  #### THIS IS SO HIGH BECAUSE DATA HAVE BEEN SCALED
columns_selected <- unlist(index_selected)
df_ts_subset <- df_ts[,columns_selected]
dim(df_ts_subset)  




############## PART 2: THEIL SEN AND TREND DETECTION ############



################################################################################################
# ## BLOCK THIS THEIL-SEN FOR NOW AND GOTO THE DECOMP/STL SECTION


################### Dealing with Seasonality ##########

# NEED TO CONSIDER ADDITIVE OR MULTIPLICATIVE;  HOW TO DETERMINE FREQUENCY AUTOMATICALLY (AUTOMATED)?

### let's analyze a specific time series
#names(df_ts)
i <- 1 

#undebug(cycle_pattern_detection)
cycle_obj<- cycle_pattern_detection(i,df_ts_data=df_ts_subset,method_opt="stl",freq_val=12,save_fig=F)
#cycle_obj<- cycle_pattern_detection(i,df_ts_data=df_ts_subset,method_opt="decompose",freq_val=12,save_fig=F)

n_col <- ncol(df_ts_subset)
cycle_obj <- lapply(i:n_col,
       FUN=cycle_pattern_detection,
       df_ts_data=df_ts_subset,
       method_opt="stl",
       #method_opt="decompose",
       freq_val=12,
       save_fig=F)


cycle_obj[[i]]$ts_zoo_obs #original data
plot(cycle_obj[[4]]$ts_zoo_obs,col="red") #original data
lines(cycle_obj[[4]]$ts_test) #de-seasonalized data
cycle_obj[[i]]$ts_test #checking

list_df_ts_removed<- lapply(1:length(cycle_obj),
                            function(i,x){as.numeric(x[[i]]$ts_test)},
                            x=cycle_obj)


mat_removed <- do.call(cbind,list_df_ts_removed)  
df_removed <- as.data.frame(mat_removed)
df_ts_removed <- zoo(df_removed,date(df_ts_subset))
names(df_ts_removed) <- names(df_ts_subset)
  
View(df_ts_removed)
dim(df_ts_removed)
class(df_ts_removed)

#####################################
##### Add code to compute de ratio: 07/24, modified 08/01
df_ts_removed[,1]
mean(df_ts_removed[,1])

range1 <- c("2011-01-01","2016-11-01") #reference  # total 2 years; 6month current
range2 <- c("2016-12-01","2017-05-01") #current
df_ts_removed_ref <- window(df_ts_removed,start=range1[1],end=range1[2])
df_ts_removed_current <- window(df_ts_removed,start=range2[1],end=range2[2])

plot(df_ts_removed[,1])

df_ts_removed_ratio <- df_ts_removed

df_ts_removed_ratio[,1] <- mean(df_ts_removed_current[,1])/mean(df_ts_removed_ref[,1])
class(df_ts_removed_current)
class(df_ts_removed_ref)

### New updated function 07/24/2017
#undebug(trend_pattern_detection)
test_ratio <- trend_pattern_detection(df_ts_removed,
                                      range1=NULL,
                                      range2=NULL,
                                      roll_window=12,
                                      align_val="right",
                                      out_suffix="",
                                      out_dir=".")


# FOR GRAPHS
# df_ts_removed_stl <- as.data.frame(df_ts_removed)
# class(df_ts_removed_stl)
# write.csv(df_ts_removed_stl, "df_ts_removed_stl.csv", row.names = FALSE)

# df_ts_removed_decom <- as.data.frame(df_ts_removed)
# class(df_ts_removed_decom)
# write.csv(df_ts_removed_decom, "df_ts_removed_decom.csv", row.names = FALSE)

# # i wanted to use this decomposed data in theilsen, so i converted it to a dataframe (many-stepped process below) 
# 
# #convert zoo object to dataframe
# write.zoo(
#   df_ts_removed,
#   index.name = "Date",
#   file = "tmp.txt",
#   sep = ",",
#   col.names = TRUE
# )
# df_ts_removed.df <- read.csv("tmp.txt", sep = ',')
# 
# #write.csv(df_ts_removed.df, "df_ts_removed_df_stl.csv", row.names = FALSE)
# write.csv(df_ts_removed.df, "df_ts_removed_df_decom.csv", row.names = FALSE)
# #df_ts_removed_df<-read.csv("df_ts_removed_df4_stl_a.csv", header=TRUE, stringsAsFactors = FALSE) 
# df_ts_removed_df<-read.csv("df_ts_removed_df4_decom_a.csv", header=TRUE, stringsAsFactors = FALSE) 
# #import as a dataframe; had to change formatting - paste co-sp column names (copied from df_ts_subset4) and removed date column
# # must change working directory to have it read from the main folder
# #for the graphs code
# # transform; add dates, make country column, clean up sp names (remove co and _ )
# #must move copy of df_ts_removed_df from the output folder to the main folder so it can be found
# 
# 
# 
# 
# # add dates back in 
# #range_dates <- names(df_ts_removed_df)[n_col_start_date:ncol(df_ts_removed_df)]
# range_dates <- names(df)[n_col_start_date:ncol(df)]
# range_dates_str <- as.character(range_dates)
# range_dates<- ymd(range_dates_str)
# 
# # df_ts_decom <- zoo(df_ts_removed_df,range_dates)
# # 
# # dim(df_ts_decom)
# # View(df_ts_decom)
# # class(df_ts_decom)
# 
# df_ts_stl <- zoo(df_ts_removed_df,range_dates)
# 
# dim(df_ts_stl)
# View(df_ts_stl)
# class(df_ts_stl)


# TRY THE SAME PROCESS WITH STL







######### try Theil-sen on decomposed data ##############

# TEST (No Rolling average with default dates (full time period))
test_decom <- trend_pattern_detection_decom(df_ts_removed,range1=NULL,range2=NULL,out_suffix="",out_dir=".")
View(test_decom)
out_filename <- paste0("theil_sen_trend_detection_test_decom_full",out_suffix,".csv")
write.csv(test_decom,out_filename,row.names = F)


range1 <- c("2015-04-01","2016-10-01") #reference  # total 2 years; 6month current
range2 <- c("2016-11-01","2017-05-01") #current
# TEST (No Rolling average)
test_decom_a <- trend_pattern_detection_decom(df_ts_decom,range1=range1,range2=range2,out_suffix="",out_dir=".")
View(test_decom_a)
out_filename <- paste0("theil_sen_trend_detection_test_decom_a_",out_suffix,".csv")
write.csv(test_decom_a,out_filename,row.names = F)


# TEST_1 (Rolling average with default dates: REF = "2011-01-01","2016-12-01" AND CURRENT = "2017-01-01","2017-05-01"
test_decom2 <- trend_pattern_detection_decom(df_ts_decom,range1=NULL,roll_window = 5,range2=NULL,out_suffix="",out_dir=".")
# does not work when roll_window >7 with current time period = 5 months b/c of the way rolling shortens the ts
View(test_decom2)
out_filename <- paste0("theil_sen_trend_detection_test_decom2_",out_suffix,".csv")
write.csv(test_decom2,out_filename,row.names = F)


#range1 <- c("2013-12-01","2016-12-01") #reference  #MIDWAY - 1/2 OF TIME OF TIME SERIES EXCLUDING LAST 5 TEST MONTHS
#range1 <- c("2011-01-01","2016-12-01") #reference  #FULL - time period
range1 <- c("2015-06-01","2016-12-01") #reference  # total 2 years

range2 <- c("2017-01-01","2017-05-01") #current


# MAKE CURRENT TIME PERIOD 6 MONTHS
#range1 <- c("2013-12-01","2016-12-01") #reference  #MIDWAY - 1/2 OF TIME OF TIME SERIES EXCLUDING LAST 5 TEST MONTHS
#range1 <- c("2011-01-01","2016-12-01") #reference  #FULL - time period
range1 <- c("2015-04-01","2016-11-01") #reference  # total 2 years; 6month current

range2 <- c("2016-12-01","2017-05-01") #current

# TEST_2  
test_decom4 <- trend_pattern_detection_decom(df_ts_decom,range1=range1,roll_window = 4,range2=range2,out_suffix="",out_dir=".")
View(test_decom4)
out_filename <- paste0("theil_sen_trend_detection_test_decom4_",out_suffix,".csv")
write.csv(test_decom4,out_filename,row.names = F)








######### try Theil-sen on STL data ##############

test_stl_roll_right <- trend_pattern_detection(df_ts_removed,
                                      range1=NULL,
                                      range2=NULL,
                                      roll_window=3,
                                      #align_val="right",
                                      out_suffix="",
                                      out_dir=".")

View(test_stl_roll_right)


# TEST (No Rolling average with default dates (full time period))
test_stl <- trend_pattern_detection_stl(df_ts_removed,range1=NULL,range2=NULL,out_suffix="",out_dir=".")
View(test_stl)
out_filename <- paste0("theil_sen_trend_detection_test_stl_full",out_suffix,".csv")
write.csv(test_stl,out_filename,row.names = F)


####################################################

# PLOT SPECIFIC COUNTRY-SPECIES COMBOS

plot(df_ts_subset$Indonesia_Bubo_scandiacus) #normalized-original
lines(df_ts_removed$Indonesia_Bubo_scandiacus,col="red") #de-seasoned

####################################################



range1 <- c("2011-01-01","2016-09-01") #reference  # total 2 years; 6month current
range2 <- c("2016-10-01","2017-05-01") #current

range1 <- c("2015-09-01","2016-10-01") #reference  #5months# 
range2 <- c("2016-11-01","2017-05-01") #current #5months

# TEST (No Rolling average)
test_stl.4 <- trend_pattern_detection_stl(df_ts_removed,range1=range1,range2=range2,out_suffix="",out_dir=".")
View(test_stl.4)
out_filename <- paste0("theil_sen_trend_detection_test_stl_a_",out_suffix,".csv")
write.csv(test_stl_a,out_filename,row.names = F)


test_stl_roll_right <- trend_pattern_detection(df_ts_removed,
                                               range1=range1,
                                               range2=range2,
                                               roll_window=3,
                                               #align_val="right",
                                               out_suffix="",
                                               out_dir=".")

View(test_stl_roll_right)






# TEST_1 (Rolling average with default dates: REF = "2011-01-01","2016-12-01" AND CURRENT = "2017-01-01","2017-05-01"
test_stl.4a <- trend_pattern_detection_stl(df_ts_removed,range1=NULL,roll_window = 3,range2=NULL,out_suffix="",out_dir=".")
# does not work when roll_window >7 with current time period = 5 months b/c of the way rolling shortens the ts
View(test_stl.4a)
out_filename <- paste0("theil_sen_trend_detection_test_stl2_",out_suffix,".csv")
write.csv(test_stl.2,out_filename,row.names = F)


#range1 <- c("2013-12-01","2016-12-01") #reference  #MIDWAY - 1/2 OF TIME OF TIME SERIES EXCLUDING LAST 5 TEST MONTHS
#range1 <- c("2011-01-01","2016-09-01") #reference  #FULL - time period
range1 <- c("2015-06-01","2016-12-01") #reference  # total 2 years; 5month current
# range2 <- c("2017-01-01","2017-05-01") #current
#range1 <- c("2015-04-01","2016-09-01") #reference  # total 2 years; 6month current
#range2 <- c("2016-12-01","2017-05-01") #current
range2 <- c("2016-10-01","2017-05-01") #current

range1 <- c("2015-09-01","2016-10-01") #reference  #5months# 
range2 <- c("2016-11-01","2017-05-01") #current #5months




# TEST_2  
test_stl.5r <- trend_pattern_detection_stl(df_ts_removed,range1=range1,roll_window =4,range2=range2,out_suffix="",out_dir=".")
View(test_stl.5r)
out_filename <- paste0("theil_sen_trend_detection_test_stl.3r_",out_suffix,".csv")
write.csv(test_stl.3r,out_filename,row.names = F)


# PLOT SPECIFIC COUNTRY-SPECIES COMBOS


plot(df_ts_subset$USA_Acanthonus_armatus) #normalized-original
lines(df_ts_removed$USA_Acanthonus_armatus,col="red") #de-seasoned



# # TRYING TO PLOT ROLLED DATA
# rollmean_ts <- rollmean(df_ts_removed_stl[,108], 12) #can change the length #108 is seasonal
# par(mfrow=c(1,1))
# plot(rollmean_ts)
# lines((df_ts_removed_stl[,108]))
# 



################### End of script ################
