############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/15/2017
## DATE MODIFIED: 07/21/2017
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

functions_time_series_analyses_script <- "time_series_functions_07112017c.R" #PARAM 1
functions_processing_data_script <- "processing_data_google_search_time_series_functions_07202017.R" #PARAM 1
script_path <- "/nfs/edaut-data/Time Series MARSS" #path to script #PARAM 2
#script_path <- "/nfs/bparmentier-data/Data/projects/animals_trade/scripts" #path to script #PARAM 2
source(file.path(script_path,functions_processing_data_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_time_series_analyses_script)) #source all functions used in this script 1.


#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/nfs/edaut-data/Time Series MARSS"
#ARGS 2
infile_name <- "hi_gst60k_species_concatenated2.csv" 
#infile_name <- "hi_gst60k_species_normalized_csv2.csv" 
#ARGS 3
start_date <- "2011-01-01"
#ARGS 4
end_date <- NULL
#ARGS 5
#scaling_factor <- 100000 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
scaling_factor <- 1000 
#ARGS 6
out_dir <- "/nfs/edaut-data/Time Series MARSS/outputs" #parent directory where the new output directory is placed
#out_dir <- "/nfs/bparmentier-data/Data/projects/animals_trade/outputs"
#ARGS 7
create_out_dir_param=TRUE #create a new ouput dir if TRUE
#ARGS 8
out_suffix <- "60k_time_series_analyses_07202017_norm"
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

df_original_norm <- read.table(data_ts_filename,sep=",",fill=T,header=T,check.names = F)

df_norm <- df_original_norm

#dates_val <- seq(as.Date(date_range[1]), as.Date(date_range[2]), by="month")
#test<- df_dat_animals * scaling_factor

dim(df_norm)


############## RUN WITHOUT "ALL" CATEGORY BECAUSE WANT COUNTRY-SPECIFIC RESULTS AND "ALL" IS DOMINATING##################
#df_norm <- filter(df_norm, !(country == "All")) # All have already been removed from in imported CSV     


########### NEED TO FILTER TO KEEP JUST THOSE ROWS WHERE THE LAST VALUE IS GREATER THAN THE SECOND TO LAST VALUE ##################
names(df_norm)
selected_last <- df_norm[,c("2017-05-01")] > df_norm[,c("2017-04-01")]  #can change dates to compare different ending time periods
#sum(selected_last)
df_last <- df_norm[selected_last,]
df_norm <- df_last
#write.csv(df, "df_check_last_values.csv", row.names = FALSE) YES - last month values are greater then previous month

#View(df_dat_animals)
range_dates <- names(df_norm)[n_col_start_date:ncol(df_norm)]
range_dates_str <- as.character(range_dates)

range_dates<- ymd(range_dates_str)

#class(range_dates)
#start_date
#names_col <- c("g_id","sci_name","country",range_dates_str)
#df <- df[-1,-n_col] #remove the first row with incomplete header and last column with incomplete data
#names(df) <- names_col
#df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])

df_subset <- df_norm

df_ts <- t(df_norm[n_col_start_date:ncol(df_subset)])
df_ts <- as.data.frame(df_ts)
dim(df_ts)
df_ts_test <- df_ts[,]*100000  ############### ADDED SCALING HERE ###################################################



df_ts <- zoo(df_ts_test,range_dates)

names_countries <- as.character(df_subset$country)
names_species <- as.character(df_subset$sci_name)
names_species <- sub(" ","_",names_species)

names_col <- paste(names_countries,names_species,sep="_")
names(df_ts) <- names_col
View(df_ts)

########
## really don't want to filter again for high-gst activity b/c has already been filter for that, 
# so set the threshold really low to get all sp

df_w_ts_ref <- window(df_ts,start=as.Date("2016-05-01"),end= as.Date("2017-05-01")) # filter low volumn gst rows based on last 12 months

#df_subset <- mclapply(df_ts[1:16],
#                      thresold_val=0.00004,
#              FUN=above_threshold,
#              mc.preschedule=FALSE,
#              mc.ores=num_cores)

index_selected <- lapply(df_w_ts_ref, #for the windowed period identified in the functions file
                      FUN=above_threshold,
                      threshold_val=0.02)  #### THIS IS SO HIGH BECAUSE DATA HAVE BEEN SCALED
columns_selected <- unlist(index_selected)
df_ts_subset <- df_ts[,columns_selected]
dim(df_ts_subset)  



#########################################
# need to get list of co_sp column names to use for the graphing csv and script

View(df_ts_subset)
class(df_ts_subset)
df_ts_subset4 <- as.data.frame(df_ts_subset)
class(df_ts_subset4)
write.csv(df_ts_subset4, "df_ts_subset4.csv", row.names = FALSE)

###########################################

plot(df_ts_subset$USA_Pantherophis_spiloides)



#undebug(plot_ts)

############### WHAT IS PLOTTING FUNCTION PLOTTING??? ##################################### 

test_plot <- plot_ts(df_ts_subset,in_dir=".",scaling=1,
                     n_col_start_date=4,start_date="2011-01-01",
                     end_date=NULL,selected_countries="Japan",
                     selected_species="Acanthopagrus schlegelii",save_fig = TRUE,out_dir=".",
                     out_suffix=out_suffix)

plot(df_ts[,"USA_Aegithina tiphia"])  ###################### script out of bounds error




# NEED TO TEST WITH SOME SP-CO WHERE THE LAST MONTH INCREASES DRAMATICALLY TO SEE WHAT THE SMOOTHING DOES

## Moving average, smoothing and rollingmean
# #an example of smoothing from zoo
 rollmean_ts <- rollmean(df_ts_subset[,108], 6) #can change the length #108 is seasonal
# 
# rollmean_ts <- rollmean(df_ts_subset[,"Malaysia_Balaenoptera_musculus"], 12) #HOW TO DO BY SP_CO
# rollmean_ts <- rollmean(df_ts_subset[,], 12) #HOW TO RUN FOR EACH ROW
# rollmean_ts <- as.data.frame(rollmean_ts)
# 
# dim(rollmean_ts) #need to fix this  ############################################### 
# 
par(mfrow=c(1,1))
plot(rollmean_ts)
plot((df_ts_subset[,108]))
# plot((df_ts_subset[,"Malaysia_Balaenoptera_musculus"]))  #HOW TO DO BY SP_CO

#df_ts_subset <- rollmean_ts  # NEED TO FIX THE ROLLMEAN_TS SO IT CAN BE SEEN... A ZOO OBJECT, BUT NOT WORKING
# dim(df_ts_subset)

#########################






############## PART 2: THEIL SEN AND TREND DETECTION ############

############ SKIP DOWN TO TEST SECTION ###################

### Theil Sen slope slope calculation

#make this a function later: example with the first country
#undebug(calculate_theil_sen_time_series)

# ###### TESTING SINGLE ROWS
# mod_mblm_test <- calculate_theil_sen_time_series(i=55,
#                                                  data_df=df_ts_subset,
#                                                  out_dir=out_dir,
#                                                  out_suffix=out_suffix)
# mod_mblm_test
# 
# #debug(calculate_theil_sen_time_series)  
# 
# ##### TESTING THE ENTIRE SUBSET- FULL TIME PERIOD
# list_mod_mblm_test <- lapply(1:ncol(df_ts_subset), # input parameter i as a list
#                              FUN=calculate_theil_sen_time_series,
#                              data_df=df_ts_subset,
#                              out_dir=out_dir,
#                              out_suffix=out_suffix)

# 
# #### TESTING A FEW ROWS OF SUBSET
# list_mod_mblm_test <- lapply(1:10, # input parameter i as a list
#                              FUN=calculate_theil_sen_time_series,
#                              data_df=df_ts_subset,
#                              out_dir=out_dir,
#                              out_suffix=out_suffix)
# 
# 
# #### NEED TO FIX THE MCL APPLY
# # list_mod_mblm_test <- mclapply(1:16, # input parameter i as a list
# #                              FUN=calculate_theil_sen_time_series,
# #                              data_df=df_ts_subset,
# #                              out_dir=out_dir,
# #                              out_suffix=out_suffix,
# #                              mc.preschedule=FALSE,
# #                              mc.cores=num_cores)
# 
# ###### PRODUCES AN OUTPUT TABLE - ORIGINAL ORDER AND DESCENDING 
# 
# class(list_mod_mblm_test[[1]])
# names(list_mod_mblm_test[[1]])
# list_df_theil_sen <- lapply(list_mod_mblm_test,FUN=function(x){x$df_theil_sen})
# 
# df_theil_sen <- do.call(rbind,list_df_theil_sen)
# 
# View(df_theil_sen)
# 
# test_df <- arrange(df_theil_sen,desc(slope)) #order data by magnitude of slope (theil sen)
# View(test_df)





################################################################################################
# ## BLOCK THIS THEIL-SEN FOR NOW AND GOTO THE DECOMP/STL SECTION
# 
# ####
# #undebug(calculate_theil_sen_time_series)
# 
# # TEST (No Rolling average with default dates (full time period))
# test <- trend_pattern_detection(df_ts_subset,range1=NULL,range2=NULL,out_suffix="",out_dir=".")
# View(test)
# out_filename <- paste0("theil_sen_trend_detection_test",out_suffix,".csv")
# write.csv(test,out_filename,row.names = F)
# 
# # TEST_1 (Rolling average with default dates: REF = "2011-01-01","2016-12-01" AND CURRENT = "2017-01-01","2017-05-01"
# test_1 <- trend_pattern_detection(df_ts_subset,range1=NULL,roll_window = 6,range2=NULL,out_suffix="",out_dir=".")
# # does not work when roll_window >7 with current time period = 5 months b/c of the way rolling shortens the ts
# View(test_1)
# out_filename <- paste0("theil_sen_trend_detection_test_1",out_suffix,".csv")
# write.csv(test_1,out_filename,row.names = F)
# 
# 
# 
# 
# # NORMALIZED TESTING #########################
# 
# # TEST2  (this is the same as 'test' above) FULL time period
# range1 <- c("2011-01-01","2016-12-01") #reference
# range2 <- c("2017-01-01","2017-05-01") #current
# test_1n <- trend_pattern_detection(df_ts_subset,range1=range1,range2=range2,out_suffix="",out_dir=".")
# View(test_1n)
# out_filename <- paste0("theil_sen_trend_detection_test_1n_",out_suffix,".csv")
# write.csv(test_1n,out_filename,row.names = F)
# 
# # TEST2b with rolling
# range1 <- c("2011-01-01","2016-11-01") #reference
# range2 <- c("2016-12-01","2017-05-01") #current
# test_2n <- trend_pattern_detection(df_ts_subset,range1=range1,roll_window = 6,range2=range2,out_suffix="",out_dir=".")
# View(test_2n)
# out_filename <- paste0("theil_sen_trend_detection_test_2n_",out_suffix,".csv")
# write.csv(test_2n,out_filename,row.names = F)
# 
# 
# 
# # TEST_n1 
# #range1 <- c("2013-12-01","2016-12-01") #reference  #MIDWAY - 1/2 OF TIME OF TIME SERIES EXCLUDING LAST 5 TEST MONTHS
# #range1 <- c("2011-01-01","2016-12-01") #reference  #FULL - time period
# range1 <- c("2015-06-01","2016-12-01") #reference  # total 2 years
# range2 <- c("2017-01-01","2017-05-01") #current
# 

############################ this plotting-option section won't work b/c the object is not created. 
# if want to be able to plot directly see Function from 7/10/17

# ### PATTERN 1 - ENTIRE TS IS ROLLED
# out_suffix_str <- paste0(c("test_3n.1"),"_",out_suffix)
# test_3n.1 <- trend_pattern_detection(df_ts_subset,
#                                     range1=range1,
#                                     range2=range2,
#                                     roll_window=6,
#                                     out_suffix=out_suffix_str,
#                                     out_dir=out_dir)
# 
# #names(test_3n.1)
# 
# df_w_ts_ref<-  test_3n.1$df_w_ts_ref
# df_w_ts_current<-  test_3n.1$df_w_ts_current
# names(test_3n.1$df_w_ts_ref)
# plot(df_w_ts_ref[,98])
# plot(df_w_ts_current[,98])
# 
# 
# plot_ts_cosp<- function(country_species_name,df_time_series){
#   
#   #country_species_name <- "USA_Ursus_americanus"
#   #country_species_name <- "Indonesia_Sus_scrofa"
#   #plot(subset(df_w_ts_ref,select=country_species_name))
#   #plot(subset(df_ts_subset,select=country_species_name))
#   plot(subset(df_time_series,select=country_species_name))
#   #plot(subset(df_ts_,select=country_species_name))
#   
# }
# 
# 
# plot_ts_cosp("Malaysia_Balaenoptera_musculus",df_ts_subset)
# plot_ts_cosp("Malaysia_Balaenoptera_musculus",df_w_ts_ref)
# plot_ts_cosp("Malaysia_Balaenoptera_musculus",df_w_ts_current)
# 
# # plot(test_3n.1$df_w_ts_ref)
# # dim(df_w_ts_ref)
# 
# #View(test_3n.2) # HOW CAN I VIEW?
# out_filename <- paste0("theil_sen_trend_detection_test_3n.1",out_suffix,".csv")
# write.csv(test_3n.1,out_filename,row.names = F)







#START AGAIN FROM HERE TO RUN WITH STL TO COMPARE TO DECOMP



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


# FOR GRAPHS
df_ts_removed_stl <- as.data.frame(df_ts_removed)
class(df_ts_removed_stl)
write.csv(df_ts_removed_stl, "df_ts_removed_stl.csv", row.names = FALSE)

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


range1 <- c("2015-04-01","2016-11-01") #reference  # total 2 years; 6month current
range2 <- c("2016-12-01","2017-05-01") #current
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

# TEST (No Rolling average with default dates (full time period))
test_stl <- trend_pattern_detection_stl(df_ts_removed,range1=NULL,range2=NULL,out_suffix="",out_dir=".")
View(test_stl)
out_filename <- paste0("theil_sen_trend_detection_test_stl_full",out_suffix,".csv")
write.csv(test_stl,out_filename,row.names = F)


####################################################

# PLOT SPECIFIC COUNTRY-SPECIES COMBOS

plot(df_ts_subset$USA_Agelasticus_thilius) #normalized-original
lines(df_ts_removed$USA_Agelasticus_thilius,col="red") #de-seasoned

####################################################




range1 <- c("2015-04-01","2016-11-01") #reference  # total 2 years; 6month current
range2 <- c("2016-12-01","2017-05-01") #current
# TEST (No Rolling average)
test_stl_a <- trend_pattern_detection_stl(df_ts_removed,range1=range1,range2=range2,out_suffix="",out_dir=".")
View(test_stl_a)
out_filename <- paste0("theil_sen_trend_detection_test_stl_a_",out_suffix,".csv")
write.csv(test_stl_a,out_filename,row.names = F)


# TEST_1 (Rolling average with default dates: REF = "2011-01-01","2016-12-01" AND CURRENT = "2017-01-01","2017-05-01"
test_stl.2 <- trend_pattern_detection_stl(df_ts_removed,range1=NULL,roll_window = 4,range2=NULL,out_suffix="",out_dir=".")
# does not work when roll_window >7 with current time period = 5 months b/c of the way rolling shortens the ts
View(test_stl.2)
out_filename <- paste0("theil_sen_trend_detection_test_stl2_",out_suffix,".csv")
write.csv(test_stl.2,out_filename,row.names = F)


#range1 <- c("2013-12-01","2016-12-01") #reference  #MIDWAY - 1/2 OF TIME OF TIME SERIES EXCLUDING LAST 5 TEST MONTHS
range1 <- c("2011-01-01","2016-07-01") #reference  #FULL - time period
# range1 <- c("2015-06-01","2016-12-01") #reference  # total 2 years; 5month current
# range2 <- c("2017-01-01","2017-05-01") #current
#range1 <- c("2015-04-01","2016-09-01") #reference  # total 2 years; 6month current
#range2 <- c("2016-12-01","2017-05-01") #current
range2 <- c("2016-08-01","2017-05-01") #current

# TEST_2  
test_stl.8r <- trend_pattern_detection_stl(df_ts_removed,range1=range1,roll_window =5,range2=range2,out_suffix="",out_dir=".")
View(test_stl.8r)
out_filename <- paste0("theil_sen_trend_detection_test_stl.3r_",out_suffix,".csv")
write.csv(test_stl.3r,out_filename,row.names = F)


# PLOT SPECIFIC COUNTRY-SPECIES COMBOS

plot(df_ts_subset$USA_Bunolagus_monticularis) #normalized-original
lines(df_ts_removed$USA_Bunolagus_monticularis,col="red") #de-seasoned



################### End of script ################
