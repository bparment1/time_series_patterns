############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/23/2017
## DATE MODIFIED: 06/26/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: quick example of time series example ##SHOULD CHANGE THIS
##

###################################################
#

###### Library used

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(graphics)

#####  Functions used in this script

functions_time_series_analyses_script <- "time_series_functions_05252017.R"  #MAY NEED TO CHANGE THE DATE IF I UPDATE THE FUNCTIONS SCRIPT
script_path <- "/nfs/edaut-data/Time Series MARSS"
source(file.path(script_path,functions_time_series_analyses_script)) #to direct it to my functions file


##########  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/nfs/edaut-data/Time Series MARSS"
#ARGS 2
infile_name <- "vert_sp_gst_original_06122017.csv" 
#ARGS 3
start_date <- "2004-01-01"
#ARGS 4
end_date <- NULL
#ARGS 5
#scaling_factor <- 100000 #MODIFY THE SCALING FACTOR - FOR NORMALIZED DATA SHOULD BE 10,000 AT LEAST
scaling_factor <- 1000 
#ARGS 6
out_dir <- "/nfs/edaut-data/Time Series MARSS/outputs" #parent directory where the new output directory is placed
#ARGS 7
create_out_dir_param=TRUE #create a new output dir if TRUE
#ARGS 8
out_suffix <- "60k_time_series_analyses_06262017"  # CHANGE THIS FOR EACH OUTPUT??
#ARGS_9
n_col_start_date <- 4


################# START SCRIPT #######################

### PART I: READ AND PREPARE DATA FOR REGRESSIONS #######
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

################# NOTES #####################

#normalize the data for those countries where i have total gst (or for all of them depending on what Michel sends me)
#filter to remove species-country combinations that have median values (excluding zero values)
      # <0.0000005 (once normalized) 
#filter full dataset for just those species i want (once i have the data from google)

###########################################

## Step 1: read in the data and generate time stamps

data_ts_filename <- import_data_ts(infile_name = infile_name,
                                   in_dir = in_dir,
                                   scaling = scaling,
                                   n_col_start_date=4,
                                   start_date = start_date,
                                   end_date=NULL,
                                   out_dir = out_dir,
                                   out_suffix = out_suffix)

df <- read.table(file.path(in_dir,infile_name),sep=",",fill=T,quote="",header=F)

#remove first row that contains names
#rename columns using:
#n_col_start_date <- 3 #may need to be changed
n_col <- ncol(df)
nt <- n_col - n_col_start_date #drop the last month because it is often incomplete

range_dates <- seq.Date(from=as.Date(start_date),by="month",length.out = nt )
range_dates_str <- as.character(range_dates)
#class(range_dates)
#start_date
names_col <- c("g_id","sci_name","country",range_dates_str)
df <- df[-1,-n_col] #remove the first row with incomplete header and last column with incomplete data
names(df) <- names_col
df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])


#remove low gst volume species-country combinations
above_threshold <- function(x) {
  #x <- x[x > 0]
  median(x, na.rm=FALSE) >0.004 #change to see various amounts of species retained;
  # max(x, na.rm=FALSE) >0.001 #remove those with 0 values
  #median(x, na.rm=FALSE) >0.0000005 #change to see various amounts of species retained; #use this if normalized
}

#select high gst volume rows for last 12 months
rows_to_accept <- apply(df[151:164],1, above_threshold)
hi_gst_species <- df[rows_to_accept, ]
# names(hi_gst_species)
df_ts <- hi_gst_species[, 4:164]#*scaling_factor ### how to multilply these columns and keep the others
#### WHY IS THE SCALING FACTOR ABOVE NOT WORKING?



# need to transpose, keep as a df, and add column name
df_ts <- gather(hi_gst_species, "date", "gst") #use this if selecting all countries
#df_ts <- gather(gst60k_c, "date", "gst", -name, -country) #use this if selecting 1 country

#make real date for R code
df_ts$date <- as.Date(df_ts$date, format = "X%m.%d.%Y")





#df_dat <- subset(df_dat_animals,) #improve later ############################
names_country <- as.character(df_ts$country) #adds country names back in after transposing
names_species <- as.character(df_ts$sci_name) #need to add underscore btwn genus_species
names_species <- sub(" ","_",names_species)

names_col <- paste(names_country,names_species,sep="_")



##############  Threshold MEDIAN value #####################

# Run function
filter_median_threshold <- median_gst_above_threshold 







###### Zoo analyses ###########

df_ts <- t(df_ts) #transpose, the result is a matrix
df_ts <- as.data.frame(df_ts) #coerce in data.frame object
names(df_ts) <- names_col #in this case, country code # when transposed, lost the names y need to add back

df_ts <- zoo(df_ts,dates_val) #making a zoo object - can select a window of time; used with ts
plot(df_ts[,1:3])  





#################


## EXAMPLE of windowing by dates using the zoo package
df_w_ts <- window(df_ts,start=as.Date("2016-01-01"),end= as.Date("2017-01-01"))
plot(df_w_ts[,1:3]) #this is plotting a zoo object (check class - should be zoo) which is different from regular plotting

df_agg_ts <- aggregate(df_ts, by=as.yearqtr) # depending on the input - in this case from zoo
#the as.yearqtr is a zoo function; the aggregate function will aggregate quarterly; could also aggregate annually
#df_agg_ts <- aggregate(df_ts, by=year) aggregates yearly
class(df_agg_ts)
dim(df_agg_ts)


#plot(df_agg_ts[,1:3]) #the default is to plot separately
plot(df_agg_ts[,1:3], plot.type="single") #this plots all 3 on the same graph #???? WOULD NEED A LEGEND Y COLORS IF DID IT THIS WAY.
  #BETTER TO DO REAL PLOTS WITH GGPLOT??




### to compare monthly and quarterly
# should be able to plot both on the same plot  #########???? THESE PLOTS ARE OF THE FIRST COLUMN, RIGHT?  NEED LEGEND
par(mfrow=c(2,1)) #plots 2 rows y 1 column
plot(df_ts[,1])
plot(df_agg_ts[,1],col="red")


#################### DON'T UNDERSTAND THIS ###############

## Moving average, smoothing and rollingmean
#an example of smoothing from zoo
rollmean_ts <- rollmean(df_ts[,1], 12) #can change the length
    # can also do the rollmedian, rollmax
    #could look at the difference from the max-mean to find peaks
    #the weight is just one
    #if don't overlap the windows
    #can also "align" to R or L; use the center one for smoothing; and the running mean for prediction
    #could have very small windows and not overlap
    #can do rollapply with a function to identify the peaks
dim(rollmean_ts) #need to fix this

par(mfrow=c(1,1))
plot(rollmean_ts)
plot((df_ts[,1]))




############ needs more work

## Fourier Analysis, Wavelet etc.
# an easy way to find seasonal variation (cycles of any length)

# test<-(fft(df_ts[,1]))
# plot(test)
# View(test)
# class(fft(df_ts[,1]))
# 
# plot(Re(fft(test))^2)





############### Theil Sen slope slope calculation ######################

#to test the fucntion without it being in the function
# mod_mblm_test <- calculate_theil_sen_time_series(i=1,
#                                                  data_df=df_ts,
#                                                  out_dir=out_dir,
#                                                  out_suffix="time_series_analyses_05232017")


#testing as a loop using lapply - a loop across a list
#the "i" from above is now the first imput = the 1:ncol
list_mod_mblm_test <- lapply(1:ncol(df_ts), #this generates a list
                             FUN=calculate_theil_sen_time_series,
                             data_df=df_ts,
                             out_dir=out_dir,
                             out_suffix=out_suffix)

mod_mblm <- list_mod_mblm_test[[3]]
mod_mblm
coef(mod_mblm)

############ NOT SURE WHAT THIS DOES ###############

list_theil_sen_obj <- list.files(path=out_dir,pattern="*.RData")
test_obj <- load_obj(list_theil_sen_obj[[1]])
test_obj

##### TESTING 

#in the mean time can use this rolling (raw) mean 

rollmean_ts <- rollmean(df_ts[,1:10], 12) 

mod_mblm_rollmean <- calculate_theil_sen_time_series(i=1,
                                                        data_df=df_w_ts_ref,
                                                        out_dir=out_dir,
                                                        out_suffix="time_series_analyses_05232017")
#this coefficient (0.004) is different than the non-smoothed the theilsen which was 0.0008
coef(mod_mblm_rollmean)

####################################

# Set up windowing to run theil-sen slope estimator to detect the 3 gst patterns
# w=window; ref=reference period (versus current)
#REALLY ANNOYING TO HAVE TO THINK ABOUT THE MONTHS TO CALCULATE... HOW TO AUTOMATE??


#1. short-term acute increase
  #compare the last 4 months (current) to the previous 8 months (ref)
df_w_ts_ref <- window(df_ts,start=as.Date("2016-05-01"),end=as.Date("2016-12-01")) #previous 8 months
df_w_ts_current <- window(df_ts,start=as.Date("2017-01-01"),end=as.Date("2017-04-01")) #last 4 months
df_w_ts_combined <- window(df_ts,start=as.Date("2016-05-01"),end=as.Date("2017-04-01")) #entire test window

  #compare the last 4 months (current) to the previous 12 months (ref)
df_w_ts_ref <- window(df_ts,start=as.Date("2016-01-01"),end=as.Date("2016-12-01")) #previous 12 months
df_w_ts_current <- window(df_ts,start=as.Date("2017-01-01"),end=as.Date("2017-04-01")) #last 4 months
df_w_ts_combined <- window(df_ts,start=as.Date("2016-01-01"),end=as.Date("2017-04-01")) #entire test window





#2. long-term steady increase
  #compare the last 4 months (current) to the previous 24 months (ref)
df_w_ts_ref <- window(df_ts,start=as.Date("2015-01-01"),end=as.Date("2016-12-01"))
df_w_ts_current <- window(df_ts,start=as.Date("2017-01-01"),end=as.Date("2017-04-01"))
df_w_ts_combined <- window(df_ts,start=as.Date("2015-01-01"),end=as.Date("2017-04-01"))

  #compare the last 4 months (current) to the previous 48 months (ref)
df_w_ts_ref <- window(df_ts,start=as.Date("2014-01-01"),end=as.Date("2016-12-01"))
df_w_ts_current <- window(df_ts,start=as.Date("2017-01-01"),end=as.Date("2017-04-01"))
df_w_ts_combined <- window(df_ts,start=as.Date("2014-01-01"),end=as.Date("2017-04-01"))

  #compare the last 4 months (current) to the previous 72 months (ref)
df_w_ts_ref <- window(df_ts,start=as.Date("2011-01-01"),end=as.Date("2016-12-01"))
df_w_ts_current <- window(df_ts,start=as.Date("2017-01-01"),end=as.Date("2017-04-01"))
df_w_ts_combined <- window(df_ts,start=as.Date("2011-01-01"),end=as.Date("2017-04-01")) #entire ts



# running column by column (must select the column)
names(df_w_ts_current)[94]
plot(df_w_ts_ref[,1])
plot(df_w_ts_current[,1])
plot(df_w_ts_combined[,1])

mod_mblm_df_w_ts_ref <- calculate_theil_sen_time_series(i=94,
                                                 data_df=df_w_ts_ref,
                                                 out_dir=out_dir,
                                                 out_suffix="time_series_analyses_05232017")


mod_mblm_df_w_ts_current <- calculate_theil_sen_time_series(i=94,
                                                        data_df=df_w_ts_current,
                                                        out_dir=out_dir,
                                                        out_suffix="time_series_analyses_05232017")

mod_mblm_df_w_ts_combined <- calculate_theil_sen_time_series(i=94,
                                                            data_df=df_w_ts_combined,
                                                            out_dir=out_dir,
                                                            out_suffix="time_series_analyses_05232017")

coef(mod_mblm_df_w_ts_ref)
coef(mod_mblm_df_w_ts_current)
coef(mod_mblm_df_w_ts_combined)
coef(mod_mblm_df_w_ts_ref)[2] #this is the slope time index (the second output component)
theil_sen_ratio <- (coef(mod_mblm_df_w_ts_current)[2]/coef(mod_mblm_df_w_ts_ref)[2])*100  ###  SHOULD THE BE REVERSED???
theil_sen_ratio

##### need to convert the zoo object to a time-series object - do it automatically

plot(df_ts[,1])
nt <- nrow(df_ts) #to keep track of the number of rows
as.numeric((df_ts[,1]),start=df_ts)




##### then decompose the ts oject first, then run the theil-sen on the smoothed data

decompose(df_ts[,1]) - #this expects a ts object (instead of a zoo object)
df_ts[1,1]
0.0008147/df_ts[1,1]*100 #we have the rate of increase (the theilsen slope)
(0.365+(0.0008147)*12)-0.365 #the first value in the series (and then the other stuff) = the value in 12 months


plot(df_w_ts[,1:3]) #this is plotting a zoo object (check class - should be zoo) which is different from regular plotting

df_agg_ts <- aggregate(df_ts, by=as.yearqtr) # depending on the input - in this case from zoo
#the as.yearqtr is a zoo function; the aggregate function will aggregate quarterly; could also aggregate annually
#df_agg_ts <- aggregate(df_ts, by=year) aggregates yearly
class(df_agg_ts)
dim(df_agg_ts)


#plot(df_agg_ts[,1:3]) #the default is to plot separately
plot(df_agg_ts[,1:3], plot.type="single") #this plots all 3 on the same graph #???? WOULD NEED A LEGEND Y COLORS IF DID IT THIS WAY.
#BETTER TO DO REAL PLOTS WITH GGPLOT??


### to compare monthly and quarterly
# should be able to plot both on the same plot  #########???? THESE PLOTS ARE OF THE FIRST COLUMN, RIGHT?  NEED LEGEND
par(mfrow=c(2,1)) #plots 2 rows y 1 column
plot(df_ts[,1])
plot(df_agg_ts[,1],col="red")






#####
#1. Mann-Kendall test with 4 to 6 month window periods to test for significant differences
#among the trends - set up the 3 categores - increasing, decreasing, neutral
# 2.then Thiel-sen


#### HOW TO SET UP SO ALERTS FOR THE POSTIIVE SPECIES-COUNTRY COMBINATIONS?


# to have the output sent to a separate fil - need to create in the functions file






############ Lowess Curve Spike Detection ###############

#Part 1: Prepare data

df_lc <- df_dat

#make real date for R code
#fb1$date <- as.Date(fb1$date, format = "X%m.%d.%Y")

#select first row to test
df_lc <- df_lc[1,2:ncol(df_lc)]

# Run function
mod_lowess_test <- spike_detection_test(i=1,
                                        data_df=df_lc,
                                        out_dir=out_dir,
                                        out_suffix="time_series_analyses_05232017")

# Run function as a loop with lapply
# list_mod_lowess_test <- lapply(1:ncol(df_lc), #this generates a list
#                              FUN=spike_detection_test,
#                              data_df=df_lc,
#                              out_dir=out_dir,
#                              out_suffix=out_suffix)


##### MAKES NEW OUTPUT FOLDER FOR EACH RUN? How to get output in useful format with positive sp-countries listed??? ##############

################### End of script ################