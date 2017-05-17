############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/15/2017
## DATE MODIFIED: 05/16/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: quick example of time series example
##

###################################################
#

###### Functions used in this script

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate)

###### Functions used in this script


#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/research-home/bparmentier/projects/animals_trade"
#ARGS 2
infile_name <- "Gekko gecko PCA.csv"
#ARGS 3
#date_range <- c("2011.01.01","2017.03.01") #
date_range <- c("2011/1/1","2017/3/1")
#ARGS 4
scaling_factor <- 1000

################# START SCRIPT ###############################

## Step 1: read in the data and generate time stamps

dates_val <- seq(as.Date(date_range[1]), as.Date(date_range[2]), by="month")
df_dat_animals <- read.table(file.path(in_dir,infile_name),sep=",",header=T)

df_dat_animals * scaling_factor

dim(df_dat_animals)
#View(df_dat_animals)

### Step 2: Subset and transpose to create
names(df_dat_animals)
dates_years <- year(dates_val)

df_dat <- df_dat_animals[,15:89]
#df_dat <- subset(df_dat_animals,) #improve later
names_col <- as.character(df_dat_animals$country)

### PCA analyses

#make a function from previous code

### Zoo analyses

df_ts <- t(df_dat) #transpose, the result is a matrix
df_ts <- as.data.frame(df_ts) #coerce in data.frame object
names(df_ts) <- names_col #in this case, country code

df_ts <- zoo(df_ts,dates_val)
plot(df_ts[,1:3])

## Example of windowing by dates
df_w_ts <- window(df_ts,start=as.Date("2016-01-01"),end= as.Date("2017-01-01"))
plot(df_w_ts[,1:3])

df_agg_ts <- aggregate(df_ts, by=as.yearqtr)
class(df_agg_ts)

plot(df_agg_ts[,1:3])


par(mfrow=c(2,1))
plot(df_ts[,1])
plot(df_agg_ts[,1],col="red")

## Moving average, smoothing and rollingmean
rollmean_ts <- rollmean(df_ts[,1], 12)
dim(rollmean_ts) #need to fix this

par(mfrow=c(1,1))

plot(rollmean_ts)
plot((df_ts[,1]))

## Fourier Analysis, Wavelet etc.

test<-(fft(df_ts[,1]))
plot(test)
View(test)
class(fft(df_ts[,1]))

plot(Re(fft(test))^2)

### Theil Sen slope slope calculation

#make this a function later: example with the first country

time_index <- 1:nrow(df_ts)
subset_name <- names(df_ts)[1]
df_mblm <- subset(df_ts,select=subset_name) #subset using relevant name

df_mblm <- as.data.frame(df_mblm ) #convert to data.frame since it was a zoo df
df_mblm$time_index <- time_index
names(df_mblm)
formula_str <- paste0(names(df_mblm)[1]," ~ ","time_index")
formula_mblm <- as.formula(formula_str)

mod_mblm<- mblm(formula_mblm,df_mblm)
#can get the conf interval or significance if wanted
plot(mod_mblm)


################### End of script ################