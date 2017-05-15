############### SESYNC Research Support 
##
## The relationship between PCA and FA is also shown.
##
## DATE CREATED: 05/15/2017
## DATE MODIFIED: 05/15/2017
## Version: 1
## PROJECT: 
##
## TO DO:
## COMMIT: quick example of time series example
##

###################################################
#

###### Functions used in this script

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate)
#install.packages("lubridate")

###### Functions used in this script


#####  Parameters and argument set up ###########

#ARGS 1
in_dir <- "/research-home/bparmentier/projects/animals_trade"
#ARGS 2
infile_name <- "Gekko gecko PCA.csv"

#date_range <- c("2011.01.01","2017.03.01") #
date_range <- c("2011/1/1","2017/3/1")


################# START SCRIPT ###############################

## Step 1: read in the data and generate time stamps

dates_val <- seq(as.Date(date_range[1]), as.Date(date_range[2]), by="month")
df_dat_animals <- read.table(file.path(in_dir,infile_name),sep=",",header=T)

dim(df_dat_animals)
View(df_dat_animals)

### Step 2: Subset and transpose to create
names(df_dat_animals)
dates_years <- year(dates_val)
#dates_year %in% names(df_dat_animals)
pattern_str <- as.character(dates_years)
grep(pattern, x, ignore.case = FALSE, perl = FALSE, value = FALSE,
     fixed = FALSE, useBytes = FALSE, invert = FALSE)

df_dat <- subset(df_dat_animals,)

### PCA analyses

### Zoo analyses

### Theil Sen slop

################### End of script ################