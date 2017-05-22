############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/22/2017
## DATE MODIFIED: 05/22/2017
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

###### Functions used in this script

script_path <- "/research-home/bparmentier/Data/projects/animals_trade/scripts"
function_time_series_analyses <- "animals_trade_time_series_analyses_functions_05222017.R" #PARAM 1
source(file.path(script_path,function_time_series_analyses ))

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
out_suffix <-"animals_trade_time_series_05222017" #output suffix for the files and ouptu folder #PARAM 8
#ARGS 8
create_out_dir_param <- TRUE 
#ARGS 9
n_pca <- 6 #number of pca to produce 
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

if(is.null(n_pca)){
  n_pca <- ncol(df_dat)
}

##### Performing PCA

## This is T-mode!!! should have option for S-mode
pca_mod <-principal(df_dat,nfactors=n_factors,rotate="none",covar = FALSE)

#Save correlation matrix

#Save loadings matrix

#Save eigenvalues

#Save scores if produce_scores==T


##### plotting of PCAs

plot_pca_loadings_time_series(1,pca_mod=pca_mod,dates_val,out_dir=out_dir, out_suffix=out_suffix)

lapply(1:10,FUN=plot_pca_loadings_time_series,pca_mod=pca_mod,dates_val,out_dir=out_dir,out_suffix=out_suffix)

plot_pca_loadings_time_series <- function(i,pca_mod,dates_val,out_dir=".", out_suffix=""){
  #This function generates loadings plots based on a pca model object from psych package in R.
  #
  #INPUTS
  #1) i: component i e.g. 1
  #2) pca_mod: pca model object from psych R package
  #3) dates_val: list of dates matching the temporal sequence
  #4) out_dir: output directory
  #5) out_suffix: output suffix
  #
  #OUTPUTS
  #1) png_filename: png files with plot of of pca loadings
  #
  #TO DO: add option for mutiple plots
  
  ########## BEGIN SCRIPT  ############
  # 
  ## Select loadings
  
  loadings_df <- as.data.frame(pca_mod$loadings[,i])
  pca_loadings_dz <- zoo(loadings_df,dates_val) #create zoo object from data.frame and date sequence object
  
  
  ##Figure components
  
  png_filename <- file.path(out_dir,paste("Figure_pc_components_",i,"_",out_suffix,".png",sep=""))
  
  res_pix<-960
  col_mfrow<- 1
  row_mfrow<- 1
  png(filename= png_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  #par(mfrow=c(1,2))
  
  plot(pca_loadings_dz,
       type="b",
       xlab="time steps",
       ylab="PC loadings",
       ylim=c(-1,1),
       col="blue")
  
  title(paste0("Loadings for component ", i))
  
  names_vals <- paste0("pc",i)
  
  legend("topright",legend=names_vals,
         pt.cex=0.8,cex=1.1,col=c("black"),
         lty=c(1,1), # set legend symbol as lines
         pch=1, #add circle symbol to line
         lwd=c(1,1),bty="n")
  
  dev.off()
  
  return(png_filename)
  
}


## Add scree plot of eigen values
png_filename <- file.path(out_dir,paste("Figure_scree_plots_eigen_values_components_",out_suffix,".png",sep=""))

res_pix<-960
col_mfrow<- 1
row_mfrow<- 1
png(filename= png_filename,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(pca_mod$values,
     main="Scree plot: Eigen values",
     type="b",
     ylab="Eigenvalue",
     xlab="components")

dev.off()

## Add scree plot of variance from components
png_filename <- file.path(out_dir,paste("Figure_scree_plots_percent_variance_components_",out_suffix,".png",sep=""))

res_pix<-960
col_mfrow<- 1
row_mfrow<- 1
png(filename= png_filename,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(pca_mod$values/sum(pca_mod$values)*100,
     main="Scree plot: Variance from component in %",
     type="b",
     ylab="% variance",
     xlab="components")

dev.off()


################### End of script ################