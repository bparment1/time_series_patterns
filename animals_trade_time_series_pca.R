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
library(plotrix)

###### Functions used in this script

script_path <- "/research-home/bparmentier/Data/projects/animals_trade/scripts"
function_time_series_analyses <- "animals_trade_time_series_analyses_functions_05222017d.R" #PARAM 1
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

if(is.null(n_pca)){
  npc <- ncol(df_dat)
}

##### Performing PCA
  
## By default we use the correlation matrix in T mode, this will be changed later on.
#debug(run_pca_fun)
principal_pca_obj <- run_pca_fun(A=df_dat,mode="T",rotation_opt="none",matrix_val=NULL,npc=npc,loadings=TRUE,scores_opt=TRUE)
pca_mod <- principal_pca_obj$pca_principal

#Save eigenvalues
df_eigen <- data.frame(variable=names(df_dat),eigenvalue=pca_mod$values)
df_eigen_filename <- file.path(out_dir,paste("df_eigen_",out_suffix_str,".txt",sep=""))
write.table(df_eigen,file=df_eigen_filename,sep=",")

#Save correlation matrix
out_suffix_str <- paste0("cor_",out_suffix)
pca_input_square_matrix_filename <- file.path(out_dir,paste("pca_input_square_matrix_",out_suffix_str,".txt",sep=""))
write.table(principal_pca_obj$matrix_val,file=pca_input_square_matrix_filename,sep=",")

#Save loadings matrix
pca_loadings_filename <- file.path(out_dir,paste("pca_loadings_",out_suffix,".txt",sep=""))
write.table(principal_pca_obj$loadings,file=pca_loadings_filename,sep=",")

#Save scores if produce_scores==T
if(produce_scores==T){
  #
  pca_scores_filename <- file.path(out_dir,paste("pca_scores_",out_suffix,".txt",sep=""))
  write.table(principal_pca_obj$scores,file=pca_scores_filename,sep=",")
}

############################
##### plotting of PCAs

#Add option for non time series
plot_pca_loadings_time_series(1,pca_mod=pca_mod,dates_val,out_dir=out_dir, out_suffix=out_suffix)

lapply(1:6,FUN=plot_pca_loadings_time_series,pca_mod=pca_mod,dates_val,out_dir=out_dir,out_suffix=out_suffix)

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

######## Plot in pc space

loadings_df <- as.data.frame(pca_mod$loadings[,pcs_selected])
#pca_loadings_dz <- zoo(loadings_df,dates_val) #create zoo object from data.frame and date sequence object

pcs_selected <- c(1,2)


png_filename <- file.path(out_dir,paste("Figure_pc_components_space_loadings",pcs_selected[1],pcs_selected[2],"_",out_suffix,".png",sep=""))

res_pix<-960
col_mfrow<- 1
row_mfrow<- 1
png(filename= png_filename,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
#par(mfrow=c(1,2))

plot(loadings_df[,1],loadings_df[,2],
     type="p",
     pch = 20,
     col ="blue",
     xlab=names(loadings_df)[1],
     ylab=names(loadings_df)[2],
     ylim=c(-1,1),
     xlim=c(-1,1),
     axes = FALSE,
     cex.lab = 1.2)
axis(1, at=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1
axis(2, las=1,at=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1

box()    #This draws a box...

title(paste0("Loadings for component ", names(loadings_df)[1]," and " ,names(loadings_df)[2] ))
draw.circle(0,0,c(1.0,1),nv=200)#,border="purple",
text(loadings_df[,1],loadings_df[,2],rownames(loadings_df),pos=1,cex=1)            
grid(2,2)

dev.off()





################### End of script ################