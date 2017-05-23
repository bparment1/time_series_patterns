#################################    General script  #######################################
#############################         PCA and EOF functionality         #######################################
#This script collects functions about PCA and EOF from multiple projects'
#
#
#AUTHOR: Benoit Parmentier                                                                #
#DATE CREATED: 05/22/2017 
#DATE MODIFIED: 05/23/2017
#Version: 1
#PROJECT: General script
#   
#COMMIT: initial code PCA and EOF, utility analyses
#TODO:

#################################################################################################

#################################################################################################

###Loading R library and packages                                                      
library(gtools)                          # loading some useful tools 
library(mgcv)                            # GAM package by Simon Wood
library(sp)                              # Spatial pacakge with class definition by Bivand et al.
library(spdep)                           # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                           # GDAL wrapper for R, spatial utilities
library(gstat)                           # Kriging and co-kriging by Pebesma et al.
library(fields)                          # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                          # Hijmans et al. package for raster processing
library(foreign)                         # Library for format exchange (e.g. dbf,spss,sas etc.)
library(gdata)                           # various tools with xls reading
library(xts)                             # basic package for time series analysis
library(zoo)                             # basic package for time series analysis
library(forecast)                        # package containing ARIMA procedures
library(rasterVis)                       # plotting raster
library(nnet)                            # Contains multinom and neural net functions
library(ggplot2)                         # plotting package
library(reshape2)                        # data wrangling
library(mlogit)                          # maximum liklihood estimation and multinomial model
library(parallel)                        # parralel programming and multi cores
library(plyr)                            # miscealleous tools
library(rgeos)                           # topology and vector spatial queries and operations
library(afex)                            # functions related to ANOVA
library(car)                             # regression analysis tools
library(MASS)                            # contains negative binomial model
library(multcomp)                        # contains Tukey comparison
library(lubridate)                       # parse and manipulate dates and date-time objects in a friendly manner
library(nlme)                            # mixed models and GLS
library(psych)                           # pca/eigenvector decomposition functionalities
library(GPArotation)                     # Rotations functionality for PCA
library(plotrix)                         # plot functionality and graphic drawing tools

###### Functions used in this script

run_pca_fun <- function(A,mode="T",rotation_opt="none",matrix_val=NULL,npc=1,loadings=TRUE,scores_opt=TRUE){
  #This function performs PCA with an input data or square matrix.
  #This is to be use for smal datasets.
  
  #INPUTS:
  #1) A: data matrix 
  #2) mode: T, s mode or other mode , note that this option is not in use at this moment, use matrix_val
  #3) rotation_opt: which option to use to rotate pc components, the default is none
  #4) npc: number of components produced
  #5) matrix_val: square matrix used in the computation of eigen values, vectors
  #
  #OUTPUTS
  #obj_principal: object in the form of list with the following items
  #1) pca_principal
  #2) loadings
  #3) scores
  #4) data_matrix
  
  #######
  
  
  if(is.null(matrix_val)){
    matrix_val <- cor(A)
  }
  
  ##Cross product PCA S modes
  
  pca_principal_obj <-principal(r=matrix_val, #
                                nfactors = npc, 
                                residuals = FALSE, 
                                covar=TRUE, #use covar option...
                                rotate = rotation_opt,
                                scores=scores_opt,
                                oblique.scores=TRUE,
                                method="regression")
  
  principal_loadings_df <- as.data.frame(unclass(pca_principal_obj$loadings)) # extract the matrix of ??
  #Add scores...
  if(scores_opt==T){
    principal_scores <- as.data.frame(predict(pca_principal_obj,A))
  }
  
  ###
  obj_principal <- list(pca_principal_obj,principal_loadings_df,principal_scores,matrix_val,A)
  names(obj_principal) <- c("pca_principal","loadings","scores","matrix_val","data_matrix")
  return(obj_principal)
  
}


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

plot_pca_components_space_loadings <- function(pcs_selected,pca_mod,var_labels=NULL,out_dir=".", out_suffix=""){
  #This function generates loadings plots in pc space based on a pca model object from psych package in R.
  #
  #INPUTS
  #1) pcs_selected: components defining the pc space
  #2) pca_mod: pca model object from psych R package
  #3) var_labels: names of input variables, if null it is taken from the row names
  #4) out_dir: output directory
  #5) out_suffix: output suffix
  #
  #OUTPUTS
  #1) png_filename: png files with plot of of pca loadings
  #
  #TO DO:
  
  ########## BEGIN SCRIPT  ############
  # 
  ## Select loadings
  
  loadings_df <- as.data.frame(pca_mod$loadings[,pcs_selected])
  
  if(is.null(var_labels)){
    var_labels <- rownames(loadings_df)
  }
  
  png_filename <- file.path(out_dir,paste("Figure_pc_components_space_loadings",pcs_selected[1],"_",pcs_selected[2],"_",out_suffix,".png",sep=""))
  
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
  axis(1, at=seq(-1,1,0.2),cex=1.2)
  #axis(1, at=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1
  axis(2, las=1,at=seq(-1,1,0.2),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1
  
  box()    #This draws a box...
  
  title(paste0("Loadings for component ", names(loadings_df)[1]," and " ,names(loadings_df)[2] ))
  draw.circle(0,0,c(1.0,1),nv=200)#,border="purple",
  text(loadings_df[,1],loadings_df[,2],var_labels,pos=1,cex=1)            
  grid(2,2)
  
  dev.off()
  
  return(png_filename)
}

run_pca_analysis <- function(data_df,matrix_val=NULL,npc=1,pcs_selected=NULL,time_series_loadings=F,var_labels=NULL,mode_val=T, rotation_opt="none",scores_opt=FALSE,out_dir=".",out_suffix=""){
  ##General function to run PCA/EOF analyses. The function includes options for T and S mode or a matrix
  ##given by the user.
  #INPUTS
  #
  #1) data_df
  #2) matrix_val
  #3) npc=1
  #4)pcs_selected=NULL
  #5)time_series_loadings=F
  #6)var_labels=NULL
  #7)mode_val=T
  #8)rotation_opt="none"
  #9)scores_opt=FALSE
  #10)out_dir="."
  #11)out_suffix=""
  #
  #OUTPUTS
  #
  #
  
  ########### Beging function ###########
  
  principal_pca_obj <- run_pca_fun(A=data_df,mode=mode_val,rotation_opt="none",
                                   matrix_val=NULL,npc=npc,loadings=TRUE,scores_opt=TRUE)
  pca_mod <- principal_pca_obj$pca_principal
  
  #Save eigenvalues
  
  df_eigen <- data.frame(variable=names(data_df),eigenvalue=pca_mod$values)
  df_eigen_filename <- file.path(out_dir,paste("df_eigen_",out_suffix,".txt",sep=""))
  write.table(df_eigen,file=df_eigen_filename,sep=",")
  
  #Save correlation matrix
  #reformat, two digits
  pca_input_square_matrix_filename <- file.path(out_dir,paste("pca_input_square_matrix_",out_suffix,".txt",sep=""))
  write.table(principal_pca_obj$matrix_val,file=pca_input_square_matrix_filename,sep=",")
  
  #Save loadings matrix
  pca_loadings_filename <- file.path(out_dir,paste("pca_loadings_",out_suffix,".txt",sep=""))
  write.table(principal_pca_obj$loadings,file=pca_loadings_filename,sep=",")
  
  #Save scores if users is asking for it
  if(scores_opt==T){
    #
    pca_scores_filename <- file.path(out_dir,paste("pca_scores_",out_suffix,".txt",sep=""))
    write.table(principal_pca_obj$scores,file=pca_scores_filename,sep=",")
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
  
  ######## Plot in pc space
  
  #pcs_selected <- c(1,2)
  ## can also run a sequence
  #list_pcs_selected <- list(c(1,2),c(2,3),c(3,4))
  
  lapply(pcs_selected, 
         FUN= plot_pca_components_space_loadings,
         pca_mod,
         var_labels=var_labels,
         out_dir=".", 
         out_suffix=out_suffix)
  
  ### Plot loadings as sequence: if time series
  if(time_series_loadings=F){
    pcs_selected_unique <- unique(unlist(pcs_selected))
    #plot_pca_loadings_time_series(1,pca_mod=pca_mod,dates_val,out_dir=out_dir, out_suffix=out_suffix)
    
    lapply(pcs_selected,FUN=plot_pca_loadings_time_series,pca_mod=pca_mod,dates_val,out_dir=out_dir,out_suffix=out_suffix)
    
  }
  
  ## generate object
  
  return()
}


############### END OF SCRIPT ###################

