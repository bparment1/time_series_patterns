############### SESYNC Research Support: Animals Trade project ########## 
## General purpose function to analyze time series data.
## Includes:
## 1) PCA analyses (both time series and other)
## 2) Theil Sen slope estimator
## 3) Mann Kenall trend test 
##
## DATE CREATED: 05/22/2017
## DATE MODIFIED: 05/23/2017
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
library(GPArotation) #Rotations functionality for PCA

###### Functions used in this script

create_dir_fun <- function(out_dir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}


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

plot_pca_copmonents_space_loadings <- function(pcs_selected,pca_mod,var_labels=NULL,out_dir=".", out_suffix=""){
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
  axis(2, las=1,at=seq(-1,1,0.2),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1
  
  box()    #This draws a box...
  
  title(paste0("Loadings for component ", names(loadings_df)[1]," and " ,names(loadings_df)[2] ))
  draw.circle(0,0,c(1.0,1),nv=200)#,border="purple",
  text(loadings_df[,1],loadings_df[,2],var_labels,pos=1,cex=1)            
  grid(2,2)
  
  dev.off()
  
  return(png_filename)
}



################### End of script ################