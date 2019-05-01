############################## Harmonic regression #################### 
##
## Research Support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 05/01/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis Managing Hurricanes
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

###Loading R library and packages                                                      
#library(gstat) #spatial interpolation and kriging methods
library(sp) # spatial/geographic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
library(readxl) #functionalities to read in excel type data
library(sf) # spatial ojbects simple feature model implementation OGC
#library(gstat)
#library(spacetime)

###### Functions used in this script and sourced from other files

#https://stats.stackexchange.com/questions/60500/how-to-find-a-good-fit-for-semi-sinusoidal-model-in-r


### This needs to be modified.
#SSTlm2 <- lm(Degrees ~ sin(2*pi*ToY)+cos(2*pi*ToY)
#             +sin(4*pi*ToY)+cos(4*pi*ToY),data=SST)
#summary(SSTlm2)

fit_harmonic <- function(p,n,y,mod_obj=F,figure=F){
  
  t <- 1:n
  omega_val=2*pi*p/n #may be more than 1
  
  cos_val <- lapply(omega_val,function(omega){cos(omega*t)})
  sin_val <- lapply(omega_val,function(omega){sin(omega*t)})
  
  cos_df <- as.data.frame(do.call(cbind,cos_val))
  names(cos_df) <- paste("cos",p,sep="")
  sin_df <- as.data.frame(do.call(cbind,sin_val))
  names(sin_df) <- paste("sin",p,sep="")
  
  in_df <- data.frame(y_var = y)
  in_df <- cbind(in_df,cos_df,sin_df)
  #View(in_df)
  #cos_val =cos(omega*t)
  #sin_val =sin(omega*t)
  
  #omega_val = lapply(p,function(p){2*pi*p/n})
  
  #plot(cos_val)
  #plot(sin_val)
  
  #y_var <- "invasion.status"
  #in_dir[[y_var]] <- as.factor(data[[y_var]]) #this is needed for randomForest to get a classification
  
  explanatory_variables <- names(in_df)[-1] #drop the first column
  
  right_side_formula <- paste(explanatory_variables,collapse = " + ")
  model_formula_str <- paste0("y_var"," ~ ",right_side_formula)
  
  #in_df <- data.frame(y=y,cos_val=cos_val,sin_val=sin_val)
  mod <- lm(model_formula_str ,data=in_df)
  summary(mod)
  #mod2 <- lm(model_formula_str ,data=in_df,na.action = na.exclude)
  #mod$coefficients
  #mod2$coefficients
  #mod2$model
  
  ### Extract coefficients for each harmonic
 
  p_val<-2
  #p
  #debug(extract_harmonic_coef)
  test_df<- extract_harmonic_coef(p_val,n,mod)
    
  harmonic_df <- lapply(p,FUN=extract_harmonic_coef,n=n,mod=mod)
  
  ### Figure
  if(figure==TRUE){
    y_range <- range(mod$fitted.values,y,na.rm = T)
    plot(mod$fitted.values,ylim=y_range)
    points(y,col="blue",pch="+")
  }
  ###
  harmonic_obj <- harmonic_df
  
  if(mod_obj==T){
    harmonic_obj <- list(harmonic_df=harmonic_df,mod=mod)
  }else{
    harmonic_obj <- list(harmonic_df=harmonic_df)
  }
  ### 
  return(harmonic_obj)
}

harmonic_regression<- function(y,n,harmonic_val=NULL,mod_obj=F,figure=F){
  ##
  # if mod_obj is True then return the model object 
  #
  ## Default to two first harmonic:
  if(is.null(harmonic_val)){
    p <- 1:2
  }else{
    p <- 1:harmonic_val
  }
  
  #harmonic_val = 1 # pr from 1 to n/2
  
  #n<-24
  
  #l_harmonic_obj <- lapply(p,
  #       FUN=fit_harmonic,
  #       n=n,
  #       y=y,
  #       mod_obj=mod_obj,
  #       figure=figure)
  
  #debug(fit_harmonic)
  l_harmonic_obj <- fit_harmonic(p,n,y,mod_obj=F,figure=F)
  
  l_df <- lapply(p,function(i){l_harmonic_obj$harmonic_df[[i]]})
  harmonic_df <- do.call(rbind,l_df)
  rownames(harmonic_df) <- NULL
  
  #View(harmonic_df)
  harmonic_results_obj  <- list(harmonic_df,l_harmonic_obj)
  names(harmonic_results_obj) <- c("harmonic_df","l_harmonic_obj")
  
  return(harmonic_results_obj)
}

split_sequence <- function(x,n,overlap=0){
  if(overlap==0){
    n_splits <- floor(length(x)/n)
    #n_modified <- n - overlap
    intervals_val <- seq(1,to=length(x),by=n)
    intervals_val <- c(intervals_val,length(x))
    n_splits
    list_intervals <- lapply(2:length(intervals_val),function(i){data.frame(start=intervals_val[[i-1]],end=intervals_val[[i]])})
    intervals_df <- do.call(rbind,list_intervals)
    intervals_df
    #lapply(2:length(intervals_val),function(i){intervals_val[[i]]-overlap})
    #length(intervals_val)
    
  }
  ##implement the other option later
  
  ## now split:
  test <- lapply(1:nrow(intervals_df),function(i){x[intervals_df[i,]$start:intervals_df[i,]$end-1]})
  
  return(test)
}


extract_harmonic_coef <- function(p_val,n,mod){
  
  summary(mod)
  coef_df <- (as.data.frame(t(mod$coefficients)))
  
  cos_term <- paste("cos",p_val,sep="")
  sin_term <- paste("sin",p_val,sep="")
  
  a <- coef_df[[sin_term]] #sine term
  b <- coef_df[[cos_term]] #cosine term
  A0 <- coef_df[["(Intercept)"]] #mean
  p_significance <- coef_df[["(Intercept)"]] #mean
  
  A = sqrt(a^2 + b^2)
  phase = atan(-b/a)
  ## Add p values later?
  #n <- nrow(mod$model)
  omega_val=2*pi*p_val/n #may be more than 1
  
  #class((summary(mod))$coefficients)
  results_df <- (as.data.frame(summary(mod)$coefficients))
  results_df <- results_df[rownames(results_df)%in%c(sin_term,cos_term,"(Intercept)"),]
  
  pr <- results_df$`Pr(>|t|)`
  pr_a <- pr[3] 
  pr_b <- pr[2]
  pr_A0 <- pr[1]
  
  harmonic_df <- data.frame(A0=A0,A=A,a=a,b=b,
                            pr_A0=pr_A0,pr_a=pr_a,pr_b=pr_b,
                            phase=phase,harmonic=p_val,omega=omega_val)
  
  return(harmonic_df)
}

################################### End of script #######################################

