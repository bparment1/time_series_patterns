############################## Harmonic regression #################### 
##
## Research Support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 04/22 /2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis Managin Hurricanes
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

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
  
  cos_df <- do.call(cbind,cos_val)
  sin_df <- do.call(cbind,sin_val)
  
  #cos_val =cos(omega*t)
  #sin_val =sin(omega*t)
  
  #omega_val = lapply(p,function(p){2*pi*p/n})
  
  #plot(cos_val)
  #plot(sin_val)
  
  y_var <- "invasion.status"
  in_dir[[y_var]] <- as.factor(data[[y_var]]) #this is needed for randomForest to get a classification
  
  explanatory_variables <- names(data)[-1] #drop the first column
  
  right_side_formula <- paste(explanatory_variables,collapse = " + ")
  model_formula_str <- paste0(y_var," ~ ",right_side_formula)
  
  in_df <- data.frame(y=y,cos_val=cos_val,sin_val=sin_val)
  mod <- lm(y~ sin_val+ cos_val ,data=in_df)
  summary(mod)
  a <- mod$coefficients[2] #sine term
  b <- mod$coefficients[3] #cosine term
  A0 <- mod$coefficients[3] #mean
  
  A = sqrt(a^2 + b^2)
  phase = atan(-b/a)
  
  ## Add p values later?
  
  harmonic_df <- data.frame(A0=A0,A=A,a=a,b=b,phase=phase,harmonic=p,omega=omega)
  
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
  
  l_df <- lapply(l_harmonic_obj,function(x){x$harmonic_df})
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
################################### End of script #######################################

