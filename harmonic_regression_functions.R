############################## Harmonic regression #################### 
##
## Research Support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 04/17/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis Managin Hurricanes
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

###### Functions used in this script and sourced from other files

fit_harmonic <- function(p,n,y,mod_obj=F,figure=F){
  t <- 1:n
  omega=2*pi*p/n
  cos_val =cos(omega*t)
  sin_val =sin(omega*t)
  
  #plot(cos_val)
  #plot(sin_val)
  
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
  }
  #harmonic_val = 1 # pr from 1 to n/2
  
  #n<-24
  
  #omega_val = lapply(p,function(p){2*pi*p/n})
  
  l_harmonic_obj <- lapply(p,
         FUN=fit_harmonic,
         n=n,
         y=y,
         mod_obj=mod_obj,
         figure=figure)
  
  l_df <- lapply(l_harmonic_obj,function(x){x$harmonic_df})
  harmonic_df <- do.call(rbind,l_df)
  rownames(harmonic_df) <- NULL
  
  #View(harmonic_df)
  harmonic_results_obj  <- list(harmonic_df,l_harmonic_obj)
  names(harmonic_results_obj) <- c("harmonic_df","l_harmonic_obj")
  
  return(harmonic_results_obj)
}

################################### End of script #######################################

