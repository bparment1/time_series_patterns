############### SESYNC Research Support: Animals Trade project ########## 
##
## 
## DATE CREATED: 05/23/2017
## DATE MODIFIED: 08/11/2017
## AUTHORS: Benoit Parmentier and Elizabeth Daut
## Version: 1
## PROJECT: Animals trade
## ISSUE: 
## TO DO:
##
## COMMIT: rolling alone not working to remove seasonality; 
          #removed trend_pattern_detection2 with NO rolling for current ts
##

###################################################
#

###### Library used

library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(mblm) #Theil Sen estimator
library(lubridate) 

###### Functions used in this script

##create an output directory
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


calculate_theil_sen_time_series <- function(i,data_df,save_opt=FALSE,out_dir=".",out_suffix=""){
  #This function generates Theil Sen slope estimate using mblm package.
  #
  #INPUTS
  #1) i: column index for input data.frame, used as y variable (trend variable)
  #2) data_df: input data containing time series
  #4) out_dir: output directory
  #5) out_suffix: output suffix appended to written output on drive
  #
  #OUTPUTS
  #1)
  #
  #TO DO: 
  
  ############# BEGIN FUNCTION #############
  
  #setting up the data input: data.frame with 2 columns
  time_index <- 1:nrow(data_df) #
  subset_name <- names(data_df)[i]
  df_mblm <- subset(data_df,select=subset_name) #subset using relevant name
  
  #stores dates for later processing
  if(inherits(data_df,"zoo")){
    dates_val <- date(df_mblm)
  }else{
    dates_val <- NULL
  }
  
  df_mblm <- as.data.frame(df_mblm) #convert to data.frame since it was a zoo df
  df_mblm$time_index <- time_index
  #names(df_mblm)
  
  ### automate formula input
  formula_str <- paste0(names(df_mblm)[1]," ~ ","time_index")
  formula_mblm <- as.formula(formula_str) #transform object into formula
  
  mod_mblm <- mblm(formula_mblm,df_mblm)
  
  ##### Extract information from model object mblm

  slope_theil_sen <- coef(mod_mblm)[2]
  intercept_theil_sen <- coef(mod_mblm)[1]
  ID_ts <- subset_name
  method_str <- "theil_sen"

  if(!is.null(dates_val)){
    nt <- length(dates_val)
    start_date <- dates_val[1]
    end_date <- dates_val[nt]
  }else{
    nt <- NA
    start_date <- NA
    end_date <- NA
  }

  duration_ts <- nt
  slope_sign <- sign(slope_theil_sen)
  
  df_theil_sen <- data.frame(ID_ts=ID_ts,
                             intercept=intercept_theil_sen,
                             slope=slope_theil_sen,
                             slope_sign=slope_sign,
                             method=method_str,
                             start_date= start_date,
                             end_date = end_date,
                             duration=nt)
  
  #### Prepare object to return
  obj_theil_sen <- list(mod_mblm,df_theil_sen)
  names(obj_theil_sen) <- c("mod_mblm","df_theil_sen")
  
  ##### save to disk
  
  if(save_opt==TRUE){
    obj_theil_sen_filename <- file.path(out_dir,paste("theil_sen_obj_",subset_name,"_",out_suffix,".RData",sep=""))
    save(obj_theil_sen,file= obj_theil_sen_filename)
  } 
  
  return(obj_theil_sen)
}

plot_ts <- function(df,in_dir=".",scaling=1,
                    n_col_start_date=4,start_date="2004-01-01",
                    end_date=NULL,selected_countries=NULL,
                    selected_species=NULL,save_fig=FALSE,out_dir=".", 
                    out_suffix=""){
  #general function to plot ts data
  ## add conditions to select all countries and/or species
  #e.g. selected_countries==ALL
  #e.g. selected_species==ALL
  
  if(is.null(selected_countries)){
    selected_countries <- c("USA")
  }
  
  if(is.null(selected_species)){
    selected_species <- df$sci_name[1]
  }
  
  df_subset <- subset(df,df$sci_name%in%selected_species 
                      & df$country%in%selected_countries)
  
  #### transform subset data into a time series zoo object
  
  #TO DO:
  #df_ts 
  #
  # need to transpose, keep as a df, and add column name
  #gst60k_1<- gather(gst60k_select, "date", "gst", -name, -country) #change if selecting 1 country
  names_countries <- as.character(df_subset$country)
  names_species <- as.character(df_subset$sci_name)
  names_species <- sub(" ","_",names_species)
  
  names_col <- paste(names_countries,names_species,sep="_")
  
  n_col <- ncol(df_subset)
  df_ts <- t(df_subset[,n_col_start_date:n_col]) #transpose, the result is a matrix
  df_ts <- as.data.frame(df_ts) #coerce in data.frame object
  names(df_ts) <- names_col #in this case, country code
  
  ### checkt that we are not dropping columns
  if(is.null(end_date)){
    
    n_col <- ncol(df_subset)
    nt <- n_col - n_col_start_date #we are dropping the last date because it is often incomplete
    range_dates <- seq.Date(from=as.Date(start_date),by="month",
                            length.out = nt )
    range_dates_str <- as.character(range_dates)
    end_date <- range_dates[nt]
  }
  
  #else{
  #  df_w_ts <- window(df_ts,start=start_date,end= end_date)
  #}
  
  ## EXAMPLE of windowing by dates using the zoo package
  #df_w_ts <- window(df_ts,start=start_date,end= end_date)
  
  df_ts <- zoo(df_ts,range_dates)
  
  plot(df_ts)
  
  if(save_fig==TRUE){
    
    res_pix<-480 #set as function argument...
    col_mfrow<-1
    #row_mfrow<-2
    row_mfrow<-1
    
    png_file_name<- paste("Figure_",selected_countries,"_",selected_species,
                          out_suffix,".png", sep="")
    
    png(filename=file.path(out_dir,png_file_name),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    par(mfrow=c(row_mfrow,col_mfrow))
    
    plot(df_ts) #plot now
    dev.off()
    
  }else{
    png_file_name<- paste("Figure_",selected_countries,"_",selected_species,
                          out_suffix,".png", sep="")
  }
  
  return(png_file_name)
}


#remove low gst volume species-country combinations
above_threshold <- function(x,threshold_val=0.004) {
  #x <- x[x > 0]
  median(x, na.rm=FALSE) > threshold_val #change to see various amounts of species retained;
  #var(x,na.rm=FALSE) !=0
  # max(x, na.rm=FALSE) >0.001 #remove those with 0 values
  #median(x, na.rm=FALSE) >0.0000005 #change to see various amounts of species retained; #use this if normalized
}


trend_pattern_detection <- function(df_ts,range1=NULL,range2=NULL,roll_window=NULL,align_val="right",out_suffix="",out_dir="."){
  
  library(zoo)
  
  #add documentation
  #
  if(is.null(range1)){
    range1 <- c("2011-01-01","2016-12-01")
  }
  if(is.null(range2)){
    range2 <- c("2017-01-01","2017-05-01")
  }
  
  if(!is.null(roll_window)){
    #smooth time series if not null
    df_ts <- rollmean(df_ts,k=roll_window,align=align_val)
    #df_ts_removed_roll <- rollmean(df_ts_removed,k=12,align="right")
    
    #df_ts_tmp <- rollmean(df_ts,roll_window)
    
    #note that the time series is cut by roll_window -1 
    #e.g. if roll_window is 12, 11 timesteps are lost, 5 at the beginning and 6 at the end.
    #df_ts <- rollmean(df_ts,roll_window,fill=TRUE) #HOW TO RUN FOR EACH ROW
  }else{
    roll_window <- NA
  }
  
  #df_w_ts_ref <- window(df_ts_subset,start=as.Date(range1[1]),end= as.Date(range1[2]))
  #df_w_ts_current <- window(df_ts_subset,start=as.Date(range2[1]),end= as.Date(range2[2]))
  
  df_w_ts_current <- window(df_ts,start=as.Date(range2[1]),end= as.Date(range2[2]))  # i changed to the lines below
  df_w_ts_ref <- window(df_ts,start=as.Date(range1[1]),end= as.Date(range1[2]))
  
  # df_w_ts_current <- window(df_ts_removed,start=as.Date(range2[1]),end= as.Date(range2[2]))
  # df_w_ts_ref <- window(df_ts_removed,start=as.Date(range1[1]),end= as.Date(range1[2]))
  
  
  #1:ncol(df_w_ts
  #df_ts_removed_ratio[,1] <- mean(df_ts_removed_current[,1])/mean(df_ts_removed_ref[,1])
         
  #df_ts_ratio <- df_w_ts
  
  mean_ratio <- lapply(1:ncol(df_w_ts_ref),
         FUN=function(i){  mean_ratio_val <- mean(df_w_ts_current[,i])/mean(df_w_ts_ref[,i])}
  )
  mean_ratio <- unlist(mean_ratio)
  
  #detach_package("quantmod", TRUE)
  #library(zoo)
  
  #df_w_ts_combined <- window(df_ts,start=as.Date("2016-05-01"),end= as.Date("2017-04-01"))
  
  #plot(df_w_ts_combined[,1]) 
  #names(df_w_ts_current)[1]
  
  #df_w_ts_ref[,1] + df_w_ts_current[,1] 
  #plot(df_w_ts_current[,1])
  #plot(df_w_ts_ref[,1])
  #out_suffix="time_series_analyses_05232017")
  mod_mblm_df_w_ts_ref <- lapply(1:ncol(df_w_ts_ref), # input parameter i as a list
                                 FUN=calculate_theil_sen_time_series,
                                 data_df=df_w_ts_ref,
                                 out_dir=out_dir,
                                 out_suffix=out_suffix)
  
  mod_mblm_df_w_ts_current<- lapply(1:ncol(df_w_ts_current), # input parameter i as a list
                                    FUN=calculate_theil_sen_time_series,
                                    data_df=df_w_ts_current,
                                    out_dir=out_dir,
                                    out_suffix=out_suffix)
  
  #class(list_mod_mblm_test[[1]])
  #names(list_mod_mblm_test[[1]])
  list_df_theil_sen_ref <- lapply(mod_mblm_df_w_ts_ref,FUN=function(x){x$df_theil_sen})
  list_df_theil_sen_current <- lapply(mod_mblm_df_w_ts_current,FUN=function(x){x$df_theil_sen})
  
  df_theil_sen_ref <- do.call(rbind,list_df_theil_sen_ref)
  df_theil_sen_current <- do.call(rbind,list_df_theil_sen_current)
  
  ## compare ratio of coef
  #mod_mblm_df_w_ts_current$
  
  ### Ordering will change
  #c("ID_ts","start_date","end_date","duration","method","intercept","slope","slope_sign")
  names(df_theil_sen_ref) <- c("ID_ts","intercept1","slope1","slope_sign1","method1","start_date1","end_date1","duration1")
  names(df_theil_sen_current) <- c("ID_ts","intercept2","slope2","slope_sign2","method2","start_date2","end_date2","duration2")
  
  #merge(df_theil_sen_ref,df_theil_sen_current,by=ID_ts) #not unique identifier
  #Use cbind
  #drop method column too
  df_theil_sen_combined <- cbind(df_theil_sen_ref,df_theil_sen_current[,-1])
  
  df_theil_sen_combined$ts_ratio <- df_theil_sen_combined$slope2/df_theil_sen_combined$slope1
  df_theil_sen_combined$mean_ratio <- mean_ratio
  
  df_theil_sen_combined <- arrange(df_theil_sen_combined,desc(mean_ratio))
  
  df_theil_sen_combined$roll_window <- roll_window
  
  out_filename <- paste0("df_theil_sen_combined_",out_suffix,".csv")
  write.csv(df_theil_sen_combined,out_filename,row.names = F)
  
  return(df_theil_sen_combined)  
}




cycle_pattern_detection <-function(i,df_ts_data,method_opt="decompose",freq_val=12,save_fig=F){
  #comments about the function
  #df_ts_data = df_ts_subset
  
  ##### Begin script ####
  
  ts_zoo_obs <- df_ts_data[,i]
 
  if(method_opt=="decompose"){
    
    ts_cycle <- decompose(ts(ts_zoo_obs,frequency=freq_val))
    
    #plot(ts_decomp)
    #ts_decomp$seasonal
    #max_peak <- which.max(ts_cycle$seasonal) # this suggests peaks of seasonality in ??? 
    
    #length(ts_decomp$seasonal) - length(ts_zoo_obs) #ok no data point lost
    
    #decomp
    ts_test <- as.numeric(ts_zoo_obs) - as.numeric(ts_cycle$seasonal)
    ts_test <- zoo(ts_test,date(ts_zoo_obs))
    
  }
  
  if(method_opt=="stl"){
    
   
    ts_cycle <- stl(ts(ts_zoo_obs,frequency=freq_val),s.window=freq_val)
    #plot(ts_stl)
    #str(ts_stl)
    
    # dim(ts_stl$time.series)
    # class(ts_stl$time.series)
    # plot(ts_stl$time.series[,1]) #seasonality
    # plot(ts_stl$time.series[,2]) #trend
    # length(ts_stl$time.series[,2]) #ok 161  (or 77 if using normalized data)
    
    
    #### Now remove periodic cycle:
   # ts_test <- as.numeric(ts_zoo_obs) - as.numeric(ts_cycle$time.series[,4]) # FOR COUNTRY-SPECIES WITH KNOWN SEASONALITY, 
    # #I'M GETTING AN ERROR Error in `[.default`(ts_stl$time.series, , 4) : subscript out of bounds
    # ts_test <- zoo(ts_test,date(ts_zoo_obs))
    # #plot(ts_test)
    # #lines(ts_zoo_obs,add=T,col="red")
    # 
    # 
   # ts_test <- as.numeric(ts_zoo_obs) - as.numeric(ts_cycle$seasonal)
    ts_test <- as.numeric(ts_zoo_obs) - (ts_cycle$time.series[,1]) #### is this right??
    #ts_test <- as.numeric(ts_zoo_obs) - (ts_cycle$seasonal)
    ts_test <- zoo(ts_test,date(ts_zoo_obs))
    
  }
  
  
  #
  if(save_fig==T){
    #
    #plot(ts_test)
    #lines(ts_zoo_obs,add=T,col="red")
    
  }
  
  
  
  ### 
  cycle_obj <-  list(ts_test,ts_zoo_obs,ts_cycle)
  names(cycle_obj) <- c("ts_test","ts_zoo_obs","ts_cycle")
  return(cycle_obj)
}


normalize_by_country_totals <- function(country_name,df_input,df_gst_totals,n_col_start_date){
  #This functions normalizes data using monthly totals by countries.
  
  #### BEGIN #####
  
  n_col <- ncol(df_input)
  
  n_select <- n_col_start_date-1
  
  

  df_input_subset <- subset(df_input,country_name==df_input$country)
  df_input_info <- df_input_subset[,1:n_select]
  df_input_subset <- df_input_subset[,n_col_start_date:n_col]
  
  df_input_subset <- t(df_input_subset)
  
  
  df_gst_totals_subset <- subset(df_gst_totals,select=country_name)
  monthly_averages <- as.numeric(df_gst_totals_subset[,1])
  
  df_input_norm <- df_input_subset/monthly_averages
  
  df_input_norm <- t(df_input_norm)
  rownames(df_input_norm)<- NULL
  
  df_input_norm <- cbind(df_input_info,df_input_norm)
  #test <- cbind(df_input_info,df_input_norm)
  
  return(df_input_norm)
}


################### End of script ################
