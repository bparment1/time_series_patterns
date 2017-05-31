

in_dir <- "/nfs/edaut-data/Time Series MARSS"
infile_name <- "gst40k_original_5262017.csv"
out_dir <- "/nfs/edaut-data/Time Series MARSS/outputs"
out_suffix <- "imported_ts_05312017"

scaling <- 1000
start_date="2004-01-01"

import_data_ts <- function(infile_name,in_dir=".",scaling=1,n_col_start_date=4,start_date="2004-01-01",end_date=NULL,out_dir=".", out_suffix=""){
  #Import data sent by google and format columns for future processing and visualization.
  
  ##### BEGIN FUNCTION ####
  
  ### Add quote="" otherwise EOF warning and error in reading
  df <- read.table(file.path(in_dir,infile_name),sep=",",fill=T,quote="",header=F)
  
  #remove first row that contains names
  #rename columns using:
  #n_col_start_date <- 3 #may need to be changed
  n_col <- ncol(df)
  nt <- n_col - n_col_start_date #we are dropping the last date because it is often incomplete
  
  range_dates <- seq.Date(from=as.Date(start_date),by="month",length.out = nt )
  range_dates_str <- as.character(range_dates)
  #class(range_dates)
  #start_date
  names_col <- c("g_id","sci_name","country",range_dates_str)
  df <- df[-1,-n_col] #remove the first row with incomplete header and last column with incomplete data
  names(df) <- names_col
  df[,n_col_start_date] <- sub('"',"",df[,n_col_start_date])
  
  out_filename_tmp <- sub(extension(infile_name),"",infile_name)#
  out_filename <- file.path(out_dir,paste0(out_filename_tmp,"_",out_suffix,".csv"))
  write.table(df,filename=out_filename)
  return(out_filename)
}
