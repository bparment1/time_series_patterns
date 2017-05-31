library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(grid)

##########  Parameters and argument set up ###########
#ARGS 1
in_dir <- "/nfs/edaut-data/Time Series MARSS"
#ARGS 2
infile_name <-"gst40k_original_5262017.csv"
#ARGS 3
start_date="2004-01-01" 
#ARGS 4
scaling_factor <- 1000 
#ARGS 5
out_dir <- "/nfs/edaut-data/Time Series MARSS/outputs" #parent directory where the new output directory is placed
#ARGS 6
create_out_dir_param=TRUE #create a new output dir if TRUE
#ARGS 7
out_suffix <- "imported_ts_05312017"  # CHANGE THIS FOR EACH OUTPUT

#testing to create the function - must give the arguments without the word function y {}
#calling the function
lf_processed <- import_data_ts(infile_name=infile_name,
               in_dir=in_dir,
               scaling=1,
               n_col_start_date=4,
               start_date="2004-01-01",
               end_date=NULL,
               out_dir=out_dir,
               out_suffix=out_suffix)

infile_name <- lf_processed[1] #for testing
  

  
####### BEGIN FUNCTION #######
#for testing you want to read the data just for testing
plot_ts <- function(infile_name,in_dir=".",scaling=1,n_col_start_date=4,start_date="2004-01-01",end_date=NULL,out_dir=".", out_suffix=""){
 #general function to plot ts data

  
  #df <- read.table(infile_name,sep=",",header=T) #just for testing
  df <- read.table(file.path(in_dir,infile_name),sep=",",header=T) 
  
############ TO DO
  #copy code from the theilsen file that subsets by country, dates, species??
  #plotting code
  
  
##################
  
  
  ## EXAMPLE of windowing by dates using the zoo package
  df_w_ts <- window(df_ts,start=as.Date("2016-01-01"),end= as.Date("2017-01-01"))
  plot(df_w_ts[,1:3]) #this is plotting a zoo object (check class - should be zoo) which is different from regular plotting
  
  df_agg_ts <- aggregate(df_ts, by=as.yearqtr) # depending on the input - in this case from zoo
  #the as.yearqtr is a zoo function; the aggregate function will aggregate quarterly; could also aggregate annually
  #df_agg_ts <- aggregate(df_ts, by=year) aggregates yearly
  class(df_agg_ts)
  dim(df_agg_ts)
  
  
  #plot(df_agg_ts[,1:3]) #the default is to plot separately
  plot(df_agg_ts[,1:3], plot.type="single") #this plots all 3 on the same graph #???? WOULD NEED A LEGEND Y COLORS IF DID IT THIS WAY.
  #BETTER TO DO REAL PLOTS WITH GGPLOT??
  
  
  
  
  ### to compare monthly and quarterly
  # should be able to plot both on the same plot  #########???? THESE PLOTS ARE OF THE FIRST COLUMN, RIGHT?  NEED LEGEND
  par(mfrow=c(2,1)) #plots 2 rows y 1 column
  plot(df_ts[,1])
  plot(df_agg_ts[,1],col="red")

  ######################
  
  
  #remove first row y replace with column names
  n_col_start_date <- 4
  nt <- ncol(df) - n_col_start_date #number of columns - named columns
  names_col <- c("g_id", "sci_name","country","")
   
  
  gst<-read.table("gst_month_w_values.csv", sep=",", header=TRUE)[,-c(1)] 
  
  # add a time-step (month 1-144) column
  import$month <- 1:144
  gst$month <- 1:144
  
  # need to transpose, keep as a df, and add column name
  import1<- gather(import, sci_name, qty, -month) #don't want to gather month
  gst1<- gather(gst, sci_name, avg, -month)    
  
  #join imports y gst
  imp_gst <- full_join(import1, gst1) # with full-join will automatically join on any identically-named columns
  # so don't need to specify "by = ""; only need the by if didn't want to join on all of the columns.
  # also use by to join by columns with different names like sci_name y scientific name
  
  # remove the stupid dot
  imp_gst <- mutate(imp_gst, sci_name = gsub(".", " ", sci_name, fixed = TRUE)) #replacing the dot with a space
  
  # write y then reimport csv file after changing months (1-144) to dates
  #write.csv(imp_gst, "imp_gst_dates.csv", row.names = FALSE)
  imp_gst<-read.csv("imp_gst_dates.csv") #change if just want to limit the number of years
  #imp_gst<-read.csv("imp_gst_dates_10-15.csv") 
  #imp_gst<-read.csv("imp_gst_dates_sugar_glider.csv") 
  
  #make real date for R code
  imp_gst$date <- as.Date(imp_gst$date, format = "%m/%d/%Y")
  
  ###############
  # select species for plots  # CHANGE name depending if import or gst
  
  #sp_select <- read.csv("graph_gst_sp_herps.txt") #sp to keep; need header = false so doesn't put 1st sp as header OR add sci_name at top of list
  #sp_select <- read.csv("graph_movie_bird.txt")
  #sp_select <- read.csv("graph_reptile_expo.txt")
  #sp_select <- read.csv("graph_500_2010.txt")
  #sp_select <- read.csv("graph_gst_sp_except_mammal.txt")
  #sp_select <- read.csv("graph_HP_movies.txt")
  sp_select <- read.csv("graph_iucn_sp.txt")
  imp_gst <- semi_join(imp_gst, sp_select) # sci_name is the only identical column name in both df
  imp_gst <- subset(imp_gst, sci_name=="Erythrura gouldiae")
  #imp_gst <- subset(imp_gst, sci_name=="Pogona henrylawsoni")
  
  # rescale the gst values so they can go on the same plot as imports
  # imp_gst <- group_by(imp_gst, sci_name) %>%
  #     mutate(avg_scaled=(avg/max(avg)* max(qty)))
  # 
  # 
  # ggplot(imp_gst, aes(x = month)) +
  #             geom_line(aes(y=qty), col="blue") +
  #             geom_line(aes(y=avg_scaled), col="green") +
  #             facet_wrap("sci_name", ncol=2, scales="free") +
  #             scale_y_continuous(sec.axis=sec_axis(~./max(.), "avg"))
  
  
  #okay for small number of graphs on a single page
  # grid.arrange(grobs=grobs, ncol=2,
  #              bottom=textGrob("month", vjust=1),
  #              left = textGrob("Imported Individuals", rot = 90, vjust = 1),
  #              right = textGrob("Average GST", rot = 270, vjust = 1))
  
  
  ##############
  
  # for species with >500% change and started in/after 2010
  grobs <- list()
  for (sp in unique(imp_gst$sci_name)) {
    imp_gst_sp <- filter(imp_gst, sci_name==sp)
    max_qty <- max(imp_gst_sp$qty)
    max_avg <- max(imp_gst_sp$avg)
    imp_gst_sp <- mutate(imp_gst_sp, avg_scaled=(avg * max_qty / max_avg))
    plt <- ggplot(imp_gst_sp, aes(date)) +
      geom_line(aes(y=qty), col='blue') +
      geom_line(aes(y=avg_scaled), col='green') +
      scale_x_date(NULL) +
      scale_y_continuous(NULL, sec.axis=sec_axis(~.*max_avg/max_qty)) +
      ggtitle(sp) +
      theme(plot.title = element_text(size = 8, face = "italic"))
    grobs <- c(grobs, list(plt))
  }
  
  
  
  return(out_dir)
}

# read files



##############

# for MOVIES
grobs <- list()
for (sp in unique(imp_gst$sci_name)) {
  imp_gst_sp <- filter(imp_gst, sci_name==sp)
  max_qty <- max(imp_gst_sp$qty)
  max_avg <- max(imp_gst_sp$avg)
  imp_gst_sp <- mutate(imp_gst_sp, avg_scaled=(avg * max_qty / max_avg))
    plt <- ggplot(imp_gst_sp, aes(date)) +
    geom_line(aes(y=qty), col='blue') +
    geom_line(aes(y=avg_scaled), col='green') +
    geom_vline(xintercept = as.numeric(as.Date("2011-07-15"))) +  # HP movie dates 
    geom_vline(xintercept = as.numeric(as.Date("2010-11-19"))) +
    geom_vline(xintercept = as.numeric(as.Date("2005-11-18"))) +   
    geom_vline(xintercept = as.numeric(as.Date("2009-07-15"))) +
    geom_vline(xintercept = as.numeric(as.Date("2007-07-11"))) +   
    geom_vline(xintercept = as.numeric(as.Date("2004-06-04"))) +
    scale_x_date(NULL) +
    scale_y_continuous(NULL, sec.axis=sec_axis(~.*max_avg/max_qty)) +
    ggtitle(sp) +
    theme(plot.title = element_text(size = 8, face = "italic"))
  grobs <- c(grobs, list(plt))
}

###############

# for US EXPOS (Daytona expo)
grobs <- list()
for (sp in unique(imp_gst$sci_name)) {
   imp_gst_sp <- filter(imp_gst, sci_name==sp)
    max_qty <- max(imp_gst_sp$qty)
    max_avg <- max(imp_gst_sp$avg)
    imp_gst_sp <- mutate(imp_gst_sp, avg_scaled=(avg * max_qty / max_avg))
    plt <- ggplot(imp_gst_sp, aes(date)) +
        geom_line(aes(y=qty), col='blue') +
        geom_line(aes(y=avg_scaled), col='green') +
       #geom_vline(xintercept = as.numeric(as.Date("2004-08-20"))) +  # reptile expo dates 
       #geom_vline(xintercept = as.numeric(as.Date("2005-08-20"))) +
       #geom_vline(xintercept = as.numeric(as.Date("2006-08-20"))) +  
       #geom_vline(xintercept = as.numeric(as.Date("2007-08-20"))) +
       #geom_vline(xintercept = as.numeric(as.Date("2008-08-20"))) +  
       #geom_vline(xintercept = as.numeric(as.Date("2009-08-20"))) +
        # geom_vline(xintercept = as.numeric(as.Date("2010-08-20"))) +   
        # geom_vline(xintercept = as.numeric(as.Date("2011-08-20"))) +
        # geom_vline(xintercept = as.numeric(as.Date("2012-08-20"))) +  
        # geom_vline(xintercept = as.numeric(as.Date("2013-08-20"))) +
        # geom_vline(xintercept = as.numeric(as.Date("2014-08-20"))) +  
        # geom_vline(xintercept = as.numeric(as.Date("2015-08-20"))) +
        scale_x_date(NULL) +
        scale_y_continuous(NULL, sec.axis=sec_axis(~.*max_avg/max_qty)) + #not sure why this isn't working???
        ggtitle(sp) +
        theme(plot.title = element_text(size = 8, face = "italic"))
    grobs <- c(grobs, list(plt))
}


###############

pdf(file='Erythrura gouldiae', width=7.5, height=10) #change name
plots_per_page <- 16
plots_this_page <- 1:plots_per_page
while (length(plots_this_page) == plots_per_page) {
  plots_this_page <- plots_this_page[plots_this_page <= length(grobs)]
  grid.arrange(grobs=grobs[plots_this_page], ncol=2, nrow=plots_per_page/2, #change # of columns if necessary (was 2)
               bottom=textGrob("Year", vjust=1),
               left = textGrob("Imported Individuals", rot = 90, vjust = 1),
               right = textGrob("Average OST", rot = 270, vjust = 1))
  plots_this_page <- plots_this_page + plots_per_page
}
dev.off()




#EXAMPLES

p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(aspect.ratio=1/1)
p1 <- p1 + theme(axis.text = element_text(size=rel(0.7)))
p1 <- p1 + labs (y = expression(Individuals~Imported)) # left y-axis
p1 <- p1 + labs (y = expression(Average~Google~Search~Trends)) # right y-axis
p1 <- p1 + labs(x = "Month (n=144)")
p1 <- p1 + theme(legend.title=element_blank())

#EXAMPLES


#p1
plot(import$Ameerega.hahneli, type="l",col="red")
plot(gst$Ameerega.hahneli, type="l",col="green")

#p2
plot(import$Calumma.oshaughnessyi, type="l",col="red")
plot(gst$Calumma.oshaughnessyi, type="l",col="green")

