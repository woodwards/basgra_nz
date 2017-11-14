## initialise_BASGRA_xxx.R ##
   cat(file=stderr(), 'Starting initialise_BASGRA_Scott.r', "\n")

## 1. GENERAL INITIALISATION ##
   dyn.load("BASGRA_WG.DLL") # use this one to use own PET
   source('scripts/initialise_BASGRA_general.R')

## 2. SITE CONDITIONS ##
   # Scott weather data is 2010-1 to 2016-366 (2557 rows)
   year_start     <- as.integer(2011)
   doy_start      <- as.integer(152) # 1 June
   year_stop      <- as.integer(2015)
   doy_stop       <- as.integer(151) # 30 May
   NDAYS          <- as.integer(difftime(
     as.Date(doy_stop -1, origin=paste(year_stop ,1,1,sep='-')), 
     as.Date(doy_start-1, origin=paste(year_start,1,1,sep='-')),
     units="days")) 
   #   NDAYS          <- as.integer(1460) # WARNING this value needs to fit within weather data!

   file_weather   <- 'model_inputs/weather_Scott.txt'
   file_params    <- 'model_inputs/parameters_Scott.txt' # can contain multiple columns
     parcol       <- 1 # which one are we going to use? (row names are ignored)
     
   # read harvest days (Simon) which are set up in make_weather1.r   
   days_harvest   <- as.matrix(read.table(file="model_inputs/harvest_Scott.txt", sep='\t'))
     
## 3. CREATE HARVEST CALENDAR AND WEATHER INPUT ##
   days_harvest   <- as.integer(days_harvest)
   matrix_weather <- read_weather_WG(year_start,doy_start,NDAYS,file_weather) # function in initialise_BASGRA_general.R
   
## 4. CREATE VECTOR "PARAMS" ##
   df_params      <- read.table(file_params,header=T,sep="\t",row.names=1)
   params         <- df_params[,parcol]
   
## 5. CREATE EMPTY MATRIX y ##
   y              <- matrix(0,NDAYS,NOUT)
   
   #
   cat(file=stderr(), 'Finished initialise_BASGRA_Scott.r', "\n")
   
