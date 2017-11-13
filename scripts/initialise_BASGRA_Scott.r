## initialise_BASGRA_Saerheim_00_early_Gri.R ##

## 1. GENERAL INITIALISATION ##
   dyn.load("BASGRA_WG.DLL") # use this one to use own PET
   source('scripts/initialise_BASGRA_general.R')

## 2. SITE CONDITIONS ##
   year_start     <- as.integer(2011)
   doy_start      <- as.integer(152) # 1 June
   NDAYS          <- as.integer(1460) # WARNING this value needs to fit within weather data!
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
