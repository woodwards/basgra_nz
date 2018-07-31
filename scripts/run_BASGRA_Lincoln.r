## 1. INITIALISATION ##
  library(tidyverse)
  cat(file=stderr(), 'Starting run_BASGRA_Lincoln.r', "\n")

  ## 1. GENERAL INITIALISATION ##
  cat(file=stderr(), 'Initialising BASGRA', "\n")
  source('scripts/initialise_BASGRA_general.R')
  BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
  dyn.load(BASGRA_DLL) 
  
  ## 2. SITE CONDITIONS ##
  cat(file=stderr(), 'Setting up simulation', "\n")
  
  # run model for this period 
  # Scott weather data is 2010-1 to 2017-318 (14 Nov)
  year_start     <- as.integer(2011)
  doy_start      <- as.integer(244) # 1 Sept
  year_stop      <- as.integer(2014)
  doy_stop       <- as.integer(151) # 31 May
  
  # calculate sim length --- WARNING this needs to fit within weather data!
  NDAYS          <- as.integer(difftime(
    as.Date(doy_stop -1, origin=paste(year_stop ,1,1,sep='-')), 
    as.Date(doy_start-1, origin=paste(year_start,1,1,sep='-')),
    units="days")) + 1L
  
  # specify weather file
  file_weather   <- 'model_inputs/weather_Lincoln.txt'

  # specify  parameter set
  # file_params    <- 'model_inputs/parameters_Scott.txt' # default parameters
  file_params    <- 'model_outputs/BASGRA_parModes_MAP1.txt' # optimised parameters
  parcol       <- 1 # which one are we going to use? (row names are ignored)

  # read harvest days (Simon) which are set up in make_weather1.r   
  days_harvest   <- as.matrix(read.table(file="model_inputs/harvest_Lincoln_3.txt", sep='\t'))
  
  ## 3. CREATE HARVEST CALENDAR AND WEATHER INPUT ##
  days_harvest   <- as.integer(days_harvest)
  matrix_weather <- read_weather_WG(year_start,doy_start,NDAYS,file_weather) 
  
  ## 4. CREATE VECTOR "PARAMS" ##
  df_params      <- read.table(file_params,header=T,sep="\t",row.names=1)
  params         <- df_params[,parcol]
  
  ## 5. CREATE EMPTY MATRIX y ##
  y              <- matrix(0,NDAYS,NOUT)
  
## 2. RUNNING ##
  cat(file=stderr(), 'Calling run_model()', "\n")
   output <- run_model()
  
## 3. OUTPUT ##
   cat(file=stderr(), 'Calling export_output()', "\n")
   export_output()

   # unload model   
   dyn.unload(BASGRA_DLL) 
   
   # prepare output for plotting
   output <- as.data.frame(output)
   names(output) <- outputNames
   x <- output$Time
   # plot(x,output$TILTOT)
   # plot(x,output$DEBUG)
   
   