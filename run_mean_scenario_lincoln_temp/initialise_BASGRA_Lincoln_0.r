## initialise_BASGRA_xxx.R ##

## 1. GENERAL INITIALISATION ##
source("scripts/initialise_BASGRA_general.R")

## 2. SITE CONDITIONS ##

# run model for this period 
# Scott weather data is 2010-1 to 2017-318 (14 Nov)
year_start     <- as.integer(2011)
# doy_start      <- as.integer(244) # 1 Sept
doy_start      <- as.integer(121) # 1 May
year_stop      <- as.integer(2017)
# doy_stop       <- as.integer(151) # 31 May
doy_stop       <- as.integer(120) # 30 April

# calculate sim length --- WARNING this needs to fit within weather data!
# NDAYS          <- as.integer(1460) 
NDAYS          <- as.integer(difftime(
  as.Date(doy_stop -1, origin=paste(year_stop ,1,1,sep="-")), 
  as.Date(doy_start-1, origin=paste(year_start,1,1,sep="-")),
  units="days")) + 1L

# use data in this period (replaced by data weights)
# year_start_data     <- as.integer(2012)
# doy_start_data      <- as.integer(1) # 1 Jan
# year_stop_data      <- as.integer(2017)
# doy_stop_data       <- as.integer(151) # 31 May

file_weather   <- paste(scenario, "/weather_Lincoln.txt", sep="")

file_params    <- paste(parameter_location, "/parameters_All.txt", sep="") # can contain multiple columns
parcol       <- 3 # which one are we going to use? (row names are ignored)

# read harvest days (Simon) which are set up in make_weather1.r
file_harvest   <- paste(scenario, "/harvest_Lincoln_0.txt", sep="")
days_harvest   <- as.matrix(read.table(file=file_harvest, sep="\t"))

## 3. CREATE HARVEST CALENDAR AND WEATHER INPUT ##
days_harvest   <- as.integer(days_harvest)
matrix_weather <- read_weather_WG(year_start,doy_start,NDAYS,file_weather) # function in initialise_BASGRA_general.R
matrix_weather[, 4:5] <- matrix_weather[, 4:5] - 0.0

## 4. CREATE VECTOR "PARAMS" ##
df_params      <- read.table(file_params,header=T,sep="\t",row.names=1)
params         <- df_params[,parcol]

## 5. CREATE EMPTY MATRIX y ##
y              <- matrix(0,NDAYS,NOUT)


